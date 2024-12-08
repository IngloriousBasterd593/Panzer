#ifndef PIPELINE
#define PIPELINE

#include "../utilities/common.h"

void lineBresenham(Mesh* mesh, unsigned int* frameColors, vec2f vertexStart, vec2f vertexEnd, unsigned int color) 
{
    int x1 = vertexStart.x;
    int y1 = vertexStart.y;
    int x2 = vertexEnd.x;
    int y2 = vertexEnd.y;

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    int stepX = (x1 < x2) ? 1 : -1;
    int stepY = (y1 < y2) ? 1 : -1;

    int error = dx - dy;
    int error2;

    while(1) 
    {
        // not out of bounds, trust me
        if(((y1 + mesh->Yposition) * SCREENWIDTH) + x1 + mesh->Xposition > SCREENSIZE || ((y1 + mesh->Yposition) * SCREENWIDTH) + x1 + mesh->Xposition < 0) 
        {
            break;
        }

        frameColors[((y1 + mesh->Yposition) * SCREENWIDTH) + x1 + mesh->Xposition] = color;

        if(x1 == x2 && y1 == y2) 
        {
            break;
        }

        error2 = 2 * error;

        if(error2 > -dy) 
        {
            error -= dy;
            x1 += stepX;
        }

        if(error2 < dx) 
        {
            error += dx;
            y1 += stepY;
        }
    }

    return;
}

void fillTriangle(Mesh* mesh, unsigned int* frameColors, vec2f baseVertex1, vec2f baseVertex2, vec2f centralVertex, int trianglePrecision, unsigned int color) 
{
    float dx1 = (centralVertex.x - baseVertex1.x) / trianglePrecision;
    float dy1 = (centralVertex.y - baseVertex1.y) / trianglePrecision;
    float dx2 = (centralVertex.x - baseVertex2.x) / trianglePrecision;
    float dy2 = (centralVertex.y - baseVertex2.y) / trianglePrecision;

    for (int i = 0; i < trianglePrecision; i++) 
    {
        baseVertex1.x += dx1;
        baseVertex1.y += dy1;

        baseVertex2.x += dx2;
        baseVertex2.y += dy2;

        lineBresenham(mesh, frameColors, baseVertex1, baseVertex2, color);
    }

    return;
}

void fillRectangle(Mesh* mesh, unsigned int* frameColors, vec2f vertexA, vec2f vertexB, vec2f vertexC, vec2f vertexD, int drawPrecision, unsigned int color) 
{
    float dxU = (vertexD.x - vertexC.x) / drawPrecision;
    float dyU = (vertexD.y - vertexC.y) / drawPrecision;
    float dxL = (vertexB.x - vertexA.x) / drawPrecision;
    float dyL = (vertexB.y - vertexA.y) / drawPrecision;

    for(int i = 0; i < drawPrecision; i++) 
    {
        vertexA.x += dxU;
        vertexA.y += dyU;

        vertexC.x += dxL;
        vertexC.y += dyL;

        lineBresenham(mesh, frameColors, vertexA, vertexC, color);
    }

    return;
}

int checkMeshBoundingBoxCollision(Mesh* mesh1, Mesh* mesh2)
{
    if(mesh1->boundingBox.xmax < mesh2->boundingBox.xmin || mesh1->boundingBox.xmin > mesh2->boundingBox.xmax) 
    {
        return false;
    }

    if(mesh1->boundingBox.ymax < mesh2->boundingBox.ymin || mesh1->boundingBox.ymin > mesh2->boundingBox.ymax) 
    {
        return false;
    }

    if(mesh1->boundingBox.zmax < mesh2->boundingBox.zmin || mesh1->boundingBox.zmin > mesh2->boundingBox.zmax) 
    {
        return false;
    }

    return true;
}

inline int checkBoundingBoxCollision(BoundingBox* b1, BoundingBox* b2)
{
   



    return true;
}

int checkIfMeshInRegion(Mesh* mesh, BoundingBox* region)
{
    if(mesh->pos.x > region.xmin && mesh->pos.x < region.xmax && mesh->pos.y > region.ymin && mesh->pos.y < region.ymax && mesh->pos.z > region.zmin && mesh->pos.z < region.zmax)
    {
        return true;
    } 

    return true;
}

int checkForMeshCollision(Mesh** meshes, int numberOfMeshes)
{









    return EXIT_SUCCESS;
}

void Mesh_draw(Mesh* mesh, vec3f* meshNormals, Camera* camera, unsigned int* frameColors, int drawPrecision) 
{
    if(drawPrecision < 1) 
    {
        fprintf(stderr, "write a correct drawprecision, bozo\n");
        return;
    }

    int index;
    int indexPlusPixelsPlusOne;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    int grayscaleRGB;
    float grayscaleCoefficient;
    unsigned int color;

    vec2f vertex1;
    vec2f vertex2;
    vec2f vertexUpper1;
    vec2f vertexUpper2;
    mat4f orthoMatrix;
    mat4f perspectiveMatrix;
    mat4f matrix;

       
/*
    perspectiveProjectionMatrix(&perspectiveMatrix, camera);

    for(int q = 0; q < PIXELS; q++) 
    {
        for(int l = 0; l < PIXELS; l++) 
        {
            index = q * PIXELS + l;

            vec4f vertex = { mesh->x[index] + mesh->Xposition, mesh->y[index] + mesh->Yposition, mesh->z[index] + mesh->Zposition, 1.0f };

            multiplyVectorByMatrix(&perspectiveMatrix, &vertex);
            vec3f projVertex = perspectiveNdcToScreen(mesh, &vertex);

            mesh->xProj[index] = projVertex.x;
            mesh->yProj[index] = projVertex.y;
            mesh->zProj[index] = projVertex.z;
        }
    } */


    mesh->pos.x += mesh->velocity.x;
    mesh->pos.y += mesh->velocity.x;
    mesh->pos.z += mesh->velocity.x;

    for(int q = 0; q < PIXELS; q++) 
    {
        for(int l = 0; l < PIXELS; l++) 
        {
            index = q * PIXELS + l;

            if(dotProduct3f(&camera->POV, &meshNormals[index]) < 0) {
                continue;
            }

            nextL = (l + 1) % PIXELS;
            nextQ = (q + 1) % PIXELS;

            indexPlusOne = q * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + l;
            indexPlusPixelsPlusOne = nextQ * PIXELS + nextL;

            grayscaleCoefficient = (dotProduct3f(&meshNormals[index], &camera->lightingDirectionVector));

            grayscaleRGB = (unsigned char) (255 - (132 * (1 - grayscaleCoefficient)));

            color = (0xFF << 24) | (grayscaleRGB << 16) | (grayscaleRGB << 8) | grayscaleRGB;

/*
            vertex1.x = mesh->xProj[index];
            vertex1.y = mesh->yProj[index];

            vertex2.x = mesh->xProj[indexPlusOne];
            vertex2.y = mesh->yProj[indexPlusOne];

            vertexUpper1.x = mesh->xProj[indexPlusPixels];
            vertexUpper1.y = mesh->yProj[indexPlusPixels];

            vertexUpper2.x = mesh->xProj[indexPlusPixelsPlusOne];
            vertexUpper2.y = mesh->yProj[indexPlusPixelsPlusOne];   */


            vertex1.x = mesh->x[index];
            vertex1.y = mesh->yProj[index];

            vertex2.x = mesh->x[indexPlusOne];
            vertex2.y = mesh->y[indexPlusOne];

            vertexUpper1.x = mesh->x[indexPlusPixels];
            vertexUpper1.y = mesh->y[indexPlusPixels];

            vertexUpper2.x = mesh->x[indexPlusPixelsPlusOne];
            vertexUpper2.y = mesh->y[indexPlusPixelsPlusOne];  

            fillTriangle(mesh, frameColors, vertex1, vertex2, vertexUpper1, drawPrecision, color);
            fillTriangle(mesh, frameColors, vertex2, vertexUpper1, vertexUpper2, drawPrecision, color);
        }
    }

    return;
}  


void Mesh_rotate(Mesh* mesh, vec3f* meshNormals, float rad, char axis) 
{
    vec3f previousVector;

    float cosRad = (float) cos(rad);
    float sinRad = (float) sin(rad);

    int index;

    switch(axis) 
    {
        case 'x':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {
                    index = j * PIXELS + k;

                    previousVector.x = mesh->x[index];
                    previousVector.y = mesh->y[index];
                    previousVector.z = mesh->z[index];

                    mesh->y[index] = previousVector.y * cosRad - previousVector.z * sinRad;
                    mesh->z[index] = previousVector.y * sinRad + previousVector.z * cosRad;

                    previousVector = meshNormals[index];

                    meshNormals[index].y = previousVector.y * cosRad - previousVector.z * sinRad;
                    meshNormals[index].z = previousVector.y * sinRad + previousVector.z * cosRad;
                }
            }
        break;

        case 'y':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {
                    index = j * PIXELS + k;

                    previousVector.x = mesh->x[index];
                    previousVector.y = mesh->y[index];
                    previousVector.z = mesh->z[index];

                    mesh->x[index] = previousVector.x * cosRad + previousVector.z * sinRad;
                    mesh->z[index] = -previousVector.x * sinRad + previousVector.z * cosRad;

                    previousVector = meshNormals[index];

                    meshNormals[index].x = previousVector.x * cosRad + previousVector.z * sinRad;
                    meshNormals[index].z = -previousVector.x * sinRad + previousVector.z * cosRad;
                }
            }
            break;

        case 'z':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {
                    index = j * PIXELS + k;

                    previousVector.x = mesh->x[index];
                    previousVector.y = mesh->y[index];
                    previousVector.z = mesh->z[index];

                    mesh->x[index] = previousVector.x * cosRad - previousVector.y * sinRad;  
                    mesh->y[index] = previousVector.x * sinRad + previousVector.y * cosRad;

                    previousVector = meshNormals[index];

                    meshNormals[index].x = previousVector.x * cosRad - previousVector.y * sinRad;  
                    meshNormals[index].y = previousVector.x * sinRad + previousVector.y * cosRad; 
                }
            }
            break;

        default:
        break;
    }

    return;
}


#endif