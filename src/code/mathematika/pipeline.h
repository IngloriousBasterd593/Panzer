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

        frameColors[((y1 + mesh->pos.y + mesh->Yposition) * SCREENWIDTH) + x1 + mesh->pos.y + mesh->Xposition] = color;

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


// check for mesh bounding box colission
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
    if (b1->xmin > b2->xmax || b1->xmax < b2->xmin) return false;
    if (b1->ymin > b2->ymax || b1->ymax < b2->ymin) return false;
    if (b1->zmin > b2->zmax || b1->zmax < b2->zmin) return false;
    
    return true;
}

// check if mesh is in region
int checkIfMeshInRegion(Mesh* mesh, BoundingBox* region)
{
    if(mesh->pos.x > region.xmin && mesh->pos.x < region.xmax && mesh->pos.y > region.ymin && mesh->pos.y < region.ymax && mesh->pos.z > region.zmin && mesh->pos.z < region.zmax)
    {
        return true;
    } 

    return true;
}

int checkForMeshCollisionAndUpdateMeshParameters(Mesh** meshes, int numberOfMeshes) 
{
    for(int i = 0; i < numberOfMeshes; i++) 
    {
        for(int j = i + 1; j < numberOfMeshes; j++) 
        {
           
            if(checkMeshBoundingBoxCollision(meshes[i], meshes[j])) 
            {
                for(int bi = 0; bi < BOUNDINGBOXCOUNT; bi++) 
                {
                    for(int bj = 0; bj < BOUNDINGBOXCOUNT; bj++) 
                    {
                        if(checkBoundingBoxCollision(&meshes[i]->boundingBoxes[bi], &meshes[j]->boundingBoxes[bj])) 
                        {
                            
                            
                            // Collision detected - update parameters

                            printf("Collision detected");
                            
                            return EXIT_SUCCESS;
                        }
                    }
                }
            }
        }
    }

    return EXIT_SUCCESS;
}

void Mesh_draw(Mesh** meshes, int numberOfMeshes, Camera* camera, unsigned int* frameColors, int drawPrecision)
{
    checkForMeshCollisionAndUpdateMeshParameters(meshes, numberOfMeshes);
    
    for (int m = 0; m < numberOfMeshes; m++)
    {
        meshes[m]->pos.x += meshes[m]->velocity.x;
        meshes[m]->pos.y += meshes[m]->velocity.x;
        meshes[m]->pos.z += meshes[m]->velocity.x;

        for(int q = 0; q < PIXELS; q++) 
        {
            for(int l = 0; l < PIXELS; l++) 
            {
                int index = q * PIXELS + l;

                if(dotProduct3f(&camera->POV, &meshes[m]->meshNormals[index]) < 0) {
                    continue;
                }

                int nextL = (l + 1) % PIXELS;
                int nextQ = (q + 1) % PIXELS;

                int indexPlusOne = q * PIXELS + nextL;
                int indexPlusPixels = nextQ  * PIXELS + l;
                int indexPlusPixelsPlusOne = nextQ * PIXELS + nextL;

                float grayscaleCoefficient = dotProduct3f(&meshes[m]->meshNormals[index], &camera->lightingDirectionVector);

                unsigned char grayscaleRGB = (unsigned char)(255 - (132 * (1 - grayscaleCoefficient)));
                unsigned int color = (0xFF << 24) | (grayscaleRGB << 16) 
                                                   | (grayscaleRGB << 8) 
                                                   | grayscaleRGB;

                /*
                // If you have projected coordinates, you could use these:
                vertex1.x = meshes[m]->xProj[index];
                vertex1.y = meshes[m]->yProj[index];

                vertex2.x = meshes[m]->xProj[indexPlusOne];
                vertex2.y = meshes[m]->yProj[indexPlusOne];

                vertexUpper1.x = meshes[m]->xProj[indexPlusPixels];
                vertexUpper1.y = meshes[m]->yProj[indexPlusPixels];

                vertexUpper2.x = meshes[m]->xProj[indexPlusPixelsPlusOne];
                vertexUpper2.y = meshes[m]->yProj[indexPlusPixelsPlusOne];
                */

                // Otherwise, using “raw” coordinates mixed with yProj:
                vec2f vertex1;
                vertex1.x = meshes[m]->x[index];
                vertex1.y = meshes[m]->yProj[index];

                vec2f vertex2;
                vertex2.x = meshes[m]->x[indexPlusOne];
                vertex2.y = meshes[m]->y[indexPlusOne];

                vec2f vertexUpper1;
                vertexUpper1.x = meshes[m]->x[indexPlusPixels];
                vertexUpper1.y = meshes[m]->y[indexPlusPixels];

                vec2f vertexUpper2;
                vertexUpper2.x = meshes[m]->x[indexPlusPixelsPlusOne];
                vertexUpper2.y = meshes[m]->y[indexPlusPixelsPlusOne];  

                fillTriangle(meshes[m], frameColors, vertex1, vertex2, vertexUpper1, drawPrecision, color);
                fillTriangle(meshes[m], frameColors, vertex2, vertexUpper1, vertexUpper2, drawPrecision, color);
            }
        }
    }

    return;
}



void Mesh_rotate(Mesh** meshes, int numberOfMeshes, float rad, char axis)
{
    float cosRad = (float)cos(rad);
    float sinRad = (float)sin(rad);

    for (int m = 0; m < numberOfMeshes; m++)
    {
        // For each mesh in the array
        Mesh* currentMesh = meshes[m]; 
        vec3f previousVector;
        int index;

        switch (axis) 
        {
            case 'x':
            {
                for (int j = 0; j < PIXELS; j++) 
                {
                    for (int k = 0; k < PIXELS; k++) 
                    {
                        index = j * PIXELS + k;

                        // Rotate vertex
                        previousVector.x = currentMesh->x[index];
                        previousVector.y = currentMesh->y[index];
                        previousVector.z = currentMesh->z[index];

                        currentMesh->y[index] = previousVector.y * cosRad - previousVector.z * sinRad;
                        currentMesh->z[index] = previousVector.y * sinRad + previousVector.z * cosRad;

                        // Rotate normal
                        previousVector = currentMesh->meshNormals[index];
                        currentMesh->meshNormals[index].y = previousVector.y * cosRad - previousVector.z * sinRad;
                        currentMesh->meshNormals[index].z = previousVector.y * sinRad + previousVector.z * cosRad;
                    }
                }
                break;
            }
            case 'y':
            {
                for (int j = 0; j < PIXELS; j++) 
                {
                    for (int k = 0; k < PIXELS; k++) 
                    {
                        index = j * PIXELS + k;

                        // Rotate vertex
                        previousVector.x = currentMesh->x[index];
                        previousVector.y = currentMesh->y[index];
                        previousVector.z = currentMesh->z[index];

                        currentMesh->x[index] = previousVector.x * cosRad + previousVector.z * sinRad;
                        currentMesh->z[index] = -previousVector.x * sinRad + previousVector.z * cosRad;

                        // Rotate normal
                        previousVector = currentMesh->meshNormals[index];
                        currentMesh->meshNormals[index].x = previousVector.x * cosRad + previousVector.z * sinRad;
                        currentMesh->meshNormals[index].z = -previousVector.x * sinRad + previousVector.z * cosRad;
                    }
                }
                break;
            }
            case 'z':
            {
                for (int j = 0; j < PIXELS; j++) 
                {
                    for (int k = 0; k < PIXELS; k++) 
                    {
                        index = j * PIXELS + k;

                        // Rotate vertex
                        previousVector.x = currentMesh->x[index];
                        previousVector.y = currentMesh->y[index];
                        previousVector.z = currentMesh->z[index];

                        currentMesh->x[index] = previousVector.x * cosRad - previousVector.y * sinRad;
                        currentMesh->y[index] = previousVector.x * sinRad + previousVector.y * cosRad;

                        // Rotate normal
                        previousVector = currentMesh->meshNormals[index];
                        currentMesh->meshNormals[index].x = previousVector.x * cosRad - previousVector.y * sinRad;  
                        currentMesh->meshNormals[index].y = previousVector.x * sinRad + previousVector.y * cosRad;
                    }
                }
                break;
            }
            default:
                // Do nothing if the axis is not x, y, or z
                break;
        }
    }
}



#endif