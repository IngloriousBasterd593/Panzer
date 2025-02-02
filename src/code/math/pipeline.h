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

void checkForMeshCollisionAndUpdateMeshParameters(Scene* sceneInstance) 
{

    
   

    return;
}

void Mesh_draw(Scene* sceneInstance)
{
    checkForMeshCollisionAndUpdateMeshParameters(sceneInstance->meshes, sceneInstance->numberOfMeshes);
    
    for(int m = 0; m < sceneInstance->numberOfMeshes; m++)
    {
        sceneInstance->meshes[m]->pos.x += sceneInstance->meshes[m]->velocity.x;
        sceneInstance->meshes[m]->pos.y += sceneInstance->meshes[m]->velocity.x;
        sceneInstance->meshes[m]->pos.z += sceneInstance->meshes[m]->velocity.x;

        for(int q = 0; q < PIXELS; q++) 
        {
            for(int l = 0; l < PIXELS; l++) 
            {
                int index = q * PIXELS + l;

                if(dotProduct3f(sceneInstance->camera->POV, sceneInstance->meshes[m]->meshNormals[index]) < 0) {
                    continue;
                }

                int nextL = (l + 1) % PIXELS;
                int nextQ = (q + 1) % PIXELS;

                int indexPlusOne = q * PIXELS + nextL;
                int indexPlusPixels = nextQ  * PIXELS + l;
                int indexPlusPixelsPlusOne = nextQ * PIXELS + nextL;

                float grayscaleCoefficient = dotProduct3f(sceneInstance->meshes[m]->meshNormals[index], sceneInstance->camera->lightingDirectionVector);

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
                vertex1.x = sceneInstance->meshes[m]->x[index];
                vertex1.y = sceneInstance->meshes[m]->yProj[index];

                vec2f vertex2;
                vertex2.x = sceneInstance->meshes[m]->x[indexPlusOne];
                vertex2.y = sceneInstance->meshes[m]->y[indexPlusOne];

                vec2f vertexUpper1;
                vertexUpper1.x = sceneInstance->meshes[m]->x[indexPlusPixels];
                vertexUpper1.y = sceneInstance->meshes[m]->y[indexPlusPixels];

                vec2f vertexUpper2;
                vertexUpper2.x = sceneInstance->meshes[m]->x[indexPlusPixelsPlusOne];
                vertexUpper2.y = sceneInstance->meshes[m]->y[indexPlusPixelsPlusOne];  

                fillTriangle(sceneInstance->meshes[m], sceneInstance->frameColors, vertex1, vertex2, vertexUpper1, sceneInstance->drawPrecision, color);
                fillTriangle(sceneInstance->meshes[m], sceneInstance->frameColors, vertex2, vertexUpper1, vertexUpper2, sceneInstance->drawPrecision, color);
            }
        }
    }

    return;
}


void Mesh_rotate(Mesh* mesh, float rad, vec3f axis)
{
    float axisLen = sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);

    axis.x /= axisLen;
    axis.y /= axisLen;
    axis.z /= axisLen;
    
    float cosRad = cos(rad);
    float sinRad = sin(rad);
    float oneMinusCos = 1.0f - cosRad;

        for (int j = 0; j < PIXELS; j++)
        {
            for (int k = 0; k < PIXELS; k++)
            {
                int index = j * PIXELS + k;
  
                float vx = mesh->x[index];
                float vy = mesh->y[index];
                float vz = mesh->z[index];
                
                float crossX = axis.y * vz - axis.z * vy;
                float crossY = axis.z * vx - axis.x * vz;
                float crossZ = axis.x * vy - axis.y * vx;
                
                float dot = axis.x * vx + axis.y * vy + axis.z * vz;
                
                float newVx = vx * cosRad + crossX * sinRad + axis.x * dot * oneMinusCos;
                float newVy = vy * cosRad + crossY * sinRad + axis.y * dot * oneMinusCos;
                float newVz = vz * cosRad + crossZ * sinRad + axis.z * dot * oneMinusCos;
                
                mesh->x[index] = newVx;
                mesh->y[index] = newVy;
                mesh->z[index] = newVz;



                
                vec3f n = mesh->meshNormals[index];

                float nx = n.x;
                float ny = n.y;
                float nz = n.z;
                
                float crossNX = axis.y * nz - axis.z * ny;
                float crossNY = axis.z * nx - axis.x * nz;
                float crossNZ = axis.x * ny - axis.y * nx;
                
                float dotN = axis.x * nx + axis.y * ny + axis.z * nz;
                
                vec3f n_rot;

                n_rot.x = nx * cosRad + crossNX * sinRad + axis.x * dotN * oneMinusCos;
                n_rot.y = ny * cosRad + crossNY * sinRad + axis.y * dotN * oneMinusCos;
                n_rot.z = nz * cosRad + crossNZ * sinRad + axis.z * dotN * oneMinusCos;
                
                mesh->meshNormals[index] = n_rot;




                traverseBVH(mesh->head, );
            }
        }
    }


#endif