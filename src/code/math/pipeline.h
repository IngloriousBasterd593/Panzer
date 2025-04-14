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

    // deprecated

    /*

    while(1) 
    {
        if(((y1 + mesh->pos.y) * SCREENWIDTH) + x1 + mesh->pos.x > SCREENSIZE || ((y1 + mesh->pos.y) * SCREENWIDTH) + x1 + mesh->pos.x < 0) 
        {
            break;
        }

        frameColors[((y1 + mesh->pos.y) * SCREENWIDTH) + x1 + mesh->pos.x] = color;

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

    */


    while(1) 
    {
        if(x1 > SCREENWIDTH || x1 < 0 || y1 > SCREENHEIGHT || y2 < 0) 
        {
            break;
        }

        frameColors[(y1 * SCREENWIDTH) + x1] = color;

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

void fillTriangle(Mesh* mesh, unsigned int* frameColors, triangle2f* t, int trianglePrecision, unsigned int color) 
{
    float dx1 = (t->p3.x - t->p1.x) / trianglePrecision;
    float dy1 = (t->p3.y - t->p1;.y) / trianglePrecision;
    float dx2 = (t->p3.x - t->p2.x) / trianglePrecision;
    float dy2 = (t->p3.y - t->p2.y) / trianglePrecision;

    for (int i = 0; i < trianglePrecision; i++) 
    {
        t->p1;.x += dx1;
        t->p1;.y += dy1;

        t->p2.x += dx2;
        t->p2.y += dy2;

        lineBresenham(mesh, frameColors, t->p1;, t->p2, color);
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
    sceneInstance->meshes[m]->pos.x += sceneInstance->meshes[m]->velocity.x;
    sceneInstance->meshes[m]->pos.y += sceneInstance->meshes[m]->velocity.y;
    sceneInstance->meshes[m]->pos.z += sceneInstance->meshes[m]->velocity.z;

    for(int i = 0; i < sceneInstance->meshCount; i++)
    {
        for(int j = 0; j < sceneInstance->meshCount; i++)
        {
            if(i == j)
            {
                continue;
            }

            if(checkBoundingBoxCollision({sceneInstance->meshes[i]->meshMin, sceneInstance->meshes[i]->meshMin}, {sceneInstance->meshes[j]->meshMax, sceneInstance->meshes[j]->meshMax}))
            {
                compareBVHAABBs(sceneInstance->meshes[i]->head, sceneInstance->meshes[j]->head);
            }
        }   
    }

    return;
}

void Mesh_draw(Scene* sceneInstance)
{
    checkForMeshCollisionAndUpdateMeshParameters(sceneInstance);

    for (int m = 0; m < sceneInstance->meshCount; m++)
    {
        Mesh* mesh = sceneInstance->meshes[m];

        for (int i = 0; i < mesh->triangle_count; i++)
        {
            if (dotProduct3f(sceneInstance->camera.POV, mesh->triangles[i].normal) < 0)
            {
                continue;
            }

            float grayscaleCoefficient = dotProduct3f(mesh->triangles[i].normal, sceneInstance->camera.lightingDirectionVector);
            unsigned char grayscaleRGB = (unsigned char)(255 - (132 * (1 - grayscaleCoefficient)));
            unsigned int color = (0xFF << 24) | (grayscaleRGB << 16) 
                                 | (grayscaleRGB << 8) | grayscaleRGB;

            mat4f proj;
            perspective(sceneInstance->camera.FOV, sceneInstance->camera.aspectRatio, sceneInstance->camera.nearPlane, sceneInstance->camera.farPlane, &proj);

            vec4f v1 = { mesh->triangles[i].p1.x, mesh->triangles[i].p1.y, mesh->triangles[i].p1.z, 1.0f };
            vec4f v2 = { mesh->triangles[i].p2.x, mesh->triangles[i].p2.y, mesh->triangles[i].p2.z, 1.0f };
            vec4f v3 = { mesh->triangles[i].p3.x, mesh->triangles[i].p3.y, mesh->triangles[i].p3.z, 1.0f };

            vec4f p1_proj = multiplyMatrixVector4f(&proj, &v1);
            vec4f p2_proj = multiplyMatrixVector4f(&proj, &v2);
            vec4f p3_proj = multiplyMatrixVector4f(&proj, &v3);

            if (p1_proj.w != 0.0f)
            {
                p1_proj.x /= p1_proj.w;
                p1_proj.y /= p1_proj.w;
                p1_proj.z /= p1_proj.w;
            }
            if (p2_proj.w != 0.0f)
            {
                p2_proj.x /= p2_proj.w;
                p2_proj.y /= p2_proj.w;
                p2_proj.z /= p2_proj.w;
            }
            if (p3_proj.w != 0.0f)
            {
                p3_proj.x /= p3_proj.w;
                p3_proj.y /= p3_proj.w;
                p3_proj.z /= p3_proj.w;
            }

            int p1_screen_x = (int)((p1_proj.x + 1.0f) * HALFWINWIDTH);
            int p1_screen_y = (int)((1.0f - p1_proj.y) * HALFWINHEIGHT);
            int p2_screen_x = (int)((p2_proj.x + 1.0f) * HALFWINWIDTH);
            int p2_screen_y = (int)((1.0f - p2_proj.y) * HALFWINHEIGHT);
            int p3_screen_x = (int)((p3_proj.x + 1.0f) * HALFWINWIDTH);
            int p3_screen_y = (int)((1.0f - p3_proj.y) * HALFWINHEIGHT);

            triangle2f t =
            {
                { p1_screen_x, p1_screen_y },
                { p2_screen_x, p2_screen_y },
                { p3_screen_x, p3_screen_y }
            };

            fillTriangle(mesh, sceneInstance->frameColors, &t, sceneInstance->drawPrecision, color);
        }
    }
}

void Mesh_rotate(Mesh* mesh, vec3f omega)
{
    
}


#endif