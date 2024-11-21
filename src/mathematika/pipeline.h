#ifndef COMMONH
#define COMMONH

#include "../utilities/common.h"



void lineBresenham(Manifold* manifold, unsigned int* frameColors, vec2f vertexStart, vec2f vertexEnd, unsigned int color) 
{
    
    short x1 = vertexStart.x;
    short y1 = vertexStart.y;
    short x2 = vertexEnd.x;
    short y2 = vertexEnd.y;

    short dx = abs(x2 - x1);
    short dy = abs(y2 - y1);

    short stepX = (x1 < x2) ? 1 : -1;
    short stepY = (y1 < y2) ? 1 : -1;

    short error = dx - dy;
    short error2;

    while(1) 
    {
        // not out of bounds, trust me
        if(((y1 + manifold->Yposition) * SCREENWIDTH) + x1 + manifold->Xposition > SCREENSIZE || ((y1 + manifold->Yposition) * SCREENWIDTH) + x1 + manifold->Xposition < 0) 
        {
            break;
        }

        frameColors[((y1 + manifold->Yposition) * SCREENWIDTH) + x1 + manifold->Xposition] = color;

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



void fillTriangle(Manifold* manifold, unsigned int* frameColors, vec2f baseVertex1, vec2f baseVertex2, vec2f centralVertex, int trianglePrecision, unsigned int color) 
{

    float dx1 = (centralVertex.x - baseVertex1.x) / trianglePrecision;
    float dy1 = (centralVertex.y - baseVertex1.y) / trianglePrecision;
    float dx2 = (centralVertex.x - baseVertex2.x) / trianglePrecision;
    float dy2 = (centralVertex.y - baseVertex2.y) / trianglePrecision;

    vec2f point1;
    vec2f point2;

    for (int i = 0; i < trianglePrecision; i++) 
    {

        point1.x = baseVertex1.x + dx1;
        point1.y = baseVertex1.y + dy1;

        point2.x = baseVertex2.x + dx2;
        point2.y = baseVertex2.y + dy2;

        lineBresenham(manifold, frameColors, point1, point2, color);

    }

    return;
}



void fillRectangle(Manifold* manifold, unsigned int* frameColors, vec2f vertexA, vec2f vertexB, vec2f vertexC, vec2f vertexD, int drawPrecision, unsigned int color) 
{
    /*
    if(drawPrecision <= 0) 
    {
        fprintf(stderr, "write a correct drawprecision, bozo\n");
        return;
    } */

    float dxU = (vertexD.x - vertexC.x) / drawPrecision;
    float dyU = (vertexD.y - vertexC.y) / drawPrecision;
    float dxL = (vertexB.x - vertexA.x) / drawPrecision;
    float dyL = (vertexB.y - vertexA.y) / drawPrecision;

    for(int i = 0; i <= drawPrecision; i++) 
    {

        vertexA.x += dxU;
        vertexA.y += dyU;

        vertexC.x += dxL;
        vertexC.y += dyL;

        lineBresenham(manifold, frameColors, vertexA, vertexC, color);
    }

    return;
}



void Manifold_draw(Manifold* manifold, vec3f* manifoldNormals, Camera* camera, unsigned int* frameColors, int drawPrecision) 
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
    float f = 1.0f / tanf(camera->FOV / 2.0f);

    vec2f vertex1;
    vec2f vertex2;
    vec2f vertexUpper1;
    vec2f vertexUpper2;

    // perform perspective projection
    for(int q = 0; q < PIXELS; q++) 
    {
        for(int l = 0; l < PIXELS; l++) 
        {
            index = q * PIXELS + l;

            /*
            if (manifold->z[index] <= 0.1f) 
            {
                // perror("lmao");
                continue; 
            }
            */

            vec3f projectedVertex = {((manifold->x[index] + manifold->Xposition) / manifold->z[index]) + HALFWINWIDTH, ((manifold->y[index] + manifold->Yposition) / (manifold->z[index] + manifold->Xposition)) + HALFWINHEIGHT, manifold->z[index]};
            // printf("%f %f %f\n", projectedVertex.x, projectedVertex.y, projectedVertex.z);



/*
            vec4f projectedVertex = {
                ((((manifold->x[index] + manifold->Xposition) * f) / (camera->aspectRatio * (manifold->z[index] + manifold->Zposition))) + 1.0f) * (SCREENWIDTH / 2.0f),
                ((1.0f - (((manifold->y[index] + manifold->Yposition) * f) / manifold->z[index])) * (SCREENHEIGHT / 2.0f)),
                ((camera->farPlane + camera->nearPlane) * (manifold->z[index] + manifold->Zposition) + 2.0f * camera->farPlane * camera->nearPlane) / (camera->nearPlane - camera->farPlane),
                1.0f
            };
*/ 
            manifold->xProj[index] = projectedVertex.x;
            manifold->yProj[index] = projectedVertex.y;
            manifold->zProj[index] = projectedVertex.z; 

/*
            manifold->xProj[index] = manifold->x[index];
            manifold->yProj[index] = manifold->y[index];
            manifold->zProj[index] = manifold->z[index]; */

            //printf("%f %f %f", manifold->x[index] + manifold->Xposition, manifold->y[index] + manifold->Yposition, manifold->z[index] + manifold->Zposition);
            //printf("\t%f %f %f\n", manifold->xProj[index], manifold->yProj[index], manifold->zProj[index]);
            //usleep(1000);

        }
    }


    for(int q = 0; q < PIXELS; q++) 
    {
        for(int l = 0; l < PIXELS; l++) 
        {
            index = q * PIXELS + l;

            if(dotproductf3(&camera->POV, &manifoldNormals[index]) < 0) {
                continue;
            }
            
            nextL = (l + 1) % PIXELS;
            nextQ = (q + 1) % PIXELS;

            indexPlusOne = q * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + l;
            indexPlusPixelsPlusOne = nextQ * PIXELS + nextL;

            grayscaleCoefficient = (dotproduct(&manifoldNormals[index], &camera->lightingDirectionVector));

            grayscaleRGB = (uint8_t) (255 - (132 * (1 - grayscaleCoefficient)));

            color = (0xFF << 24) | (grayscaleRGB << 16) | (grayscaleRGB << 8) | grayscaleRGB;

            vertex1.x = manifold->xProj[index];
            vertex1.y = manifold->yProj[index];

            vertex2.x = manifold->xProj[indexPlusOne];
            vertex2.y = manifold->yProj[indexPlusOne];

            vertexUpper1.x = manifold->xProj[indexPlusPixels];
            vertexUpper1.y = manifold->yProj[indexPlusPixels];

            vertexUpper2.x = manifold->xProj[indexPlusPixelsPlusOne];
            vertexUpper2.y = manifold->yProj[indexPlusPixelsPlusOne];   
             
/*
            vertex1.x = manifold->x[index];
            vertex1.y = manifold->y[index];

            vertex2.x = manifold->x[indexPlusOne];
            vertex2.y = manifold->y[indexPlusOne];

            vertexUpper1.x = manifold->x[indexPlusPixels];
            vertexUpper1.y = manifold->y[indexPlusPixels];

            vertexUpper2.x = manifold->x[indexPlusPixelsPlusOne];
            vertexUpper2.y = manifold->y[indexPlusPixelsPlusOne]; 
*/

            fillRectangle(manifold, frameColors, vertex1, vertex2, vertexUpper1, vertexUpper2, drawPrecision, color);

            // fillTriangle(frameColors, vertex1, vertex2, vertexUpper1, drawPrecision, color);
            // fillTriangle(frameColors, vertex2, vertexUpper1, vertexUpper2, drawPrecision, color);

        }
    }

    return;
}  



void Manifold_rotate(Manifold* manifold, vec3f* manifoldNormals, float rad, char axis) 
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
                    // Rotate manifold VERTICES
                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->y[index] = previousVector.y * cosRad - previousVector.z * sinRad;
                    manifold->z[index] = previousVector.y * sinRad + previousVector.z * cosRad;

                    // Rotate normal vectors

                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].y = previousVector.y * cosRad - previousVector.z * sinRad;
                    manifoldNormals[index].z = previousVector.y * sinRad + previousVector.z * cosRad;

                }
            }
        break;

        case 'y':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {
                    
                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->x[index] = previousVector.x * cosRad + previousVector.z * sinRad;
                    manifold->z[index] = -previousVector.x * sinRad + previousVector.z * cosRad;


                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].x = previousVector.x * cosRad + previousVector.z * sinRad;
                    manifoldNormals[index].z = -previousVector.x * sinRad + previousVector.z * cosRad;
                }
            }
            break;

        case 'z':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {

                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->x[index] = previousVector.x * cosRad - previousVector.y * sinRad;  
                    manifold->y[index] = previousVector.x * sinRad + previousVector.y * cosRad;


                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].x = previousVector.x * cosRad - previousVector.y * sinRad;  
                    manifoldNormals[index].y = previousVector.x * sinRad + previousVector.y * cosRad; 

                }
            }
            break;

        default:
        break;
    }

    return;
}

#endif