#ifndef MATHCORE
#define MATHCORE

#include "../utilities/common.h"

vec3f unit(vec3f* v) 
{

    float magnitude = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
    vec3f normalized;

    if(magnitude == 0.0f) 
    {
        // fprintf(stderr, "magnitude zero\n");
        normalized.x = normalized.y = normalized.z = 0;
        return normalized;
    } 

    normalized.x = v->x / magnitude;
    normalized.y = v->y / magnitude;
    normalized.z = v->z / magnitude;

    return normalized;
}



vec3f crossproduct(vec3f* v1, vec3f* v2) 
{
    return (vec3f) {
                    v1->y * v2->z - v1->z * v2->y, 
                    v1->z * v2->x - v1->x * v2->z, 
                    v1->x * v2->y - v1->y * v2->x 
    };
}



float dotproduct(vec3f* v1, vec3f* v2) 
{
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}



mat4f* OrthographicProjectionMatrix4f(float left, float right, float bottom, float top, float nearZ, float farZ, mat4f* result) // applies to the entire matrix
{
    result->vecRows[0] = (Vector4f){ 2.0/(right-left), 0.0, 0.0, -(right+left)/(right-left) };
    result->vecRows[1] = (Vector4f){ 0.0, 2.0/(top-bottom), 0.0, -(top+bottom)/(top-bottom) };
    result->vecRows[2] = (Vector4f){ 0.0, 0.0, -2.0/(farZ-nearZ), -(farZ+nearZ)/(farZ-nearZ) };
    result->vecRows[3] = (Vector4f){ 0.0, 0.0, 0.0, 1.0 };
    return result;
}


vec3f perspectiveProject(vec3f* vertex, Camera* camera) 
{
    float f = 1.0f / tanf(camera->FOV / 2.0f);
    
    vec4f projectedVector = {
                            (((vertex->x * f) / (camera->aspectRatio * vertex->z)) + 1.0f) * (SCREENWIDTH / 2.0f),
                            (1.0f - ((vertex->y * f) / vertex->z)) * (SCREENHEIGHT / 2.0f),
                            ((camera->farPlane + camera->nearPlane) * vertex->z + 2.0f * camera->farPlane * camera->nearPlane) / (camera->nearPlane - camera->farPlane), 
                            vertex->z
    };

    return (vec3f) {projectedVector.x, projectedVector.y, projectedVector.z};
}

void sphere_init(Manifold* manifold, vec3f* manifoldNormals, int radius, int offsetX, int offsetY, int offsetZ) 
{

    if(radius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    manifold->Xposition = HALFWINWIDTH;
    manifold->Yposition = HALFWINHEIGHT;
    manifold->Zposition = 2000;
    
    int index;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    vec3f normalVector;
    vec3f partialDerivativeU;
    vec3f partialDerivativeV;
    
    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        index = j * PIXELS + k; 
        
        manifold->x[index] = ((radius * cos(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j)));
        manifold->y[index] = ((radius * sin(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j)));
        manifold->z[index] = ((radius * cos(PIOVERPIXELS * j)));

        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        nextL = (k + 1) % PIXELS;
        nextQ = (j + 1) % PIXELS;

        index = j * PIXELS + k; 
        indexPlusOne = j * PIXELS + nextL;
        indexPlusPixels = nextQ * PIXELS + k;

        partialDerivativeU.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusOne]);
        partialDerivativeU.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusOne]);
        partialDerivativeU.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusOne]);

        partialDerivativeV.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusPixels]);
        partialDerivativeV.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusPixels]);
        partialDerivativeV.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusPixels]);
        
        normalVector = crossproduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit(&normalVector);

        }
    }

    return;
}  



void torus_init(Manifold* manifold, vec3f* manifoldNormals, int innerRadius, int outerRadius, int offsetX, int offsetY, int offsetZ) 
{

    if(innerRadius <= 10 || outerRadius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    manifold->Xposition = HALFWINWIDTH;
    manifold->Yposition = HALFWINHEIGHT;
    manifold->Zposition = 700;

    vec3f normalVector;
    vec3f partialDerivativeU;
    vec3f partialDerivativeV;

    int index;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        index = j * PIXELS + k; 

        manifold->x[index] = (((outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j)) + offsetX);
        manifold->y[index] = (((outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j)) + offsetY);
        manifold->z[index] = ((innerRadius * sin(TWOPIOVERPIXELS * k)) + offsetZ);

        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        nextL = (k + 1) % PIXELS;
        nextQ = (j + 1) % PIXELS;

        index = j * PIXELS + k; 
        indexPlusOne = j * PIXELS + nextL;
        indexPlusPixels = nextQ * PIXELS + k;

        partialDerivativeU.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusOne]);
        partialDerivativeU.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusOne]);
        partialDerivativeU.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusOne]);

        partialDerivativeV.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusPixels]);
        partialDerivativeV.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusPixels]);
        partialDerivativeV.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusPixels]);
        
        normalVector = crossproduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit(&normalVector);

        }
    }

    return;
} 

#endif