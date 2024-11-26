#ifndef MATHCORE
#define MATHCORE

#include "../utilities/common.h"
#include "../utilities/shared.h"


vec3f unit3f(vec3f* v) 
{

    float magnitude = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);

    if(magnitude == 0.0f) 
    {
        // fprintf(stderr, "magnitude zero\n");
        return (vec3f) {0, 0, 0};
    } 

    return (vec3f) {v->x / magnitude, v->y / magnitude, v->z / magnitude};
}



vec3f crossProduct(vec3f* v1, vec3f* v2) 
{
    return (vec3f) {
                    v1->y * v2->z - v1->z * v2->y, 
                    v1->z * v2->x - v1->x * v2->z, 
                    v1->x * v2->y - v1->y * v2->x 
    };
}



float dotProduct3f(vec3f* v1, vec3f* v2) 
{
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}



float dotProduct2f(vec2f* v1, vec2f* v2) 
{
    return v1->x * v2->x + v1->y * v2->y;
}


/*
void OrthographicProjectionMatrix4f(mat4f* result, Camera* camera) 
{
    result->column[0] = (vec4f) { 2.0f / (right - left), 0.0f, 0.0f, 0.0f };
    result->column[1] = (vec4f) { 0.0f, 2.0f / (top - bottom), 0.0f, 0.0f };
    result->column[2] = (vec4f) { 0.0f, 0.0f, -2.0/(farZ-nearZ), 0.0f };
    result->column[3] = (vec4f) { -(right + left) / (right - left), -(top + bottom) / (top - bottom), -(farZ + nearZ) / (farZ - nearZ), 1.0f };

    return;
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
*/




void orthographicProjectionMatrix(mat4f* result, Camera* camera) 
{
    memset(result->raw, 0, 64);

    result->column[0] = (vec4f) { 2.0f / (camera->right - camera->left), 0.0f, 0.0f, 0.0f };
    result->column[1] = (vec4f) { 0.0f, 2.0f / (camera->top - camera->bottom), 0.0f, 0.0f };
    result->column[2] = (vec4f) { 0.0f, 0.0f, -2.0f / (camera->farPlane - camera->nearPlane), 0.0f };
    result->column[3] = (vec4f) { -(camera->right + camera->left) / (camera->right - camera->left), -(camera->top + camera->bottom) / (camera->top - camera->bottom), -(camera->farPlane + camera->nearPlane) / (camera->farPlane - camera->nearPlane), 1.0f };

    return;
}




void perspectiveProjectionMatrix(mat4f* result, Camera* camera) 
{

    memset(result->raw, 0, 64);
    
    float f = 1.0f / tanf(camera->FOV / 2.0f);
    float nf = 1.0f / (camera->nearPlane - camera->farPlane);

    result->column[0].x = f / camera->aspectRatio;
    result->column[1].y = f;
    result->column[2].z = (camera->farPlane + camera->nearPlane) * nf;
    result->column[2].w = 1.0f;
    result->column[3].z = (2.0f * camera->farPlane * camera->nearPlane) * nf;

    return;
}




void multiplyVectorByMatrix(mat4f* matrix, vec4f* vertex) 
{

    float x = vertex->x;
    float y = vertex->y;
    float z = vertex->z;
    float w = 1.0f;

    vertex->x = matrix->column[0].x * x + matrix->column[1].x * y + matrix->column[2].x * z + matrix->column[3].x * w;
    vertex->y = matrix->column[0].y * x + matrix->column[1].y * y + matrix->column[2].y * z + matrix->column[3].y * w;
    vertex->z = matrix->column[0].z * x + matrix->column[1].z * y + matrix->column[2].z * z + matrix->column[3].z * w;
    vertex->w = matrix->column[0].w * x + matrix->column[1].w * y + matrix->column[2].w * z + matrix->column[3].w * w;

    return;
}




void multiplyMatrixByMatrix(mat4f* m1, mat4f* m2, mat4f* resultMatrix)
{
   
    resultMatrix->column[0].x = m1->column[0].x * m2->column[0].x + m1->column[1].x * m2->column[0].y + m1->column[2].x * m2->column[0].z + m1->column[3].x * m2->column[0].w;
    resultMatrix->column[0].y = m1->column[0].y * m2->column[0].x + m1->column[1].y * m2->column[0].y + m1->column[2].y * m2->column[0].z + m1->column[3].y * m2->column[0].w;
    resultMatrix->column[0].z = m1->column[0].z * m2->column[0].x + m1->column[1].z * m2->column[0].y + m1->column[2].z * m2->column[0].z + m1->column[3].z * m2->column[0].w;
    resultMatrix->column[0].w = m1->column[0].w * m2->column[0].x + m1->column[1].w * m2->column[0].y + m1->column[2].w * m2->column[0].z + m1->column[3].w * m2->column[0].w;

    resultMatrix->column[1].x = m1->column[0].x * m2->column[1].x + m1->column[1].x * m2->column[1].y + m1->column[2].x * m2->column[1].z + m1->column[3].x * m2->column[1].w;
    resultMatrix->column[1].y = m1->column[0].y * m2->column[1].x + m1->column[1].y * m2->column[1].y + m1->column[2].y * m2->column[1].z + m1->column[3].y * m2->column[1].w;
    resultMatrix->column[1].z = m1->column[0].z * m2->column[1].x + m1->column[1].z * m2->column[1].y + m1->column[2].z * m2->column[1].z + m1->column[3].z * m2->column[1].w;
    resultMatrix->column[1].w = m1->column[0].w * m2->column[1].x + m1->column[1].w * m2->column[1].y + m1->column[2].w * m2->column[1].z + m1->column[3].w * m2->column[1].w;

    resultMatrix->column[2].x = m1->column[0].x * m2->column[2].x + m1->column[1].x * m2->column[2].y + m1->column[2].x * m2->column[2].z + m1->column[3].x * m2->column[2].w;
    resultMatrix->column[2].y = m1->column[0].y * m2->column[2].x + m1->column[1].y * m2->column[2].y + m1->column[2].y * m2->column[2].z + m1->column[3].y * m2->column[2].w;
    resultMatrix->column[2].z = m1->column[0].z * m2->column[2].x + m1->column[1].z * m2->column[2].y + m1->column[2].z * m2->column[2].z + m1->column[3].z * m2->column[2].w;
    resultMatrix->column[2].w = m1->column[0].w * m2->column[2].x + m1->column[1].w * m2->column[2].y + m1->column[2].w * m2->column[2].z + m1->column[3].w * m2->column[2].w;

    resultMatrix->column[3].x = m1->column[0].x * m2->column[3].x + m1->column[1].x * m2->column[3].y + m1->column[2].x * m2->column[3].z + m1->column[3].x * m2->column[3].w;
    resultMatrix->column[3].y = m1->column[0].y * m2->column[3].x + m1->column[1].y * m2->column[3].y + m1->column[2].y * m2->column[3].z + m1->column[3].y * m2->column[3].w;
    resultMatrix->column[3].z = m1->column[0].z * m2->column[3].x + m1->column[1].z * m2->column[3].y + m1->column[2].z * m2->column[3].z + m1->column[3].z * m2->column[3].w;
    resultMatrix->column[3].w = m1->column[0].w * m2->column[3].x + m1->column[1].w * m2->column[3].y + m1->column[2].w * m2->column[3].z + m1->column[3].w * m2->column[3].w;

    return;
}




void sphere_init(Manifold* manifold, vec3f* manifoldNormals, int radius, int offsetX, int offsetY, int offsetZ) 
{

    if(radius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    manifold->Xposition = offsetX;
    manifold->Yposition = offsetY;
    manifold->Zposition = offsetZ;
    
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
        
        normalVector = crossProduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit3f(&normalVector);

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

    manifold->Xposition = offsetX;
    manifold->Yposition = offsetY;
    manifold->Zposition = offsetZ;

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
        
        normalVector = crossProduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit3f(&normalVector);

        }
    }

    return;
} 

#endif