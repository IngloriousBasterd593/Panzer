#ifndef MATHCORE
#define MATHCORE

#include "../utilities/common.h"
#include "../utilities/shared.h"

vec3f unit3f(vec3f* v) 
{
    float magnitude = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);

    if(magnitude == 0.0f) 
    {
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

float maxOfArray(float* array) 
{
    float maxValue = array[0];
    for (int i = 1; i < VERTICES; i++) 
    {
        if (array[i] > maxValue) 
        {
            maxValue = array[i];
        }
    }
    return maxValue;
}

float maxOfTwoF(float num1, float num2)
{
    if(num1 > num2)
    {
        return num1;
    } else if(num1 < num2)
    {
        return num2;
    } else
    {
        return num1;
    }
}

float minOfTwoF(float num1, float num2)
{
    if(num1 > num2)
    {
        return num2;
    } else if(num1 < num2)
    {
        return num1;
    } else
    {
        return num1;
    }
}

void generateBoundingBox(Mesh* mesh)
{
    mesh->boundingBox.xmax = mesh->x[0] + mesh->Xposition;
    mesh->boundingBox.ymax = mesh->y[0] + mesh->Yposition;
    mesh->boundingBox.zmax = mesh->z[0] + mesh->Zposition;
    mesh->boundingBox.xmin = mesh->x[0] + mesh->Xposition;
    mesh->boundingBox.ymin = mesh->y[0] + mesh->Yposition;
    mesh->boundingBox.zmin = mesh->z[0] + mesh->Zposition;

    for(int i = 1; i < VERTICES; i++)
    {
        if(mesh->boundingBox.xmax < mesh->x[i] + mesh->Xposition)
        {
            mesh->boundingBox.xmax = mesh->x[i] + mesh->Xposition;
        }
        if(mesh->boundingBox.xmin > mesh->x[i] + mesh->Xposition)
        {
            mesh->boundingBox.xmin = mesh->x[i] + mesh->Xposition;
        }

        if(mesh->boundingBox.ymax < mesh->y[i] + mesh->Yposition)
        {
            mesh->boundingBox.ymax = mesh->y[i] + mesh->Yposition;
        }
        if(mesh->boundingBox.ymin > mesh->y[i] + mesh->Yposition)
        {
            mesh->boundingBox.ymin = mesh->y[i] + mesh->Yposition;
        }

        if(mesh->boundingBox.zmax < mesh->z[i] + mesh->Zposition)
        {
            mesh->boundingBox.zmax = mesh->z[i] + mesh->Zposition;
        }
        if(mesh->boundingBox.zmin > mesh->z[i] + mesh->Zposition)
        {
            mesh->boundingBox.zmin = mesh->z[i] + mesh->Zposition;
        }
    }

    int index;
    int indexPlusPixelsPlusOne;
    int nextJ;
    int nextI;

    for(int i = 0; i < BOUNDINGBOXCOUNT; i += BOUNDINGBOXSTEP)
    {
        for(int j = 0; j < BOUNDINGBOXCOUNT; j += BOUNDINGBOXSTEP)
        {
            index = i * PIXELS + j;
            nextJ = (j + BOUNDINGBOXSTEP) % PIXELS;
            nextI = (i + BOUNDINGBOXSTEP) % PIXELS;
            indexPlusPixelsPlusOne = nextI * PIXELS + nextJ;

            mesh->boundingBoxes[index].xmax = maxOfTwoF(mesh->x[index], mesh->x[indexPlusPixelsPlusOne]);
            mesh->boundingBoxes[index].xmin = minOfTwoF(mesh->x[index], mesh->x[indexPlusPixelsPlusOne]);   
            mesh->boundingBoxes[index].ymax = maxOfTwoF(mesh->y[index], mesh->y[indexPlusPixelsPlusOne]);
            mesh->boundingBoxes[index].ymin = minOfTwoF(mesh->y[index], mesh->y[indexPlusPixelsPlusOne]);   
            mesh->boundingBoxes[index].zmax = maxOfTwoF(mesh->z[index], mesh->z[indexPlusPixelsPlusOne]);
            mesh->boundingBoxes[index].zmin = minOfTwoF(mesh->z[index], mesh->z[indexPlusPixelsPlusOne]);
        }
    }

    return;
}

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
    for (int i = 0; i < 4; i++) 
    {
        for (int j = 0; j < 4; j++) 
        {
            resultMatrix->raw[i * 4 + j] = 0;
            for (int k = 0; k < 4; k++) 
            {
                resultMatrix->raw[i * 4 + j] += m1->raw[i * 4 + k] * m2->raw[k * 4 + j];
            }
        }
    }
}

vec3f perspectiveNdcToScreen(Mesh* mesh, vec4f* vertex)
{
    return (vec3f) { 
        (((vertex->x / vertex->w) + 1.0f) * HALFWINWIDTH) + mesh->Xposition,
        ((1.0f - (vertex->y / vertex->w)) * HALFWINHEIGHT) + mesh->Yposition,
        vertex->z / vertex->w
    };
}

void sphere_init(Mesh* mesh, vec3f* meshNormals, int radius, int offsetX, int offsetY, int offsetZ) 
{
    if(radius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    mesh->Xposition = offsetX;
    mesh->Yposition = offsetY;
    mesh->Zposition = offsetZ;
    
    int index, nextL, nextQ;

    vec3f normalVector, partialDerivativeU, partialDerivativeV;
    
    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            mesh->x[index] = radius * cos(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j);
            mesh->y[index] = radius * sin(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j);
            mesh->z[index] = radius * cos(PIOVERPIXELS * j);
        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;

            index = j * PIXELS + k; 

            partialDerivativeU.x = mesh->x[index] - mesh->x[index + nextL];
            partialDerivativeU.y = mesh->y[index] - mesh->y[index + nextL];
            partialDerivativeU.z = mesh->z[index] - mesh->z[index + nextL];

            partialDerivativeV.x = mesh->x[index] - mesh->x[nextQ * PIXELS + k];
            partialDerivativeV.y = mesh->y[index] - mesh->y[nextQ * PIXELS + k];
            partialDerivativeV.z = mesh->z[index] - mesh->z[nextQ * PIXELS + k];
            
            normalVector = crossProduct(&partialDerivativeU, &partialDerivativeV);

            meshNormals[index] = unit3f(&normalVector);
        }
    }

    generateBoundingBox(mesh);
    mesh->velocity = zerovector3i;

    return;
}

void torus_init(Mesh* mesh, vec3f* meshNormals, int innerRadius, int outerRadius, int offsetX, int offsetY, int offsetZ) 
{
    if(innerRadius <= 10 || outerRadius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    mesh->Xposition = offsetX;
    mesh->Yposition = offsetY;
    mesh->Zposition = offsetZ;

    vec3f normalVector, partialDerivativeU, partialDerivativeV;

    int index, nextL, nextQ;

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 

            mesh->x[index] = (outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j) + offsetX;
            mesh->y[index] = (outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j) + offsetY;
            mesh->z[index] = innerRadius * sin(TWOPIOVERPIXELS * k) + offsetZ;
        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;

            index = j * PIXELS + k; 

            partialDerivativeU.x = mesh->x[index] - mesh->x[index + nextL];
            partialDerivativeU.y = mesh->y[index] - mesh->y[index + nextL];
            partialDerivativeU.z = mesh->z[index] - mesh->z[index + nextL];

            partialDerivativeV.x = mesh->x[index] - mesh->x[nextQ * PIXELS + k];
            partialDerivativeV.y = mesh->y[index] - mesh->y[nextQ * PIXELS + k];
            partialDerivativeV.z = mesh->z[index] - mesh->z[nextQ * PIXELS + k];
            
            normalVector = crossProduct(&partialDerivativeU, &partialDerivativeV);

            meshNormals[index] = unit3f(&normalVector);
        }
    }

    generateBoundingBox(mesh);

    mesh->velocity = zerovector3i;

    return;
} 

#endif
