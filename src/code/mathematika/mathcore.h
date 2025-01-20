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

float maxOfArrayF(float* array) 
{
    float maxValue = array[0];
    for(int i = 1; i < VERTICES; i++) 
    {
        if(array[i] > maxValue) 
        {
            maxValue = array[i];
        }
    }
    return maxValue;
}

float minOfArrayF(float* array) 
{
    float minValue = array[0];
    for(int i = 1; i < VERTICES; i++) 
    {
        if(array[i] < minValue) 
        {
            minValue = array[i];
        }
    }
    return minValue;
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
    mesh->boundingBox.xmax = maxOfArrayF(mesh->x);
    mesh->boundingBox.ymax = minOfArrayF(mesh->x);
    mesh->boundingBox.zmax = maxOfArrayF(mesh->y);
    mesh->boundingBox.xmin = minOfArrayF(mesh->y);
    mesh->boundingBox.ymin = maxOfArrayF(mesh->z);
    mesh->boundingBox.zmin = minOfArrayF(mesh->z);

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
        (((vertex->x / vertex->w) + 1.0f) * HALFWINWIDTH) + mesh->pos.x,
        ((1.0f - (vertex->y / vertex->w)) * HALFWINHEIGHT) + mesh->pos.y,
        vertex->z / vertex->w
    };
}

// a function to traverse an octree
void traverseOctree(OctreeNode* head) 
{
    if(head == NULL) 
    {
        return;
    }

    for(int i = 0; i < 8; i++) 
    {
        partitionBoundingBox(head->children[i]->boundingBoxes, head->boundingBox);
        traverseOctree(head->children[i]);
    }
}

void partitionBoundingBox(AABB* boundingBoxes[8], AABB* boundingBox) 
{
    float xmid = (boundingBox->xmax + boundingBox->xmin) / 2;
    float ymid = (boundingBox->ymax + boundingBox->ymin) / 2;
    float zmid = (boundingBox->zmax + boundingBox->zmin) / 2;

    for (int i = 0; i < 8; i++)
    {
        boundingBoxes[i] = (AABB*) malloc(sizeof(AABB));
    }

    boundingBoxes[0]->xmin = boundingBox->xmin;
    boundingBoxes[0]->xmax = xmid;
    boundingBoxes[0]->ymin = boundingBox->ymin;
    boundingBoxes[0]->ymax = ymid;
    boundingBoxes[0]->zmin = boundingBox->zmin;
    boundingBoxes[0]->zmax = zmid;

    boundingBoxes[1]->xmin = xmid;
    boundingBoxes[1]->xmax = boundingBox->xmax;
    boundingBoxes[1]->ymin = boundingBox->ymin;
    boundingBoxes[1]->ymax = ymid;
    boundingBoxes[1]->zmin = boundingBox->zmin;
    boundingBoxes[1]->zmax = zmid;

    boundingBoxes[2]->xmin = boundingBox->xmin;
    boundingBoxes[2]->xmax = xmid;
    boundingBoxes[2]->ymin = ymid;
    boundingBoxes[2]->ymax = boundingBox->ymax;
    boundingBoxes[2]->zmin = boundingBox->zmin;
    boundingBoxes[2]->zmax = zmid;

    boundingBoxes[3]->xmin = xmid;
    boundingBoxes[3]->xmax = boundingBox->xmax;
    boundingBoxes[3]->ymin = ymid;
    boundingBoxes[3]->ymax = boundingBox->ymax;
    boundingBoxes[3]->zmin = boundingBox->zmin;
    boundingBoxes[3]->zmax = zmid;
    
    boundingBoxes[4]->xmin = boundingBox->xmin;
    boundingBoxes[4]->xmax = xmid;
    boundingBoxes[4]->ymin = boundingBox->ymin;
    boundingBoxes[4]->ymax = ymid;
    boundingBoxes[4]->zmin = zmid;
    boundingBoxes[4]->zmax = boundingBox->zmax;

    boundingBoxes[5]->xmin = xmid;
    boundingBoxes[5]->xmax = boundingBox->xmax;
    boundingBoxes[5]->ymin = boundingBox->ymin;
    boundingBoxes[5]->ymax = ymid;
    boundingBoxes[5]->zmin = zmid;
    boundingBoxes[5]->zmax = boundingBox->zmax;

    boundingBoxes[6]->xmin = boundingBox->xmin;
    boundingBoxes[6]->xmax = xmid;
    boundingBoxes[6]->ymin = ymid;
    boundingBoxes[6]->ymax = boundingBox->ymax;
    boundingBoxes[6]->zmin = zmid;
    boundingBoxes[6]->zmax = boundingBox->zmax;

    boundingBoxes[7]->xmin = xmid;
    boundingBoxes[7]->xmax = boundingBox->xmax;
    boundingBoxes[7]->ymin = ymid;
    boundingBoxes[7]->ymax = boundingBox->ymax;
    boundingBoxes[7]->zmin = zmid;
    boundingBoxes[7]->zmax = boundingBox->zmax;

    return;
}



void sphere_init(Mesh* mesh, int radius, int offsetX, int offsetY, int offsetZ) 
{
    if(radius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;
    
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

            mesh->meshNormals[index] = unit3f(&normalVector);
        }
    }

    generateBoundingBox(mesh);
    mesh->velocity = zerovector3i;

    return;
}

void torus_init(Mesh* mesh, int innerRadius, int outerRadius, int offsetX, int offsetY, int offsetZ) 
{
    if(innerRadius <= 10 || outerRadius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;

    vec3f normalVector, partialDerivativeU, partialDerivativeV;

    int index, nextL, nextQ;

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 

            mesh->x[index] = (outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j);
            mesh->y[index] = (outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j);
            mesh->z[index] = innerRadius * sin(TWOPIOVERPIXELS * k);
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

            mesh->meshNormals[index] = unit3f(&normalVector);
        }
    }

    generateBoundingBox(mesh);

    mesh->velocity = zerovector3i;

    return;
} 

#endif
