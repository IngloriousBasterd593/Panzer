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

float magnitude3fS(vec3f* v) 
{
    return v->x * v->x + v->y * v->y + v->z * v->z;
}

float maxOfArrayF(float* array, int lower, int upper) 
{
    float maxValue = array[lower];
    for(int i = lower + 1; i < upper; i++) 
    {
        if(array[i] > maxValue) 
        {
            maxValue = array[i];
        }
    }
    return maxValue;
}

float minOfArrayF(float* array, int lower, int upper) 
{
    float minValue = array[lower];
    for(int i = lower + 1; i < upper; i++) 
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

AABB generateBoundingBox3D(int* x, int* y, int* z, int l, int u)
{
    AABB AABB;
    AABB.xmax = maxOfArrayF(x, l, u);
    AABB.ymax = minOfArrayF(x, l, u);
    AABB.zmax = maxOfArrayF(y, l, u);
    AABB.xmin = minOfArrayF(y, l, u);
    AABB.ymin = maxOfArrayF(z, l, u);
    AABB.zmin = minOfArrayF(z, l, u);

    return AABB;
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

unsigned int UDFS(vec3f* v1, vec3f* v2)
{
    return (v2->x - v1->x) * (v2->x - v1->x) + (v2->y - v1->y) * (v2->y - v1->y) + (v2->z - v1->z) * (v2->z - v1->z);
}

int checkBoundingBoxCollision(AABB* b1, AABB* b2) 
{
    if (b1->xmin > b2->xmax || b1->xmax < b2->xmin)
        return false;
    if (b1->ymin > b2->ymax || b1->ymax < b2->ymin) 
        return false;
    if (b1->zmin > b2->zmax || b1->zmax < b2->zmin) 
        return false;
    
    return true;
}

/*
void traverseOctree(Mesh* mesh, (void (*callback) (Mesh*))) 
{
    if(mesh->head == NULL) 
    {
        return;
    }

    for(int i = 0; i < 8; i++) 
    {
        callback(mesh);
        traverseOctree(head->children[i]);
    }
}

void createMeshOctree(Mesh* mesh)
{
    for(int i = 0; i < 8; i++)
    {
        head->children[i] = (OctreeNode*) malloc(sizeof(OctreeNode));
        if(head->children[i] == NULL)
        {
            for(int j = 0; j < i; j++)
            {
                free(head->children[j]);
            }

            fprintf(stderr, "Failed to allocate memory for octree child\n");
            return;
        }

        head->children[i]->children = NULL;
        head->children[i]->boundingBox = NULL;
    }

    return;
}

void partitionMeshBoundingBox(Mesh* mesh) 
{
    if (mesh->head == NULL)
    {
        return;
    }

    float xmid = (mesh->head->boundingBox->xmax + mesh->head->boundingBox->xmin) / 2;
    float ymid = (mesh->head->boundingBox->ymax + mesh->head->boundingBox->ymin) / 2;
    float zmid = (mesh->head->boundingBox->zmax + mesh->head->boundingBox->zmin) / 2;

    for (int i = 0; i < 8; i++)
    {
        mesh->head->children[i]->boundingBox = (AABB*) malloc(sizeof(AABB));
        if (mesh->head->children[i]->boundingBox == NULL)
        {
            for(int j = 0; j < i; j++)
            {
                free(mesh->head->children[j]->boundingBox);
            }
            fprintf(stderr, "Failed to allocate memory for boundingBox\n");
            return;
        }
    }

    mesh->head->children[0]->boundingBox->xmin = mesh->head->boundingBox->xmin;
    mesh->head->children[0]->boundingBox->xmax = xmid;
    mesh->head->children[0]->boundingBox->ymin = mesh->head->boundingBox->ymin;
    mesh->head->children[0]->boundingBox->ymax = ymid;
    mesh->head->children[0]->boundingBox->zmin = mesh->head->boundingBox->zmin;
    mesh->head->children[0]->boundingBox->zmax = zmid;

    mesh->head->children[1]->boundingBox->xmin = xmid;
    mesh->head->children[1]->boundingBox->xmax = mesh->head->boundingBox->xmax;
    mesh->head->children[1]->boundingBox->ymin = mesh->head->boundingBox->ymin;
    mesh->head->children[1]->boundingBox->ymax = ymid;
    mesh->head->children[1]->boundingBox->zmin = mesh->head->boundingBox->zmin;
    mesh->head->children[1]->boundingBox->zmax = zmid;

    mesh->head->children[2]->boundingBox->xmin = mesh->head->boundingBox->xmin;
    mesh->head->children[2]->boundingBox->xmax = xmid;
    mesh->head->children[2]->boundingBox->ymin = ymid;
    mesh->head->children[2]->boundingBox->ymax = mesh->head->boundingBox->ymax;
    mesh->head->children[2]->boundingBox->zmin = mesh->head->boundingBox->zmin;
    mesh->head->children[2]->boundingBox->zmax = zmid;

    mesh->head->children[3]->boundingBox->xmin = xmid;
    mesh->head->children[3]->boundingBox->xmax = mesh->head->boundingBox->xmax;
    mesh->head->children[3]->boundingBox->ymin = ymid;
    mesh->head->children[3]->boundingBox->ymax = mesh->head->boundingBox->ymax;
    mesh->head->children[3]->boundingBox->zmin = mesh->head->boundingBox->zmin;
    mesh->head->children[3]->boundingBox->zmax = zmid;
    
    mesh->head->children[4]->boundingBox->xmin = mesh->head->boundingBox->xmin;
    mesh->head->children[4]->boundingBox->xmax = xmid;
    mesh->head->children[4]->boundingBox->ymin = mesh->head->boundingBox->ymin;
    mesh->head->children[4]->boundingBox->ymax = ymid;
    mesh->head->children[4]->boundingBox->zmin = zmid;
    mesh->head->children[4]->boundingBox->zmax = mesh->head->boundingBox->zmax;

    mesh->head->children[5]->boundingBox->xmin = xmid;
    mesh->head->children[5]->boundingBox->xmax = mesh->head->boundingBox->xmax;
    mesh->head->children[5]->boundingBox->ymin = mesh->head->boundingBox->ymin;
    mesh->head->children[5]->boundingBox->ymax = ymid;
    mesh->head->children[5]->boundingBox->zmin = zmid;
    mesh->head->children[5]->boundingBox->zmax = mesh->head->boundingBox->zmax;

    mesh->head->children[6]->boundingBox->xmin = mesh->head->boundingBox->xmin;
    mesh->head->children[6]->boundingBox->xmax = xmid;
    mesh->head->children[6]->boundingBox->ymin = ymid;
    mesh->head->children[6]->boundingBox->ymax = mesh->head->boundingBox->ymax;
    mesh->head->children[6]->boundingBox->zmin = zmid;
    mesh->head->children[6]->boundingBox->zmax = mesh->head->boundingBox->zmax;

    mesh->head->children[7]->boundingBox->xmin = xmid;
    mesh->head->children[7]->boundingBox->xmax = mesh->head->boundingBox->xmax;
    mesh->head->children[7]->boundingBox->ymin = ymid;
    mesh->head->children[7]->boundingBox->ymax = mesh->head->boundingBox->ymax;
    mesh->head->children[7]->boundingBox->zmin = zmid;
    mesh->head->children[7]->boundingBox->zmax = mesh->head->boundingBox->zmax;

    return;
}

void compareOctreeAABBs(OctreeNode* n1, OctreeNode* n2) 
{
    for(int i = 0; i < 8; i++)
    {
        for(int j = 0; j < 8; j++)
        {
            if(n1->children[i] == NULL || n2->children[j] == NULL)
            {
                
            }

            if(checkBoundingBoxCollision(n1->children[i]->boundingBox, n2->children[j]->boundingBox))
            {
                compareOctreeAABBs(n1->children[i], n2->children[j]);
            }
        }
    }

    return;
}

void partitionBoundingBox(AABB b, AABB* result[8]) 
{
    for(int i = 0; i < 8; i++)
    {
        if(result[i])
        {
            fprintf(stderr, "nullptr\n");
            return;
        }
    }

    float xmid = (b->xmax + b->xmin) / 2;
    float ymid = (b->ymax + b->ymin) / 2;
    float zmid = (b->zmax + b->zmin) / 2;

    result[1]->xmin = mesh->head->boundingBox->xmin;
    result[0]->xmax = xmid;
    result[0]->ymin = mesh->head->boundingBox->ymin;
    result[0]->ymax = ymid;
    result[0]->zmin = mesh->head->boundingBox->zmin;
    result[0]->zmax = zmid;

    result[1]->xmin = xmid;
    result[1]->xmax = mesh->head->boundingBox->xmax;
    result[1]->ymin = mesh->head->boundingBox->ymin;
    result[1]->ymax = ymid;
    result[1]->zmin = mesh->head->boundingBox->zmin;
    result[1]->zmax = zmid;

    result[2]->xmin = mesh->head->boundingBox->xmin;
    result[2]->xmax = xmid;
    result[2]->ymin = ymid;
    result[2]->ymax = mesh->head->boundingBox->ymax;
    result[2]->zmin = mesh->head->boundingBox->zmin;
    result[2]->zmax = zmid;

    result[3]->xmin = xmid;
    result[3]->xmax = mesh->head->boundingBox->xmax;
    result[3]->ymin = ymid;
    result[3]->ymax = mesh->head->boundingBox->ymax;
    result[3]->zmin = mesh->head->boundingBox->zmin;
    result[3]->zmax = zmid;
    
    result[4]->xmin = mesh->head->boundingBox->xmin;
    result[4]->xmax = xmid;
    result[4]->ymin = mesh->head->boundingBox->ymin;
    result[4]->ymax = ymid;
    result[4]->zmin = zmid;
    result[4]->zmax = mesh->head->boundingBox->zmax;

    result[5]->xmin = xmid;
    result[5]->xmax = mesh->head->boundingBox->xmax;
    result[5]->ymin = mesh->head->boundingBox->ymin;
    result[5]->ymax = ymid;
    result[5]->zmin = zmid;
    result[5]->zmax = mesh->head->boundingBox->zmax;

    result[6]->xmin = mesh->head->boundingBox->xmin;
    result[6]->xmax = xmid;
    result[6]->ymin = ymid;
    result[6]->ymax = mesh->head->boundingBox->ymax;
    result[6]->zmin = zmid;
    result[6]->zmax = mesh->head->boundingBox->zmax;

    result[7]->xmin = xmid;
    result[7]->xmax = mesh->head->boundingBox->xmax;
    result[7]->ymin = ymid;
    result[7]->ymax = mesh->head->boundingBox->ymax;
    result[7]->zmin = zmid;
    result[7]->zmax = mesh->head->boundingBox->zmax;

    return;
}

*/

void traverseBVH(TreeNode* node, void (*callback) (Mesh*))
{
    if(node == NULL) 
    {
        return;
    }

    callback(node);

    traverseBVH(node->children[0], callback);
    traverseBVH(node->children[1], callback);
}

void createBoundingBoxForNode(Mesh* mesh, TreeNode* node, int start, int end)
{
    int mid = (start + end) / 2;

    int l1 = start;
    int u1 = mid;

    int l2 = mid;
    int u2 = end;

    memcpy(node->children[0], generateBoundingBox3D(mesh->x, mesh->y, mesh->z, l1, u1), sizeof(AABB));
    memcpy(node->children[1], generateBoundingBox3D(mesh->x, mesh->y, mesh->z, l2, u2), sizeof(AABB));

    createBoundingBoxForNode(mesh, node->children[0], l1, u1);
    createBoundingBoxForNode(mesh, node->children[1], l2, u2);

    return;
}

void initializeBVHTree(TreeNode* head, int currentDepth)
{
    head->boundingBox = (AABB*) malloc(sizeof(AABB));
    if (head->boundingBox[0] == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory for bounding box\n");
        return;
    }

    if(currentDepth < BVH_DEPTH) 
    {
        for (int i = 0; i < 2; i++) 
        {
            head->children[i] = (TreeNode*)malloc(sizeof(TreeNode));
            if (!head->children[i]) 
            {
                for (int j = 0; j < i; j++) 
                {
                    free(head->children[j]);
                    head->children[j] = NULL;
                }

                fprintf(stderr, "Failed to allocate memory for child\n");
                return;
            }

            initializeBVHTree(head->children[i], currentDepth + 1);
        }
    } 
    else 
    {
        head->children[0] = NULL;
        head->children[1] = NULL;
    }
}


int compareBVHAABBs(TreeNode* n1, TreeNode* n2, int depth)
{
    if(n1 == NULL || n2 == NULL || n1->boundingBox == NULL || n2->boundingBox == NULL) 
    {
        return false;
    }

    if(checkBoundingBoxCollision(n1->boundingBox, n2->boundingBox)) 
    {
        if(depth == BVH_DEPTH) 
        {
            // collision detected
            return true;
        }

        for(int i = 0; i < 2; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                compareBVHAABBs(n1->children[i], n2->children[j], ++depth);
            }
        }
    }

    return false;
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

    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh, mesh->head, 0, VERTICES);
    
    mesh->velocity = zerovector3i;

    return;
} 

#endif
