#ifndef MATHCORE
#define MATHCORE

#include "../utilities/common.h"
#include "../utilities/shared.h"

vec3f add3f(vec3f* v1, vec3f* v2)
{
    return (vec3f) {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

vec3f sub3f(vec3f* v1, vec3f* v2)
{
    return (vec3f) {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

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

unsigned int UDFS(vec3f* v1, vec3f* v2)
{
    return (v2->x - v1->x) * (v2->x - v1->x) + (v2->y - v1->y) * (v2->y - v1->y) + (v2->z - v1->z) * (v2->z - v1->z);
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

void multiplyVectorByMatrix4f(mat4f* matrix, vec4f* vertex) 
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

void multiplyVectorByMatrix3f(mat3f* matrix, vec3f* vertex) 
{
    float x = vertex->x;
    float y = vertex->y;
    float z = vertex->z;

    vertex->x = matrix->column[0].x * x + matrix->column[1].x * y + matrix->column[2].x * z + matrix->column[3].x * w;
    vertex->y = matrix->column[0].y * x + matrix->column[1].y * y + matrix->column[2].y * z + matrix->column[3].y * w;
    vertex->z = matrix->column[0].z * x + matrix->column[1].z * y + matrix->column[2].z * z + matrix->column[3].z * w;

    return;
}

void multiplyMatrixByMatrix4f(mat4f* m1, mat4f* m2, mat4f* resultMatrix)
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





/*    
*
*
*
*
*
*   brow wtf
*   not needed
*
*
*
*/


void traverseOctree(Mesh* mesh) 
{
    if(mesh->spatialPartitioningHead == NULL) 
    {
        return;
    }

    for(int i = 0; i < 8; i++) 
    {
        traverseOctree(head->children[i]);
    }
}

void traverseAndRotateBVH(TreeNode* node)
{
    if(node == NULL) 
    {
        return;
    }

    

    traverseBVH(node->children[0]);
    traverseBVH(node->children[1]);
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

/*    
*
*
*
*
*
*   not needeed
*
*
*
*
*/

// generates bounding box from binary partitioned arrays for BVH inirialisation 
AABB generateBoundingBoxFromPartitionedArray(int* x, int* y, int* z, int l, int u)
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

// chech if two bounding boxes are colliding
int checkBoundingBoxCollision(AABB* b1, AABB* b2) 
{
    if (b1->xmin > b2->xmax || b1->xmax < b2->xmin)
    { 
    return false;
    }
    if (b1->ymin > b2->ymax || b1->ymax < b2->ymin) 
    {
    return false;
    }
    if (b1->zmin > b2->zmax || b1->zmax < b2->zmin) 
    {
    return false;
    }
    
    return true;
}


// generates and populates an octree with partitioned bounding boxes for broad phase collision detection
void generateSimulationSpaceOctree(OctreeNode* node, int currentDepth, int maxDepth) 
{
    if(currentDepth == maxDepth) 
    {
        return;
    }

    AABB partitionedBoundingBoxes[8];

    float xmid = (node->boundingBox->xmax + node->boundingBox->xmin) / 2;
    float ymid = (node->boundingBox->ymax + node->boundingBox->ymin) / 2;
    float zmid = (node->boundingBox->zmax + node->boundingBox->zmin) / 2;

    partitionedBoundingBoxes[0]->xmin = node->boundingBox->xmin;
    partitionedBoundingBoxes[0]->xmax = xmid;
    partitionedBoundingBoxes[0]->ymin = node->boundingBox->ymin;
    partitionedBoundingBoxes[0]->ymax = ymid;
    partitionedBoundingBoxes[0]->zmin = node->boundingBox->zmin;
    partitionedBoundingBoxes[0]->zmax = zmid;

    partitionedBoundingBoxes[1]->xmin = xmid;
    partitionedBoundingBoxes[1]->xmax = node->boundingBox->xmax;
    partitionedBoundingBoxes[1]->ymin = node->boundingBox->ymin;
    partitionedBoundingBoxes[1]->ymax = ymid;
    partitionedBoundingBoxes[1]->zmin = node->boundingBox->zmin;
    partitionedBoundingBoxes[1]->zmax = zmid;

    partitionedBoundingBoxes[2]->xmin = node->boundingBox->xmin;
    partitionedBoundingBoxes[2]->xmax = xmid;
    partitionedBoundingBoxes[2]->ymin = ymid;
    partitionedBoundingBoxes[2]->ymax = node->boundingBox->ymax;
    partitionedBoundingBoxes[2]->zmin = node->boundingBox->zmin;
    partitionedBoundingBoxes[2]->zmax = zmid;

    partitionedBoundingBoxes[3]->xmin = xmid;
    partitionedBoundingBoxes[3]->xmax = node->boundingBox->xmax;
    partitionedBoundingBoxes[3]->ymin = ymid;
    partitionedBoundingBoxes[3]->ymax = node->boundingBox->ymax;
    partitionedBoundingBoxes[3]->zmin = node->boundingBox->zmin;
    partitionedBoundingBoxes[3]->zmax = zmid;
        
    partitionedBoundingBoxes[4]->xmin = node->boundingBox->xmin;
    partitionedBoundingBoxes[4]->xmax = xmid;
    partitionedBoundingBoxes[4]->ymin = node->boundingBox->ymin;
    partitionedBoundingBoxes[4]->ymax = ymid;
    partitionedBoundingBoxes[4]->zmin = zmid;
    partitionedBoundingBoxes[4]->zmax = node->boundingBox->zmax;

    partitionedBoundingBoxes[5]->xmin = xmid;
    partitionedBoundingBoxes[5]->xmax = node->boundingBox->xmax;
    partitionedBoundingBoxes[5]->ymin = node->boundingBox->ymin;
    partitionedBoundingBoxes[5]->ymax = ymid;
    partitionedBoundingBoxes[5]->zmin = zmid;
    partitionedBoundingBoxes[5]->zmax = node->boundingBox->zmax;

    partitionedBoundingBoxes[6]->xmin = node->boundingBox->xmin;
    partitionedBoundingBoxes[6]->xmax = xmid;
    partitionedBoundingBoxes[6]->ymin = ymid;
    partitionedBoundingBoxes[6]->ymax = node->boundingBox->ymax;
    partitionedBoundingBoxes[6]->zmin = zmid;
    partitionedBoundingBoxes[6]->zmax = node->boundingBox->zmax;

    partitionedBoundingBoxes[7]->xmin = xmid;
    partitionedBoundingBoxes[7]->xmax = node->boundingBox->xmax;
    partitionedBoundingBoxes[7]->ymin = ymid;
    partitionedBoundingBoxes[7]->ymax = node->boundingBox->ymax;
    partitionedBoundingBoxes[7]->zmin = zmid;
    partitionedBoundingBoxes[7]->zmax = node->boundingBox->zmax;

    for(int i = 0; i < 8; i++) 
    {
        node->children[i] = (OctreeNode*) malloc(sizeof(OctreeNode));
        if(node->children[i] == NULL) 
        {
            for(int j = 0; j < i; j++) 
            {
                free(node->children[j]);
            }

            fprintf(stderr, "Failed to allocate memory for children\n");
            return;
        }

        node->children[i]->boundingBox = (AABB*) malloc(sizeof(AABB));
        if(node->boundingBox == NULL) 
        {
            fprintf(stderr, "Failed to allocate memory for bounding box\n");
            return;
        }

        memcpy(node->children[i]->boundingBox, &partitionedBoundingBoxes[i], sizeof(AABB));
    }

    for(int i = 0; i < 8; i++) 
    {
        generateSimulationSpaceOctree(node->children[i], ++currentDepth, maxDepth);
    }

    return;
}

// traverses a mesh BVH tree and generates bounding boxes for each node with partitioned arrays
// non trivial to implement box generation and initialisation in a single function
// therefore this function is a helper function called after initializeBVHTree
void createBoundingBoxForNode(vertex_array* v, TreeNode* node, int start, int end)
{
    int mid = (start + end) / 2;

    int l1 = start;
    int u1 = mid;

    int l2 = mid;
    int u2 = end;

    memcpy(node->children[0], generateBoundingBox3D(v->x, v->y, v->z, l1, u1), sizeof(AABB));
    memcpy(node->children[1], generateBoundingBox3D(v->x, v->y, v->z, l2, u2), sizeof(AABB));

    createBoundingBoxForNode(v, node->children[0], l1, u1);
    createBoundingBoxForNode(v, node->children[1], l2, u2);

    return;
}

// initializes a BVH tree for mesh
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
            head->children[i] = (TreeNode*) malloc(sizeof(TreeNode));
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

// compares two BVH AABBs for collision detection
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

// helper function to partition the scene simulation space into an octree for broad phase collision detection
void partitionSceneSimulationSpace(Scene* sceneInstance) 
{
    vec3f min = {sceneInstance->simulationSpace.xmin, sceneInstance->simulationSpace.ymin, sceneInstance->simulationSpace.zmin};
    vec3f max = {sceneInstance->simulationSpace.xmax, sceneInstance->simulationSpace.ymax, sceneInstance->simulationSpace.zmax};

    int SPATIALPARTITIONINGDEPTH = (int) (log2f((3 * sqrt(UDFS(&min, &max))) / (MAXMESHDISTANCE))); 

    sceneInstance->spatialPartitioningHead = (OctreeNode*) malloc(sizeof(OctreeNode));
    if(sceneInstance->spatialPartitioningHead == NULL) 
    {
        fprintf(stderr, "couldn't get space for octree");
        return;
    }

    sceneInstance->spatialPartitioningHead->boundingBox = (AABB*) malloc(sizeof(AABB));
    if(sceneInstance->spatialPartitioningHead->boundingBox == NULL) 
    {
        fprintf(stderr, "couldn't get space for bounding box");
        return;
    }

    memcpy(sceneInstance->spatialPartitioningHead->boundingBox, &sceneInstance->simulationSpace, sizeof(AABB));
    generateSimulationSpaceOctree(sceneInstance->spatialPartitioningHead, 0, SPATIALPARTITIONINGDEPTH);

    return;
}


void traverseBinaryTreeRoateAABBs(TreeNode* node, mat3f* rotationMatrix)
{
    if(!node)
    {
        return;
    }

    vec3f corners[6] =
    {
        {node->boundingBox->xmin, node->boundingBox->ymin, node->boundingBox->zmin},
        {node->boundingBox->xmax, node->boundingBox->ymin, node->boundingBox->zmin},
        {node->boundingBox->xmin, node->boundingBox->ymax, node->boundingBox->zmin},
        {node->boundingBox->xmax, node->boundingBox->ymax, node->boundingBox->zmin},
        {node->boundingBox->xmin, node->boundingBox->ymin, node->boundingBox->zmax},
        {node->boundingBox->xmax, node->boundingBox->ymax, node->boundingBox->zmax}
    };

    for(int i = 0; i < 6; i++)
    {
        vec3f corner = { corners[i].x, corners[i].y, corners[i].z };

        multiplyVectorByMatrix3f(rotationMatrix, &corner);

        corners[i].x = corner.x;
        corners[i].y = corner.y;
        corners[i].z = corner.z;
    }

    traverseBinaryTree(node->left, rotationMatrix);
    traverseBinaryTree(node->right, rotationMatrix);
}

void traverseBinaryTree

/*    
*
*
*
*
*
*   end
*
*
*
*
*/






void torus_init(Mesh* mesh, int innerRadius, int outerRadius, int offsetX, int offsetY, int offsetZ) 
{
    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;

    vec3f normalVector, partialDerivativeU, partialDerivativeV;

    int index, nextL, nextQ, i = 0, t_index = 0;

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            mesh->vert_array->x[index] = (outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j);
            mesh->vert_array->y[index] = (outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j);
            mesh->vert_array->z[index] = innerRadius * sin(TWOPIOVERPIXELS * k);  
        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;

            index = j * PIXELS + k; 

            partialDerivativeU.x = mesh->vert_array->x[index] - mesh->vert_array->x[index + nextL];
            partialDerivativeU.y = mesh->vert_array->y[index] - mesh->vert_array->y[index + nextL];
            partialDerivativeU.z = mesh->vert_array->z[index] - mesh->vert_array->z[index + nextL];

            partialDerivativeV.x = mesh->vert_array->x[index] - mesh->vert_array->x[nextQ * PIXELS + k];
            partialDerivativeV.y = mesh->vert_array->y[index] - mesh->vert_array->y[nextQ * PIXELS + k];
            partialDerivativeV.z = mesh->vert_array->z[index] - mesh->vert_array->z[nextQ * PIXELS + k];
            
            normalVector = crossProduct(&partialDerivativeU, &partialDerivativeV);

            mesh->triangles.normal[index] = unit3f(&normalVector);
        }
    }

    for(int i = 0; i < VERTICES; i += 3) 
    {
        mesh->vert_array.x[i + 2] = mesh->triangles[i].p3.z;
        mesh->vert_array.y[i + 2] = mesh->triangles[i].p3.y;
        mesh->vert_array.y[i + 2] = mesh->triangles[i].p3.x;
        mesh->vert_array.x[i + 1] = mesh->triangles[i].p2.z;
        mesh->vert_array.y[i + 1] = mesh->triangles[i].p2.y;
        mesh->vert_array.y[i + 1] = mesh->triangles[i].p2.x;
        mesh->vert_array.x[i] = mesh->triangles[i].p1.z;
        mesh->vert_array.y[i] = mesh->triangles[i].p1.y;
        mesh->vert_array.x[i] = mesh->triangles[i].p1.x;
    }

    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh->triangles, mesh->head, 0, VERTICES);
    
    mesh->velocity = zerovector3i;

    return;
} 

void sphere_init(Mesh* mesh, int radius, int offsetX, int offsetY, int offsetZ) 
{
    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;
    
    int index, nextL, nextQ, indexPlusOne, indexPlusPixels;
    vec3f normalVector, partialDerivativeU, partialDerivativeV;
    
    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            mesh->vert_array->x[index] = radius * cos(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j);
            mesh->vert_array->y[index] = radius * sin(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j);
            mesh->vert_array->z[index] = radius * cos(PIOVERPIXELS * j);
        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;

            index = j * PIXELS + k; 

            partialDerivativeU.x = mesh->vert_array->x[index] - mesh->vert_array->x[index + nextL];
            partialDerivativeU.y = mesh->vert_array->y[index] - mesh->vert_array->y[index + nextL];
            partialDerivativeU.z = mesh->vert_array->z[index] - mesh->vert_array->z[index + nextL];

            partialDerivativeV.x = mesh->vert_array->x[index] - mesh->vert_array->x[nextQ * PIXELS + k];
            partialDerivativeV.y = mesh->vert_array->y[index] - mesh->vert_array->y[nextQ * PIXELS + k];
            partialDerivativeV.z = mesh->vert_array->z[index] - mesh->vert_array->z[nextQ * PIXELS + k];
            
            normalVector = crossProduct(&partialDerivativeU, &partialDerivativeV);

            mesh->triangles.normal[index] = unit3f(&normalVector);
        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            indexPlusOne = (index + 1) % PIXELS;
            indexPlusOne = (index + PIXELS);
        }
    }

    for(int i = 0; i < VERTICES; i += 3) 
    {
        mesh->vert_array.x[i + 2] = mesh->triangles[i].p3.z;
        mesh->vert_array.y[i + 2] = mesh->triangles[i].p3.y;
        mesh->vert_array.y[i + 2] = mesh->triangles[i].p3.x;
        mesh->vert_array.x[i + 1] = mesh->triangles[i].p2.z;
        mesh->vert_array.y[i + 1] = mesh->triangles[i].p2.y;
        mesh->vert_array.y[i + 1] = mesh->triangles[i].p2.x;
        mesh->vert_array.x[i] = mesh->triangles[i].p1.z;
        mesh->vert_array.y[i] = mesh->triangles[i].p1.y;
        mesh->vert_array.x[i] = mesh->triangles[i].p1.x;
    }

    
    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh->triangles, mesh->head, 0, VERTICES);

    mesh->velocity = zerovector3i;

    return;
}

void mesh_init(Mesh* mesh, triangle* outTriangles, uint32_t outTriangleCount, int offsetX, int offsetY, int offsetZ) 
{
    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;
    
    for(int i = 0; i < VERTICES; i++) 
    {
        mesh->triangles[i] = outTriangles[i];   
    }

    for(int i = 0; i < VERTICES; i += 3) 
    {
        mesh->vert_array.x[i] = outTriangles[i].p1.x;
        mesh->vert_array.y[i] = outTriangles[i].p1.y;
        mesh->vert_array.x[i] = outTriangles[i].p1.z;
        mesh->vert_array.y[i + 1] = outTriangles[i].p2.x;
        mesh->vert_array.y[i + 1] = outTriangles[i].p2.y;
        mesh->vert_array.x[i + 1] = outTriangles[i].p2.z;
        mesh->vert_array.y[i + 2] = outTriangles[i].p3.x;
        mesh->vert_array.y[i + 2] = outTriangles[i].p3.y;
        mesh->vert_array.x[i + 2] = outTriangles[i].p3.z;
    }
    
    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh->triangles, mesh->head, 0, VERTICES);

    mesh->velocity = zerovector3i;

    return;
}

#endif
