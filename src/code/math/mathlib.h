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

void perspective(float fov, float aspect, float near, float far, mat4f* m) 
{
    float t = tanf(fov / 2.0f);
    m->mat[0][0] = 1.0f / (aspect * t);
    m->mat[0][1] = 0.0f;
    m->mat[0][2] = 0.0f;
    m->mat[0][3] = 0.0f;
    m->mat[1][0] = 0.0f;
    m->mat[1][1] = 1.0f / t;
    m->mat[1][2] = 0.0f;
    m->mat[1][3] = 0.0f;
    m->mat[2][0] = 0.0f;
    m->mat[2][1] = 0.0f;
    m->mat[2][2] = -(far + near) / (far - near);
    m->mat[2][3] = -1.0f;
    m->mat[3][0] = 0.0f;
    m->mat[3][1] = 0.0f;
    m->mat[3][2] = -(2.0f * far * near) / (far - near);
    m->mat[3][3] = 0.0f;
}

vec4f multiplyMatrixVector4f(const mat4f* m, const vec4f* v) 
{
    vec4f result;

    result.x = m->mat[0][0] * v->x + m->mat[0][1] * v->y + m->mat[0][2] * v->z + m->mat[0][3] * v->w;
    result.y = m->mat[1][0] * v->x + m->mat[1][1] * v->y + m->mat[1][2] * v->z + m->mat[1][3] * v->w;
    result.z = m->mat[2][0] * v->x + m->mat[2][1] * v->y + m->mat[2][2] * v->z + m->mat[2][3] * v->w;
    result.w = m->mat[3][0] * v->x + m->mat[3][1] * v->y + m->mat[3][2] * v->z + m->mat[3][3] * v->w;

    return result;
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
    return (vec3f) 
    { 
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
int compareBVHAABBs(TreeNode* n1, TreeNode* n2, Mesh* m1, Mesh* m2, int depth)
{
    if(n1 == NULL || n2 == NULL || n1->boundingBox == NULL || n2->boundingBox == NULL) 
    {
        return false;
    }

    if (checkBoundingBoxCollision({n1->boundingBox.xmin + m1->pos.x, n1->boundingBox.ymin + m1->pos.y, n1->boundingBox.zmin + m1->pos.z,
                                    n1->boundingBox.xmax + m1->pos.x, n1->boundingBox.ymax + m1->pos.y, n1->boundingBox.zmax + m1->pos.z},

                                {n2->boundingBox.xmin + m2->pos.x, n2->boundingBox.ymin + m2->pos.y, n2->boundingBox.zmin + m2->pos.z,
                                    n2->boundingBox.xmax + m2->pos.x, n2->boundingBox.ymax + m2->pos.y, n2->boundingBox.zmax + m2->pos.z})) 
    {

        if(depth == BVH_DEPTH) 
        {
            resolveElasticCollision(m1, m2);

            // collision detected
            return true;
        }

        for(int i = 0; i < 2; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                compareBVHAABBs(n1->children[i], n2->children[j], m1, m2, depth + 1);
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



void resolveElasticCollision(Mesh* m1, Mesh* m2)
{






    return;
}

void mesh_inertia(Mesh* mesh, double rho, mat3f *I)
{
    vec3f g = {0, 0, 0};

    for(int i = 0; i < VERTICES; i++) 
    {
        g.x += mesh->vert_array[i].x;
        g.y += mesh->vert_array[i].y;
        g.z += mesh->vert_array[i].z;
    }

    g = v_scale(g, 1.0 / VERTICES);

    mat3f I_origin = m_zero();
    double M = 0.0;
    vec3f Mc = {0, 0, 0};

    for(int t = 0; t < mesh->triangle_count; t++)
    {
        vec3f p1 = mesh->triangles[t].p1;
        vec3f p2 = mesh->triangles[t].p2;
        vec3f p3 = mesh->triangles[t].p3;

        vec3f a = v_sub(p2, p1);
        vec3f b = v_sub(p3, p1);
        vec3f c = v_sub(g , p1);

        double V = fabs(v_dot(a,v_cross(b,c))) / 6.0;
        double m = rho * V;  
        M += m;

        vec3f cent = v_scale(v_add(v_add(p1,p2),v_add(p3, g)), 0.25);
        Mc = v_add(Mc, v_scale(cent,m));

        double x[4] = {p1.x, p2.x, p3.x, g.x};
        double y[4] = {p1.y, p2.y, p3.y, g.y};
        double z[4] = {p1.z, p2.z, p3.z, g.z};

        double Ixx = m * 0.1 * (S2(y)+S2(z));
        double Iyy =  m * 0.1 * (S2(z) + S2(x));
        double Izz =  m * 0.1 * (S2(x) + S2(y));
        double Ixy = -m * 0.05 * SP(x, y);
        double Iyz = -m * 0.05 * SP(y, z);
        double Izx = -m * 0.05 * SP(z, x);

        mat3f Io = {{{ Ixx, Ixy, Izx },
                    { Ixy, Iyy, Iyz },
                    { Izx, Iyz, Izz }}};

        I_origin = m_add(I_origin, Io);

        vec3f C = v_scale( Mc, 1.0 / M );

        double C2 = v_dot( C, C );
        mat3f shift = m_add(m_identity(M * C2), m_scale(outer(C, C), -M));

        *I = m_add(I_origin, m_scale(shift, -1));
    }
}






void torus_init(Mesh* mesh, int innerRadius, int outerRadius, int offsetX, int offsetY, int offsetZ) 
{
    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;

    int index, nextL, nextQ, indexPlusOne, indexPlusPixels, indexPlusPixelsPlusOne, triangle_index;
    vec3f normalVector, partialDerivativeU, partialDerivativeV;

    mesh->triangle_count = VERTICES * 2;

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

            mesh->triangles[triangle_index].normal = unit3f(&normalVector);
            triangle_index++;
            mesh->triangles[triangle_index].normal = unit3f(&normalVector);
            triangle_index++;
        }
    }

    triangle_index = 0;

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;
            indexPlusOne = j * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + k;
            indexPlusPixelsPlusOne = nextQ * PIXELS + k + 1;

            mesh->triangles[triangle_index].p1 = { mesh->vert_array.x[index],  mesh->vert_array.y[index],  mesh->vert_array.z[index]};
            mesh->triangles[triangle_index].p2 = { mesh->vert_array.x[indexPlusOne],  mesh->vert_array.y[indexPlusOne],  mesh->vert_array.z[indexPlusOne]};
            mesh->triangles[triangle_index].p3 = { mesh->vert_array.x[indexPlusPixels],  mesh->vert_array.y[indexPlusPixels],  mesh->vert_array.z[indexPlusPixels]};

            triangle_index++;

            mesh->triangles[triangle_index].p1 = { mesh->vert_array.x[indexPlusOne],  mesh->vert_array.y[indexPlusOne],  mesh->vert_array.z[indexPlusOne]};
            mesh->triangles[triangle_index].p2 = { mesh->vert_array.x[indexPlusPixels],  mesh->vert_array.y[indexPlusPixels],  mesh->vert_array.z[indexPlusPixels]};
            mesh->triangles[triangle_index].p3 = { mesh->vert_array.x[indexPlusPixelsPlusOne],  mesh->vert_array.y[indexPlusPixelsPlusOne],  mesh->vert_array.z[indexPlusPixelsPlusOne]};

            triangle_index++;
        }
    }

    vec3f geometric_center = {0, 0, 0};
    vec3f COM = {0, 0, 0};
    vec3f weightedCOMaccumulator = {0, 0, 0};
    double volume = 0;
    double mass;

    for(int i = 0; i < VERTICES; i++)
    {
        geometric_center.x += mesh->vert_array->x[i];
        geometric_center.y += mesh->vert_array->y[i];
        geometric_center.z += mesh->vert_array->z[i];
    }

    geometric_center.x /= VERTICES;
    geometric_center.y /= VERTICES;
    geometric_center.z /= VERTICES;

    for(int i = 0; i < mesh->triangle_count; i++) 
    {
        vec3f v2 = {mesh->triangles.p2.x - mesh->triangles.p1.x, mesh->triangles.p2.y - mesh->triangles.p1.y, mesh->triangles.p2.z - mesh->triangles.p1.z};
        vec3f v3 = {mesh->triangles.p3.x - mesh->triangles.p1.x, mesh->triangles.p3.y - mesh->triangles.p1.y, mesh->triangles.p3.z - mesh->triangles.p1.z};
        vec3f v4 = {geometric_center.x - mesh->triangles.p1.x, geometric_center.y - mesh->triangles.p1.y, geometric_center.z - mesh->triangles.p1.z};

        double curr_volume = fabs(dotProduct3f(&(crossProduct(&v2, &v3)), &v4)) / 6;

        vec3f tetrahedronCOM = {mesh->triangles.p1.x + mesh->mesh->triangles.p2.x + mesh->triangles.p3.x + geometric_center.x,
                                mesh->triangles.p1.y + mesh->mesh->triangles.p2.y + mesh->triangles.p3.y + geometric_center.y,
                                mesh->triangles.p1.z + mesh->mesh->triangles.p2.z + mesh->triangles.p3.z + geometric_center.z}; 

        weightedCOMaccumulator.x += curr_volume * tetrahedronCOM.x;
        weightedCOMaccumulator.y += curr_volume * tetrahedronCOM.y;
        weightedCOMaccumulator.z += curr_volume * tetrahedronCOM.z;


        double curr_mass = curr_volume * uniform_density;
        volume += curr_volume;
        mass += curr_mass;
    }

    weightedCOMaccumulator.x /= volume;
    weightedCOMaccumulator.y /= volume;
    weightedCOMaccumulator.z /= volume;

    mesh->COM = weightedCOMaccumulator;

    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh->vert_array, mesh->head, 0, VERTICES);
    
    mesh->velocity = zerovector3i;

    return;
} 

void sphere_init(Mesh* mesh, int radius, int offsetX, int offsetY, int offsetZ) 
{
    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;
    
    int index, nextL, nextQ, indexPlusOne, indexPlusPixels, indexPlusPixelsPlusOne, triangle_index;
    vec3f normalVector, partialDerivativeU, partialDerivativeV;

    mesh->triangle_count = 2 * VERTICES;
    
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
            // set a different normal vector for an adjacent triangle

            mesh->triangles[triangle_index].normal = unit3f(&normalVector);
            triangle_index++;
            mesh->triangles[triangle_index].normal = unit3f(&normalVector);
            triangle_index++;
        }
    }

    triangle_index = 0;
 
    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;
            indexPlusOne = j * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + k;
            indexPlusPixelsPlusOne = nextQ * PIXELS + k + 1;

            mesh->triangles[triangle_index].p1 = { mesh->vert_array.x[index],  mesh->vert_array.y[index],  mesh->vert_array.z[index]};
            mesh->triangles[triangle_index].p2 = { mesh->vert_array.x[indexPlusOne],  mesh->vert_array.y[indexPlusOne],  mesh->vert_array.z[indexPlusOne]};
            mesh->triangles[triangle_index].p3 = { mesh->vert_array.x[indexPlusPixels],  mesh->vert_array.y[indexPlusPixels],  mesh->vert_array.z[indexPlusPixels]};

            triangle_index++;

            mesh->triangles[triangle_index].p1 = { mesh->vert_array.x[indexPlusOne],  mesh->vert_array.y[indexPlusOne],  mesh->vert_array.z[indexPlusOne]};
            mesh->triangles[triangle_index].p2 = { mesh->vert_array.x[indexPlusPixels],  mesh->vert_array.y[indexPlusPixels],  mesh->vert_array.z[indexPlusPixels]};
            mesh->triangles[triangle_index].p3 = { mesh->vert_array.x[indexPlusPixelsPlusOne],  mesh->vert_array.y[indexPlusPixelsPlusOne],  mesh->vert_array.z[indexPlusPixelsPlusOne]};

            triangle_index++;
        }
    }

    vec3f geometric_center = {0, 0, 0};
    vec3f COM = {0, 0, 0};
    vec3f weightedCOMaccumulator = {0, 0, 0};
    double volume = 0;
    double mass;

    for(int i = 0; i < VERTICES; i++)
    {
        geometric_center.x += mesh->vert_array->x[i];
        geometric_center.y += mesh->vert_array->y[i];
        geometric_center.z += mesh->vert_array->z[i];
    }

    geometric_center.x /= VERTICES;
    geometric_center.y /= VERTICES;
    geometric_center.z /= VERTICES;

    for(int i = 0; i < mesh->triangle_count; i++) 
    {
        vec3f v2 = {mesh->triangles.p2.x - mesh->triangles.p1.x, mesh->triangles.p2.y - mesh->triangles.p1.y, mesh->triangles.p2.z - mesh->triangles.p1.z};
        vec3f v3 = {mesh->triangles.p3.x - mesh->triangles.p1.x, mesh->triangles.p3.y - mesh->triangles.p1.y, mesh->triangles.p3.z - mesh->triangles.p1.z};
        vec3f v4 = {geometric_center.x - mesh->triangles.p1.x, geometric_center.y - mesh->triangles.p1.y, geometric_center.z - mesh->triangles.p1.z};

        double curr_volume = fabs(dotProduct3f(&(crossProduct(&v2, &v3)), &v4)) / 6;

        vec3f tetrahedronCOM = {mesh->triangles.p1.x + mesh->mesh->triangles.p2.x + mesh->triangles.p3.x + geometric_center.x,
                                mesh->triangles.p1.y + mesh->mesh->triangles.p2.y + mesh->triangles.p3.y + geometric_center.y,
                                mesh->triangles.p1.z + mesh->mesh->triangles.p2.z + mesh->triangles.p3.z + geometric_center.z}; 

        weightedCOMaccumulator.x += curr_volume * tetrahedronCOM.x;
        weightedCOMaccumulator.y += curr_volume * tetrahedronCOM.y;
        weightedCOMaccumulator.z += curr_volume * tetrahedronCOM.z;


        double curr_mass = curr_volume * uniform_density;
        volume += curr_volume;
        mass += curr_mass;
    }

    mesh->V = volume;
    mesh->M = mass;

    weightedCOMaccumulator.x /= volume;
    weightedCOMaccumulator.y /= volume;
    weightedCOMaccumulator.z /= volume;

    mesh->COM = weightedCOMaccumulator;
    
    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh->vert_array, mesh->head, 0, VERTICES);

    mesh->velocity = zerovector3i;

    return;
}

void mesh_init(Mesh* mesh, triangle* outTriangles, uint32_t outTriangleCount, int offsetX, int offsetY, int offsetZ, double uniform_density) 
{
    mesh->pos.x = offsetX;
    mesh->pos.y = offsetY;
    mesh->pos.z = offsetZ;

    int triangle_index = 0;
    
    for(int i = 0; i < outTriangleCount; i++) 
    {
        mesh->triangles[i] = outTriangles[i];   
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {
            index = j * PIXELS + k; 
            nextL = (k + 1) % PIXELS;
            nextQ = (j + 1) % PIXELS;
            indexPlusOne = j * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + k;

            mesh->vert_array->x[index] = mesh->triangles[triangle_index].p1.x;
            mesh->vert_array->y[index] = mesh->triangles[triangle_index].p1.y;
            mesh->vert_array->z[index] = mesh->triangles[triangle_index].p1.z;      

            triangle_index++;
        }
    }

    // calculate center of mass

    vec3f geometric_center = {0, 0, 0};
    vec3f COM = {0, 0, 0};
    vec3f weightedCOMaccumulator = {0, 0, 0};
    double volume = 0;
    double mass;

    for(int i = 0; i < VERTICES; i++)
    {
        geometric_center.x += mesh->vert_array->x[i];
        geometric_center.y += mesh->vert_array->y[i];
        geometric_center.z += mesh->vert_array->z[i];
    }

    geometric_center.x /= VERTICES;
    geometric_center.y /= VERTICES;
    geometric_center.z /= VERTICES;

    for(int i = 0; i < outTriangles; i++) 
    {
        vec3f v2 = {mesh->triangles.p2.x - mesh->triangles.p1.x, mesh->triangles.p2.y - mesh->triangles.p1.y, mesh->triangles.p2.z - mesh->triangles.p1.z};
        vec3f v3 = {mesh->triangles.p3.x - mesh->triangles.p1.x, mesh->triangles.p3.y - mesh->triangles.p1.y, mesh->triangles.p3.z - mesh->triangles.p1.z};
        vec3f v4 = {geometric_center.x - mesh->triangles.p1.x, geometric_center.y - mesh->triangles.p1.y, geometric_center.z - mesh->triangles.p1.z};

        double curr_volume = fabs(dotProduct3f(&(crossProduct(&v2, &v3)), &v4)) / 6;

        vec3f tetrahedronCOM = {(mesh->triangles.p1.x + mesh->mesh->triangles.p2.x + mesh->triangles.p3.x + geometric_center.x) / 4,
                                (mesh->triangles.p1.y + mesh->mesh->triangles.p2.y + mesh->triangles.p3.y + geometric_center.y) / 4,
                                (mesh->triangles.p1.z + mesh->mesh->triangles.p2.z + mesh->triangles.p3.z + geometric_center.z) / 4}; 

        weightedCOMaccumulator.x += curr_volume * tetrahedronCOM.x;
        weightedCOMaccumulator.y += curr_volume * tetrahedronCOM.y;
        weightedCOMaccumulator.z += curr_volume * tetrahedronCOM.z;


        double curr_mass = curr_volume * uniform_density;
        volume += curr_volume;
        mass += curr_mass;
    }

    weightedCOMaccumulator.x /= volume;
    weightedCOMaccumulator.y /= volume;
    weightedCOMaccumulator.z /= volume;

    mesh->COM = weightedCOMaccumulator;

    // inertia tensor

    mat3f I_origin = mat3f_zero();
    double mass_total = 0.0;

    for(int t = 0; t < mesh->triangle_count; t++) 
    {
        vec3f p1 = mesh->triangles[t].p1;
        vec3f p2 = mesh->triangles[t].p2;
        vec3f p3 = mesh->triangles[t].p3;
        vec3f p4 = geometric_center;

        vec3f a = { p2.x - p1.x, p2.y - p1.y, p2.z - p1.z };
        vec3f b = { p3.x - p1.x, p3.y - p1.y, p3.z - p1.z };
        vec3f c = { p4.x - p1.x, p4.y - p1.y, p4.z - p1.z };

        double V = fabs(dotProduct3f(a, crossProduct(b, c))) / 6.0;
        double m = V * uniform_density;
        mass_total += m;

        double x[4] = { p1.x, p2.x, p3.x, p4.x };
        double y[4] = { p1.y, p2.y, p3.y, p4.y };
        double z[4] = { p1.z, p2.z, p3.z, p4.z };

        double Ixx =  m * 0.1 * (S2(y) + S2(z));
        double Iyy =  m * 0.1 * (S2(z) + S2(x));
        double Izz =  m * 0.1 * (S2(x) + S2(y));
        double Ixy = -m * 0.05 * SP(x, y);
        double Iyz = -m * 0.05 * SP(y, z);
        double Izx = -m * 0.05 * SP(z, x);

        mat3f Io = {{{ Ixx, Ixy, Izx },
                    { Ixy, Iyy, Iyz },
                    { Izx, Iyz, Izz }}};

        I_origin = mat3f_add(I_origin, Io);
    }

    double COM2 = dotProduct3f(COM, COM);
    mat3f shift = mat3f_add(mat3f_identity(mass_total * COM2), mat3f_scale(mat3f_outer(COM, COM), -mass_total));
    mat3f I_body = mat3f_add(I_origin, mat3f_scale(shift, -1));

    mesh->I_body = I_body;

    memcpy(mesh->head->boundingBox, generateBoundingBox(mesh));
    initializeBVHTree(mesh->head, 0);
    createBoundingBoxForNode(mesh->vert_array, mesh->head, 0, VERTICES);

    mesh->velocity = zerovector3i;

    return;
}

#endif
