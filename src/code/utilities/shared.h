#ifndef SDLUTILITIES
#define SDLUTILITIES

#include "common.h"

// include SDL

#ifndef UNITY_BUILD
#define UNITY_BUILD 1
#endif       
#ifdef _WIN64
 #include <SDL.h>
 #include <SDL_image.h>
#else
 #include <SDL2/SDL.h>
 #include <SDL2/SDL_image.h>
#endif



int readSTL(const char* filename, triangle** outTriangles, uint32_t* outTriangleCount)
{
    FILE* file = fopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        return false;
    }

    char header[80];
    if(fread(header, 1, 80, file) != 80) 
    {
        fprintf(stderr, "Failed to read STL header\n");
        fclose(file);
        return false;
    }

    // Read the number of triangles (4 bytes)
    uint32_t triangleCount = 0;
    if (fread(&triangleCount, sizeof(uint32_t), 1, file) != 1) 
    {
        fprintf(stderr, "Failed to read triangle count\n");
        fclose(file);
        return false;
    }

    triangle* triangles = (triangle*) malloc(triangleCount * sizeof(triangle));
    if(!triangles) 
    {
        fprintf(stderr, "Failed to allocate memory for triangles\n");
        fclose(file);
        return false;
    }

    for(uint32_t i = 0; i < triangleCount; i++) 
    {
        if(fread(&triangles[i].normal, sizeof(vec3f), 1, file) != 1) 
        {
            fprintf(stderr, "Failed to read triangle normal\n");
            free(triangles);
            fclose(file);
            return false;
        }

        if(fread(&triangles[i].p1, sizeof(vec3f), 1, file) != 1 ||
           fread(&triangles[i].p2, sizeof(vec3f), 1, file) != 1 ||
           fread(&triangles[i].p3, sizeof(vec3f), 1, file) != 1) 
        {
            fprintf(stderr, "Failed to read triangle vertices\n");
            free(triangles);
            fclose(file);
            return false;
        }

        fseek(file, 2, SEEK_CUR);
    }

    fclose(file);

    *outTriangles = triangles;
    *outTriangleCount = triangleCount;

    return true;
}




int SDL_init(SDL_Window** window, SDL_Renderer** renderer, SDL_Texture** texture, const char* name, int width, int height) 
{
    if(SDL_Init(SDL_INIT_VIDEO) < 0) 
    {
        fprintf(stderr, "couldn't init SDL\n");
        return 1;
    }

    *window = SDL_CreateWindow(name, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
    if(*window == NULL) 
    {
        fprintf(stderr, "couldn't init window\n");
        SDL_Quit();
        return 1;
    }

    *renderer = SDL_CreateRenderer(*window, -1, SDL_RENDERER_ACCELERATED);
    if(*renderer == NULL) 
    {
        fprintf(stderr, "couldn't init renderer\n");
        SDL_DestroyWindow(*window);
        SDL_Quit();
        return 1;
    }

    *texture = SDL_CreateTexture(*renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, width, height);
    if(*texture == NULL) 
    {
        fprintf(stderr, "couldn't init texture\n");
        SDL_DestroyWindow(*window);
        SDL_DestroyRenderer(*renderer);
        SDL_Quit();
        return 1;
    }

    SDL_SetRenderDrawColor(*renderer, 255, 255, 255, 255);
    SDL_RenderClear(*renderer);
    SDL_RenderPresent(*renderer);
    SDL_SetRenderTarget(*renderer, *texture);
    SDL_RenderPresent(*renderer);

    return 0;
}

int get_space(Mesh* meshes) 
{
    mesh->triangle = malloc(VERTICES * sizeof(triangle));
    if(mesh->triangle == NULL) 
    {
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->AABBes = (AABB*) malloc(sizeof(AABB) * AABBCOUNT);
    if(mesh->AABBes == NULL)
    {
        fprintf(stderr, "Failed to initialize heap\n");
        return 1;
    }

    return 0;
}

void free_space(Scene* sceneInstance) 
{
    for(int i = 0; i < sceneInstance->meshCount; i++) 
    {
        free(sceneInstance->meshes[i]->triangles);
        free(sceneInstance->meshes[i]->AABBes);
        free(sceneInstance->meshes[i]);
    }

    free(sceneInstance->meshes);
    free(sceneInstance->shaderProgram);
    free(sceneInstance->frameColors);
    free(sceneInstance);

    return;
}

// scene constructor
void initializeScene(Scene* sceneInstance, int meshCount, int drawPrecision, AABB simulationSpace)
{
    sceneInstance = malloc(sizeof(Scene));
    if(sceneInstance == NULL) 
    {
        fprintf(stderr, "couldn't get space for scene");
        return;
    }

    sceneInstance->meshCount = meshCount;
    sceneInstance->drawPrecision = drawPrecision;

    sceneInstance->Mesh** meshes = malloc(sceneInstance->meshCount * sizeof(Mesh*));
    if(sceneInstance->meshes == NULL) 
    {
        fprintf(stderr, "couldn't get space for meshes");
    }

    for(int i = 0; i < sceneInstance->meshCount; i++) 
    {
        sceneInstance->meshes[i] = malloc(sizeof(Mesh));
        if(sceneInstance->meshes[i] == NULL) 
        {
            fprintf(stderr, "couldn't get space for mesh");
            
        }

        if(get_space(sceneInstance->meshes[i]) == 1) 
        {
            
        }
    }

    sceneInstance->shaderProgram = malloc(BUFFLEN * sizeof(char));
    if(sceneInstance->shaderProgram == NULL) 
    {
        fprintf(stderr, "couldn't get space for shader program");
        return 1;
    }

    vec3f lightPerspectiveVector = {0, 0, 1};
    vec3f viewVector = {0, 0, 1};

    sceneInstance->camera = { viewVector, PI / 3, viewVector, SCREENWIDTH / SCREENHEIGHT, 1, 100, 10, 20, 10, 20 };

    sceneInstance->frameColors = malloc(SCREENSIZE * sizeof(unsigned int));
    if(sceneInstance->frameColors == NULL) 
    {
        perror("bruh");
    } 

    memcpy(sceneInstance->simulationSpace, &simulationSpace, sizeof(AABB));
    partitionSceneSimulationSpace(sceneInstance);

    return;
}

void printVector3f(vec3f vector) 
{
    printf("%f %f %f\n", vector.x, vector.y, vector.z);

    return;
}



void printVector4f(vec4f vector) 
{
    printf("%f %f %f %f\n", vector.x, vector.y, vector.z, vector.w);

    return;
}



void printMatrix4f(mat4f mat) 
{
    for(int i = 0; i < 4; i++) 
    {
        for(int j = 0; j < 4; j++) {
            printf("[ %f ]", mat.mat[j][i]);
        }
        printf("\n");
    }

    printf("\n");

    return;
}



#endif