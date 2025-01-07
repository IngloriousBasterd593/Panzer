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

int get_space(Mesh* mesh) 
{
    mesh->x = malloc(VERTICES * sizeof(float));
    if(mesh->x == NULL) 
    {
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->y = malloc(VERTICES * sizeof(float));
    if(mesh->y == NULL) 
    {
        free(mesh->x);
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->z = malloc(VERTICES * sizeof(float));
    if(mesh->z == NULL) 
    {
        free(mesh->x);
        free(mesh->y);
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->xProj = malloc(VERTICES * sizeof(float));
    if(mesh->xProj == NULL) 
    {
        free(mesh->x);
        free(mesh->y);
        free(mesh->z);
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->yProj = malloc(VERTICES * sizeof(float));
    if(mesh->yProj == NULL) 
    {
        free(mesh->x);
        free(mesh->y);
        free(mesh->z);
        free(mesh->xProj);
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->zProj = malloc(VERTICES * sizeof(float));
    if(mesh->zProj == NULL) 
    {
        free(mesh->x);
        free(mesh->y);
        free(mesh->z);
        free(mesh->xProj);
        free(mesh->yProj);
        fprintf(stderr, "couldn't init mesh\n");
        return 1;
    }

    mesh->boundingBoxes = (BoundingBox*) malloc(sizeof(BoundingBox) * BOUNDINGBOXCOUNT);
    if(mesh->boundingBoxes == NULL)
    {
        fprintf(stderr, "Failed to initialize heap\n");
        return 1;
    }

    return 0;
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