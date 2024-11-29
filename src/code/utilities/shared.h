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
        fprintf(stderr, "couldnt init sdl\n");
        return 1;
    }

    *window = SDL_CreateWindow(name, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
    if(*window == NULL) 
    {
        fprintf(stderr, "couldnt init window\n");
        SDL_Quit();
        return 1;
    }

    *renderer = SDL_CreateRenderer(*window, -1, SDL_RENDERER_ACCELERATED);
    if(*renderer == NULL) 
    {
        fprintf(stderr, "couldnt init renderer\n");
        SDL_DestroyWindow(*window);
        SDL_Quit();
        return 1;
    }

    *texture = SDL_CreateTexture(*renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, width, height);
    if(*texture == NULL) 
    {
        fprintf(stderr, "couldnt init texture\n");
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



int get_space(Manifold* manifold) 
{

    manifold->x = malloc(VERTICES * sizeof(float));
    if(manifold->x == NULL) 
    {
        fprintf(stderr, "couldnt init manifold\n");
        return 1;
    }

    manifold->y = malloc(VERTICES * sizeof(float));
    if(manifold->y == NULL) 
    {
        free(manifold->x);
        fprintf(stderr, "couldnt init manifold\n");
        return 1;
    }

    manifold->z = malloc(VERTICES * sizeof(float));
    if(manifold->z == NULL) 
    {
        free(manifold->x);
        free(manifold->y);
        fprintf(stderr, "couldnt init manifold\n");
        return 1;
    }

    manifold->xProj = malloc(VERTICES * sizeof(float));
    if(manifold->xProj == NULL) 
    {
        free(manifold->x);
        free(manifold->y);
        free(manifold->z);
        fprintf(stderr, "couldnt init manifold\n");
        return 1;
    }

    manifold->yProj = malloc(VERTICES * sizeof(float));
    if(manifold->yProj == NULL) 
    {
        free(manifold->x);
        free(manifold->y);
        free(manifold->z);
        free(manifold->xProj);
        fprintf(stderr, "couldnt init manifold\n");
        return 1;
    }

    manifold->zProj = malloc(VERTICES * sizeof(float));
    if(manifold->zProj == NULL) 
    {
        free(manifold->x);
        free(manifold->y);
        free(manifold->z);
        free(manifold->xProj);
        free(manifold->yProj);
        fprintf(stderr, "couldnt init manifold\n");
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