#ifndef SDL_GRAPHICS_LIB
#define SDL_GRAPHICS_LIB


// external libs

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h> 
#include <math.h>
#include <unistd.h>
#include <time.h>


// SDL

#define UNITY_BUILD 1             
#ifdef _WIN64
 #include <SDL.h>
 #include <SDL_image.h>
#else
 #include <SDL2/SDL.h>
 #include <SDL2/SDL_image.h>
#endif
#include "C:\sdl-vscode-c-cpp\common.h"
#include "C:\sdl-vscode-c-cpp\sdl_utils.h"


// SDL addition




// macros

#define launadainais 1
#define PI 3.141592653589793
#define PIXELS 50      // 50 for 30 fps
#define POINTS PIXELS * PIXELS
#define S_WIDTH 1280
#define S_HEIGHT 720
#define C_winX S_WIDTH / 2
#define C_winY S_HEIGHT / 2
#define MAX_LOOP_LEN 4096


// structs

typedef struct{
    float x;
    float y;
    float z;
} Vector3;


typedef struct{
    float x;
    float y;
} Vector2;


typedef struct{
    float* x;
    float* y;
    float* z;
} Manifold;



// function prototypes 


SDL_Renderer* SDL_INIT(SDL_Window** window, const char* name, int width, int height);


void get_space(Manifold* manifold);


Vector3 normal(Vector3 vector1, Vector3 vector2);


Vector3 unit(Vector3 vector);


float dotproduct(Vector3 vector1, Vector3 vector2);


void sphere_init(Manifold* manifold, float radius);


void torus_init(Manifold* manifold, float Rinner, float Router);


void line(SDL_Renderer* renderer, Vector2 vector1, Vector2 vector2, int R, int G, int B, int Xoffset, int Yoffset);


void fillTriangle(SDL_Renderer* renderer, Vector2 vertex1, Vector2 vertex2, Vector2 vertexC, float precision, int R, int G, int B, int Xoffset, int Yoffset);


void sphere_draw(SDL_Renderer* renderer, Manifold manifold, int Xoffset, int Yoffset);


void torus_draw(SDL_Renderer* renderer, Manifold* manifold, int Xoffset, int Yoffset, int precision);


void M_rotate(Manifold* manifold, float rad, char axis);


#endif 