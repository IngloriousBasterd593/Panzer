#ifndef SDL_GRAPHICS_LIB
#define SDL_GRAPHICS_LIB

// external libs

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h> 
#include <math.h>
#include <unistd.h>
#include <time.h>


// SDL && OpenGL

#include "C:\glfw-3.4.bin.WIN64\include\GLFW\glfw3.h"   

#define UNITY_BUILD 1             
#ifdef _WIN64
 #include <SDL.h>
 #include <SDL_image.h>
#else
 #include <SDL2/SDL.h>
 #include <SDL2/SDL_image.h>
#endif
#include "common.h"
#include "sdl_utils.h"


// macros

#define PI 3.141592653589793
#define TWOPI 2 * PI
#define TWOPIOVERPIXELS TWOPI / (PIXELS - 1)
#define PIOVERPIXELS PI / (PIXELS - 1)
#define PIXELS 30       // 25 pixels - 35 precision     
#define POINTS PIXELS * PIXELS
#define launadainais 1
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


int get_space(Manifold* manifold);


Vector3 crossproduct(Vector3* vector1, Vector3* vector2);


Vector3 unit(Vector3* v);


float dotproduct(Vector3* vector1, Vector3* vector2);


void sphere_init(Manifold* manifold, Vector3 manifoldNormals[], float radius);


void torus_init(Manifold* manifold, Vector3 manifoldNormals[], float Rinner, float Router);


void line(SDL_Renderer* renderer, Vector2 vector1, Vector2 vector2, int R, int G, int B, int Xoffset, int Yoffset);


void fillTriangle(SDL_Renderer* renderer, Vector2 baseVertex1, Vector2 baseVertex2, Vector2 centralVertex, float trianglePrecision, int R, int G, int B, int Xoffset, int Yoffset);


void Manifold_draw(SDL_Renderer* renderer, Manifold* manifold, Vector3 manifoldNormals[], int Xoffset, int Yoffset, int trianglePrecision);


void Manifold_rotate(Manifold* manifold, Vector3 manifoldNormals[], float rad, char axis);


#endif 