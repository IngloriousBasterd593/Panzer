#ifndef COMMONH
#define COMMONH

/*
*
*       This is supposed to be a 3d physics simulation.
*       The implementation is inherently bloated and cant be fixed, at this point it's easier to rewrite the entire thing than to fix the existing codebase
*       Do not use any of these functions for any purpose, as they lack proper error handling (designed to make them faster) and are non refactorable, 
*       lacking proper documentation and inhereting implementation issues from the person who wrote them.
*       
*/



// external libs

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>


// macros

#define PIXELS 160                                  // loop counter for 3 dimensional shape access    
#define VERTICES (PIXELS * PIXELS)                  // number of vertices for each object
#define PI M_PI
#define TWOPI (2 * PI)
#define TWOPIOVERPIXELS ((TWOPI) / (PIXELS - 1))    // step size for initialization
#define PIOVERPIXELS ((PI) / (PIXELS - 1))          // step size for initialization
#define launadainais 1
#define SCREENWIDTH 1280
#define SCREENHEIGHT 720
#define SCREENSIZE (SCREENWIDTH * SCREENHEIGHT)
#define HALFWINWIDTH (SCREENWIDTH / 2)
#define HALFWINHEIGHT (SCREENHEIGHT / 2)
#define MAXLOOPLEN 4096                             // lmao
#define BUFFLEN 1024                                // safety
#define FAR 2                                     // far plane distance for perspective projection
#define NEAR 1                                    // near plane distance for perspective projection

// structs

typedef struct {
    float x;
    float y;
    float z;
} vec3f;

typedef struct {
    float x;
    float y;
    float z;
    float w;
} vec4f;

typedef union {
    float mat[4][4];
    float raw[16];
    vec4f column[4];
} mat4f;

typedef struct {
    float x;
    float y;
} vec2f;

typedef struct {
    int x;
    int y;
    int z;
} vec3i;

typedef struct {
    int x;
    int y;
} vec2i;

typedef struct {
    float* x;
    float* y;
    float* z;
    float* xProj;
    float* yProj;
    float* zProj;
    int Xposition;
    int Yposition;
    int Zposition;
} Manifold;

typedef struct {
    vec3f POV;
    float FOV;
    vec3f lightingDirectionVector;
    float aspectRatio;
    float nearPlane;
    float farPlane;
    float far;
    float near;
} Camera;


#endif
