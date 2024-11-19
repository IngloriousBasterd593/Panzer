#ifndef SDL_GRAPHICS_LIB
#define SDL_GRAPHICS_LIB



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


// OpenGL

#include <glad/glad.h>  
#include <GLFW/glfw3.h>


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

typedef struct {
    float mat[4][4];
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
} Camera;


// Function Implementations

void copyShaderSource(const char* source, char* destination) 
{
    if(source == NULL || destination == NULL) 
    {
        fprintf(stderr, "source or destination dont exist\n");
    }

    FILE* shaderFile = fopen(source, "r");
    if(shaderFile == NULL) 
    {
        fprintf(stderr, "failed to open file\n");
        return;
    }

    int i = 0; 
    char ch;

    while((ch = fgetc(shaderFile)) != EOF) 
    {
        if(i > BUFFLEN - 2) 
        {
            destination[i] = '\0';
            break;
        }

        destination[i] = ch;
        i++;
    }

    fclose(shaderFile);

    return;
}

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



vec3f unit(vec3f* v) 
{

    float magnitude = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
    vec3f normalized;

    if(magnitude == 0.0f) 
    {
        // fprintf(stderr, "magnitude zero\n");
        normalized.x = normalized.y = normalized.z = 0;
        return normalized;
    } 

    normalized.x = v->x / magnitude;
    normalized.y = v->y / magnitude;
    normalized.z = v->z / magnitude;

    return normalized;
}



vec3f crossproduct(vec3f* v1, vec3f* v2) 
{
    return (vec3f) {
                    v1->y * v2->z - v1->z * v2->y, 
                    v1->z * v2->x - v1->x * v2->z, 
                    v1->x * v2->y - v1->y * v2->x 
    };
}



float dotproduct(vec3f* v1, vec3f* v2) 
{
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}



void sphere_init(Manifold* manifold, vec3f* manifoldNormals, int radius, int offsetX, int offsetY, int offsetZ) 
{

    if(radius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    manifold->Xposition = HALFWINWIDTH;
    manifold->Yposition = HALFWINHEIGHT;
    manifold->Zposition = 2000;
    
    int index;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    vec3f normalVector;
    vec3f partialDerivativeU;
    vec3f partialDerivativeV;
    
    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        index = j * PIXELS + k; 
        
        manifold->x[index] = ((radius * cos(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j)));
        manifold->y[index] = ((radius * sin(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j)));
        manifold->z[index] = ((radius * cos(PIOVERPIXELS * j)));

        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        nextL = (k + 1) % PIXELS;
        nextQ = (j + 1) % PIXELS;

        index = j * PIXELS + k; 
        indexPlusOne = j * PIXELS + nextL;
        indexPlusPixels = nextQ * PIXELS + k;

        partialDerivativeU.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusOne]);
        partialDerivativeU.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusOne]);
        partialDerivativeU.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusOne]);

        partialDerivativeV.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusPixels]);
        partialDerivativeV.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusPixels]);
        partialDerivativeV.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusPixels]);
        
        normalVector = crossproduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit(&normalVector);

        }
    }

    return;
}  



void torus_init(Manifold* manifold, vec3f* manifoldNormals, int innerRadius, int outerRadius, int offsetX, int offsetY, int offsetZ) 
{

    if(innerRadius <= 10 || outerRadius <= 10) 
    {
        fprintf(stderr, "write a correct radius, bozo\n");
        return;
    }

    manifold->Xposition = HALFWINWIDTH;
    manifold->Yposition = HALFWINHEIGHT;
    manifold->Zposition = 700;

    vec3f normalVector;
    vec3f partialDerivativeU;
    vec3f partialDerivativeV;

    int index;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        index = j * PIXELS + k; 

        manifold->x[index] = (((outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j)) + offsetX);
        manifold->y[index] = (((outerRadius + innerRadius * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j)) + offsetY);
        manifold->z[index] = ((innerRadius * sin(TWOPIOVERPIXELS * k)) + offsetZ);

        }
    }

    for(int j = 0; j < PIXELS; j++) 
    {
        for(int k = 0; k < PIXELS; k++) 
        {

        nextL = (k + 1) % PIXELS;
        nextQ = (j + 1) % PIXELS;

        index = j * PIXELS + k; 
        indexPlusOne = j * PIXELS + nextL;
        indexPlusPixels = nextQ * PIXELS + k;

        partialDerivativeU.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusOne]);
        partialDerivativeU.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusOne]);
        partialDerivativeU.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusOne]);

        partialDerivativeV.x = (int)(10000 * manifold->x[index]) - (int)(10000 * manifold->x[indexPlusPixels]);
        partialDerivativeV.y = (int)(10000 * manifold->y[index]) - (int)(10000 * manifold->y[indexPlusPixels]);
        partialDerivativeV.z = (int)(10000 * manifold->z[index]) - (int)(10000 * manifold->z[indexPlusPixels]);
        
        normalVector = crossproduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit(&normalVector);

        }
    }

    return;
} 



void lineBresenham(Manifold* manifold, unsigned int* frameColors, vec2f vertexStart, vec2f vertexEnd, unsigned int color) 
{
    
    short x1 = vertexStart.x;
    short y1 = vertexStart.y;
    short x2 = vertexEnd.x;
    short y2 = vertexEnd.y;

    short dx = abs(x2 - x1);
    short dy = abs(y2 - y1);

    short stepX = (x1 < x2) ? 1 : -1;
    short stepY = (y1 < y2) ? 1 : -1;

    short error = dx - dy;
    short error2;

    while(1) 
    {
        // not out of bounds, trust me
        if(((y1 + manifold->Yposition) * SCREENWIDTH) + x1 + manifold->Xposition > SCREENSIZE || ((y1 + manifold->Yposition) * SCREENWIDTH) + x1 + manifold->Xposition < 0) 
        {
            break;
        }

        frameColors[((y1 + manifold->Yposition) * SCREENWIDTH) + x1 + manifold->Xposition] = color;

        if(x1 == x2 && y1 == y2) 
        {
            break;
        }

        error2 = 2 * error;

        if(error2 > -dy) 
        {
            error -= dy;
            x1 += stepX;
        }

        if(error2 < dx) 
        {
            error += dx;
            y1 += stepY;
        }
    }

    return;
}



void fillTriangle(Manifold* manifold, unsigned int* frameColors, vec2f baseVertex1, vec2f baseVertex2, vec2f centralVertex, int trianglePrecision, unsigned int color) 
{

    float dx1 = (centralVertex.x - baseVertex1.x) / trianglePrecision;
    float dy1 = (centralVertex.y - baseVertex1.y) / trianglePrecision;
    float dx2 = (centralVertex.x - baseVertex2.x) / trianglePrecision;
    float dy2 = (centralVertex.y - baseVertex2.y) / trianglePrecision;

    vec2f point1;
    vec2f point2;

    for (int i = 0; i < trianglePrecision; i++) 
    {

        point1.x = baseVertex1.x + dx1;
        point1.y = baseVertex1.y + dy1;

        point2.x = baseVertex2.x + dx2;
        point2.y = baseVertex2.y + dy2;

        lineBresenham(manifold, frameColors, point1, point2, color);

    }

    return;
}



void fillRectangle(Manifold* manifold, unsigned int* frameColors, vec2f vertexA, vec2f vertexB, vec2f vertexC, vec2f vertexD, int drawPrecision, unsigned int color) 
{
    /*
    if(drawPrecision <= 0) 
    {
        fprintf(stderr, "write a correct drawprecision, bozo\n");
        return;
    } */

    float dxU = (vertexD.x - vertexC.x) / drawPrecision;
    float dyU = (vertexD.y - vertexC.y) / drawPrecision;
    float dxL = (vertexB.x - vertexA.x) / drawPrecision;
    float dyL = (vertexB.y - vertexA.y) / drawPrecision;

    for(int i = 0; i <= drawPrecision; i++) 
    {

        vertexA.x += dxU;
        vertexA.y += dyU;

        vertexC.x += dxL;
        vertexC.y += dyL;

        lineBresenham(manifold, frameColors, vertexA, vertexC, color);
    }

    return;
}




mat4f* OrthographicProjectionMatrix4f(float left, float right, float bottom, float top, float nearZ, float farZ, Matrix4f* result) // applies to the entire matrix
{
    result->vecRows[0] = (Vector4f){ 2.0/(right-left), 0.0, 0.0, -(right+left)/(right-left) };
    result->vecRows[1] = (Vector4f){ 0.0, 2.0/(top-bottom), 0.0, -(top+bottom)/(top-bottom) };
    result->vecRows[2] = (Vector4f){ 0.0, 0.0, -2.0/(farZ-nearZ), -(farZ+nearZ)/(farZ-nearZ) };
    result->vecRows[3] = (Vector4f){ 0.0, 0.0, 0.0, 1.0 };
    return result;
}


/*
mat4f* OrthographicProjectionMatrix4f(float left, float right, float bottom, float top, float nearZ, float farZ, mat4f* result) 
{
    result->vecRows[0] = { 2.0/(right-left), 0.0, 0.0, -(right+left)/(right-left) };
    result->vecRows[1] = { 0.0, 2.0/(top-bottom), 0.0, -(top+bottom)/(top-bottom) };
    result->vecRows[2] = { 0.0, 0.0, -2.0/(farZ-nearZ), -(farZ+nearZ)/(farZ-nearZ) };
    result->vecRows[3] = { 0.0, 0.0, 0.0, 1.0 };
    return result;
}

*/

/*
vec3f perspectiveProject(vec3f* vertex, Camera* camera) 
{
    float f = 1.0f / tanf(camera->FOV / 2.0f);
    
    vec4f projectedVector = {
                            (((vertex->x * f) / (camera->aspectRatio * vertex->z)) + 1.0f) * (SCREENWIDTH / 2.0f),
                            (1.0f - ((vertex->y * f) / vertex->z)) * (SCREENHEIGHT / 2.0f),
                            ((camera->farPlane + camera->nearPlane) * vertex->z + 2.0f * camera->farPlane * camera->nearPlane) / (camera->nearPlane - camera->farPlane), 
                            vertex->z
    };

    return (vec3f) {projectedVector.x, projectedVector.y, projectedVector.z};
}
*/


void Manifold_draw(Manifold* manifold, vec3f* manifoldNormals, Camera* camera, unsigned int* frameColors, int drawPrecision) 
{

    if(drawPrecision < 1) 
    {
        fprintf(stderr, "write a correct drawprecision, bozo\n");
        return;
    }

    int index;
    int indexPlusPixelsPlusOne;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    int grayscaleRGB;
    float grayscaleCoefficient;
    unsigned int color;
    float f = 1.0f / tanf(camera->FOV / 2.0f);

    vec2f vertex1;
    vec2f vertex2;
    vec2f vertexUpper1;
    vec2f vertexUpper2;

    // perform perspective projection
    for(int q = 0; q < PIXELS; q++) 
    {
        for(int l = 0; l < PIXELS; l++) 
        {
            index = q * PIXELS + l;

            /*
            if (manifold->z[index] <= 0.1f) 
            {
                // perror("lmao");
                continue; 
            }
            */

            vec3f projectedVertex = {((manifold->x[index] + manifold->Xposition) / manifold->z[index]) + HALFWINWIDTH, ((manifold->y[index] + manifold->Yposition) / (manifold->z[index] + manifold->Xposition)) + HALFWINHEIGHT, manifold->z[index]};
            // printf("%f %f %f\n", projectedVertex.x, projectedVertex.y, projectedVertex.z);



/*
            vec4f projectedVertex = {
                ((((manifold->x[index] + manifold->Xposition) * f) / (camera->aspectRatio * (manifold->z[index] + manifold->Zposition))) + 1.0f) * (SCREENWIDTH / 2.0f),
                ((1.0f - (((manifold->y[index] + manifold->Yposition) * f) / manifold->z[index])) * (SCREENHEIGHT / 2.0f)),
                ((camera->farPlane + camera->nearPlane) * (manifold->z[index] + manifold->Zposition) + 2.0f * camera->farPlane * camera->nearPlane) / (camera->nearPlane - camera->farPlane),
                1.0f
            };
*/ 
            manifold->xProj[index] = projectedVertex.x;
            manifold->yProj[index] = projectedVertex.y;
            manifold->zProj[index] = projectedVertex.z; 

/*
            manifold->xProj[index] = manifold->x[index];
            manifold->yProj[index] = manifold->y[index];
            manifold->zProj[index] = manifold->z[index]; */

            //printf("%f %f %f", manifold->x[index] + manifold->Xposition, manifold->y[index] + manifold->Yposition, manifold->z[index] + manifold->Zposition);
            //printf("\t%f %f %f\n", manifold->xProj[index], manifold->yProj[index], manifold->zProj[index]);
            //usleep(1000);

        }
    }


    for(int q = 0; q < PIXELS; q++) 
    {
        for(int l = 0; l < PIXELS; l++) 
        {
            index = q * PIXELS + l;

            if(dotproduct(&camera->POV, &manifoldNormals[index]) < 0) {
                continue;
            }
            
            nextL = (l + 1) % PIXELS;
            nextQ = (q + 1) % PIXELS;

            indexPlusOne = q * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + l;
            indexPlusPixelsPlusOne = nextQ * PIXELS + nextL;

            grayscaleCoefficient = (dotproduct(&manifoldNormals[index], &camera->lightingDirectionVector));

            grayscaleRGB = (uint8_t) (255 - (132 * (1 - grayscaleCoefficient)));

            color = (0xFF << 24) | (grayscaleRGB << 16) | (grayscaleRGB << 8) | grayscaleRGB;

            vertex1.x = manifold->xProj[index];
            vertex1.y = manifold->yProj[index];

            vertex2.x = manifold->xProj[indexPlusOne];
            vertex2.y = manifold->yProj[indexPlusOne];

            vertexUpper1.x = manifold->xProj[indexPlusPixels];
            vertexUpper1.y = manifold->yProj[indexPlusPixels];

            vertexUpper2.x = manifold->xProj[indexPlusPixelsPlusOne];
            vertexUpper2.y = manifold->yProj[indexPlusPixelsPlusOne];   
             
/*
            vertex1.x = manifold->x[index];
            vertex1.y = manifold->y[index];

            vertex2.x = manifold->x[indexPlusOne];
            vertex2.y = manifold->y[indexPlusOne];

            vertexUpper1.x = manifold->x[indexPlusPixels];
            vertexUpper1.y = manifold->y[indexPlusPixels];

            vertexUpper2.x = manifold->x[indexPlusPixelsPlusOne];
            vertexUpper2.y = manifold->y[indexPlusPixelsPlusOne]; 
*/

            fillRectangle(manifold, frameColors, vertex1, vertex2, vertexUpper1, vertexUpper2, drawPrecision, color);

            // fillTriangle(frameColors, vertex1, vertex2, vertexUpper1, drawPrecision, color);
            // fillTriangle(frameColors, vertex2, vertexUpper1, vertexUpper2, drawPrecision, color);

        }
    }

    return;
}  



void Manifold_rotate(Manifold* manifold, vec3f* manifoldNormals, float rad, char axis) 
{

    vec3f previousVector;

    float cosRad = (float) cos(rad);
    float sinRad = (float) sin(rad);

    int index;
    
    switch(axis) 
    {
        case 'x':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {
                    // Rotate manifold VERTICES
                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->y[index] = previousVector.y * cosRad - previousVector.z * sinRad;
                    manifold->z[index] = previousVector.y * sinRad + previousVector.z * cosRad;

                    // Rotate normal vectors

                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].y = previousVector.y * cosRad - previousVector.z * sinRad;
                    manifoldNormals[index].z = previousVector.y * sinRad + previousVector.z * cosRad;

                }
            }
        break;

        case 'y':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {
                    
                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->x[index] = previousVector.x * cosRad + previousVector.z * sinRad;
                    manifold->z[index] = -previousVector.x * sinRad + previousVector.z * cosRad;


                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].x = previousVector.x * cosRad + previousVector.z * sinRad;
                    manifoldNormals[index].z = -previousVector.x * sinRad + previousVector.z * cosRad;
                }
            }
            break;

        case 'z':
            for (int j = 0; j < PIXELS; j++) 
            {
                for (int k = 0; k < PIXELS; k++) 
                {

                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->x[index] = previousVector.x * cosRad - previousVector.y * sinRad;  
                    manifold->y[index] = previousVector.x * sinRad + previousVector.y * cosRad;


                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].x = previousVector.x * cosRad - previousVector.y * sinRad;  
                    manifoldNormals[index].y = previousVector.x * sinRad + previousVector.y * cosRad; 

                }
            }
            break;

        default:
        break;
    }

    return;
}


#endif
