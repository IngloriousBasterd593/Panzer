#ifndef SDL_GRAPHICS_LIB
#define SDL_GRAPHICS_LIB

// external libs

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
// #include <unistd.h>
#include <time.h>


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

#define PI 3.141592653589793
#define PIXELS 40           
#define POINTS PIXELS * PIXELS
#define TWOPI 2 * PI
#define TWOPIOVERPIXELS TWOPI / (PIXELS - 1)
#define PIOVERPIXELS PI / (PIXELS - 1)
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


void SDL_init(SDL_Window** window, SDL_Renderer** renderer, SDL_Texture** texture, const char* name, int width, int height) {

    if(SDL_Init(SDL_INIT_VIDEO) < 0) {
        perror("SDL lmao could not initialize!");
        return;
    }

    *window = SDL_CreateWindow(name, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
    if(*window == NULL) {
        perror("window null lmfao");
        SDL_Quit();
        return;
    }

    *renderer = SDL_CreateRenderer(*window, -1, SDL_RENDERER_ACCELERATED);
    if(*renderer == NULL) {
        perror("renderer null bruh");
        SDL_DestroyWindow(*window);
        SDL_Quit();
        return;
    }

    *texture = SDL_CreateTexture(*renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_TARGET, width, height);
    if(*texture == NULL) {
        perror("texture null bruh");
        SDL_DestroyWindow(*window);
        SDL_DestroyRenderer(*renderer);
        SDL_Quit();
        return;
    }


    SDL_SetRenderDrawColor(*renderer, 255, 255, 255, 255);
    SDL_RenderClear(*renderer);
    SDL_RenderPresent(*renderer);
    SDL_SetRenderTarget(*renderer, *texture);
    SDL_RenderPresent(*renderer);

    return;
}




int get_space(Manifold* manifold, unsigned int** colors) {

    manifold->x = malloc(POINTS * sizeof(float));
    if(manifold->x == NULL) {
        perror("bruh");
        return 1;
    }

    manifold->y = malloc(POINTS * sizeof(float));
    if(manifold->y == NULL) {
        free(manifold->x);
        perror("bruh");
        return 1;
    }

    manifold->z = malloc(POINTS * sizeof(float));
    if(manifold->z == NULL) {
        free(manifold->x);
        free(manifold->y);
        perror("bruh");
        return 1;
    }

    
    *colors = malloc(S_WIDTH * S_HEIGHT * sizeof(unsigned int));
    if(*colors == NULL) {
        free(manifold->x);
        free(manifold->y);
        free(manifold->z);
        perror("bruh");
        return 1;
    }


    for(int i = 0; i < S_WIDTH * S_HEIGHT; i++) {
        (*colors)[i] = (255 << 24) | (255 << 16) | (255 << 8) | 255; 
    }
    


    return 0;
}




Vector3 unit(Vector3* v) {

    double magnitude = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
    Vector3 normalized;

    if (magnitude == 0) {
        normalized.x = normalized.y = normalized.z = 0;
        return normalized;
    } 

    normalized.x = v->x / magnitude;
    normalized.y = v->y / magnitude;
    normalized.z = v->z / magnitude;

    return normalized;
}




Vector3 crossproduct(Vector3* vector1, Vector3* vector2) {
    Vector3 Vnormal;

    Vnormal.x = vector1->y * vector2->z - vector1->z * vector2->y;
    Vnormal.y = vector1->z * vector2->x - vector1->x * vector2->z;
    Vnormal.z = vector1->x * vector2->y - vector1->y * vector2->x;

    return Vnormal;
}




float dotproduct(Vector3* vector1, Vector3* vector2) {
    return vector1->x * vector2->x + vector1->y * vector2->y + vector1->z * vector2->z;
}




void sphere_init(Manifold* manifold, Vector3 manifoldNormals[], float radius) {

    if(radius <= 0) {
        perror("bozo");
        return;
    }
    
    int r = radius * 25;
    int index = 0;

    Vector3 normalVector;
    Vector3 partialDerivativeU;
    Vector3 partialDerivativeV;
    
    for(int j = 0; j < PIXELS; j++) {
        for(int k = 0; k < PIXELS; k++) {
        index = j * PIXELS + k;  

        manifold->x[index] = (r * cos(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j));
        manifold->y[index] = (r * sin(TWOPIOVERPIXELS * k) * sin(PIOVERPIXELS * j));
        manifold->z[index] = (r * cos(PIOVERPIXELS * j));
        
        
        partialDerivativeU.x = -r * sin(j * TWOPIOVERPIXELS) * sin(k * PIOVERPIXELS);
        partialDerivativeU.y = r * cos(j * TWOPIOVERPIXELS) * sin(k * PIOVERPIXELS);
        partialDerivativeU.z = 0;

        partialDerivativeV.x = r * cos(j * TWOPIOVERPIXELS) * cos(k * PIOVERPIXELS);
        partialDerivativeV.y = r * sin(j * TWOPIOVERPIXELS) * cos(k * PIOVERPIXELS);
        partialDerivativeV.z = -r * sin(k * PIOVERPIXELS);

        normalVector = crossproduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit(&normalVector);

        }
    }

    return;
}  




void torus_init(Manifold* manifold, Vector3 manifoldNormals[], float radiusInner, float radiusOuter) {

    if(radiusInner <= 0 || radiusOuter <= 0) {
        perror("bozoo");
        return;
    }

    int Rinner = radiusInner * 25;
    int Router = radiusOuter * 50;

    Vector3 normalVector;
    Vector3 partialDerivativeU;
    Vector3 partialDerivativeV;

    int index = 0;

    for(int j = 0; j < PIXELS; j++) {
        for(int k = 0; k < PIXELS; k++) {
        index = j * PIXELS + k; 

        manifold->x[index] = ((Router + Rinner * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j));
        manifold->y[index] = ((Router + Rinner * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j));
        manifold->z[index] = (Rinner * sin(TWOPIOVERPIXELS * k));


        partialDerivativeU.x = ((Router + Rinner * cos(TWOPIOVERPIXELS * k)) * cos(TWOPIOVERPIXELS * j));
        partialDerivativeU.y = (-(Router + Rinner * cos(TWOPIOVERPIXELS * k)) * sin(TWOPIOVERPIXELS * j));
        partialDerivativeU.z = 0;

        partialDerivativeV.x = (-Rinner * sin(TWOPIOVERPIXELS * k) * cos(TWOPIOVERPIXELS * j));
        partialDerivativeV.y = (-Rinner * sin(TWOPIOVERPIXELS * k) * sin(TWOPIOVERPIXELS * j));
        partialDerivativeV.z = (Rinner * cos(TWOPIOVERPIXELS * k));

        normalVector = crossproduct(&partialDerivativeU, &partialDerivativeV);

        manifoldNormals[index] = unit(&normalVector);

        }
    }

    return;
} 




void lineBresenham(SDL_Renderer* renderer, unsigned int* frameColors, Vector2 vertexStart, Vector2 vertexEnd, int deltaX, int deltaY, uint8_t R, uint8_t G, uint8_t B, uint8_t A) {
    
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

    while(1) {

        frameColors[(y1 + deltaY) * S_WIDTH + x1 + deltaX] = (R << 24) | (G << 16) | (B << 8) | A;

        // SDL_RenderDrawPoint(renderer, x1 + deltaX, y1 + deltaY);

        if (x1 == x2 && y1 == y2) {
            break;
        }

        error2 = 2 * error;

        if (error2 > -dy) {
            error -= dy;
            x1 += stepX;
        }
        if (error2 < dx) {
            error += dx;
            y1 += stepY;
        }
    }
}




void fillTriangle(SDL_Renderer* renderer, unsigned int* frameColors, Vector2 baseVertex1, Vector2 baseVertex2, Vector2 centralVertex, int trianglePrecision, uint8_t R, uint8_t G, uint8_t B, uint8_t A, int Xoffset, int Yoffset) {

    int deltaX = C_winX + Xoffset;
    int deltaY = C_winY + Yoffset;

    float dx1 = (centralVertex.x - baseVertex1.x) / trianglePrecision;
    float dy1 = (centralVertex.y - baseVertex1.y) / trianglePrecision;
    float dx2 = (centralVertex.x - baseVertex2.x) / trianglePrecision;
    float dy2 = (centralVertex.y - baseVertex2.y) / trianglePrecision;

    Vector2 point1;
    Vector2 point2;

    for (int i = 0; i < trianglePrecision; i++) {

        point1.x = baseVertex1.x + dx1;
        point1.y = baseVertex1.y + dy1;

        point2.x = baseVertex2.x + dx2;
        point2.y = baseVertex2.y + dy2;

        lineBresenham(renderer, frameColors, point1, point2, deltaX, deltaY, R, G, B, A);

    }

    return;
}




void fillRectangle(SDL_Renderer* renderer, unsigned int* frameColors, Vector2 vertexA, Vector2 vertexB, Vector2 vertexC, Vector2 vertexD, int drawPrecision, uint8_t R, uint8_t G, uint8_t B, uint8_t A, int deltaX, int deltaY) {
 
    if(drawPrecision <= 0) {
        perror("rip bozo");
        return;
    }

    SDL_SetRenderDrawColor(renderer, R, G, B, 255);

    float dxU = (vertexD.x - vertexC.x) / drawPrecision;
    float dyU = (vertexD.y - vertexC.y) / drawPrecision;
    float dxL = (vertexB.x - vertexA.x) / drawPrecision;
    float dyL = (vertexB.y - vertexA.y) / drawPrecision;

    for(int i = 0; i <= drawPrecision; i++) {

        vertexA.x += dxU;
        vertexA.y += dyU;

        vertexC.x += dxL;
        vertexC.y += dyL;

        lineBresenham(renderer, frameColors, vertexA, vertexC, deltaX, deltaY, R, G, B, A);
    }

    return;
}




void Manifold_draw(SDL_Renderer* renderer, Manifold* manifold, Vector3 manifoldNormals[], unsigned int* frameColors, int Xoffset, int Yoffset, int drawPrecision) {

    if(drawPrecision < 1) {
        perror("you twisted sack of shit");
        return;
    }

    int index;
    int indexPlusPixelsPlusOne;
    int indexPlusPixels;
    int indexPlusOne;
    int nextL;
    int nextQ;

    int deltaX = C_winX + Xoffset;
    int deltaY = C_winY + Yoffset;

    int grayscaleRGB;
    float grayscaleCoefficient;

    Vector3 lightPerspectiveVector = {1, 0, 0};

    Vector2 vertex1;
    Vector2 vertex2;
    Vector2 vertexUpper1;
    Vector2 vertexUpper2;

    for(int q = 0; q < PIXELS; q++) {
        for(int l = 0; l < PIXELS; l++) {
            
            nextL = (l + 1) % PIXELS;
            nextQ = (q + 1) % PIXELS;

            index = q * PIXELS + l;
            indexPlusOne = q * PIXELS + nextL;
            indexPlusPixels = nextQ * PIXELS + l;
            indexPlusPixelsPlusOne = nextQ * PIXELS + nextL;

            grayscaleCoefficient = fabs(dotproduct(&manifoldNormals[index], &lightPerspectiveVector));

            grayscaleRGB = (uint8_t) (255 - (90 * (1 - grayscaleCoefficient)));

            vertex1.x = manifold->y[index];
            vertex1.y = manifold->z[index];

            vertex2.x = manifold->y[indexPlusOne];
            vertex2.y = manifold->z[indexPlusOne];

            vertexUpper1.x = manifold->y[indexPlusPixels];
            vertexUpper1.y = manifold->z[indexPlusPixels];

            vertexUpper2.x = manifold->y[indexPlusPixelsPlusOne];
            vertexUpper2.y = manifold->z[indexPlusPixelsPlusOne];   

            fillRectangle(renderer, frameColors, vertex1, vertex2, vertexUpper1, vertexUpper2, drawPrecision, grayscaleRGB, grayscaleRGB, grayscaleRGB, 255, deltaX, deltaY);

        }
    }

    return;
}  




void Manifold_rotate(Manifold* manifold, Vector3 manifoldNormals[], float rad, char axis) {

    Vector3 previousVector;

    float cosRad = (float) cos(rad);
    float sinRad = (float) sin(rad);

    int index = 0;
    
    switch(axis) {
        case 'x':
            for (int j = 0; j < PIXELS; j++) {
                for (int k = 0; k < PIXELS; k++) {
                    // rotēt manifolda punktus
                    index = j * PIXELS + k;

                    previousVector.x = manifold->x[index];
                    previousVector.y = manifold->y[index];
                    previousVector.z = manifold->z[index];

                    manifold->y[index] = previousVector.y * cosRad - previousVector.z * sinRad;
                    manifold->z[index] = previousVector.y * sinRad + previousVector.z * cosRad;

                    // rotēt normālvektorus

                    previousVector = manifoldNormals[index];

                    manifoldNormals[index].y = previousVector.y * cosRad - previousVector.z * sinRad;
                    manifoldNormals[index].z = previousVector.y * sinRad + previousVector.z * cosRad;

                }
            }
        break;

        case 'y':
            for (int j = 0; j < PIXELS; j++) {
                for (int k = 0; k < PIXELS; k++) {
                    
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
            for (int j = 0; j < PIXELS; j++) {
                for (int k = 0; k < PIXELS; k++) {

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