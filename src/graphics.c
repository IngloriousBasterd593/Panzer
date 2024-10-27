#include "sdl_lib.h"



SDL_Renderer* SDL_INIT(SDL_Window** window, const char* name, int width, int height) {

    if(SDL_Init(SDL_INIT_VIDEO) < 0) {
        perror("SDL lmao could not initialize!");
        return NULL;
    }

    *window = SDL_CreateWindow(name, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
    if(*window == NULL) {
        perror("window null lmfao");
        SDL_Quit();
        return NULL;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(*window, -1, SDL_RENDERER_ACCELERATED);
    if(renderer == NULL) {
        perror("renderer null bruh");
        SDL_DestroyWindow(*window);
        SDL_Quit();
        return NULL;
    }

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);
    SDL_RenderPresent(renderer);

    return renderer;
}




int get_space(Manifold* manifold) {

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




void line(SDL_Renderer* renderer, Vector2 vertexStart, Vector2 vertexEnd, int deltaX, int deltaY, int linePrecision) {
    // slow
    float dx = (vertexEnd.x - vertexStart.x) / linePrecision;
    float dy = (vertexEnd.y - vertexStart.y) / linePrecision;

    for (int i = 0; i <= linePrecision; i++) {

        vertexStart.x += dx;
        vertexStart.y += dy;
        
        SDL_RenderDrawPoint(renderer, (int) vertexStart.x + deltaX, (int) vertexStart.y + deltaY);
    }

    return;
}




void lineNew(SDL_Renderer* renderer, Vector2 vertexStart, Vector2 vertexEnd, int deltaX, int deltaY) {

    short delta = floor((vertexEnd.y - vertexStart.y) / (vertexEnd.x - vertexStart.x));
    short stepX = 0;
    short stepY = 0;

    if(delta < 0) {
        delta--;

        for(int i = 0; i < vertexEnd.y - vertexStart.y; i++) {
        SDL_RenderDrawPoint(renderer, stepX + deltaX, stepY + deltaY);

        stepY++;
        stepY % delta == 0 ? stepX--;

        }

    } else {
        delta++;

        for(int i = 0; i < vertexEnd.y - vertexStart.y; i++) 
        SDL_RenderDrawPoint(renderer, stepX + deltaX, stepY + deltaY);  

        stepY++;
        stepY % delta == 0 ? stepX++;

        }

    return;
}




void fillTriangle(SDL_Renderer* renderer, Vector2 baseVertex1, Vector2 baseVertex2, Vector2 centralVertex, int trianglePrecision, int linePrecision, int R, int G, int B, int Xoffset, int Yoffset) {

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

        line(renderer, point1, point2, deltaX, deltaY, linePrecision);

    }

    return;
}




void fillRectangle(SDL_Renderer* renderer, Vector2 vertexA, Vector2 vertexB, Vector2 vertexC, Vector2 vertexD, int drawPrecision, int linePrecision, int R, int G, int B, int Xoffset, int Yoffset) {

    if(drawPrecision <= 0) {
        perror("rip bozo");
        return;
    }

    SDL_SetRenderDrawColor(renderer, R, G, B, 255);

    int deltaX = C_winX + Xoffset;
    int deltaY = C_winY + Yoffset;

    float dxU = (vertexD.x - vertexC.x) / drawPrecision;
    float dyU = (vertexD.y - vertexC.y) / drawPrecision;
    float dxL = (vertexB.x - vertexA.x) / drawPrecision;
    float dyL = (vertexB.y - vertexA.y) / drawPrecision;

    for(int i = 0; i <= drawPrecision; i++) {

        vertexA.x += dxU;
        vertexA.y += dyU;

        vertexC.x += dxL;
        vertexC.y += dyL;

        line(renderer, vertexA, vertexC, deltaX, deltaY, 10);
    }

    return;
}




void Manifold_draw(SDL_Renderer* renderer, Manifold* manifold, Vector3 manifoldNormals[], int Xoffset, int Yoffset, int drawPrecision, int linePrecision) {

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

            grayscaleRGB = (255 - (90 * (1 - grayscaleCoefficient)));

            vertex1.x = manifold->y[index];
            vertex1.y = manifold->z[index];

            vertex2.x = manifold->y[indexPlusOne];
            vertex2.y = manifold->z[indexPlusOne];

            vertexUpper1.x = manifold->y[indexPlusPixels];
            vertexUpper1.y = manifold->z[indexPlusPixels];

            vertexUpper2.x = manifold->y[indexPlusPixelsPlusOne];
            vertexUpper2.y = manifold->z[indexPlusPixelsPlusOne];   

            fillRectangle(renderer, vertex1, vertex2, vertexUpper1, vertexUpper2, drawPrecision, linePrecision, grayscaleRGB, grayscaleRGB, grayscaleRGB, Xoffset, Yoffset);

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