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




void get_space(Manifold* manifold) {

    manifold->x = malloc(POINTS * sizeof(float));
    if(manifold->x == NULL) {
        perror("bruh");
        return;
    }

    manifold->y = malloc(POINTS * sizeof(float));
    if(manifold->y == NULL) {
        perror("bruh");
        return;
    }

    manifold->z = malloc(POINTS * sizeof(float));
    if(manifold->z == NULL) {
        perror("bruh");
        return;
    }

    return;
}




Vector3 crossproduct(Vector3 vector1, Vector3 vector2) {

    Vector3 Vnormal;

    Vnormal.x = vector1.y * vector2.z - vector1.z * vector2.y;
    Vnormal.y = vector1.z * vector2.x - vector1.x * vector2.z;
    Vnormal.z = vector1.x * vector2.y - vector1.y * vector2.x;

    return Vnormal;
}




Vector3 unit(Vector3 vector) {
    float magnitude = sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);

    if(magnitude == 0) {
        // perror("divided by zero kys");
        return vector;
    }

    vector.x /= magnitude;
    vector.y /= magnitude;
    vector.z /= magnitude;

    return vector;
}




float dotproduct(Vector3 vector1, Vector3 vector2) {
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}




void sphere_init(Manifold* manifold, Vector3 manifoldNormals[], float radius) {
    
    int r = radius * 25;
    int index = 0;

    Vector3 partialDerivativeU;
    Vector3 partialDerivativeV;
    
    for(int j = 0; j < PIXELS; j++) {
        for(int k = 0; k < PIXELS; k++) {
        index = j * PIXELS + k;  

        manifold->x[index] = (r * cos(twopiOverPixels * k) * sin(piOverPixels * j));
        manifold->y[index] = (r * sin(twopiOverPixels * k) * sin(piOverPixels * j));
        manifold->z[index] = (r * cos(piOverPixels * j));
        
        partialDerivativeU.x = -sin(j * twopiOverPixels) * sin(k * piOverPixels);
        partialDerivativeU.y = cos(j * twopiOverPixels) * sin(k * piOverPixels);
        partialDerivativeU.z = 0;

        partialDerivativeV.x = cos(j * twopiOverPixels) * cos(k * piOverPixels);
        partialDerivativeV.y = sin(j * twopiOverPixels) * sin(k * piOverPixels);
        partialDerivativeV.z = -sin(k * piOverPixels);

        manifoldNormals[index] = unit(crossproduct(partialDerivativeU, partialDerivativeV));

        }
    }

    return;
}  



/*
void torus_init(Manifold* manifold, float radiusInner, float radiusOuter) {

    int Rinner = radiusInner * 25;
    int Router = radiusOuter * 50;

    int index = 0;

    for(int j = 0; j < PIXELS; j++) {
        for(int k = 0; k < PIXELS; k++) {
        index = j * PIXELS + k;  

        manifold->x[index] = (int) ((Router + Rinner * cos(twopiOverPixels * k)) * sin(twopiOverPixels * j));
        manifold->y[index] = (int) ((Router + Rinner * cos(twopiOverPixels * k)) * cos(twopiOverPixels * j));
        manifold->z[index] = (int) (Rinner * sin(twopiOverPixels * k));

        }
    }

    return;
} */




void line(SDL_Renderer* renderer, Vector2 vector1, Vector2 vector2, int R, int G, int B, int Xoffset, int Yoffset) {

    double dx = vector2.x - vector1.x;
    double dy = vector2.y - vector1.y;

    float s = sqrt((dx * dx) + (dy * dy));
    double theta = atan2(dy, dx);

    float cosT = (float) cos(theta);
    float sinT = (float) sin(theta);

    SDL_SetRenderDrawColor(renderer, R, G, B, 255);

    for (int i = 0; i <= (int) s; i++) {
        SDL_RenderDrawPoint(renderer, (int) (vector1.x + i * cosT) + C_winX + Xoffset, (int) (vector1.y + i * sinT) + C_winY + Yoffset);
    }

    return;
}




void fillTriangle(SDL_Renderer* renderer, Vector2 baseVertex1, Vector2 baseVertex2, Vector2 centralVertex, float trianglePrecision, int R, int G, int B, int Xoffset, int Yoffset) {

    float t = 0;

    Vector2 point1;
    Vector2 point2;

    for (int i = 0; i < trianglePrecision; i++) {

        t = (float)i / trianglePrecision;

        point1.x = baseVertex1.x + t * (centralVertex.x - baseVertex1.x);
        point1.y = baseVertex1.y + t * (centralVertex.y - baseVertex1.y);

        point1.x = baseVertex2.x + t * (centralVertex.x - baseVertex2.x);
        point1.y = baseVertex2.y + t * (centralVertex.y - baseVertex2.y);

        line(renderer, point1, point2, R, G, B, Xoffset, Yoffset);

    }

    return;
}




void fillRectangle(SDL_Renderer* renderer, Vector2 vertexA, Vector2 vertexB, Vector2 vertexC, Vector2 vertexD, float drawPrecision, int R, int G, int B) {






    return;
}




void sphere_draw(SDL_Renderer* renderer, Manifold* manifold, Vector3 manifoldNormals[], int Xoffset, int Yoffset, int trianglePrecision) {

    int index = 0;
    int grayscaleRGB = 0;
    float grayscaleCoefficient = 0;

    Vector3 surfaceNormalVector;
    Vector3 lightPerspectiveVector = {1, 0, 0};

    Vector2 vertex1;
    Vector2 vertex2;
    Vector2 vertexUpper1;
    Vector2 vertexUpper2;

    for(int q = 0; q < PIXELS; q++) {
        for(int l = 0; l < PIXELS; l++) {
            index = q * PIXELS + l;

            if(index + PIXELS > POINTS) {
                continue;
            }   // broken


            surfaceNormalVector = manifoldNormals[index];

            grayscaleCoefficient = dotproduct(surfaceNormalVector, lightPerspectiveVector);

            grayscaleRGB = (255 - (90 * (1 - grayscaleCoefficient)));

            vertex1.x = manifold->y[index];
            vertex1.y = manifold->z[index];

            vertex2.x = manifold->y[index + 1];
            vertex2.y = manifold->z[index + 1];

            vertexUpper1.x = manifold->y[index + PIXELS];
            vertexUpper1.y = manifold->z[index + PIXELS];

            vertexUpper2.x = manifold->y[index + PIXELS + 1];
            vertexUpper2.y = manifold->z[index + PIXELS + 1];   // out of bounds - salabot

            fillTriangle(renderer, vertex1, vertex2, vertexUpper1, trianglePrecision, grayscaleRGB, grayscaleRGB, grayscaleRGB, Xoffset, Yoffset);
            fillTriangle(renderer, vertex2, vertexUpper1, vertexUpper2, trianglePrecision, grayscaleRGB, grayscaleRGB, grayscaleRGB, Xoffset, Yoffset);

        }
    }

    return;
}  



/*
void torus_draw(SDL_Renderer* renderer, Manifold* manifold, float Rinner, float Router, int Xoffset, int Yoffset, int precision) { 

    int index = 0;
    int adaskrasa = 0;
    float melnums = 0;

    Rinner = Rinner * 25;
    Router = Router * 50;

    Vector3 Snormal;
    Vector3 v1;
    Vector3 v2;
    Vector3 lightPerspective = {0, 1, 0};

    Vector2 vertex1;
    Vector2 vertex2;
    Vector2 vertexU1;
    Vector2 vertexU2;

    for(int q = 0; q < PIXELS; q++) {
        for(int l = 0; l < PIXELS; l++) {
            index = q * PIXELS + l;

            if(index + PIXELS > POINTS) {
                continue;
            }   // broken

            v1.x =  ((Router + Rinner * cos(twopiOverPixels * l)) * cos(twopiOverPixels * q));
            v1.y =  (-(Router + Rinner * cos(twopiOverPixels * l)) * sin(twopiOverPixels * q));
            v1.z = 0;  

            v2.x =  (-Rinner * sin(twopiOverPixels * l) * cos(twopiOverPixels * q));
            v2.y =  (-Rinner * sin(twopiOverPixels * l) * sin(twopiOverPixels * q));
            v2.z =  (Rinner * cos(twopiOverPixels * l));

            Snormal = normal(v1, v2);
            Snormal = unit(Snormal);

            melnums = fabs(dotproduct(Snormal, lightPerspective));
            

            adaskrasa = (255 - (62 * (1 - melnums))); 

            vertex1.x = manifold->x[index];
            vertex1.y = manifold->y[index];

            vertex2.x = manifold->x[index + 1];
            vertex2.y = manifold->y[index + 1];

            vertexU1.x = manifold->x[index + PIXELS];
            vertexU1.y = manifold->y[index + PIXELS];

            vertexU2.x = manifold->x[index + PIXELS + 1];
            vertexU2.y = manifold->y[index + PIXELS + 1];   // out of bounds - salabot

            fillTriangle(renderer, vertex1, vertex2, vertexU1, precision, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset);
            fillTriangle(renderer, vertex2, vertexU1, vertexU2, precision, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset);
            
        }
    }

    return;
} */




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









