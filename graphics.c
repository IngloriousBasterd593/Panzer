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




Vector3 normal(Vector3 vector1, Vector3 vector2) {
    Vector3 normal;
    normal.x = vector1.y * vector2.z - vector1.z * vector2.y;
    normal.y = vector1.z * vector2.x - vector1.x * vector2.z;
    normal.z = vector1.x * vector2.y - vector1.y * vector2.x;
    return normal;
}




Vector3 unit(Vector3 vector) {
    float magnitude = sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    vector.x = vector.x / magnitude;
    vector.y = vector.y / magnitude;
    vector.z = vector.z / magnitude;

    return vector;
}




float dotproduct(Vector3 vector1, Vector3 vector2) {
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}




void sphere_init(Manifold* manifold, float radius) {
    
    int r = radius * 25;
    int index = 0;

    for(int j = 0; j < PIXELS; j++) {
        for(int k = 0; k < PIXELS; k++) {
        index = j * PIXELS + k;  

        manifold->x[index] = (int) (r * cos(2 * PI * k / PIXELS) * sin(PI * j / PIXELS));
        manifold->y[index] = (int) (r * sin(2 * PI * k / PIXELS) * sin(PI * j / PIXELS));
        manifold->z[index] = (int) (r * cos(PI * j / PIXELS));
        }
    }

    return;
}




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
}




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




void fillTriangle(SDL_Renderer* renderer, Vector2 vertex1, Vector2 vertex2, Vector2 vertexC, float precision, int R, int G, int B, int Xoffset, int Yoffset) {

    float t = 0;

    for (int i = 0; i < precision; i++) {

        t = (float)i / precision;

        Vector2 v1 = {vertex1.x + t * (vertexC.x - vertex1.x), vertex1.y + t * (vertexC.y - vertex1.y)};
        Vector2 v2 = {vertex2.x + t * (vertexC.x - vertex2.x), vertex2.y + t * (vertexC.y - vertex2.y)};

        line(renderer, v1, v2, R, G, B, Xoffset, Yoffset);

    }

    return;
}



/*
void sphere_draw(SDL_Renderer* renderer, Manifold manifold, int Xoffset, int Yoffset) {

    int index = 0;

    Vector3 Snormal;
    Vector3 v1;
    Vector3 v2;
    Vector3 lightPerspective = {1, 0, 0};

    for(int q = 0; q < PIXELS; q++) {
        for(int l = 0; l < PIXELS; l++) {
            index = q * PIXELS + l;

            v1.x = -sin(l * 2 * PI / PIXELS) * sin(q * PI / PIXELS);
            v1.y = cos(l * 2 * PI / PIXELS) * sin(q * PI / PIXELS);
            v1.z = 0;
            v2.x = cos(l * 2 * PI / PIXELS) * cos(q * PI / PIXELS);
            v2.y = sin(l * 2 * PI / PIXELS) * sin(q * PI / PIXELS);
            v2.z = -sin(q * PI / PIXELS);

            Snormal = normal(v1, v2);
            Snormal = unit(Snormal);
            float melnums = dotproduct(Snormal, lightPerspective);

            if(melnums < 0) {
                melnums = -melnums;
            }

            int adaskrasa = 255 - (196 * (1 - melnums)); 

            Vector2 p1 = {manifold.x[index], manifold.y[index]};
            Vector2 p2 = {manifold.x[index + 1], manifold.y[index + 1]};
            Vector2 pC = {manifold.x[index + PIXELS], manifold.y[index + PIXELS]};
            Vector2 pA = {manifold.x[index + PIXELS + 1], manifold.y[index + PIXELS + 1]};

            fillTriangle(renderer, p1, p2, pC, 5, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset); // 20 precision

        }
    }

    SDL_RenderPresent(renderer);
    return;
}  */




void torus_draw(SDL_Renderer* renderer, Manifold* manifold, int Xoffset, int Yoffset, int precision) { 

    int index = 0;
    int adaskrasa = 0;

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
            }

            v1.x = -(1 + cos(q * piOverPixels)) * sin(l * twopiOverPixels);
            v1.y = (1 + cos(q * piOverPixels)) * cos(l * twopiOverPixels);
            v1.z = 0;

            v2.x = -sin(q * piOverPixels) * cos(l * twopiOverPixels);
            v2.y = -sin(q * piOverPixels) * sin(l * twopiOverPixels);
            v2.z = cos(q * piOverPixels); 

            Snormal = normal(v1, v2);
            Snormal = unit(Snormal);
            float melnums = dotproduct(Snormal, lightPerspective);

            if(melnums < 0) {
                melnums = -melnums;
            } 

            adaskrasa = 255 - (196 * (1 - melnums)); 

            vertex1.x = manifold->x[index];
            vertex1.y = manifold->y[index];
            vertex2.x = manifold->x[index + 1];
            vertex2.y = manifold->y[index + 1];
            vertexU1.x = manifold->x[index + PIXELS];
            vertexU1.y = manifold->y[index + PIXELS];
            vertexU2.x = manifold->x[index + PIXELS + 1];
            vertexU2.y = manifold->y[index + PIXELS + 1];

            fillTriangle(renderer, vertex1, vertex2, vertexU1, precision, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset);
            fillTriangle(renderer, vertex2, vertexU1, vertexU2, precision, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset);
            
        }
    }

    return;
} 




void M_rotate(Manifold* manifold, float rad, char axis) {
    float X = 0;
    float Y = 0;
    float Z = 0;

    float cosRad = (float) cos(rad);
    float sinRad = (float) sin(rad);

    int index = 0;
    switch(axis) {
        case 'x':
            for(int j = 0; j < PIXELS; j++) {
                for(int k = 0; k < PIXELS; k++) {
                    index = j * PIXELS + k;
                    X = manifold->x[index];
                    Y = manifold->y[index];
                    Z = manifold->z[index];

                    manifold->x[index] = X;
                    manifold->y[index] = Y * cosRad - Z * sinRad;
                    manifold->z[index] = Y * sinRad + Z * cosRad;

                }
            }
        break;

        case 'y':
            for(int j = 0; j < PIXELS; j++) {
                for(int k = 0; k < PIXELS; k++) {
                    index = j * PIXELS + k;
                    X = manifold->x[index];
                    Y = manifold->y[index];
                    Z = manifold->z[index];

                    manifold->x[index] = X * cosRad + Z * sinRad;
                    manifold->y[index] = Y;
                    manifold->z[index] = -X * sinRad + Z * cosRad;
                }   
            }
        break;

        case 'z':
            for(int j = 0; j < PIXELS; j++) {
                for(int k = 0; k < PIXELS; k++) {
                    X = manifold->x[index];
                    Y = manifold->y[index];
                    Z = manifold->z[index];
                    index = j * PIXELS + k;

                    manifold->x[index] = X * cosRad - Y * sinRad;
                    manifold->y[index] = X * sinRad + Y * cosRad;
                    manifold->z[index] = Z;
                }   
            }
        break;
        default:
        break;
    }

    return;
}


