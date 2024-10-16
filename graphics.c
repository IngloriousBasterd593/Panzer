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




void sphere_init(Manifold* manifold, SurfaceNormals* manifoldNormals, float radius) {
    
    int r = radius * 25;
    int index = 0;

    for(int j = 0; j < PIXELS; j++) {
        for(int k = 0; k < PIXELS; k++) {
        index = j * PIXELS + k;  

        manifold->x[index] = (int) (r * cos(twopiOverPixels * k) * sin(piOverPixels * j));
        manifold->y[index] = (int) (r * sin(twopiOverPixels * k) * sin(piOverPixels * j));
        manifold->z[index] = (int) (r * cos(piOverPixels * j));
        
        manifoldNormals->u[index].x = -sin(j * twopiOverPixels) * sin(k * piOverPixels);
        manifoldNormals->u[index].y = cos(j * twopiOverPixels) * sin(k * piOverPixels);
        manifoldNormals->u[index].z = 0;

        manifoldNormals->u[index].x = cos(j * twopiOverPixels) * cos(k * piOverPixels);
        manifoldNormals->u[index].y = sin(j * twopiOverPixels) * sin(k * piOverPixels);
        manifoldNormals->u[index].z = -sin(k * piOverPixels);

        printf("%f\n", manifoldNormals->u[index].x);
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




void fillRectangle(SDL_Renderer* renderer, Vector2 vertexA, Vector2 vertexB, Vector2 vertexC, Vector2 vertexD, float precision, int R, int G, int B) {






    return;
}




void sphere_draw(SDL_Renderer* renderer, Manifold* manifold, SurfaceNormals* manifoldNormals, int Xoffset, int Yoffset, int precision) {

    int index = 0;
    int adaskrasa = 0;
    float melnums = 0;

    Vector3 Snormal;
    Vector3 v1;
    Vector3 v2;
    Vector3 lightPerspective = {0, 0, 1};

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

            /*

            v1.x = -sin(l * twopiOverPixels) * sin(q * piOverPixels);
            v1.y = cos(l * twopiOverPixels) * sin(q * piOverPixels);
            v1.z = 0;

            v2.x = cos(l * twopiOverPixels) * cos(q * piOverPixels);
            v2.y = sin(l * twopiOverPixels) * sin(q * piOverPixels);
            v2.z = -sin(q * piOverPixels);

            Snormal = normal(v1, v2);
            Snormal = unit(Snormal);

            melnums = dotproduct(Snormal, lightPerspective);

            adaskrasa = (255 - (126 * (1 - melnums))); */

            v1 = manifoldNormals->u[index];
            v2 = manifoldNormals->v[index];

            Snormal = normal(v1, v2);
            Snormal = unit(Snormal);

            melnums = dotproduct(Snormal, lightPerspective);

            adaskrasa = (255 - (126 * (1 - melnums)));

            vertex1.x = manifold->y[index];
            vertex1.y = manifold->z[index];

            vertex2.x = manifold->y[index + 1];
            vertex2.y = manifold->z[index + 1];

            vertexU1.x = manifold->y[index + PIXELS];
            vertexU1.y = manifold->z[index + PIXELS];

            vertexU2.x = manifold->y[index + PIXELS + 1];
            vertexU2.y = manifold->z[index + PIXELS + 1];   // out of bounds - salabot

            fillTriangle(renderer, vertex1, vertex2, vertexU1, precision, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset);
            fillTriangle(renderer, vertex2, vertexU1, vertexU2, precision, adaskrasa, adaskrasa, adaskrasa, Xoffset, Yoffset);

        }
    }

    return;
}  




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

            melnums = dotproduct(Snormal, lightPerspective);

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




void V_rotate(SurfaceNormals* manifoldNormals, float rad, char axis) {

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

                    X = manifoldNormals->u[index].x;
                    Y = manifoldNormals->u[index].y;
                    Z = manifoldNormals->u[index].z;

                    manifoldNormals->u[index].x = X;
                    manifoldNormals->u[index].y = Y * cosRad - Z * sinRad;
                    manifoldNormals->u[index].z = Y * sinRad + Z * cosRad;

                    X = manifoldNormals->v[index].x;
                    Y = manifoldNormals->v[index].y;
                    Z = manifoldNormals->v[index].z;

                    manifoldNormals->v[index].x = X;
                    manifoldNormals->v[index].y = Y * cosRad - Z * sinRad;
                    manifoldNormals->v[index].z = Y * sinRad + Z * cosRad;

                }
            }
        break;

        case 'y':
            for(int j = 0; j < PIXELS; j++) {
                for(int k = 0; k < PIXELS; k++) {
                    index = j * PIXELS + k;

                    X = manifoldNormals->u[index].x;
                    Y = manifoldNormals->u[index].y;
                    Z = manifoldNormals->u[index].z;

                    manifoldNormals->u[index].x = X * cosRad + Z * sinRad;
                    manifoldNormals->u[index].y = Y;
                    manifoldNormals->u[index].z = -X * sinRad + Z * cosRad;

                    X = manifoldNormals->v[index].x;
                    Y = manifoldNormals->v[index].y;
                    Z = manifoldNormals->v[index].z;

                    manifoldNormals->v[index].x = X * cosRad + Z * sinRad;
                    manifoldNormals->v[index].y = Y;
                    manifoldNormals->v[index].z = -X * sinRad + Z * cosRad;

                }   
            }
        break;

        case 'z':
            for(int j = 0; j < PIXELS; j++) {
                for(int k = 0; k < PIXELS; k++) {
                    
                    X = manifoldNormals->u[index].x;
                    Y = manifoldNormals->u[index].y;
                    Z = manifoldNormals->u[index].z;

                    manifoldNormals->u[index].x = X * cosRad - Y * sinRad;
                    manifoldNormals->u[index].y = X * sinRad + Y * cosRad;
                    manifoldNormals->u[index].z = Z;

                    X = manifoldNormals->v[index].x;
                    Y = manifoldNormals->v[index].y;
                    Z = manifoldNormals->v[index].z;

                    manifoldNormals->v[index].x = X * cosRad - Y * sinRad;
                    manifoldNormals->v[index].y = X * sinRad + Y * cosRad;
                    manifoldNormals->v[index].z = Z;

                }   
            }
        break;

        default:
        break;
    }

    return;
}






