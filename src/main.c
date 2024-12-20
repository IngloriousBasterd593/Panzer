#include "code/mathematika/mathcore.h"
#include "code/mathematika/pipeline.h"
#include "code/utilities/shared.h"
#include "code/utilities/glutilities.h"


int main(int argc, char** argv) {

    SDL_Window* window = NULL;
    SDL_Texture* texture = NULL;
    SDL_Renderer* renderer = NULL;



    if(SDL_init(&window, &renderer, &texture, "3D", SCREENWIDTH, SCREENHEIGHT) == 1) {
        goto quit;
    }

    clock_t start, end;
    double cpu_time_used;

    char* shaderProgram = malloc(BUFFLEN * sizeof(char));
    if(shaderProgram == NULL) {
        fprintf(stderr, "couldn't get space for shader program");
        return 1;
    }

    start = clock();







    Mesh torus;
    Mesh sphere;

    vec3f lightPerspectiveVector = {0, 0, 1};
    vec3f viewVector = {0, 0, 1};

    Camera camera = { viewVector, PI / 3, viewVector, SCREENWIDTH / SCREENHEIGHT, 1, 100, 10, 20, 10, 20 };

    unsigned int* frameColors = NULL;

    // vec3f torusNormals[VERTICES];
    // vec3f sphereNormals[VERTICES];

    vec3f* torusNormals;
    vec3f* sphereNormals;




    if(get_space(&torus) == 1) {
        goto SDLquit;
    }

    if(get_space(&sphere) == 1) {
        goto SDLquit;
    }

    frameColors = malloc(SCREENSIZE * sizeof(unsigned int));
    if(frameColors == NULL) {
        perror("bruh");
        goto free;
    }



    // torus_init(&torus, torusNormals, 100, 50, HALFWINWIDTH, HALFWINHEIGHT, 0);
    sphere_init(&sphere, sphereNormals, 30, HALFWINWIDTH, HALFWINHEIGHT, 600);

    // Mesh_draw(&torus, torusNormals, &camera, frameColors, 0, 0, 20);
    Mesh_draw(&sphere, sphereNormals, &camera, frameColors, 20);

    SDL_RenderPresent(renderer);

    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f milliseconds to init\n", cpu_time_used);



    printf("%f %f %f\n", sphere.velocity.x, sphere.velocity.x, sphere.velocity.x);








    double deltaTime = 0;
    float deltaRad = PI / 540;
    float Rad = 0;
    int radius = 150;
    float theta = 0;
    int drawPrecision = 20;

    int quit = 0;
    SDL_Event e;




    while(!quit) {
        while(SDL_PollEvent(&e) != 0) {
            switch(e.type) {
                case SDL_QUIT:
                    quit = 1; 
                    break;
            }
        } 

        start = clock();

        memset(frameColors, 0xFF, SCREENWIDTH * SCREENHEIGHT * sizeof(unsigned int));

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        Mesh_draw(&sphere, sphereNormals, &camera, frameColors, drawPrecision);
        // Mesh_draw(&torus, torusNormals, frameColors, 200, 130 * sin(Rad), drawPrecision);

        SDL_UpdateTexture(texture, NULL, frameColors, SCREENWIDTH * sizeof(unsigned int));
        SDL_RenderCopy(renderer, texture, NULL, NULL);

        SDL_RenderPresent(renderer);

        // Mesh_rotate(&torus, torusNormals, deltaRad, 'x');
        // Mesh_rotate(&torus, torusNormals, deltaRad, 'y');
        // Mesh_rotate(&torus, torusNormals, deltaRad, 'z');

        Mesh_rotate(&sphere, sphereNormals, deltaRad, 'x');
        Mesh_rotate(&sphere, sphereNormals, deltaRad, 'y');
        Mesh_rotate(&sphere, sphereNormals, deltaRad, 'z');

        theta += 3 * deltaRad;
        Rad += deltaRad;

        end = clock();
        cpu_time_used = ((double) (end - start));

        if(cpu_time_used < 16) {
            usleep(1000 * (16 - cpu_time_used));
        }

        printf("%4.0f milliseconds per frame\n", cpu_time_used);  
    }

    free:

    free(torus.x);
    free(torus.y);
    free(torus.z);

    free(sphere.x);
    free(sphere.y);
    free(sphere.z);

    free(torus.xProj);
    free(torus.yProj);
    free(torus.zProj);

    free(sphere.xProj);
    free(sphere.yProj);
    free(sphere.zProj);

    free(frameColors);
    free(shaderProgram);

    SDLquit: 

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    quit:
    
    return 0;
}
