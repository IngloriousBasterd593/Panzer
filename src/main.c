#include "code/mathematika/mathcore.h"
#include "code/mathematika/pipeline.h"
#include "code/utilities/shared.h"
#include "code/utilities/glutilities.h"


int main(int argc, char** argv) {

    SDL_Window** window = NULL;
    SDL_Texture** texture = NULL;
    SDL_Renderer** renderer = NULL;



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

    int meshCount = 1;

    Mesh** meshes = malloc(meshCount * sizeof(Mesh*));
    if(meshes == NULL) 
    {
        fprintf(stderr, "couldn't get space for meshes");
        goto SDLquit;
    }

    for(int i = 0; i < meshCount; i++) 
    {
        meshes[i] = malloc(sizeof(Mesh));
        if(meshes[i] == NULL) 
        {
            fprintf(stderr, "couldn't get space for mesh");
            goto free;
        }

        if(get_space(meshes[i]) == 1) 
        {
            goto free;
        }
    }

    vec3f lightPerspectiveVector = {0, 0, 1};
    vec3f viewVector = {0, 0, 1};

    Camera camera = { viewVector, PI / 3, viewVector, SCREENWIDTH / SCREENHEIGHT, 1, 100, 10, 20, 10, 20 };

    unsigned int* frameColors = NULL;

    frameColors = malloc(SCREENSIZE * sizeof(unsigned int));
    if(frameColors == NULL) {
        perror("bruh");
        goto free;
    }

    torus_init(meshes[0], 100, 50, HALFWINWIDTH, HALFWINHEIGHT, 0);

    Mesh_draw(meshes[0], &camera, frameColors, 0, 0, 20);
    
    SDL_RenderPresent(renderer);

    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f milliseconds to init\n", cpu_time_used);



    printf("%f %f %f\n", sphere.velocity.x, sphere.velocity.y, sphere.velocity.z);








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

        

        SDL_UpdateTexture(texture, NULL, frameColors, SCREENWIDTH * sizeof(unsigned int));
        SDL_RenderCopy(renderer, texture, NULL, NULL);

        SDL_RenderPresent(renderer);


       
            Mesh_draw(meshes, &camera, frameColors, drawPrecision);

            Mesh_rotate(meshes, deltaRad, 'x');
            Mesh_rotate(meshes, deltaRad, 'y');
            Mesh_rotate(meshes, deltaRad, 'z');
        
  

        theta += 3 * deltaRad;
        Rad += deltaRad;

        end = clock();
        cpu_time_used = ((double) (end - start));

        if(cpu_time_used < 16) {
            usleep(1000 * (16 - cpu_time_used));
        }

        printf("%4.0f milliseconds per frame\n", cpu_time_used);



        /*
        pipiline:
        1. clear frame
        2. draw 
        3. present frame
        4. rotate 
        5. rotate 
        6. check for collisions
            





        */  
    }

    free:

    free_space(meshes, 1);

    free(frameColors);
    free(shaderProgram);

    SDLquit: 

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    quit:
    
    return 0;
}
