


#include "code/mathematika/mathcore.h"
#include "code/mathematika/pipeline.h"
#include "code/utilities/shared.h"
#include "code/utilities/glutilities.h"


int main(int argc, char** argv) 
{
    // SDL
    SDL_Window** window = NULL;
    SDL_Texture** texture = NULL;
    SDL_Renderer** renderer = NULL;

    if(SDL_init(&window, &renderer, &texture, "3D", SCREENWIDTH, SCREENHEIGHT) == 1) {
        goto quit;
    }

    clock_t start, end;
    double cpu_time_used;


    start = clock();

    Scene sceneInstance;

    initializeScene(&sceneInstance);

    

    

    
    torus_init(sceneInstance->meshes[0], 100, 50, HALFWINWIDTH, HALFWINHEIGHT, 0);

    Mesh_draw(sceneInstance->meshes[0], sceneInstance->camera, sceneInstance->frameColors, 0, 0, 20);
    
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

        memset(sceneInstance->frameColors, 0xFF, SCREENWIDTH * SCREENHEIGHT * sizeof(unsigned int));

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        

        SDL_UpdateTexture(texture, NULL, sceneInstance->frameColors, SCREENWIDTH * sizeof(unsigned int));
        SDL_RenderCopy(renderer, texture, NULL, NULL);

        SDL_RenderPresent(renderer);


       
            Mesh_draw(sceneInstance->meshes, sceneInstance->camera, sceneInstance->frameColors, sceneInstance->drawPrecision);

            Mesh_rotate(sceneInstance->meshes, deltaRad, 'x');
            Mesh_rotate(sceneInstance->meshes, deltaRad, 'y');
            Mesh_rotate(sceneInstance->meshes, deltaRad, 'z');
        
  

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

 
    free_space(sceneInstance);

    
    
    return 0;
}
