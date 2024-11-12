#include "graphics.h"


int main(int argc, char** argv) {

    SDL_Window* window = NULL;
    SDL_Texture* texture = NULL;
    SDL_Renderer* renderer = NULL;

    if(SDL_init(&window, &renderer, &texture, "3D", S_WIDTH, S_HEIGHT) == 1) {
        goto quit;
    }


    clock_t start, end;
    double cpu_time_used;


   






    // SDL_Delay(3000);


    start = clock();



    

    Manifold torus;
    Manifold sphere;

    unsigned int* frameColors = NULL;

    Vector3f torusNormals[POINTS];
    Vector3f sphereNormals[POINTS];


    if(get_space(&torus) == 1) {
        goto exit;
    }

    if(get_space(&sphere) == 1) {
        goto exit;
    }

    frameColors = malloc(S_WIDTH * S_HEIGHT * sizeof(unsigned int));
    if(frameColors == NULL) {
        perror("bruh");
        goto free;
    }
        



  
  
    torus_init(&torus, torusNormals, 2.5, 2);
    sphere_init(&sphere, sphereNormals, 6.5);


    Manifold_draw(renderer, &torus, torusNormals, frameColors, 0, 0, 20);
    Manifold_draw(renderer, &sphere, sphereNormals, frameColors, 0, 0, 20);



    SDL_RenderPresent(renderer);








    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f miliseconds to init\n", cpu_time_used);
    



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

        memset(frameColors, 0xFF, S_WIDTH * S_HEIGHT * sizeof(unsigned int));


        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);



        Manifold_draw(renderer, &sphere, sphereNormals, frameColors, -200, 0, drawPrecision);
        Manifold_draw(renderer, &torus, torusNormals, frameColors, 200, 130 * sin(Rad), drawPrecision);


        SDL_UpdateTexture(texture, NULL, frameColors, S_WIDTH * sizeof(unsigned int));
        


        SDL_RenderCopy(renderer, texture, NULL, NULL);

        SDL_RenderPresent(renderer);



        Manifold_rotate(&torus, torusNormals, deltaRad, 'x');
        Manifold_rotate(&torus, torusNormals, deltaRad, 'y');
        Manifold_rotate(&torus, torusNormals, deltaRad, 'z');


        Manifold_rotate(&sphere, sphereNormals, deltaRad, 'x');
        Manifold_rotate(&sphere, sphereNormals, deltaRad, 'y');
        Manifold_rotate(&sphere, sphereNormals, deltaRad, 'z');
      



        theta += 3 * deltaRad;
        Rad += deltaRad;




        end = clock();

        cpu_time_used = ((double) (end - start));

        if(cpu_time_used < 16) {
            usleep(1000 * (16 - cpu_time_used));
        }

        printf("%4.0f miliseconds per frame\n", cpu_time_used);  


    

        
    
    }


    free:


    free(torus.x);
    free(torus.y);
    free(torus.z);

    free(sphere.x);
    free(sphere.y);
    free(sphere.z);

    free(frameColors);

    exit: 

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    quit:
    
    return 0;
}
