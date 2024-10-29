#include "graphics.h"



int main(int argc, char** argv) {

    SDL_Window* window = NULL;
    SDL_Texture* texture = NULL;
    SDL_Renderer* renderer = NULL;

    SDL_init(&window, &renderer, &texture, "3D", S_WIDTH, S_HEIGHT);


    clock_t start, end;
    double cpu_time_used;





    start = clock();



    

    Manifold torus;
    Manifold sphere;

    unsigned int* frameColors = NULL;
    unsigned int* IDK_butthecodeisbroken = NULL;

    Vector3 torusNormals[POINTS] = {0};


    if(get_space(&torus, &frameColors) == 1) {
        goto exit;
    }

    if(get_space(&sphere, &IDK_butthecodeisbroken) == 1) {
        goto exit;
    }

        



  
  
    torus_init(&torus, torusNormals, 2.5, 2);


    Manifold_draw(renderer, &torus, torusNormals, frameColors, 0, 0, 20);



    SDL_RenderPresent(renderer);








    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f miliseconds to init\n", cpu_time_used);
    



    double deltaTime = 0;
    float deltaRad = PI / 120;
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

        

        


        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);




        Manifold_draw(renderer, &torus, torusNormals, frameColors, 0, 120 * sin(2 * Rad), drawPrecision);

        SDL_UpdateTexture(texture, NULL, frameColors, S_WIDTH * S_HEIGHT * sizeof(unsigned int));


        SDL_RenderCopy(renderer, texture, NULL, NULL);

        SDL_RenderPresent(renderer);



        Manifold_rotate(&torus, torusNormals, deltaRad, 'x');
        Manifold_rotate(&torus, torusNormals, deltaRad, 'y');
        Manifold_rotate(&torus, torusNormals, deltaRad, 'z');

      



        theta += 3 * deltaRad;
        Rad += deltaRad;




        end = clock();


        cpu_time_used = ((double) (end - start));
        printf("%4.0f FPS\n", 1000 / cpu_time_used);  


    

        
    
    }


    free(torus.x);
    free(torus.y);
    free(torus.z);

    free(sphere.x);
    free(sphere.y);
    free(sphere.z);

    free(frameColors);
    free(IDK_butthecodeisbroken);

    exit: 

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}
