#include "sdl_lib.h"


int main(int argc, char** argv) {

    SDL_Window* window = NULL;
    SDL_Renderer* renderer = SDL_INIT(&window, "3D", 1980, 1080);
    if(renderer == NULL) {
        perror("renderer null bruh");
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }


    clock_t start, end;
    double cpu_time_used;

    /*
    Vector2 vec1 = {50, 50};
    Vector2 vec2 = {200, 0};
    Vector2 vecC = {100, 100};

    fillTriangle(renderer, vec1, vec2, vecC, 50, 128, 128, 128);

    SDL_RenderPresent(renderer); */
   





    
    

    start = clock();

    Manifold torus = get_space(torus); 
    // Manifold torus1 = get_space(torus1);
  
    torus = torus_init(torus, 2.5, 2);
    // torus1 = torus_init(torus, 3, 2);
 

    //torus = M_rotate(torus, PI / 6, 'x');
    //torus = M_rotate(torus, PI / 4, 'y');

    torus_draw(renderer, torus, 0, 0);
    // torus_draw(renderer, torus1, 400, 0);


    
    //sphere_init(sphere, 4);
    //draw(renderer, sphere, 0, 0); 



    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f miliseconds to init\n\n", cpu_time_used);
    

    double deltaTime = 0;
    float deltaRad = PI / 1440;
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
    

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);
        //SDL_RenderPresent(renderer);


        torus_draw(renderer, torus, 0, 100 * sin(theta));
        // torus_draw(renderer, torus1, 400, 120 * sin(theta));
       

        SDL_RenderPresent(renderer);


        torus = M_rotate(torus, deltaRad, 'x');
        torus = M_rotate(torus, deltaRad, 'y');

        // torus1 = M_rotate(torus1, deltaRad, 'x');
        // torus1 = M_rotate(torus1, deltaRad, 'y');
               

        theta += 3 * deltaRad;
        Rad += deltaRad;

        //usleep(50000); 

        end = clock();

        cpu_time_used = ((double) (end - start));
        printf("%4.0f FPS\n", 1000 / cpu_time_used); 

    

    
        

    }

    free(torus.x);
    free(torus.y);
    free(torus.z);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}