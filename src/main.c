#include "sdl_lib.h"



int main(int argc, char** argv) {

    SDL_Window* window = NULL;
    SDL_Renderer* renderer = SDL_INIT(&window, "3D", S_WIDTH, S_HEIGHT);
    if(renderer == NULL) {
        perror("renderer null bruh");
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }


    clock_t start, end;
    double cpu_time_used;

    start = clock();





    Manifold torus;
    Manifold sphere;

    Vector3 sphereNormals[POINTS] = {0};

    get_space(&torus); 
    get_space(&sphere); 
  
    sphere_init(&sphere, sphereNormals, 6);



    sphere_draw(renderer, &sphere, sphereNormals, 0, 0, 10);

    // printf("why slow");

    SDL_RenderPresent(renderer);

    // sphere_draw(renderer, &sphere, sphereNormals, 0, 0, 10);









    





    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f miliseconds to init\n", cpu_time_used);
    



    double deltaTime = 0;
    float deltaRad = PI / 120;
    float Rad = 0;
    int radius = 150;
    float theta = 0;
    int trianglePrecision = 10;

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



        sphere_draw(renderer, &sphere, sphereNormals, 0, 0, trianglePrecision);

        SDL_RenderPresent(renderer);
    

       
        Manifold_rotate(&sphere, sphereNormals, deltaRad, 'x');




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

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}
