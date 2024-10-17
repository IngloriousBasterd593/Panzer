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

    /*
    Vector2 vec1 = {50, 50};
    Vector2 vec2 = {200, 0};
    Vector2 vecC = {100, 100};

    fillTriangle(renderer, vec1, vec2, vecC, 50, 128, 128, 128);

    SDL_RenderPresent(renderer); */

    
    //Vector3 vec1 = {1, 1, 1};
    //Vector3 vec2 = {3, 2, 8};
    //vec1 = unit(vec1);
    //vec2 = unit(vec2);

    // Vector3 vector = normal(vec1, vec2);

    // printf("%f\n", vector.x);







    // usleep(1000000); 












    start = clock();

    Manifold torus;
    Manifold sphere;

    SurfaceNormals sphereNormals;



    // SDL_Delay(5000);

    get_space(&torus); 
    get_space(&sphere); 
  
    sphere_init(&sphere, &sphereNormals, 6);
    // torus_init(&torus, 2.5, 2);

    // torus_draw(renderer, &torus, 2.5, 2, 0, 0, 10);
    sphere_draw(renderer, &sphere, &sphereNormals, 0, 0, 10);


    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f miliseconds to init\n", cpu_time_used);
    

    double deltaTime = 0;
    float deltaRad = twopiOverPixels / 12;
    float Rad = 0;
    int radius = 150;
    float theta = 0;
    int precision = 25;
    // long iterations = 0;


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


        // torus_draw(renderer, &torus, 2.5, 2, 2, 100 * sin(theta), precision);
        sphere_draw(renderer, &sphere, &sphereNormals, 0, 0, precision);
       

        SDL_RenderPresent(renderer);

        M_rotate(&sphere, deltaRad, 'y');
        M_rotate(&sphere, deltaRad, 'x');

        V_rotate(&sphereNormals, deltaRad, 'y');
        V_rotate(&sphereNormals, deltaRad, 'x');

        printf("%f\n", sphereNormals.u[0].x);


        // M_rotate(&torus, deltaRad, 'x');
        // M_rotate(&torus, deltaRad, 'y');
               

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