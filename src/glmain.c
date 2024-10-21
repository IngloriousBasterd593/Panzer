#include "C:\glfw-3.4.bin.WIN64\include\GLFW\glfw3.h"
#include <unistd.h>
#include <time.h>
#include <stdio.h>

int main(void)


{

    clock_t start, end;
    double cpu_time_used;
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {

        start = clock();





        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();






    usleep(50000);



    end = clock();
    cpu_time_used = ((double) (end - start));
    printf("%f FPS\n", 1000 / cpu_time_used);

    }

    glfwTerminate();
    return 0;
}