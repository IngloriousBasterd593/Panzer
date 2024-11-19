#ifndef GLUTILITIES
#define GLUTILITIES

#include "common.h"

#include <glad/glad.h>  
#include <GLFW/glfw3.h>



void copyShaderSource(const char* source, char* destination) 
{
    if(source == NULL || destination == NULL) 
    {
        fprintf(stderr, "source or destination dont exist\n");
    }

    FILE* shaderFile = fopen(source, "r");
    if(shaderFile == NULL) 
    {
        fprintf(stderr, "failed to open file\n");
        return;
    }

    int i = 0; 
    char ch;

    while((ch = fgetc(shaderFile)) != EOF) 
    {
        if(i > BUFFLEN - 2) 
        {
            destination[i] = '\0';
            break;
        }

        destination[i] = ch;
        i++;
    }

    fclose(shaderFile);

    return;
}

#endif