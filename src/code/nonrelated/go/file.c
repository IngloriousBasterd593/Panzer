#include "utils.h"

int main()
{
    Input** inputs = initUserInputBuffers(10);

    getUserInput(inputs[0]);


    printf("%d\n", kempelaUzdevums(inputs[0]->parsedToInt, inputs[0]->len));

























    for(int i = 0; i < 10; i++)
    {
        freeInputObject(inputs[i]);
    }

    free(inputs);
   
    return EXIT_SUCCESS;
}