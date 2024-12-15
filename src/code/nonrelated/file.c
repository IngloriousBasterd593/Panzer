#include "utils.h"

int main()
{
    setvbuf(stdout, NULL, _IOFBF, 1073741824);

    // 4294967296
    // 65536
    // 16384
    // 1048576
    // 4194304
    // 67108864
    // 268435456
    // 1073741824

    int n;
    scanf("%d", &n);

    int lines = n;
    int stars = 1;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < lines; j++)
        {
            printf(" ");
        }

        for(int k = 0; k < stars; k++)
        {
            printf("*");
        }

        printf("\n");

        stars += 2;
        lines--;
    } 

   
    fflush(stdout);

    return 0;
}
