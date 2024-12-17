// #include "utils.h"


#include <stdio.h>
#include <limits.h> 
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define BUFFLEN 512
#define MAXNUMCOUNT 256



int main()
{
    // setvbuf(stdout, NULL, _IOFBF, 1073741824);

    // 4294967296
    // 65536
    // 16384
    // 1048576
    // 4194304
    // 67108864
    // 268435456
    // 1073741824

/*

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
    */
/*
    int n;

    printf("Ievadiet balonu skaitu: ");
    scanf("%d", &n);

    int balloons[n];

    printf("Ievadiet balonu augstumus:\n");
    for(int i = 0; i < n; i++) 
    {
        scanf("%d", &balloons[i]);
    }

    int arrows[1000] = {0}; 
    int arrow_count = 0;    

    for(int i = 0; i < n; i++) 
    {
        int height = balloons[i];
        int shot = 0;

        for(int j = 0; j < arrow_count; j++) 
        {
            if(arrows[j] == height) 
            {
                arrows[j]--; 
                shot = 1;    
                break;
            }
        }

        if(!shot) 
        {
            arrows[arrow_count++] = height - 1;
        }
    }

    printf("Nepieciešamais bultu skaits: %d\n", arrow_count);

*/


    Node* list = NULL;


    return 0;
}
