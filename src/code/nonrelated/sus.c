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
    int t;
    scanf("%d", &t);

    while(t--)
    {
        int n;
        scanf("%d", &n);

        if(n % 2 == 0)
        {
            printf("%d\n", n);
        } else
        {
            printf("%d\n", n - 1);
        }

    }

    return 0;
}