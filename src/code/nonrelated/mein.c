#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int parse(const char *str, int *arr, int max_size) 
{
    int count = 0;  
    const char *ptr = str;

    while(*ptr && count < max_size) 
    {
        while(*ptr == ' ') 
        {
            ptr++;
        }
     
        if(*ptr != '\0') 
        {
            arr[count++] = strtol(ptr, (char **) &ptr, 10);
        }
    }

    return count;  
}

int main() 
{
    int t, n, k;
    char buff[512];
    memset(buff, '\0', 512);
    scanf("%d", &t);

    int arr[2];
    memset(arr, 0, 8);

    while(t >= 0)
    {
        fgets(buff, 512, stdin);
        buff[strlen(buff) - 1]= '\0';

        parse(buff, arr, 2);

        int a[arr[0]];
        memset(a, 0, sizeof(int) * arr[0]);


    







        t--;
    }

    return 0;
}
