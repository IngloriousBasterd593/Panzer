#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* int maxOfArray(int arr[]) 
{   
    int i = 0, lenght = 0;
    while((i * sizeof(int)) < sizeof(arr)) {
        lenght++;
        i++;
    }
    int max = arr[0];
    
    for(i = 0; i < lenght; i++) {
        if((i + 1) < lenght) {
            if(max < arr[i + 1]) {
                max = arr[i + 1];
            }
        }
    }
    
    return max;
} */

int main(int argc, char** argv) 
{
    int t;
    int arr[] = {5, 4, 5};
    int i = 0, lenght = 0;
    while((i * sizeof(int)) < sizeof(arr)) {
        lenght++;
        i++;
    }
    int max = arr[0];
    
    for(i = 0; i < lenght; i++) {
        if((i + 1) < lenght) {
            if(max < arr[i + 1]) {
                max = arr[i + 1];
            }
        }
    }
    
    int OddOrEven;
    for(i = 0; i < lenght; i++) {
        if(arr[i] == max) {
            break;
        }
    }

    if(i % 2 == 0) {
        OddOrEven = 1;
    } else {
        OddOrEven = 0;
    }

    printf("%d", OddOrEven);

    int sum;

    if(OddOrEven == 0 && lenght % 2 != 0) {
        sum = max + (lenght / 2) + 1;
    }
    if(OddOrEven == 1 && lenght % 2 != 0) {
        sum = max + (lenght / 2);
    }
    if(OddOrEven == 0 && lenght % 2 == 0) {
        sum = max + (lenght / 2);
    }
    if(OddOrEven == 1 && lenght % 2 != 0) {
        sum = max + (lenght / 2);
    }

    printf("%d", sum);

    return 0;
}