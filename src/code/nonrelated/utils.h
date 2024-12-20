#ifndef UTILS
#define UTILS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define MAXBUFFLEN 512
#define MAXNUMCOUNT 256
#define INTLIMITDIGITS 10


typedef struct {
    char* buff;
    int* parsedToInt;
    int len;
} Input;


typedef struct {
    void* data;
    Node* next;
} Node;


static enum Datatypes {
    int,
    long, 
    float,
    double,
    char,
    char*,
    char**,
    int*,
    float*
}

int bubbleSort(int* arr, int len)
{
    if(arr == NULL)
    {
        return EXIT_FAILURE;
    }

    int temp;
    for(int j = 0; j < len; j++)
    {
        for(int i = 0; i < len; i++)
        {
            if(i + 1 >= len) 
            {
                break;
            }
            if(arr[i] > arr[i + 1])
            {
                temp = arr[i + 1];
                arr[i + 1] = arr[i];
                arr[i] = temp;
            }
        }
    }

    return EXIT_SUCCESS;
}


int IntArrayMax(int* arr, int len)
{
    int max = arr[0];
    
    for(int i = 1; i < len; i++) 
    {
        if(arr[i] > max)
        {
            max = arr[i];
        }
    }

    return max;
}


int intArraySum(int* arr, int len)
{
    int sum = 0;

    for(int i = 0; i < len; i++)
    {
        sum += arr[i];
    }

    return sum;
}

int isCharInRange(char c, int lower, int upper)
{
    if(c >= lower && c <= upper)
    {
        return true;
    } else
    {
        return false;
    }
}


int isIntInRange(int num, int lower, int upper)
{
    if(num >= lower && num <= upper)
    {
        return true;
    } else
    {
        return false;
    }
}


int parseStringToIntArray(char* buff, int* arr, char delimeter, int base)
{
    if(buff == NULL || arr == NULL || !isIntInRange(base, 2, 10))
    {
        return EXIT_FAILURE;
    }

    int lenght = 0;
    int len = strlen(buff);
    char num[INTLIMITDIGITS + 1];
    int numOfInts = 1;
    int number = 0;
    int numCount = 0;
    int numSize = 0;
    int power = 0;

    memset(num, '\0', INTLIMITDIGITS + 1);

    

    for(int i = 0; i < len; i++)
    {
        if(lenght >= MAXNUMCOUNT)
        {
            return EXIT_FAILURE;
        }

        if(i + 1 >= len)
        {
            lenght++;
            for(int j = numSize - 1; j >= 0; j--)
            {
                number += num[j] * ((int)pow(base, power));
                power++;
            }

            arr[numCount] = number;

            return lenght;
        }

        if(buff[i] != delimeter)
        {
            if(isCharInRange(buff[i], 48, 57))
            {
                num[numSize] = buff[i] - 48;
            } else if(isCharInRange(buff[i], 65, 70))
            {
                num[numSize] = buff[i] - 55;
            } else
            {
                fprintf(stderr, "Enter a valid symbol!");
                return EXIT_FAILURE;
            }

            
            numSize++;
        }

        if(buff[i] == delimeter)
        {
            for(int j = numSize - 1; j >= 0; j--)
            {
                number += num[j] * ((int)pow(base, power));
                power++;
            }


            memset(num, '\0', INTLIMITDIGITS + 1);
            arr[numCount] = number;
            power = 0;
            number = 0;
            numCount++;
            numSize = 0;
            lenght++;
        }
    }

    return EXIT_FAILURE;
}


int greedyMaxSubarraySum(int* arr, int len) 
{
    int maxSum;
    int depth = 0;
    int currArr[len];
    int sums[MAXNUMCOUNT];
    int counter = 0;

    memset(sums, 0, sizeof(int) * MAXNUMCOUNT);
    memset(currArr, 0, sizeof(int) * len);

    for(int i = 0; i < len; i++)
    {
        for(int j = 0; j < len; j++)
        {
            for(int k = depth; k < len; k++)
            {
                currArr[k - depth] = arr[k];
                sums[counter] = intArraySum(currArr, len);
                counter++;

                printf("%d ", intArraySum(currArr, len));
            }

            memset(currArr, 0, sizeof(int) * len);
            depth++;
        }
    }

    return IntArrayMax(sums, counter);
}


int maxSubarraySum(int* arr, int len) 
{
    int maxSum = INT_MIN; 
    int currSum = 0;   

    for(int i = 0; i < len; i++) 
    {
        currSum += arr[i];
        if(currSum > maxSum) 
        {
            maxSum = currSum; 
        }
        if(currSum < 0) 
        {
            currSum = 0; 
        }
    }

    return maxSum;
}


Input** initUserInputBuffers(int count)
{
    Input** inputs = (Input**) malloc(count * sizeof(Input*));
    if(inputs == NULL)
    {
        return NULL;
    }

    for(int i = 0; i < count; i++)
    {
        inputs[i] = (Input*) malloc(count * sizeof(Input));
        if(inputs[i] == NULL)
        {
            while(i--) 
            {
                free(inputs[i]);
            }

            return NULL;
        } 
    }

    for(int i = 0; i < count; i++)
    {
        inputs[i]->buff = malloc(sizeof(char) * MAXBUFFLEN);
        if(inputs[i]->buff == NULL)
        {
            while(i--) 
            {
                free(inputs[i]->buff);
            }

            return NULL;
        }

        inputs[i]->parsedToInt = malloc(sizeof(char) * MAXNUMCOUNT);
        if(inputs[i]->buff == NULL)
        {
            while(i--) 
            {
                free(inputs[i]->parsedToInt);
            }

            return NULL;
        }
    }

    return inputs;
}


void getUserInput(Input* inputInstance)
{
    int base;
    memset(inputInstance->buff, '\0', MAXBUFFLEN);
    memset(inputInstance->parsedToInt, 0, sizeof(int) * MAXNUMCOUNT);

    printf("Enter a base: ");
    scanf("%d", &base);
    getchar();

    printf("Enter a masivs: ");
    fgets(inputInstance->buff, MAXBUFFLEN, stdin);
    inputInstance->buff[MAXBUFFLEN - 2] = '\0';

    inputInstance->len = parseStringToIntArray(inputInstance->buff, inputInstance->parsedToInt, ' ', base);

    return; 
}


void freeInputObject(Input* input)
{
    free(input->buff);
    free(input->parsedToInt);
    free(input);
}


int isPrime(int n)
{
    if(n <= 2)
    {
        return true;
    }

    for(int i = 2; i < n; i++)
    {
        if(n % i == 0)
        {
            return false;
        }
    }

    return true;
}


int stalinSort(int* arr, int* res, int len)
{
    int j = 0;
    for(int i = 0; i < len - 1; i++)
    {
    
        if(arr[i] < arr[i + 1])
        {
            res[j] = arr[i];
            j++;
            if(i >= len - 1)
            {
                printf("sniga");
                res[j] = arr[i + 1];
                break;
            }
        } else if(arr[i] > arr[i + 1]) 
        { 
            continue;
        } else 
        {
            res[j] = arr[i];
            j++;
        }
    }

    return EXIT_SUCCESS;
}


int initLinkedList(Node* list, int len)
{
    list = (Node*) malloc(len * sizeof(Node));
    if(list == NULL)
    {
        return EXIT_FAILURE;
    }

    for(int i = 0; i < len - 1; i++)
    {
        list[i].next = &list[i + 1];
    }

    list[len - 1].next = NULL;

    return EXIT_SUCCESS;
}


int addNewNode(Node* list, int len)
{


    return EXIT_SUCCESS;
}


int deleteNode(Node* list, int len)
{


    return EXIT_SUCCESS;
}

int appendToNode(Node* list, void* data, enum Datatypes type)
{
    switch(type) 
    {
        case: int*

        break;

        case: long*

        break;

        case: float*

        break;

        case: double*

        break;

        case: char*

        break;

        case: char**,

        break;

        case: int**

        break;

        default: 
        break;
    }


    return EXIT_SUCCESS;
}

static enum Datatypes {
    int,
    long, 
    float,
    double,
    char,
    char*,
    char**,
    int*,
    float*
}


#endif
