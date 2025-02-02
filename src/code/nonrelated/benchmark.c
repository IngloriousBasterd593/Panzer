
#include <stdint.h>
#include <time.h>
#include <stdio.h>
#define DX 100
#define FD 10
#define STEP (2 * 3.14159) / DX


typedef struct {
    uint32_t a, b;
    uint32_t c, d;
} Matrix;


static Matrix multiply(Matrix x, Matrix y)
{
    Matrix result;


    result.a = x.a * y.a + x.b * y.c;
    result.b = x.a * y.b + x.b * y.d;
    result.c = x.c * y.a + x.d * y.c;
    result.d = x.c * y.b + x.d * y.d;


    return result;
}


static Matrix matrix_power(Matrix base, uint16_t exponent)
{
    Matrix result = {1, 0, 0, 1};


    while(exponent > 0)
    {
        if(exponent & 1)
        {
            result = multiply(result, base);
        }


        base = multiply(base, base);
        exponent >>= 1;
    }


    return result;
}


static uint32_t fibonacci(uint16_t n)
{
    if(n == 0)
    {
        return 0;
    }
    else if(n == 1)
    {
        return 1;
    }


    Matrix base = {1, 1, 1, 0};
    Matrix result = matrix_power(base, n - 1);
   
    return result.a;
}


static void bubble_sort(uint16_t *arr, uint8_t n)
{
    uint32_t pass, i, j;
   
    for(pass = 0; pass < n; pass++)
    {
        for(i = 0; i < n - 1; i++)
        {
            for(j = 0; j < n - i - 1; j++)
            {
                if(arr[j] > arr[j+1])
                {
                    uint16_t temp = arr[j];
                    arr[j] = arr[j+1];
                    arr[j + 1] = temp;
                }
            }
        }
    }
}


int main()
{
    uint32_t n = 100000000;
    // 100 biti
    uint32_t result = fibonacci(n);


   




    uint16_t arr[50] = {
    148, 67, 253, 11, 36, 245, 0, 255, 127, 64,
    89, 50, 15, 99, 234, 128, 1, 200, 33, 90,
    55, 77, 123, 222, 187, 11, 67, 241, 2, 69,
    211, 45, 95, 13, 197, 38, 112, 76, 45, 157,
    3, 208, 46, 57, 69, 1, 250, 222, 99, 74
    };


    // 126 biti
    bubble_sort(arr, 50);








    return 0;
}


















static uint32_t factorial(uint32_t n)
{
    uint32_t result = 1;




    for(uint32_t i = 1; i <= n; i++)
    {
        result *= i;
    }




    return result;
}




static float power(float base, uint32_t exponent)
{
    float result = 1.0;




    for (uint32_t i = 0; i < exponent; i++)
    {
        result *= base;
    }




    return result;
}




static float sint(float x)
{
    float sum = 0.0;




    for (uint32_t n = 0; n < FD; n++)
    {
        float sign = (n % 2 == 0) ? 1.0 : -1.0;
        float numerator = power(x, 2 * n + 1);
        float denom = factorial(2 * n + 1);
        sum += sign * (numerator / denom);
    }




    return sum;
}




// 56 biti
static float cost(float x)
{
    float sum = 0.0;




    for (uint32_t n = 0; n < FD; n++)
    {
        float sign = (n % 2 == 0) ? 1.0 : -1.0;
        float numerator = power(x, 2 * n);
        float denom = factorial(2 * n);
        sum += sign * (numerator / denom);
    }




    return sum;
}




static void benchmark_torus(uint32_t R, uint32_t r)
{
    for(int j = 0; j < DX; j++)
    {
        for(int k = 0; k < DX; k++)
        {
            if(j * STEP > 2 * 3.14159)
            {
                break;
            }


            float garbage = (R + r * cost(STEP * k)) * sint(STEP * j);
            garbage = (R + r * cost(STEP * k)) * cost(STEP * j);
            garbage = r * sint(STEP * k);
        }
    }
}
     












