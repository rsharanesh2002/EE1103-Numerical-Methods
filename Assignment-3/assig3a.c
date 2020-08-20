/********************************************************************************************************************************************
Purpose:Take the values of x and y=tan(x), for x = 1, 1.1, 1.2 and 1.3 and use these values to interpolate and find the value for y at x=1.15  and x=1.25 using linear interpolation,lagranges second degree polynomial and cubic splines.
Authors:Sharansh R(EE19B116),Sakthi Harish D.T(EE19B054),TMVS Ganesh(EE19B124)
Date : 4th September,2019
Output:The results of interpolation using various interpolating methods
 "The value of tan(1.15) using Lagrange multiplier is : 2.229543
         Error in tan(1.15) is : 0.004954
 The value of tan(1.35) using Lagrange multiplier is : 4.345078
         Error in tan(1.35) is : 0.110143
 The value of tan(1.15) using newton's method is: 2.229543
         Error in tan(1.15) is : 0.004954
 The value of tan(1.35) using newton's method is: 4.345076
         Error in tan(1.35) is : 0.110146
 The value of tan(1.15) using Cubic spline is : 2.219990
         Error in tan(1.15) is : 0.014507
 The value of tan(1.35) using Cubic spline is : 4.156694
         Error in tan(1.35) is : 0.298527"
**********************************************************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void lagrange(float*, float*); // Function for interpolating using lagrange polynomial method
void newton(float*); // Function for interpolating using Newton method
void cspline(float*, float*); // Function for interpolating using Cubic Spline method

float x[4] = { 1.0, 1.1, 1.2, 1.3 }; // Input x values
float y[4] = { tan(1.0), tan(1.1), tan(1.2), tan(1.3) }; // Input function values
float z[2] = { 1.15, 1.35 }; //
float value[2] = { tan(1.15), tan(1.35) }; // Actual function values
float error[2][3]; // Array to store the error

int i, j, k, n = 4;

int main()
{

    lagrange(x, y); // interpolating using Lagrange Polynomial method
    newton(x); // interpolating using Newton method
    cspline(x, y); //Function for interpolating using Cubic Spline method

    return 0;

}


void lagrange(float x[], float y[]) // Function for interpolating using Lagrange Polynomial method
{
    float term[5], out[5] = { 0.0 };
    for (k = 0; k < 2; k++)
    {
        for (int i = 0; i < 4; i++)
        {
            term[k] = y[i];

            for (int j = 0; j < 4; j++)
            {
                if (j != i)
                {
                    term[k] = term[k] * (z[k] - x[j]) / (x[i] - x[j]);
                }
            }

            out[k] += term[k];
        }
        error[k][0] = (value[k] - out[k]); // calculating error
        printf("\n The value of tan(%.2lf) using Lagrange multiplier is : %lf", z[k], out[k]);
        printf("\n \t Error in tan(%.2lf) is : %lf ", z[k], error[k][0]);
    }
}


void newton(float x[]) // Function for interpolating using Newton method
{
    float nt[10][10], sum, pro;

    for (i = 0; i < n; i++) // initialize the newton array
        nt[i][0] = tan(x[i]);


    for (i = 1; i < n; i++)
    {
        for (j = 0; j < n - i; j++)
        {
            nt[j][i] = (nt[j][i - 1] - nt[j + 1][i - 1]) / (x[j] - x[i + j]);
        }
    }
    for (i = 0; i < 2; i++)
    {
        sum = nt[0][0];

        for (int k = 1; k < n; k++)
        {
            pro = 1;
            for (int j = 0; j < k; j++)
            {
                pro = pro * (z[i] - x[j]);
            }
            sum = sum + (pro * nt[0][k]);
        }
        error[i][1] = value[i] - sum; // calculating error
        printf("\n The value of tan(%.2lf) using newton's method is: %lf ", z[i], sum);
        printf("\n \t Error in tan(%.2lf) is : %lf ", z[i], error[i][1]);
    }
}


void cspline(float x[], float y[]) // Function for interpolating using Cubic Spline method
{
    float h[3], b[3], v[3], u[3], w[4], out[2];
    for (k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            h[i] = x[i + 1] - x[i];
            b[i] = (y[i + 1] - y[i]) * (1.0 / h[i]);
        }

        for (int i = 1; i < 3; i++)
        {
            v[i] = 2 * (h[i - 1] + h[i]);
            u[i] = 6 * (b[i] - b[i - 1]);
        }

        w[0] = w[3] = 0;

        // gaussian elimination for tridiagonal systems

        w[2] = u[2] / v[2];
        w[1] = (u[1] - h[1] * w[2]) / v[1];

        for (int i = 0; i < 4; i++)
        {
            if (z[k] == x[i])
            {
                out[k] = y[i];
                break;
            }
            else if (z[k] < x[i])
            {

                out[k] = w[i] * pow((z[k] - x[i - 1]), 3) / (6 * h[i - 1])
                    + w[i - 1] * pow((x[i] - z[k]), 3) / (6 * h[i - 1])
                    + ((y[i] / h[i - 1]) - ((w[i] * h[i - 1]) / 6)) * (z[k] - x[i - 1])
                    + ((y[i - 1] / h[i - 1]) - ((w[i - 1] * h[i - 1]) / 6)) * (x[i] - z[k]);
                break;
            }
            else if (z[k] > x[3])
            {
                out[k] = w[3] * pow((z[k] - x[2]), 3) / (6 * h[2])
                    + w[2] * pow((x[3] - z[k]), 3) / (6 * h[2])
                    + ((y[3] / h[2]) - ((w[3] * h[2]) / 6)) * (z[k] - x[2])
                    + ((y[2] / h[2]) - ((w[2] * h[2]) / 6)) * (x[3] - z[k]);
                break;
            }
        }
        error[k][2] = (value[k] - out[k]); // calculating error
        printf("\n The value of tan(%.2lf) using Cubic spline is : %lf", z[k], out[k]);
        printf("\n \t Error in tan(%.2lf) is : %lf ", z[k], error[k][2]);
    }
}

