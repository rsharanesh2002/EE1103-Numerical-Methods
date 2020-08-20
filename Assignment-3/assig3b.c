/******************************************************************************************************************************************
Purpose: To downsample the data in the file 'out1_test0.csv' by taking every 50th point and then
regenerating the intermediate points by interpolation using cubic splines and and comparing the root
mean square error between the actual data set and the interpolated points.
Authors:Sharansh R(EE19B116),Sakthi Harish D.T(EE19B054),TMVS Ganesh(EE19B124)
Date   : 4th September,2019
Input  : The data from the file 'out1_test0.csv'.
Output : The value of the root mean square error.
          "The RMS error is : 0.000463"
*******************************************************************************************************************************************/

#include <stdio.h>
#include <math.h>

void tridiagonalCubicSplineGen(int n, double h[n], double a[n - 1][n], double y[n + 1]);
void gaussEliminationLS(int m, int n, double a[m][n], double x[n - 1]); // Function to perform Gauss elimination
void cSCoeffCalc(int n, double h[n], double sig[n + 1], double y[n + 1], double a[n], double b[n],double c[n], double d[n]);
void rmserror(); // Function to find RMS error

void main()
{

    double y[51], h[50], result,inp; // y array to store the y-axis points, h array to store the successive interval widths
    double a[50], b[50], c[50], d[50]; // arrays to store the a,b,c,d values
    double sig[51]; // array to store Si's
    double sigTemp[49]; // array to store the Si's except S0 and Sn
    double tri[49][50]; // matrix to store the tridiagonal system of equations that will solve for Si's

    int x[51], i1 = 0, c1 = 1, count = 0; // x array to store the x-axis points

    FILE *in, *out; // File pointers
    in = fopen("out1_test0.csv", "r"); // Input file
    out = fopen("outdata.txt", "w"); // Output file


/*******************************************************************
To downsample the data from 250ps to 25ps by taking every 10th point
*******************************************************************/

    while (!feof(in))
    {
        fscanf(in, "%lf", &inp);
        count++;
        if ((count % 10) == 1)
        {
            x[i1] = c1;
            y[i1] = inp;
            i1++;
            c1 += 10;
        }
    }
    fclose(in);


    for (int i = 0; i < 50; i++) // computing the h matrix
        h[i] = x[i + 1] - x[i];


/************************
To initialize tri[n-1][n]
************************/

    for (int i = 0; i < 49; i++)
    {
        for (int j = 0; j < 50; j++)
        {
            tri[i][j] = 0;
        }
    }

    sig[0] = 0; // Initialize Si at 0
    sig[50] = 0; // Initialize Si at n


    tridiagonalCubicSplineGen(50, h, tri, y); // To generate the matrix
    gaussEliminationLS(49, 50, tri, sigTemp); // Perform Gauss Elimination


/********************
To compute Si values
********************/

    for (int i = 1; i < 50; i++)
    {

        sig[i] = sigTemp[i - 1];
    }


    cSCoeffCalc(50, h, sig, y, a, b, c, d); // calculate the values of ai's, bi's, ci's, and di's


    for (int i = 0; i < 50; i++)
    {
        for (int j = (i * 10); j < ((i + 1) * 10); j++)
        {
            result
                = a[i] * pow((j - x[i]), 3) + b[i] * pow((j - x[i]), 2) + c[i] * (j - x[i]) + d[i];

            fprintf(out, "%lf\n", result); // Writing the interpolated values to a file
        }
    }


    rmserror(); // calculate the RMS error

}

/*****************************
Function to calculte RMS error
******************************/

void rmserror()
{
    FILE *ptr1, *ptr2;
    double x1, x2, err = 0;
    ptr1 = fopen("result.txt", "r");
    ptr2 = fopen("out1_test0.csv", "r");
    while (fscanf(ptr1, "%lf", &x1) != EOF)
    {
        fscanf(ptr2, "%lf", &x2);
        err += ((x1 - x2) * (x1 - x2));//Calculating the RMS error
    }
    err /= 500;
    err = pow(err, 0.5);
    printf("The RMS error is : %lf", err); // Printing the RMS error
}

/************************************
Function to perform Gauss elimination
************************************/

void gaussEliminationLS(int m, int n, double a[m][n], double x[n - 1])
{
    double term;
    int i, j, k;
    for (i = 0; i < m - 1; i++)
    {

        // Begin Gauss Elimination
        for (k = i + 1; k < m; k++)
        {
            term = a[k][i] / a[i][i];
            for (j = 0; j < n; j++)
            {
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }
    }

    // Begin Back-substitution
    for (i = m - 1; i >= 0; i--)
    {
        x[i] = a[i][n - 1];
        for (j = i + 1; j < n - 1; j++)
        {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }
}

/**********************************************************
Cubic Spline coefficients calculator
Function that calculates the values of ai, bi, ci, and di's
***********************************************************/

void cSCoeffCalc(int n, double h[n], double sig[n + 1], double y[n + 1], double a[n], double b[n],
    double c[n], double d[n])
{
    int i;
    for (i = 0; i < n; i++)
    {
        d[i] = y[i];
        b[i] = sig[i] / 2.0;
        a[i] = (sig[i + 1] - sig[i]) / (h[i] * 6.0);
        c[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * sig[i] + sig[i + 1]) / 6.0;
    }
}

/****************************************************
Function to generate the tridiagonal aurgmented matrix
for cubic spline
****************************************************/

void tridiagonalCubicSplineGen(int n, double h[n], double a[n - 1][n], double y[n + 1])
{
    int i;
    for (i = 0; i < n - 1; i++)
    {
        a[i][i] = 2 * (h[i] + h[i + 1]);
    }
    for (i = 0; i < n - 2; i++)
    {
        a[i][i + 1] = h[i + 1];
        a[i + 1][i] = h[i + 1];
    }
    for (i = 1; i < n; i++)
    {
        a[i - 1][n - 1]
            = (y[i + 1] - y[i]) * 6.0 / (double)h[i] - (y[i] - y[i - 1]) * (6.0 / (double)h[i - 1]);
    }
}

