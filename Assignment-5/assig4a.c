/**********************************************************************************************************************************************
Purpose: To compute the LU Decomposition and also compute the solution using the above LU Decomposition for the cubic spline interpolation.
Authors: Sharansh R(EE19B116),Vibhhu Sharma(EE19B128),Adithyan CG(EE19B003)
Date   : 18th September,2019
Input  :                                                      PART-A
         The number of equation and the corresponding prefactors of variables and the coefficient
         matrix.

                                                              PART-B
         The number of data-points and the corresponding data-points and the particuar vale of x to
         find the interpolated value.


Output :                                                      PART-A
         The solution of the given input equations.

         "!!!!!  PART-A !!!!! (To compute the LU Decomposition and also compute the solution using
Gaussian Elimintion)

          Enter the number of equations : 2

          Enter the coefficients of variables seperated by space (in order) :
                1 2 1 1

          Enter the coeffient matrix seperated by space (in order) : 4 2

          The Matrix obtained by Multiplication of Lower and Upper triangular is same as the
coefficient matrix!!!

          The Solution of the given input is given below :

                x1 = 0.000000
                x2 = 2.000000 "



                                                              PART-B
         The interpolated value of the function at a particular value specified by the user.
         "!!!!!  PART-B !!!!! (Using the above method of LU Deomposition in Cubic Spline
Interpolation)

          Enter the number of points : 4

	  Enter the x value and the y value of the function for the number of points mentioned above(x  and y format): 
          1 2 
	  2 3
	  3 4
	  5 9
	  Enter the value of x at which the value of y is needed : 4
	  UPPER TRIANGULAR MATRIX
	  4.00 1.00 
	  0.00 5.75 
	  LOWER TRIANGULAR MATRIX
	  1.00 0.00 
	  0.25 1.00 
	  value of the function at x = 4.00 is y =  6.1087 
 "

***********************************************************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/************************************************************ !!!PART-A!!!*******************************************************************/

void parta()
{

    double a[5][5]; // coeficient of equations matrix

    double a1[5][5]; // to verify the matrix multiplication of lower and upper triangular matrix

    double u[5][5], l[5][5]; // to find the lower and upper matrix

    double b[5]; // coeficient matrix

    double x[5], y[5]; // to store the intermediate values

    double sum = 0.0, s = 0, diff = 0.0;

    int n, i, j, p, q, r, k; // iteration variables

    printf("\nEnter the number of equations : "); /// input the number of equations from user
    scanf("%d", &n);

    printf("\nEnter the coefficients of variables seperated by space (in order) :\n \t"); /// input the coefients of equations from user
                                                                                       

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            scanf("%lf", &a[i][j]);

    printf("\nEnter the coeffient matrix seperated by space (in order) : "); /// input the coefficient matrix
                                                                            
    for (i = 0; i < n; i++)
        scanf("%lf", &b[i]);

    /// Dolittle's algorithm for LU Decomposition
    for (i = 0; i < n; i++)
    {
        for (k = i; k < n; k++) /// computing upper triangular matrix
        {
            sum = 0.0;

            for (j = 0; j < i; j++)
                sum += l[i][j] * u[j][k];

            u[i][k] = a[i][k] - sum;
        }

        for (k = i; k < n; k++) /// computing lower triangular matrix
        {
            if (i == k)
                l[i][i] = 1.0;

            else
            {
                sum = 0.0;

                for (j = 0; j < i; j++)
                    sum += l[k][j] * u[j][i];

                l[k][i] = (a[k][i] - sum) / u[i][i];
            }
        }
    }

    /// computing the matrix multiplication of lower and upper triangular matrices
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                s = s + l[i][k] * u[k][j];
            }

            a1[i][j] = s;
            s = 0;
        }
    }

    /// computing the differences between the original and multiplied matrices
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            diff += (a[i][j] - a1[i][j]);

    if (diff == 0.0)
        printf("\nThe Matrix obtained by Multiplication of Lower and Upper triangular is same as "
               "the coefficient matrix!!! \n");


    // solve Ly=b by forward substitution
    y[0] = b[0] / l[0][0];

    for (i = 1; i < n; i++)
    {
        y[i] = b[i];

        for (j = 0; j < i; j++)
        {
            y[i] -= (l[i][j] * y[j]);
        }

        y[i] /= l[i][i];
    }


    // solve Ux=y by backward substitution
    x[n - 1] = y[n - 1] / u[n - 1][n - 1];

    for (i = n - 2; i >= 0; i--)
    {
        x[i] = y[i];

        for (j = i + 1; j <= n - 1; j++)
        {
            x[i] -= (u[i][j] * x[j]);
        }

        x[i] /= u[i][i];
    }

    printf("\nThe Solution of the given input is given below : \n"); // printing the solution

    for (i = 0; i < n; i++)
        printf("\n\t x%d = %lf ", i + 1, x[i]);
}


/************************************************************ !!!PART-B!!!*******************************************************************/

void partb()
{

    void cubicspline(float x[100], float y[100], int k, float xi)
    {
        int n;
        n = k - 1;
        float h[n], b[n], v[n - 1], u[n - 1], z[n + 1];
        for (int i = 0; i < n; i++)
        {
            h[i] = x[i + 1] - x[i];
        }

        for (int i = 0; i < n; i++)
        {
            b[i] = (y[i + 1] - y[i]) / h[i];
        }

        for (int i = 0; i < n - 1; i++)
        {
            v[i] = 2 * (h[i + 1] + h[i]);
        }

        for (int i = 0; i < n - 1; i++)
        {
            u[i] = 6 * (b[i + 1] - b[i]);
        }

        z[0] = z[n] = 0;

        float mat[n - 1][n - 1];
        for (int i = 0; i < n - 1; i++)
        {
            mat[i][i] = v[i];
        }

        for (int i = 0; i < n - 2; i++)
        {
            mat[i][i + 1] = h[i + 1];
            mat[i + 1][i] = h[i + 1];
        }


        // using by LU decomposition method to solve the matrix
        int row, column, ctr;
        float lowertriangle[n - 1][n - 1]; // matrix for lower triangle
        float uppertriangle[n - 1][n - 1]; // matrix for upper triangle
        for (row = 0; row < n - 1; row++)
        {
            for (column = 0; column < n - 1; column++)
            {
                if (row > column) // initialise all elements of upper triangle in upper triangular
                                  // matrix as 0
                    uppertriangle[row][column] = 0.0;
                if (row < column) // initialise all elements of lower triangle in lower triangular
                                  // matix as 0
                    lowertriangle[row][column] = 0.0;
                if (row
                    == column) // initialise all diagonal elements of lower triangular matrix as 1
                    lowertriangle[row][column] = 1.0;
            }
        }
        int ab = 0;
        int m = 1;
        while (ab < n - 1)
        {
            for (ctr = m; ctr < n - 1; ctr++)
            {
                float mul[n - 1];
                mul[ctr] = mat[ctr][ab] / mat[ab][ab];
                lowertriangle[ctr][ab] = mul[ctr];
                row = ctr;
                for (column = ab; column < n - 1; column++)
                {
                    mat[row][column] = mat[row][column] - mul[ctr] * mat[ab][column];
                }
            }
            m++;
            ab++;
        }

        printf("UPPER TRIANGULAR MATRIX\n");
        for (row = 0; row < n - 1; row++)
        { // displays upper triangular matrix
            for (column = 0; column < n - 1; column++)
            {
                uppertriangle[row][column] = mat[row][column];
                printf("%.2f ", uppertriangle[row][column]);
            }
            printf("\n");
        }
        printf("LOWER TRIANGULAR MATRIX\n");
        for (row = 0; row < n - 1; row++)
        { // displays lower triangular matrix
            for (column = 0; column < n - 1; column++)
                printf("%.2f ", lowertriangle[row][column]);
            printf("\n");
        }

        float Y[n - 1];

        // Generating intermediate matrix Y
        for (int i = 0; i < n - 1; i++)
        {
            Y[i] = u[i];
            for (int j = 0; j < i; j++)
            {
                Y[i] = Y[i] - lowertriangle[i][j] * Y[j];
            }
        }

        for (int i = n - 1; i >= 1; i--)
        {
            z[i] = Y[i - 1];
            for (int j = i + 1; j < n; j++)
            {
                z[i] = z[i] - uppertriangle[i - 1][j - 1] * z[j];
            }
            z[i] = z[i] / uppertriangle[i - 1][i - 1];
        }

        float result;


        for (int i = 1; i < n + 1; i++)
        {
            if ((xi > x[i - 1]) && (xi < x[i]))
            {

               result = z[i] * pow((xi - x[i - 1]), 3) / (6 * h[i - 1])
                  + z[i - 1] * pow((x[i] - xi), 3) / (6 * h[i - 1])
                  + ((y[i] / h[i - 1]) - ((z[i] * h[i - 1]) / 6)) * (xi - x[i - 1])
                  + ((y[i - 1] / h[i - 1]) - ((z[i - 1] * h[i - 1]) / 6)) * (x[i] - xi);//computing the result using cubic spline interpolation


                printf("value of the function at x = %.2f is y =  %.4f \n", xi, result);
            }
        }
    }


    int num;
    float x[100], y[100];
    printf("Enter the number of points : ");
    scanf("%d", &num);


    printf("\nEnter the x value and the y value of the function for the number of points mentioned above(x  and y format): \n");
    for (int i = 0; i < num; i++)
    {
        scanf("%f %f", &x[i], &y[i]);
    }
    float xi;
    printf("Enter the value of x at which the value of y is needed(x within max and min of values entered) : ");
    scanf("%f", &xi);

    cubicspline(x, y, num, xi);
}

int main()
{
    printf("\n!!!!!  PART-A !!!!! (To compute the LU Decomposition and also compute the solution "
           "using Gaussian Elimintion)\n");
    parta();

    printf("\n\n!!!!!  PART-B !!!!! (Using the above method of LU Deomposition in Cubic Spline "
           "Interpolation)\n");
    partb();
    return 0;
}
