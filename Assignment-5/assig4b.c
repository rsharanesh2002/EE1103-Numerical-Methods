/**********************************************************************************************************************************************
Purpose: To create a model equation satisfying the given set of data points and to find out the output for particular given point.
Authors: Sharansh R(EE19B116),Vibhhu Sharma(EE19B128),Adithyan CG(EE19B003)
Date   : 17th September,2019
Input  : The dataset and the input values for computing a particular case.
Output : The generated model and the output for a particular case.
          "The model created by using the given set of data points is given below :
 	   Area = 0.013659 * h^( 0.521803 ) * w^( 0.519018 )
	   Enter height and weight
	   187 78

	   The area for the given h= 187.000000 cm and w=78.000000 kg is :2.008646"

***********************************************************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int i, j, k;// iteration variables
    int r = 2; // number of dependent quantities
    int q = 3; // q=r+1
    int n = 9; // number of inut values

    double coefarr[r + 1][r + 1], lowarr[r + 1][r + 1], upparr[r + 1][r + 1], coef[r + 1],
        coefpre[r + 1], sum = 0.0, arr[r + 1]; //matrices to hold the lu decomposition

    double  pref, hpow, wpow;//variables to compute the equation

    double h[10] = { log10(182), log10(180), log10(179), log10(187), log10(189), log10(194),
        log10(195), log10(193), log10(200) }; // input values of height

    double w[10] = { log10(74), log10(88), log10(94), log10(78), log10(84), log10(98), log10(76),
        log10(86), log10(96) }; // input values of weight

    double a[10] = { log10(1.92), log10(2.11), log10(2.15), log10(2.02), log10(2.09), log10(2.31),
        log10(2.02), log10(2.16), log10(2.31) }; // input values of area

    double aout; 

    // initialising the coefficent array
    coefarr[0][0] = n;
    for (i = 0; i < n; i++)
    {

        coefarr[0][1] += h[i];
        coefarr[0][2] += w[i];
        coefarr[1][0] += h[i];
        coefarr[1][1] += h[i] * h[i];
        coefarr[1][2] += h[i] * w[i];
        coefarr[2][0] += w[i];
        coefarr[2][1] += h[i] * w[i];
        coefarr[2][2] += w[i] * w[i];
        coef[0] += a[i];
        coef[1] += h[i] * a[i];
        coef[2] += a[i] * w[i];

    }

    /// Dolittle's algorithm for LU Decomposition
    for (i = 0; i < q; i++)
    {

        for (k = i; k < q; k++)/// computing upper triangular matrix
        {
            sum = 0.0;

            for (j = 0; j < i; j++)
                sum += lowarr[i][j] * upparr[j][k];

            upparr[i][k] = coefarr[i][k] - sum;

        }

        for (k = i; k < q; k++)/// computing lower triangular matrix
        {

            if (i == k)
                lowarr[i][i] = 1.0;

            else
            {
                sum = 0.0;

                for (j = 0; j < i; j++)
                    sum += lowarr[k][j] * upparr[j][i];

                lowarr[k][i] = (coefarr[k][i] - sum) / upparr[i][i];
            }

        }
    }

    // solving by forward substitution
    coefpre[0] = coef[0] / lowarr[0][0];

    for (i = 1; i < q; i++)
    {

        coefpre[i] = coef[i];

        for (j = 0; j < i; j++)
        {
            coefpre[i] -= (lowarr[i][j] * coefpre[j]);
        }

        coefpre[i] /= lowarr[i][i];

    }


    // solving by backward substitution
    arr[q - 1] = coefpre[q - 1] / upparr[q - 1][q - 1];

    for (i = q - 2; i >= 0; i--)
    {
        arr[i] = coefpre[i];

        for (j = i + 1; j <= q - 1; j++)
            arr[i] -= (upparr[i][j] * arr[j]);

        arr[i] /= upparr[i][i];
    }

    /// assigning the values for the variables in the equation
    pref = pow(10, arr[0]);
    hpow = arr[1];
    wpow = arr[2];

    /// printing the required output
    printf("The model created by using the given set of data points is given below :\n \t");

    printf("Area = %lf * h^( %lf ) * w^( %lf )", pref, arr[1], arr[2]);
    double height,weight;
    printf("\nEnter height and weight\n");
    scanf("%lf %lf",&height ,&weight);	
    aout = pref * pow(height, arr[1]) * pow(weight, arr[2]);
    printf("\nThe area for the given h= %lf cm and w=%lf kg is :%lf\n",height,weight,aout);

    return 0;
}

