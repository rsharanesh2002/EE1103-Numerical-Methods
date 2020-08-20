/**********************************************************************************************************************************************

                                             *********______QUIZ-2______**********
GROUP   : 1J
PURPOSE : To write a C-program to simulate the brownian motion of particles at room temperature:
           1.That takes in 3 argufloatments on the command line (D(size of particle),T(Temperature),N(number of steps))
           2.Find the motion of the particle by solving the langevin equation.
           3.Plot the RMS distance of the particle as a function of tempearture.

AUTHORS: Vibhhu Sharma(EE19B128),Sharansh R(EE19B116),Adityan CG(EE19B003)
DATE   : October 31st, 2019
INPUTS : Five inputs required:
        --->Size of the Particle(D) (in nanometer)
        --->Temperature(T)
        --->Temperature_start
        --->Temperature_stop
        --->Number of steps(N)
SAMPLE INPUTS
OUTPUTS: The Corresponding plots in GNUPLOT.

***********************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

/*************  UNITS Taken  **************
   Length in nanometer(nm)
   Mass in milligram(mg)
*******************************************/
#define PI 3.141
#define KB 13.8 // boltzman constant
#define KAPPA 0.1 // stiffness
#define ETA 8.9E-7 // viscosity of water(i.e the medium taken)


double uni_rand() // genrating a uniform distribution
{
    return (rand() + 1.0) / (RAND_MAX + 1.0);
}

double randforce(double y) // function to generate a random white noise floor force
{
    double x = uni_rand();

    if (x > 0.5)
        return y * sqrt(-2 * log(uni_rand())) * cos(2 * PI * uni_floatrand());

    else
        return (-1 * y) * sqrt(-2 * log(uni_rand())) * cos(2 * PI * uni_rand());
}


float retb(float x, float T, float g) // function to return the value of differential equation
{
    double y = 2 * g * KB * T;

    return (randforce(y) - KAPPA * x) / g;
}


float solve_runge(int N, float temp, float gamma,int p) // function that solves the differntial equation using runge-kutta method
{

    float* b = (float*)malloc(N * sizeof(float)); // array to store the solutions of the differntial equation in x

    b[0] = 0; // initial b value

    float h = 0.1; // step size of time
    float t = 1; // initial value of time
    float k1, k2, k3, k4; // constants in runge kutta's solution for b
    float rms_sum = 0, rms; // variables to calculate the rms value

    FILE* gnuplotfp = popen("gnuplot -persistent", "w"); // file pointer

    if (p == 0) // plotting the trajectory only in the first case for a fixed value of temperature
    {
        fprintf(gnuplotfp,"set title 'X vs t for an optically trapped particle' \n set xlabel 't(s)' \n "
        "set ylabel 'X(nm)'\n "); // setting the format of plot
        fprintf(gnuplotfp, "plot [][]'-' u 1:2 ls 1 title ' Trajectory' w l\n");
    }

    for (int i = 1; i < N; i++)
    {

        /*Calculating the values of all the constants required in runge kutta's final equation*/
        k1 = retb(b[i - 1], temp, gamma);
        k2 = retb(b[i - 1] + k1 * 0.5 * h, temp, gamma);
        k3 = retb(b[i - 1] + k2 * 0.5 * h, temp, gamma);
        k4 = retb(b[i - 1] + k3 * h, temp, gamma);

        b[i] = b[i - 1] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0; // runge-kutta fourth order formulae
                                                               
        if (p == 0) // plotting the trajectory only in the first case for a fixed value of temperature
            fprintf(gnuplotfp, "%f %f \n", t, b[i]);

        t = t + h; // incrementing the time
        rms_sum += b[i] * b[i]; // summing up the squares to calculate the rms value

    }

    rms = sqrt((rms_sum * h) / t); // calculating rms value

    free(b); // freeing up the allocated memory
    return rms; // returns the rms value

}

int main(int argc, char* argv[])
{
    srand(time(0));
    /* Reading in data from command line */
    if (argc != 6)
    {
        fprintf(stderr, "Five inpfloatuts required:\n\t--->Size of the Particle(enter in nanometer "
                        "\n\t--->Temperature"
                        "\n\t--->Temperature_start"
                        "\n\t--->Temperature_stop"
                        "\n\t--->Number of steps\n");
        exit(1);
    }

    float psize = atof(argv[1]); // size of the particle
    float temp = atof(argv[2]); // temperature
    float temp_start = atof(argv[3]); // start value of temperature
    float temp_stop = atof(argv[4]); // end value of temperature
    int N = atoi(argv[5]); // total number of steps to be taken

    float gamma = 3 * PI * ETA * psize; // dragging coefficient

    float rungeret; // to get the rms value

    rungeret = solve_runge(N, temp, gamma, 0); // passing with zero to get the trajectory

    FILE* gnuplotfile = popen("gnuplot -persistent", "w"); // file pointer
    fprintf(gnuplotfile, "set title 'Brownian motion in 1D' \n set xlabel 'T(K)' \n "
                         "set ylabel 'X-RMS(nm)'\n "); // setting the format of plot
    fprintf(gnuplotfile, "plot [][]'-'u 1:2 ls 1 title 'RMS vs T' w l\n");

    while (temp_start < temp_stop)
    {
        fprintf(gnuplotfile, "%f %f\n", temp_start,solve_runge(N, temp_start, gamma,1)); // writing the rms value of 'b' by varying the 											 // temperature
        temp_start += 0.5; // increasing the temperature
    }
    return 0;
}

