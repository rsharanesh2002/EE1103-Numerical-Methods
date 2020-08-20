/**********************************************************************************************************************************************
Purpose: To plot the trajectory of foucault's pendulum by solving differential equations using
		1)Euler's method
		2)Heun's method
		3)Runge-Kutta Fourth Order method
		4)RK45
Authors: Vibhhu Sharma(EE19B128),Sharansh R(EE19B116),Adithyan CG(EE19B003)
Date   : 21st October,2019
Input  : Three inputs required:
	--->Method : 1 (for) EULER
		     2 (for) EULER-HEUN
		     3 (for) RUNGE-KUTTA
		     4 (for) RK45
	--->Total time
	--->Stepsize of time

Output :Corresponding plot in gnuplot 
***********************************************************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float solve1(float omega_earth, float omegasq_shm, float x, float v)//function for the differntial equation in x
{
    return ((2 * omega_earth * v) - (omegasq_shm * x));
}

float solve2(float omega_earth, float omegasq_shm, float y, float u)//function for the differntial equation in y
{
    return (-(2 * omega_earth * u) - (omegasq_shm * y));
}

void euler(int N, float h)
{
    FILE* gnuplotfile = popen("gnuplot -persistent", "w"); // file pointer
    fprintf(gnuplotfile, "set title 'Foucalt Pendulum using EULER-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n "); // setting the format of plot
    fprintf(gnuplotfile, "splot [][]'-'u 1:2:3 ls 1 title ' Trajectory' w l\n");
    FILE* gnuplotfile1 = popen("gnuplot -persistent1", "w"); // file pointer

    fprintf(gnuplotfile1, "set title 'Foucalt Pendulum using EULER-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n "); // setting the format of plot
    fprintf(gnuplotfile1, "plot [][]'-' ls 1 title ' Trajectory' w l\n");


    float x = 0, y = 0;//initial x and y
    float u = 0.00019, v = 0.00019;//initial x an y velocities
    float omega_earth = 0.0000727;//angular velocity of rotating surface(in case of foucault pendulum, the earth)
    float omegasq_shm = 3.14;//(w0)^2 term in differential equation

    int t; //initialising variable for the time
    float time_step = h; //stepsize for time

    for (int i = 0; i < N; i++)
    {
        float utemp = u;
        float vtemp = v;
        float xtemp = x;
        float ytemp = y;

        x = xtemp + u * time_step;//solving differential equation in x
        y = ytemp + v * time_step;//solving differential equation in y

        fprintf(gnuplotfile, "%f %f %f\n", x, y,i*time_step); // writing x and y co-ordinates to a file
	fprintf(gnuplotfile1, "%f %f\n", x, y);
        u = utemp + solve1(omega_earth, omegasq_shm, x, vtemp) * time_step;//solving differential equation in u
        v = vtemp + solve2(omega_earth, omegasq_shm, y, utemp) * time_step;//solving differential equation in v
    }
}

void heun(int N, float h1)
{
    int n = N;
    float* x = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in x
    float* y = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in y
    float* u = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in u
    float* v = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in v

    float omega_earth = 0.0000727;//angular velocity of rotating surface(in case of foucault pendulum, the earth)
    float omegasq_shm = 3.14;//(w0)^2 term in differential equation
    x[0] = 0, y[0] = 0;//initial x and y
    u[0] = 0.00019, v[0] = 0.00019;//initial x an y velocities
    
    float h = h1;       // timestep
    float k1, k2;	//constants in heun's solution for u
    float k11, k22;	//constants in heun's solution for v
    float kx1, kx2;	//constants in heun's solution for x
    float ky1, ky2;	//constants in heun's solution for y

    FILE* gnuplotfile = popen("gnuplot -persistent", "w"); // file pointer
    fprintf(gnuplotfile, "set title 'Foucalt Pendulum using HEUNS-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n set zlabel 'T'\n "); // setting the format of plot
    fprintf(gnuplotfile, "splot [][]'-' u 1:2:3 ls 1 title ' Trajectory' w l\n");
    FILE* gnuplotfile1 = popen("gnuplot -persistent1", "w"); // file pointer
     fprintf(gnuplotfile1, "set title 'Foucalt Pendulum using HEUNS-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n "); // setting the format of plot
    fprintf(gnuplotfile1, "plot [][]'-' ls 1 title ' Trajectory' w l\n");


    for (int i = 1; i < n; i++)
    {
        //assigning the values to the respective constants for solving those equations
        k1 = ((2 * omega_earth * v[i - 1]) - (omegasq_shm * x[i - 1]));
        k11 = (-(2 * omega_earth * u[i - 1]) - (omegasq_shm * y[i - 1]));
        k2 = ((2 * omega_earth * (v[i - 1] + k1 * h)) - (omegasq_shm * (x[i - 1] + h)));
        k22 = (-(2 * omega_earth * (u[i - 1] + k11 * h)) - (omegasq_shm * (y[i - 1] + h)));

        u[i] = u[i - 1] + (h * (k1 + k2) * 0.5);//solving differential equation in u
        v[i] = v[i - 1] + (h * (k11 + k22) * 0.5);//solving differential equation in v

	//assigning the values to the respective constants for solving those equations
        kx1 = u[i];
        kx2 = u[i] + h;
        ky1 = v[i];
        ky2 = v[i] + ky1 * h;

        x[i] = x[i - 1] + h * 0.5 * (kx1 + kx2);//solving differential equation in x
        y[i] = y[i - 1] + h * 0.5 * (ky1 + ky2);//solving differential equation in y
    }

    for (int i = 0; i < n; i++)
    {
        fprintf(gnuplotfile, "%f %f %f\n", x[i], y[i],i*h1); // writing x and y co-ordinates to a file
        fprintf(gnuplotfile1, "%f %f\n", x[i], y[i]); // writing x and y co-ordinates to a file
    }

    free(x);
    free(y);
    free(u);
    free(v);
}

void runge(int N, float h1)
{
    int n = N;
    float* x = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in x
    float* y = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in y
    float* u = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in u
    float* v = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in v

    float omega_earth = 0.0000727;//angular velocity of rotating surface(in case of foucault pendulum, the earth)
    float omegasq_shm = 3.14;//(w0)^2 term in differential equation
    x[0] = 0, y[0] = 0;//initial x and y
    u[0] = 0.00019, v[0] = 0.00019;//initial x an y velocities

    float h = h1; // timestep
    float k1, k2, k3, k4;//constants in runge kutta's solution for u
    float k11, k22, k33, k44;//constants in runge kutta's solution for v
    float kx1, kx2, kx3, kx4;//constants in runge kutta's solution for x
    float ky1, ky2, ky3, ky4;//constants in runge kutta's solution for y

    FILE* gnuplotfile = popen("gnuplot -persistent", "w"); // file pointer
    fprintf(gnuplotfile,
        "set title 'Foucalt Pendulum using RUNGE-KUTTA-METHOD' \n set xlabel 'X' \n "
        "set ylabel 'Y'\n set zlabel 'T'\n "); // setting the format of plot
    fprintf(gnuplotfile, "splot [][]'-' u 1:2:3 ls 1 title ' Trajectory' w l\n");
    FILE* gnuplotfile1 = popen("gnuplot -persistent1", "w"); // file pointer
    fprintf(gnuplotfile1, "set title 'Foucalt Pendulum using RUNGE-KUTTA-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n "); // setting the format of plot
    fprintf(gnuplotfile1, "plot [][]'-' ls 1 title ' Trajectory' w l\n");


    for (int i = 1; i < n; i++)
    {
    
    	/*Calculating the values of all the constants required in runge kutta's final equation*/
        k1 = ((2 * omega_earth * v[i - 1]) - (omegasq_shm * x[i - 1]));
        k11 = (-(2 * omega_earth * u[i - 1]) - (omegasq_shm * y[i - 1]));
        k2 = ((2 * omega_earth * (v[i - 1] + 0.5 * k1 * h)) - (omegasq_shm * (x[i - 1] + 0.5 * h)));
        k22 = (-(2 * omega_earth * (u[i - 1] + 0.5 * k11 * h))
            - (omegasq_shm * (y[i - 1] + 0.5 * h)));
        k3 = ((2 * omega_earth * (v[i - 1] + 0.5 * k2 * h)) - (omegasq_shm * (x[i - 1] + 0.5 * h)));
        k33 = (-(2 * omega_earth * (u[i - 1] + 0.5 * k22 * h))
            - (omegasq_shm * (y[i - 1] + 0.5 * h)));
        k4 = ((2 * omega_earth * (v[i - 1] + k3 * h)) - (omegasq_shm * (x[i - 1] + h)));
        k44 = (-(2 * omega_earth * (u[i - 1] + k33 * h)) - (omegasq_shm * (y[i - 1] + h)));
        u[i] = u[i - 1] + (h * (k1 + 2 * k2 + 2 * k3 + k4) / 6);
        v[i] = v[i - 1] + (h * (k11 + 2 * k22 + 2 * k33 + k44) / 6);
        kx1 = u[i];
        kx2 = u[i] + 0.5 * h;
        kx3 = u[i] + 0.5 * h;
        kx4 = u[i] + h;
        ky1 = v[i];
        ky2 = v[i] + ky1 * 0.5 * h;
        ky3 = v[i] + ky2 * 0.5 * h;
        ky4 = v[i] + ky3 * h;
        x[i] = x[i - 1] + h * (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6.0;
        y[i] = y[i - 1] + h * (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6.0;
    }

    for (int i = 0; i < n; i++)
    {
        fprintf(gnuplotfile, "%f %f %f\n", x[i], y[i],i*h1); // writing x and y co-ordinates to a file
        fprintf(gnuplotfile1, "%f %f\n", x[i], y[i]); // writing x and y co-ordinates to a file
    }

    free(x);
    free(y);
    free(u);
    free(v);
}

void fehlberg(int N, float h1)
{
    int n = N;
    float* x = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in x
    float* y = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in y
    float* u = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in u
    float* v = (float*)malloc(n * sizeof(float));//array to store the solutions of the differntial equation in v

    float omega_earth = 0.0000727;//angular velocity of rotating surface(in case of foucault pendulum, the earth)
    float omegasq_shm = 3.14;//(w0)^2 term in differential equation
    x[0] = 0, y[0] = 0;//initial x and y
    u[0] = 0.00019, v[0] = 0.00019;//initial x an y velocities

    float h = h1; // timestep
    float k1, k2, k3, k4, k5;
    float k11, k22, k33, k44, k55;
    float kx1, kx2, kx3, kx4, kx5;
    float ky1, ky2, ky3, ky4, ky5;
    FILE* gnuplotfile = popen("gnuplot -persistent", "w"); // file pointer

    fprintf(gnuplotfile, "set title 'Foucalt Pendulum using RK45-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n set zlabel 'T'\n "); // setting the format of plot
    fprintf(gnuplotfile, "splot [][]'-' u 1:2:3 ls 1 title ' Trajectory' w l\n");
    FILE* gnuplotfile1 = popen("gnuplot -persistent1", "w"); // file pointer
    fprintf(gnuplotfile1, "set title 'Foucalt Pendulum using RK45-METHOD' \n set xlabel 'X' \n "
                         "set ylabel 'Y'\n "); // setting the format of plot
    fprintf(gnuplotfile1, "plot [][]'-' ls 1 title ' Trajectory' w l\n");


    for (int i = 1; i < n; i++)
    {
    
    	/*Calculating the values of all constants involved in the final fehlberg equation*/
        k1 = ((2 * omega_earth * v[i - 1]) - (omegasq_shm * x[i - 1]));
        k11 = (-(2 * omega_earth * u[i - 1]) - (omegasq_shm * y[i - 1]));
        k2 = ((2 * omega_earth * (v[i - 1] + 0.25 * k1 * h))
            - (omegasq_shm * (x[i - 1] + 0.25 * h)));
        k22 = (-(2 * omega_earth * (u[i - 1] + 0.25 * k11 * h))
            - (omegasq_shm * (y[i - 1] + 0.25 * h)));
        k3 = ((2 * omega_earth * (v[i - 1] + 0.09375 * k1 * h + 0.28125 * k2 * h))
            - (omegasq_shm * (x[i - 1] + 0.375 * h)));
        k33 = (-(2 * omega_earth * (u[i - 1] + 0.09375 * k11 * h + 0.28125 * k22 * h))
            - (omegasq_shm * (y[i - 1] + 0.375 * h)));
        k4 = ((2 * omega_earth * (v[i - 1] + 3.32 * k3 * h + 0.87938 * k1 * h - 3.2772 * k2 * h))
            - (omegasq_shm * (x[i - 1] + 0.923 * h)));
        k44 = (-(2 * omega_earth
                   * (u[i - 1] + 3.32 * k33 * h + 0.87938 * k11 * h - 3.2772 * k22 * h))
            - (omegasq_shm * (y[i - 1] + 0.923 * h)));
        k5 = ((2 * omega_earth * (v[i - 1] + 2.0324 * k1 * h + 7.17348928 * k3 * h - 8 * k2 * h
                                     - 0.20589669 * k4 * h))
            - (omegasq_shm * (x[i - 1] + h)));
        k55 = (-(2 * omega_earth * (u[i - 1] + 2.0324 * k11 * h + 7.17348928 * k33 * h - 8 * k22 * h
                                       - 0.20589669 * k44 * h))
            - (omegasq_shm * (y[i - 1] + h)));
        u[i] = u[i - 1] + (h * (0.11574 * k1 + 0.54892788 * k3 + 0.535723 * k4 - 0.2 * k5));
        v[i] = v[i - 1] + (h * (0.11574 * k11 + 0.54892788 * k33 + 0.535723 * k44 - 0.2 * k55));
        kx1 = u[i];
        kx2 = u[i] + 0.25 * h;
        kx3 = u[i] + 0.375 * h;
        kx4 = u[i] + 0.923 * h;
        kx5 = u[i] + h;
        ky1 = v[i];
        ky2 = v[i] + 0.25 * ky1 * h;
        ky3 = v[i] + 0.09375 * ky1 * h + 0.28125 * ky2 * h;
        ky4 = v[i] + 3.32 * ky3 * h + 0.87938 * ky1 * h - 3.2772 * ky2 * h;
        ky5 = v[i] + 2.0324 * ky1 * h + 7.17348928 * ky3 * h - 8 * ky2 * h - 0.20589669 * ky4 * h;
        x[i] = x[i - 1] + (h * (0.11574 * kx1 + 0.54892788 * kx3 + 0.535723 * kx4 - 0.2 * kx5));
        y[i] = y[i - 1] + (h * (0.11574 * ky1 + 0.54892788 * ky3 + 0.535723 * ky4 - 0.2 * ky5));
    }

    for (int i = 0; i < n; i++)
    {
        fprintf(gnuplotfile, "%f %f %f\n", x[i], y[i],i*h1); // writing x and y co-ordinates to a file
        fprintf(gnuplotfile1, "%f %f\n", x[i], y[i]); // writing x and y co-ordinates to a file
    }

    free(x);
    free(y);
    free(u);
    free(v);
}

int main(int argc, char* argv[])
{
    /* Reading in data from command line */
    if (argc != 4)
    {
        fprintf(stderr, "Three inputs required:\n\t--->Method : 1 (for) EULER\n\t\t     2 (for) "
                        "EULER-HEUN\n\t\t     3 (for) RUNGE-KUTTA\n\t\t     4 (for) RK45"
                        "\n\t--->Total time"
                        "\n\t--->Stepsize of time\n");
        exit(1);
    }
    int n = atoi(argv[1]);//getting the user choice
    int N = atoi(argv[2]);//total time
    float delta = atof(argv[3]);//stepsize for time

    if (n == 1)
        euler(N, delta);//using euler method
    if (n == 2)
        heun(N, delta);//using euler-heuns method
    if (n == 3)
        runge(N, delta);//using Runge-Kutta Fourth Order method
    if (n == 4)
        fehlberg(N, delta);//using RK45 method
    return 0;
}
