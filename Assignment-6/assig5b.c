/*********************************************************************************************************************************************
Purpose: To create a phase plot of predator_prey model and dynamic popultation of predator and prey with time using euler's method.
Authors: Vibhhu Sharma(EE19B128),Sharansh R(EE19B116),Adithyan CG(EE19B003)
Date   : 20th September,2019
Input  : The required inputs from user are:
         ./a.out {alpha} {beta} {gamma} {delta} {initial prey population} {initial predator population} {toatal time} {delta_time}
         Sample : ./a.out 1 2 3 4 1 1 50 0.01
         Note : The given input must be in a optimal range.(Also the delta_time must be small enough to get a optimal plot).
Output : The phase plot of predator_prey model and dynamic popultation of predator and prey with time.
**********************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>

float f(float x, float y, float a, float b)// function to compute derivative
   {return ((a*x)-(b*x*y));}
   
float g(float x, float y, float c, float d)// function to compute derivative
   {return ((d*x*y)-(c*y));}
  
int main(int argc, char* argv[])
{

    /* Reading in data from command line */
    if (argc != 9)
    {
        fprintf(stderr, "Eight inputs required (alpha, beta, gamma, delta, initial prey population, initial "
                        "predator population, total time, change in time)\n");
        exit(1);
    }

    float a, b, c, d;
    float x0, y0, y, x, h, tn;
    float mx, my; // derivatives
    a = strtod(argv[1], NULL);//alpha
    b = strtod(argv[2], NULL);//beta
    c = strtod(argv[3], NULL);//gamma
    d = strtod(argv[4], NULL);//delta
    x0 = strtod(argv[5], NULL);//initial prey population
    y0 = strtod(argv[6], NULL);//initial predator population
    tn = strtod(argv[7], NULL);//total time
    h = strtod(argv[8], NULL);//change in time

    FILE* gnuplotphaseplot = popen("gnuplot -persistent", "w");//file to plot the phase diagram in gnuplot
    FILE* gnuplottime = popen("gnuplot -persistent", "w");//file to plot time dependence in gnuplot

    FILE* lvdata = fopen("lvdata.txt", "w"); // storing the data in a file named lvdata.txt

    fprintf(gnuplotphaseplot, "set title 'Predator versus Prey Phaseplot' \n set xlabel 'Prey' \n "
                              "set ylabel 'Predator'\n ");//setting the format of plot 

    fprintf(gnuplottime, "set title 'Dynamic population of Prey and Predator' \n set xlabel "
                         "'Time' \n set ylabel 'Predator (or) Prey'\n ");//setting the format of plot

    fprintf(gnuplotphaseplot, "plot '-' with d ls 1 title 'Predator Prey Phaseplot'\n");

    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

    x = x0; // initial number of prey
    y = y0; // initial number of predators
    float t = 0; // time

    
    //iteration starts for euler's method
    
    while (t < tn)
    {
        mx = f(x, y, a, b);
        my = g(x, y, c, d);

        y = y + my * h;
        x = x + mx * h;

        if (x <= 0)
            x = 0;
        if (y <= 0)
            y = 0;

        t = t + h;//incrementing time

        fprintf(gnuplotphaseplot, "%f %f\n", x, y); //writing to file
        fprintf(lvdata, "%f %f %f\n", t, x, y); //writing to file
    }

    fprintf(gnuplotphaseplot, "e\n");
    fprintf(gnuplottime, "plot 'lvdata.txt' using 1:2 w d ls 2 title 'Prey', 'lvdata.txt' using 1:3 w d ls 4 title 'Predator'\n");

    // making the gnuplot window open without terminating the C program
    fflush(gnuplottime);
    fflush(gnuplotphaseplot);

    fclose(lvdata);//closing the textfile

}

