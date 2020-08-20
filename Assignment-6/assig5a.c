/*********************************************************************************************************************************************
Purpose: To create a logistic map in an intresting range of r entered by user.
Authors: Sharansh R(EE19B116),Vibhhu Sharma(EE19B128),Adithyan CG(EE19B003)
Date   : 20th September,2019
Input  : The initial value and final value of r along with the delta_r 
         ./a.out {initial r} {final_r} {delta_r}
Output : The logistic map in the range of r as entered by the user.
**********************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    /* Reading in data from command line */
    if (argc != 4)
    {
        fprintf(stderr, "Three inputs required (rstart, rend, delta_r).\n");
        exit(1);
    }

        double r_start = strtod(argv[1],NULL);//inputing starting value of r

        double r_end = strtod(argv[2],NULL);//inputing ending value of r

        double r_delta = strtod(argv[3],NULL);//inputing value of delta r

        double r; //for iteration purposes

        FILE* gnuplotfile = popen("gnuplot -persistent", "w"); //file pointer
	
	fprintf(gnuplotfile, "set title 'Logistic Map' \n set xlabel 'r' \n "
                              "set ylabel 'x'\n ");//setting the format of plot 
        fprintf(gnuplotfile, "plot [][-0.25:1]'-' with d ls 1 title 'Logistic map'\n");

        FILE* gnuplotPipe = popen("gnuplot -persistent", "w");

        double x, x0 = 0.5;// initial values of x

        int m, M = 10000;// variabe to iterate number of times
        

        for (r = r_start; r < r_end; r = r + r_delta) //starting iteration
        {

            if (r <= 1)// if r<=1 then the population will eventually die at some time
            {

                x = 0;
                fprintf(gnuplotfile, "%lf %lf \n", r, x); //writing output to file

            }

            else if ((r > 1) && (r < 3))// between r=3 and r=4 the population converges towards (r-1)/r
            {

                x = (r - 1) / r;
                fprintf(gnuplotfile, "%lf %lf \n", r, x); //writing output to file

            }

            else //for other values of r
            {

                x = x0; //initialising values of x

                for (m = 0; m < M; ++m) // iterating m number of times
                {

                    x = r * x * (1 - x);
                    fprintf(gnuplotfile, "%lf %lf \n", r, x); //writing output to file

                }

            }

        }
        
        fflush(gnuplotfile);
        

        return 0;
    }
