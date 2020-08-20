/**********************************************************************************************************************************************

                                             *********______QUIZ-1______**********           
GROUP   : 1J                                         
PURPOSE : To write a C-program to:
           1.That takes in 3 arguments on the command line (after a.out. N mean std)
           2.Converts each of the arguments into integers and floats (as needed)
           3.Generates N random numbers with given mean and stdev
           4.Writes the random numbers in to a file
           5.Reads the random numbers from the file
           6.Calculates the mean and stdev and returns the difference in values between that expect and that calculated
AUTHORS: SHARANESH R(EE19B116);SAKTHI HARISH D T(EE19B054);TMVS GANESH (EE19B124)
DATE   : SEPTEMBER 12th,2019
INPUTS : *Number of values* *ACTUAL MEAN* *ACTUAL STANDARD DEVIATION*
OUTPUTS: The differnce between the actual and expected mean and standard deviation
SAMPLE OUTPUT:
              
             "The difference between the actual mean and expected mean is : 0.027927
              The difference between the actual standard deviation and expected deviation is : 0.028846"

***********************************************************************************************************************************************/    

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#define PI 3.14159265358979323846


/**********************************
Function to generate random numbers
**********************************/
float drand()
{

    double temp = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2*log(temp)) * cos(2*PI*temp);

}

int main(int argc, char* argv[])
{
    srand(time(0));  //to seed the random number generator

    int N;   //number of random numbers

    float std, mean;  //expexted mean and standrad deviation
    float temp,temp1;  //temporary variables for file purposes
    float sum = 0, sumsq = 0;  //sum and sum of squares of random numbers
    float am, as;  //actual mean and standrad deviation

    FILE *fp1,*fp2;   //file pointers
    fp1 = fopen("rand.txt", "w");  //opening file in write mode

    N = atoi(argv[1]);  //getting number from user
    mean = atof(argv[2]);  //getting expected mean from user
    std = atof(argv[3]);  //getting expected standard deviation from user
    
    for (int i = 0; i < N; i++)  //"for" loop for iterating
    {
        temp = std*drand() + mean;  //generating random numbers with particular mean and standard deviation
        fprintf(fp1, "%f \n", temp);  //printing to the file
       
    }

    fclose(fp1);  //closing file

    fp2= fopen("rand.txt", "r");  //opening file in read mode

    while(fscanf(fp2, "%f", &temp1)!=EOF)  //reading the random number again from the file
    {  
        sum += temp1;  //calculating sum of random numbers
        sumsq += temp1 * temp1;  //calculating sum of squares of random numbers
    }  

    fclose(fp1);  //closing file 

    am = sum / N;  //calculating actual mean
    as = sqrt((sumsq / N) - (am * am));  //calculating actual standrad deviation

    //printing the required output
    printf("\nThe difference between the actual mean and expected mean is : %f", fabs(am-mean));
    printf("\nThe difference between the actual standard deviation and expected deviation is : %f", fabs(as-std));
    return 0;
}

