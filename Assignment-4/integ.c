/*********************************************************************************************************************************************
Purpose: To compute the integration of given datafile("out1_test0.csv") and compare the various method(Box, Trapezium, Simpson's, Romberg's) of integration and to plot a graph to compare the various methods of integration.
Authors: Sharansh R(EE19B116),Sakthi Harish D.T(EE19B054),TMVS Ganesh(EE19B124)
Date   : 12th September,2019
Output : The value of the integrals using various method and the graph comparing the various methods of integration.
         "The value of the integral using romberg method is 7.148271 

	  The value of the integral using Box method is 7.148276 
	  	Error in Box method is : 0.000005 

	  The value of the integral using Trapezoidal method is 7.148260 
	        Error in Trapezoidal method is : 0.000010 

          The value of the integral using Simpson's method is 7.143318 
		Error in Simpsons method is : 0.004953"
**********************************************************************************************************************************************/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double func(double x,double *arr) //Function for calculating the interpolated value using lagrange polynomial
    {                               
	double fu=0;
	int y=(int)x;

        if(y==x)       //fitting a straight line for the last two points
	   return arr[y];

	if(y==500)        
		fu=arr[500]+((arr[501]-arr[500])*(x-500));

	else           //fitting quadraric function for every three consecutive poins using lagrange's formula 
	{
		double px;
		int x1=(y%2==0)?y:y-1; 

		for(int i=0;i<3;i++)
		{
			px=arr[x1+i];

			for(int j=0;j<3;j++)
			{
			   if(i==j)
				continue;
			   px=px*(x-(x1+j))/(i-j);
			}

			fu=fu+px;
		}
		
	}
	return fu; //returning the function value ater interpolating
}


/***************************************************
Function to calculate integral using Romberg's metod
***************************************************/
double romberg(double (*f), double  a, double  b, size_t max_steps)
{
   double R1[max_steps], R2[max_steps];  
   double *Pr = &R1[0], *Cr = &R2[0]; //Pr is previous row, Cr is current row
   double h = (b-a),ac=0.000001; //Step size and accuracy

   Pr[0] = (f[(int)a] + f[(int)b])*h*.5; //first trapezoidal step


   for(size_t i = 1; i < max_steps; ++i)
   {
      h /= 2.0;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)

      for(size_t j = 1; j <= ep; ++j)
      {
         c += func(a+(2*j-1)*h,f);
      }
      Cr[0] = h*c + 0.5*Pr[0]; //R(i,0)

      for(size_t j = 1; j <= i; ++j)
      {
         double n_k = pow(4, j);
         Cr[j] = (n_k*Cr[j-1] - Pr[j-1])/(n_k-1); //compute R(i,j)
      }

      if(i > 1 && fabs(Pr[i-1]-Cr[i]) < ac)
      {
         return Cr[i-1];
      }

      //swap Pr and Cr as we only need the last row
      double *tr = Pr;

      Pr = Cr;
      Cr = tr;
   }

   return Pr[max_steps-1]; //return integral value

}


/**************************************************
Function to calculate integral using Simpson's rule
**************************************************/
double simpsons(double n,double arr[])
{
	double o=0,e=0;
	int c=1;
	double integ=arr[0]+arr[501];
	for(double i=n;i<501;i=i+n)
	{
		if(c%2==0)
		 e=e+func(i,arr);
		else
		 o=o+func(i,arr);
		c=c+1;
	}
	integ=n*(integ+4*o+2*e)/3;
	return integ;
}


/****************************************************
Function to calculate integral using Trapezoidal rule
****************************************************/
double trapezoidal(double n,double arr[])  
{
	double integ=arr[0]+arr[501];
	double s=0;
	for(double i=n;i<501;i=i+n)
		s=s+func(i,arr);
	integ=n*(integ+2*s)/2;
	return integ;
}


/**********************************************
Function to calculate integral using Box method
***********************************************/
double box(double n,double arr[]) 
{
	double integ=0;
	for(double i=0;i<501;i=i+n)
		integ=integ+(func(i+(n/2),arr));
	integ=n*integ;
	return integ;
}


int main()
{
	int n=502;
	double arr[502];

	FILE *ptr;
	ptr=fopen("out1_test0.csv","r"); //to copy the test file into an array 

	for(int i=0;i<n;i++)
		fscanf(ptr,"%lf",&arr[i]);
        
	printf("The value of the integral using romberg method is %lf \n",romberg(arr,0,501,100)); //integrating using romberg method 

	printf("\nThe value of the integral using Box method is %lf \n",box(1,arr)); //integrating using Box method 
        printf("\tError in Box method is : %lf \n",fabs(romberg(arr,0,501,100)-box(1,arr))); //error in box method

	printf("\nThe value of the integral using Trapezoidal method is %lf \n",trapezoidal(1,arr)); //integrating using Trapezoidal method 
        printf("\tError in Trapezoidal method is : %lf \n",fabs(romberg(arr,0,501,100)-trapezoidal(1,arr))); //error in trapezium 

	printf("\nThe value of the integral using Simpson's method is %lf \n",simpsons(1,arr)); //integrating using Simpsons method 
        printf("\tError in Simpsons method is : %lf \n",fabs(romberg(arr,0,501,100)-simpsons(1,arr))); //error in simpsons method

	fclose(ptr);

	FILE *ptr1,*ptr2,*ptr3; //file pointers 

	ptr1=fopen("boxg.txt","w");
	ptr2=fopen("trag.txt","w");
	ptr3=fopen("simg.txt","w");

	for(double j=0.1;j<=100;j+=0.1) //to find the value of integral for different step sizes
	{
		fprintf(ptr1,"%lf\t%lf\n",j,box(j,arr));
		fprintf(ptr2,"%lf\t%lf\n",j,trapezoidal(j,arr));
		fprintf(ptr3,"%lf\t%lf\n",j,simpsons(j,arr));
	}
	fclose(ptr1);
	fclose(ptr2);
	fclose(ptr3);		
return 0;
}
