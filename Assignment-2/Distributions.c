/***********************************************************************************************************
Purpose	    :To generate various types of distribution and plot ther graphs
Authors     :SHARANESH R(EE19B116),SAKTHI HARISH D T(EE19B054),TMVS GANESH(EE19B124)
Date        :27/08/2019
Input       :Number ofrandom numbers,Type of distribution
Output      :Graph of the particular distribution
Input Format:./a.out <N>(Number of random numbers) <Type_of_distribution>{normal,rayleigh,cauchy,poissonian}
************************************************************************************************************/
#include <stdlib.h>
#include <math.h>
#include<stdio.h>
#include<time.h>
#include<string.h>
#include<limits.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif

void normal(double);// generating normal distribution
void rayleigh(double);// generating rayleigh distribution
void cauchy(double);// generating cauchy distribution
void poissonian(double);//generating poissonian distribution
double random_normal();
double random_rayleigh();
double random_cauchy();
double fac(int);
double drand();


int main(int argc,char* argv[])
{
   double N;

   N=atof(argv[1]);

	if(strncasecmp(argv[2],"normal", 20)==0)
		normal(N);

  	if(strncasecmp(argv[2],"rayleigh", 20)==0)
		rayleigh(N);

  	if(strncasecmp(argv[2],"poissonian", 20)==0)
		poissonian(N);

  	if(strncasecmp(argv[2],"cauchy", 20)==0)
	cauchy(N);

   return 0;

}

double drand()
{
   return (rand()+1.0)/(RAND_MAX+1.0);
}

double random_normal()
{
   return sqrt(-2*log(drand())) * cos(2*PI*drand());
}

double random_rayleigh()
{
   return sqrt(-2.0*log(drand()));
}

double random_cauchy()
{
   return tan(PI*(drand()-0.5));
}

double fac(int p)// computing factorial
{
   int t;
   double fac=1.0;

	for(t=1;t<=p;t++)
		fac=fac*t;

   return fac;
}


void normal(double N)
{
   srand(time(0)); //time seeder

   int i,j;

   double min=INT_MAX,max=0,nbin,binsize;

   double *n,*cbin,*bin,*mbin;

   FILE *fp;
   fp=fopen("nrand.txt","w");

   n=(double*)malloc(N*sizeof(double)); //dynamic memory allocation
   cbin=(double*)malloc(N*sizeof(double));
   bin=(double*)malloc(N*sizeof(double));
   mbin=(double*)malloc(N*sizeof(double));

  
	for (i=0; i<N; i++)
		n[i] = 0.5*random_normal(); // generating normal distribution

	for(i=0;i<N;i++) //calculating min and max
	   {

		if(n[i]<min)
		min=n[i];
		if(n[i]>max)
		max=n[i];

	   }

  nbin=1+(3.322*(log10(N)));// calculating bin sizes
  binsize=(max-min)/nbin;

	for(j=0;j<nbin;j++)  //assigning values into the bins
		bin[j]=min+(j*binsize);
	
	for(j=0;j<nbin;j++)//generating a histogram
	   {

	        cbin[j]=0;

   		for(i=0;i<N;i++)
		   {

	            if((n[i]>bin[j])&&(n[i]<bin[j+1]))
	               cbin[j]++;
	       	   }
	    }

	for(j=0;j<nbin;j++)
	   mbin[j]=bin[j] + (0.5*binsize);

        for(j=0;j<nbin;j++)
	   fprintf(fp,"%.2lf  %.5lf \n",mbin[j],cbin[j]/N);

       	free(n);
	free(cbin);
	free(bin);
	free(mbin);
}


void rayleigh(double N)
{
  srand(time(0)); //time seeder

  int i,j;

  double min=INT_MAX,max=0,nbin,binsize;
  double *r,*cbin,*bin,*mbin;

  FILE *fp;
  fp=fopen("rrand.txt","w");
    
  r=(double*)malloc(N*sizeof(double)); //dynamic memory allocation
  cbin=(double*)malloc(N*sizeof(double));
  bin=(double*)malloc(N*sizeof(double));
  mbin=(double*)malloc(N*sizeof(double));

  
	for (i=0; i<N; i++)
	   r[i] = random_rayleigh(); // generating rayleigh distribution
	 
	for(i=0;i<N;i++) //calculating min and max
	   {
	      if(r[i]<min)
              min=r[i];
	      if(r[i]>max)
	      max=r[i];
	   }

   nbin=1+(3.322*(log10(N))); // calculating bin sizes
   binsize=(max-min)/nbin; 

	for(j=0;j<nbin;j++) //assigning values into the bins
	    bin[j]=min+(j*binsize);
	
        for(j=0;j<nbin;j++) //generating a histogram
  	   {

           cbin[j]=0;

           for(i=0;i<N;i++)
  	      {

	        if((r[i]>bin[j])&&(r[i]<bin[j+1]))
 	        	cbin[j]++;
	      
       	      }
	    }

	for(j=0;j<nbin;j++)
 	   mbin[j]=bin[j] + (0.5*binsize);


	for(j=0;j<nbin;j++)
	    fprintf(fp,"%.2lf  %.5lf \n",mbin[j],cbin[j]/N);
	
	
	free(r);
	free(cbin);
	free(bin);
	free(mbin);
   
}



void cauchy(double N)
{
  srand(time(0)); //time seeder

  int i,j;

  double min=INT_MAX,max=0,nbin,binsize;
  double *c,*cbin,*bin,*mbin;

  FILE *fp;
  fp=fopen("crand.txt","w");    

  c=(double*)malloc(N*sizeof(double)); //dynamic memory allocation
  cbin=(double*)malloc(N*sizeof(double));
  bin=(double*)malloc(N*sizeof(double));
  mbin=(double*)malloc(N*sizeof(double));
 
 
   for (i=0; i<N; i++)
	c[i] = random_cauchy(); // generating cauchy distribution
	
   for(i=0;i<N;i++) //calculating min and max
	{
	   if(c[i]<min)
           min=c[i];
	   if(c[i]>max)
	   max=c[i];
	}

  nbin=1+(3.322*(log10(N))); // calculating bin sizes
  binsize=(max-min)/nbin; 

    for(j=0;j<nbin;j++)//assigning values into the bins
	bin[j]=min+(j*binsize);
	
    for(j=0;j<nbin;j++) //generating a histogram
  	{
           cbin[j]=0;

        for(i=0;i<N;i++)
  	   {

 	      if((c[i]>bin[j])&&(c[i]<bin[j+1]))
 	      	  cbin[j]++;
	      
	    }
	}

    for(j=0;j<nbin;j++)
         mbin[j]=bin[j] + (0.5*binsize);

    for(j=0;j<nbin;j++)
 	 fprintf(fp,"%.2lf  %.5lf \n",mbin[j],cbin[j]/N);
	

 	free(c);
	free(cbin);
	free(bin);
	free(mbin);
}


void poissonian(double N)
{

  int i;
  double *fact,*p;
 
  FILE *fp1,*fp2,*fp3,*fp10;

  fp1=fopen("pois1.txt","w");
  fp2=fopen("pois2.txt","w");
  fp3=fopen("pois3.txt","w");
  fp10=fopen("pois10.txt","w");

  fact=(double*)malloc(N*sizeof(double));
  p=(double*)malloc(N*sizeof(double));

  	for(i=0;i<N;i++)
	   {

		p[i]=pow(1,i)*exp(-1)/fac(i);
		fprintf(fp1,"%d %.3lf \n",i+1,p[i]);

           }

        for(i=0;i<N;i++)
  	   {

		p[i]=pow(2,i)*exp(-2)/fac(i);
		fprintf(fp2,"%d %.3lf \n",i+1,p[i]);

            }
	for(i=0;i<N;i++)
	    {

		p[i]=pow(3,i)*exp(-3)/fac(i);
		fprintf(fp3,"%d %.3lf \n",i+1,p[i]);

 	    }
	for(i=0;i<N;i++)
 	    {

		p[i]=pow(10,i)*exp(-10)/fac(i);
		fprintf(fp10,"%d %.3lf \n",i+1,p[i]);

	     }

	free(p);
	free(fact);
}

