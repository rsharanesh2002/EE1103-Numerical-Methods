/**********************************************************************************************************************************************
Purpose:To find the factorial of a given number
Authors:Sakthi Harish(EE19B054),Sharanesh R(EE19B116)
Date   :August 1,2019
Inputs :Positive integers
Outputs:Factorial of the input
***********************************************************************************************************************************************/
#include<stdio.h>
int main()
{
  float n; //input number
  int i; //variable for iteration
  unsigned long long nfac=1; //initialising the value of factorial
  printf("\n Enter the Integer:");
  scanf("%f",&n);
  if((int)n==n) //check for integer
   {
    if(n<0) //check for negative numbers
    printf("\n ERROR!!Factorial doesn't exist for negative numbers\n");
    else
     { 
      for(i=1;i<=n;i++)
       {
        nfac=nfac*i; //computing factorial
       }
     printf("\n Factorial of %d is %llu ",(int)n,nfac);
     }
   }
  else
  printf("\n Its time you got to know what an integer is,isn't it?\n"); //printing the valueof factorial
  return 0;
}
