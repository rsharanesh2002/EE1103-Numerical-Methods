/**********************************************************************************************************************************************
 Purpose: Mimic a real time system that analyzes fluorescence data
 
 Authors: Vibhhu Sharma(EE19B128),Sharansh R(EE19B116),Adithyan CG(EE19B003)
 
 Date   : 27th October,2019
 
 Input: Read one value at a time, from a very large file "CompleteData_fsc_ssc_10um.txt"
 
 Output: No. of peaks seen common in both graphs, mean width of these peaks, mean and median distance between consecutive peaks
 	Sample output:
	Median distance between maxima=2114.000000units(time)
	No of peaks=298
	Average width of maxima=16.104027 units(time)
	Average distance between maxima=3728.000000 units(time)
	Rate of input=1 ms
	In these units, average and median distance=3.728863 and 2.114000 seconds
***********************************************************************************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX_ITEMS   20
#define horizontal_threshold 300
#define vertical_threshold 150
typedef struct circularQueue_s//defining basic structure of the circular queue
{
    int     first;
    int     last;
    int     validItems;
    int     data[MAX_ITEMS];
} circularQueue_t;

void initializeQueue(circularQueue_t *theQueue);//function to initialise the queue
int isEmpty(circularQueue_t *theQueue);//checks if the queue is empty
int putItem(circularQueue_t *theQueue, float theItemValue);//used to insert an item into the queue
int getItem(circularQueue_t *theQueue, float theItemValue);//used to obtain an item from the queue
void printQueue(circularQueue_t *theQueue);//prints the elements stored within the queue

int main(){
  FILE *fp = fopen("CompleteData_fsc_ssc_10um.txt","r");
  FILE *fp1=fopen("last1.txt","w");//file to store data about forward deflection
  FILE *fp2=fopen("last2.txt","w");//file to store data about sideways deflection
  int N = 1116000;//total number of lines in the file
  int* res = (int*)malloc(N * sizeof(int));
  int W = MAX_ITEMS;//No. of items in the queue at one time
  
  			/*The circular queue always has exactly 20 items. Whenever we nedd to ad an item, we remove the first item from the 				queue and each element is then shifted to the left*/
  float sum = 0,sum1=0;
  int i;
  float ret,dat,bat;
  int counter=0;//counter variable
  circularQueue_t dataQ;
  circularQueue_t dataQ1;
  initializeQueue(&dataQ);
  initializeQueue(&dataQ1);
  for(i=0;i<W;i++)
  {
    fscanf(fp,"%f %f",&dat,&bat);
    fprintf(fp1,"%d %f\n",i,dat);
    fprintf(fp2,"%d %f\n",i,bat);
    sum += abs((int)dat-2675);//the variable sum stores the sum of absolute values of difference between each number and 2675(Near about mean 					of non extrema data)(vertical)
    sum1+=abs((int)bat);	//the variable sum1 stores the sum of every 20 elements(horizontal)
    putItem(&dataQ, dat);
    putItem(&dataQ1, dat);
  }//till now I have just put W elements into the queue

  do
  {
    getItem(&dataQ,dat);
    getItem(&dataQ,bat);     // get first item in circular queue
    sum = sum - abs((int)dat-2675);//removing the element of the queue from sum
    sum1=sum1-abs((int)bat);  
    ret=fscanf(fp,"%f %f",&dat,&bat); // read next item from file
    sum = sum + abs((int)dat-2675);//using the next element from the file in place of the recently removed element
    sum1=sum1+abs((int)bat);
    if(sum1>=horizontal_threshold && sum >= vertical_threshold)//condition for maxima to exist on both curves
		{		
		res[counter]=i;//array res stores the location all the extrema
		counter++;//counts the number of maxima(even consecutive maxima are counted separate)
		}
    putItem(&dataQ,dat);        // push last item into circuar queue
    putItem(&dataQ,bat);
        i++;
    fprintf(fp1,"%d %f\n",i,dat);
    fprintf(fp2,"%d %f\n",i,bat);
  }  
  while(ret!=EOF);//loop executes till the end of the file is reached
  
  int c=0;
  int d=1,ct=0;
  int* width = (int*)malloc(counter * sizeof(int));//array to store the distance between consecutive extrema
  int* start = (int*)malloc(counter * sizeof(int));//array to store the starting position of an extrema
  int finalsum=0;
  printf("Starting points of various peaks in ms\n");
  for(int j=0;j<counter;j=j+d)
  {
  d=1;
	  while(1)
	  {
		  if(res[j+d]-res[j+d-1]>=10 && d>1)//condition to locate immediately adjacent block extrema
		  {
		  break;
		  }
	  d++;
	  }
	  if((j+d)>counter)//to avoid array index out of bounds error
	  break;
  width[ct]=res[j+d]-res[j];//distance between consecutive extrema
  finalsum=finalsum+width[ct];
  start[ct]=res[j];//starting point of an extrema, point from where extrema first crosses threshold value
  printf("%d\n",start[ct]);
  ct++;
  }
  /*Sorting the widths in descending order*/
  int temp=0;
  for(i=0 ; i<ct ; i++)
    {
        for(int j=0 ; j<ct-1 ; j++)
        {
            if(width[j]>width[j+1])
            {
                temp        = width[j];
                width[j]    = width[j+1];
                width[j+1]  = temp;
            }
        }
    }
    float median=0;
    /*Two cases for median*/
    if(ct%2 == 0)
        median = (width[(ct-1)/2] + width[ct/2])/2.0;
    // if number of elements are odd
    else
        median = width[ct/2];
  printf("Median distance between maxima=%funits(time)\n",median);
  printf("No of peaks=%d\n",ct);
  printf("Average width of maxima=%f units(time)\n",(float)(1116000-finalsum)/(ct));
  printf("Average distance between maxima=%f units(time)\n",(float)(finalsum/(ct)));
  printf("Rate of input=1 ms\n");
  printf("In these units, average and median distance=%f and %f seconds\n",(float)((float)finalsum/(1000*ct)),(median/1000.0));
 free(start);
 free(width);
 free(res);
  return 0;
}


void initializeQueue(circularQueue_t *theQueue)//initialises the queue with default values
{
    int i;
    theQueue->validItems  =  0;
    theQueue->first       =  0;
    theQueue->last        =  0;
    for(i=0; i<MAX_ITEMS; i++)
    {
        theQueue->data[i] = 0;
    }        
    return;
}

int isEmpty(circularQueue_t *theQueue)
{
    if(theQueue->validItems==0)
        return(1);
    else
        return(0);
}

int putItem(circularQueue_t *theQueue, float theItemValue)
{
    if(theQueue->validItems>=MAX_ITEMS)
    {
        printf("The queue is full\n");
        printf("You cannot add items\n");
        return(-1);
    }
    else
    {
        theQueue->validItems++;
        theQueue->data[theQueue->last] = theItemValue;
        theQueue->last = (theQueue->last+1)%MAX_ITEMS;//implementing circular ques and pushing each element one unit to its left
	return(0);
    }
}

int getItem(circularQueue_t *theQueue, float theItemValue)
{
    if(isEmpty(theQueue))
    {
        printf("isempty\n");
        return(-1);
    }
    else
    {
        theItemValue=theQueue->data[theQueue->first];
        theQueue->first=(theQueue->first+1)%MAX_ITEMS;//implementing circular ques and pushing each element one unit to its left
        theQueue->validItems--;
        return(0);
    }
}

void printQueue(circularQueue_t *theQueue)
{
    int aux, aux1;
    aux  = theQueue->first;
    aux1 = theQueue->validItems;
    while(aux1>0)
    {
        printf("Element #%d = %d\n", aux, theQueue->data[aux]);
        aux=(aux+1)%MAX_ITEMS;
        aux1--;
    }
    return;
}


