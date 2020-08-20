// Mimic a real time system that analyzes fluorescence data
// Input: read one value at a time, from a very large file
// Output: location and width of each peak, median distance between peaks

#include<stdio.h>
#include<stdlib.h>

#define MAX_ITEMS    16
typedef struct circularQueue_s
{
    int     first;
    int     last;
    int     validItems;
    int     data[MAX_ITEMS];
} circularQueue_t;

void initializeQueue(circularQueue_t *theQueue);
int isEmpty(circularQueue_t *theQueue);
int putItem(circularQueue_t *theQueue, int theItemValue);
int getItem(circularQueue_t *theQueue, int *theItemValue);
void printQueue(circularQueue_t *theQueue);

int main(){
  FILE *fp = fopen("FLdata.txt","r");
  int N = 15276000;
  int W = MAX_ITEMS;          // set a moving average window
  int thresh = 2*W;
  
  int sum = 0;
  int i,ret, dat;
  circularQueue_t dataQ;
  
  initializeQueue(&dataQ);

  for(i=0;i<W;i++){
    fscanf(fp,"%d",&dat);
    sum += dat;
    putItem(&dataQ, dat);
  }

  // printf("Sum of %d pts = %d\n",W,sum);
  do{
    getItem(&dataQ,&dat);        // get first item in circular queue
    sum = sum - dat;
    ret = fscanf(fp,"%d",&dat);  // read next item from file
    sum = sum + dat;
    putItem(&dataQ,dat);        // push last item into circuar queue

    if(sum > thresh)
      printf("%d\t%d\n",i,sum);
    i++;
  }  while(ret!=EOF);
  // printf("Read %d data points\n",i);
}


void initializeQueue(circularQueue_t *theQueue)
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

int putItem(circularQueue_t *theQueue, int theItemValue)
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
        theQueue->last = (theQueue->last+1)%MAX_ITEMS;
	return(0);
    }
}

int getItem(circularQueue_t *theQueue, int *theItemValue)
{
    if(isEmpty(theQueue))
    {
        printf("isempty\n");
        return(-1);
    }
    else
    {
        *theItemValue=theQueue->data[theQueue->first];
        theQueue->first=(theQueue->first+1)%MAX_ITEMS;
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


