//medfilt_new.c    Median filter
//http://www.embedded.com/2000/0011/0011feat3.htm
//thoroughly modified by DST 032508
// NOTE: assumes odd TAU?

#include <stdio.h>
#include <math.h>
//window size
//#define TAU 9

//Marker for chain end; Smaller than any datum 
#define STOPPER -10000




//Initialize circular linked list buffer
/*
    head->[________]->tail-|
            ^
            oldest datapoint
            
    Scan from head; add new datum in sorted location.
    Scanning is done at twice the speed of median pointer iteration;
    which is why the code is repeated...
*/
struct node
{
    struct node   *next;
    double  value;
};

struct node* buffer;             //datapoint buffer
struct node  tail={NULL,STOPPER};     //list tail (virtual node)
struct node  head={&tail,0};          //list head (virtual node)
struct node* oldest;           //p* to oldest datum

int TAU;  //window size





double mfilt(double datum);

double* medianfilter( double* data, const int Npoints, const int szWindow )
{
    int i;
	TAU = szWindow;
    const int lag = int(TAU/2);
    
    //initialize the buffer
	buffer = new node[TAU];
	oldest = buffer;
    head.next = buffer;

    for(i=0; i<TAU; ++i)
    {
        buffer[i].value = data[0];
        buffer[i].next = buffer+i+1;
    }
    buffer[TAU-1].next = &tail;
    
    //fill with data until it reaches the median position
    for(i=0; i<lag; ++i)
        mfilt( data[i] );
    
    //Collect the results; last datapoints using zero padding
    for(i=lag; i<Npoints; ++i)
        data[i-lag] = mfilt( data[i] );
        
    for(i=Npoints; i<Npoints+lag; ++i)
        data[i-lag] = mfilt( data[Npoints-1] );
    
	delete [] buffer;
    return data;
}



//median filter function.  returns median of last TAU datapoints sent
double mfilt(double datum)
{

    //variables declared for this call
    struct node *successor;     //p* to successor of replaced data item
    struct node *scan;          //p* used to scan down the sorted list
    struct node *prev;          //previous value of scan


    //No stoppers allowed.
    if(datum <= STOPPER) datum = STOPPER + 1;
    
    //Add new datum into buffer 
    if( (++oldest - buffer) >= TAU) oldest=buffer;
    oldest->value = datum;
    
    //save pointer to old value's successor          
    successor = oldest->next;
    
    //Handle chain-out of first item in chain as special case
    scan = &head;
    if( scan->next == oldest ) scan->next = successor;
    
    //Scan through linked list from largest to smallest value
    do
    {
        prev = scan;
        scan = scan->next;
        
        //Remove old datapoint by linking its neighbor nodes
        if( scan->next == oldest )
            scan->next = successor;
        
        //Add new datum here if it fits in sorted list
        if( scan->value < datum )
        {
            oldest->next = prev->next;
            prev->next = oldest;
            datum = STOPPER;  //prevent adding datum in future
        }
        
    } while( scan != &tail );
    
    
    //find median value
    scan = &head;
    for(int i=0; i<=int(TAU/2); ++i)
        scan = scan->next;

    return( scan->value );
}









