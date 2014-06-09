//medfilt_new.c    Median filter
//http://www.embedded.com/2000/0011/0011feat3.htm
//thoroughly modified by DST 032508
// NOTE: assumes odd TAU?

#include <stdio.h>
#include <math.h>
#include <assert.h>

//
#ifdef MEX_FCN
    #include "mex.h"
    #include "matrix.h"
    #define printf mexPrintf
    #define INF mxGetInf()
#endif

// float math.
/*
#define MIN(a,b) ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) )
#define SWAP(a,b)  do { double temp=a; a=b; b=temp; } while(0)
*/


// This implementation stores the current window (TAU in size) of the data
// as a sorted linked list of data elements. "buffer" is the linked list
// of data, *oldest is the least recently added node that is removed when
// a new element is added, *head is the most recently added element, and 
// tail simply marks the end of the list.
// 
// Traversing the linked-list at any time gives the SORTED list of values
// in the current window of the data. Selecting the middle one gives us
// the median (if TAU is odd).
//
// Because the buffer is circular, we can change the value of the oldest
// element to the new/to-insert element and it is automatically in exactly
// the place it needs to be without any additional work.
//
// The elements are stored in the buffer in the order they are added and
// are to be removed. The linked list defines the sorted order.
//
// This isn't very different from removing the "oldest" element, adding the
// new one, and running insertion sort on the list.

//Initialize circular linked list buffer
/*
    head->[________]->tail
            ^
            oldest datapoint
            
    Scan from head; add new datum in sorted location.
*/
struct node
{
    struct node *next;
    struct node *prev;
    double value;
};

struct node *buffer;        //datapoint buffer (window)
struct node  head={0,0,0};  //list head virtual node
struct node  tail={0,0,0};  //list head virtual node
struct node *oldest;        //p* to oldest datum to be removed next

int TAU;  //window size





double mfilt(double datum);
void printList();



double* medianfilter( double* data, const int nPoints, const int szWindow )
{
    int i;
    TAU = szWindow;
    int middle = int(TAU/2); //median position in window.
    
    //printf("Points=%d\n", nPoints );
    
    // Fill the buffer with the first data point in the input array. We do 
    // this because there are not yet TAU/2 points to find the median of.
	//node buffer[szWindow];
    buffer = new node[szWindow];
	oldest    = &buffer[0];
    
    head.next = &buffer[0];
    head.prev = 0;
    head.value = INF;

    for(i=0; i<TAU; ++i)
    {
        buffer[i].value = data[0];
        buffer[i].prev  = buffer+i-1;
        buffer[i].next  = buffer+i+1;
    }
    
    buffer[0].prev     = &head;
    buffer[TAU-1].next = &tail;
    
    tail.prev  = &buffer[TAU-1];
    tail.next  = 0;
    tail.value = -INF;
    
    
    // Finish filling the buffer before we can take anything out.
    for(i=0; i<middle; ++i)
       mfilt( data[i] );
    
    //printf("\nmiddle\n");
    
    // Filter, scanning the window over the data one element at a time.
    for(i=middle; i<nPoints; ++i)
       data[i-middle] = mfilt( data[i] );
    
    //printf("\nend\n");
    
    // For the last TAU datapoints, there is no more data to take the
    // median of, so we extend the last data point.
    for(i=nPoints; i<nPoints+middle; ++i)
       data[i-middle] = mfilt( data[nPoints-1] );
    
    //printf("\ndone\n");
    
    delete [] buffer;
    return data;
}



//median filter function.  returns median of last TAU datapoints sent
//we assume the circular buffer linked list is already initialized!
double mfilt(double datum)
{
    struct node *successor;   //p* to oldest node's successor (next node in sorted list)
    struct node *itr;         //iterater into the sorted list
    struct node *prev;        //previous value during iteration
    int i=0;
    
    //printf( "oldest=%ld (%ld)\n", oldest-buffer, long(oldest) );
    assert( (oldest-buffer) < TAU );
    
    
    //printf("List before removal:\n");
    //printList();
    
    // 1) Remove oldest datapoint from list.
    oldest->prev->next = oldest->next;
    oldest->next->prev = oldest->prev;
    //printf("\nRemoved %ld\n",oldest-buffer);
    
    //printf("List after removal:\n");
    //printList();
    
    
    // 2) Insert new datum at correct point. "oldest" now points to the
    // location where we will add the new node. This will iterate to the
    // tail and insert in front of it. Since the tail has a NULL next
    // pointer, it ends there.
    prev = &head;
    itr  = head.next;
    i=0;
    
    while( itr )
    {
        assert( i<=TAU );
        
        // If the new value fits here (sorted), insert it.
        // "oldest" will now point to the newly inserted node.
        if( itr->value < datum )
        {
            //printf("\nInsert %.1f (%d)\n",datum,i);
            prev->next    = oldest;
            itr->prev     = oldest;
            oldest->next  = itr;
            oldest->prev  = prev;
            oldest->value = datum;
            break;
        }
        
        // Iterate to next link.
        prev = itr;
        itr  = itr->next;
        ++i;
    }
    
    
    // Now that we've taken care of inserting the new node and removing the
    // old one, increment the pointer to the next oldest node in line.
    oldest = oldest+1;
    if( (oldest-buffer) >= TAU)  oldest=&buffer[0];  //wrap around circular buffer
    
    
    // Find and return the median value of the (now) sorted list.
    itr = head.next;
    
    for( i=0; i<int(TAU/2); ++i )
    {
        itr = itr->next;
        assert( itr>=buffer );
    }
    
    return itr->value;
}


void printList()
{
    node* itr  = head.next;
    int i = 0;
    while( itr )
    {
        printf("%ld-> %ld ->%ld  =  %.1f\n", (itr->prev)-buffer, itr-buffer, (itr->next)-buffer, itr->value  );
        itr  = itr->next;
        ++i;
        
        assert( i<12 );
    }
}



