#ifndef BINARY_HEAPS_H_
#define BINARY_HEAPS_H_

#include <vector>
#include "elements.h"
using namespace std;
// BinaryHeap class
//
// CONSTRUCTION: with an optional capacity (that defaults to 100)
//
// ******************PUBLIC OPERATIONS*********************
// void insert( x )       --> Insert x
// deleteMin( minItem )   --> Remove (and optionally return) smallest item
// Type findMin( )  --> Return smallest item
// bool isEmpty( )        --> Return true if empty; else false
// bool isFull( )         --> Return true if full; else false
// void makeEmpty( )      --> Remove all items
// ******************ERRORS********************************
// Throws Underflow and Overflow as warranted

template <class Type>
class BinaryHeap
{
  public:
    explicit BinaryHeap( int capacity);
    bool isEmpty( ) const;
    bool isFull( ) const;
    const Type & findMin( ) const;

    void insert( const Type & x );
    Type deleteMin( );
    void deleteMin( Type & minItem );
    void makeEmpty( );
	
  private:
    int  currentSize;  // Number of elements in heap
	vector<Type> Array;        // The heap Array

    void buildHeap( );
    void percolateDown( int hole );

};
template <class Type>
BinaryHeap<Type>::BinaryHeap( int capacity )
{
	Array.resize(capacity + 1);
	currentSize=0;
}

/**
 * Insert item x into the priority queue, maintaining heap order.
 * Duplicates are allowed.
 * Throw Overflow if container is full.
 */
template <class Type>
void BinaryHeap<Type>::insert( const Type & x )
{
    if( isFull( ) )
	{
		std::cout<<"Resize the heap "<<Array.size()<<endl;

		Array.resize(2*Array.size()+1);
	}
        // Percolate up
    int hole = ++currentSize;
    for( ; hole > 1 && x < Array[ hole / 2 ]; hole /= 2 )
        Array[ hole ] = Array[ hole / 2 ];
    Array[ hole ] = x;
}

/**
 * Find the smallest item in the priority queue.
 * Return the smallest item, or throw Underflow if empty.
 */
template <class Type>
const Type & BinaryHeap<Type>::findMin( ) const
{
    if( isEmpty( ) )
       std::cout<<"The heap is empty"<<endl;
    return Array[ 1 ];
}

/**
 * Remove the smallest item from the priority queue.
 * Throw Underflow if empty.
 */
template <class Type>
Type BinaryHeap<Type>::deleteMin( )
{	
    if( isEmpty( ) )
        std::cout<<"The heap is empty"<<endl;
	Type minItem = Array[ 1 ];
    Array[ 1 ] = Array[ currentSize-- ];
    percolateDown( 1 );
	return minItem;
}

/**
 * Remove the smallest item from the priority queue
 * and place it in minItem. Throw Underflow if empty.
 */
template <class Type>
void BinaryHeap<Type>::deleteMin( Type & minItem )
{
    if( isEmpty( ) )
        std::cout<<"The heap is empty"<<endl;

    minItem = Array[ 1 ];
    Array[ 1 ] = Array[ currentSize-- ];
    percolateDown( 1 );
}

/**
 * Establish heap order property from an arbitrary
 * arrangement of items. Runs in linear time.
 */
template <class Type>
void BinaryHeap<Type>::buildHeap( )
{
    for( int i = currentSize / 2; i > 0; i-- )
        percolateDown( i );
}

/**
 * Test if the priority queue is logically empty.
 * Return true if empty, false otherwise.
 */
template <class Type>
bool BinaryHeap<Type>::isEmpty( ) const
{
    return currentSize == 0;
}

/**
 * Test if the priority queue is logically full.
 * Return true if full, false otherwise.
 */
template <class Type>
bool BinaryHeap<Type>::isFull( ) const
{
    return currentSize == Array.size( ) - 1;
}

/**
 * Make the priority queue logically empty.
 */
template <class Type>
void BinaryHeap<Type>::makeEmpty( )
{
    currentSize = 0;
}

/**
 * Internal method to percolate down in the heap.
 * hole is the index at which the percolate begins.
 */
template <class Type>
void BinaryHeap<Type>::percolateDown( int hole )
{
	int child;
	Type tmp = Array[ hole ];

	for( ; hole * 2 <= currentSize; hole = child )
	{
	  child = hole * 2;
	  if( child != currentSize && Array[ child + 1 ] < Array[ child ] )
		  child++;
	  if( Array[ child ] < tmp )
		 Array[ hole ] = Array[ child ];
	else
		  break;
	}
	Array[ hole ] = tmp;
}

#endif
