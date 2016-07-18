#include "binaryheaps.h"

/**
 * Construct the binary heap.
 * capacity is the capacity of the binary heap.
 */
template <class Type>
BinaryHeap<Type>::BinaryHeap( int capacity )
{
	Array=new Type[capacity + 1];
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
		std::cout<<"Resize the heap"<<endl;
		Array.resize(2*Array.max_size()+1);
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
void BinaryHeap<Type>::deleteMin( )
{
    if( isEmpty( ) )
        std::cout<<"The heap is empty"<<endl;

    Array[ 1 ] = Array[ currentSize-- ];
    percolateDown( 1 );
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
