// BinaryHeap class interface
//
// Etype: must have zero-parameter constructor and operator=;
//     must have operator<
// CONSTRUCTION: with (a) Etype representing negative infinity
// Copy construction of BinaryHeap objects is DISALLOWED
// Deep copy is supported
//
// ******************PUBLIC OPERATIONS******************
// void Inesrt( Etype X ) --> Insert X
// Etype FindMin( )       --> Return smallest item
// void DeleteMin( )      --> Remove smallest item
// void DeleteMin( Etype & X ) --> Same, but put it in X
// int IsEmpty( )         --> Return 1 if empty; else return 0
// int IsFull( )          --> Return 1 if full; else return 0
// void MakeEmpty( )      --> Remove all items
// void Toss( Etype X )   --> Insert X (lazily)
// void FixHeap( )        --> Reestablish heap order property
// ******************ERRORS*****************************
// Predefined exception is propogated if new fails
// EXCEPTION is called for FindMin or DeleteMin when empty

#ifndef __BinaryHeap
#define __BinaryHeap
#include <stdlib.h>
template <class Etype>
class BinaryHeap
{
  public:
        // Constructor, destructor, and copy assignment
    BinaryHeap( const Etype & MinVal );
    ~BinaryHeap( ) { delete [ ] Array; }

    const BinaryHeap & operator=( const BinaryHeap & Rhs );

        // Add an item maintaining heap order
    void Insert( const Etype & X );

        // Add an item but do not maintain order
    void Toss( const Etype & X );

        // Return minimum item in heap
    const Etype & FindMin( );

        // Delete minimum item in heap
    void DeleteMin( );
    void DeleteMin( Etype & X );

        // Reestablish heap order (linear-time function)
    void FixHeap( );

    int IsEmpty( ) const { return CurrentSize == 0; }
    int IsFull( ) const  { return 0; }
    void MakeEmpty( )    { CurrentSize = 0; }
  private:
    enum { DefaultSize = 10 };

    int MaxSize;         // Number of elements that can be stored
    int CurrentSize;     // Number of elements currently stored
    int OrderOK;         // Zero if heap order is not guaranteed
    Etype *Array;        // Dynamically allocated array

    BinaryHeap( const BinaryHeap & ); // Disable copy constructor

        // Internal routines
    void PercolateDown( int Index );
    void GetArray( int NewMaxSize );  // Allocate array
    void CheckSize( );   // Used for Toss and Insert
};
#endif


//#include "Exception.h"

// Routine to allocate the heap array

template <class Etype>
void
BinaryHeap<Etype>::GetArray( int NewMaxSize )
{
    Array = new Etype [ NewMaxSize + 1 ];
}

// Constructor for BinaryHeap

template <class Etype>
BinaryHeap<Etype>::BinaryHeap( const Etype & MinVal ) :
        MaxSize( DefaultSize ), CurrentSize( 0 ), OrderOK( 1 )
{
    GetArray( MaxSize );
    Array[ 0 ] = MinVal;
}

// If heap is full, double heap array

template <class Etype>
void
BinaryHeap<Etype>::CheckSize( )
{
    if( CurrentSize == MaxSize )
    {
        Etype *Old = Array;
        GetArray( MaxSize * 2 );
        for( int i = 0; i <= MaxSize; i++ )
            Array[ i ] = Old[ i ];
        delete [ ] Old;
        MaxSize *= 2;
    }
}

// Add X into the heap without maintaining order

template <class Etype>
void
BinaryHeap<Etype>::Toss( const Etype & X )
{
    CheckSize( );
    Array[ ++CurrentSize ] = X;
    if( X < Array[ CurrentSize / 2 ] )
        OrderOK = 0;
}

// Insert X into heap and if heap order is being maintained,
// percolate X up as needed

template <class Etype>
void
BinaryHeap<Etype>::Insert( const Etype & X )
{
    if( OrderOK == 0 )
    {
        Toss( X );
        return;
    }

    CheckSize( );

        // Percolate up
    int Hole = ++CurrentSize;
    for( ; X < Array[ Hole / 2 ]; Hole /= 2 )
        Array[ Hole ] = Array[ Hole / 2 ];
    Array[ Hole ] = X;
}

// Return minimum item in the heap
// Call FixHeap first if necessary

template <class Etype>
const Etype &
BinaryHeap<Etype>::FindMin( )
{
    //EXCEPTION( IsEmpty( ), "Binary heap is empty" );

    if( OrderOK == 0 )
        FixHeap( );
    return Array[ 1 ];
}

// Delete the minimum item and place it in X

template <class Etype>
void
BinaryHeap<Etype>::DeleteMin( Etype & X )
{
    X = FindMin( );
    Array[ 1 ] = Array[ CurrentSize-- ];
    PercolateDown( 1 );
}

// Delete the minimum item; throw it away
// NOTE: It would be better to write an additional
// private member to consolidate the common work between
// the two forms of DeleteMin.

template <class Etype>
void
BinaryHeap<Etype>::DeleteMin( )
{
    Etype X;
    DeleteMin( X );
}

// Private member to percolate down in the heap

template <class Etype>
void
BinaryHeap<Etype>::PercolateDown( int Hole )
{
    int Child;
    Etype Tmp = Array[ Hole ];

    for( ; Hole * 2 <= CurrentSize; Hole = Child )
    {
        Child = Hole * 2;
        if( Child != CurrentSize &&
                Array[ Child + 1 ] < Array[ Child ] )
            Child++;
        if( Array[ Child ] < Tmp )
            Array[ Hole ] = Array[ Child ];
        else
            break;
    }
    Array[ Hole ] = Tmp;
}

// Linear time FixHeap member

template <class Etype>
void
BinaryHeap<Etype>::FixHeap( )
{
    for( int i = CurrentSize / 2; i > 0; i-- )
        PercolateDown( i );
    OrderOK = 1;
}

