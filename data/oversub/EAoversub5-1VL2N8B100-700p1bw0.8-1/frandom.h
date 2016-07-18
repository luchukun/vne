/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 */  

#include <iostream>
#include <math.h>
#ifndef frandom_H
#define frandom_H
//#include <time.h>
/* Period parameters */  
const int RD_N=624;
const int RD_M=397;
const int MATRIX_A=0x9908b0dfUL;   /* constant vector a */
const int UPPER_MASK=0x80000000UL; /* most significant w-r bits */
const int LOWER_MASK=0x7fffffffUL;/* least significant r bits */

static unsigned long mt[RD_N]; /* the array for the state vector  */
static int mti=RD_N+1; /* mti==RD_N+1 means mt[RD_N] is not initialized */

/* initializes mt[RD_N] with a seed */
void init_genrand(unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);
/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);
/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void);
/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);
/* generates a random number on [0,1)-real-interval */
double genrand_real2(void);
/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);
/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* generates a random integer 1+[0,n-1] from [0,1)*/
int unif_int(int n);
/* generates a random integer [x,y] */
int unif_int(int x,int y);
/* generates a random integer [x,y] */
float unif_float(float x,float y);
/* generate bip0=1-p,p1=p; from [0,1)*/
bool rand_b01(float p);
	
float exprnd(float mean);
#endif