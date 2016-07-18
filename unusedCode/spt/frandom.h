/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 */  

#include <iostream>
#include <math.h>
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
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<RD_N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (RD_N>key_length ? RD_N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=RD_N) { mt[0] = mt[RD_N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=RD_N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=RD_N) { mt[0] = mt[RD_N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= RD_N) { /* generate RD_N words at one time */
        int kk;

        if (mti == RD_N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<RD_N-RD_M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+RD_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<RD_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(RD_M-RD_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[RD_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[RD_N-1] = mt[RD_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* generates a random integer 1+[0,n-1] from [0,1)*/
int unif_int(int n)
{	
	int i;
	double f;
   /* Seed the random-number generator with current time so that
    * the numbers will be different every time we run.  */	
	//init_genrand( (unsigned long)time( NULL ));
	/* generates a random integer 1+[0,n-1] from [0,1) */
	f=1+n*genrand_real2();
	i=(int) (floor(f));
	return i;
}



/* generate bip0=1-p,p1=p; from [0,1)*/
bool rand_b01(float p)
{
	double f; 	
	//init_genrand( (unsigned long)time( NULL ));
	/* generates a random integer [0,n-1]+1 from [0,1) */
	f=genrand_real2();	
	if (f < (1-p))
		return 0;
	else
		return 1;
}	
	
float exprnd(float p)
{
	double f;
	f=genrand_real3();
	f=1/(1-f);
	return (float) log(f)/p;
}

