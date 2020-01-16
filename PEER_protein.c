#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include <time.h>
#define MAX_CARDIN 0x80000
//#ifdef DIST_PRINT
//#define NUM_COL 200
//#endif
//unsigned char dna[MAX];
unsigned char *dna;
unsigned int dna_len;
unsigned int cardinality;
static int dist[MAX_CARDIN][2];
//long int **dist;
double total_time;
long avg_seq_len;
//#ifdef DIST_PRINT
//static unsigned int dist_mat[MAX_CARDIN][10];//[NUM_COL];
//#endif
double fip1[MAX_CARDIN];// feature vector to store sequence 1
double fip2[MAX_CARDIN]; // feature vector to store sequence 2
double log2_table[MAX];
void parse( double * );
double euclidian(void);
//void entropy();
unsigned int lookup[26]={0,-1,1,2,3,4,5,6,7,-1,8,9,10,11,-1,12,13,14,16,-1,17,18,-1,19,-1};


int main (int argc, char *argv[]) 
{
   // declaring necessary variables to read raw dna file
   FILE *fp1, *fp2;
   int i=0; double ed = 0.0;
   long seconds, ns; 
   struct timespec start, finish;
   void gen_log2_table(void);
   // opening the data files

   fp1 = fopen(argv[1], "rb");
   if( fp1 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
   fp2 = fopen(argv[2], "rb");
   if( fp2 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
  
   gen_log2_table();
   dna = (unsigned char *)malloc(MAX);
   fread(dna, MAX , 1, fp1);   
   dna_len = strlen(dna)-1;
   avg_seq_len = dna_len;
   fflush(stdout);
   //dist = (long int **)malloc(MAX_CARDIN*sizeof(long int *));
    //if (dist == NULL) {
     //printf("Error in allocating dist \n");
    //return 0;
   //}

  // for(i=0; i<MAX_CARDIN; i++) { 
   //dist[i] = (long int *)malloc(9*sizeof(long int));
   //if (dist[i] == NULL) {
   //   printf("Error in allocating dist6[%llu] \n", i);
   //   return 0;
   //}
 //}
   parse(fip1);
   
   fread(dna, MAX , 1, fp2);   
   dna_len = strlen(dna)-1;
   avg_seq_len = avg_seq_len + dna_len;
   avg_seq_len = avg_seq_len/2; 
   fflush(stdout);
   parse(fip2);
   clock_gettime(CLOCK_REALTIME, &start);
   ed = euclidian();  
   clock_gettime(CLOCK_REALTIME, &finish);
   seconds = finish.tv_sec - start.tv_sec;
   ns = finish.tv_nsec - start.tv_nsec;
   total_time = total_time + (double)seconds + (double)ns/(double)100000000;
   printf("\n\nTotal time needed to compute ED between two wait vetors : %e \n and the Distance is %f and average sequence length %ld\n", total_time, ed,avg_seq_len);
   
   return(0);
}

void gen_log2_table(void)
{
  unsigned int i;

  for(i=0; i<MAX; i++)
    log2_table[i] = log((double)i)/log(2.0);
}

unsigned int GetIdx(unsigned int i, unsigned int kopt)
{
  int j;
  static int idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<kopt; j++)
    {
#ifdef DEBUG_PARSING
      printf("%c", dna[j]);
#endif
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    #ifdef DEBUG_PARSING
    printf(": ");
    printf("Idx = %d", idx);
    printf(": loc = %d\n",i);
#endif
    return idx;  
  }
  #ifdef DEBUG_PARSING
    for(j=0; j<kopt; j++)
      printf("%c", dna[i+j]);
#endif
    idx = (idx*4+lookup[dna[i+kopt-1]-'a'])%cardinality;
#ifdef DEBUG_PARSING
    printf(": ");
    printf("Idx = %d", idx);
    printf(": loc = %d\n",i);
#endif
    return idx;   
}


double my_log2(unsigned int d)
{
  return log2_table[d];
}

void parse( double * fip)
{
  unsigned int kopt, i, j, idx, d; 
  long seconds, ns; 
  struct timespec start, finish;
  int num_col =0;

  kopt = (unsigned int)(round(log((double)dna_len)/log(20.0)));
  cardinality = pow(4,kopt);
  
  
  printf("\nInside parse() : \n dna len = %d, kopt = %d\n", dna_len, kopt);
  for(i=0; i<cardinality; i++) dist[i][0] = -1;// initializing first location to -1
  
  clock_gettime(CLOCK_REALTIME, &start);


  for(i=0; i<(dna_len-kopt+1); i++)
  {
    idx = GetIdx(i, kopt);
#ifdef DIST_PRINT
    printf("\nIdx = %d: ", idx);
#endif
    // setting first location of the idx-th kopt-word  
    if (dist[idx][0] == -1)
    {
      dist[idx][0] = i;
      dist[idx][1] = i;
      //printf(" Initializing %d-th kopt-word with first location %d", idx,i);
    }
    else
    {
      d = i - dist[idx][1];
      //printf("Re-entering at %d-th kopt-word",idx);
      //printf(" Previous location was %d -- current location is %d",dist[idx][1], i);
      //printf(" interval %d  ", d);
//#ifdef DIST_PRINT
      //dist_mat[idx][dist[idx][2]] = d;
//#endif
      dist[idx][1] = i;
      fip[idx] = fip[idx] + (double)d * my_log2(d); // (double)(d)/dna_len  * log2( dna_len/(double)(d) );
      //dist[idx][2] ++; 
      //printf(" my_log2(d)  %e\n",  my_log2(d) ) ;
    }
  }

  for (i =0 ; i< cardinality; i++)
  {
    if ( dist[i][0] == -1 )
    ;
    else {
     d = dna_len - dist[i][1] + dist[i][0];
     //printf("\n\n\n\n Adjusting circular entropy:  \n");
     //printf("Re-entering at %d-th kopt-word",i);
     //printf(" Previous location was %d ",dist[i][1]);
     //printf(" interval %d  ", d);
//#ifdef DIST_PRINT
      //dist_mat[i][dist[i][2]] = d;
//#endif
     //dist[i][2] ++;
     fip[i] = fip[i] + (double)d * my_log2(d); // (double)(d)/dna_len  * log2( dna_len/(double)(d) );
     fip[i] = (double)my_log2(dna_len) - fip[i] / (double)dna_len;
     //printf(" my_log2(d) %e\n", my_log2(d) );
    } 
   //printf(" First localtion of %d-th kopt-word is %d\n",i,dist[i][0]); 
  }
 
  clock_gettime(CLOCK_REALTIME, &finish);
  seconds = finish.tv_sec - start.tv_sec;
  ns = finish.tv_nsec - start.tv_nsec;
  total_time = total_time + (double)seconds + (double)ns/(double)100000000;
  fflush(stdout);
  printf( "%e seconds\n", (double)seconds + (double)ns/(double)100000000);
  fflush(stdout);
  
  for (i=0;i<cardinality;i++)
  {
    printf("%f \n",fip[i]);
  }
  
 
}

double euclidian()

{
  int i=0; double distance =0.0, similarity =0.0, s=0.0, t=0.0;
  struct timespec start, finish;
  clock_gettime(CLOCK_REALTIME, &start);
  long seconds, ns;
  for (i =0; i < cardinality; i++) 
  {
    distance = distance +  (fip1[i] - fip2[i]) * (fip1[i] - fip2[i]); 
    s += fip1[i]*fip1[i];
    t += fip2[i]*fip2[i];
  }
  distance = sqrt(distance); s = sqrt(s); t = sqrt(t);
  if (s < t)
  similarity = (1 - distance/s) * 100;
  else
  similarity = (1 - distance/t) * 100;
  printf("%lf\n",similarity);
  clock_gettime(CLOCK_REALTIME, &finish);
  seconds = finish.tv_sec - start.tv_sec;
  ns = finish.tv_nsec - start.tv_nsec;
  total_time = total_time + (double)seconds + (double)ns/(double)100000000;
  fflush(stdout);
  //printf( "%e seconds\n", (double)seconds + (double)ns/(double)100000000);
  fflush(stdout);
  return(distance);
}


