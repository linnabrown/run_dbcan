#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double **dmatrix(int num_row, int num_col);
double *dvector(int nh);
int **imatrix(int num_row, int num_col);
int *ivector(int nh);

void free_dvector(double *v);
void free_dmatrix(double **m,int num_row);
void free_ivector(int *v);
void free_imatrix(int **m,int num_row);

int tr2int (char *nt);
int nt2int (char nt);
int nt2int_rc (char nt);


int trinucleotide (char a, char b, char c);
double log2(double a);
void get_protein(char *dna, char *protein, int strand, int whole_genome);
void print_usage();
  
