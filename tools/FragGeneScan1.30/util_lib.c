#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

double log2(double a){
  return log(a)/log(2);
}


double **dmatrix(int num_row, int num_col){

  int i, j;
  double **m;

  m=(double **) malloc(num_row * sizeof(double*));
  if (!m) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in dmatrix()");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<num_row; i++) {
    m[i]=(double *) malloc(num_col * sizeof(double));
    if (!m[i]) {
      fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i, " in dmatrix()");
      exit(EXIT_FAILURE);
    }

    for(j=0; j<num_col; j++){
      m[i][j] = 0.0;
    }
  }
  return m;
}


int **imatrix(int num_row, int num_col){

  int i,j;
  int **m;

  m=(int **) malloc(num_row * sizeof(int*));
  if (!m) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in imatrix()");
    exit(EXIT_FAILURE);
  }
  
  for(i=0; i<num_row; i++) {
    m[i]=(int *) malloc(num_col * sizeof(int));
    if (!m[i]) {
      fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i ," in imatrix()");
      exit(EXIT_FAILURE);
    }

    for(j=0; j<num_col; j++){
      m[i][j] = 0;
    }
  }
  return m;
}


double *dvector(int nh){
  
  int j;
  double *v;  

  v=(double *)malloc(nh * sizeof(double));

  if (!v) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure in dvector()");
    exit(EXIT_FAILURE);
  }

  for(j=0; j<nh; j++){
    v[j] = 0.0;
  }
  return v;
}


int *ivector(int nh){
  
  int j;
  int *v;

  v=(int *)malloc(nh * sizeof(int));

  if (!v) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure in ivector()");
    exit(EXIT_FAILURE);
  }

  for(j=0; j<nh; j++){
    v[j] = 0;
  }
  return v;
}


void free_dvector(double *v){
  
  free(v);
}


void free_ivector(int *v){
  
  free(v);
}


void free_dmatrix(double **m, int num_row){

  int i;

  for(i=num_row-1; i>=0; i--) 
    free(m[i]);

  free(m);
}


void free_imatrix(int **m,int num_row){  

  int i;

  for(i=num_row-1; i>=0; i--) 
    free(m[i]);

  free(m);
}


int tr2int (char *tr){

  int result;

  if      (strcmp(tr, "MM")==0){   result = 0; }
  else if (strcmp(tr, "MI")==0){   result = 1; }
  else if (strcmp(tr, "MD")==0){   result = 2; }
  else if (strcmp(tr, "II")==0){   result = 3; }
  else if (strcmp(tr, "IM")==0){   result = 4; }
  else if (strcmp(tr, "DD")==0){   result = 5; }
  else if (strcmp(tr, "DM")==0){   result = 6; }
  else if (strcmp(tr, "GE")==0){   result = 7; }
  else if (strcmp(tr, "GG")==0){   result = 8; }
  else if (strcmp(tr, "ER")==0){   result = 9; }
  else if (strcmp(tr, "RS")==0){   result = 10;}
  else if (strcmp(tr, "RR")==0){   result = 11;}
  else if (strcmp(tr, "ES")==0){   result = 12;}    /* ES: E+ -> S+, E- -> S- */
  else if (strcmp(tr, "ES1")==0){   result = 13;}   /* ES1: E+ -> S-, E- -> S+ */
  
  return result;
}


int nt2int (char nt){

  int result;

  if      (nt == 'A' || nt == 'a'){  result = 0; }
  else if (nt == 'C' || nt == 'c'){  result = 1; }
  else if (nt == 'G' || nt == 'g'){  result = 2; }
  else if (nt == 'T' || nt == 't'){  result = 3; }
  else                            {  result = 4; }

  return result;
}


int nt2int_rc (char nt){

  int result;

  if      (nt == 'A' || nt == 'a'){  result = 3; }
  else if (nt == 'C' || nt == 'c'){  result = 2; }
  else if (nt == 'G' || nt == 'g'){  result = 1; }
  else if (nt == 'T' || nt == 't'){  result = 0; }
  else                            {  result = 4; }

  return result;
}

int nt2int_rc_indel (char nt){

  int result;

  if      (nt == 'A' ){  result = 3; }
  else if (nt == 'C' ){  result = 2; }
  else if (nt == 'G' ){  result = 1; }
  else if (nt == 'T' ){  result = 0; }
  else if (nt == 'a' ){  result = 8; }
  else if (nt == 'c' ){  result = 7; }
  else if (nt == 'g' ){  result = 6; }
  else if (nt == 't' ){  result = 5; }
  else if (nt == 'n' ){  result = 9; }
  else if (nt == 'x' ){  result = 10;}
  else                {  result = 4; }

  return result;
}


int trinucleotide (char a, char b, char c){

  int freq_id;

  if      (a == 'A' || a == 'a'){  freq_id = 0;}
  else if (a == 'C' || a == 'c'){  freq_id = 16;}
  else if (a == 'G' || a == 'g'){  freq_id = 32;}
  else if (a == 'T' || a == 't'){  freq_id = 48;}
  else { freq_id = 0;}

  if      (b == 'A' || b == 'a'){  freq_id += 0;}
  else if (b == 'C' || b == 'c'){  freq_id += 4;}
  else if (b == 'G' || b == 'g'){  freq_id += 8;}
  else if (b == 'T' || b == 't'){  freq_id += 12;}
  else {freq_id = 0;}

  if      (c == 'A' || c == 'a'){  freq_id += 0;}
  else if (c == 'C' || c == 'c'){  freq_id += 1;}
  else if (c == 'G' || c == 'g'){  freq_id += 2;}
  else if (c == 'T' || c == 't'){  freq_id += 3;}
  else {freq_id = 0;}

  return freq_id;
}

int trinucleotide_pep (char a, char b, char c){

  int freq_id;

  if      (a == 'A' || a == 'a'){  freq_id = 0;}
  else if (a == 'C' || a == 'c'){  freq_id = 16;}
  else if (a == 'G' || a == 'g'){  freq_id = 32;}
  else if (a == 'T' || a == 't'){  freq_id = 48;}
  else { freq_id = 64;}

  if (freq_id <64){
    if      (b == 'A' || b == 'a'){  freq_id += 0;}
    else if (b == 'C' || b == 'c'){  freq_id += 4;}
    else if (b == 'G' || b == 'g'){  freq_id += 8;}
    else if (b == 'T' || b == 't'){  freq_id += 12;}
    else {freq_id = 64;}
  }

  if (freq_id < 64){
    if      (c == 'A' || c == 'a'){  freq_id += 0;}
    else if (c == 'C' || c == 'c'){  freq_id += 1;}
    else if (c == 'G' || c == 'g'){  freq_id += 2;}
    else if (c == 'T' || c == 't'){  freq_id += 3;}
    else {freq_id = 64;}
  }
  return freq_id;
}

void get_rc_dna(char *dna, char *dna1){

  char codon[5] = {'A', 'C', 'G', 'T', 'N'};
  int i;
  int dna_len = strlen(dna);
  for (i=0; i<dna_len; i++){
    
    dna1[dna_len-i-1] =codon[nt2int_rc(dna[i])];
  }
}

void get_rc_dna_indel(char *dna, char *dna1){

  char codon[11] = {'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', 'x'};
  int i;
  int dna_len = strlen(dna);
  for (i=0; i<dna_len; i++){
    
    dna1[dna_len-i-1] =codon[nt2int_rc_indel(dna[i])];
  }
}


void get_protein(char *dna, char *protein,  int strand, int whole_genome){

  int i;
  char codon_code[65] = {'K','N','K','N',
			 'T','T','T','T',
			 'R','S','R','S',
			 'I','I','M','I',
			 'Q','H','Q','H',
			 'P','P','P','P',
			 'R','R','R','R',
			 'L','L','L','L',
			 'E','D','E','D',
			 'A','A','A','A',
			 'G','G','G','G',
			 'V','V','V','V',
			 '*','Y','*','Y',
			 'S','S','S','S',
			 '*','C','W','C',
			 'L','F','L','F', 'X'};

  char anti_codon_code[65] = {'F','V','L','I',  
			 'C','G','R','S',
			 'S','A','P','T',
			 'Y','D','H','N',
			 'L','V','L','M',
			 'W','G','R','R',
			 'S','A','P','T',
			 '*','E','Q','K',
			 'F','V','L','I',
			 'C','G','R','S',
			 'S','A','P','T',
			 'Y','D','H','N',
			 'L','V','L','I',
			 '*','G','R','R',
			 'S','A','P','T',
			 '*','E','Q','K','X'};
  int dna_len = strlen(dna);

  if (strand ==1){
    for (i=0; i<dna_len; i+=3){
      protein[i/3] = codon_code[trinucleotide_pep(dna[i], dna[i+1], dna[i+2])];
    }
  }else{
    int protein_len = dna_len/3;
    if (dna_len % 3 == 2){
      dna_len -= 2;
    }else if (dna_len % 3 == 1){
      dna_len -= 1;
    }
    for (i=0; i<dna_len; i+=3){
      protein[(dna_len-i)/3-1] = anti_codon_code[trinucleotide_pep(dna[i], dna[i+1], dna[i+2])];
      protein_len --;
    }
  }

  if(protein[strlen(protein) - 1] == '*') { //remove the ending *
    protein[strlen(protein) - 1] = 0;
  }

  //alternative start codons still encode for Met
  //E. coli uses 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG and one or two others (e.g., an AUU and possibly a CUG)
  //only consider two major alternative ones, GTG and TTG
  if(whole_genome == 0) return; //short reads, skip
  if(strand == 1) {
  	int s = trinucleotide_pep(dna[0], dna[1], dna[2]); 
	if(s == trinucleotide_pep('G', 'T', 'G') || s == trinucleotide_pep('T', 'T', 'G')){
		protein[0] = 'M';
	}
  }
  else {
  	int s = trinucleotide_pep(dna[dna_len - 3], dna[dna_len - 2], dna[dna_len - 1]); 
	if(s == trinucleotide_pep('C', 'A', 'C') || s == trinucleotide_pep('C', 'A', 'A')) {
		protein[0] = 'M';
	}
  } 
}

void print_usage(){

  printf("%s", "USAGE: ./FragGeneScan.pl -s [seq_file_name] -o [output_file_name] -w [1 or 0] -t [train_file_name] (-p [thread_num])\n\n");
  printf("%s", "       Mandatory parameters\n");
  printf("%s", "       [seq_file_name]:    sequence file name including the full path\n");
  printf("%s", "       [output_file_name]: output file name including the full path\n");
  printf("%s", "       [1 or 0]:           1 if the sequence file has complete genomic sequences\n");
  printf("%s", "                           0 if the sequence file has short sequence reads\n");
  printf("%s", "       [train_file_name]:  file name that contains model parameters; this file should be in the \"train\" directory\n");
  printf("%s", "                           Note that four files containing model parameters already exist in the \"train\" directory\n");
  printf("%s", "                           [complete] for complete genomic sequences or short sequence reads without sequencing error\n");
  printf("%s", "                           [sanger_5] for Sanger sequencing reads with about 0.5% error rate\n");
  printf("%s", "                           [sanger_10] for Sanger sequencing reads with about 1% error rate\n");
  printf("%s", "                           [454_5] for 454 pyrosequencing reads with about 0.5% error rate\n");
  printf("%s", "                           [454_10] for 454 pyrosequencing reads with about 1% error rate\n");
  printf("%s", "                           [454_30] for 454 pyrosequencing reads with about 3% error rate\n");
  printf("%s", "                           [illumina_5] for Illumina sequencing reads with about 0.5% error rate\n");
  printf("%s", "                           [illumina_10] for Illumina sequencing reads with about 1% error rate\n\n");
  printf("%s", "       Optional parameter\n");
  printf("%s", "       [thread_num]:       the number of threads used by FragGeneScan; default is 1 thread.\n");
}
