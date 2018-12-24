#include <stdio.h>
#include <math.h>
#include <string.h>
#include <values.h>
#include <stdlib.h>
#include <malloc.h>
#include <ctype.h>
#include "hmm.h"
#include "util_lib.h"

void dump_memory(void *p, int size);

void viterbi(HMM *hmm_ptr, TRAIN *train_ptr, char *O, FILE *fp_out, FILE *fp_aa, FILE *fp_dna,
char *head, int whole_genome, int cg, int format){
  double max_dbl = 10000000000.0;
  int debug=0;

  int *vpath;                          /* optimal path after backtracking */
  double **alpha;                      /* viterbi prob array */
  int **path;                          /* viterbi path array */
  int i,j,l,m,n, x,y,z,t,jj,kk;
  int temp_t;
  int ag_num;

  int orf_start=0;
  int orf_strand=0;
  int print_save;
  int frame_save;
  int orf_save;
  double prob_save, start_freq;

  int best;
  int from, from0, to;   /*from0: i-2 position, from: i-1 position */
  int from2;             /* from2: i-2, i-1 for condition in emission
probability */
  double temp_alpha, temp, prob, prob2;
  int len_seq;
  int gene_len;
  int count_ag;

  int temp_strand;
  int num_d;          /* the number of delete */
  int freq_id;
  double h_kd, r_kd, p_kd;

  int codon_start;
  char dna_tmp[300000];
  char dna[300000];
  char dna1[300000];
  char dna_f[300000];
  char dna_f1[300000];
  char protein[100000];
  int dna_id=0;
  int dna_f_id=0;
  int out_nt;
  int start_t, dna_start_t;
  int end_t;
  int prev_match;
  int start_orf;
  int frame;
  double final_score;

  int insert[100];
  int delete[100];
  int insert_id, delete_id;
  char *head_short=NULL;
  char delimi[] = " ";

  int temp_i[6] = {0,0,0,0,0,0};
  int temp_i_1[6] = {0,0,0,0,0,0};

  int num_N=0;
  /***************************************************************/
  /* initialize                                                  */
  /***************************************************************/
  int refine = 0;
  if (whole_genome==1){
    gene_len = 120;
    refine = 1;
  }else{
    gene_len = 60;
  }

  len_seq = strlen(O);
  alpha = (double **)dmatrix(hmm_ptr->N, len_seq);
  path = (int **)imatrix(hmm_ptr->N, len_seq);
  vpath = (int *)ivector(len_seq);

  for (i=0; i<hmm_ptr->N; i++){
    alpha[i][0] = -1 * hmm_ptr->pi[i];
  }

  double log53 = log(0.53);
  double log16 = log(0.16);
  double log30 = log(0.30);
  double log25 = log(0.25);
  double log95 = log(0.95);
  double log54 = log(0.54);
  double log83 = log(0.83);
  double log07 = log(0.07);
  /* stop state */
  if ((O[0] == 'T'|| O[0] == 't')  &&
      (((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'A'|| O[2] == 'a')) ||
       ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'G'|| O[2] == 'g')) ||
       ((O[1] == 'G'|| O[1] == 'g') && (O[2] == 'A'|| O[2] == 'a')))) {

    alpha[E_STATE][0] = max_dbl;
    alpha[E_STATE][1] = max_dbl;
    path[E_STATE][1] = E_STATE;
    path[E_STATE][2] = E_STATE;

    alpha[M6_STATE][2] = max_dbl;
    alpha[M5_STATE][1] = max_dbl;
    alpha[M4_STATE][0] = max_dbl;
    alpha[M3_STATE][2] = max_dbl;
    alpha[M2_STATE][1] = max_dbl;
    alpha[M1_STATE][0] = max_dbl;

    if ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'A'|| O[2] == 'a')){
      alpha[E_STATE][2] = alpha[E_STATE][2] - log53;
    }else if ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'G'|| O[2] == 'g')){
      alpha[E_STATE][2] = alpha[E_STATE][2] - log16;
    }else if((O[1] == 'G'|| O[1] == 'g') && (O[2] == 'A'|| O[2] == 'a')){
      alpha[E_STATE][2] = alpha[E_STATE][2] - log30;
    }
  }

  if ((O[2] == 'A'|| O[0] == 'a')  &&
      (((O[0] == 'T'|| O[0] == 't') && (O[1] == 'T'|| O[1] == 't')) ||
       ((O[0] == 'C'|| O[0] == 'c') && (O[1] == 'T'|| O[1] == 't')) ||
       ((O[0] == 'T'|| O[0] == 't') && (O[1] == 'C'|| O[1] == 'c')))) {

    alpha[S_STATE_1][0] = max_dbl;
    alpha[S_STATE_1][1] = max_dbl;
    alpha[S_STATE_1][2] = alpha[S_STATE][0];
    path[S_STATE_1][1] = S_STATE_1;
    path[S_STATE_1][2] = S_STATE_1;

    alpha[M3_STATE_1][2] = max_dbl;
    alpha[M6_STATE_1][2] = max_dbl;

    if ((O[0] == 'T'|| O[0] == 't') && (O[1] == 'T'|| O[1] == 't')){
      alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log53;
    }else if ((O[0] == 'C'|| O[0] == 'c') && (O[1] == 'T'|| O[1] == 't')){
      alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log16;
    }else if((O[0] == 'T'|| O[0] == 't') && (O[1] == 'C'|| O[1] == 'c')){
      alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log30;
    }
  }

  /******************************************************************/
  /*  fill out the rest of the columns                              */
  /******************************************************************/
  for (t = 1; t < len_seq; t++) {
    from = nt2int(O[t-1]);
    if (t>1){
      from0 = nt2int(O[t-2]);
    }else{
      from0 = 2;
    }
    to = nt2int(O[t]);

    /* if DNA is other than ACGT, do it later */
    if (from==4){ from=2; }
    if (from0==4){ from0=2; }
    if (to==4){
      to=2;
      num_N += 1;
    }else{
      num_N=0;
    }
    from2 = from0*4+from;

    /******************/
    /* M state        */
    /******************/

    for (i=M1_STATE; i<=M6_STATE; i++)   {

      if (alpha[i][t]<max_dbl){

	if (t==0){

	}else{

	  if (i==M1_STATE){

	    /* from M state */
	    j = M6_STATE;
	    alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
	    path[i][t] = j;

	    /* from D state */
	    if (whole_genome==0){
	      for (j=M5_STATE; j>=M1_STATE; j--){
		if (j >= i ){
		  num_d = i-j+6;
		}else if (j+1<i){
		  num_d = i-j;
		}else{
		  num_d = -10;
		}
		if(num_d>0){
		  temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to]
		    - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
		  if ( temp_alpha < alpha[i][t]){
		    alpha[i][t] = temp_alpha;
		    path[i][t] = j;
		  }
		}
	      }
	    }

	    /* from Start state */
	    temp_alpha = alpha[S_STATE][t-1] - hmm_ptr->e_M[0][from2][to];
	    if ( temp_alpha < alpha[i][t] ){
	      alpha[i][t] = temp_alpha;
	      path[i][t] = S_STATE;
	    }

	  }else{   /*i ==M2-M6*/

	    /* from M state */
	    j = i - 1;
	    alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i-M1_STATE][from2][to];
	    path[i][t] = j;


	    /* from D state */
	    if (whole_genome==0){
	      for (j=M6_STATE; j>=M1_STATE; j--){
		if (j >= i ){
		  num_d = i-j+6;
		}else if (j+1 < i){
		  num_d = i-j;
		}else{
		  num_d = -10;
		}
		if (num_d > 0){


		  temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i-M1_STATE][from2][to]
		    - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
		  if ( temp_alpha < alpha[i][t]){
		    alpha[i][t] = temp_alpha;
		    path[i][t] = j;
		  }
		}
	      }
	    }
	  }

	  /* from I state */
	  if (i==M1_STATE) { j = I6_STATE;
	  }else{ j = I1_STATE + (i - M1_STATE -1); }


	  /* to aviod stop codon */
          if (t<2){
          }else if((i==M2_STATE || i==M5_STATE) && (O[temp_i[j-I1_STATE]] == 'T'||O[temp_i[j-I1_STATE]] =='t') &&
		   (((O[t] == 'A'||O[t] == 'a') && (O[t+1] =='A'||O[t+1] =='a')) ||
		    ((O[t] == 'A'||O[t] == 'a') && (O[t+1] =='G'||O[t+1] =='g')) ||
		    ((O[t] == 'G'||O[t] == 'g') && (O[t+1] =='A'||O[t+1] =='a')))){

	  }else if ((i==M3_STATE || i==M6_STATE) && (O[temp_i[j-I1_STATE]-1] == 'T'||O[temp_i[j-I1_STATE]-1] =='t') &&
		    (((O[temp_i[j-I1_STATE]] == 'A'||O[temp_i[j-I1_STATE]] == 'a') && (O[t] =='A'||O[t] == 'a')) ||
		     ((O[temp_i[j-I1_STATE]] == 'A'||O[temp_i[j-I1_STATE]] == 'a') && (O[t] =='G'||O[t] == 'g')) ||
		     ((O[temp_i[j-I1_STATE]] == 'G'||O[temp_i[j-I1_STATE]] == 'g') && (O[t] =='A'||O[t] == 'a')))){
	  }else{
	    temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_IM] - log25;
	    if ( temp_alpha < alpha[i][t]){
	      alpha[i][t] = temp_alpha;
	      path[i][t] = j;
	    }
	  }
	}
      }
    }

    /******************/
    /* I state        */
    /******************/
    for (i=I1_STATE; i<=I6_STATE; i++) {

      if (t==0){
      }else{

	/* from I state */
	j = i;
	alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
	path[i][t] = j;

	/* from M state */
	j = i - I1_STATE + M1_STATE ;
	if (i==I6_STATE){
	  temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
	}else{
	  temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
	}
	if (temp_alpha < alpha[i][t]){
	  alpha[i][t] = temp_alpha;
	  path[i][t] = j;

	  temp_i[i-I1_STATE] = t-1;
	}
      }
    }

    /******************/
    /* M' state        */
    /******************/

    for (i=M1_STATE_1; i<=M6_STATE_1; i++)   {
      if  ((i==M1_STATE_1 || i==M4_STATE_1)&& t>=3 &&
	   (((O[t-3] == 'T'||O[t-3] == 't') && (O[t-2] == 'T'||O[t-2] == 't') && (O[t-1] == 'A'||O[t-1] =='a')) ||
	    ((O[t-3] == 'C'||O[t-3] == 'c') && (O[t-2] == 'T'||O[t-2] == 't') && (O[t-1] == 'A'||O[t-1] =='a')) ||
	    ((O[t-3] == 'T'||O[t-3] == 't') && (O[t-2] == 'C'||O[t-2] == 'c') && (O[t-1] == 'A'||O[t-1] =='a')))){

	/* from Start state  since this is actually stop codon in minus strand */
	alpha[i][t] = alpha[S_STATE_1][t-1] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to];
	path[i][t] = S_STATE_1;

      }else{

	if (t==0){
	}else{

	  if (i==M1_STATE_1 ){

	    /* from M state */
	    j = M6_STATE_1;
	    alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[0][from2][to];
	    path[i][t] = j;

	    /* from D state */
	    if (whole_genome==0){
	      for (j=M5_STATE_1; j>=M1_STATE_1; j--){
		if (j >= i){
		  num_d = i-j+6;
		}else if (j+1 <i){
		  num_d = i-j;
		}else{
		  num_d = -10;
		}
		if (num_d > 0){
		  temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[0][from2][to]
		    - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
		  if ( temp_alpha < alpha[i][t]){
		    alpha[i][t] = temp_alpha;
		    path[i][t] = j;
		  }
		}
	      }
	    }

	  }else{

	    /* from M state */
	    j = i - 1;
	    alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to];
	    path[i][t] = j;

	    /* from D state */
	    if (whole_genome==0){
	      for (j=M6_STATE_1; j>=M1_STATE_1; j--){
		if (j >= i ){
		  num_d = i-j+6;
		}else if (j+1 < i){
		  num_d = i-j;
		}else{
		  num_d = -10;
		}
		if (num_d>0){
		  temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to]
		    - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
		  if ( temp_alpha < alpha[i][t]){
		    alpha[i][t] = temp_alpha;
		    path[i][t] = j;
		  }
		}
	      }
	    }
	  }

	  /* from I state */
	  if (i==M1_STATE_1) { j = I6_STATE_1;
	  }else{ j = I1_STATE_1 + (i - M1_STATE_1 -1); }


	  /* to aviod stop codon */
	  	  /* to aviod stop codon */
          if (t<2){
          }else  if((i==M2_STATE_1 || i==M5_STATE_1) && (O[t+1] == 'A'||O[t+1] == 'a') &&
		    (((O[temp_i_1[j-I1_STATE_1]] == 'T'|| O[temp_i_1[j-I1_STATE_1]] == 't') && (O[t] =='T'|| O[t] =='t')) ||
		     ((O[temp_i_1[j-I1_STATE_1]] == 'C'|| O[temp_i_1[j-I1_STATE_1]] == 'c') && (O[t] =='T'|| O[t] =='t')) ||
		     ((O[temp_i_1[j-I1_STATE_1]] == 'T'|| O[temp_i_1[j-I1_STATE_1]] == 't') && (O[t] =='C'|| O[t] =='c')))){

	  }else if ((i==M3_STATE_1 || i==M6_STATE_1) && (O[t] == 'A'||O[t] == 'a') &&
		    (((O[temp_i_1[j-I1_STATE_1]-1] == 'T'|| O[temp_i_1[j-I1_STATE_1]-1]=='t') &&
		      (O[temp_i_1[j-I1_STATE_1]] =='T'|| O[temp_i_1[j-I1_STATE_1]] =='t')) || 
		     ((O[temp_i_1[j-I1_STATE_1]-1] == 'C'|| O[temp_i_1[j-I1_STATE_1]-1]=='c') &&
		      (O[temp_i_1[j-I1_STATE_1]] =='T'|| O[temp_i_1[j-I1_STATE_1]] =='t')) || 
		     ((O[temp_i_1[j-I1_STATE_1]-1] == 'T'|| O[temp_i_1[j-I1_STATE_1]-1]=='t') &&
		      (O[temp_i_1[j-I1_STATE_1]] =='C'|| O[temp_i_1[j-I1_STATE_1]] =='c')))){
	  }else{

	    temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_IM] - log25;
	    if ( temp_alpha < alpha[i][t]){
	      alpha[i][t] = temp_alpha;
	      path[i][t] = j;
	    }
	  }
	}
      }
    }

    /******************/
    /* I' state        */
    /******************/
    for (i=I1_STATE_1; i<=I6_STATE_1; i++) {

      if (t==0){
      }else{
	/* from I state */
	j = i;
	alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
	path[i][t] = j;

	/* from M state */
	if (path[S_STATE_1][t-3]!= R_STATE && path[S_STATE_1][t-4] !=R_STATE && path[S_STATE_1][t-5] !=R_STATE){
	  j = i - I1_STATE_1 + M1_STATE_1;
	  if (i==I6_STATE_1){
	    temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
	  }else{
	    temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
	  }
	  if (temp_alpha < alpha[i][t]){
	    alpha[i][t] = temp_alpha;
	    path[i][t] = j;

	    temp_i_1[i-I1_STATE_1] = t-1;
	  }
	}
      }
    }

    /***********************/
    /* Non_coding state    */
    /***********************/

    if (t==0){
    }else{
      alpha[R_STATE][t] = alpha[R_STATE][t-1] - hmm_ptr->tr_R_R[from][to] - hmm_ptr->tr[TR_RR];
      path[R_STATE][t] = R_STATE;

      temp_alpha = alpha[E_STATE][t-1]  - hmm_ptr->tr[TR_ER];
      if (temp_alpha < alpha[R_STATE][t] ){
	alpha[R_STATE][t] = temp_alpha;
	path[R_STATE][t] = E_STATE;
      }

      temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ER] ;
      if (temp_alpha < alpha[R_STATE][t] ){
        alpha[R_STATE][t] = temp_alpha;
        path[R_STATE][t] = E_STATE_1;
      }
      alpha[R_STATE][t] -= log95;
    }

    /******************/
    /* END state      */
    /******************/
    if (alpha[E_STATE][t] == 0){

      alpha[E_STATE][t] = max_dbl;
      path[E_STATE][t] = NOSTATE;

      if (t < len_seq -2 && (O[t] == 'T'||O[t] == 't')  &&
	(((O[t+1] == 'A'||O[t+1] == 'a') && (O[t+2] == 'A'||O[t+2] =='a')) ||
	 ((O[t+1] == 'A'||O[t+1] == 'a') && (O[t+2] == 'G'||O[t+2] =='g')) ||
	 ((O[t+1] == 'G'||O[t+1] == 'g') && (O[t+2] == 'A'||O[t+2] =='a')))) {

	alpha[E_STATE][t+2] = max_dbl;
	/* transition from frame4,frame5,and frame6 */
	temp_alpha = alpha[M6_STATE][t-1] - hmm_ptr->tr[TR_GE];
	if (temp_alpha < alpha[E_STATE][t+2]){
	  alpha[E_STATE][t+2] = temp_alpha;
	  path[E_STATE][t] = M6_STATE;
	}

	/* transition from frame1,frame2,and frame3 */
	temp_alpha  = alpha[M3_STATE][t-1] - hmm_ptr->tr[TR_GE];
	if (temp_alpha < alpha[E_STATE][t+2]){
	  alpha[E_STATE][t+2] = temp_alpha;
	  path[E_STATE][t] = M3_STATE;
	}

	alpha[E_STATE][t] = max_dbl;
	alpha[E_STATE][t+1] = max_dbl;
	path[E_STATE][t+1] = E_STATE;
	path[E_STATE][t+2] = E_STATE;

	alpha[M6_STATE][t+2] = max_dbl;
	alpha[M5_STATE][t+1] = max_dbl;
	alpha[M4_STATE][t] = max_dbl;
	alpha[M3_STATE][t+2] = max_dbl;
	alpha[M2_STATE][t+1] = max_dbl;
	alpha[M1_STATE][t] = max_dbl;

	if ((O[t+1] == 'A'||O[t+1] =='a') && (O[t+2] == 'A'||O[t+2] =='a')){
	  alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log54;
	}else if ((O[t+1] == 'A'||O[t+1] =='a') && (O[t+2] == 'G'||O[t+2] =='g')){
	  alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log16;
	}else if((O[t+1] == 'G'||O[t+1] == 'g') && (O[t+2] == 'A'||O[t+2] =='a')) {
	  alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log30;
	}

	/* adjustment based on probability distribution */
	start_freq=0;
	freq_id = 0;

        double sub_sum = 0;
        int sub_count = 0;

        if (t>=60){ /* bug reported by Yu-Wei */
                for(i=-60; i<=-3; i++){
			if (t+i+2 < len_seq)
			{
                        	start_freq -= hmm_ptr->tr_E[i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
			}
                }
        }else{
                for(i=(-1*t); i<=-3; i++){
			if (t+i+2 < len_seq)
			{
                        	sub_sum += hmm_ptr->tr_E[i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
			}
                }
                sub_sum = sub_sum * 58 / (-3 + t + 1);
                start_freq -= sub_sum;
        }

	h_kd = hmm_ptr->E_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E_dist[1],2)/(2*pow(hmm_ptr->E_dist[0],2)));
	r_kd = hmm_ptr->E_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E_dist[4],2)/(2*pow(hmm_ptr->E_dist[3],2)));
	p_kd = h_kd / (h_kd + r_kd);
	if (p_kd<0.01){
	  p_kd=0.01;
	}else if (p_kd>0.99){
	  p_kd=0.99;
	}
	alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log(p_kd);
      }
    }

    /*************************************************/
    /* START' state                                  */
    /* origianlly stop codon of genes in - strand    */
    /*************************************************/
    if (alpha[S_STATE_1][t] == 0){

      alpha[S_STATE_1][t] = max_dbl;
      path[S_STATE_1][t] = NOSTATE;


      if (t<len_seq-2 && (O[t+2] == 'A'||O[t+2] == 'a') &&
	  (((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'T'|| O[t+1] == 't')) ||
	   ((O[t] == 'C'||O[t] =='c') && (O[t+1] == 'T'|| O[t+1] == 't')) ||
	   ((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'C'|| O[t+1] == 'c')))) {

	alpha[S_STATE_1][t] = max_dbl;
	path[S_STATE_1][t] = R_STATE;
	alpha[S_STATE_1][t+1] = max_dbl;
	alpha[S_STATE_1][t+2] = alpha[R_STATE][t-1] - hmm_ptr->tr[TR_RS];
	path[S_STATE_1][t+1] = S_STATE_1;
	path[S_STATE_1][t+2] = S_STATE_1;

	temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ES];
	if (temp_alpha < alpha[S_STATE_1][t+2]){
	  alpha[S_STATE_1][t+2] = temp_alpha;
	  path[S_STATE_1][t] = E_STATE_1;
	}

	temp_alpha = alpha[E_STATE][t-1] - hmm_ptr->tr[TR_ES1];
	if (temp_alpha < alpha[S_STATE_1][t+2]){
	  alpha[S_STATE_1][t+2] = temp_alpha;
	  path[S_STATE_1][t] = E_STATE;
	}

        alpha[M3_STATE_1][t+2] = max_dbl;
        alpha[M6_STATE_1][t+2] = max_dbl;

	if ((O[t] == 'T'||O[t] == 't') && (O[t+1] == 'T'||O[t+1] == 't')){
	  alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log54;
	}else if ((O[t] == 'C'||O[t] =='c') && (O[t+1] == 'T'||O[t+1]=='t')){
	  alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log16;
	}else if((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'C'||O[t+1] =='c')){
	  alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log30;
	}

	/* adjustment based on probability distribution */
	start_freq=0;
	freq_id = 0;
	for(i=3; i<=60; i++){
	  if (t+i+2 < len_seq)
	  {
	    start_freq -= hmm_ptr->tr_S_1[i-3][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
	  }
	}
	h_kd = hmm_ptr->S1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[1],2)/(2*pow(hmm_ptr->S1_dist[0],2)));
	r_kd = hmm_ptr->S1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[4],2)/(2*pow(hmm_ptr->S1_dist[3],2)));
	p_kd = h_kd / (h_kd + r_kd);
	if (p_kd<0.01){
	  p_kd=0.01;
	}else if (p_kd>0.99){
	  p_kd=0.99;
	}
	alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log(p_kd);
      }
    }

    /************************/
    /* START state          */
    /************************/
    if (alpha[S_STATE][t] == 0){

      alpha[S_STATE][t] = max_dbl;
      path[S_STATE][t] = NOSTATE;

      if (t<len_seq-2 &&  (O[t+1] == 'T'||O[t+1] =='t') && (O[t+2] == 'G'||O[t+2] =='g')&&
	  ((O[t] == 'A'||O[t] =='a') || (O[t] == 'G'||O[t] =='g') ||  (O[t] == 'T'||O[t] =='t'))) {

	alpha[S_STATE][t] = max_dbl;
	alpha[S_STATE][t+1] = max_dbl;
	alpha[S_STATE][t+2] = alpha[R_STATE][t-1] - hmm_ptr->tr[TR_RS];
	path[S_STATE][t] = R_STATE;
	path[S_STATE][t+1] = S_STATE;
	path[S_STATE][t+2] = S_STATE;

	temp_alpha = alpha[E_STATE][t-1] - hmm_ptr->tr[TR_ES];
	if (temp_alpha < alpha[S_STATE][t+2]){
	  alpha[S_STATE][t+2] = temp_alpha;
	  path[S_STATE][t] = E_STATE;
	}

	temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ES1];
	if (temp_alpha < alpha[S_STATE][t+2]){
	  alpha[S_STATE][t+2] = temp_alpha;
	  path[S_STATE][t] = E_STATE_1;
	}


	if ((O[t] == 'A'||O[t] =='a') ){
	  alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log83;
	}else if ((O[t] == 'G'||O[t] =='g')){
	  alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(0.10);
	}else if((O[t] == 'T'||O[t] == 't')) {
	  alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log07;
	}

	/* adjustment based on probability distribution */
	start_freq=0;
	freq_id = 0;
        double sub_sum = 0;
        int sub_count = 0;

        if (t>=30){
                for(i=-30; i<=30; i++){
			if (t+i+2 < len_seq)
			{	
				start_freq -= hmm_ptr->tr_S[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
			}
                }
        }else{
                for(i=(-1*t); i<=30; i++){
			if (t+i+2 < len_seq)
			{	
				sub_sum += hmm_ptr->tr_S[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
			}
                }
                sub_sum = sub_sum * 61 / (30 + t + 1);
                start_freq -= sub_sum;
        }

	h_kd = hmm_ptr->S_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S_dist[1],2)/(2*pow(hmm_ptr->S_dist[0],2)));
	r_kd = hmm_ptr->S_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S_dist[4],2)/(2*pow(hmm_ptr->S_dist[3],2)));
	p_kd = h_kd / (h_kd + r_kd);
	if (p_kd<0.01){
	  p_kd=0.01;
	}else if (p_kd>0.99){
	  p_kd=0.99;
	}
	alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(p_kd);

      }
    }

    /**********************************************/
    /* END' state                                 */
    /* originally start codon of genes in - strand */
    /**********************************************/
    if (alpha[E_STATE_1][t] == 0){

      alpha[E_STATE_1][t] = max_dbl;
      path[E_STATE_1][t] = NOSTATE;

      if (t < len_seq - 2 && (O[t] == 'C'||O[t] =='c') && (O[t+1] == 'A'||O[t+1] == 'a') &&
	  ((O[t+2] == 'T'||O[t+2] =='t') || (O[t+2] == 'C'||O[t+2] =='c') || (O[t+2] == 'A'||O[t+2] =='a'))) {

	/* transition from frame6 */
	alpha[E_STATE_1][t+2] = alpha[M6_STATE_1][t-1] - hmm_ptr->tr[TR_GE];
	path[E_STATE_1][t] = M6_STATE_1;
	alpha[E_STATE_1][t] = max_dbl;
	alpha[E_STATE_1][t+1] = max_dbl;
	path[E_STATE_1][t+1] = E_STATE_1;
	path[E_STATE_1][t+2] = E_STATE_1;

	if ((O[t+2] == 'T'||O[t+2] == 't') ){
	  alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log83;
	}else if ((O[t+2] == 'C'||O[t+2] =='c') ){
	  alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(0.10);
	}else if((O[t+2] == 'A'||O[t+2] =='a') ){
	  alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log07;
	}

	/* adjustment based on probability distribution */
	start_freq=0;
	freq_id = 0;

        double sub_sum = 0;
        int sub_count = 0;

        if (t>=30){
                for(i=-30; i<=30; i++){
			if (t+i+2 < len_seq)
			{	
				start_freq -= hmm_ptr->tr_E_1[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
			}
                }
        }else{
                for(i=(-1*t); i<=30; i++){
			if (t+i+2 < len_seq)
			{	
				sub_sum += hmm_ptr->tr_E_1[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
			}
                }
                sub_sum = sub_sum * 61 / (30 + t + 1);
                start_freq -= sub_sum;
        }

	h_kd = hmm_ptr->E1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)));
	r_kd = hmm_ptr->E1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)));
	p_kd = h_kd / (h_kd + r_kd);

	if (p_kd<0.01){
	  p_kd=0.01;
	}else if (p_kd>0.99){
	  p_kd=0.99;
	}
	alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(p_kd);
      }
    }
    if (num_N>9){

      for (i=0; i<NUM_STATE; i++){
	if (i!=R_STATE){
	  alpha[i][t] = max_dbl;
	  path[i][t] = R_STATE;
	}
      }
    }
  }

  /***********************************************************/
  /* backtrack array to find the optimal path                */
  /***********************************************************/

  head_short = strtok(head, delimi);
  fprintf(fp_out, "%s\n", head_short); //use head_short, Ye, April 22, 2016

  /* find the state for O[N] with the highest probability */
  prob = max_dbl;
  for (i = 0; i < hmm_ptr->N; i++){

    if (alpha[i][len_seq-1] < prob){
      prob = alpha[i][len_seq-1];
      vpath[len_seq-1] = i;
    }
  }

  /* backtrack the optimal path */
  for(t=len_seq-2; t>=0; t--){
    vpath[t] = path[vpath[t+1]][t+1];
  }

  print_save = 0;
  codon_start=0;
  start_t=-1;

  int glen = strlen(O);
  char codon[4], utr[65];
  for (t=0; t<len_seq; t++){

    if (codon_start==0 && start_t < 0 &&
	((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
	 (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1) ||
	 vpath[t] == S_STATE || vpath[t] == S_STATE_1 )){
      start_t=t+1;
      //printf("Note assign start_t %d (t=%d vpath %d)\n", start_t, t, vpath[t]);
    }

    if (codon_start==0 &&
	(vpath[t]==M1_STATE || vpath[t]==M4_STATE ||
	 vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)){

      memset(dna,0,300000);
      memset(dna1,0,300000);
      memset(dna_f,0,300000);
      memset(dna_f1,0,300000);
      memset(protein,0, 100000);
      memset(insert,0,100);
      memset(delete,0,100);

      insert_id = 0;
      delete_id = 0;
      dna_id = 0;
      dna_f_id = 0;
      dna[dna_id]=O[t];
      dna_start_t = t + 1; //Ye April 21, 2016
      dna_f[dna_f_id]=O[t];
      //printf("Note start dna: t = %d, dna_id %d, dna_f_id %d, add %c\n", t, dna_id, dna_f_id, O[t]);
      start_orf=t+1;
      prev_match = vpath[t];

      if (vpath[t] < M6_STATE){
	codon_start=1;
      }else{
	codon_start=-1;
      }

    }else if (codon_start!=0 && (vpath[t]==E_STATE || vpath[t]==E_STATE_1 || t==len_seq-1)){

      if (vpath[t]==E_STATE || vpath[t]==E_STATE_1){
	end_t=t+3;
      }else{
	end_t=t+1;

	/* FGS1.12 start: remove incomplete codon */
	temp_t = t;
	while(vpath[temp_t] != M1_STATE && vpath[temp_t] != M4_STATE  &&
vpath[temp_t] != M1_STATE_1  && vpath[temp_t] != M4_STATE_1){
	  dna_f[dna_f_id] = '\0';
	  dna_f_id--;

	  dna[dna_id] = '\0';
	  dna_id--;

	  temp_t--;
	}
	/* FGS1.12 end: remove incomplete codon */
      }
      final_score = (alpha[vpath[end_t-4]][end_t-4]- alpha[vpath[start_t+2]][start_t+2] )/(end_t-start_t-5);
      frame = start_orf%3;
      if (frame==0){
	frame=3;
      }

      if (dna_id > gene_len  ){
	if (codon_start==1){
          if(start_t == dna_start_t - 3) { //add complete start codon to dna, Ye April 21, 2016
		strcpy(dna_tmp, dna);
		sprintf(dna, "%c%c%c%s", O[start_t-1], O[start_t], O[start_t+1], dna_tmp);
		//printf("add start codon to dna: %c%c%c\n", O[start_t-1], O[start_t], O[start_t+1]);
		//printf("old dna %d %s\n", strlen(dna_tmp), dna_tmp);
		//printf("new dna %d %s\n", strlen(dna), dna);
	  }
	  if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
  	    int start_old = start_t;
	    codon[0] = 0;
            strncpy(codon, O + start_old-1, 3);
	    codon[3] = 0;
	    int s = 0;
	    //find the optimal start codon within 30bp up- and downstream of start codon
	    double e_save;
            int s_save;
	    while((!(!strcmp(codon, "TAA") || !strcmp(codon, "TAG") || !strcmp(codon, "TGA"))) && (start_old-1-s-35>=0)) {
	      if(!strcmp(codon, "ATG") || !strcmp(codon, "GTG") || !strcmp(codon, "TTG")) {
		utr[0] = 0;
		strncpy(utr, O+start_old-1-s-30,63);
		utr[63] = 0;
		//printf("check s=%d, codon %s\n", s, codon);
		double freq_sum = 0;
		for(j = 0; j < strlen(utr) - 2; j ++) {
		   int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]); 
		   freq_sum -= train_ptr->start[cg][j][idx];
		   //printf("j=%d, key=%c%c%c %d, start %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->start[cg][j][idx]);
		}
		if(s == 0) { e_save = freq_sum; s_save = s; }
		else if(freq_sum < e_save) { e_save = freq_sum; s_save = -1 * s; }
		//printf("s=%d freq_sum %lf\n", s, freq_sum);
		//getchar();
	      }
	      s += 3;
	      codon[0] = 0;
	      strncpy(codon, O+start_old-1-s, 3);
	      codon[3] = 0;
	    }
	    start_t = start_old+s_save;
	    //update dna
	    if(s_save != 0) {
              //printf("start refined + %d -> %d\n", start_old, start_t);
  	      dna[0] = 0;
	      strncpy(dna, O + start_t - 1, end_t - start_t + 1);
	      dna[end_t - start_t + 1] = 0;
            }
 	  }
	  fprintf(fp_out, "%d\t%d\t+\t%d\t%lf\t", start_t, end_t, frame, final_score);

	  fprintf(fp_out, "I:");
	  for (i=0; i<insert_id; i++){
	    fprintf(fp_out, "%d,", insert[i]);
	  }
	  fprintf(fp_out, "\tD:");
	  for (i=0; i<delete_id; i++){
	    fprintf(fp_out, "%d,", delete[i]);
	  }
	  fprintf(fp_out, "\n");

	  fprintf(fp_aa, "%s_%d_%d_+\n", head_short, start_t, end_t);
	  fprintf(fp_dna, "%s_%d_%d_+\n", head_short, start_t, end_t);

	  //printf("dna-start %c%c%c exp %c%c%c (%c%c%c)\n", dna[0], dna[1], dna[2], O[start_t-1], O[start_t], O[start_t+1], O[start_t+2], O[start_t+3], O[start_t+4]);
	  //printf("dna-len %d start_t %d end_t %d exp-len %d diff %d\n", strlen(dna), start_t, end_t, end_t - start_t + 1, end_t - start_t + 1 - strlen(dna));
	  get_protein(dna,protein,1, whole_genome);
	  fprintf(fp_aa, "%s\n", protein);
	  if (format==0){
	    fprintf(fp_dna, "%s\n", dna);
	  }else if (format==1){
	    fprintf(fp_dna, "%s\n", dna_f);
	  }
	}else if (codon_start==-1){
          //printf("reverse strand dna-len %d, start_t %d, dna_start %d, add-up %d, end_t %d\n", strlen(dna), start_t, dna_start_t, dna_start_t + strlen(dna), end_t);
	  //getchar();
          if(dna_start_t + strlen(dna) == end_t - 2) { //add complete start codon (on reverse strand) to dna, Ye April 21, 2016
		strcpy(dna_tmp, dna);
		sprintf(dna, "%s%c%c%c", dna_tmp, O[end_t-3], O[end_t-2], O[end_t-1]);
		//printf("add start codon on the reverse strand to dna: %c%c%c\n", O[end_t-3], O[end_t-2], O[end_t-1]);
	  }
	  if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
  	    int end_old = end_t; //reverse
	    codon[0] = 0;
            strncpy(codon, O + end_t-1-2, 3);
	    codon[3] = 0;
	    int s = 0;
	    //find the optimal start codon within 30bp up- and downstream of start codon
	    double e_save;
            int s_save;
	    while((!(!strcmp(codon, "TTA") || !strcmp(codon, "CTA") || !strcmp(codon, "TCA"))) && (end_old-2+s+35 < glen)) {
	      if(!strcmp(codon, "CAT") || !strcmp(codon, "CAC") || !strcmp(codon, "CAA")) {
		utr[0] = 0;
		strncpy(utr, O+end_old-1-2+s-30,63);
		utr[63] = 0;
		//printf("check s=%d, codon %s\n", s, codon);
		double freq_sum = 0;
		for(j = 0; j < strlen(utr) - 2; j ++) {
		   int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]); 
		   freq_sum -= train_ptr->stop1[cg][j][idx]; //stop1?? Ye, April 18, 2016
		   //printf("j=%d, key=%c%c%c %d, stop1 %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->stop1[cg][j][idx]);
		}
		if(s == 0) { e_save = freq_sum; s_save = s; }
		else if(freq_sum < e_save) { e_save = freq_sum; s_save = s; }
		//printf("s=%d freq_sum %lf\n", s, freq_sum);
		//getchar();
	      }
	      s += 3;
	      codon[0] = 0;
	      strncpy(codon, O+end_old-1-2+s, 3);
	      codon[3] = 0;
	    }
	    end_t = end_old+s_save;
	    //update dna
	    if(s_save != 0) {
              //printf("start refined - end %d -> %d\n", end_old, end_t);
   	      dna[0] = 0;
	      strncpy(dna, O + start_t - 1, end_t - start_t + 1);
	      dna[end_t - start_t + 1] = 0;
            }
 	  }

	  fprintf(fp_out, "%d\t%d\t-\t%d\t%lf\t", start_t, end_t, frame, final_score);
	  fprintf(fp_out, "I:");
	  for (i=0; i<insert_id; i++){
	    fprintf(fp_out, "%d,", insert[i]);
	  }
	  fprintf(fp_out, "\tD:");
	  for (i=0; i<delete_id; i++){
	    fprintf(fp_out, "%d,", delete[i]);
	  }
	  fprintf(fp_out, "\n");

	  fprintf(fp_aa, "%s_%d_%d_-\n", head_short, start_t, end_t);
	  fprintf(fp_dna, "%s_%d_%d_-\n", head_short, start_t, end_t);

	  get_protein(dna,protein,-1, whole_genome);
	  get_rc_dna(dna, dna1);
	  get_rc_dna_indel(dna_f, dna_f1);
	  fprintf(fp_aa, "%s\n", protein);
	  if (format==0){
	    fprintf(fp_dna, "%s\n", dna1);
	  }else if (format==1){
	    fprintf(fp_dna, "%s\n", dna_f1);
	  }
	}
      }
      codon_start=0;
      start_t = -1;
      end_t = -1;
      dna_id = 0;
      dna_f_id = 0;

    }else if (codon_start!=0 &&
	      ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
	       (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1)) &&
	      vpath[t]-prev_match<6){

      if (vpath[t] < prev_match){
	out_nt = vpath[t]+6-prev_match;
      }else{
	out_nt = vpath[t]-prev_match;
      }
      for (kk=0; kk<out_nt; kk++){   /* for deleted nt in reads */
	dna_id ++;
	dna[dna_id] = 'N';
        //printf("dna_id %d, dna-len %d\n", dna_id, strlen(dna)); 
	dna_f_id ++;
	dna_f[dna_f_id] = 'x';
	if (kk>0){
	  delete[delete_id]=t+1;
	  delete_id++;
	}
      }
      dna[dna_id]=O[t];
      //printf("dna_id %d, add %d %c dna-len %d\n", dna_id, t, O[t], strlen(dna)); 
      dna_f[dna_f_id]=O[t];
      prev_match = vpath[t];

    }else if (codon_start!=0 &&
	      ((vpath[t]>=I1_STATE && vpath[t]<=I6_STATE) ||
	       (vpath[t]>=I1_STATE_1 && vpath[t]<=I6_STATE_1))){
      dna_f_id ++;
      dna_f[dna_f_id] = tolower(O[t]);
      insert[insert_id]=t+1;
      insert_id++;

    }
    else if (codon_start!=0 && vpath[t]==R_STATE){
      /* for long NNNNNNNNN, pretend R state */
      codon_start=0;
      start_t=-1;
      end_t = -1;
      dna_id=0;
      dna_f_id=0;

    }
  }

  free_dmatrix(alpha, hmm_ptr->N);
  free_imatrix(path, hmm_ptr->N);
  free_ivector(vpath);
}

int get_prob_from_cg(HMM *hmm_ptr, TRAIN *train_ptr, char *O){ //change from void to int, Ye, April 18, 2016
  int cg_id = -1;
  int cg_count=0;
  int len_seq;
  int i,j,k;

  len_seq = strlen(O);
  for (i=0; i<len_seq; i++){
    if ((O[i] == 'C'||O[i] =='c') || (O[i] == 'G'||O[i] == 'g') ){
      cg_count++;
    }
  }
   cg_count = floor((cg_count*1.0/len_seq)*100)-26;
   if (cg_count < 0){
     cg_count = 0;
   }else if (cg_count > 43){
     cg_count = 43;
   }

  memcpy(hmm_ptr->e_M, train_ptr->trans[cg_count], sizeof(hmm_ptr->e_M));
  memcpy(hmm_ptr->e_M_1, train_ptr->rtrans[cg_count],
sizeof(hmm_ptr->e_M_1));
  memcpy(hmm_ptr->tr_R_R, train_ptr->noncoding[cg_count],
sizeof(hmm_ptr->tr_R_R));
  memcpy(hmm_ptr->tr_S, train_ptr->start[cg_count], sizeof(hmm_ptr->tr_S));
  memcpy(hmm_ptr->tr_E, train_ptr->stop[cg_count], sizeof(hmm_ptr->tr_E));
  memcpy(hmm_ptr->tr_S_1, train_ptr->start1[cg_count],
sizeof(hmm_ptr->tr_S_1));
  memcpy(hmm_ptr->tr_E_1, train_ptr->stop1[cg_count],
sizeof(hmm_ptr->tr_E_1));
  memcpy(hmm_ptr->S_dist, train_ptr->S_dist[cg_count],
sizeof(hmm_ptr->S_dist));
  memcpy(hmm_ptr->E_dist, train_ptr->E_dist[cg_count],
sizeof(hmm_ptr->E_dist));
  memcpy(hmm_ptr->S1_dist, train_ptr->S1_dist[cg_count],
sizeof(hmm_ptr->S1_dist));
  memcpy(hmm_ptr->E1_dist, train_ptr->E1_dist[cg_count],
sizeof(hmm_ptr->E1_dist));
 
  return cg_count;
}



void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename,
			 char *sfilename,char *pfilename,char *s1filename,char *p1filename,char *dfilename, TRAIN *train_ptr){

  int i, j, k, p;
  double prob;
  FILE *fp, *fpm, *fpm1, *fpn, *fps, *fpp, *fps1, *fpp1, *fpd;

  char name[10];
  char head[20];
  char start[10];
  char end[10];

  /* probabilities saved in log Ye April 18, 2016 */
  /* start <- ./train/start (start in forward)
     stop <- ./train/stop (stop in forward)
     start1 <- ./train/stop1 (start in reverse)
     stop1 <- ./train/start1 (stop in reverse)
  */

  /****************************************************/
  /* transition                                       */
  /****************************************************/
  fp = fopen (filename , "r");

  /* Transition */
  fscanf(fp, "%s", head);
  for (i=0; i<14; i++){
    fscanf(fp, "%s %lf", name, &prob);
    hmm_ptr->tr[tr2int(name)] = log(prob);
  }

  /* TransitionMI */
  fscanf(fp, "%s", head);
  for (i=0; i<16; i++){
    fscanf(fp, "%s %s %lf\n", start, end, &prob);
    hmm_ptr->tr_M_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
  }

  /* TransitionII */
  fscanf(fp, "%s", head);
  for (i=0; i<16; i++){
    fscanf(fp, "%s %s %lf", start, end, &prob);
    hmm_ptr->tr_I_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
 }

  /* PI */
  fscanf(fp, "%s", head);
  for (i=0; i<NUM_STATE; i++){
    fscanf(fp, "%s %lf", name, &prob);
    hmm_ptr->pi[i] = log(prob);
  }
  fclose(fp);


  /****************************************************/
  /* M state transition                               */
  /****************************************************/
  fpm = fopen (mfilename , "r");
  for (p=0; p<44; p++){                        /* cg */
    fscanf(fpm, "%s", head);
    for (i=0; i<6; i++){                       /* period */
      for (j=0; j<16; j++){                    /* condition */
	for (k=0; k<4; k++){                   /* emission */
	  fscanf(fpm, "%lf", &prob);
	  train_ptr->trans[p][i][j][k] = log(prob);
	}
      }
    }
  }
  fclose(fpm);


  /****************************************************/
  /* M state_1 transition                             */
  /****************************************************/
  fpm1 = fopen (mfilename1 , "r");
  for (p=0; p<44; p++){
    fscanf(fpm1, "%s", head);
    for (i=0; i<6; i++){
      for (j=0; j<16; j++){
	for (k=0; k<4; k++){
	  fscanf(fpm1, "%lf", &prob);
	  train_ptr->rtrans[p][i][j][k] = log(prob);
	}
      }
    }
  }
  fclose(fpm1);


  /****************************************************/
  /* noncoding state  transition                      */
  /****************************************************/
  fpn = fopen (nfilename, "r");
  for (p=0; p<44; p++){
    fscanf(fpn, "%s", head);
    for (j=0; j<4; j++){
      for (k=0; k<4; k++){
	fscanf(fpn, "%lf", &prob);
	train_ptr->noncoding[p][j][k] = log(prob);
      }
    }
  }
  fclose(fpn);


  /****************************************************/
  /* start                                            */
  /****************************************************/
  fps = fopen (sfilename, "r");
  for (p=0; p<44; p++){
    fscanf(fps, "%s", head);
    for (j=0; j<61; j++){
      for (k=0; k<64; k++){
	fscanf(fps, "%lf", &prob);
	train_ptr->start[p][j][k] = log(prob);
      }
    }
  }
  fclose(fps);


  /****************************************************/
  /* stop                                             */
  /****************************************************/
  fpp = fopen (pfilename, "r"); //sfilename->pfilename, Ye, April 18, 2016
  for (p=0; p<44; p++){
    fscanf(fpp, "%s", head);
    //for (j=0; j<58; j++){ //58->61, Ye, April 18, 2016
    for (j=0; j<61; j++){
      for (k=0; k<64; k++){
	fscanf(fpp, "%lf", &prob);
	train_ptr->stop[p][j][k] = log(prob);
      }
    }
  }
  fclose(fpp);


  /****************************************************/
  /* start1                                           */
  /****************************************************/
  fps1 = fopen (s1filename, "r");
  for (p=0; p<44; p++){
    fscanf(fps1, "%s", head);
    for (j=0; j<61; j++){ //58->61 Ye, April 18, 2016
      for (k=0; k<64; k++){
	fscanf(fps1, "%lf", &prob);
	train_ptr->start1[p][j][k] = log(prob);
      }
    }
  }
  fclose(fps1);


  /****************************************************/
  /* stop1                                            */
  /****************************************************/
  fpp1 = fopen (p1filename, "r");
  for (p=0; p<44; p++){
    fscanf(fpp1, "%s", head);
    for (j=0; j<61; j++){
      for (k=0; k<64; k++){
	fscanf(fpp1, "%lf", &prob);
	train_ptr->stop1[p][j][k] = log(prob);
      }
    }
  }
  fclose(fpp1);


  /****************************************************/
  /* pwm distribution                                 */
  /* S_dist, E_dist, S1_dist, E1_dist NOT in log      */
  /****************************************************/
  fpd = fopen (dfilename, "r");
  for (p=0; p<44; p++){
    fscanf(fpd, "%s", head);
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->S_dist[p][k] = prob;
    }
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->E_dist[p][k] = prob;
    }
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->S1_dist[p][k] = prob;
    }
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->E1_dist[p][k] = prob;
    }
  }
  fclose(fpd);

}

void free_hmm(HMM *hmm_ptr){

  free_dvector(hmm_ptr->pi);
}

void dump_memory(void *p, int size)
{
	int i, s;
	double *c;
	c = (double*)p;
	s = size / sizeof(double);

	printf("Dump size %d\n", size);
	for (i = 0; i < s; i++)
	{
		if (i > 0 && i % 10 == 0)
		{
			printf("\n");
		}
		printf("%f ", *c);
		c++;
	}
	printf("\n");
}

