#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define A 0
#define C 1
#define G 2
#define T 3

#define NUM_STATE 29

#define NOSTATE -1
#define S_STATE 0
#define E_STATE 1
#define R_STATE 2
#define S_STATE_1 3
#define E_STATE_1 4
#define M1_STATE 5
#define M2_STATE 6
#define M3_STATE 7
#define M4_STATE 8
#define M5_STATE 9
#define M6_STATE 10
#define M1_STATE_1 11
#define M2_STATE_1 12
#define M3_STATE_1 13
#define M4_STATE_1 14
#define M5_STATE_1 15
#define M6_STATE_1 16
#define I1_STATE 17
#define I2_STATE 18
#define I3_STATE 19
#define I4_STATE 20
#define I5_STATE 21
#define I6_STATE 22
#define I1_STATE_1 23
#define I2_STATE_1 24
#define I3_STATE_1 25
#define I4_STATE_1 26
#define I5_STATE_1 27
#define I6_STATE_1 28


#define TR_MM 0
#define TR_MI 1
#define TR_MD 2
#define TR_II 3
#define TR_IM 4
#define TR_DD 5
#define TR_DM 6
#define TR_GE 7
#define TR_GG 8
#define TR_ER 9
#define TR_RS 10
#define TR_RR 11
#define TR_ES 12
#define TR_ES1 13


typedef struct {

  double  pi[29];    /* pi[1..N] pi[i] is the initial state distribution. */
  int N;               /* number of state */

  double tr[14];                 /* transition probability from a (delete/insert/match) state to a state */

  double e_M_1[6][16][4];      /* transition probability from a lowest-level state  to a  lowest-level state*/
  double e_M[6][16][4];

  double tr_R_R[4][4];
  double tr_I_I[4][4];
  double tr_M_I[4][4];

  double tr_S[61][64];
  double tr_E[61][64];
  double tr_S_1[61][64];
  double tr_E_1[61][64];

  double S_dist[6];  /*sigma, mu,alpha, sigma_r, mu_r, alpha_r */
  double E_dist[6];
  double S1_dist[6];
  double E1_dist[6];
} HMM;



typedef struct {

  double trans[44][6][16][4];
  double rtrans[44][6][16][4];
  double noncoding[44][4][4];
  double start[44][61][64];
  double stop[44][61][64];
  double start1[44][61][64];
  double stop1[44][61][64];

  double S_dist[44][6];
  double E_dist[44][6];
  double S1_dist[44][6];
  double E1_dist[44][6];

} TRAIN;



int get_prob_from_cg(HMM *hmm, TRAIN *train, char *O); //return cg - 26 Ye April 16, 2016
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename, char *sfilename,char *pfilename,char *s1filename,char *p1filename, char *dfilename, TRAIN *train_ptr);
void viterbi(HMM *hmm_ptr, TRAIN *train_ptr, char *O, FILE *out_filename, FILE *log_filename,FILE *dna_filename, char *head, int metagene, int cg, int format);
void free_hmm(HMM *hmm);
void get_protein(char *dna, char *protein, int strand, int whole_genome);
void get_rc_dna(char *dna, char *dna1);
void get_corrected_dna(char *dna, char *dna_f);
