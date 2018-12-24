#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "hmm.h"

#include <pthread.h>

#define ADD_LEN 1024
#define STRINGLEN 4096

typedef struct thread_data
{
	FILE *out;
	FILE *aa;
	FILE *dna;
	char *obs_head;
	char *obs_seq;
	int wholegenome;
	int cg;
	int format;
	HMM *hmm;
	TRAIN *train;
} thread_data;

void* thread_func(void *threadarr);

int main (int argc, char **argv)
{
  clock_t start = clock();
  int i, j, c, max;
  HMM hmm;
  char *obs_seq, *obs_head;
  TRAIN train;
  int wholegenome;
  int format=0;
  FILE *fp_out, *fp_aa, *fp_dna, *fp;
  char hmm_file[STRINGLEN] = "";
  char out_header[STRINGLEN] = "";
  char aa_file[STRINGLEN] = "";
  char seq_file[STRINGLEN] = "";
  char out_file[STRINGLEN] = "";
  char dna_file[STRINGLEN] = ""; 
  char train_file[STRINGLEN] = "";
  char mstate_file[STRINGLEN] = "";
  char rstate_file[STRINGLEN] = "";
  char nstate_file[STRINGLEN] = "";
  char sstate_file[STRINGLEN] = "";
  char pstate_file[STRINGLEN] = "";
  char s1state_file[STRINGLEN] = "";     /* stop codon of gene in - stand */
  char p1state_file[STRINGLEN] = "";
  char dstate_file[STRINGLEN] = "";
  char train_dir[STRINGLEN] = "";
  int count=0;
  int currcount = 0;
  int total = 0;
  char mystring[STRINGLEN] = "";
  int *obs_seq_len;
  int bp_count;  /* count the length of each line in input file */

  int threadnum = 1;
  int rc;

  thread_data *threadarr;
  char **lastline, **currline;

  strncpy(train_dir, argv[0], strlen(argv[0])-12);
  strcat(train_dir, "train/");
  strcpy(mstate_file, train_dir);
  strcat(mstate_file, "gene");
  strcpy(rstate_file, train_dir);
  strcat(rstate_file, "rgene");
  strcpy(nstate_file, train_dir);
  strcat(nstate_file, "noncoding");
  strcpy(sstate_file, train_dir);
  strcat(sstate_file, "start");
  strcpy(pstate_file, train_dir);
  strcat(pstate_file, "stop");
  strcpy(s1state_file, train_dir);
  strcat(s1state_file, "stop1");
  strcpy(p1state_file, train_dir);
  strcat(p1state_file, "start1");
  strcpy(dstate_file, train_dir);
  strcat(dstate_file, "pwm");


  /* read command line argument */
  if (argc <= 8){    
    fprintf(stderr, "ERROR: You missed some parameters for input\n");
    print_usage();
    exit(EXIT_FAILURE);
  }

  while ((c=getopt(argc, argv, "fs:o:w:t:p:")) != -1){
    switch (c){
    case 's':
      strcpy(seq_file, optarg);
      if (access(seq_file, F_OK)==-1){
	fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
	print_usage();
	exit(EXIT_FAILURE);
      }
      break;  
    case 'w':
      wholegenome = atoi(optarg);
      if (wholegenome != 0 && wholegenome != 1){
	fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
	print_usage();
	exit(EXIT_FAILURE);
      }
      break;
    case 'p':
      threadnum = atoi(optarg);
      if (threadnum < 1){
	fprintf(stderr, "ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
	print_usage();
	exit(EXIT_FAILURE);
      }
      printf("Using %d threads.\n", threadnum);
      break;
    case 'o':
      strcpy(out_header, optarg);
      break;
    case 't':
      strcpy(train_file, optarg);
      strcpy(hmm_file, train_dir);
      strcat(hmm_file, train_file);

      if (access(hmm_file, F_OK)==-1){
	fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
	print_usage();
	exit(EXIT_FAILURE);
      }
      break;
    case 'f':
      format = 1;
      break;
    }
  }

  
  /* check whether the specified files exist */
  if (access(mstate_file, F_OK)==-1){
    fprintf(stderr, "Forward prob. file [%s] does not exist\n", mstate_file);
    exit(1);
  }
  if (access(rstate_file, F_OK)==-1){
    fprintf(stderr, "Backward prob. file [%s] does not exist\n", rstate_file);
    exit(1);
  }
  if (access(nstate_file, F_OK)==-1){
    fprintf(stderr, "noncoding prob. file [%s] does not exist\n", nstate_file);
    exit(1);
  }
  if (access(sstate_file, F_OK)==-1){
    fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
    exit(1);
  }
  if (access(pstate_file, F_OK)==-1){
    fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
    exit(1);
  }
  if (access(s1state_file, F_OK)==-1){
    fprintf(stderr, "start1 prob. file [%s] does not exist\n", s1state_file);
    exit(1);
  }
  if (access(p1state_file, F_OK)==-1){
    fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
    exit(1);
  }
  if (access(dstate_file, F_OK)==-1){
    fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
    exit(1);
  }
  if (access(hmm_file, F_OK)==-1){
    fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
    exit(1);
  }
  
  /* read all initial model */
  hmm.N=NUM_STATE;
  get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

  // Initialize thread data structure
  threadarr = (thread_data*)malloc(sizeof(thread_data) * threadnum);
  memset(threadarr, '\0', sizeof(thread_data) * threadnum);
  for (i = 0; i < threadnum; i++)
  {
    if(threadnum > 1) sprintf(mystring, "%s.out.tmp.%d", out_header, i);
    else sprintf(mystring, "%s.out", out_header); 
    threadarr[i].out = fopen(mystring, "w");
    if(threadnum > 1) sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
    else sprintf(mystring, "%s.faa", out_header);
    threadarr[i].aa = fopen(mystring, "w");
    if(threadnum > 1) sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
    else sprintf(mystring, "%s.ffn", out_header);
    threadarr[i].dna = fopen(mystring, "w");

    threadarr[i].hmm = (HMM*)malloc(sizeof(HMM));
    memcpy(threadarr[i].hmm, &hmm, sizeof(HMM));
    threadarr[i].train = (TRAIN*)malloc(sizeof(TRAIN));
    memcpy(threadarr[i].train, &train, sizeof(TRAIN));

    //threadarr[i].hmm->N=NUM_STATE;
    //get_train_from_file(hmm_file, threadarr[i].hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, threadarr[i].train);

    threadarr[i].wholegenome = wholegenome;
    threadarr[i].format = format;
  }

  pthread_t *thread;
  thread = (pthread_t*)malloc(sizeof(thread) * threadnum);
  memset(thread, '\0', sizeof(thread) * threadnum);
  void *status;
  fp = fopen (seq_file, "r");
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
      count++;
    }
  }
  obs_seq_len = (int *)malloc(count * sizeof(int));
  printf("no. of seqs: %d\n", count);  

  i = 0;
  count = 0;
  rewind(fp);
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
      if (i>0){
        obs_seq_len[count] = i;
        count++;
      }
      i = 0;
    }else{
      bp_count = strlen(mystring);
      while(mystring[bp_count-1] == 10 || mystring[bp_count-1]==13){
	bp_count --;
      }

      i += bp_count;
    }
  }
  obs_seq_len[count] = i;

  rewind(fp);
  total = 0;
  count = 0;
  j = 0;

  while (!(feof(fp)))
  {
    memset(mystring, '\0', sizeof mystring);
    fgets (mystring , sizeof mystring  , fp);
    bp_count = strlen(mystring);
    while(mystring[bp_count - 1] == 10 || mystring[bp_count - 1]==13){
      //mystring[bp_count - 1] = 0;
      bp_count --;
    }

    if (mystring[0] == '>' || feof(fp)){
      if (feof(fp))
      {
        memcpy(threadarr[currcount].obs_seq + j, mystring, bp_count);
        j += bp_count;
        //max = appendSeq(mystring, &(threadarr[currcount].obs_seq), max);
      }
      if ((count > 0 && count % threadnum == 0) || feof(fp))
      {
        // Deal with the thread
	for (i = 0; i < count; i++)
	{
	  rc = pthread_create(&thread[i], NULL, thread_func, (void*)&threadarr[i]);
	  if (rc)
	  {
	    printf("Error: Unable to create thread, %d\n", rc);
	    exit(-1);
	  }
        }
	for (i = 0; i < count; i++)
	{
	  rc = pthread_join(thread[i], &status);
	  if (rc)
	  {
	    printf("Error: Unable to join threads, %d\n", rc);
	    exit(-1);
	  }
	}
	for (i = 0; i < count; i++)
	{
	  free(threadarr[i].obs_head);
	  free(threadarr[i].obs_seq);
          threadarr[i].obs_head = NULL;
          threadarr[i].obs_seq = NULL;
	}

	count = 0;
      }

      if (!(feof(fp)))
      {
        threadarr[count].obs_head = (char *)malloc((bp_count+1) * sizeof(char));
        memset(threadarr[count].obs_head, 0, (bp_count+1) * sizeof(char));
        memcpy(threadarr[count].obs_head, mystring, bp_count);
        //threadarr[count].obs_seq = NULL;
        threadarr[count].obs_seq = (char*)malloc((obs_seq_len[total] + 1) * sizeof(char));
        memset(threadarr[count].obs_seq, '\0', (obs_seq_len[total] + 1) * sizeof(char));
        total++;
        currcount = count;
        count++;
        j = 0;
        max = 0;
      }

    }else{
      memcpy(threadarr[currcount].obs_seq + j, mystring, bp_count);
      j += bp_count;
      //max = appendSeq(mystring, &(threadarr[currcount].obs_seq), max);
    }
    if (feof(fp))
    {
      break;
    }
  }
  for (i = 0; i < threadnum; i++)
  {
    fclose(threadarr[i].out);
    fclose(threadarr[i].aa);
    fclose(threadarr[i].dna);
  }

  if(threadnum > 1) {
    /* create output file name */
    strcpy(aa_file, out_header);
    strcat(aa_file, ".faa");
    strcpy(dna_file, out_header);
    strcat(dna_file, ".ffn");
    strcpy(out_file, out_header);
    strcat(out_file, ".out");

    remove (out_file);
    remove (aa_file);
    remove (dna_file);

    fp_aa = fopen (aa_file , "w");
    fp_out = fopen (out_file , "w");
    fp_dna = fopen (dna_file , "w");

    lastline = (char**)malloc(sizeof(char*) * threadnum);
    memset(lastline, '\0', sizeof(char*) * threadnum);
    currline = (char**)malloc(sizeof(char*) * threadnum);
    memset(currline, '\0', sizeof(char*) * threadnum);
    for (i = 0; i < threadnum; i++)
    {
      sprintf(mystring, "%s.out.tmp.%d", out_header, i);
      threadarr[i].out = fopen(mystring, "r");
      sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
      threadarr[i].aa = fopen(mystring, "r");
      sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
      threadarr[i].dna = fopen(mystring, "r");

      lastline[i] = (char*)malloc(sizeof(char) * (STRINGLEN + 1));
      memset(lastline[i], '\0', sizeof(char) * (STRINGLEN + 1));
      currline[i] = (char*)malloc(sizeof(char) * (STRINGLEN + 1));
      memset(currline[i], '\0', sizeof(char) * (STRINGLEN + 1));
    }

    // Organize out file
    while (1)
    {
      j = 0;
      for (i = 0; i < threadnum; i++)
      {
        if (lastline[i][0] != '\0')
        {
          fputs(lastline[i], fp_out);
          lastline[i][0] = '\0';
        }
        while(fgets(currline[i], STRINGLEN, threadarr[i].out))
        {
          if (currline[i][0] == '>')
          {
            memcpy(lastline[i], currline[i], strlen(currline[i]) + 1);
            break;
          }
          else
          {
            fputs(currline[i], fp_out);
          }
        }
        if (feof(threadarr[i].out))
        {
          j++;
        }
      }
      if (j == threadnum)
      {
        break;
      }
    }
    // Organize faa file
    for (i = 0; i < threadnum; i++)
    {
      lastline[i][0] = '\0';
    }
    while (1)
    {
      j = 0;
      for (i = 0; i < threadnum; i++)
      {
        if (lastline[i][0] != '\0')
        {
          fputs(lastline[i], fp_aa);
          lastline[i][0] = '\0';
        }
        while(fgets(currline[i], STRINGLEN, threadarr[i].aa))
        {
          if (currline[i][0] == '>')
          {
            memcpy(lastline[i], currline[i], strlen(currline[i]) + 1);
            break;
          }
          else
          {
            fputs(currline[i], fp_aa);
          }
        }
        if (feof(threadarr[i].aa))
        {
          j++;
        }
      }
      if (j == threadnum)
      {
        break;
      }
    }

    // Organize dna file
    for (i = 0; i < threadnum; i++)
    {
      lastline[i][0] = '\0';
    }
    while (1)
    {
      j = 0;
      for (i = 0; i < threadnum; i++)
      {
        if (lastline[i][0] != '\0')
        {
          fputs(lastline[i], fp_dna);
          lastline[i][0] = '\0';
        }
        while(fgets(currline[i], STRINGLEN, threadarr[i].dna))
        {
          if (currline[i][0] == '>')
          {
            memcpy(lastline[i], currline[i], strlen(currline[i]) + 1);
            break;
          }
          else
          {
            fputs(currline[i], fp_dna);
          }
        }
        if (feof(threadarr[i].dna))
        {
          j++;
        }
      }
      if (j == threadnum)
      {
        break;
      }
    }

    for (i = 0; i < threadnum; i++)
    {
      fclose(threadarr[i].out);
      fclose(threadarr[i].aa);
      fclose(threadarr[i].dna);
      sprintf(mystring, "%s.out.tmp.%d", out_header, i);
      remove(mystring);
      sprintf(mystring, "%s.faa.tmp.%d", out_header, i);
      remove(mystring);
      sprintf(mystring, "%s.ffn.tmp.%d", out_header, i);
      remove(mystring);
      free(threadarr[i].hmm);
      free(threadarr[i].train);
      free(lastline[i]);
      free(currline[i]);
    }
    free(threadarr);
    free(lastline);
    free(currline);
  
    free(obs_seq_len);
    //free(obs_head);
    fclose(fp_out);
    fclose(fp_aa);
    fclose(fp_dna);
    fclose(fp);
  }
  clock_t end = clock();
  printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (end - start) / (60.0 * CLOCKS_PER_SEC));
}


void* thread_func(void *threadarr)
{
  thread_data *d;
  d = (thread_data*)threadarr;
  d->cg = get_prob_from_cg(d->hmm, d->train, d->obs_seq); //cg - 26 Ye April 16, 2016
  if (strlen(d->obs_seq)>70){
    viterbi(d->hmm, d->train, d->obs_seq, d->out, d->aa, d->dna, d->obs_head, d->wholegenome, d->cg, d->format);
  }
}

int appendSeq(char *input, char **seq, int input_max)
{
	int len, inputlen, max;
	char *tmp;

	max = input_max;
	if (*seq != NULL)
	{
		len = strlen(*seq);
	}
	else
	{
		len = 0;
	}
	inputlen = strlen(input);
	if ((len + inputlen) >= max)
	{
		while ((len + inputlen) >= max)
		{
			max += ADD_LEN;
		}
		tmp = (char*)malloc(sizeof(char) * max);
		memset(tmp, '\0', sizeof(char) * max);
		if (*seq != NULL)
		{
			memcpy(tmp, *seq, len);
		}
		free(*seq);
		*seq = tmp;
	}
	strcat(*seq, input);
	return max;
}

