CC=	gcc
CFLAG= -O3
SRCS=	util_lib.c hmm_lib.c run_hmm.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o 

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC)  $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o  -lm -lpthread

clean:
	rm -rf *.o FragGeneScan* *~
