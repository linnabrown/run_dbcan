FROM conda/miniconda3

# install conda packages
RUN conda install -c anaconda conda-build \
  && conda install wget \
  && conda install -c bioconda -y diamond hmmer=3.1b2 prodigal fraggenescan \
  && conda clean -a

# install dbcan via pip
RUN pip install dbcan

# download and make the databases
RUN mkdir -p /app/db && cd /app/db \
  && wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.09242021.fa && diamond makedb --in CAZyDB.09242021.fa -d CAZy \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt && mv dbCAN-HMMdb-V10.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm

WORKDIR /app

CMD [ "run_dbcan -h" ]
