FROM continuumio/miniconda3

# set work dir
WORKDIR /app

# install conda packages
RUN conda create -y -n run_dbcan python=3.8 dbcan -c conda-forge -c bioconda && conda clean -ya

# add run_dbcan environment to path
RUN echo "source activate run_dbcan" > ~/.bashrc
ENV PATH /opt/conda/envs/run_dbcan/bin:$PATH

# download and make the databases
RUN test -d db || mkdir db
RUN cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa && diamond makedb --in CAZyDB.08062022.fa -d CAZy \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt && mv dbCAN-HMMdb-V11.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm && hmmpress tf-2.hmm \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff


ENTRYPOINT ["run_dbcan"]

CMD [ "run_dbcan -h" ]
