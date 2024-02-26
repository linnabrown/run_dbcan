FROM continuumio/miniconda3

# set work dir
WORKDIR /app
COPY dbcan.yml .

# install conda packages
RUN conda env create -f dbcan.yml && conda clean -ya

# add run_dbcan environment to path
RUN echo "source activate CAZyme_annotation" > ~/.bashrc
ENV PATH /opt/conda/envs/CAZyme_annotation/bin:$PATH

# download and make the databases
RUN dbcan_build --cpus 1 --db-dir db --clean

# define the entrypoint
ENTRYPOINT ["run_dbcan"]

# define the default command
CMD [ "run_dbcan -h" ]
