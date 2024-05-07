FROM continuumio/miniconda3

# set work dir
WORKDIR /app
COPY . /app/

# install conda packages
RUN conda env create -f envs/docker_build.yml && conda clean -ya

# add run_dbcan environment to path
RUN echo "source activate run_dbcan" > ~/.bashrc
ENV PATH /opt/conda/envs/run_dbcan/bin:$PATH

# install the dbcan package
RUN pip install . && rm -rf .git

# download and make the databases
RUN dbcan_build --cpus 1 --db-dir db --clean

# define the entrypoint
ENTRYPOINT ["run_dbcan"]

# define the default command
CMD [ "run_dbcan -h" ]
