# From https://github.com/kaust-vislab/python-data-science-project
FROM ubuntu:20.04

LABEL maintainer="johnlees <john@johnlees.me>"

SHELL [ "/bin/bash", "--login", "-c" ]

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create a non-root user
ARG username=poppipe-usr
ARG uid=1000
ARG gid=100
ENV USER $username
ENV UID $uid
ENV GID $gid
ENV HOME /home/$USER

RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --gid $GID \
    --home $HOME \
    $USER

COPY environment.yml /tmp/
RUN chown $UID:$GID /tmp/environment.yml

COPY docker/entrypoint.sh /usr/local/bin/
RUN chown $UID:$GID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

USER $USER

# install miniforge
# NB need to set homedir explicitly, or use `WORKDIR $HOME`
ENV CONDA_DIR $HOME/miniforge3
RUN wget -O ~/Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    chmod +x ~/Miniforge3.sh && \
    bash ~/Miniforge3.sh -b -p $CONDA_DIR && \
    rm ~/Miniforge3.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
RUN echo ". $CONDA_DIR/etc/profile.d/mamba.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interactive shells
RUN mamba shell init
RUN conda shell init

# create a project directory inside user home
# (this isn't used with a clone running snakemake)
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
# copy the code in
COPY . $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build the conda environment
ENV ENV_PREFIX $PROJECT_DIR/env
RUN mamba env create -r $ENV_PREFIX --prefix $ENV_PREFIX --file /tmp/environment.yml && \
    mamba clean --all --yes

# use an entrypoint script to insure conda environment is properly activated at runtime
USER root
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# default of running shell is fine
#CMD [ "bash" ]