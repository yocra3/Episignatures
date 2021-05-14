FROM nfcore/base:1.13.3
LABEL authors="Carlos Ruiz-Arenas" \
      description="Docker image containing the python packages for Episignatures project"

# Install mamba to manage conda
RUN conda install mamba -c conda-forge

# Install the conda environment
COPY environment.yml /
RUN mamba env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/episignatures/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name episignatures > episignatures.yml
