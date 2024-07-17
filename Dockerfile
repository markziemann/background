FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Update apt-get
RUN apt-get update \
        && apt-get upgrade -y \
        && apt-get install -y nano git  libncurses-dev \
        ## Install the python package magic wormhole to send files
        && pip install magic-wormhole           \
        ## Remove packages in '/var/cache/' and 'var/lib'
        ## to remove side-effects of apt-get update
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e 'install.packages(c("eulerr","gplots","kableExtra","parallel","RColorBrewer","stringi","HGNChelper","beeswarm","gridExtra","png","RhpcBLASctl","tictoc"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("biomaRt","clusterProfiler","DESeq2","edgeR","fgsea","getDEE2","limma","mitch"))'

# get a clone of the codes
RUN git clone https://github.com/markziemann/background.git

# Set the container working directory
ENV DIRPATH /background
WORKDIR $DIRPATH
