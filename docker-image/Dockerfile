FROM rocker/r-ver:4.4.1

# Install system dependencies for Python
RUN apt-get update && \
    apt-get install -y \
    python3 \
    python3-pip

# Set working directory
WORKDIR /app

# Copy requirements.txt for Python dependencies
COPY requirements.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt

# Install R packages (if needed)
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('optparse')"
RUN R -e "BiocManager::install('R.utils')"
RUN R -e "BiocManager::install('sangerseqR',version='3.20',dependencies = TRUE, force = TRUE)"
RUN R -e "BiocManager::install('sangeranalyseR',version='3.20',dependencies = TRUE, force = TRUE)"