FROM python:3.12-slim

LABEL org.opencontainers.image.title="bio-oracle-replication"
LABEL org.opencontainers.image.description="Replication of Rutterford et al. (2023): SST as primary driver of fish community structure, tested in the Mediterranean using ClimateFish + Bio-ORACLE"
LABEL org.opencontainers.image.authors="Anne Fouilloux <anne.fouilloux@gmail.com>"
LABEL org.opencontainers.image.source="https://github.com/annefou/bio-oracle"
LABEL org.opencontainers.image.licenses="MIT"

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY 01_download_climatefish.py .
COPY 02_download_biooracle.py .
COPY 03_spatial_join.py .
COPY 04_community_analysis.py .
COPY 05_variable_importance.py .

RUN mkdir -p data results

CMD ["sh", "-c", "python 01_download_climatefish.py && python 02_download_biooracle.py && python 03_spatial_join.py && python 04_community_analysis.py && python 05_variable_importance.py"]
