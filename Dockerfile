FROM python:3.11-alpine

RUN apk add make gcc musl-dev zlib-dev

RUN mkdir -p /bcftools && cd /bcftools && wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && tar -jxvf bcftools-1.21.tar.bz2 && cd bcftools-1.21 && make

COPY ./requirements.txt /requirements.txt
RUN python -m pip install -r /requirements.txt

COPY ./gwas-indexing /gwas-indexing

ENTRYPOINT ["python", "/gwas-indexing/main.py"]
