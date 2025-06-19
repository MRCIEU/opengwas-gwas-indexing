FROM python:3.11-alpine

RUN apk add make gcc musl-dev zlib-dev xz-dev bzip2-dev curl-dev

RUN mkdir -p /bcftools && cd /bcftools && wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && tar -jxvf bcftools-1.21.tar.bz2 && cd bcftools-1.21 && make

RUN mkdir -p /plink && cd /plink && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip && unzip -q plink_linux_x86_64_20241022.zip plink

RUN mkdir -p /ref && cd /ref && wget http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz && tar -zxvf data_maf0.01_rs_ref.tgz && rm data_maf0.01_rs_ref.tgz

COPY ./requirements.txt /requirements.txt
RUN python -m pip install -r /requirements.txt

RUN apk add autossh openssh-client

COPY ./gwas-indexing _oci.py start.sh /gwas-indexing/

RUN chmod +x /gwas-indexing/start.sh
CMD ["/gwas-indexing/start.sh"]
