FROM gnina/gnina:v1.1@sha256:acf03094a27a7904d75e5c70f5e37b07d7f70372b2f56cfd4989b47b9c310fc0

RUN pip3 install rdkit==2023.9.6 pandas==2.1.4 posebusters==0.2.7

RUN git clone https://github.com/rlabduke/reduce.git && cd reduce && make && make install
