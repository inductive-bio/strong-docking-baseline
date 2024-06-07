FROM gnina/gnina:v1.1

RUN pip3 install rdkit==2023.9.6 pandas==2.1.4 posebusters==0.2.7

RUN git clone https://github.com/rlabduke/reduce.git && cd reduce && make && make install
