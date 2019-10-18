FROM ubuntu:16.04
RUN apt-get -m update && apt-get install -y libgsl0-dev python-pip git
RUN pip install numpy scipy
RUN pip install ete2
RUN git init . && git remote add -t \* -f origin https://github.com/morrislab/phylowgs.git && git checkout smchet5
RUN g++ -o mh.o -O3 mh.cpp  util.cpp `gsl-config --cflags --libs`
