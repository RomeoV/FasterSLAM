FROM gcc:9

RUN apt-get update
RUN apt-get install -y cmake valgrind python-pip
RUN pip install gcovr
