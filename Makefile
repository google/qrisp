
#$HOME/google_libs/bin/protoc *.proto --cpp_out=./
VPATH=third_party:third_party/proto

LIB_PATHS=/usr/local/google/home/fdb/google_libs/lib
LIBS=-lglog -lprotobuf -lgflags
INCL=/usr/local/google/home/fdb/google_libs/include

CC=g++
CPPFLAGS=--std=c++11

SRC=cluster.cc\
dataset-utils.cc\
learning-utils.cc\
model.cc\
performance.cc\
plif.cc\
recurrences-nbest.cc\
recurrences.cc\
rna-structure.cc\
sgd.cc\
utils.cc\
config.pb.cc\
parameters.pb.cc\
structure.pb.cc

OBJS=cluster.o\
dataset-utils.o\
learning-utils.o\
model.o\
performance.o\
plif.o\
recurrences-nbest.o\
recurrences.o\
rna-structure.o\
sgd.o\
utils.o\
config.pb.o\
parameters.pb.o\
structure.pb.o

qrisp: ${OBJS}
	${CC} ${CPPFLAGS} third_party/qrisp-train.cc -o qrisp-train ${OBJS} -I${INCL} -L${LIB_PATHS} ${LIBS}

.cc.o:
	$(CC) -c $(CPPFLAGS) -I${INCL} -o $@ $<

clean:
	rm -f *.o
