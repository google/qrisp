# Makefile for the QRisp project.
#VPATH=third_party

LIBS=-lglog -lprotobuf -lgflags

PROTOC=protoc
PROTOCFLAGS=--cpp_out=.

CXX=g++
CXXFLAGS=--std=c++11

SOURCES=$(wildcard *.cc)
HEADERS=$(wildcard *.h)

PROTOS=$(wildcard *.proto)
PROTO_OBJS=$(PROTOS:.proto=.pb.o)

OBJS=$(SOURCES:%.cc=%.o)

PBSRCS   := $(wildcard *.proto)
PBOBJS   := $(PROTOS:.proto=.pb.o)
PBGENS   := $(PROTOS:.proto=.pb.cc) $(PROTOS:.proto=.pb.h)

.PRECIOUS: ${PBGENS}

all: qrisp

%.pb.cc: %.proto
	${PROTOC} ${PROTOCFLAGS} $<

%.pb.o : %.pb.cc
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

qrisp: ${PBOBJS} ${OBJECTS}
	${CXX} ${CXXFLAGS} third_party/qrisp.cc -o qrisp ${OBJECTS} ${LIBS}

clean:
	rm -f *.o
