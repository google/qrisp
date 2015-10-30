# Makefile for the QRisp project.
VPATH=third_party

LIBS=-lglog -lprotobuf -lgflags

PROTOC=protoc
PROTOCFLAGS=--cpp_out=third_party

CXX=g++
CXXFLAGS=--std=c++11

SOURCES=$(wildcard *.cc)
HEADERS=$(wildcard *.h)
PROTOS=$(wildcard *.proto)

OBJECTS=$(SOURCES:%.cc=%.o)

PBS=$(PROTOS:%.proto=%.pb)

qrisp: protos ${OBJECTS}
	${CXX} ${CXXFLAGS} third_party/qrisp.cc -o qrisp ${OBJECTS} ${LIBS}

protos: ${PBS}
	@ echo ${PROTOS}

.cc.o:
	$(CXX) -c $(CXXFLAGS) -o $@ $<

%.pb: %.proto
	${PROTOC} ${PROTOCFLAGS} $*.proto
	${CXX} ${CXXFLAGS} -c -o $*.pb.o $*.pb.cc

clean:
	rm -f *.o
