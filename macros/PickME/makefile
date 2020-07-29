
#CXX = g++
CXX = $(shell root-config --cxx)
ROOTFLAGS = $(shell root-config --cflags)

CFLAGS = ${ROOTFLAGS} -g -IStRoot/MyAnalysis

ROOTLIBS = $(shell root-config --libs)


all: pickME.bin

pickME.bin: pickME.o dict.o mcTrackEvent.o mcTrackParticle.o 
	${CXX} ${CFLAGS} $^ -o$@ ${ROOTLIBS}

dict.C: StRoot/MyAnalysis/mcTrackEvent.h StRoot/MyAnalysis/mcTrackParticle.h
	rootcint -f $@ -c -p StRoot/MyAnalysis/mcTrackEvent.h StRoot/MyAnalysis/mcTrackParticle.h

dict.o: dict.C
	${CXX} -fPIC  -c ${CFLAGS} $< -o $@

pickME.o: pickME.C
	${CXX} -fPIC  -c ${CFLAGS} $< -o $@

%.o: StRoot/MyAnalysis/%.cxx StRoot/MyAnalysis/%.h
	${CXX} -fPIC -c ${CFLAGS} $< -o $@

%.so: StRoot/MyAnalysis/%.cxx StRoot/MyAnalysis/%.h
	rootcint -f 

clean:
	rm -vf *.so *.o *.d pickME.bin dict.C

#g++ -o pickME.bin pickME.C StRoot/MyAnalysis/ITrack.cxx StRoot/MyAnalysis/JEvent.cxx StRoot/MyAnalysis/mcTrackParticle.cxx -I/usr/local/root/include -I/StRoot/MyAnalysis -L/usr/local/root/lib -lCore -lRIO -lTree -lMathCore -lPhysics -lHist
#	echo "Starting to compile."
#g++ -c StRoot/MyAnalysis/ITrack.cxx  -I/usr/local/root/include -I/StRoot/MyAnalysis `root-config --libs`
#echo "Finished compiling "
#g++ -c StRoot/MyAnalysis/JEvent.cxx -I/usr/local/root/include -I/StRoot/MyAnalysis `root-config --libs`
#echo "Finished compiling JEvent.o"
#g++ -c StRoot/MyAnalysis/mcTrackParticle.cxx -I/usr/local/root/include -I/StRoot/MyAnalysis `root-config --libs`
#echo "Finished compiling mcTrackParticle.o"
#g++ -g -o pickME.bin pickME.C JEvent.o mcTrackParticle.o  -I/usr/local/root/include -I/StRoot/MyAnalysis `root-config --libs`
#echo "Everything is compiled now???"





