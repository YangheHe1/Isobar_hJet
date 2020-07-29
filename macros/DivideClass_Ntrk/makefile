OBJS=DivideClass.o
EXE=analysis

ROOTFLAGS       = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

FASTJETINCFLAGS  = $(shell fastjet-config --cxxflags)
FASTJETLDFLAGS   = $(shell fastjet-config --libs)

FASTJETFLAGS  = $(shell fastjet-config --cxxflags)
FASTJETFLAGS += $(shell fastjet-config --libs)
#FASTJETFLAGS    = `/star/u/nie/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
#FASTJETINCFLAGS = -I/star/u/nie/fastjet-install/include
#FASTJETLDFLAGS  = -L/star/u/nie/fastjet-install/lib #-L/lib/x86_64-linux-gnu/

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CPP = g++ -m32
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS) $(FASTJETINCFLAGS) $(FASTJETLDFLAGS)
#FLAGS = -Wall -g $(FASTJETINCFLAGS) $(FASTJETLDFLAGS)

COMPILE = $(CPP) $(FLAGS) -c

$(EXE): $(OBJS)
	$(CPP) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS) $(FASTJETFLAGS)
	#$(CPP) -o $(EXE) $(OBJS) $(FASTJETFLAGS)

%.o: %.C
	$(COMPILE)  $<
clean :
	rm $(PWD)/$(OBJS)
	rm $(PWD)/$(EXE)

