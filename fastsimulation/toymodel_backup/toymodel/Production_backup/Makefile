DEPENDROOT    = $(shell [ -n "${ROOTSYS}" ] && echo ${ROOTSYS} || echo "ROOTSYS_NOTSET" )

#include $(DEPENDROOT)/test/Makefile
include ../dev/Makefile.arch

DEPENDFJ      = $(shell [ -n "${FASTJETDIR}" ] && echo ${FASTJETDIR} || echo "FASTJET_NOTSET" )

DEBUG         = -g

MODULE	      = Thrm

#DEPENDENCIES  =	$(DEPENDFJ)/include/FjPseudoJet.hh \
#					$(DEPENDFJ)/include/fastjet/SISConePlugin.hh \

DEPENDENCIES  =	$(DEPENDFJ)/include/fastjet/PseudoJet.hh \
					$(DEPENDFJ)/include/fastjet/ClusterSequence.hh \
					$(DEPENDFJ)/include/fastjet/ClusterSequenceArea.hh \
					$(DEPENDFJ)/include/fastjet/ClusterSequenceActiveArea.hh \
					$(DEPENDFJ)/include/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh \
					$(DEPENDFJ)/include/fastjet/Selector.hh \
					$(DEPENDFJ)/include/fastjet/config.h \
					$(DEPENDFJ)/include/fastjet/tools/Subtractor.hh \
					$(DEPENDFJ)/include/fastjet/tools/JetMedianBackgroundEstimator.hh \

INCLUDES      = $(shell ${DEPENDFJ}/bin/fastjet-config --cxxflags) \

CXXFLAGS     += $(INCLUDES) -Wall -c $(DEBUG)

EXTRAROOTLIBS = -lEG -lEGPythia6 -lProof -lProofPlayer -lMinuit -lGeom -lVMC

LDFLAGS      += $(ROOTLIBS) $(EXTRAROOTLIBS) \
		$(shell ${DEPENDFJ}/bin/fastjet-config --libs --plugins) \
		-Wall $(DEBUG)

MODDICTSBASE  = $(MODULE)RootDict
MODDICTS      = $(MODDICTSBASE).$(SrcSuf)
MODDICTO      = $(MODDICTSBASE).$(ObjSuf)

LINKDEFFILE   = $(MODULE)LinkDef.h

# here put all the files for which cint should be called (the generation of dictionaries
# remember to add a correct line to the $(MODULE)LinkDef.h file for each class or funtion or namespace
# use the wild card but then the file name should start with $(MODULE)
SRCTMP        = $(wildcard $(MODULE)*.h)

SRCTMPC       = $(patsubst $(MODDICTSBASE).h,,$(SRCTMP))
MODHEADERS    = $(patsubst $(MODULE)LinkDef.h,,$(SRCTMPC))
CINTHEADERS   = $(patsubst $(MODULE)LinkDef.h,,$(SRCTMPC))

#here the files for which we will NOT call cint
MODHEADERS   += dNdpT.h JetReconstruction.h #fjwrapper.h 

MODSRCS	      = $(patsubst %.h,%.$(SrcSuf),$(MODHEADERS)) \
		$(MODDICTS)

MODO          = $(patsubst %.h,%.$(ObjSuf),$(MODHEADERS)) \
                $(MODDICTO)

OPT	      = -g
OPT2	      = -g

MODSO         = lib$(MODULE).$(DllSuf)

ifeq ($(PLATFORM),win32)
MODLIB        = lib$(MODULE).lib
else
MODLIB        = $(shell pwd)/$(MODSO)
endif

OBJS          = $(MODO)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(DEPENDENCIES) $(MODSO)

$(MODSO):	$(MODO)

ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(MODLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
		$(MT_DLL)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "$@ done"

clean:
		@rm -rf $(OBJS) $(MODSO) $(MODDICTS) $(MODDICTSBASE).h \
		core DataOut *~

distclean:      clean
		-@mv -f linearIO.root linearIO.roott
		@rm -f $(MODSO) $(MODLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations

.SUFFIXES: .$(SrcSuf)

$(MODDICTS): $(CINTHEADERS) $(LINKDEFFILE)
	@echo "[i] Generating dictionary $@..."
	@rootcint -f $@ -c -I$(DEPENDPICO) -I$(DEPENDPICOPY) $^

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

