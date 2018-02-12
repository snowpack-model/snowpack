#/***********************************************************************************/
#/*  Swiss Federal Institute for Snow and Avalanche Research    SLF-DAVOS           */
# -----------------------------------------------------------------------------
#   MAKEFILE AUTHOR: Mathias Bavay    SLF-DAVOS
# -----------------------------------------------------------------------------
####### USER CONFIGURATION
# build mode: release or debug
MODE	 = release

####### DIRECTORIES
MIODIR = ${HOME}/usr
TARGET = runoff_prevah
RUNOFFDIR = .

####### COMPILERS AND OPTIONS
CC		= colorgcc -DGNU -std=c99
CXX		= g++
LINKER		= g++ -DGNU

DEBUG		= -g -O0 -D__DEBUG #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC #-dr
debug_arithm	= -DDEBUG_ARITHM -g
WARNING		=  -Wall  -Wextra #-Weffc++ #-pedantic-errors #-Werror
ARCH		=  -march=native
OPTIM		=  -g -O3 -ftree-vectorize
MAKE_OPTIM	=  -pipe

FF		= gfortran
FFWARNING	=  -Wall
FFARCH		=  -march=native
FFOPTIM		=  -fomit-frame-pointer -O3
FFMAKE_OPTIM	=  -pipe

ifeq ($(MODE),debug)		#for debug
	CFLAGS   = $(debug_arithm) $(MAKE_OPTIM) $(ARCH) $(DEBUG) $(WARNING)
	CCFLAGS  = $(debug_arithm) $(MAKE_OPTIM) $(ARCH) $(DEBUG) $(WARNING)
	FFLAGS	 = $(FFMAKE_OPTIM) $(FFARCH) $(DEBUG) $(FFWARNING)
else				#for release
	CFLAGS   = $(debug_arithm) $(MAKE_OPTIM) $(ARCH) $(OPTIM) $(WARNING)
	CCFLAGS  = $(debug_arithm) $(MAKE_OPTIM) $(ARCH) $(OPTIM) $(WARNING)
	FFLAGS	 = $(FFMAKE_OPTIM) $(FFARCH) $(FFOPTIM) $(FFWARNING)
endif

CFLAGS += -DSVNREV="\"$(SVNREV)\""
CCFLAGS += -DSVNREV="\"$(SVNREV)\""

####### TOOLS
DEL_FILE = rm -f
CP_FILE  = cp -d

####### Libraries
#keep in mind that the order of enumeration of libraries IS important for the linking
#if two libraries depend on each other, the only solution is to list the libraries twice
#(see http://gcc.gnu.org/ml/gcc-help/2004-04/msg00110.html (this is why the popc LDFLAGS are so ugly...)
LIBS          = -lm -lc -ldl -lstdc++ -L/usr/lib
FLIBS         = -lgfortran -lm
LDFLAGS_SEQ   = -rdynamic -lmeteoio $(FLIBS)
LDFLAGS       = -L$(RUNOFFDIR) -L$(MIODIR)/lib $(LIBS)
INCPATH       = -I$(RUNOFFDIR) -I$(MIODIR)/include
CFLAGS   += $(INCPATH)
CCFLAGS  += $(INCPATH)
FFLAGS   += $(INCPATH)

######## Sources, objects, headers
RUNOFF_OBJ = 	$(RUNOFFDIR)/modulerunoff.mod \
		$(RUNOFFDIR)/initrunoff.o \
		$(RUNOFFDIR)/coreralpha.o \
		$(RUNOFFDIR)/helproutines.o \
		$(RUNOFFDIR)/posnam.o \
		$(RUNOFFDIR)/isfile.o \
		$(RUNOFFDIR)/runoffalpha.o \
		$(RUNOFFDIR)/Runoff.o

####### IMPLICIT RULES
%.o : %.cc
	$(CXX) $(CCFLAGS) -c $< -o $@
%.o : %.for
	$(FF) $(FFLAGS) -c $< -o $@
%.mod : %.for
	$(FF) $(FFLAGS) -c $< -o $@

####### GENERAL BUILD RULES

$(TARGET): $(RUNOFF_OBJ) runoff_prevah.o
	$(LINKER) -o $(TARGET) runoff_prevah.o $(LDFLAGS_SEQ) $(LDFLAGS) $(RUNOFF_OBJ) $(FLIBS)

clean:
	-$(DEL_FILE) $(RUNOFFDIR)/*~ $(RUNOFFDIR)/core $(RUNOFFDIR)/*.core $(RUNOFFDIR)/*.a $(RUNOFFDIR)/*.o $(RUNOFFDIR)/*.module $(RUNOFFDIR)/*.mod
	-$(DEL_FILE) $(RUNOFFDIR)/$(TARGET) $(RUNOFFDIR)/$(TARGET).exe
