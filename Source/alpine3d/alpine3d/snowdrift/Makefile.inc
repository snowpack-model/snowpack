#/***********************************************************************************/ 
#/*  Swiss Federal Institute for Snow and Avalanche Research    SLF-DAVOS           */
# -----------------------------------------------------------------------------
#   GENERAL INFO   AUTHORS: PERRY BARTELT and MICHAEL LEHNING    SLF-DAVOS      
# -----------------------------------------------------------------------------
#/**********************************************************************************/ 
#/*                             RESEARCH VERSION 3.0                               */
#/*                    Corresponds to OPERATIONAL VERSION 4.1                      */
#/**********************************************************************************/

######## Sources, objects, headers
SNOWDRIFT_OBJ = $(SNOWDRIFTDIR)/SnowDrift.o \
		$(SNOWDRIFTDIR)/Saltation.o \
		$(SNOWDRIFTDIR)/Suspension.o \
		$(SNOWDRIFTDIR)/Sublimation.o \
		$(SNOWDRIFTDIR)/checksum.o \
		$(SNOWDRIFTDIR)/LegacyIO.o \
		$(SNOWDRIFTDIR)/SnowDriftFEInit.o \
		$(SNOWDRIFTDIR)/SnowDriftFENumerics.o \
		$(SNOWDRIFTDIR)/SnowDriftFEControl.o \
		$(SNOWDRIFTDIR)/Cell.o
SNOWDRIFT_OBJ_POPC  = $(SNOWDRIFTDIR)/SnowDrift.stub.o \
		$(SNOWDRIFTDIR)/PackSnowDrift_par.o \
		$(SNOWDRIFTDIR)/SnowDrift_par.o \
		$(SNOWDRIFTDIR)/Saltation_par.o \
		$(SNOWDRIFTDIR)/Suspension_par.o \
		$(SNOWDRIFTDIR)/Sublimation_par.o \
		$(SNOWDRIFTDIR)/checksum_par.o \
		$(SNOWDRIFTDIR)/LegacyIO_par.o \
		$(SNOWDRIFTDIR)/SnowDriftFEInit_par.o \
		$(SNOWDRIFTDIR)/SnowDriftFENumerics_par.o \
		$(SNOWDRIFTDIR)/SnowDriftFEControl_par.o \
		$(SNOWDRIFTDIR)/Cell_par.o

####### Build rules

snowdrift: $(SNOWDRIFTDIR)/libsnowdrift.a

snowdrift_par: snowdrift_lib_par snowdrift_module_par

snowdrift_lib_par: $(SNOWDRIFTDIR)/libsnowdriftpopc.a

snowdrift_module_par: $(SNOWDRIFTDIR)/snowdrift.module

snowdrift_clean:
	-$(DEL_FILE) $(SNOWDRIFTDIR)/*~ $(SNOWDRIFTDIR)/core $(SNOWDRIFTDIR)/*.core $(SNOWDRIFTDIR)/*.a $(SNOWDRIFTDIR)/*.o $(SNOWDRIFTDIR)/*.module
	-$(DEL_FILE) $(SNOWDRIFTDIR)/_par*

####### Compile
$(SNOWDRIFTDIR)/libsnowdrift.a: $(SNOWDRIFT_OBJ)
	ar -r $@ $(SNOWDRIFT_OBJ)
	ranlib $@

$(SNOWDRIFTDIR)/snowdrift.module: $(SNOWDRIFTDIR)/PackSnowDrift_par.o $(LIBS_POPC)
	$(POPCC) -object -popcld=$(LINKER) -o $@   $(SNOWDRIFTDIR)/PackSnowDrift_par.o $(LDFLAGS) $(LDFLAGS_POPC)

$(SNOWDRIFTDIR)/libsnowdriftpopc.a:  $(SNOWDRIFT_OBJ_POPC)
	ar -r $@ $(SNOWDRIFT_OBJ_POPC)
	ranlib $@

