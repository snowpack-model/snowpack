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
EB_OBJ 	  =	$(EBDIR)/EnergyBalance.o \
		$(EBDIR)/radiationfield.o

EB_OBJ_POPC  = $(EBDIR)/EnergyBalance_par.o \
		$(EBDIR)/radiationfield.stub.o \
		$(EBDIR)/radiationfield_par.o \
		$(EBDIR)/EnergyBalance.stub.o \
		$(EBDIR)/PackEB_par.o
#		$(EBDIR)/EnergyBalance.stub.o \
#		$(SNOWPACKDIR)/PackSnowPack_par.o \

####### Build rules

ebalance: $(EBDIR)/libebalance.a

ebalance_par: ebalance_lib_par ebalance_module_par

ebalance_lib_par: $(EBDIR)/libebalancepopc.a

ebalance_module_par: $(EBDIR)/ebalance.module

ebalance_clean:
	-$(DEL_FILE) $(EBDIR)/*~ $(EBDIR)/core $(EBDIR)/*.core $(EBDIR)/*.a $(EBDIR)/*.o $(EBDIR)/*.module
	-$(DEL_FILE) $(EBDIR)/_par*

####### Compile
$(EBDIR)/libebalance.a: $(EB_OBJ)
	ar -r $@ $(EB_OBJ)
	ranlib $@

$(EBDIR)/libebalancepopc.a: $(EB_OBJ_POPC)
	ar -r $@ $(EB_OBJ_POPC)
	ranlib $@

$(EBDIR)/ebalance.module: $(EBDIR)/PackEB_par.o $(LIBS_POPC)
	$(POPCC) -object -popcld=$(LINKER) $(CCFLAGS) -o $@ $(EBDIR)/PackEB_par.o $(LDFLAGS_POPC) $(LDFLAGS)
