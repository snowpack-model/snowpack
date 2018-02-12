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
RUNOFF_OBJ 	  =	$(RUNOFFDIR)/Runoff.o
#RUNOFF_OBJ_POPC  = $(COMMONDIR)/marshal_Alpine3D_par.o

####### Build rules

runoff: $(RUNOFFDIR)/librunoff.a

#runoff_par:  $(RUNOFF_OBJ_POPC) runoff_lib_par

#runoff_lib_par: $(RUNOFFDIR)/librunoffpopc.a

#runoff_module_par: $(RUNOFFDIR)/libcommonpopc.a

runoff_clean:
	-$(DEL_FILE) $(RUNOFFDIR)/*~ $(RUNOFFDIR)/core $(RUNOFFDIR)/*.core $(RUNOFFDIR)/*.a $(RUNOFFDIR)/*.o $(RUNOFFDIR)/*.module
	-$(DEL_FILE) $(RUNOFFDIR)/_par*

####### Compile
$(RUNOFFDIR)/librunoff.a: $(RUNOFF_OBJ)
	ar -r $@ $(RUNOFF_OBJ)
	ranlib $@

#$(RUNOFFDIR)/libcommonpopc.a: $(RUNOFF_OBJ_POPC) $(RUNOFF_OBJ)
#	ar -r $@  $(RUNOFF_OBJ_POPC) $(RUNOFF_OBJ)
#	ranlib $@

#$(RUNOFF_OBJ_POPC): %_par.o : %.cc
#	$(POPCC) $(CCFLAGS) -c $< -o $@

