IN = M_DEF
include $(IN)

# METEOROLOGY_LIB = -L$(LIB) -lmeteorology
# CHEMISTRY_LIB = -L$(LIB) -lchemistry
# CHEMISTRY_FILES = \
# 	$(CHEMISTRY_CASE)_Precision.f90 \
# 	$(CHEMISTRY_CASE)_Parameters.f90 \
# 	$(CHEMISTRY_CASE)_Monitor.f90 \
# 	$(CHEMISTRY_CASE)_Global.f90 \
# 	$(CHEMISTRY_CASE)_Initialize.f90 \
# 	$(CHEMISTRY_CASE)_Function.f90 \
# 	$(CHEMISTRY_CASE)_Stochastic.f90 \
# 	$(CHEMISTRY_CASE)_JacobianSP.f90 \
# 	$(CHEMISTRY_CASE)_Jacobian.f90 \
# 	$(CHEMISTRY_CASE)_LinearAlgebra.f90 \
# 	$(CHEMISTRY_CASE)_Rates.f90 \
# 	$(CHEMISTRY_CASE)_HessianSP.f90 \
# 	$(CHEMISTRY_CASE)_Hessian.f90 \
# 	$(CHEMISTRY_CASE)_StoichiomSP.f90 \
# 	$(CHEMISTRY_CASE)_Stoichiom.f90 \
# 	$(CHEMISTRY_CASE)_Util.f90 \
# 	$(CHEMISTRY_CASE)_Model.f90 \
# 	$(CHEMISTRY_CASE)_Integrator.f90 \
# 	Chemistry_Mod.f90 \
# 
# CHEMISTRY_SRC = $(addprefix $(CHEMISTRY_DIR), $(CHEMISTRY_FILES))
# CHEMISTRY_OBJ = $($(addprefix $(LIB), $(CHEMISTRY_FILES)):.f90=.o)





# OBJ = $(CHEMISTRY_LIB)

all: main.exe

main.exe: chemistry_obj
# 	$(FC) $(FFLAGS) -I$(LIB) MAIN/main.f90 -o $@
# 	$(FC) $(FFLAGA) -I$(LIB) $(LIB)/main.o
	@ar rv $(LIB)/liball.a $(LIB)/*.o
	$(FC) $(FFLAGS) -I$(LIB) MAIN/main.f90 -L$(LIB) -lall -o $@

chemistry_obj:
	@cd CHEMISTRY; make "IN=$(IN)"

# 
# # meteorology:
# # 	@cd Meteorology; make "IN=$(IN)"
# 
# chemistry_obj:
# # 	@cd Chemistry; make "IN=$(IN)"
# 
clean:
	rm -f $(LIB)/*.o
	rm -f $(LIB)/*.mod
# 
# cleanall:
# 	rm -f $(LIB)/*.a
# 	rm -f $(LIB)/*.mod
# 	@cd Chemistry; make "IN=$(IN)" clean
# 
# 
# all: main.exe
# 
# main.exe: $(OBJS)
# 
# 
# 
# $(OBJDIR)/%.o: $(CHEMISTRY_DIR)/%.f90
# 	$(FC) $(FFLAGS) -c 
