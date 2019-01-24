include ../$(IN)

# MAKEFILE = Makefile_$(CASE)

#  $(LIB)/libchemistry.a: *.f90
#  # 	make -f $^ clean
#  	make -f $(MAKEFILE)
#  	$(LINK) $(OPT) -I. $(INCL) -c Chemistry_Mod.f90
#  	@ar rv $@ *.o
#  	@ar ts $@
#  	@rm -f *.o
#  	@for i in *.mod; do \
#  		mv $$i $(LIB); \
#  	done

SRC = \
	$(CHEMISTRY_CASE)_Precision.f90 \
	$(CHEMISTRY_CASE)_Parameters.f90 \
	$(CHEMISTRY_CASE)_Monitor.f90 \
	$(CHEMISTRY_CASE)_Global.f90 \
	$(CHEMISTRY_CASE)_Initialize.f90 \
	$(CHEMISTRY_CASE)_Function.f90 \
	$(CHEMISTRY_CASE)_JacobianSP.f90 \
	$(CHEMISTRY_CASE)_Jacobian.f90 \
	$(CHEMISTRY_CASE)_LinearAlgebra.f90 \
	$(CHEMISTRY_CASE)_Rates.f90 \
	$(CHEMISTRY_CASE)_Util.f90 \
	$(CHEMISTRY_CASE)_Integrator.f90 \
	$(CHEMISTRY_CASE)_Model.f90 \
	$(CHEMISTRY_CASE)_Main.f90 \
	chemistry_mod.f90

# Object files
OBJ = $(patsubst %.f90, $(LIB)/%.o, $(SRC))

all: $(OBJ)

# -J: where to put '.mod' files
$(LIB)/%.o: $(CHEMISTRY_CASE)/%.f90
	$(FC) $(FFLAGS) -I$(LIB) -J$(LIB) -c $^ -o $@
