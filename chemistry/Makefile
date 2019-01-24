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
	$(CHEMCASE)_Precision.f90 \
	$(CHEMCASE)_Parameters.f90 \
	$(CHEMCASE)_Monitor.f90 \
	$(CHEMCASE)_Global.f90 \
	$(CHEMCASE)_Initialize.f90 \
	$(CHEMCASE)_Function.f90 \
	$(CHEMCASE)_JacobianSP.f90 \
	$(CHEMCASE)_Jacobian.f90 \
	$(CHEMCASE)_LinearAlgebra.f90 \
	$(CHEMCASE)_Rates.f90 \
	$(CHEMCASE)_Util.f90 \
	$(CHEMCASE)_Integrator.f90 \
	$(CHEMCASE)_Model.f90 \
	$(CHEMCASE)_Main.f90 \
	chemistry_mod.f90

# Object files
OBJ = $(patsubst %.f90, $(LIB)/%.o, $(SRC))

all: $(OBJ)

# -J: where to put '.mod' files
$(LIB)/%.o: $(CHEMDIR)/$(CHEMCASE)/%.f90
	$(FC) $(FFLAGS) -I$(LIB) -J$(LIB) -c $^ -o $@
