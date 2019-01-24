IN = M_DEF
include $(IN)

# OBJ = $(CHEMISTRY_LIB)

all: main.exe

main.exe: chemistry_obj
	@ar rv $(LIB)/liball.a $(LIB)/*.o
	$(FC) $(FFLAGS) -I$(LIB) main/main.f90 -L$(LIB) -lall -o $@

chemistry_obj:
	@cd chemistry; make "IN=$(IN)"

clean:
	rm -f $(LIB)/*.o
	rm -f $(LIB)/*.mod
# 
# cleanall:
# 	rm -f $(LIB)/*.a
# 	rm -f $(LIB)/*.mod
# 	@cd Chemistry; make "IN=$(IN)" clean
