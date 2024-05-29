FC=gfortran
FFLAGS= -O3 -Wall -fopenmp

.SUFFIXES: .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

OBJECTS=\
   su2main.o\
   su2lattice.o\
   su2thermalize.o\
   su2measures.o\
   su2center.o\
   su2gaugefix.o\

SRC=\
   su2main.f90\
   su2lattice.f90\
   su2thermalize.f90\
   su2measures.f90\
   su2center.f90\
   su2gaugefix.f90\

su2main: $(OBJECTS)

	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm *\.o su2main $(OBJECTS)

ctags: $(SRC)
	ctags $(SRC)

su2main.o: su2lattice.o su2thermalize.o su2measures.o su2gaugefix.o su2center.o su2main.f90
	$(FC) $(FFLAGS) -c su2main.f90 -o su2main.o 
	
su2center.o: su2lattice.o su2measures.o su2center.f90
	$(FC) $(FFLAGS) -c su2center.f90 -o su2center.o 
	
su2thermalize.o: su2lattice.o su2thermalize.f90
	$(FC) $(FFLAGS) -c su2thermalize.f90 -o su2thermalize.o 

su2gaugefix.o: su2lattice.o su2gaugefix.f90
	$(FC) $(FFLAGS) -c su2gaugefix.f90 -o su2gaugefix.o 

su2measures.o: su2lattice.o su2measures.f90
	$(FC) $(FFLAGS) -c su2measures.f90 -o su2measures.o
	
su2lattice.o: su2lattice.f90
	$(FC) $(FFLAGS) -c su2lattice.f90 -o su2lattice.o 