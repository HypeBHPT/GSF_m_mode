BASE = $(shell /bin/pwd)
EXED = $(BASE)/$(EXECTBL)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src
INCD = $(BASE)/header


#Compiler of C
CC=gcc

#Flags -I tells which directory to look for the .h files 
CFLAGS =-Wall -I$(INCD) -O3 -I /opt/local/include/libomp -I /opt/local/include/gsl  -L /opt/local/lib/libomp -L/opt/local/lib -I/opt/local/include

#Libraries (Math)
LIBS=-Xpreprocessor -fopenmp -lomp -lgsl -lgslcblas -lmpfr -lgmp

#set of headers .h files
DEPS = $(INCD)/$(wildcard *.h)


#set of source .c files
SRC = kerr-circular.c load_data.c test_routines.c hyperboloidal_functions.c SelfForce_Functions.c parameters.c Coordinates.c IndexRoutines.c solve_equations.c derivatives.c FieldEqs.c Allocate.c get_InitialGuess.c debug.c newton.c newton_direct.c output.c Solve_PDE.c spectral_utilities.c utilities.c FuncAndJacobian.c
  
#set of object .o files
OBJS = $(SRC:.c=.o)

EXECTBL  = SolveODE

#rule to say that .o files depends on .c and files on DEPS. All those
#files are going to be checked on the execution for the up-to-date version


all:$(EXECTBL)


$(EXECTBL): $(OBJS)
	@echo '' 
	@echo '============ ...building the executable ============'
	@echo ''
	cd $(OBJD); $(CC) $(CFLAGS) -o ../$(EXECTBL) $(OBJS) $(LIBS)
	@echo ''
	@echo '============ ...Done! ============'


vpath %.o $(OBJD)
vpath %.c $(SRCD)
%.o:%.c $(DEPS)
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	cd $(OBJD); $(CC) -o $@ -c $(SRCD)/$(*).c $(LIBS) $(CFLAGS)



flush:
	rm -f $(EXECTBL) *.dat fort.* $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc*

clean:
	rm -f *.dat *.yg fort.* $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc* $(SRCD)/*.o $(EXECTBL)
	rm -f *~ $(OBJD)/*~ $(SRCD)/*~ header/*~
