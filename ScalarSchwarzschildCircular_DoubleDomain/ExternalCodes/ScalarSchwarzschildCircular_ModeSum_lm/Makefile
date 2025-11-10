BASE = $(shell /bin/pwd)
EXED = $(BASE)/$(EXECTBL)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src
INCD = $(BASE)/header


#Compiler of C
CC=gcc

#Flags -I tells which directory to look for the .h files 
CFLAGS =-Wall -I$(INCD) -O3 -I /opt/local/include/libomp  -L /opt/local/lib/libomp

#Libraries (Math)
LIBS=-Xpreprocessor -fopenmp -lomp

#set of headers .h files
DEPS = $(INCD)/$(wildcard *.h)


#set of source .c files
SRC = Field_Eqs.c parameters.c hyperboloidal_functions.c SelfForce_Functions.c Solve_ODE.c FuncAndJacobian.c Coordinates.c Allocate.c derivatives.c InitialGuess.c IndexRoutines.c debug.c    newton.c  output.c  spectral_utilities.c  utilities.c
  
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
