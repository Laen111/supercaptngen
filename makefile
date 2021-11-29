FC=gfortran
FOPT= -O3 -fPIC -std=legacy# -Wall -fbounds-check -g  #legacy is required if you are running gcc 10 or later 
NUMDIR = ./numerical
QAGDIR = ./numerical/dqag
WDIR = ./Wfunctions
RDIR = ./Rfunctions

MAIN = main.o
SUPER = supermain.o

SHAREDOBJ = sharedcap.o
GENOBJ = gencap.o
OPEROBJ = opercap.o
SUPEROBJ = supercap.o

NUMFOBJ =  dgamic.o d1mach.o
QAG = dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o

WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o


gencaplib.so: $(SHAREDOBJ) $(GENOBJ) $(OPEROBJ) $(NUMFOBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $^

# -L tells the linker where to look for shared libraries
# -rpath puts the location of the libraries in the executable so the load can find them at runtime
# -Wl lets us send options to the linker (which are comma seperated)
gentest.x: $(MAIN) gencaplib.so
	${FC} $(FOPT) -L. -Wl,-rpath,. -o $@ $^

supertest.x: $(SUPER) gencaplib.so
	${FC} $(FOPT) -L. -Wl,-rpath,. -o $@ $^

$(NUMFOBJ): %.o : $(NUMDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(SHAREDOBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(GENOBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(OPEROBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(SUPEROBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(SUPER): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(NUMOBJ): %.o: $(NUMDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(QAG): %.o: $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(WFUNC): %.o: $(WDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(RFUNC): %.o: $(RDIR)/%.f
	$(FC) $(FOPT) -c  $<


clean:
	rm -f *.o *.mod *.so *.x
