
# Uncomment (i.e. remove the `#' character at the begining of) 
# the two lines below if you want to use MPI.
# Comment the two lines below if you don't want to use MPI.

# CC=mpicc
# DEFS=-DUNIX -D$(PROG) -DDEBUG -DMPI

# Comment the line below if you want to use MPI.
# Uncomment the line below if you don't want to use MPI.


DEFS=-DUNIX -D$(PROG) -DDEBUG

bin_PROGRAMS = phyml
EXTRA_PROGRAMS = m4 mg mc
PROG = PHYML


phyml_SOURCES = main.c \
     utilities.c utilities.h \
     optimiz.c optimiz.h \
     lk.c lk.h \
     bionj.c bionj.h \
     models.c models.h \
     free.c free.h \
     options.c options.h \
     simu.c simu.h \
     eigen.c eigen.h \
     pars.c pars.h \
     alrt.c alrt.h \
     interface.c interface.h \
     cl.c cl.h \
     spr.c spr.h \
     draw.c draw.h \
     mpi_boot.c mpi_boot.h \
     numeric.c numeric.h

phyml_LDADD = -lm


all-am:	intro phyml$(EXEEXT)
	@echo ""
	@echo "Done."

intro:	
	@echo ""
	@echo ""
	@echo ":: Building PhyML version $(VERSION) ::"
	@echo ""
	@echo ""


