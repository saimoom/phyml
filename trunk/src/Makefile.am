
# Uncomment (i.e. remove the `#' character at the begining of) 
# the two lines below if you want to use MPI.
# Comment the two lines below if you don't want to use MPI.

# CC=mpicc
# DEFS=-DUNIX -D$(PROG) -DDEBUG -DMPI

# Comment the line below if you want to use MPI.
# Uncomment the line below if you don't want to use MPI.





if WANT_MC
bin_PROGRAMS = mc
PROG = MC
else
bin_PROGRAMS = phyml
PROG = PHYML
endif

DEFS=-DUNIX -D$(PROG) -DDEBUG

if WANT_MC
mc_SOURCES = main.c utilities.c optimiz.c lk.c bionj.c models.c free.c\
options.c simu.c eigen.c pars.c alrt.c interface.c cl.c spr.c mc.c m4.c draw.c rates.c mcmc.c numeric.c
mc_LDADD = -lm
else
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
endif

# all-am:	intro bin_PROGRAMS$(EXEEXT)
# 	@echo ""
# 	@echo "Done."

# intro:	
# 	@echo ""
# 	@echo ""
# 	@echo ":: Building PhyML version $(VERSION) ::"
# 	@echo ""
# 	@echo ""


