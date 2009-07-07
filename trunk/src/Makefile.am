# Uncomment (i.e. remove the `#' character at the begining of) 
# the two lines below if you want to use MPI.

# CC=mpicc
# DEFS=-DUNIX -D$(PROG) -DDEBUG -DMPI

# Comment the line below if you want to use MPI.
# Uncomment the line below if you don't want to use MPI.

DEFS=-DUNIX -D$(PROG) -DDEBUG

if WANT_MC
bin_PROGRAMS = mc
PROG = MC
else
if WANT_PART
bin_PROGRAMS = part
PROG = PART
else
bin_PROGRAMS = phyml
PROG = PHYML
endif
endif

if WANT_MC
mc_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
options.c options.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
spr.c spr.h\
mc.c mc.h\
m4.c m4.h\
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h
mc_LDADD = -lm
else 
if WANT_PART
part_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
options.c options.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
spr.c spr.h\
mc.c mc.h\
m4.c m4.h\
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h
part_LDADD = -lm
else
phyml_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
options.c  options.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
mpi_boot.c  mpi_boot.h\
stats.c stats.h
phyml_LDADD = -lm
endif
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


