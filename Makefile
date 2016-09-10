.PHONY: all

#  These objects do not all compile from the source provided by Allen and
#  Tildesley since some files contain duplicate definitions of the same
#  subroutine. These samples were originally crafted as demonstrative
#  snippets, not standalone code.

OBJECTS=alpha-fcc.o avoid-sqrt.o brownian.o clustering.o constraint.o ewald.o \
        fft.o gear.o hardsphere.o init-velocity.o leapfrog-rotation.o \
        leapfrog.o link-cell-sheared.o link-cell.o lj.o mc-dumbbell.o \
        mc-hardlines.o mc-muvt-indices.o mc-muvt.o mc-npt.o mc-nvt.o \
        md-multitimestep.o md-nph-extended.o md-npt-constraint.o \
        md-nvt-constraint.o md-nvt-extended.o order-parameter.o pbc.o \
        quaternion.o random-rotate.o rattle.o shake.o time-correlation.o \
        unfold-pbc.o velocityverlet.o verlet-list.o voronai.o

all: $(OBJECTS)

%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@
