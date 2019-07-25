FC1 = gfortran ./source/prep.f -o run && ./run && rm run

FC2 = gfortran -O3 -mcmodel=large -funroll-loops

OBJ =                        \
./source/dirhbz.f            \
./source/mait_multipole.f    \
./source/init_spurious.f     \
./source/init_coulomb.f      \
./source/init_pairing.f      \
./source/start_fam.f         \
./source/iter_fam.f          \
./source/fam_drhodkappa.f    \
./source/fam_ddensdcurr.f    \
./source/fam_dpotentials.f   \
./source/fam_dcoulomb.f      \
./source/fam_dh1.f           \
./source/fam_ddelta.f        \
./source/fam_broyden.f       \
./source/fam_strength.f      \
./source/printout.f          \
./source/utility_functions.f \

prep:
	$(FC1)

run: $(OBJ)
	$(FC2) $(OBJ) -lblas -llapack -o run
