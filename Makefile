FC = gfortran

FFLAGS = -O3 -fopenmp
LFLAGS = -fopenmp ~/ARPACK/libarpack_LINUX.a -llapack -lblas

SRCF = 	src/dneupd2.f src/hp_sort2.f
SRCF90 = src/grid_construct.f90 \
		src/potential_module.f90 \
		src/mapping_parameters.f90 \
		src/index_generator.f90 \
		src/index_mapping.f90 \
		src/matrix_vector_products.f90 \
		src/basis_module.f90 \
		src/matrix_construction.f90 \
		src/arpack_module.f90 \
		src/output_module.f90 \
		src/coll.f90

OBJ = $(patsubst src/%.f,build/%.o,$(SRCF)) \
      $(patsubst src/%.f90,build/%.o,$(SRCF90))


build/%.o: src/%.f
	$(FC) $(FFLAGS) -c $< -o $@

build/%.o: src/%.f90
	$(FC) $(FFLAGS) -Jbuild -Ibuild -c $< -o $@


coll: $(OBJ)
	$(FC) $(OBJ) -o $@ $(LFLAGS)

clean:
	rm -f build/*.o build/*.mod coll


!OMP_NUM_THREADS=n ./coll