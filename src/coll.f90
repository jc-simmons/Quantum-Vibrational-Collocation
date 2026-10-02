program main

    use iso_fortran_env, only: real64
    use grid_module
    use potential_module
    use mapping_parameters
    use index_generator_module
    use index_mapping_module
    use mvp_module
    use basis_keo_module
    use matrix_module
    use arpack_module
    use output_module

    implicit none

    !$use omp_lib

    integer :: ndim, keo_terms
    integer :: nb, np
    integer :: info
    integer :: basis_max_ind, grid_max_ind
    integer :: config_unit
    integer :: time1, time2, rate
    integer :: ido
    integer :: iparam(11)

    integer, allocatable :: basis_limits(:)
    integer, allocatable :: grid_limits(:)
    integer, allocatable :: basis_inds(:,:)
    integer, allocatable :: grid_inds(:,:)
    integer, allocatable :: limits(:)
    integer, allocatable :: n(:)
    integer, allocatable :: mappings(:,:,:,:,:)
    integer, allocatable :: sortlens(:,:)
    integer, allocatable :: sortlims(:,:,:,:)

    real(real64), allocatable :: Bz(:,:,:,:,:,:)
    real(real64), allocatable :: Binv(:,:,:,:)
    real(real64), allocatable :: potp(:)
    real(real64), allocatable :: grid(:,:)
    real(real64), allocatable :: dr(:)
    real(real64), allocatable :: di(:)
    real(real64), allocatable :: zeig(:,:)

    real(real64) :: ran

    namelist /system/ ndim, keo_terms
    namelist /indices/ basis_limits, grid_limits

    !$ call omp_set_num_threads(1)

    call system_clock(time1)
    call system_clock(count_rate=rate)

    call random_seed()
    call random_number(ran)

    open(newunit=config_unit, file='config.nml', status='old', action='read')

    read(config_unit, nml=system)

    allocate(basis_limits(ndim))
    allocate(grid_limits(ndim))
    allocate(n(ndim))
  
    rewind(config_unit)
    read(config_unit, nml=indices)
    rewind(config_unit)
    read(config_unit, nml=basis_config)
    rewind(config_unit)
    read(config_unit, nml=grid_config)
    rewind(config_unit)
    read(config_unit, nml=potential_config)
    rewind(config_unit)
    read(config_unit, nml=output_config)

    close(config_unit)

    call system_clock(time1, rate)

    basis_max_ind = maxval(basis_limits)
    grid_max_ind  = maxval(grid_limits)

    allocate(limits(ndim))

    limits = basis_max_ind
    call generate_indices(limits, basis_inds)
    nb = size(basis_inds, 2)

    limits = grid_max_ind
    call generate_indices(limits, grid_inds)
    np = size(grid_inds, 2)

    n = grid_max_ind

    print*, 'basis/grid size: ', nb, np

    basis_max_ind = basis_max_ind + 1
    grid_max_ind  = grid_max_ind + 1


    call compute_grid(grid, ndim, grid_max_ind)

    call compute_potential(grid, grid_inds, potp)


    allocate( &
        Binv(ndim,NUM_LU_FACTORS,grid_max_ind,grid_max_ind), &
        Bz(NUM_B_MATRICES, NUM_LU_FACTORS,ndim,keo_terms, &
        grid_max_ind,grid_max_ind) &
    )

    call construct_matrices(grid, n, ndim, hob, hob_keo, Bz, Binv)

    allocate(mappings(NUM_MAP_TYPES,NUM_STATES,ndim,np,grid_max_ind))
    allocate(sortlens(NUM_MAP_TYPES,ndim))
    allocate(sortlims(NUM_MAP_TYPES,NUM_STATES,ndim,np))

    mappings = 0
    sortlens = 0
    sortlims = 0

    call build_space_map(basis_inds, basis_max_ind, grid_max_ind, &
        .false., mappings(BASIS_TO_POINTS,:,:,:,:), sortlens, sortlims, &
        BASIS_TO_POINTS, ndim)

    call build_space_map(grid_inds, grid_max_ind, basis_max_ind, &
        .false., mappings(POINTS_TO_BASIS,:,:,:,:), sortlens, sortlims, &
        POINTS_TO_BASIS, ndim)

    call build_space_map(grid_inds, grid_max_ind, grid_max_ind, &
        .false., mappings(POINTS_TO_POINTS,:,:,:,:), sortlens, sortlims, &
        POINTS_TO_POINTS, ndim)

    call build_space_map(basis_inds, basis_max_ind, basis_max_ind, &
        .false., mappings(BASIS_TO_BASIS,:,:,:,:), sortlens, sortlims, &
        BASIS_TO_BASIS, ndim)

    deallocate(basis_inds, grid_inds)
 

   call initialize_mvp(Bz, Binv, potp, sortlens, sortlims, mappings)
 
   print *, 'starting ARPACK..'

   call solve_arpack(nb, mvp, dr, di, zeig, iparam, ido, info)

   call system_clock(time2)

   call write_output(ndim, basis_max_ind, grid_max_ind, nb, np, size(dr), &
   iparam(5), iparam(3), info, dr, di, &
   real(time2 - time1, real64) / real(rate, real64)) 

end program main 
  
!some legacy utils 
!  used in the sort() subroutine 
integer function findminimum(x, Start, sto)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                :: Start, sto
      real*8 :: x(start-sto+1)
      real*8                             :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)       ! assume the first is the min
      Location = Start          ! record its position
      DO i = Start+1, sto       ! start with next elements
         IF (x(i) < Minimum) THEN   !   if x(i) less than the min?
            Minimum  = x(i)     !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real*8, INTENT(INOUT) :: a, b
      real*8                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------


   SUBROUTINE  Sort(x, len)
      IMPLICIT  NONE
      real*8 :: x(len)
      INTEGER, INTENT(IN)                   :: len
      INTEGER                               :: i
      INTEGER                               :: Location
      INTEGER :: FindMinimum


      DO i = 1, len-1           ! except for the last
         Location = FindMinimum(x, i, len)  ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
   
   

   
