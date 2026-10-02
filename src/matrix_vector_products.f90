module mvp_module

    use iso_fortran_env, only: real64
    use mapping_parameters

    implicit none

    private

    public :: initialize_mvp, mvp, mapped_multiply

    real(real64), allocatable :: Bz(:,:,:,:,:,:)
    real(real64), allocatable :: Binv(:,:,:,:)
    real(real64), allocatable :: potp(:)

    integer, allocatable :: sortlens(:,:)
    integer, allocatable :: sortlims(:,:,:,:)
    integer, allocatable :: mappings(:,:,:,:,:)

    integer :: np
    integer :: nx
    integer :: d
    integer :: terms

    logical :: initialized = .false.

contains

    subroutine initialize_mvp(Bz_in, Binv_in, potp_in, sortlens_in, &
                              sortlims_in, mappings_in)

        real(real64), intent(in) :: Bz_in(:,:,:,:,:,:)
        real(real64), intent(in) :: Binv_in(:,:,:,:)
        real(real64), intent(in) :: potp_in(:)

        integer, intent(in) :: sortlens_in(:,:)
        integer, intent(in) :: sortlims_in(:,:,:,:)
        integer, intent(in) :: mappings_in(:,:,:,:,:)

        Bz       = Bz_in
        Binv     = Binv_in
        potp     = potp_in
        sortlens = sortlens_in
        sortlims = sortlims_in
        mappings = mappings_in

        np    = size(potp)
        nx    = size(Bz, 5)
        d     = size(Bz, 3)
        terms = size(Bz, 4)

        initialized = .true.

    end subroutine initialize_mvp


function mvp(x) result(y)

    real(real64), intent(in) :: x(:)
    real(real64) :: y(size(x))

    real(real64) :: v(np)
    real(real64) :: res(np + 1)
    real(real64) :: uu(np)

    integer :: i
    integer :: l
    integer :: nb

    if (.not. initialized) then
        error stop 'mvp_module: initialize_mvp must be called first'
    end if

    nb = size(x)

    if (nb > np) then
        error stop 'mvp_module: input vector is larger than the potential grid'
    end if

    uu  = 0.0_real64
    v   = 0.0_real64
    res = 0.0_real64

    res(2:nb + 1) = x

        ! BU
        call mapped_multiply(res, v, &
            mappings(BASIS_TO_BASIS,INPUT_STATE,:,:,:), &
            mappings(BASIS_TO_BASIS,INPUT_STATE,:,:,:), &
            sortlims(BASIS_TO_BASIS,INPUT_STATE,:,:), &
            sortlims(BASIS_TO_BASIS,OUTPUT_STATE,:,:), &
            sortlens(BASIS_TO_BASIS,:), &
            Bz(B_MATRIX,U_FACTOR,:,1,:,:), d)

        ! BL
        call mapped_multiply(res, v, &
            mappings(BASIS_TO_POINTS,INPUT_STATE,:,:,:), &
            mappings(BASIS_TO_POINTS,OUTPUT_STATE,:,:,:), &
            sortlims(BASIS_TO_POINTS,INPUT_STATE,:,:), &
            sortlims(BASIS_TO_POINTS,OUTPUT_STATE,:,:), &
            sortlens(BASIS_TO_POINTS,:), &
            Bz(B_MATRIX,L_FACTOR,:,1,:,:), d)

        ! Potential
        do i = 1, np
            uu(i) = uu(i) + res(i + 1) * potp(i)
        end do

        ! B''
        do l = 1, terms

                res = 0.0_real64
                res(2:nb + 1) = x

            call mapped_multiply(res, v, &
                mappings(BASIS_TO_BASIS,INPUT_STATE,:,:,:), &
                mappings(BASIS_TO_BASIS,INPUT_STATE,:,:,:), &
                sortlims(BASIS_TO_BASIS,INPUT_STATE,:,:), &
                sortlims(BASIS_TO_BASIS,OUTPUT_STATE,:,:), &
                sortlens(BASIS_TO_BASIS,:), &
                Bz(B2_MATRIX,U_FACTOR,:,terms-l+1,:,:), d)

            call mapped_multiply(res, v, &
                mappings(BASIS_TO_POINTS,INPUT_STATE,:,:,:), &
                mappings(BASIS_TO_POINTS,OUTPUT_STATE,:,:,:), &
                sortlims(BASIS_TO_POINTS,INPUT_STATE,:,:), &
                sortlims(BASIS_TO_POINTS,OUTPUT_STATE,:,:), &
                sortlens(BASIS_TO_POINTS,:), &
                Bz(B2_MATRIX,L_FACTOR,:,terms-l+1,:,:), d)

            do i = 1, np
                uu(i) = uu(i) + res(i + 1)
            end do

        end do

        ! inv(BL)
            res = 0.0_real64
            res(2:) = uu

        call mapped_multiply(res, v, &
            mappings(POINTS_TO_POINTS,INPUT_STATE,:,:,:), &
            mappings(POINTS_TO_POINTS,INPUT_STATE,:,:,:), &
            sortlims(POINTS_TO_POINTS,INPUT_STATE,:,:), &
            sortlims(POINTS_TO_POINTS,OUTPUT_STATE,:,:), &
            sortlens(POINTS_TO_POINTS,:), &
            Binv(:,L_FACTOR,:,:), d)

        ! inv(BU)
        call mapped_multiply(res, v, &
            mappings(POINTS_TO_BASIS,INPUT_STATE,:,:,:), &
            mappings(POINTS_TO_BASIS,OUTPUT_STATE,:,:,:), &
            sortlims(POINTS_TO_BASIS,INPUT_STATE,:,:), &
            sortlims(POINTS_TO_BASIS,OUTPUT_STATE,:,:), &
            sortlens(POINTS_TO_BASIS,:), &
            Binv(:,U_FACTOR,:,:), d)

    y = v(:nb)

end function mvp


subroutine mapped_multiply(input, output, mapin, mapout, limin, limout, &
                           lens, matrix, d)

    implicit none

    integer, intent(in) :: d

    integer, intent(in) :: lens(:)
    integer, intent(in) :: mapin(:,:,:)
    integer, intent(in) :: mapout(:,:,:)
    integer, intent(in) :: limin(:,:)
    integer, intent(in) :: limout(:,:)

    real(real64), intent(inout) :: input(:)
    real(real64), intent(out) :: output(:)
    real(real64), intent(in) :: matrix(:,:,:)
    real(real64) :: sum

    integer :: i
    integer :: j
    integer :: k
    integer :: l

    sum = 0.0_real64

    do i = 1, d

        output = 0.0_real64

        do j = 1, lens(i)
            do k = 1, limout(i,j)
                sum = 0.0_real64
                do l = 1, limin(i,j)
                    sum = sum + matrix(d-i+1,k,l) * input(mapin(i,j,l) + 1)
                end do
                output(mapout(i,j,k)) = sum
            end do
        end do

        input(2:) = output

    end do

end subroutine mapped_multiply

end module mvp_module