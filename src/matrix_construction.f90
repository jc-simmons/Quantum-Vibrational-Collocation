module matrix_module

    use iso_fortran_env, only: real64

    implicit none

    private

    public :: construct_matrices
    public :: construct_basis_matrix
    public :: construct_keo_matrix

!===============================================================
    ! Interface for a basis function
    !
    ! basis(k,x,dim) returns the one-dimensional factor of the
    ! basis function with index k at coordinate x in dimension dim.
    !===============================================================

    abstract interface

        function basis_function(k, x, dim) result(value)

            import :: real64

            integer, intent(in) :: k, dim
            real(real64), intent(in) :: x
            real(real64) :: value

        end function basis_function

    end interface


    !===============================================================
    ! Interface for a KEO function
    !
    ! keo(k,x,dim,nterm) returns the one-dimensional factor
    ! associated with dimension dim after application of KEO term
    ! nterm to the multidimensional basis function.
    !===============================================================

    abstract interface

        function keo_function(k, x, dim, nterm) result(value)

            import :: real64

            integer, intent(in) :: k, dim, nterm
            real(real64), intent(in) :: x
            real(real64) :: value

        end function keo_function

    end interface


contains


    !===============================================================
    ! Construct all matrix-derived quantities required by the
    ! calculation.
    !
    ! B and KEO-applied matrices are constructed internally, then
    ! factorized into L and U and stored in Bz.
    !
    ! Bz(1,1,dim,1,:,:) = L for B
    ! Bz(1,2,dim,1,:,:) = U for B
    !
    ! Bz(2,1,dim,nterm,:,:) = L for KEO term nterm
    ! Bz(2,2,dim,nterm,:,:) = U for KEO term nterm
    !
    ! Binv(dim,1,:,:) = inverse of L for B
    ! Binv(dim,2,:,:) = inverse of U for B
    !
    ! The matrix size for each dimension is n(dim).
    !===============================================================

    subroutine construct_matrices(grid, n, kterms, basis, keo, Bz, Binv)

        real(real64), intent(in) :: grid(:,:)
        integer, intent(in) :: n(:)
        integer, intent(in) :: kterms

        procedure(basis_function) :: basis
        procedure(keo_function) :: keo

        real(real64), intent(out) :: Bz(:,:,:,:,:,:)
        real(real64), intent(out) :: Binv(:,:,:,:)

        real(real64), allocatable :: B(:,:)
        real(real64), allocatable :: Bkeo(:,:)
        real(real64), allocatable :: L(:,:)
        real(real64), allocatable :: U(:,:)

        integer :: ndim
        integer :: nx
        integer :: dim
        integer :: nterm
        integer :: ni

        ndim = size(n)
        nx = size(Bz,5)

        if (size(grid,2) /= ndim) then
            error stop "construct_matrices: grid has incorrect number of dimensions."
        end if

        if (size(Bz,1) /= 2 .or. size(Bz,2) /= 2) then
            error stop "construct_matrices: incorrect first two dimensions of Bz."
        end if

        if (size(Bz,3) /= ndim) then
            error stop "construct_matrices: incorrect dimension dimension of Bz."
        end if

        if (size(Bz,4) < kterms) then
            error stop "construct_matrices: Bz does not have enough KEO terms."
        end if

        if (size(Binv,1) /= ndim .or. size(Binv,2) /= 2) then
            error stop "construct_matrices: incorrect dimensions of Binv."
        end if

        Bz = 0.0_real64
        Binv = 0.0_real64

        do dim = 1, ndim

            ni = n(dim)

            if (ni > nx) then
                error stop "construct_matrices: n(dim) exceeds matrix storage size."
            end if

            allocate(B(ni,ni))
            allocate(L(ni,ni))
            allocate(U(ni,ni))

            !-------------------------------------------------------
            ! Basis matrix for this dimension
            !-------------------------------------------------------

            call construct_basis_matrix( &
                grid(1:ni,dim), basis, dim, B)

            call lu_nopivot(B, L, U)

            Bz(1,1,dim,1,1:ni,1:ni) = L
            Bz(1,2,dim,1,1:ni,1:ni) = U

            !-------------------------------------------------------
            ! Inverses of L and U for the ordinary basis matrix
            !-------------------------------------------------------

            Binv(dim,1,1:ni,1:ni) = inverse(L)
            Binv(dim,2,1:ni,1:ni) = inverse(U)

            !-------------------------------------------------------
            ! KEO terms for this dimension
            !-------------------------------------------------------

            do nterm = 1, kterms

                allocate(Bkeo(ni,ni))

                call construct_keo_matrix( &
                    grid(1:ni,dim), keo, dim, nterm, Bkeo)

                call lu_nopivot(Bkeo, L, U)

                Bz(2,1,dim,nterm,1:ni,1:ni) = L
                Bz(2,2,dim,nterm,1:ni,1:ni) = U

                deallocate(Bkeo)

            end do

            deallocate(B)
            deallocate(L)
            deallocate(U)

        end do

    end subroutine construct_matrices



    !===============================================================
    ! Construct 1D basis collocation matrix
    !
    ! B(j,k) = basis(k-1, points(j), dim)
    !
    ! The basis index passed to the basis routine is zero-based,
    ! while the matrix indices are one-based.
    !===============================================================

    subroutine construct_basis_matrix(points, basis, dim, B)

        real(real64), intent(in) :: points(:)
        procedure(basis_function) :: basis
        integer, intent(in) :: dim
        real(real64), intent(out) :: B(:,:)

        integer :: j,k

        if (size(B,1) /= size(points)) then
            error stop "construct_basis_matrix: incorrect number of rows in B."
        end if

        
        do j = 1, size(points)

            do k = 1, size(B,2)

                B(j,k) = basis(k-1, points(j), dim)

            end do

        end do

    end subroutine construct_basis_matrix


    !===============================================================
    ! Construct 1D KEO-applied collocation matrix
    !
    ! Bkeo(j,k) = keo(k-1, points(j), dim, nterm)
    !
    ! The returned matrix contains the one-dimensional factor
    ! associated with dimension dim after application of KEO term
    ! nterm to the multidimensional basis function.
    !===============================================================

    subroutine construct_keo_matrix(points, keo, dim, nterm, Bkeo)

        real(real64), intent(in) :: points(:)
        procedure(keo_function) :: keo
        integer, intent(in) :: dim
        integer, intent(in) :: nterm
        real(real64), intent(out) :: Bkeo(:,:)

        integer :: j
        integer :: k

        if (size(Bkeo,1) /= size(points)) then
            error stop "construct_keo_matrix: incorrect number of rows in Bkeo."
        end if

        do j = 1, size(points)

            do k = 1, size(Bkeo,2)

                Bkeo(j,k) = keo(k-1, points(j), dim, nterm)

            end do

        end do

    end subroutine construct_keo_matrix


    !===============================================================
    ! LU factorization without pivoting
    !
    ! A = L * U
    !
    ! This is the existing no-pivot algorithm.
    !===============================================================
    subroutine lu_nopivot(A, L, U)

        real(real64), intent(in) :: A(:,:)
        real(real64), intent(out) :: L(:,:)
        real(real64), intent(out) :: U(:,:)

        real(real64), allocatable :: work(:,:)

        integer :: n
        integer :: i
        integer :: j
        integer :: k

        n = size(A,1)

        if (size(A,2) /= n) then
            error stop "lu_nopivot: A must be square."
        end if

        if (size(L,1) /= n .or. size(L,2) /= n) then
            error stop "lu_nopivot: incorrect size for L."
        end if

        if (size(U,1) /= n .or. size(U,2) /= n) then
            error stop "lu_nopivot: incorrect size for U."
        end if

        allocate(work(n,n))

        work = A
        L = 0.0_real64

        do i = 1, n
            L(i,i) = 1.0_real64
        end do

        do k = 1, n

            if (abs(work(k,k)) < 1.0e-8_real64) then
                error stop "lu_nopivot: zero or near-zero pivot."
            end if

            if (k < n) then

                L(k+1:n,k) = work(k+1:n,k) / work(k,k)

                do j = k + 1, n
                    work(j,:) = work(j,:) &
                              - L(j,k) * work(k,:)
                end do

            end if

        end do

        U = work

        deallocate(work)

    end subroutine lu_nopivot


    !===============================================================
    ! Matrix inverse using LAPACK
    ! modified from http://fortranwiki.org/fortran/show/Matrix+inversion
    !===============================================================
    function inverse(A) result(Ainv)

        real(real64), intent(in) :: A(:,:)
        real(real64) :: Ainv(size(A,1),size(A,2))

        integer :: n
        integer :: info
        integer, allocatable :: ipiv(:)
        real(real64), allocatable :: work(:)

        external :: DGETRF
        external :: DGETRI

        n = size(A,1)

        if (size(A,2) /= n) then
            error stop "inverse: A must be square."
        end if

        allocate(ipiv(n))
        allocate(work(n))

        Ainv = A

        call dgetrf(n, n, Ainv, n, ipiv, info)

        if (info < 0) then
            error stop "inverse: invalid argument in DGETRF."
        else if (info > 0) then
            error stop "inverse: matrix is singular."
        end if

        call dgetri(n, Ainv, n, ipiv, work, n, info)

        if (info < 0) then
            error stop "inverse: invalid argument in DGETRI."
        else if (info > 0) then
            error stop "inverse: matrix inversion failed."
        end if

        deallocate(ipiv)
        deallocate(work)

    end function inverse

end module matrix_module