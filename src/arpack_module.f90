module arpack_module

    use iso_fortran_env, only: real64
    use mvp_module, only: mvp

    implicit none
    private

    public :: solve_arpack

contains

    subroutine solve_arpack(nx, op, dr, di, zeig, iparam, ido, info)

        integer, intent(in) :: nx
        procedure(mvp) :: op

        real(real64), allocatable, intent(out) :: dr(:)
        real(real64), allocatable, intent(out) :: di(:)
        real(real64), allocatable, intent(out) :: zeig(:,:)
        integer, intent(out) :: iparam(11)
        integer, intent(out) :: ido
        integer, intent(out) :: info

        integer :: nev
        integer :: ncv
        integer :: maxiter
        integer :: lworkl

        integer :: ipntr(14)

        integer, allocatable :: iselect(:)

        real(real64) :: tol
        real(real64) :: sigmar
        real(real64) :: sigmai

        real(real64), allocatable :: resid(:)
        real(real64), allocatable :: workd(:)
        real(real64), allocatable :: workl(:)
        real(real64), allocatable :: v(:,:)
        real(real64), allocatable :: workev(:)

        character(len=1) :: bmat
        character(len=2) :: which
        character(len=1) :: howmny

        logical :: rvec

        maxiter = 150
        nev = min(nx - 2, 50)
        ncv = min(nx, 120)
        tol = 0.0_real64

        iparam = 0
        iparam(1) = 1
        iparam(3) = 2 * maxiter
        iparam(4) = 1
        iparam(7) = 1

        ido = 0
        bmat = 'I'
        which = 'SR'
        info = 0

        allocate(resid(nx))
        allocate(workd(3 * nx))
        allocate(v(nx, ncv))

        lworkl = 3 * ncv**2 + 8 * ncv
        allocate(workl(lworkl))

        v = 0.0_real64

        ! Reverse-communication iteration
        do while (ido < 99)

            call dnaupd(ido, bmat, nx, which, nev, tol, resid, ncv, v, nx, &
                        iparam, ipntr, workd, workl, lworkl, info)

            if (ido == -1 .or. ido == 1) then

                workd(ipntr(2):ipntr(2) + nx - 1) = &
                    op(workd(ipntr(1):ipntr(1) + nx - 1))

            end if

        end do

        ! Extract eigenvalues and eigenvectors
        rvec = .true.
        howmny = 'A'

        allocate(dr(nev + 1))
        allocate(di(nev + 1))
        allocate(zeig(nx, nev + 1))
        allocate(workev(3 * ncv))
        allocate(iselect(ncv))

        iselect = .true.
        zeig = v(:, 1:nev + 1)
        dr = 0.0_real64
        di = 0.0_real64

        call dneupd2(rvec, howmny, iselect, dr, di, zeig, nx, sigmar, sigmai, &
                     workev, bmat, nx, which, nev, tol, resid, ncv, v, nx, &
                     iparam, ipntr, workd, workl, lworkl, info)

    end subroutine solve_arpack

end module arpack_module