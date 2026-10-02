module grid_module

    use iso_fortran_env, only: real64

    implicit none

    integer :: grid_id(100) = 0
    character(len=32) :: grid_routine(100)
    character(len=256) :: grid_source(100)

    namelist /grid_config/ grid_id, grid_routine, grid_source

contains

    subroutine compute_grid(grid, ndim, pmax)

        integer, intent(in) :: ndim
        integer, intent(in) :: pmax
        real(real64), allocatable, intent(out) :: grid(:,:)

        integer :: i, j, unit, ios
        integer :: ngrid
        real(real64) :: work_grid(pmax)

        allocate(grid(pmax, ndim))

        grid = 0.0_real64
        ngrid = maxval(grid_id)

        do i = 1, ngrid

            select case (trim(grid_routine(i)))
            case ('LOAD')

                open(newunit=unit, file=trim(grid_source(i)), &
                     status='old', action='read', iostat=ios)

                if (ios /= 0) then
                    error stop 'Unable to open grid file'
                end if

                read(unit, *, iostat=ios) work_grid

                if (ios /= 0) then
                    close(unit)
                    error stop 'Unable to read grid file'
                end if

                close(unit)

            case ('GEN')

                ! Generate grid using source(i)

            case default

                error stop 'Unknown grid routine'

            end select

            do j = 1, ndim
                if (grid_id(j) == i) then
                    grid(:,j) = work_grid
                end if
            end do

        end do

    end subroutine compute_grid

end module grid_module