module potential_module

    use iso_fortran_env, only: real64

    implicit none

    character(len=256) :: potential_routine
    character(len=256) :: potential_source
    integer :: potential_terms

    namelist /potential_config/ potential_routine, &
    potential_source, potential_terms

contains

    subroutine compute_potential(grid, grid_indices, potential_vals)

        integer, intent(in) :: grid_indices(:,:)
        real(8), intent(in)  :: grid(:,:)
        real(8), allocatable, intent(out) :: potential_vals(:)
        real(8), allocatable :: point(:)
        integer :: i, j
        integer :: grid_length, grid_dims

        grid_dims = size(grid_indices, 1)
        grid_length = size(grid_indices, 2)
    
        allocate(point(grid_dims))
        allocate(potential_vals(grid_length))

        select case (trim(potential_routine))

        case ('SOP')

            do i = 1, grid_length
                do j = 1, grid_dims 
                    point(j) = grid(grid_indices(j,i),j)
                end do 
                call sop_potential(point, potential_vals(i), potential_source, potential_terms)
            end do

        case default !can extend with more cases

            error stop 'Unknown potential routine'

        end select

    end subroutine compute_potential


    subroutine sop_potential(x, val, potential_file, npotential)

        real(real64), intent(in) :: x(:)
        real(real64), intent(out) :: val
        character(len=*), intent(in) :: potential_file
        integer, intent(in) :: npotential

        integer, allocatable, save :: fexp(:,:)
        real(real64), allocatable, save :: fcoef(:)
        logical, save :: initialized = .false.

        integer :: i, j, d
        integer :: unit
        real(real64) :: prod


        if (.not. initialized) then

            d = size(x)

            open(newunit=unit, file=trim(potential_file), &
                 status='old', action='read')

            allocate(fexp(npotential,d))
            allocate(fcoef(npotential))

            do i = 1, npotential
                read(unit, *) fexp(i,:), fcoef(i)
            end do

            close(unit)

            initialized = .true.

        end if

        val = 0.0_real64

        do i = 1, size(fcoef)
            prod = 1.0_real64

            do j = 1, size(x)
                prod = prod * x(j)**fexp(i,j)
            end do

            val = val + prod * fcoef(i)
        end do

    end subroutine sop_potential

end module potential_module