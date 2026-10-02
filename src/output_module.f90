! DATE_AND_TIME
! https://gcc.gnu.org/onlinedocs/gfortran/DATE_005fAND_005f_TIME.html
!
! IOSTAT
! https://support.nag.com/nagware/np/r70_doc/manual/compiler_8_3.html

module output_module

    use iso_fortran_env, only : real64, output_unit

    implicit none

    logical :: print_output = .true.
    logical :: save_output  = .false.

    !default but can be saved to another directory 
    character(len=256) :: output_dir  = 'results/'
    character(len=256) :: run_name = ''
    character(len=256) :: output_file = ''

    namelist /output_config/ print_output, save_output, output_dir, run_name

contains


    subroutine write_output(ndim, basis_max, point_max, total_basis, &
                            total_points, &
                            n_requested, n_converged, iterations, info, &
                            dr, di, run_time)

        integer, intent(in) :: ndim
        integer, intent(in) :: basis_max
        integer, intent(in) :: point_max
        integer, intent(in) :: total_basis
        integer, intent(in) :: total_points
        integer, intent(in) :: n_requested
        integer, intent(in) :: n_converged
        integer, intent(in) :: iterations
        integer, intent(in) :: info

        real(real64), intent(in) :: dr(:)
        real(real64), intent(in) :: di(:)
        real(real64), intent(in) :: run_time

        integer :: n_values
        integer :: unit
        integer :: ios
        character(len=8)  :: date
        character(len=10) :: time

        n_values = min(n_converged, size(dr), size(di))

        call date_and_time(date=date, time=time)

        if (len_trim(run_name) == 0) then
            run_name = 'run_' // trim(date) // '_' // time(1:6)
        end if

        output_file = trim(output_dir) // trim(run_name) // '.dat'

        if (save_output) then

            open(newunit=unit, file=trim(output_file), status='replace', &
                 action='write', iostat=ios)

            if (ios /= 0) then
                write(*, '(A)') &
                    'Output error: could not open output file: ' // trim(output_file)
                stop
            end if

            call log_report(unit, .true., ndim, basis_max, point_max, &
                            total_basis, total_points,  &
                            n_requested, n_converged, iterations, info, &
                            dr, di, run_time)

            close(unit)

        end if

        if (print_output) then

            call log_report(output_unit, .false., ndim, basis_max, point_max, &
                            total_basis, total_points, &
                            n_requested, n_converged, iterations, info, &
                            dr, di, run_time)

        end if

    end subroutine write_output


    subroutine log_report(unit, file_format, ndim, basis_max, point_max, &
                          total_basis, total_points,  &
                          n_requested, n_converged, iterations, info, &
                          dr, di, run_time)

        integer, intent(in) :: unit
        logical, intent(in) :: file_format
        integer, intent(in) :: ndim
        integer, intent(in) :: basis_max
        integer, intent(in) :: point_max
        integer, intent(in) :: total_basis
        integer, intent(in) :: total_points
        integer, intent(in) :: n_requested
        integer, intent(in) :: n_converged
        integer, intent(in) :: iterations
        integer, intent(in) :: info

        real(real64), intent(in) :: dr(:)
        real(real64), intent(in) :: di(:)
        real(real64), intent(in) :: run_time

        integer :: i, n_values

        n_values = min(n_converged, size(dr), size(di))

        if (file_format) then

            write(unit, '(A)') '# Rectangular Collocation Eigenvalue Solver'
            write(unit, '(A)') '# Run: ' // trim(run_name)
            write(unit, '(A)') '#'
            write(unit, '(A)') '# System'
            write(unit, '(A,I10)') '#   Dimensions:             ', ndim
            write(unit, '(A,I10)') '#   Basis maximum:          ', basis_max
            write(unit, '(A,I10)') '#   Point maximum:          ', point_max
            write(unit, '(A,I10)') '#   Basis functions:        ', total_basis
            write(unit, '(A,I10)') '#   Grid points:            ', total_points
            write(unit, '(A)') '#'
            write(unit, '(A)') '# ARPACK'
            write(unit, '(A,I10)') '#   Requested eigenvalues:  ', n_requested
            write(unit, '(A,I10)') '#   Converged eigenvalues:  ', n_converged
            write(unit, '(A,I10)') '#   Iterations:             ', iterations
            write(unit, '(A,I10)') '#   INFO:                   ', info
            write(unit, '(A,F12.3,A)') '#   Run time:               ', run_time, ' s'
            write(unit, '(A)') '#'
            write(unit, '(A)') '# Eigenvalues'
            write(unit, '(A)') '#   Real                    Imaginary'

        else

            write(unit, '(A)') 'Rectangular Collocation Eigenvalue Solver'
            write(unit, '(A)') 'Run: ' // trim(run_name)
            write(unit, '(A)') ''
            write(unit, '(A)') 'System'
            write(unit, '(A,I10)') '  Dimensions:             ', ndim
            write(unit, '(A,I10)') '  Basis maximum:          ', basis_max
            write(unit, '(A,I10)') '  Point maximum:          ', point_max
            write(unit, '(A,I10)') '  Basis functions:        ', total_basis
            write(unit, '(A,I10)') '  Grid points:            ', total_points
            write(unit, '(A)') ''
            write(unit, '(A)') 'ARPACK'
            write(unit, '(A,I10)') '  Requested eigenvalues:  ', n_requested
            write(unit, '(A,I10)') '  Converged eigenvalues:  ', n_converged
            write(unit, '(A,I10)') '  Iterations:             ', iterations
            write(unit, '(A,I10)') '  INFO:                   ', info
            write(unit, '(A,F12.3,A)') '  Run time:               ', run_time, ' s'
            write(unit, '(A)') ''
            write(unit, '(A)') 'Eigenvalues'
            write(unit, '(A)') '    Real                    Imaginary'

        end if

        do i = 1, n_values
            write(unit, *) dr(i), di(i)
        end do

    end subroutine log_report

end module output_module
