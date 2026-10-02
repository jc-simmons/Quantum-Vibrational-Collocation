module basis_keo_module

    use iso_fortran_env, only: real64

    implicit none

    public :: hob, hob_keo, beta, omega

    real(real64), parameter :: pi = 3.14159265358979323846_real64
    real(real64) :: beta = 1.0_real64
    real(real64) :: omega(100) = 1.0_real64

    namelist /basis_config/ omega

contains

    !===============================================================
    ! One-dimensional harmonic oscillator basis function
    !
    ! k = basis-function index, starting from 0
    ! x = coordinate value
    ! dim = coordinate dimension
    !
    ! Returns the one-dimensional factor phi_k(x) associated with
    ! dimension dim in the multidimensional product basis function:
    !
    !
    ! dim is included so that the basis-function interface is
    ! dimension-aware, although it does not currently affect the
    ! harmonic oscillator basis function.
    !===============================================================

    real(real64) function hob(k, x, dim) result(value)

        integer, intent(in) :: k, dim
        real(real64), intent(in) :: x

        real(real64) :: scaled_x

        scaled_x = sqrt(beta) * x

        value = beta**0.25_real64 &
              / sqrt(2.0_real64**k * fact(k)) &
              / pi**0.25_real64 &
              * exp(-beta * x**2 / 2.0_real64) &
              * hermite(k, scaled_x)

    end function hob


    !===============================================================
    ! KEO applied to harmonic oscillator basis function
    !
    ! k = basis-function index, starting from 0
    ! x = coordinate value
    ! dim = coordinate dimension
    ! nterm = KEO term
    !
    ! Returns the one-dimensional factor associated with dimension
    ! dim after application of KEO term nterm to the multidimensional
    ! harmonic oscillator basis function.
    !===============================================================

real(real64) function hob_keo(k, x, dim, nterm) result(value)

        integer, intent(in) :: k, dim, nterm
        real(real64), intent(in) :: x
        
        real(real64) :: prefactor

        if (dim == nterm) then
            
            prefactor = omega(dim) * beta / 4.0_real64

            value = prefactor * ( &
                (2.0_real64 * real(k, real64) + 1.0_real64) * hob(k, x, dim) &
                - sqrt(real((k + 1) * (k + 2), real64)) * hob(k + 2, x, dim) &
            )

            if (k >= 2) then
                value = value - prefactor &
                    * sqrt(real(k * (k - 1), real64)) * hob(k - 2, x, dim)
            end if

        else
            value = hob(k, x, dim)
        end if

    end function hob_keo


    !===============================================================
    ! Hermite polynomial H_k(x)
    !===============================================================
    recursive real(real64) function hermite(k, x) result(value)

        integer, intent(in) :: k
        real(real64), intent(in) :: x

        if (k == 0) then

            value = 1.0_real64

        else if (k == 1) then

            value = 2.0_real64 * x

        else

            value = 2.0_real64 * x * hermite(k - 1, x) &
                  - 2.0_real64 * real(k - 1, real64) &
                  * hermite(k - 2, x)

        end if

    end function hermite


    !===============================================================
    ! Factorial
    !===============================================================
    real(real64) function fact(k) result(value)

        integer, intent(in) :: k
        integer :: i

        value = 1.0_real64

        do i = 1, k
            value = value * real(i, real64)
        end do

    end function fact

end module basis_keo_module