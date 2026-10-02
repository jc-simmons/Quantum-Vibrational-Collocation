module mapping_parameters

    implicit none

    integer, parameter :: NUM_MAP_TYPES = 4

    integer, parameter :: BASIS_TO_POINTS  = 1
    integer, parameter :: POINTS_TO_BASIS  = 2
    integer, parameter :: POINTS_TO_POINTS = 3
    integer, parameter :: BASIS_TO_BASIS  = 4

    integer, parameter :: NUM_STATES = 2

    integer, parameter :: INPUT_STATE  = 1
    integer, parameter :: OUTPUT_STATE = 2

    integer, parameter :: NUM_B_MATRICES = 2

    integer, parameter :: B_MATRIX  = 1
    integer, parameter :: B2_MATRIX = 2

    integer, parameter :: NUM_LU_FACTORS = 2

    integer, parameter :: L_FACTOR = 1
    integer, parameter :: U_FACTOR = 2

end module mapping_parameters