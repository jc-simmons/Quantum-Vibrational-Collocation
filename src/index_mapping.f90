module index_mapping_module

use index_generator_module, only : generate_indices

implicit none


!======================================================================
! Index mappings for sequential matrix-vector products
!
! Constructs index mappings used to evaluate multidimensional matrix-
! vector products as a sequence of one-dimensional operations.
!
! For each active dimension, the relevant multidimensional indices are
! sorted and grouped by their values in all non-active dimensions.
! The resulting groups identify the contiguous one-dimensional
! operations required for the active dimension.
!
! Mappings are constructed between basis and point/grid index spaces,
! as well as within a single index space. For mappings between different
! spaces, mappings are constructed both before and after the active
! dimension is changed.
!
! mapping contains the resulting index permutations.
! group_lims contains the number of indices in each group.
! mapping_lengths stores the index and group counts associated with
! each mapping, state, and active dimension.
!======================================================================

!----------------------------------------------------------------------
! Position relative to the change of the active dimension
!----------------------------------------------------------------------

integer, parameter :: BEFORE_CHANGE = 1
integer, parameter :: AFTER_CHANGE  = 2


!----------------------------------------------------------------------
! mapping_lengths quantities
!
! mapping_lengths(mapping_type, dimension)
!----------------------------------------------------------------------


contains

subroutine build_space_map(input_indices, input_max, output_max, &
same_space, mapping, &
mapping_lengths, group_lims, &
mapping_type, ndim)


! Construct index mappings for one input/output index-space pair.
!
! The input indices are sorted for each active dimension and partitioned
! into groups with identical values in all non-active dimensions.
! For mappings between different spaces, the active dimension is then
! changed in the index set, and the resulting indices are sorted and
! grouped again.
!
! The mappings and associated group information are stored in the
! arrays supplied by the caller. current_indices is maintained
! internally and is sized to the current set of indices.
!
! same_space indicates that the input and output index spaces are the
! same, in which case no new index set may need to be generated, though
! is not fully implemented 
!
! input_max and output_max define the bounds used when generating
! indices for the input and output spaces, respectively.

integer, intent(in) :: input_indices(:,:)
integer, intent(in) :: input_max
integer, intent(in) :: output_max
logical, intent(in) :: same_space
integer, intent(inout) :: mapping(:,:,:,:)
integer, intent(inout) :: mapping_lengths(:,:)
integer, intent(inout) :: group_lims(:,:,:,:)
integer, intent(in) :: mapping_type
integer, intent(in) :: ndim

integer, allocatable :: current_indices(:,:)
integer :: index_limits(ndim)
integer :: sort_permutation(maxval(shape(mapping)))

integer :: number_of_input_indices
integer :: dimension


number_of_input_indices = size(input_indices,2)

current_indices = input_indices

index_limits = input_max - 1


do dimension = 1, ndim

  current_indices([ndim,ndim-dimension+1],:) = &
      current_indices([ndim-dimension+1,ndim],:)

  call hpsort2(current_indices, &
      ndim, size(current_indices,2), &
      size(current_indices,2), sort_permutation)

  if (same_space) then

    call group_indices( &
        current_indices, &
        mapping(BEFORE_CHANGE,dimension,:,:), &
        group_lims(mapping_type,BEFORE_CHANGE,dimension,:), &
        ndim, &
        mapping_lengths(mapping_type, dimension), &
        sort_permutation)

    current_indices([ndim,ndim-dimension+1],:) = &
        current_indices([ndim-dimension+1,ndim],:)

  else

    call group_indices( &
        current_indices, &
        mapping(BEFORE_CHANGE,dimension,:,:), &
        group_lims(mapping_type,BEFORE_CHANGE,dimension,:), &
        ndim, &
        mapping_lengths(mapping_type, dimension), &
        sort_permutation)

    index_limits(ndim-dimension+1) = output_max - 1

    call generate_indices(index_limits, current_indices)

    mapping_lengths(mapping_type,dimension) = &
        size(current_indices,2)


    current_indices([ndim,ndim-dimension+1],:) = &
        current_indices([ndim-dimension+1,ndim],:)

    call hpsort2( &
        current_indices, ndim, size(current_indices,2), &
        size(current_indices,2), sort_permutation)

    call group_indices( &
        current_indices, &
        mapping(AFTER_CHANGE,dimension,:,:), &
        group_lims(mapping_type,AFTER_CHANGE,dimension,:), &
        ndim, &
        mapping_lengths(mapping_type, dimension), &
        sort_permutation)

  end if

end do


end subroutine build_space_map

!======================================================================
! group_indices
!
! Groups a set of already-sorted multidimensional indices according to
! the values of their non-active dimensions.
!
! Consecutive indices whose first ndim-1 components are identical are
! assigned to the same group. The sorting permutation is retained so
! that each group records the corresponding original index positions.
!
! group_lims records the number of indices in each group.
! mapping_lengths records the total index count and number of groups.
!======================================================================

subroutine group_indices(presort, postsort, group_lims, &
        ndim, mapping_lengths, sort_permutation)

integer :: ndim
integer :: i
integer :: group

integer :: presort(:,:)
integer :: postsort(:,:)
integer :: group_lims(:)
integer :: mapping_lengths
integer :: sort_permutation(:)

integer :: previous_index(ndim-1)
integer :: index_count

group_lims = 0
postsort = 0

previous_index = presort(1:ndim-1,1)
group = 1

index_count = size(presort, 2)

do i=1, index_count

if (all(previous_index == presort(1:ndim-1,i))) then

  group_lims(group) = group_lims(group) + 1
  postsort(group,group_lims(group)) = sort_permutation(i)

else

  group = group + 1

  group_lims(group) = group_lims(group) + 1
  postsort(group,group_lims(group)) = sort_permutation(i)

  previous_index = presort(1:ndim-1,i)

end if


end do

mapping_lengths = group

end subroutine group_indices

end module index_mapping_module
