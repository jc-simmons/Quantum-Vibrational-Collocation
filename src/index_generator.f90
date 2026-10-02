module index_generator_module

  implicit none

contains


subroutine generate_indices(limits, indices)

  integer, intent(in) :: limits(:)
  integer, allocatable, intent(inout) :: indices(:,:)

  integer :: ndim
  integer :: nindices
  integer :: dim
  integer :: max_index
  integer :: prefix_sum
  integer, allocatable :: current(:)

  ndim = size(limits)
  nindices = count_indices(limits)

  if (allocated(indices)) deallocate(indices)
  allocate(indices(ndim,nindices))

  allocate(current(ndim))

  current = 0
  dim = 1
  nindices = 0

  do

    prefix_sum = sum(current(1:dim-1))
    max_index = max(limits(dim)-prefix_sum,0)

    if (current(dim) <= max_index) then

      if (dim == ndim) then

        nindices = nindices + 1
        indices(:,nindices) = current + 1

        current(dim) = current(dim) + 1

      else

        dim = dim + 1
        current(dim) = 0

      end if

    else

      if (dim == 1) exit

      current(dim) = 0
      dim = dim - 1
      current(dim) = current(dim) + 1

    end if

  end do

  deallocate(current)

end subroutine generate_indices


  integer function count_indices(limits)

    integer, intent(in) :: limits(:)

    integer :: ndim
    integer :: max_sum
    integer :: dim
    integer :: old_sum
    integer :: new_value
    integer :: new_sum

    integer, allocatable :: counts(:)
    integer, allocatable :: new_counts(:)

    ndim = size(limits)
    max_sum = maxval(limits)

    allocate(counts(0:max_sum))
    allocate(new_counts(0:max_sum))

    counts = 0
    counts(0) = 1

    do dim = 1, ndim

      new_counts = 0

      do old_sum = 0, max_sum

        if (counts(old_sum) == 0) cycle

        do new_value = 0, max(limits(dim)-old_sum,0)

          new_sum = old_sum + new_value

          new_counts(new_sum) = new_counts(new_sum) + &
                                counts(old_sum)

        end do

      end do

      counts = new_counts

    end do

    count_indices = sum(counts)

    deallocate(counts)
    deallocate(new_counts)

  end function count_indices

end module index_generator_module

