module artemis__misc
  !! Module containing ARTEMIS-specific miscellaneous utilities.
  !!
  !! Common utilities (sort, swap, set, string ops, file ops) are provided
  !! by the coreutils library and re-exported here for backward compatibility.
  !!
  !! ARTEMIS-specific procedures:
  !! - sort2D   : sort a dim×3 array by cycling through columns
  !! - sort_col : sort a 2D array by a specified column
  !! - loadbar  : write a loading bar to the terminal
  use coreutils__kind, only: real32
  use coreutils__array, only: swap, sort1D, set
  use coreutils__string, only: icount, flagmaker, to_upper, to_lower, strip_null
  use coreutils__file, only: grep, jump, file_check
  implicit none

  !> @deprecated Use coreutils__array directly for swap, sort1D, set
  public :: swap, sort1D, set
  !> @deprecated Use coreutils__string directly for icount, flagmaker, to_upper, to_lower, strip_null
  public :: icount, flagmaker, to_upper, to_lower, strip_null
  !> @deprecated Use coreutils__file directly for grep, jump, file_check
  public :: grep, jump, file_check
  public :: sort2D, sort_col, loadbar


contains

!###############################################################################
  subroutine sort2D(arr,dim)
    !! Sort a dim×3 array by cycling through columns.
    !!
    !! ARTEMIS-specific: sorts by minimum absolute value across cycling
    !! column order. Different from coreutils sort2D.
    implicit none

    ! Arguments
    integer :: dim
    !! Number of rows in the array.
    real(real32), dimension(dim,3) :: arr
    !! Array to sort in-place.

    ! Local variables
    integer :: i, j, loc, istart
    !! Loop indices and location buffer.
    integer, dimension(3) :: a123
    !! Column cycling index.
    real(real32), dimension(3) :: buff
    !! Row buffer for swapping.

    a123(:) = [ 1, 2, 3 ]
    istart=1
    do j = 1, 3
       do i = j, dim
          loc=minloc(abs(arr(i:dim,a123(1))),dim=1,mask=(abs(arr(i:dim,a123(1))).gt.1.E-5_real32))+i-1
          buff(:)=arr(i,:)
          arr(i,:)=arr(loc,:)
          arr(loc,:)=buff(:)
       end do

       scndrow: do i = j, dim
          if(abs(arr(j,a123(1))).ne.abs(arr(i,a123(1)))) exit scndrow
          loc=minloc(abs(arr(i:dim,a123(2)))+abs(arr(i:dim,a123(3))),dim=1,&
               mask=(abs(arr(j,a123(1))).eq.abs(arr(i:dim,a123(1)))))+i-1
          buff(:)=arr(i,:)
          arr(i,:)=arr(loc,:)
          arr(loc,:)=buff(:)
       end do scndrow

       a123=cshift(a123,1)
    end do

    return
  end subroutine sort2D
!###############################################################################


!###############################################################################
  subroutine sort_col(arr1,col,reverse)
    !! Sort a 2D array by a specified column.
    implicit none

    ! Arguments
    real(real32), dimension(:,:), intent(inout) :: arr1
    !! Array to sort in-place.
    integer, intent(in) :: col
    !! Column index to sort by.
    logical, optional, intent(in) :: reverse
    !! If true, sort in descending order.

    ! Local variables
    integer :: i, dim, loc
    !! Loop index, array dimension, swap location.
    logical :: udef_reverse
    !! Local copy of reverse flag.
    real(real32), allocatable, dimension(:) :: dbuff
    !! Row buffer for swapping.


    if(present(reverse))then
       udef_reverse=reverse
    else
       udef_reverse=.false.
    end if

    allocate(dbuff(size(arr1,dim=2)))

    dim=size(arr1,dim=1)
    do i=1,dim
       if(udef_reverse)then
          loc=maxloc(arr1(i:dim,col),dim=1)+i-1          
       else
          loc=minloc(arr1(i:dim,col),dim=1)+i-1
       end if
       dbuff=arr1(i,:)
       arr1(i,:)=arr1(loc,:)
       arr1(loc,:)=dbuff

    end do

    return
  end subroutine sort_col
!###############################################################################


!###############################################################################
  subroutine loadbar(count,div,loaded)
    !! Write a loading bar to the terminal.
    implicit none

    ! Arguments
    integer, intent(in) :: count
    !! Current iteration count.
    integer, intent(in) :: div
    !! Division interval for dot printing.
    character(1), optional, intent(in) :: loaded
    !! If 'l' or 'y', clear the loading bar.

    ! Local variables
    real(real32) :: tiny = 1.E-5_real32
    !! Small threshold for mod comparison.
    character(1) :: yn
    !! Local copy of loaded flag.
    character(1) :: creturn = achar(13)
    !! Carriage return character.

    if(.not.present(loaded)) then
       yn='n'
    else
       yn=loaded
    end if

    if(yn.eq.'l'.or.yn.eq.'y') then
       write(*,'(A,20X,A)',advance='no') achar(13),achar(13)
       return
    end if

    if((real(count)/real(4*div)-floor(real(count)/real(4*div))).lt.tiny) then
       write(*,'(A,20X,A,"CALCULATING")',advance='no') creturn,creturn
    else if((real(count)/real(div)-floor(real(count)/real(div))).lt.tiny) then
       write(*,'(".")',advance='no')
    end if

    return
  end subroutine loadbar
!###############################################################################

end module artemis__misc
