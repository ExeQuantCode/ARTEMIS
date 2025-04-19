!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! module contains various miscellaneous functions and subroutines.
!!! module includes the following functions and subroutines:
!!! sort1D           (sort 1st col of array by size. Opt:sort 2nd array wrt 1st)
!!! sort2D           (sort 1st two columns of an array by size)
!!! set              (return the sorted set of unique elements)
!!! sort_col         (sort array with respect to col column)
!!! swap_i           (swap two integers around)
!!! swap_d           (swap two doubles around)
!!! swap_vec         (swap two vectors around)
!!!##################
!!! Icount           (counts words on line)
!!! grep             (finds 1st line containing the pattern)
!!! flagmaker        (read flag inputs supplied and stores variable if present)
!!! loadbar          (writes out a loading bar to the terminal)
!!! jump             (moves file to specified line number)
!!! file_check       (checks whether file exists and prompts user otherwise)
!!! to_upper         (converts all characters in string to upper case)
!!! to_lower         (converts all characters in string to lower case)
!!!#############################################################################
module artemis__misc
  use artemis__constants, only: real32
  implicit none


  interface swap
     procedure iswap, rswap, rswap_vec
  end interface swap

  interface sort1D
     procedure isort1D,rsort1D
  end interface sort1D

  interface set
     procedure iset,rset
  end interface set



contains

!!!#####################################################
!!! sorts two arrays from min to max
!!! sorts the optional second array wrt the first array
!!!#####################################################
  subroutine isort1D(arr1,arr2,reverse)
    implicit none
    integer :: i,dim,loc
    integer :: ibuff
    logical :: udef_reverse
    integer, dimension(:) :: arr1
    integer, dimension(:),intent(inout),optional :: arr2
    logical, optional, intent(in) :: reverse

    if(present(reverse))then
       udef_reverse=reverse
    else
       udef_reverse=.false.
    end if

    dim=size(arr1,dim=1)
    do i=1,dim
       if(udef_reverse)then
          loc=maxloc(arr1(i:dim),dim=1)+i-1          
       else
          loc=minloc(arr1(i:dim),dim=1)+i-1
       end if
       ibuff=arr1(i)
       arr1(i)=arr1(loc)
       arr1(loc)=ibuff

       if(present(arr2)) then
          ibuff=arr2(i)
          arr2(i)=arr2(loc)
          arr2(loc)=ibuff
       end if
    end do

    return
  end subroutine isort1D
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rsort1D(arr1,arr2,reverse)
    implicit none
    integer :: i,dim,loc,ibuff
    real(real32) :: dbuff
    logical :: udef_reverse
    real(real32), dimension(:) :: arr1
    integer, dimension(:),intent(inout),optional :: arr2
    logical, optional, intent(in) :: reverse

    if(present(reverse))then
       udef_reverse=reverse
    else
       udef_reverse=.false.
    end if

    dim=size(arr1,dim=1)
    do i=1,dim
       if(udef_reverse)then
          loc=maxloc(arr1(i:dim),dim=1)+i-1          
       else
          loc=minloc(arr1(i:dim),dim=1)+i-1
       end if
       dbuff=arr1(i)
       arr1(i)=arr1(loc)
       arr1(loc)=dbuff

       if(present(arr2)) then
          ibuff=arr2(i)
          arr2(i)=arr2(loc)
          arr2(loc)=ibuff
       end if
    end do

    return
  end subroutine rsort1D
!!!#####################################################


!!!#####################################################
!!! sort an array from min to max
!!!#####################################################
  subroutine sort2D(arr,dim)
    implicit none
    integer :: i,j,dim,loc,istart
    integer, dimension(3) :: a123
    real(real32), dimension(3) :: buff
    real(real32), dimension(dim,3) :: arr

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
!!!#####################################################


!!!#####################################################
!!! return the sorted set of unique elements
!!!#####################################################
  subroutine iset(arr)
    implicit none
    integer :: i,n
    integer, allocatable, dimension(:) :: tmp_arr
    
    integer, allocatable, dimension(:) :: arr

    call sort1D(arr)
    allocate(tmp_arr(size(arr)))

    tmp_arr(1) = arr(1)
    n=1
    do i=2,size(arr)
       if(arr(i)==tmp_arr(n)) cycle
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    call move_alloc(tmp_arr, arr)
    
  end subroutine iset
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rset(arr, tol)
    implicit none
    integer :: i,n
    real(real32) :: tiny
    real(real32), allocatable, dimension(:) :: tmp_arr
    
    real(real32), allocatable, dimension(:) :: arr
    real(real32), optional :: tol

    if(present(tol))then
       tiny = tol
    else
       tiny = 1.E-4_real32
    end if
    
    call sort1D(arr)
    allocate(tmp_arr(size(arr)))

    tmp_arr(1) = arr(1)
    n=1
    do i=2,size(arr)
       if(abs(arr(i)-tmp_arr(n)).lt.tiny) cycle
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    call move_alloc(tmp_arr, arr)
    
  end subroutine rset
!!!#####################################################


!!!#####################################################
!!! sort an array over specified column
!!!#####################################################
!!! Have it optionally take in an integer vector that ...
!!! ... lists the order of imporance of columns
  subroutine sort_col(arr1,col,reverse)
    implicit none
    integer :: i,dim,loc
    logical :: udef_reverse
    real(real32), allocatable, dimension(:) :: dbuff
    real(real32), dimension(:,:) :: arr1

    integer, intent(in) :: col
    logical, optional, intent(in) :: reverse


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
!!!#####################################################


!!!#####################################################
!!! swap two ints
!!!#####################################################
  subroutine iswap(i1,i2)
    implicit none
    integer :: i1,i2,itmp

    itmp=i1
    i1=i2
    i2=itmp
  end subroutine iswap
!!!#####################################################


!!!#####################################################
!!! swap two doubles
!!!#####################################################
  subroutine rswap(d1,d2)
    implicit none
    real(real32) :: d1,d2,dtmp

    dtmp=d1
    d1=d2
    d2=dtmp
  end subroutine rswap
!!!#####################################################


!!!#####################################################
!!! swap two vectors
!!!#####################################################
  subroutine rswap_vec(vec1,vec2)
    implicit none
    real(real32),dimension(:)::vec1,vec2
    real(real32),allocatable,dimension(:)::tvec

    allocate(tvec(size(vec1)))
    tvec=vec1(:)
    vec1(:)=vec2(:)
    vec2(:)=tvec
  end subroutine rswap_vec
!!!#####################################################


!!!#############################################################################
!!!#############################################################################
!!!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!!!#############################################################################
!!!#############################################################################


!!!#####################################################
!!! counts the number of words on a line
!!!#####################################################
  integer function Icount(full_line,tmpchar)
    character(*) :: full_line
    !ONLY WORKS WITH IFORT COMPILER
    !      character(1) :: fs
    character(len=:),allocatable :: fs
    character(100),optional :: tmpchar
    integer ::items,pos,k,length
    items=0
    pos=1

    length=1
    if(present(tmpchar)) length=len(trim(tmpchar))
    allocate(character(len=length) :: fs)
    if(present(tmpchar)) then
       fs=trim(tmpchar)
    else
       fs=" "
    end if

    loop: do
       k=verify(full_line(pos:),fs)
       if (k.eq.0) exit loop
       items=items+1
       pos=k+pos-1
       k=scan(full_line(pos:),fs)
       if (k.eq.0) exit loop
       pos=k+pos-1
    end do loop
    Icount=items
  end function Icount
!!!#####################################################


!!!#####################################################
!!! grep 
!!!#####################################################
!!! searches a file untill it finds the mattching patern
  subroutine grep(unit,input)
    integer :: unit,Reason
    character(*) :: input
    character(1024) :: buffer
    !  character(1024), intent(out), optional :: linechar
    rewind(unit)
    greploop: do
       read(unit,'(A100)',iostat=Reason) buffer
       if(Reason.lt.0) return
       if(index(trim(buffer),trim(input)).ne.0) exit greploop
    end do greploop
  end subroutine grep
!!!#####################################################


!!!#####################################################
!!! Assigns variables of flags from getarg
!!!#####################################################
!!! SHOULD MAKE THIS A FUNCTION INSTEAD !!!
  subroutine flagmaker(buffer,flag,i,skip,empty)
    integer :: i
    logical :: skip,empty
    character(*) :: flag,buffer

    if(len(trim(buffer)).eq.len(trim(flag))) then
       call getarg(i+1,buffer)
       if(scan(buffer,'-').eq.1.or.buffer.eq.'') then
          buffer=""
          empty=.true.
       else
          skip=.true.
       end if
    else
       buffer=buffer(len(trim(flag))+1:)
    end if

    return
  end subroutine flagmaker
!!!#####################################################


!!!#####################################################
!!! Writes out a loading bar to the terminal
!!!#####################################################
  subroutine loadbar(count,div,loaded)
    implicit none
    integer :: count,div !div=10
    real(real32) :: tiny=1.E-5
    character(1) :: yn,creturn = achar(13)
    character(1), optional :: loaded

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
!!!#####################################################


!!!#####################################################
!!! Jumps UNIT to input line number
!!!#####################################################
  subroutine jump(unit,linenum)
    integer :: unit, linenum, move
    rewind(unit)
    do move=1,(linenum)
       read(unit,*)
    end do
    return
  end subroutine jump
!!!#####################################################


!!!#####################################################
!!! File checker
!!!#####################################################
  subroutine file_check(UNIT,FILENAME,ACTION)
    implicit none
    integer :: i,UNIT,Reason
    character(len=*) :: FILENAME
    character(20) :: udef_action
    character(20), optional :: ACTION
    logical :: filefound

    udef_action="READWRITE"
    if(present(ACTION)) udef_action=ACTION
    udef_action=to_upper(udef_action)
    do i=1,5
       inquire(file=trim(FILENAME),exist=filefound)
       if(.not.filefound) then
          write(*,'("File name ",A," not found.")')&
               "'"//trim(FILENAME)//"'"
          write(*,'("Supply another filename: ")')
          read(*,*) FILENAME
       else
          write(*,'("Using file ",A)')  &
               "'"//trim(FILENAME)//"'"
          exit
       end if
       if(i.ge.4) then
          write(*,*) "Nope"
          call exit()
       end if
    end do
    if(trim(adjustl(udef_action)).eq.'NONE')then
       write(*,*) "File found, but not opened."
    else
       open(unit=UNIT,file=trim(FILENAME),action=trim(udef_action),iostat=Reason)
    end if


    return
  end subroutine file_check
!!!#####################################################


!!!#####################################################
!!! converts all characters in string to upper case
!!!#####################################################
  function to_upper(buffer) result(upper)
    implicit none
    integer :: i,j
    character(*) :: buffer
    character(len=:),allocatable :: upper


    allocate(character(len=len(buffer)) :: upper)
    do i=1,len(buffer)
       j=iachar(buffer(i:i))
       if(j.ge.iachar("a").and.j.le.iachar("z"))then
          upper(i:i)=achar(j-32)
       else
          upper(i:i)=buffer(i:i)
       end if
    end do

    return
  end function to_upper
!!!#####################################################


!!!#####################################################
!!! converts all characters in string to lower case
!!!#####################################################
  function to_lower(buffer) result(lower)
    implicit none
    integer :: i,j
    character(*) :: buffer
    character(len=:),allocatable :: lower


    allocate(character(len=len(buffer)) :: lower)
    do i=1,len(buffer)
       j=iachar(buffer(i:i))
       if(j.ge.iachar("A").and.j.le.iachar("Z"))then
          lower(i:i)=achar(j+32)
       else
          lower(i:i)=buffer(i:i)
       end if
    end do

    return
  end function to_lower
!!!#####################################################


!###############################################################################
  function strip_null(buffer) result(stripped)
    !! Strip null characters from a string.
    !!
    !! This is meant for handling strings passed from Python, which gain
    !! null characters at the end. The procedure finds the first null
    !! character and truncates the string at that point.
    !! Null characters are represented by ASCII code 0.
    implicit none

    ! Arguments
    character(*), intent(in) :: buffer
    !! String to be stripped.
    character(len=len(buffer)) :: stripped
    !! Stripped string.

    ! Local variables
    integer :: i
    !! Loop index.

    stripped = ""
    do i = 1, len(buffer)
       if(iachar(buffer(i:i)).ne.0)then
          stripped(i:i)=buffer(i:i)
       else
          exit
       end if
    end do

  end function strip_null
!###############################################################################

end module artemis__misc
