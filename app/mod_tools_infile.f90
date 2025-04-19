!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! module contains customn input file reading functions and subroutines.
!!! module includes the following functionsand subroutines:
!!! val              (outputs contents of string occuring after "=")
!!! getline          (gets the line using grep and goes back to start of line)
!!! assignI          (assign an integer to variable)
!!! assignD          (assign a real(real32) to variable)
!!! assignIvec       (assign an arbitrary length vector of integers to variable)
!!! assignDvec       (assign an arbitrary length vector of DP to variable)
!!! assignS          (assign a string to variable)
!!! assignL          (assign a logical to variable (T/t/1 and F/f/0 accepted)
!!! assign_list      (assign output to a user-defined value in a list)
!!! assign_listvec   (assign output vector to a user-defined value in a list)
!!! rm_comments      (remove comments from a string (anything after ! or #))
!!! cat              (cat lines until user-defined end string is encountered)
!!!#############################################################################
module infile_tools
  use artemis__constants, only: real32
  use artemis__misc, only: grep,icount
  implicit none

  interface assign
     procedure assignI,assignR,assignS,assignL
  end interface assign
  interface assign_vec
     procedure assignIvec,assignRvec
  end interface assign_vec


!!!updated 2020/01/22


contains
!!!#############################################################################
!!! val outputs the section of buffer that occurs after an "="
!!!#############################################################################
  function val(buffer)
    character(*) :: buffer
    character(100) :: val
    val=trim(adjustl(buffer((scan(buffer,"=",back=.false.)+1):)))
    return
  end function val
!!!#############################################################################


!!!#############################################################################
!!! gets the line from a grep and assigns it to buffer
!!!#############################################################################
  subroutine getline(unit,pattern,buffer)
    integer :: unit,Reason
    character(*) :: pattern,buffer
    
    call grep(unit,pattern)
    backspace(unit);read(unit,'(A)',iostat=Reason) buffer
    
  end subroutine getline
!!!#############################################################################


!!!#############################################################################
!!! assigns an integer to variable
!!!#############################################################################
  subroutine assignI(buffer,variable,found)
    integer :: found
    character(1024) :: buffer1,buffer2
    character(*) :: buffer
    integer :: variable
    buffer1=buffer(:scan(buffer,"=")-1)
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       read(buffer2,*) variable
    end if
  end subroutine assignI
!!!#############################################################################


!!!#############################################################################
!!! assigns a DP value to variable
!!!#############################################################################
  subroutine assignIvec(buffer,variable,found)
    integer :: found,i
    character(1024) :: buffer1,buffer2
    character(*) :: buffer
    integer, dimension(:) :: variable
    buffer1=buffer(:scan(buffer,"=")-1)
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       read(buffer2,*) (variable(i),i=1,size(variable))
    end if
  end subroutine assignIvec
!!!#############################################################################


!!!#############################################################################
!!! assigns a DP value to variable if the line contains the right keyword
!!!#############################################################################
  subroutine assignR(buffer,variable,found)
    integer :: found
    character(1024) :: buffer1,buffer2
    character(*) :: buffer
    real(real32) :: variable
    buffer1=buffer(:scan(buffer,"=")-1)
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       read(buffer2,*) variable
    end if
  end subroutine assignR
!!!#############################################################################


!!!#############################################################################
!!! assigns a DP value to variable
!!!#############################################################################
  subroutine assignRvec(buffer,variable,found)
    integer :: found,i
    character(1024) :: buffer1,buffer2
    character(*) :: buffer
    real(real32), dimension(:) :: variable
    buffer1=buffer(:scan(buffer,"=")-1)
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       read(buffer2,*) (variable(i),i=1,size(variable))
    end if
  end subroutine assignRvec
!!!#############################################################################


!!!#############################################################################
!!! assigns a string
!!!#############################################################################
  subroutine assignS(buffer,variable,found)
    integer::found
    character(1024)::buffer1,buffer2
    character(*) :: buffer
    character(*) :: variable
    buffer1=buffer(:scan(buffer,"=")-1)
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       read(buffer2,'(A)') variable
    end if
  end subroutine assignS
!!!#############################################################################


!!!#############################################################################
!!! assigns a logical
!!!#############################################################################
  subroutine assignL(buffer,variable,found)
    integer::found
    character(1024)::buffer1,buffer2
    character(*)::buffer
    logical::variable
    buffer1=buffer(:scan(buffer,"=")-1)
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       if(index(buffer2,"T").ne.0.or.&
            index(buffer2,"t").ne.0.or.&
            index(buffer2,"1").ne.0) then
          variable=.TRUE.
       end if
       if(index(buffer2,"F").ne.0.or.&
            index(buffer2,"f").ne.0.or.&
            index(buffer2,"0").ne.0) then
          variable=.FALSE.
       end if
    end if
  end subroutine assignL
!!!#############################################################################


!!!#############################################################################
!!! assigns a set
!!!#############################################################################
  function assign_list(buffer,tag_list,num) result(var)
    implicit none
    integer :: nlist,loc2,i
    real(real32) :: var
    character(len=1024) :: new_buffer
    integer, allocatable, dimension(:) :: loc_list
    integer, intent(in) :: num
    character(len=*), intent(in) :: buffer
    character(len=*), dimension(:), intent(in) :: tag_list


    nlist=size(tag_list,dim=1)
    allocate(loc_list(nlist))
    do i=1,nlist
       loc_list(i)=index(buffer,trim(adjustl(tag_list(i))))
    end do
    loc2 = minloc(&
         loc_list(:)-loc_list(num),dim=1,&
         mask=loc_list(:)-loc_list(num).gt.0)
    if(loc2.eq.0)then
       loc2=len(buffer)
    else
       loc2=loc_list(loc2)-1
    end if
    new_buffer = buffer(loc_list(num):loc2)

    new_buffer = trim(new_buffer(scan(new_buffer,'=')+1:))
    read(new_buffer,*) var


  end function assign_list
!!!#############################################################################


!!!#############################################################################
!!! assigns a set
!!!#############################################################################
  function assign_listvec(buffer,tag_list,num) result(var)
    implicit none
    integer :: nlist,loc2,i
    real(real32), allocatable, dimension(:) :: var
    character(len=1024) :: new_buffer
    integer, allocatable, dimension(:) :: loc_list
    integer, intent(in) :: num
    character(len=*), intent(in) :: buffer
    character(len=*), dimension(:), intent(in) :: tag_list


    nlist=size(tag_list,dim=1)
    allocate(loc_list(nlist))
    do i=1,nlist
       loc_list=index(buffer,trim(adjustl(tag_list(i))))
    end do
    loc2 = minloc(&
         loc_list(:)-loc_list(num),dim=1,&
         mask=loc_list(:)-loc_list(num).gt.0)-1
    new_buffer = buffer(loc_list(num):loc_list(loc2))

    new_buffer = trim(new_buffer(scan(new_buffer,'=')+1:))
    
    allocate(var(icount(new_buffer)))
    read(new_buffer,*) var


  end function assign_listvec
!!!#############################################################################


!!!#############################################################################
!!! remove comment from a string (anything after ! or #) 
!!!#############################################################################
  subroutine rm_comments(buffer,itmp)
    implicit none
    integer :: lbracket,rbracket,iline
    character(*) :: buffer
    integer, optional :: itmp

    iline=0
    if(present(itmp)) iline=itmp

    if(scan(buffer,'!').ne.0) buffer=buffer(:(scan(buffer,'!')-1))
    if(scan(buffer,'#').ne.0) buffer=buffer(:(scan(buffer,'#')-1))
    do while(scan(buffer,'(').ne.0.or.scan(buffer,')').ne.0)
       lbracket=scan(buffer,'(',back=.true.)
       rbracket=scan(buffer(lbracket:),')')
       if(lbracket.eq.0.or.rbracket.eq.0)then
          write(*,'(A,I0)') &
               ' NOTE: a bracketing error was encountered on line ',iline
          buffer=""
          return
       end if
       rbracket=rbracket+lbracket-1
       buffer=buffer(:(lbracket-1))//buffer((rbracket+1):)
    end do


    return
  end subroutine rm_comments
!!!#############################################################################


!!!#############################################################################
!!! cat lines
!!!#############################################################################
  subroutine cat(UNIT,end_string,string,line,end_string2,rm_cmt)
    implicit none
    integer :: line,Reason
    character(len=1024) :: buffer
    logical :: lfound_end
    integer, intent(in) :: UNIT
    character(len=*), intent(in) :: end_string
    character(len=*), optional, intent(in) ::end_string2
    character(len=*), intent(inout) :: string
    logical, optional, intent(in) :: rm_cmt

    lfound_end=.false.
    cat_loop: do
       line = line + 1
       read(UNIT,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit cat_loop
       if(present(rm_cmt).and.rm_cmt) call rm_comments(buffer,line)
       if(trim(buffer).eq.'') cycle cat_loop
       endstring_if: if(index(trim(adjustl(buffer)),end_string).eq.1)then
          if(present(end_string2))then
             if(index(trim(buffer),end_string2).eq.0) exit endstring_if
          end if
          lfound_end=.true.
          exit cat_loop
       end if endstring_if
       string=trim(string)//' '//trim(buffer)
    end do cat_loop
    if(end_string.eq.'') lfound_end=.true.
    if(.not.lfound_end)then
       write(0,'("ERROR: Error in cat subroutine in mod_tools_infile.f90")')
       write(0,'(1X,"No END tag found for cat")')
       stop
    end if
    
    
  end subroutine cat
!!!#############################################################################
end module infile_tools
