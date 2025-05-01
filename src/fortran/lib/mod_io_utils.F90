!!!#############################################################################
!!! Module to define all global variables
!!! Code written by:
!!!    Ned Thaddeus Taylor
!!! Code part of the ARTEMIS group
!!!#############################################################################
module artemis__io_utils
  use artemis__constants, only: real32
  use artemis__misc
  implicit none
  

  private

  public :: write_fmtd
  public :: err_abort,print_warning, stop_program
  public :: io_print_help
  public :: print_header
  public :: artemis__version__


  logical :: test_error_handling = .false.
  logical :: suppress_warnings = .false.
  character(len=*), parameter :: artemis__version__ = "2.0.0"
  !character(30), public, parameter :: &
  !     author(3) = [&
  !     "N. T. Taylor",&
  !     "F. H. Davies",&
  !     "I. E. M. Rudkin",&
  !     "S. P. Hepplestone"&
  !     ]
  !character(30), public, parameter :: &
  !     contributor(4) = [&
  !     "C. J. Price",&
  !     "T. H. Chan"&
  !     "J. Pitfield",&
  !     "E. A. D. Baker",&
  !     "S. G. Davies"&
  !     ]


  type, public :: tag_type
     character(25) :: name
     character(1)  :: type
     character(50) :: summary
     character(60) :: allowed
     character(60) :: default
     character(1024) :: description
     logical :: is_deprecated = .false.
     logical :: to_be_deprecated = .false.
     character(25) :: deprecated_name = ''
     character(20) :: deprecated_version
  end type tag_type



contains

!###############################################################################
  subroutine stop_program(message, exit_code, block_stop)
    !! Stop the program and print an error message.
    implicit none
    character(len=*), intent(in) :: message
    integer, intent(in), optional :: exit_code
    logical, intent(in), optional :: block_stop

    integer :: exit_code_
    logical :: block_stop_

    if(present(exit_code)) then
       exit_code_ = exit_code
    else
       exit_code_ = 1
    end if
    if(present(block_stop)) then
       block_stop_ = block_stop
    else
       block_stop_ = .false.
    end if

    write(0,*) 'ERROR: ', trim(message)
    if(.not.block_stop_)then
       if(.not.test_error_handling) then
          stop exit_code_
       end if
    end if
  end subroutine stop_program
!###############################################################################


!!!#############################################################################
!!! prints the ARTEMIS logo and author list
!!!#############################################################################
  subroutine print_header(unit)
    implicit none
    integer :: unit

    write(unit,'(A)') repeat("#",50)
    write(unit,'(A)') repeat("#",50)
    write(unit,*)
    write(unit,'(A)') "                    █████████████████████████████"
    write(unit,'(A)') "  ██   ███   █████      ███  ███  ██     ███    █"
    write(unit,'(A)') " █  █  █  █    █     ██████ █ █ █ ████ ████  ████"
    write(unit,'(A)') " ████  ████    █       ████ ██ ██ ████ █████   ██"
    write(unit,'(A)') " █  █  █ █     █     ██████ █████ ████ ███████  █"
    write(unit,'(A)') " █  █  █  █    █        ███ █████ ██     ██    ██"
    write(unit,'(A)') "                    █████████████████████████████"
    write(unit,*)
    write(unit,'(A)') repeat("#",50)
    write(unit,'(A)') repeat("#",50)
    write(unit,'(A)') "           Ab Initio Restructuring Tool           "
    write(unit,'(A)') "    Enabling Modelling of Interface Structures    "
    write(unit,*)
    write(unit,'(A,A)') " Welcome to ARTEMIS version ", artemis__version__
    write(unit,'(A,A,1X,A,A)') " (build ",__DATE__,__TIME__,")"
    write(unit,*)
    write(unit,'(A)') " Authors:"
    write(unit,'(A)') " N. T. Taylor, F. H. Davies, I. E. M. Rudkin, S. P. Hepplestone"
    write(unit,*)
    !write(unit,'(1X,A,", ")',advance="no") (author(i)(:),i=1,size(author(:)))
    !write(unit,*)
    write(unit,'(A)') " Contributors:"
    write(unit,'(A)') " C. J. Price, T. H. Chan, J. Pitfield, E. A. D. Baker, S. G. Davies"
    write(unit,*)
    write(unit,'(A)') " Artistic advisors:"
    write(unit,'(A)') " E. L. Martin"
    write(unit,*)
    write(unit,'(A)') " LICENSE:"
    write(unit,'(A)') " This work is licensed under a &
         &General Public License 3.0 (GPLv3)"
    write(unit,'(A)') " https://www.gnu.org/licenses/gpl-3.0.en.html"
    write(unit,*)
    write(unit,'(A)') repeat("#",50)

 

  end subroutine print_header
!!!#############################################################################


!!!#############################################################################
!!! customised print formatting
!!!#############################################################################
  subroutine write_fmtd(unit,message)
    implicit none
    integer :: istart,iend,itmp1
    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    
    istart=0
    iend=0
    itmp1=0
    newline_loop: do
       itmp1=itmp1+1
       if(itmp1.gt.30) call err_abort("ERROR: Internal error in write_fmtd. Too many newlines")
       istart=iend+1
       iend=index(message(istart:),'\n')+istart-2
       if(iend.lt.istart) exit newline_loop
       write(unit,'(A)') message(istart:iend)
       iend=iend+2
    end do newline_loop
    write(unit,'(A)') message(istart:)


  end subroutine write_fmtd
!!!#############################################################################


!!!#############################################################################
!!! Prints warning
!!!#############################################################################
  subroutine print_warning(message,width,fmtd)
    implicit none
    integer :: unit=6
    integer :: ipos,iend,inewline
    integer :: whitespacel,whitespacer,length,nwidth
    character(len=13) :: warning
    character(len=200) :: fmt
    character(len=*) :: message
    logical :: finished,lpresent
    character(len=:), allocatable :: line
    integer, optional, intent(in) :: width
    logical, optional, intent(in) :: fmtd


!!!-----------------------------------------------------------------------------
!!! Initialise variables and allocate line length
!!!-----------------------------------------------------------------------------
    ipos=0
    iend=0
    nwidth=50
    finished=.false.
    if(present(width)) nwidth=width
    allocate(character(len=nwidth) :: line)


!!!-----------------------------------------------------------------------------
!!! prints warning 
!!!-----------------------------------------------------------------------------
    warning="W A R N I N G"
    length=len(warning)
    whitespacel=(nwidth-length)/2-1
    whitespacer=whitespacel
    if(whitespacel+whitespacer.ne.nwidth-length-2) whitespacer=whitespacer+1
    write(fmt,'("(","""|""",",",I0,"X,A",I0,",",I0,"X,","""|""",")")') &
         whitespacel,length,whitespacer
    write(line,trim(fmt)) warning
    

    write(unit,'("+",A,"+")') repeat('-',nwidth-2)
    write(unit,'(A)') trim(line)
    write(unit,'("|",A,"|")') repeat(' ',nwidth-2)


!!!-----------------------------------------------------------------------------
!!! prints the message
!!!-----------------------------------------------------------------------------
    newline_loop: do
       ipos=iend+1
       length=len(trim(adjustl(message(ipos:))))

       if(length.le.nwidth-4)then
          finished=.true.
       else
          length=nwidth-4
       end if
       iend=ipos+length-1


       inewline=index(message(ipos:iend),'\n')
       if(inewline.eq.1)then
          iend=ipos+1
          cycle newline_loop
       elseif(inewline.ne.0)then
          finished=.false.
          iend=ipos+inewline-2
          length=inewline-1
       end if

       whitespacel=(nwidth-length)/2-1
       whitespacer=whitespacel
       if(whitespacel+whitespacer.ne.nwidth-length-2) whitespacer=whitespacer+1
       write(fmt,'("(","""|""",",",I0,"X,A",I0,",",I0,"X,","""|""",")")') &
            whitespacel,length,whitespacer
       write(line,trim(fmt)) trim(adjustl(message(ipos:iend)))

       lpresent=.false.
       if(present(fmtd))then
          if(fmtd)then
             call write_fmtd(unit,trim(line))
             lpresent=.true.
          end if
       end if
       if(.not.lpresent) write(unit,'(A)') trim(line)

       if(finished) exit newline_loop
       if(inewline.ne.0) iend=iend+2


    end do newline_loop
    write(unit,'("+",A,"+")') repeat('-',nwidth-2)


  end subroutine print_warning
!!!#############################################################################


!!!#############################################################################
!!! Prints to stderr and stops
!!!#############################################################################
  subroutine err_abort(message,fmtd)
    implicit none
    integer :: unit=0
    logical :: lpresent
    character(len=*) :: message
    logical, optional, intent(in) :: fmtd

    lpresent=.false.
    if(present(fmtd))then
       if(fmtd)then
          call write_fmtd(unit,"ERROR: "//trim(message))
          lpresent=.true.
       end if
    end if
    if(.not.lpresent) write(unit,'(A)') trim(message)
    stop

  end subroutine err_abort
!!!#############################################################################



!!!#############################################################################
!!! help and search
!!!#############################################################################
  subroutine io_print_help(unit, helpword, tags, search)
    implicit none
    integer :: i,ntags
    integer, intent(in) :: unit
    character(len=15) :: type,fmt
    character(len=*), intent(in) :: helpword
    character(len=:), allocatable :: checkword
    character(len=200) :: title
    logical :: found,lpresent
    logical, optional :: search
    type(tag_type), dimension(:), intent(in) :: tags
    

    ntags=size(tags)
    allocate(character(len=len(trim(adjustl(helpword)))) ::  checkword)
    checkword = trim(adjustl(to_upper(helpword)))


!!!-----------------------------------------------------------------------------
!!! checks that no tagname is duplicated
!!!-----------------------------------------------------------------------------
    if(count(tags(:)%name.eq.checkword).gt.1)then
       call err_abort('Error: helper: tagname entry duplicated')
    end if


!!!-----------------------------------------------------------------------------
!!! search function
!!!-----------------------------------------------------------------------------
    lpresent=.false.
    if(present(search))then
       if(search)then
          lpresent=.true.
          tagloop1: do i=1,ntags
             if(index(tags(i)%name,checkword).ne.0)then
                found=.true.

                if(tags(i)%to_be_deprecated)then
                   write(unit,'(A,T33,A)') &
                        trim(tags(i)%name),&
                        'To be deprecated ('//trim(tags(i)%deprecated_version)//')'
                elseif(tags(i)%is_deprecated)then
                   write(unit,'(A,T33,A)') &
                         trim(tags(i)%name),&
                         'Deprecated ('//trim(tags(i)%deprecated_version)//')'
                else
                   write(unit,'(A,T33,A)') &
                        trim(tags(i)%name),trim(tags(i)%summary)
                end if

             end if
          end do tagloop1
          if(.not.found) write(unit,'(3X,A)') 'No tag found'
          return
       end if
    end if
!!!-----------------------------------------------------------------------------
!!! help all function
!!!-----------------------------------------------------------------------------
    if(.not.lpresent.and.checkword.eq.'ALL')then
       tagloop2: do i=1,ntags
          write(unit,'(A,T33,A)') &
               trim(tags(i)%name),trim(tags(i)%summary)
          if(len(trim(tags(i)%summary)).gt.40)then
             write(0,'("WARNING: Internal error in io_print_help")')
             write(0,'(2X,"io_print_help in io.f90 has been supplied a&
                  & tag summary exceeding 40 characters")')
             cycle tagloop2
          end if
       end do tagloop2
       return

    end if

!!!-----------------------------------------------------------------------------
!!! finds requested tag and prints its help
!!!-----------------------------------------------------------------------------
    found=.false.
    tagloop3: do i=1,ntags
       if(trim(tags(i)%name).eq.checkword)then

          found=.true.

          title=trim(tags(i)%name)//"   --"//trim(tags(i)%summary)//"--"
          write(fmt,'("(",I0,"X,A)")') max(40-len(trim(title)),1)
          write(unit,*)
          write(unit,fmt) trim(title)
          write(unit,*)
          if(tags(i)%is_deprecated)then
             write(unit,'("DEPRECATED AS OF ",A)') &
                  trim(tags(i)%deprecated_version)
          elseif(tags(i)%to_be_deprecated)then
             write(unit,'("TO BE DEPRECATED AS OF ",A)') &
                  trim(tags(i)%deprecated_version)
          end if
          if(trim(tags(i)%deprecated_name).ne.'')then
             write(unit,'("New tag name: ",A)') trim(tags(i)%deprecated_name)
          end if
          if(tags(i)%is_deprecated.or.tags(i)%to_be_deprecated)then
             write(unit,*)
          end if

          select case(tags(i)%type)
          case('I'); type = 'Integer'
          case('R'); type = 'Real'
          case('S'); type = 'String'
          case('L'); type = 'Boolean/Logical'
          case('U'); type = 'Integer Vector'
          case('V'); type = 'Real Vector'
          case('B'); type = 'Block'
          end select

          write(unit,'("Type: ",A)') trim(type)
          write(unit,*)
          call write_fmtd(unit,trim(tags(i)%description))
          if(trim(tags(i)%allowed).ne.'') &
               write(unit,'("Allowed values: ",A)') trim(tags(i)%allowed)
          if(trim(tags(i)%default).ne.'') & 
               write(unit,'("Default value: ",A)') trim(tags(i)%default)

          exit tagloop3

       end if
    end do tagloop3
    if(.not.found) write(unit,'(3X,A)') 'No tag found'





  end subroutine io_print_help
!!!#############################################################################

end module artemis__io_utils
