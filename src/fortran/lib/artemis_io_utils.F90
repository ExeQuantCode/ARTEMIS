module artemis__io_utils
  !! Module for I/O utilities, error handling, and version information.
  !!
  !! Provides formatted output, warning/error printing with rich box formatting,
  !! help system for parameter tags, and the ARTEMIS version string.
  use artemis__constants, only: real32
  use coreutils__string, only: to_upper
  use coreutils__error, only: stop_program
  implicit none
  

  private

  public :: write_fmtd
  public :: err_abort, print_warning, stop_program
  public :: artemis__suppress_warnings
  public :: io_print_help
  public :: print_header
  public :: artemis__version__


  logical :: artemis__suppress_warnings = .false.
  !! If true, suppress all warning messages.
  character(len=*), parameter :: artemis__version__ = "2.0.0"
  !! ARTEMIS version string.


  type, public :: tag_type
     !! Type for storing parameter tag metadata for the help system.
     character(25) :: name
     !! Tag name.
     character(1)  :: type
     !! Tag data type code (I, R, S, L, U, V, B).
     character(50) :: summary
     !! Short summary of the tag.
     character(60) :: allowed
     !! Allowed values string.
     character(60) :: default
     !! Default value string.
     character(1024) :: description
     !! Full description of the tag.
     logical :: is_deprecated = .false.
     !! Whether the tag is deprecated.
     logical :: to_be_deprecated = .false.
     !! Whether the tag will be deprecated in a future version.
     character(25) :: deprecated_name = ''
     !! New tag name replacing the deprecated one.
     character(20) :: deprecated_version
     !! Version in which the tag was/will be deprecated.
  end type tag_type



contains


!###############################################################################
  subroutine print_header(unit)
    !! Print the ARTEMIS logo and author list.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! Output unit number.

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
!###############################################################################


!###############################################################################
  subroutine write_fmtd(unit,message)
    !! Write a message with embedded \n newline formatting.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! Output unit number.
    character(len=*), intent(in) :: message
    !! Message string (may contain \n for newlines).

    ! Local variables
    integer :: istart, iend, itmp1
    !! Parsing indices and safety counter.
    
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
!###############################################################################


!###############################################################################
  subroutine print_warning(message,width,fmtd)
    !! Print a warning message in a formatted box.
    implicit none

    ! Arguments
    character(len=*), intent(in) :: message
    !! Warning message text.
    integer, optional, intent(in) :: width
    !! Box width (default: 50).
    logical, optional, intent(in) :: fmtd
    !! If true, use write_fmtd for the message lines.

    ! Local variables
    integer :: unit = 0
    !! Output unit (stderr).
    integer :: ipos, iend, inewline
    !! Parsing positions.
    integer :: whitespacel, whitespacer, length, nwidth
    !! Formatting widths.
    character(len=13) :: warning
    !! Warning header string.
    character(len=200) :: fmt
    !! Format string buffer.
    logical :: finished, lpresent
    !! Loop control flags.
    character(len=:), allocatable :: line
    !! Formatted line buffer.


    if(artemis__suppress_warnings) return

    !---------------------------------------------------------------------------
    ! Initialise variables and allocate line length
    !---------------------------------------------------------------------------
    ipos=0
    iend=0
    nwidth=50
    finished=.false.
    if(present(width)) nwidth=width
    allocate(character(len=nwidth) :: line)


    !---------------------------------------------------------------------------
    ! Print warning header
    !---------------------------------------------------------------------------
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


    !---------------------------------------------------------------------------
    ! Print the message body
    !---------------------------------------------------------------------------
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
!###############################################################################


!###############################################################################
  subroutine err_abort(message,fmtd)
    !! Print an error message to stderr and stop execution.
    implicit none

    ! Arguments
    character(len=*), intent(in) :: message
    !! Error message to print.
    logical, optional, intent(in) :: fmtd
    !! If true, use write_fmtd formatting.

    ! Local variables
    integer :: unit = 0
    !! Output unit (stderr).
    logical :: lpresent
    !! Format flag.

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
!###############################################################################


!###############################################################################
  subroutine io_print_help(unit, helpword, tags, search)
    !! Print help information for parameter tags, with search support.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! Output unit number.
    character(len=*), intent(in) :: helpword
    !! Tag name or search term to look up.
    type(tag_type), dimension(:), intent(in) :: tags
    !! Array of tag definitions.
    logical, optional, intent(in) :: search
    !! If true, perform substring search across all tags.

    ! Local variables
    integer :: i, ntags
    !! Loop index and tag count.
    character(len=15) :: type, fmt
    !! Type label and format buffer.
    character(len=:), allocatable :: checkword
    !! Upper-cased version of helpword.
    character(len=200) :: title
    !! Formatted title string.
    logical :: found, lpresent
    !! Search result and format flags.
    

    ntags=size(tags)
    allocate(character(len=len(trim(adjustl(helpword)))) ::  checkword)
    checkword = trim(adjustl(to_upper(helpword)))


    !---------------------------------------------------------------------------
    ! Check that no tagname is duplicated
    !---------------------------------------------------------------------------
    if(count(tags(:)%name.eq.checkword).gt.1)then
       call err_abort('Error: helper: tagname entry duplicated')
    end if


    !---------------------------------------------------------------------------
    ! Search function
    !---------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------
    ! Help all function
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    ! Find requested tag and print its help
    !---------------------------------------------------------------------------
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
!###############################################################################

end module artemis__io_utils
