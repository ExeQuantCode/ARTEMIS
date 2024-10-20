!!!#############################################################################
!!! Module to define all global variables
!!! Code written by:
!!!    Ned Thaddeus Taylor
!!! Code part of the ARTEMIS group
!!!#############################################################################
module io
  use misc
  implicit none
  
  private   !everything is private unless explicitly defined as public

  character(25), public, parameter :: version="development version 1.0.3a"
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
     character(40) :: summary
     character(60) :: allowed
     character(60) :: default
     character(1024) :: description
  end type tag_type

  public :: write_fmtd
  public :: err_abort,print_warning
  public :: err_abort_print_struc
  public :: io_print_help
  public :: print_header
  public :: setup_input_fmt,setup_output_fmt


!!!updated 2021/11/11


contains
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
    write(unit,'(A,A)') " Welcome to ARTEMIS version ",version
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
    write(unit,'(A)') " LICENCE:"
    write(unit,'(A)') " This work is licensed under a &
         &Creative Commons Attribution-NonCommercial 3.0 &
         &Unported (CC BY-NC 3.0) License."
    write(unit,'(A)') " https://creativecommons.org/licenses/by-nc/3.0/"
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
          call write_fmtd(unit,trim(message))
          lpresent=.true.
       end if
    end if
    if(.not.lpresent) write(unit,'(A)') trim(message)
    stop

  end subroutine err_abort
!!!#############################################################################


!!!#############################################################################
!!! Prints to stderr, prints structure and stops
!!!#############################################################################
  subroutine err_abort_print_struc(in_lat,in_bas,name,message,lstop)
    use rw_geom
    implicit none
    integer :: unit=0
    character(len=*) :: name,message
    type(bas_type) :: in_bas
    double precision, dimension(3,3) :: in_lat
    logical, optional :: lstop

    
    open(100,file=name)
    call geom_write(100,in_lat,in_bas)
    close(100)
    if(message.ne.'') write(unit,'(A)') trim(message)
    if(.not.present(lstop).or.lstop) stop

  end subroutine err_abort_print_struc
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

                write(unit,'(A,T33,A)') &
                     trim(tags(i)%name),trim(tags(i)%summary)

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
  

!!!#############################################################################
!!! sets up the name of output files and subroutines to read files
!!!#############################################################################
  subroutine setup_input_fmt(fmt)
    use rw_geom, only : igeom_input
    implicit none
    character(len=*), intent(in) :: fmt
    character(len=:), allocatable :: form
    

    allocate(character(len=len(trim(adjustl(fmt)))) ::  form)
    form = trim(adjustl(to_upper(fmt)))
    
    select case(form)
    case("VASP")
       write(6,*) "Input files will be VASP formatted"
       igeom_input=1
    case("CASTEP")
       write(6,*) "Input files will be CASTEP formatted"
       igeom_input=2
       !call err_abort('ERROR: ARTEMIS not yet set up for CASTEP')
    case("QE","QUANTUMESPRESSO")
       write(6,*) "Input files will be QuantumEspresso formatted"
       igeom_input=3
       !call err_abort('ERROR: ARTEMIS not yet set up for Quantum Espresso')
    case("CRYSTAL")
       write(6,*) "Input files will be CRYSTAL formatted"
       igeom_input=4
       call err_abort('ERROR: ARTEMIS not yet set up for CRYSTAL')
    end select

    
  end subroutine setup_input_fmt
!!!#############################################################################
  

!!!#############################################################################
!!! sets up the name of output files and subroutines to read files
!!!#############################################################################
  subroutine setup_output_fmt(fmt,out_filename)
    use rw_geom, only : igeom_output
    implicit none
    character(len=*) :: out_filename
    character(len=*), intent(in) :: fmt
    character(len=:), allocatable :: form
    

    allocate(character(len=len(trim(adjustl(fmt)))) ::  form)
    form = trim(adjustl(to_upper(fmt)))
    
    select case(form)
    case("VASP")
       write(6,*) "Output files will be VASP formatted"
       if(out_filename.eq.'') out_filename="POSCAR"
       igeom_output=1
    case("CASTEP")
       write(6,*) "Output files will be CASTEP formatted"
       if(out_filename.eq.'') out_filename="struc.cell"
       igeom_output=2
       !call err_abort('ERROR: ARTEMIS not yet set up for CASTEP')
    case("QE","QUANTUMESPRESSO")
       write(6,*) "Output files will be QuantumEspresso formatted"
       if(out_filename.eq.'') out_filename="struc.geom"
       igeom_output=3
       !call err_abort('ERROR: ARTEMIS not yet set up for Quantum Espresso')
    case("CRYSTAL")
       write(6,*) "Output files will be CRYSTAL formatted"
       if(out_filename.eq.'') out_filename="INPUT_geom"
       igeom_output=4
       call err_abort('ERROR: ARTEMIS not yet set up for CRYSTAL')
    end select

    
  end subroutine setup_output_fmt
!!!#############################################################################

end module io
