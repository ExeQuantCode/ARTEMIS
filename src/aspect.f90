!!!#############################################################################
!!! Module to define all global variables
!!! Code written by:
!!!    Ned Thaddeus Taylor
!!! Code part of the ARTEMIS group
!!!#############################################################################
module aspect
  use io
  use rw_geom, only: bas_type,clone_bas
  use edit_geom
  implicit none

  private

  integer, parameter :: nopt_edits=4

  integer, parameter :: ishift_index=1
  integer, parameter :: ivacuum_index=2
  integer, parameter :: itransform_index=3
  integer, parameter :: ishift_region_index=4
  integer, parameter :: islab_index=5


  type aspect_type
     integer :: nedits
     integer, dimension(nopt_edits) :: list  !! lists order of edits to perform
     integer, dimension(nopt_edits) :: axis
     double precision, dimension(nopt_edits) :: val  !BOUNDS FOR EACH?
     double precision, dimension(nopt_edits,2) :: bounds  !BOUNDS FOR EACH?
     double precision, dimension(3,3) :: tfmat    
  end type aspect_type



  public :: aspect_type
  public :: edit_structure


!!!updated  2019/11/17


contains
!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine edit_structure(lat,bas,ofile,edits)
    implicit none
    integer :: GEOMunit,i
    type(bas_type) :: edited_bas
    double precision, dimension(3,3) :: edited_lat
    character(len=*), intent(in) :: ofile
    type(bas_type), intent(in) :: bas
    type(aspect_type), intent(in) :: edits
    double precision, dimension(3,3), intent(in) :: lat


!!! TAKE ORDER OF TASKS FROM THE ARTEMIS USER INPUT
!!! MAKE IT A CUSTOM STRUCTURE, OR JUST A LIST?
!!! THEN TAKE TFMAT,VACUUM,SHIFT,SHIFT REGION, (ROTATE?)
!!! YEAH, STORE THIS AS AN ASPECT CUSTOM STRUCTURE, CONTAINS LIST AND ALL OF THIS
!!! PUT CUSTOM STRUCTURE ELSEWHERE, BUT WRITING IT HERE FOR NOW

    call clone_bas(&
         inbas=bas,outbas=edited_bas,&
         inlat=lat,outlat=edited_lat)

    do i=1,edits%nedits
       
       select case(edits%list(i))
       case(ishift_index)
          call shifter(edited_bas,edits%axis(i),edits%val(i))
       case(ishift_region_index)
          call err_abort('ERROR: SHIFT REGION NOT YET SET UP. ISSUE WITH BOUNDS')          
       case(ivacuum_index)
          call vacuumer(edited_lat,edited_bas,&
               edits%axis(i),edits%bounds(i,1),edits%val(i))
       case(itransform_index)
          call transformer(lat=edited_lat,bas=edited_bas,tfmat=edits%tfmat)
       case(islab_index)
          call err_abort('ERROR: SLAB PRINTER NOT YET SET UP')
       end select

    end do
    
    GEOMunit=101
    open(unit=GEOMunit,file=trim(ofile))
    call geom_write(GEOMunit,edited_lat,edited_bas)
    close(GEOMunit)
    

  end subroutine edit_structure
!!!#############################################################################

end module aspect
