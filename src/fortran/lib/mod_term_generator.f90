!!!#############################################################################
!!! INTERFACES CARD SUBROUTINES
!!! Code written by Ned Thaddeus Taylor and Isiah Edward Mikel Rudkin
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
module artemis__termination_generator
  use artemis__constants,     only: real32, ierror
  use artemis__misc_types,    only: abstract_artemis_generator_type
  use artemis__geom_rw,       only: basis_type
  use artemis__io_utils,      only: err_abort, print_warning
  use artemis__io_utils_extd, only: err_abort_print_struc
  use misc_linalg,            only: modu
  use edit_geom,              only: planecutter, transformer, reducer, &
       MATNORM, compare_stoichiometry
  use artemis__sym,              only: confine_type, gldfnd
  use artemis__terminations, only: get_termination_info, term_arr_type, set_slab_height, set_layer_tol, build_slab
  implicit none


  private

  public :: artemis_termination_generator_type


  type, extends(abstract_artemis_generator_type) :: artemis_termination_generator_type

    real(real32) :: layer_separation_cutoff = 1._real32

   contains
     procedure, pass(this) :: generate => generate_terminations
  end type artemis_termination_generator_type



contains

!###############################################################################
  subroutine generate_terminations( &
       this, basis, miller_plane, axis, surface, num_layers, thickness, &
       orthogonalise, normalise, break_on_fail &
  )
    !! Generate and prints terminations parallel to the supplied miller plane
    implicit none

    ! Arguments
    class(artemis_termination_generator_type), intent(inout) :: this
    !! Instance of artemis generator type
    type(basis_type), intent(in) :: basis
    !! Atomic structure data
    integer, dimension(3), intent(in) :: miller_plane
    !! Miller plane
    integer, intent(in) :: axis
    !! Axis along which to align the slab
    integer, dimension(:), intent(in), optional :: surface
    !! Surface termination indices
    integer, intent(in), optional :: num_layers
    !! Number of layers in the slab
    real(real32), intent(in), optional :: thickness
    !! Thickness of the slab (in Å)
    logical, intent(in), optional :: orthogonalise
    !! Boolean whether to orthogonalise the lattice
    logical, intent(in), optional :: normalise
    !! Boolean whether to normalise the lattice and basis
    logical, intent(in), optional :: break_on_fail
    !! Boolean whether to break on failure

    type(basis_type), dimension(:), allocatable :: output
    !! Output structures

    ! Local variables
    integer :: itmp1, iterm, term_start, term_end, iterm_step, i
    !! Termination loop variables
    integer :: ncells, ntrans
    !! Number of cells in the slab
    integer :: num_structures
    !! Number of structures to be generated
    integer, dimension(2) :: surface_
    !! Surface termination indices
    integer :: num_layers_
    !! Number of layers in the slab
    real(real32) :: height
    !! Height of the slab
    logical :: lcycle
    !! Boolean whether to cycle through the slab
    type(basis_type) :: tmp_bas1,tmp_bas2
    !! Temporary basis structures
    type(confine_type) :: confine
    !! Confine structure along the specified axis
    type(term_arr_type) :: term
    !! List of terminations
    real(real32), dimension(3,3) :: tfmat
    !! Transformation matrix
    logical :: orthogonalise_
    !! Boolean whether to orthogonalise the lattice
    logical :: normalise_
    !! Boolean whether to normalise the lattice
    logical :: break_on_fail_
    !! Boolean whether to break on failure


    character(len=256) :: warn_msg

    integer, allocatable, dimension(:,:,:) :: bas_map,t1bas_map
    real(real32), allocatable, dimension(:,:) :: trans


    orthogonalise_ = .true.
    if(present(orthogonalise)) orthogonalise_ = orthogonalise
    break_on_fail_ = .true.
    if(present(break_on_fail)) break_on_fail_ = break_on_fail
    normalise_ = .true.
    if(present(normalise)) normalise_ = normalise
    surface_ = 0
    if(present(surface))then
       select case(size(surface,dim=1))
       case(1)
          surface_(:) = surface(1)
       case(2)
          surface_ = surface
       case default
          write(0,'(1X,"ERROR: Internal error in generate_terminations")')
          write(0,'(2X,"The surface termination indices are not of the correct size")')
          return
       end select
    end if

    !! copy lattice and basis for manipulating
    call tmp_bas1%copy(basis)
    allocate(bas_map(tmp_bas1%nspec,maxval(tmp_bas1%spec(:)%num,dim=1),2))
    bas_map = -1


    write(6,'(1X,"Using supplied plane...")')
    tfmat = planecutter(tmp_bas1%lat,real(miller_plane,real32))
    call transformer(tmp_bas1,tfmat,bas_map)
    !call err_abort_print_struc(bas,"check.vasp","stop")


    !---------------------------------------------------------------------------
    ! Finds smallest thickness of the slab and increases to ...
    ! ... user-defined thickness
    !---------------------------------------------------------------------------
    confine%l = .false.
    confine%axis = this%axis
    confine%laxis = .false.
    confine%laxis(this%axis) = .true.
    if(allocated(trans)) deallocate(trans)
    allocate(trans(minval(tmp_bas1%spec(:)%num+2),3))
    call gldfnd(confine, tmp_bas1, tmp_bas1, trans, ntrans)
    tfmat(:,:) = 0._real32
    tfmat(1,1) = 1._real32
    tfmat(2,2) = 1._real32
    if(ntrans.eq.0)then
       tfmat(3,3)=1._real32
    else
       itmp1=minloc(abs(trans(:ntrans,this%axis)),dim=1,&
            mask=abs(trans(:ntrans,this%axis)).gt.1.D-3/modu(tmp_bas1%lat(this%axis,:)))
       tfmat(3,:)=trans(itmp1,:)
    end if
    if(all(abs(tfmat(3,:)).lt.1.E-5_real32)) tfmat(3,3) = 1._real32
    call transformer(tmp_bas1,tfmat,bas_map)
    if(.not.compare_stoichiometry(tmp_bas1,basis))then
       write(0,'(1X,"ERROR: Internal error in generate_terminations")')
       write(0,'(2X,"The gldfnd subroutine could not reproduce a valid primitive cell for the material")')
       if(ierror.eq.1)then
          call err_abort_print_struc(tmp_bas1, "broken_primitive.vasp", &
           "Code exiting due to IPRINT = 1")
       end if
       write(0,'(2X,"Skipping this lattice match")')
       return
    end if

    ! get the terminations
    term = get_termination_info( &
         tmp_bas1, this%axis, &
         lprint = .true., layer_sep = this%layer_separation_cutoff, &
         break_on_fail = break_on_fail_ &
    )
    if(term%nterm .eq. 0)then
       write(warn_msg, '(A,I0,1X,I0,1X,I0,A)') &
            "No terminations found for Miller plane (",miller_plane,")"
       call print_warning(trim(warn_msg))
       return
    end if

    ! set thickness if provided by user
    if(present(num_layers))then
       num_layers_ = num_layers
    else
       num_layers_ = 1
    end if

    ! determine tolerance for layer separations (termination tolerance)
    ! ... this is different from layer_sep
    call set_layer_tol(term)

    ! determine required extension and perform that
    call set_slab_height(tmp_bas1,bas_map,term,surface_,&
         height,num_layers_, thickness, ncells,&
         term_start,term_end,iterm_step &
    )
    
    !---------------------------------------------------------------------------
    ! Normalise lattice
    !---------------------------------------------------------------------------
    if(normalise_)then
       call reducer(tmp_bas1)
       tmp_bas1%lat = MATNORM(tmp_bas1%lat)
    end if
    

    !---------------------------------------------------------------------------
    ! loop over terminations and write them
    !---------------------------------------------------------------------------
    num_structures = ( term_end - term_start ) / iterm_step + 1
    allocate(output(num_structures))
    do iterm = term_start, term_end, iterm_step
       i = ( iterm - term_start ) / iterm_step + 1 
       call output(i)%copy(tmp_bas1)
       if(allocated(t1bas_map)) deallocate(t1bas_map)
       allocate(t1bas_map,source=bas_map)
       call build_slab(output(i),bas_map,term,[iterm,surface_(2)],&
            thickness, ncells, num_layers_, height,&
            "lw", lcycle, orthogonalise_, this%vacuum_gap &
       )
    end do
    if(.not.allocated(this%structures))then
       call move_alloc(output,this%structures)
    else
       this%structures = [ this%structures, output ]
    end if

   end subroutine generate_terminations
!###############################################################################

end module artemis__termination_generator
