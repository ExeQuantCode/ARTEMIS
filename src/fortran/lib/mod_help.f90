module mod_help
  use artemis__io_utils, only: err_abort, tag_type, io_print_help
  implicit none


  private

  public :: settings_help
  public :: cell_edits_help
  public :: interface_help


  !  logical, save :: ltag_present(ntags)  !!!REPLACE READVAR WITH THIS

  ! Settings number of tags
  integer, parameter :: ntags_settings=11
  ! Settings tags
  integer, parameter :: itask_tag=1
  integer, parameter :: irestart_tag=2
  integer, parameter :: is1file_tag=3
  integer, parameter :: is2file_tag=4
  integer, parameter :: iprint_tag=5
  integer, parameter :: iclock_tag=6
  integer, parameter :: imdir_tag=7
  integer, parameter :: isdir_tag=8
  integer, parameter :: iifmt_tag=9
  integer, parameter :: iofmt_tag=10
  integer, parameter :: itol_sym_tag=11


  ! Cell_edits number of tags
  integer, parameter :: ntags_cell_edits=15
  ! Cell_edits tags
  integer, parameter :: iout_file_tag=1
  integer, parameter :: ilsurf_gen_CE_tag=2
  integer, parameter :: imiller_tag=3
  integer, parameter :: inum_layers_tag=4
  integer, parameter :: ishift_tag=5
  integer, parameter :: ishift_region_tag=6
  integer, parameter :: ivacuum_tag=7
  integer, parameter :: itfmat_tag=8
  integer, parameter :: ilayer_sep_CE_tag=9
  integer, parameter :: ilortho_CE_tag=10
  integer, parameter :: isurf_tag=11
  integer, parameter :: ilnorm_lat_tag=12
  integer, parameter :: imin_thick_tag=13
  integer, parameter :: iuse_pricel_tag=14
  integer, parameter :: irequire_stoich_tag=15

  integer, parameter :: ntags_depr_cell_edits=1
  ! Cell_edits deprecated tags
  integer, parameter :: islab_thick_tag=1


  ! Interface number of tags
  integer, parameter :: ntags_interface=59
  ! Interface tags
  integer, parameter :: inintf_tag=1
  integer, parameter :: iimatch_tag=2
  integer, parameter :: inmatch_tag=3
  integer, parameter :: igen_intfs_tag=4
  integer, parameter :: iaxis_tag=5
  integer, parameter :: ilw_miller_tag=6
  integer, parameter :: iup_miller_tag=7
  integer, parameter :: ilw_num_layers_tag=8
  integer, parameter :: iup_num_layers_tag=9
  integer, parameter :: ishiftdir_tag=10
  integer, parameter :: iishift_tag=11
  integer, parameter :: inshift_tag=12
  integer, parameter :: ic_scale_tag=13
  integer, parameter :: ishift_udef_tag=14
  integer, parameter :: iidepth_tag=15
  integer, parameter :: intf_depth_tag=16
  integer, parameter :: iswapdir_tag=17
  integer, parameter :: iiswap_tag=18
  integer, parameter :: iswap_den_tag=19
  integer, parameter :: inswap_tag=20
  integer, parameter :: inmiller_tag=21
  integer, parameter :: ilprint_matches_tag=22
  integer, parameter :: ilprint_terms_tag=23
  integer, parameter :: ilprint_shifts_tag=24
  integer, parameter :: ilw_layered_tag=25
  integer, parameter :: iup_layered_tag=26
  integer, parameter :: ilw_surf_tag=27
  integer, parameter :: iup_surf_tag=28
  integer, parameter :: itol_maxlen_tag=29
  integer, parameter :: itol_maxarea_tag=30
  integer, parameter :: itol_vec_tag=31
  integer, parameter :: itol_ang_tag=32
  integer, parameter :: itol_area_tag=33
  integer, parameter :: itol_maxfind_tag=34
  integer, parameter :: itol_maxsize_tag=35
  integer, parameter :: ilreduce_tag=36
  integer, parameter :: iicheck_tag=37
  integer, parameter :: ilsurf_gen_tag=38
  integer, parameter :: iiintf_tag=39
  integer, parameter :: ilayer_sep_tag=40
  integer, parameter :: ilw_layer_sep_tag=41
  integer, parameter :: iup_layer_sep_tag=42
  integer, parameter :: interm_tag=43
  integer, parameter :: imbond_maxlen_tag=44
  integer, parameter :: iswap_sigma_tag=45
  integer, parameter :: iswap_depth_tag=46
  integer, parameter :: iintf_loc_tag=47
  integer, parameter :: ilmirror_tag=48
  integer, parameter :: ilortho_tag=49
  integer, parameter :: ilw_use_pricel_tag=50
  integer, parameter :: iup_use_pricel_tag=51
  integer, parameter :: ilw_bulk_modulus_tag=52
  integer, parameter :: iup_bulk_modulus_tag=53
  integer, parameter :: ilc_fix_tag=54
  integer, parameter :: ilbreak_on_no_term_tag=55
  integer, parameter :: ilw_min_thick_tag=56
  integer, parameter :: iup_min_thick_tag=57
  integer, parameter :: ilw_require_stoich_tag=58
  integer, parameter :: iup_require_stoich_tag=59

  integer, parameter :: ntags_depr_interface=2
  ! Cell_edits deprecated tags
  integer, parameter :: ilw_slab_thick_tag=1
  integer, parameter :: iup_slab_thick_tag=2



contains

!!!#############################################################################
!!! setup settings tag descriptions
!!!#############################################################################
  function setup_settings_tags() result(tag)
    implicit none
    type(tag_type), dimension(ntags_settings) :: tag


    tag(itask_tag)%name    = 'TASK'
    tag(itask_tag)%type    = 'I'
    tag(itask_tag)%summary = 'Job type'
    tag(itask_tag)%allowed = '0, 1, 2'
    tag(itask_tag)%default = '1'
    tag(itask_tag)%description = &
         'Defines the type of job to be run\n&
         &ARTEMIS will ignore any cards not consistent with job type\n&
         & 0 = CELL_EDITS\n&
         & 1 = INTERFACES\n&
         & 2 = DEFECTS     # NOT YET IMPLEMENTED\n'

    tag(irestart_tag)%name    = 'RESTART'
    tag(irestart_tag)%type    = 'I'
    tag(irestart_tag)%summary = 'Whether job is a restart'
    tag(irestart_tag)%allowed = '0, 1'
    tag(irestart_tag)%default = '0'
    tag(irestart_tag)%description = &
         'Specifies ARTEMIS to run as restart job\n&
         &Use STRUC1_FILE as pregenerated interface structure'

    tag(is1file_tag)%name    = 'STRUC1_FILE'
    tag(is1file_tag)%type    = 'S'
    tag(is1file_tag)%summary = 'Input structure 1 filename'
    tag(is1file_tag)%allowed = 'Up to 200 characters'
    tag(is1file_tag)%default = 'POSCAR'
    tag(is1file_tag)%description = &
         'Name of the input structure file.\n&
         &If generating an interface, it is be taken as the lower crystal'

    tag(is2file_tag)%name    = 'STRUC2_FILE'
    tag(is2file_tag)%type    = 'S'
    tag(is2file_tag)%summary = 'Input structure 2 filename'
    tag(is2file_tag)%allowed = 'Up to 200 characters'
    tag(is2file_tag)%default = '(empty)'
    tag(is2file_tag)%description = &
         'Name of the input structure file.\nIf generating an interface, &
         &it is be taken as the upper crystal'

    tag(iprint_tag)%name    = 'IPRINT'
    tag(iprint_tag)%type    = 'I'
    tag(iprint_tag)%summary = 'Verbosity of output printing'
    tag(iprint_tag)%allowed = '-1, 0, 1, 2, 3'
    tag(iprint_tag)%default = '-1'
    tag(iprint_tag)%description = &
         'Level of printing verbosity when running the code'

    tag(iclock_tag)%name    = 'CLOCK'
    tag(iclock_tag)%type    = 'I'
    tag(iclock_tag)%summary = 'Random seed clock number'
    tag(iclock_tag)%allowed = 'Any integer number'
    tag(iclock_tag)%default = '(filled by system clock)'
    tag(iclock_tag)%description = &
         'Random number generator seed value.\n&
         &Allows the user to regenerate the same structures that were made using random numbers'

    tag(imdir_tag)%name    = 'MASTER_DIR'
    tag(imdir_tag)%type    = 'S'
    tag(imdir_tag)%summary = 'Output structure directory'
    tag(imdir_tag)%allowed = 'Up to 100 characters'
    tag(imdir_tag)%default = 'DINTERFACES'
    tag(imdir_tag)%description = &
         'Directory in which to write all generated structures.\n&
         &If not present, the directory will first be made.\n&
         &This directory will be populated with other directories'

    tag(isdir_tag)%name    = 'SUBDIR_PREFIX'
    tag(isdir_tag)%type    = 'S'
    tag(isdir_tag)%summary = 'Subdirectory prefix'
    tag(isdir_tag)%allowed = 'Up to 100 characters'
    tag(isdir_tag)%default = 'D'
    tag(isdir_tag)%description = &
         'Prefix for the subdirectories containing structure files.\n&
         &Not necessary to have'

    tag(iifmt_tag)%name    = 'INPUT_FMT'
    tag(iifmt_tag)%type    = 'S'
    tag(iifmt_tag)%summary = 'Input structure data format'
    tag(iifmt_tag)%allowed = 'VASP, CASTEP, QuantumEspresso (QE), CRYSTAL'
    tag(iifmt_tag)%default = 'VASP'
    tag(iifmt_tag)%description = &
         'Specifies the format to use to read in the structure data'

    tag(iofmt_tag)%name    = 'OUTPUT_FMT'
    tag(iofmt_tag)%type    = 'S'
    tag(iofmt_tag)%summary = 'Output structure data format'
    tag(iofmt_tag)%allowed = 'VASP, CASTEP, QuantumEspresso (QE), CRYSTAL'
    tag(iofmt_tag)%default = 'VASP'
    tag(iofmt_tag)%description = &
         'Specifies the format in which to print the generated structure data'

    tag(itol_sym_tag)%name    = 'TOL_SYM'
    tag(itol_sym_tag)%type    = 'R'
    tag(itol_sym_tag)%summary = 'Symmetry tolerance/precision'
    tag(itol_sym_tag)%allowed = 'Any positive real number'
    tag(itol_sym_tag)%default = '1.D-6'
    tag(itol_sym_tag)%description = &
         'Sets the precision to which symmetries are determined'


  end function setup_settings_tags
!!!#############################################################################


!!!#############################################################################
!!! setup settings tag descriptions
!!!#############################################################################
  function setup_cell_edits_tags() result(tag)
    implicit none
    type(tag_type), dimension(ntags_cell_edits) :: tag

    
    tag(ishift_tag)%name    = 'SHIFT'
    tag(ishift_tag)%type    = 'B'
    tag(ishift_tag)%summary = 'Shift the basis in a cell'
    tag(ishift_tag)%allowed = ''
    tag(ishift_tag)%default = ''
    tag(ishift_tag)%description = &
         'NEED TO SET UP SHIFT REGION (ISSUE WITH ITS TWO AXES!!!&
         &Example:\n&
         &  SHIFT\n&
         &    axis = 3\n&
         &    val = 0.5\n&
         &  ENDSHIFT'

    tag(ivacuum_tag)%name    = 'VACUUM'
    tag(ivacuum_tag)%type    = 'B'
    tag(ivacuum_tag)%summary = 'Add/remove vacuum at a point in a cell'
    tag(ivacuum_tag)%allowed = ''
    tag(ivacuum_tag)%default = '14'
    tag(ivacuum_tag)%description = &
         'Vacuum gap (only applied for surface generation or when &
         &user-defined).\n&
         &Example:\n&
         &  VACUUM = 14 !(Å)\n&
         &  VACUUM\n&
         &    axis = 3\n&
         &    loc = 0.8\n&
         &    val = 0.5\n&
         &  ENDVACUUM'

    tag(itfmat_tag)%name    = 'TFMAT'
    tag(itfmat_tag)%type    = 'B'
    tag(itfmat_tag)%summary = 'Apply a transformation matrix'
    tag(itfmat_tag)%allowed = ''
    tag(itfmat_tag)%default = ''
    tag(itfmat_tag)%description = &
         'Example:\n&
         &  TFMAT\n&
         &    2.0  1.0  0 0\n&
         &    0.0  1.0  0.0\n&
         &    0.0  0.0  1.0\n&
         &  ENDTFMAT'

    tag(ishift_region_tag)%name    = 'SHIFT_REGION'
    tag(ishift_region_tag)%type    = 'B'
    tag(ishift_region_tag)%summary = 'Shift a region of the basis in a cell'
    tag(ishift_region_tag)%allowed = ''
    tag(ishift_region_tag)%default = ''
    tag(ishift_region_tag)%description = &
         'ISSUE AXIS ALONG WHICH BOUNDS LIE!\n&
         &Example:\n&
         &  SHIFT\n&
         &    axis = 1\n&
         &    val = 0.7\n&
         &    bounds = 0.3 0.5\n&
         &  ENDSHIFT'

    tag(iout_file_tag)%name    = 'OUTPUT_FILE'
    tag(iout_file_tag)%type    = 'S'
    tag(iout_file_tag)%summary = 'Output structure filename'
    tag(iout_file_tag)%allowed = 'Up to 200 characters'
    tag(iout_file_tag)%default = 'POSCAR, struc.cell, struc.geom'
    tag(iout_file_tag)%description = &
         'Name of the output structure file.\n&
         &Default output filename depends on the "OUTPUT_FMT" tag.\n&
         &If the input file has the same name as the default output filename, &
         &then "_out" will be appended to the end of the output filename.'

    tag(ilsurf_gen_CE_tag)%name    = 'LSURF_GEN'
    tag(ilsurf_gen_CE_tag)%type    = 'L'
    tag(ilsurf_gen_CE_tag)%summary = 'Generate surface information'
    tag(ilsurf_gen_CE_tag)%allowed = 'TRUE or FALSE'
    tag(ilsurf_gen_CE_tag)%default = 'FALSE'
    tag(ilsurf_gen_CE_tag)%description = &
         'Prints the surface terminations of a Miller plane into DTERMINATIONS &
         &directory'

    tag(inum_layers_tag)%name    = 'NUM_LAYERS'
    tag(inum_layers_tag)%type    = 'I'
    tag(inum_layers_tag)%summary = 'Number of layers of crystal'
    tag(inum_layers_tag)%allowed = 'Any positive integer number'
    tag(inum_layers_tag)%default = '(empty)'
    tag(inum_layers_tag)%description = &
         'Defines the number of primitive layers to use for the slab'

    tag(imin_thick_tag)%name    = 'MIN_THICKNESS'
    tag(imin_thick_tag)%type    = 'R'
    tag(imin_thick_tag)%summary = 'Minimum thickness of slab'
    tag(imin_thick_tag)%allowed = 'Any positive real number'
    tag(imin_thick_tag)%default = '10.0'
    tag(imin_thick_tag)%description = &
         'Defines the minimum thickness of the lower crystal (in Å).\n&
         &The generated slab will be the smallest possible thickness equal to &
         &or greater than this value.'

    tag(iuse_pricel_tag)%name    = 'USE_PRICEL'
    tag(iuse_pricel_tag)%type    = 'L'
    tag(iuse_pricel_tag)%summary = 'Use primitive cell'
    tag(iuse_pricel_tag)%allowed = 'TRUE or FALSE'
    tag(iuse_pricel_tag)%default = 'TRUE'
    tag(iuse_pricel_tag)%description = &
         'Defines whether to generate and use the primitive unit cell &
         &for the crystal'

    tag(irequire_stoich_tag)%name    = 'REQUIRE_STOICH'
    tag(irequire_stoich_tag)%type    = 'L'
    tag(irequire_stoich_tag)%summary = 'Maintain stoichiometry for terminations'
    tag(irequire_stoich_tag)%allowed = 'TRUE or FALSE'
    tag(irequire_stoich_tag)%default = 'FALSE'
    tag(irequire_stoich_tag)%description = &
         'Defines whether to maintain stoichiometry for the terminations.\n&
         &If TRUE, ARTEMIS will only generate terminations that are consistent &
         &with the stoichiometry of the bulk crystal.\n&
         &If FALSE, ARTEMIS will generate all possible terminations.'

    tag(imiller_tag)%name    = 'MILLER_PLANE'
    tag(imiller_tag)%type    = 'U'
    tag(imiller_tag)%summary = 'Crystal Miller plane'
    tag(imiller_tag)%allowed = 'Three integer numbers'
    tag(imiller_tag)%default = '(empty)'
    tag(imiller_tag)%description = &
         'Generate surfaces using the defined Miller plane'

    tag(ilortho_CE_tag)%name    = 'LORTHO'
    tag(ilortho_CE_tag)%type    = 'L'
    tag(ilortho_CE_tag)%summary = 'Surface axis perpendicular to surface'
    tag(ilortho_CE_tag)%allowed = 'TRUE or FALSE'
    tag(ilortho_CE_tag)%default = 'TRUE'
    tag(ilortho_CE_tag)%description = &
         'Defines whether to generate surfaces with the surface axis &
         &perpendicular to the surface'
    
    tag(ilayer_sep_CE_tag)%name    = 'LAYER_SEP'
    tag(ilayer_sep_CE_tag)%type    = 'R'
    tag(ilayer_sep_CE_tag)%summary = 'Min size of gap between layers'
    tag(ilayer_sep_CE_tag)%allowed = 'Any number greater than or equal to zero'
    tag(ilayer_sep_CE_tag)%default = '1.0'
    tag(ilayer_sep_CE_tag)%description = &
         'Defines the minimum size of gaps along the Miller direction that &
         &distinguish between separate layers (in Å)'

    tag(isurf_tag)%name    = 'SURFACE'
    tag(isurf_tag)%type    = 'U'
    tag(isurf_tag)%summary = 'Crystal surface terminations'
    tag(isurf_tag)%allowed = 'One or two integer numbers (<INT> [INT])'
    tag(isurf_tag)%default = '(empty)'
    tag(isurf_tag)%description = &
         'Defines the bottom and top surface terminations of the crystal'

    tag(ilnorm_lat_tag)%name    = 'LNORM_LAT'
    tag(ilnorm_lat_tag)%type    = 'L'
    tag(ilnorm_lat_tag)%summary = 'Apply Buerger cell renormalisation'
    tag(ilnorm_lat_tag)%allowed = 'TRUE or FALSE'
    tag(ilnorm_lat_tag)%default = 'TRUE'
    tag(ilnorm_lat_tag)%description = &
         'Defines whether the lattice (cell) is printed out as-is, or whether &
         &the Buerger reduction is performed to make lattice more manageable'


!!! SET UP A CUTTER FUNCTION


  end function setup_cell_edits_tags
!!!#############################################################################


!!!#############################################################################
!!! setup interace  tag descriptions
!!!#############################################################################
  function setup_interface_tags() result(tag)
    implicit none
    type(tag_type), dimension(ntags_interface) :: tag


    tag(ilprint_matches_tag)%name    = 'LPRINT_MATCHES'
    tag(ilprint_matches_tag)%type    = 'L'
    tag(ilprint_matches_tag)%summary = 'Print lattice match information'
    tag(ilprint_matches_tag)%allowed = 'TRUE or FALSE'
    tag(ilprint_matches_tag)%default = 'FALSE'
    tag(ilprint_matches_tag)%description = &
         'Determines whether to print information of the identified lattice matches'

    tag(igen_intfs_tag)%name    = 'LGEN_INTERFACES'
    tag(igen_intfs_tag)%type    = 'L'
    tag(igen_intfs_tag)%summary = 'Generate interfaces/just print matches'
    tag(igen_intfs_tag)%allowed = 'TRUE or FALSE'
    tag(igen_intfs_tag)%default = 'TRUE'
    tag(igen_intfs_tag)%description = &
         &'Defines whether transformation matrices of lattice matches are just &
         &printed, or also used to generate interface structures.&
         &\nTRUE generates interfaces'

    tag(inintf_tag)%name    = 'NINTF'
    tag(inintf_tag)%type    = 'I'
    tag(inintf_tag)%summary = 'Max number of interfaces generated'
    tag(inintf_tag)%allowed = 'Any positive integer number'
    tag(inintf_tag)%default = '100'
    tag(inintf_tag)%description = &
         &'Defines the maximum number of interface structures to be generated (excludes swaps)\n&
         & e.g. ARTEMIS generates min(NINTF,NMATCH*NSHIFT) structure files.'

    tag(iiintf_tag)%name    = 'IINTF'
    tag(iiintf_tag)%type    = 'I'
    tag(iiintf_tag)%summary = 'Specify which interface to generate'
    tag(iiintf_tag)%allowed = 'Any positive integer number <= NINTF'
    tag(iiintf_tag)%default = '(empty)'
    tag(iiintf_tag)%description = &
         &'Specifies that only 1 unique interface is to be generated and which one'

    tag(iimatch_tag)%name    = 'IMATCH'
    tag(iimatch_tag)%type    = 'I'
    tag(iimatch_tag)%summary = 'Lattice matching method'
    tag(iimatch_tag)%allowed = '0, 1, 2'
    tag(iimatch_tag)%default = '0'
    tag(iimatch_tag)%description = &
         &'Defines the method used to match the two crystals\n \n&
         & 0 = Full lattice matching, but can be constrained to specific Miller planes&
         & for lower and upper crystals via LW_MILLER and UP_MILLER\n&
         & 1 = Cycle over all transformation matrices for lower crystal (within tolerances)&
         & and find the closest transformation matrix to map the upper crystal to the lower\n&
         & 2 = Cycle over all transformation matrices for lower and upper crystals (within tolerances)&
         & and find the closest matches between the two crystals\n \n&
         & NOTE: It is always recommended to stick to IMATCH=0 for best results.\n&
         & Other methods are experimental and may not work as expected.\n&
         & Only method 0 works with LW_MILLER and UP_MILLER.'

    tag(inmatch_tag)%name    = 'NMATCH'
    tag(inmatch_tag)%type    = 'I'
    tag(inmatch_tag)%summary = 'Max number of unique lattice matches'
    tag(inmatch_tag)%allowed = 'Any positive integer number'
    tag(inmatch_tag)%default = '5'
    tag(inmatch_tag)%description = &
         'Defines the maximum number of unique lattice matches, within tolerances.\n&
         &Identifies all matches within tolerances and outputs the NMATCH best matched'

    tag(iaxis_tag)%name    = 'AXIS'
    tag(iaxis_tag)%type    = 'I'
    tag(iaxis_tag)%summary = 'Axis along which to output interface'
    tag(iaxis_tag)%allowed = '1, 2, 3'
    tag(iaxis_tag)%default = '3'
    tag(iaxis_tag)%description = &
         'NOT YET FULLY IMPLEMENTED! Defines the axis along which to print &
         &interfaces along.\n&
         &NOTE: this does not change the interfaces generated, simply whether &
         &a generated interface will lie along a, b or c in the generated &
         &output structure file'

    tag(iintf_loc_tag)%name    = 'INTF_LOC'
    tag(iintf_loc_tag)%type    = 'V'
    tag(iintf_loc_tag)%summary = 'RESTART! Interface location (direct)'
    tag(iintf_loc_tag)%allowed = 'Two real numbers between 0 and 1'
    tag(iintf_loc_tag)%default = '(empty)'
    tag(iintf_loc_tag)%description = &
         'For a restart job, user can define the interface location &
         &if they prefer ARTEMIS not identify it itself (direct coords)\n&
         &(useful for large interfaces where ARTEMIS has already been used &
         &to determine the interface location before)'

    tag(ilprint_terms_tag)%name    = 'LPRINT_TERMS'
    tag(ilprint_terms_tag)%type    = 'L'
    tag(ilprint_terms_tag)%summary = 'Print termination information'
    tag(ilprint_terms_tag)%allowed = 'TRUE or FALSE'
    tag(ilprint_terms_tag)%default = 'FALSE'
    tag(ilprint_terms_tag)%description = &
         'Determines whether to print information of the identified terminations'

    tag(interm_tag)%name    = 'NTERM'
    tag(interm_tag)%type    = 'I'
    tag(interm_tag)%summary = 'Max terminations per Miller plane'
    tag(interm_tag)%allowed = 'Any positive integer number'
    tag(interm_tag)%default = '5'
    tag(interm_tag)%description = &
         'Defines the maximum number of possible unique terminations to consider per material per Miller plane'

    tag(ilw_layered_tag)%name    = 'LW_LAYERED'
    tag(ilw_layered_tag)%type    = 'L'
    tag(ilw_layered_tag)%summary = 'Is lower crystal layered'
    tag(ilw_layered_tag)%allowed = 'TRUE or FALSE'
    tag(ilw_layered_tag)%default = 'FALSE'
    tag(ilw_layered_tag)%description = &
         'Defines whether the lower crystal is a layered material'

    tag(iup_layered_tag)%name    = 'UP_LAYERED'
    tag(iup_layered_tag)%type    = 'L'
    tag(iup_layered_tag)%summary = 'Is upper crystal layered'
    tag(iup_layered_tag)%allowed = 'TRUE or FALSE'
    tag(iup_layered_tag)%default = 'FALSE'
    tag(iup_layered_tag)%description = &
         'Defines whether the upper crystal is a layered material'

    tag(ilw_miller_tag)%name    = 'LW_MILLER'
    tag(ilw_miller_tag)%type    = 'U'
    tag(ilw_miller_tag)%summary = 'Lower crystal Miller plane'
    tag(ilw_miller_tag)%allowed = 'Three integer numbers'
    tag(ilw_miller_tag)%default = '(empty)'
    tag(ilw_miller_tag)%description = &
         &'Confines the lower crystal to this Miller plane for lattice matching.\n\n&
         &NOTE: Can only be used with IMATCH=0.\n\n&
         &NOTE: Miller indices used in ARTEMIS are defined for the cell in &
         &use. Experimental Miller indices are presented with respect to the &
         & primitive cell. To use proper Miller indices, ensure LW_USE_PRICEL=T.'

    tag(iup_miller_tag)%name    = 'UP_MILLER'
    tag(iup_miller_tag)%type    = 'U'
    tag(iup_miller_tag)%summary = 'Upper crystal Miller plane'
    tag(iup_miller_tag)%allowed = 'Three integer numbers'
    tag(iup_miller_tag)%default = '(empty)'
    tag(iup_miller_tag)%description = &
         &'Confines the upper crystal to this Miller plane for lattice matching.\n\n&
         &NOTE: Can only be used with IMATCH=0.\n\n&
         &NOTE: Miller indices used in ARTEMIS are defined for the cell in &
         &use. Experimental Miller indices are presented with respect to the &
         & primitive cell. To use proper Miller indices, ensure UP_USE_PRICEL=T.'

    tag(inmiller_tag)%name    = 'NMILLER'
    tag(inmiller_tag)%type    = 'I'
    tag(inmiller_tag)%summary = 'Number of Miller planes'
    tag(inmiller_tag)%allowed = 'Any positive integer number'
    tag(inmiller_tag)%default = '10'
    tag(inmiller_tag)%description = &
         'Defines the number of Miller planes to search over for each crystal.'

    tag(ilw_num_layers_tag)%name    = 'LW_NUM_LAYERS'
    tag(ilw_num_layers_tag)%type    = 'I'
    tag(ilw_num_layers_tag)%summary = 'Number of layers of lower crystal'
    tag(ilw_num_layers_tag)%allowed = 'Any positive integer number'
    tag(ilw_num_layers_tag)%default = '(empty)'
    tag(ilw_num_layers_tag)%description = &
         'Defines the number of primitive layers to use for the lower crystal'

    tag(iup_num_layers_tag)%name    = 'UP_NUM_LAYERS'
    tag(iup_num_layers_tag)%type    = 'I'
    tag(iup_num_layers_tag)%summary = 'Number of layers of upper crystal'
    tag(iup_num_layers_tag)%allowed = 'Any positive integer number'
    tag(iup_num_layers_tag)%default = '(empty)'
    tag(iup_num_layers_tag)%description = &
         'Defines the number of primitive layers to use for the upper crystal'

    tag(ilw_min_thick_tag)%name    = 'LW_MIN_THICKNESS'
    tag(ilw_min_thick_tag)%type    = 'R'
    tag(ilw_min_thick_tag)%summary = 'Minimum thickness of lower crystal'
    tag(ilw_min_thick_tag)%allowed = 'Any positive real number'
    tag(ilw_min_thick_tag)%default = '10.0'
    tag(ilw_min_thick_tag)%description = &
         'Defines the minimum thickness of the lower crystal (in Å).\n&
         &The generated slab will be the smallest possible thickness equal to &
         &or greater than this value.'

    tag(iup_min_thick_tag)%name    = 'UP_MIN_THICKNESS'
    tag(iup_min_thick_tag)%type    = 'R'
    tag(iup_min_thick_tag)%summary = 'Minimum thickness of upper crystal'
    tag(iup_min_thick_tag)%allowed = 'Any positive real number'
    tag(iup_min_thick_tag)%default = '10.0'
    tag(iup_min_thick_tag)%description = &
         'Defines the minimum thickness of the upper crystal (in Å).\n&
         &The generated slab will be the smallest possible thickness equal to &
         &or greater than this value.'

    tag(ilw_surf_tag)%name    = 'LW_SURFACE'
    tag(ilw_surf_tag)%type    = 'U'
    tag(ilw_surf_tag)%summary = 'Lower crystal surface terminations'
    tag(ilw_surf_tag)%allowed = 'One or two integer numbers (<INT> [INT])'
    tag(ilw_surf_tag)%default = '(empty)'
    tag(ilw_surf_tag)%description = &
         'Defines the bottom and top surface terminations of the lower crystal'

    tag(iup_surf_tag)%name    = 'UP_SURFACE'
    tag(iup_surf_tag)%type    = 'U'
    tag(iup_surf_tag)%summary = 'Upper crystal surface terminations'
    tag(iup_surf_tag)%allowed = 'One or two integer numbers (<INT> [INT])'
    tag(iup_surf_tag)%default = '(empty)'
    tag(iup_surf_tag)%description = &
         'Defines the bottom and top surface terminations of the upper crystal'

    tag(ilayer_sep_tag)%name    = 'LAYER_SEP'
    tag(ilayer_sep_tag)%type    = 'R'
    tag(ilayer_sep_tag)%summary = 'Min size of gap between layers'
    tag(ilayer_sep_tag)%allowed = 'Any number greater than or equal to zero'
    tag(ilayer_sep_tag)%default = '1.0'
    tag(ilayer_sep_tag)%description = &
         'Defines the minimum size of gaps along the Miller direction that &
         &distinguish between separate layers (in Å).\n&
         &This tag is ignored if LW_LAYER_SEP and UP_LAYER_SEP are defined.\n&
         &If either LW_LAYER_SEP and UP_LAYER_SEP are undefined, then this value is used.'

    tag(ilw_layer_sep_tag)%name    = 'LW_LAYER_SEP'
    tag(ilw_layer_sep_tag)%type    = 'R'
    tag(ilw_layer_sep_tag)%summary = 'Min size of gap between layers for lower &
         &structure'
    tag(ilw_layer_sep_tag)%allowed = 'Any number greater than or equal to zero'
    tag(ilw_layer_sep_tag)%default = '1.0'
    tag(ilw_layer_sep_tag)%description = &
         'Defines the minimum size of gaps along the Miller direction that &
         &distinguish between separate layers (in Å) for the lower structure'

    tag(iup_layer_sep_tag)%name    = 'UP_LAYER_SEP'
    tag(iup_layer_sep_tag)%type    = 'R'
    tag(iup_layer_sep_tag)%summary = 'Min size of gap between layers for upper &
         &structure'
    tag(iup_layer_sep_tag)%allowed = 'Any number greater than or equal to zero'
    tag(iup_layer_sep_tag)%default = '1.0'
    tag(iup_layer_sep_tag)%description = &
         'Defines the minimum size of gaps along the Miller direction that &
         &distinguish between separate layers (in Å) for the upper structure'

    tag(ilbreak_on_no_term_tag)%name = 'LBREAK_ON_NO_TERM'
    tag(ilbreak_on_no_term_tag)%type = 'L'
    tag(ilbreak_on_no_term_tag)%summary = 'Stop on no termination'
    tag(ilbreak_on_no_term_tag)%allowed = 'TRUE or FALSE'
    tag(ilbreak_on_no_term_tag)%default = 'TRUE'
    tag(ilbreak_on_no_term_tag)%description = &
         'Defines whether to stop the code if no terminations are found for a &
         &given Miller plane'

    tag(ilprint_shifts_tag)%name    = 'LPRINT_SHIFTS'
    tag(ilprint_shifts_tag)%type    = 'L'
    tag(ilprint_shifts_tag)%summary = 'Print shift information'
    tag(ilprint_shifts_tag)%allowed = 'TRUE or FALSE'
    tag(ilprint_shifts_tag)%default = 'FALSE'
    tag(ilprint_shifts_tag)%description = &
         'Determines whether to print information of the shifts generated'

    tag(ishiftdir_tag)%name    = 'SHIFTDIR'
    tag(ishiftdir_tag)%type    = 'S'
    tag(ishiftdir_tag)%summary = 'Shift directory name'
    tag(ishiftdir_tag)%allowed = 'Up to 100 characters'
    tag(ishiftdir_tag)%default = 'DSHIFT'
    tag(ishiftdir_tag)%description = &
         'Name of directory in which to generate shifts'

    tag(iishift_tag)%name    = 'ISHIFT'
    tag(iishift_tag)%type    = 'I'
    tag(iishift_tag)%summary = 'Shifting method'
    tag(iishift_tag)%allowed = '0, 1, 2, 3, 4'
    tag(iishift_tag)%default = '3'
    tag(iishift_tag)%description = &
         'Defines the method to generate unique offsets of lower crystal from the upper one.\n&
         & 0 = User defined shift value\n&
         & 1 = Random shifts\n&
         & 2 = Match bond to average of two crystals'' bulk bonds\n&
         & 3 = Descriptive shifts (best and worst) using best separation\n&
         & 4 = Density of neighbours matching'

    tag(inshift_tag)%name    = 'NSHIFT'
    tag(inshift_tag)%type    = 'I'
    tag(inshift_tag)%summary = 'Number of shifts per lattice match'
    tag(inshift_tag)%allowed = 'Any positive integer number'
    tag(inshift_tag)%default = '5'
    tag(inshift_tag)%description = &
         'Defines the number of unique shifts to generate per lattice match'

    tag(ic_scale_tag)%name    = 'C_SCALE'
    tag(ic_scale_tag)%type    = 'R'
    tag(ic_scale_tag)%summary = 'Interface separation scaling factor'
    tag(ic_scale_tag)%allowed = 'Any number greater than or equal to zero'
    tag(ic_scale_tag)%default = '1.5'
    tag(ic_scale_tag)%description = &
         'Defines the amount to scale the shift perpendicular to the interface'

    tag(imbond_maxlen_tag)%name    = 'MBOND_MAXLEN'
    tag(imbond_maxlen_tag)%type    = 'R'
    tag(imbond_maxlen_tag)%summary = 'Maximum considered missing bondlength'
    tag(imbond_maxlen_tag)%allowed = 'Any positive real number'
    tag(imbond_maxlen_tag)%default = '4.0 (Å)'
    tag(imbond_maxlen_tag)%description = &
         'ONLY USED IN ISHIFT = 4\n&
         & Defines the maximum length of a bond to consider as missing (Å)'

    tag(ishift_udef_tag)%name    = 'SHIFT'
    tag(ishift_udef_tag)%type    = 'B'
    tag(ishift_udef_tag)%summary = 'User-defined interface shifts'
    tag(ishift_udef_tag)%allowed = ''
    tag(ishift_udef_tag)%default = '(empty)'
    tag(ishift_udef_tag)%description = &
         'User defined shifts.&
         &\nA negative number in the 2nd example frees that axis to generated shifts.\n&
         &Example:\n&
         &  SHIFT = 2 !(Å)\n&
         &  SHIFT = 0.5 0.5 2.0 !(direct, direct, Å)\n&
         &  SHIFT\n&
         &    0.0 0.0 1.0 !(direct, direct, Å)\n&
         &    0.0 0.0 2.0 !(direct, direct, Å)\n&
         &    0.0 0.0 3.0 !(direct, direct, Å)\n&
         &  ENDSHIFT\n'

    tag(iidepth_tag)%name    = 'IDEPTH'
    tag(iidepth_tag)%type    = 'I'
    tag(iidepth_tag)%summary = 'Interface depth method'
    tag(iidepth_tag)%allowed = '0, 1'
    tag(iidepth_tag)%default = '0'
    tag(iidepth_tag)%description = &
         'Defines the method used to determine the number of atoms to use when considering the interface shifting'

    tag(intf_depth_tag)%name    = 'INTF_DEPTH'
    tag(intf_depth_tag)%type    = 'R'
    tag(intf_depth_tag)%summary = 'Definition of interface region'
    tag(intf_depth_tag)%allowed = 'Any number greater than or equal to zero'
    tag(intf_depth_tag)%default = '1.5'
    tag(intf_depth_tag)%description = &
         'The distance from the interface to consider atoms for the shifting.\n&
         &Only used for IDEPTH > 1'

    tag(iswapdir_tag)%name    = 'SWAPDIR'
    tag(iswapdir_tag)%type    = 'S'
    tag(iswapdir_tag)%summary = 'Swap directory name'
    tag(iswapdir_tag)%allowed = 'Up to 100 characters'
    tag(iswapdir_tag)%default = 'DSWAP'
    tag(iswapdir_tag)%description = &
         'Name of directory in which to generate swaps'

    tag(iiswap_tag)%name    = 'ISWAP'
    tag(iiswap_tag)%type    = 'I'
    tag(iiswap_tag)%summary = 'Swapping method'
    tag(iiswap_tag)%allowed = '0, 1'
    tag(iiswap_tag)%default = '0'
    tag(iiswap_tag)%description = &
         'Determines the method used to generate swaps.\n&
         & 0 = No swaps generated\n&
         & 1 = Random swaps'

    tag(iswap_den_tag)%name    = 'SWAP_DENSITY'
    tag(iswap_den_tag)%type    = 'R'
    tag(iswap_den_tag)%summary = 'Swap density across an interface'
    tag(iswap_den_tag)%allowed = 'Any positive real number'
    tag(iswap_den_tag)%default = '5'
    tag(iswap_den_tag)%description = &
         'Defines the number of swaps per unit area (atoms/Å²) to populate a structure with. e.g. the intermixing area density'

    tag(iswap_sigma_tag)%name    = 'SWAP_SIGMA'
    tag(iswap_sigma_tag)%type    = 'R'
    tag(iswap_sigma_tag)%summary = 'Sigma to define gaussian of interface'
    tag(iswap_sigma_tag)%allowed = 'Any positive real number'
    tag(iswap_sigma_tag)%default = '(empty)'
    tag(iswap_sigma_tag)%description = &
         'Defines a sigma that is used in a gaussian to define the &
         &reducing chance of atoms further from the interface being swapped'

    tag(iswap_depth_tag)%name    = 'SWAP_DEPTH'
    tag(iswap_depth_tag)%type    = 'R'
    tag(iswap_depth_tag)%summary = 'Max depth from the interface to swap'
    tag(iswap_depth_tag)%allowed = 'Any positive real number'
    tag(iswap_depth_tag)%default = '3.0'
    tag(iswap_depth_tag)%description = &
         'Defines the maximum distance from the interface to consider &
         &atoms for swaps'

    tag(ilmirror_tag)%name    = 'LMIRROR'
    tag(ilmirror_tag)%type    = 'L'
    tag(ilmirror_tag)%summary = 'Require mirror symmetry for swaps'
    tag(ilmirror_tag)%allowed = 'TRUE or FALSE'
    tag(ilmirror_tag)%default = 'TRUE'
    tag(ilmirror_tag)%description = &
         'Defines whether to maintain symmetry of the two interfaces in the structure by &
         &performing equivalent swaps on both interfaces.\n&
         &TRUE  = perform equivalent swaps on both surfaces\n&
         &FALSE = perform swaps on only one surface'

    tag(inswap_tag)%name    = 'NSWAP'
    tag(inswap_tag)%type    = 'I'
    tag(inswap_tag)%summary = 'Unique swap structures per shift'
    tag(inswap_tag)%allowed = 'Any positive integer number'
    tag(inswap_tag)%default = '5'
    tag(inswap_tag)%description = &
         'Defines the number of swap structures to generate per shift'

    tag(itol_maxlen_tag)%name    = 'MAXLEN'
    tag(itol_maxlen_tag)%type    = 'R'
    tag(itol_maxlen_tag)%summary = 'Max length of a lattice vector'
    tag(itol_maxlen_tag)%allowed = 'Any positive real number'
    tag(itol_maxlen_tag)%default = '20.0 (in Å)'
    tag(itol_maxlen_tag)%description = &
         'Defines the maximum possible length of a single lattice vector in a lattice match (in Å)'

    tag(itol_maxarea_tag)%name    = 'MAXAREA'
    tag(itol_maxarea_tag)%type    = 'R'
    tag(itol_maxarea_tag)%summary = 'Max area of a lattice match'
    tag(itol_maxarea_tag)%allowed = 'Any positive real number'
    tag(itol_maxarea_tag)%default = '400.0 (in Å)'
    tag(itol_maxarea_tag)%description = &
         'Defines the maximum possible area of a lattice match (in Å)'

    tag(itol_vec_tag)%name    = 'TOL_VEC'
    tag(itol_vec_tag)%type    = 'R'
    tag(itol_vec_tag)%summary = 'Vector tolerance for lattice match'
    tag(itol_vec_tag)%allowed = 'Any positive real number'
    tag(itol_vec_tag)%default = '5 (in °)'
    tag(itol_vec_tag)%description = &
         'The allowed tolerance of lattice vector fitting during lattice matching (in %)'

    tag(itol_ang_tag)%name    = 'TOL_ANG'
    tag(itol_ang_tag)%type    = 'R'
    tag(itol_ang_tag)%summary = 'Angle tolerance for lattice match'
    tag(itol_ang_tag)%allowed = 'Any positive real number'
    tag(itol_ang_tag)%default = '1 (in °)'
    tag(itol_ang_tag)%description = &
         'The allowed tolerance of angle fitting during lattice matching (in °)'

    tag(itol_area_tag)%name    = 'TOL_AREA'
    tag(itol_area_tag)%type    = 'R'
    tag(itol_area_tag)%summary = 'Area tolerance for lattice match'
    tag(itol_area_tag)%allowed = 'Any positive real number'
    tag(itol_area_tag)%default = '10 (in %)'
    tag(itol_area_tag)%description = &
         'The allowed tolerance of area fitting during lattice matching (in %)'

    tag(itol_maxfind_tag)%name    = 'TOL_MAXFIND'
    tag(itol_maxfind_tag)%type    = 'R'
    tag(itol_maxfind_tag)%summary = 'Max number of matches per unique area'
    tag(itol_maxfind_tag)%allowed = 'Any positive real number'
    tag(itol_maxfind_tag)%default = '100'
    tag(itol_maxfind_tag)%description = &
         'The maximum number of matches for each individual area'

    tag(itol_maxsize_tag)%name    = 'TOL_MAXSIZE'
    tag(itol_maxsize_tag)%type    = 'R'
    tag(itol_maxsize_tag)%summary = 'Max lattice vector extension'
    tag(itol_maxsize_tag)%allowed = 'Any positive real number'
    tag(itol_maxsize_tag)%default = '10'
    tag(itol_maxsize_tag)%description = &
         'The maximum possible lattice vector extension in lattice matching'

    tag(ilreduce_tag)%name    = 'LREDUCE'
    tag(ilreduce_tag)%type    = 'L'
    tag(ilreduce_tag)%summary = 'Reduce lattice matches, if possible'
    tag(ilreduce_tag)%allowed = 'TRUE or FALSE'
    tag(ilreduce_tag)%default = 'FALSE'
    tag(ilreduce_tag)%description = &
         'WARNING! CURRENTLY UNSTABLE AND MOSTLY DOES NOT WORK!\n&
         &If possible, reduce any lattice match to its smallest size.'

    tag(iicheck_tag)%name    = 'ICHECK'
    tag(iicheck_tag)%type    = 'I'
    tag(iicheck_tag)%summary = 'DEVELOPER: check on an interface'
    tag(iicheck_tag)%allowed = 'Any positive integer number'
    tag(iicheck_tag)%default = '(empty)'
    tag(iicheck_tag)%description = &
         'DEVELOPER TAG\nCheck on a specific interface to help bug check'

    tag(ilsurf_gen_tag)%name    = 'LSURF_GEN'
    tag(ilsurf_gen_tag)%type    = 'L'
    tag(ilsurf_gen_tag)%summary = 'Generate surface information'
    tag(ilsurf_gen_tag)%allowed = 'TRUE or FALSE'
    tag(ilsurf_gen_tag)%default = 'FALSE'
    tag(ilsurf_gen_tag)%description = &
         'Prints the surface terminations of a Miller plane in DTERMINATIONS &
         &directory.\n&
         &Prints surfaces for crystals that have had their Miller planes supplied using the "LW_MILLER" and "UP_MILLER" tags\n&
         &Inside DTERMINATIONS, populates directory DLW_TERMS with lower &
         &parent structure surfaces.\n&
         &Inside DTERMINATIONS, populates directory DUP_TERMS with upper &
         &parent structure surfaces.'

    tag(ilortho_tag)%name    = 'LORTHO'
    tag(ilortho_tag)%type    = 'L'
    tag(ilortho_tag)%summary = 'Surface axis perpendicular to surface'
    tag(ilortho_tag)%allowed = 'TRUE or FALSE'
    tag(ilortho_tag)%default = 'TRUE'
    tag(ilortho_tag)%description = &
         'Defines whether to generate surfaces with the surface axis &
         &perpendicular to the surface'

    tag(ilw_use_pricel_tag)%name    = 'LW_USE_PRICEL'
    tag(ilw_use_pricel_tag)%type    = 'L'
    tag(ilw_use_pricel_tag)%summary = 'Use lower primitive cell'
    tag(ilw_use_pricel_tag)%allowed = 'TRUE or FALSE'
    tag(ilw_use_pricel_tag)%default = 'TRUE'
    tag(ilw_use_pricel_tag)%description = &
         'Defines whether to generate and use the primitive unit cell &
         &for the lower crystal'

    tag(iup_use_pricel_tag)%name    = 'UP_USE_PRICEL'
    tag(iup_use_pricel_tag)%type    = 'L'
    tag(iup_use_pricel_tag)%summary = 'Use upper primitive cell'
    tag(iup_use_pricel_tag)%allowed = 'TRUE or FALSE'
    tag(iup_use_pricel_tag)%default = 'TRUE'
    tag(iup_use_pricel_tag)%description = &
         'Defines whether to generate and use the primitive unit cell &
         &for the upper crystal'

    tag(ilw_bulk_modulus_tag)%name    = 'LW_BULK_MODULUS'
    tag(ilw_bulk_modulus_tag)%type    = 'R'
    tag(ilw_bulk_modulus_tag)%summary = 'Bulk modulus of lower material'
    tag(ilw_bulk_modulus_tag)%allowed = 'Any positive real number'
    tag(ilw_bulk_modulus_tag)%default = '(empty)'
    tag(ilw_bulk_modulus_tag)%description = &
         'The bulk modulus of the upper material (units=GPa).\n&
         &If LW_ and UP_ not defined, bulk modulus is ignored and upper &
         &material takes all of the strain.\n&
         &NOTE: Units are not important, as long as LW_ and UP_ have the same units.'

    tag(iup_bulk_modulus_tag)%name    = 'UP_BULK_MODULUS'
    tag(iup_bulk_modulus_tag)%type    = 'R'
    tag(iup_bulk_modulus_tag)%summary = 'Bulk modulus of upper material'
    tag(iup_bulk_modulus_tag)%allowed = 'Any positive real number'
    tag(iup_bulk_modulus_tag)%default = '(empty)'
    tag(iup_bulk_modulus_tag)%description = &
         'The bulk modulus of the upper material (units=GPa).\n&
         &If LW_ and UP_ not defined, bulk modulus is ignored and upper &
         &material takes all of the strain.\n&
         &NOTE: Units are not important, as long as LW_ and UP_ have the same units.'

    tag(ilc_fix_tag)%name    = 'LC_FIX'
    tag(ilc_fix_tag)%type    = 'L'
    tag(ilc_fix_tag)%summary = 'Fix the c axis of each material'
    tag(ilc_fix_tag)%allowed = 'TRUE or FALSE'
    tag(ilc_fix_tag)%default = 'TRUE'
    tag(ilc_fix_tag)%description = &
         'The c axis (axis perpendicular to the interface plane) &
         & can be fixed (strained) or changed (unstrained) to compensate for interfacial strains.\n&
         &  TRUE  = fix the c axis\n&
         &  FALSE = extend/compress c axis to compensate for strain.'

    tag(ilw_require_stoich_tag)%name    = 'LW_REQUIRE_STOICH'
    tag(ilw_require_stoich_tag)%type    = 'L'
    tag(ilw_require_stoich_tag)%summary = 'Maintain stoichiometry for lower terminations'
    tag(ilw_require_stoich_tag)%allowed = 'TRUE or FALSE'
    tag(ilw_require_stoich_tag)%default = 'FALSE'
    tag(ilw_require_stoich_tag)%description = &
         'Defines whether to maintain stoichiometry for the terminations of the lower structure.\n&
         &If TRUE, ARTEMIS will only generate terminations that are consistent &
         &with the stoichiometry of the bulk crystal.\n&
         &If FALSE, ARTEMIS will generate all possible terminations.'

    tag(iup_require_stoich_tag)%name    = 'UP_REQUIRE_STOICH'
    tag(iup_require_stoich_tag)%type    = 'L'
    tag(iup_require_stoich_tag)%summary = 'Maintain stoichiometry for upper terminations'
    tag(iup_require_stoich_tag)%allowed = 'TRUE or FALSE'
    tag(iup_require_stoich_tag)%default = 'FALSE'
    tag(iup_require_stoich_tag)%description = &
         'Defines whether to maintain stoichiometry for the terminations of the upper structure.\n&
         &If TRUE, ARTEMIS will only generate terminations that are consistent &
         &with the stoichiometry of the bulk crystal.\n&
         &If FALSE, ARTEMIS will generate all possible terminations.'


  end function setup_interface_tags
!!!#############################################################################


!!!#############################################################################
!!! setup deprecated interface tag descriptions
!!!#############################################################################
  function setup_depr_cell_edits_tags() result(tag)
    implicit none
    type(tag_type), dimension(ntags_depr_cell_edits) :: tag

    tag(islab_thick_tag)%name    = 'SLAB_THICKNESS'
    tag(islab_thick_tag)%type    = 'I'
    tag(islab_thick_tag)%summary = 'Number of layers of crystal'
    tag(islab_thick_tag)%allowed = 'Any positive integer number'
    tag(islab_thick_tag)%default = '(empty)'
    tag(islab_thick_tag)%is_deprecated = .false.
    tag(islab_thick_tag)%to_be_deprecated = .true.
    tag(islab_thick_tag)%deprecated_version = '2.0.0'
    tag(islab_thick_tag)%deprecated_name = 'NUM_LAYERS'
    tag(islab_thick_tag)%description = &
         'Defines the number of primitive layers to use for the lower crystal'

  end function setup_depr_cell_edits_tags
!-------------------------------------------------------------------------------
  function setup_depr_interface_tags() result(tag)
    implicit none
    type(tag_type), dimension(ntags_depr_interface) :: tag

    tag(ilw_slab_thick_tag)%name    = 'LW_SLAB_THICKNESS'
    tag(ilw_slab_thick_tag)%type    = 'I'
    tag(ilw_slab_thick_tag)%summary = 'Number of layers of lower crystal'
    tag(ilw_slab_thick_tag)%allowed = 'Any positive integer number'
    tag(ilw_slab_thick_tag)%default = '(empty)'
    tag(ilw_slab_thick_tag)%is_deprecated = .false.
    tag(ilw_slab_thick_tag)%to_be_deprecated = .true.
    tag(ilw_slab_thick_tag)%deprecated_version = '2.0.0'
    tag(ilw_slab_thick_tag)%deprecated_name = 'LW_NUM_LAYERS'
    tag(ilw_slab_thick_tag)%description = &
         'Defines the number of primitive layers to use for the lower crystal'


    tag(iup_slab_thick_tag)%name    = 'UP_SLAB_THICKNESS'
    tag(iup_slab_thick_tag)%type    = 'I'
    tag(iup_slab_thick_tag)%summary = 'Number of layers of upper crystal'
    tag(iup_slab_thick_tag)%allowed = 'Any positive integer number'
    tag(iup_slab_thick_tag)%default = '(empty)'
    tag(iup_slab_thick_tag)%is_deprecated = .false.
    tag(iup_slab_thick_tag)%to_be_deprecated = .true.
    tag(iup_slab_thick_tag)%deprecated_version = '2.0.0'
    tag(iup_slab_thick_tag)%deprecated_name = 'UP_NUM_LAYERS'
    tag(iup_slab_thick_tag)%description = &
         'Defines the number of primitive layers to use for the upper crystal'

  end function setup_depr_interface_tags
!!!#############################################################################



!!!#############################################################################
!!! settings card help
!!!#############################################################################
  subroutine settings_help(unit, helpword, search)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: helpword
    type(tag_type), dimension(ntags_settings) :: tag
    logical :: lsearch
    logical, optional :: search
    
    lsearch=.false.
    if(present(search)) lsearch=search

    tag=setup_settings_tags()

    write(unit,'("======================================")')
    write(unit,'("Help information in SETTINGS card:")')
    call io_print_help(unit,helpword,tag,lsearch)
    write(unit,*)

  end subroutine settings_help
!!!#############################################################################


!!!#############################################################################
!!! cell_edits card help
!!!#############################################################################
  subroutine cell_edits_help(unit, helpword, search)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: helpword
    type(tag_type), dimension(ntags_cell_edits + ntags_depr_cell_edits) :: tag
    logical :: lsearch
    logical, optional :: search
    
    lsearch=.false.
    if(present(search)) lsearch=search

    tag = [ setup_cell_edits_tags(), setup_depr_cell_edits_tags() ]

    write(unit,'("======================================")')
    write(unit,'("Help information in CELL_EDITS card:")')
    call io_print_help(unit,helpword,tag,lsearch)
    write(unit,*)

  end subroutine cell_edits_help
!!!#############################################################################


!!!#############################################################################
!!! interface card help
!!!#############################################################################
  subroutine interface_help(unit, helpword, search)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: helpword
    type(tag_type), dimension(ntags_interface + ntags_depr_interface) :: tag
    logical :: lsearch
    logical, optional :: search
    
    lsearch=.false.
    if(present(search)) lsearch=search

    tag = [ setup_interface_tags(), setup_depr_interface_tags() ]

    write(unit,'("======================================")')
    write(unit,'("Help information in INTERFACE card:")')
    call io_print_help(unit,helpword,tag,lsearch)
    write(unit,*)

  end subroutine interface_help
!!!#############################################################################

end module mod_help
