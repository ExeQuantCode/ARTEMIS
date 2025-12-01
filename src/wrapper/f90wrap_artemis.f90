! Module artemis defined in file ../fortran/artemis.f90

subroutine f90wrap_get_suppress_warnings(suppress_warnings)
    use artemis, only: artemis__suppress_warnings
    implicit none
    logical, intent(out) :: suppress_warnings
    
    suppress_warnings = artemis__suppress_warnings
end subroutine f90wrap_get_suppress_warnings

subroutine f90wrap_set_suppress_warnings(suppress_warnings)
    use artemis, only: artemis__suppress_warnings
    implicit none
    logical, intent(in) :: suppress_warnings
    
    artemis__suppress_warnings = suppress_warnings
end subroutine f90wrap_set_suppress_warnings

! End of module artemis defined in file ../fortran/artemis.f90

