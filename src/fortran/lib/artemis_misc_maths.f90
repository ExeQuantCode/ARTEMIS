module artemis__misc_maths
  !! Module containing miscellaneous maths functions and subroutines.
  !!
  !! Includes functions for array multiplication, Gaussian evaluation,
  !! factorial, log summation, safe inverse trigonometry, overlap computation,
  !! convolution, cross-correlation, running averages, statistical measures,
  !! turning point detection, plane identification, and distribution functions.
  use artemis__constants, only: real32
  implicit none
  integer, parameter :: QuadInt_K = selected_int_kind (16)
  !! Kind parameter for quadruple-precision integers.




contains

!###############################################################################
  function gauss(pos,centre,sigma,tol) result(output)
    !! Evaluate a Gaussian at a point
    implicit none

    ! Arguments
    real(real32) :: pos
    !! Position to evaluate the Gaussian at
    real(real32) :: centre
    !! Centre of the Gaussian
    real(real32) :: sigma
    !! Width of the Gaussian
    real(real32), intent(in), optional :: tol
    !! Tolerance for the Gaussian

    real(real32) :: output
    !! Output value of the Gaussian

    ! Local variables
    real(real32) :: x
    !! Squared distance from the centre
    real(real32) :: tol_
    !! Tolerance for the Gaussian

    tol_ = 38._real32
    if(present(tol)) tol_ = tol

    x = ( pos - centre ) ** 2._real32 / ( 2._real32 * sigma )
    if( abs(x) .lt. tol_ ) then
       output = exp( -x )
    else
       output = 0._real32
    end if

  end function gauss
!###############################################################################


!###############################################################################
  real(real32) function lnsum(n) 
    !! Return the sum of log(i) for i from 1 to n.
    implicit none
    integer :: i,n
    !! Loop index and upper limit.

    lnsum=0
    do i=1,n
       lnsum=lnsum+log(real(i))
    end do

    return
  end function lnsum
!###############################################################################


!###############################################################################
!###############################################################################
!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!###############################################################################
!###############################################################################


!###############################################################################
  function overlap_indiv_points(f,g) result(overlap)
    !! Compute the element-wise overlap between two arrays.
    implicit none
    integer :: n
    !! Loop index.
    integer :: datsize_f, datsize_g
    !! Sizes of input arrays.
    real(real32), dimension(:) :: f, g
    !! Input arrays.
    real(real32), dimension(:), allocatable :: overlap, y
    !! Result overlap array and temporary work array.
    
    datsize_f = size(f)
    datsize_g = size(g)

    allocate(y(datsize_f))
    if(allocated(overlap)) deallocate(overlap)
    allocate(overlap(datsize_f))

    do n=1,datsize_f
       y(n) = min(f(n),g(n))
    end do

    overlap = y

    
  end function overlap_indiv_points
!###############################################################################


!###############################################################################
!###############################################################################
!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!###############################################################################
!###############################################################################


!###############################################################################
  function running_avg(in_array,window,lperiodic) result(out_array)
    !! Smooth a function using a running average.
    implicit none
    integer :: i,lw,up,nstep
    !! Loop index, lower/upper window bounds, and number of steps.
    integer, intent(in) :: window
    !! Window size for the running average.
    real(real32), dimension(:), intent(in) :: in_array
    !! Input array to smooth.
    real(real32), dimension(size(in_array,dim=1)) :: out_array
    !! Smoothed output array.
    logical, optional :: lperiodic
    !! Whether the input is periodic.
    
    nstep=size(in_array)
    if(mod(dble(window),2.0).eq.0.0)then
       lw = nint(dble(window)/2.0)-1
       up = nint(dble(window)/2.0)
    else 
       lw = (floor(dble(window)/2.0))
       up = (floor(dble(window)/2.0))
    end if

    out_array=0.0
    if(present(lperiodic))then
       if(lperiodic)then
          do i=1,lw
             out_array(i)=sum(in_array(nstep-lw+i:nstep))+&
                  sum(in_array(1:i+up))
          end do
          do i=lw+1,nstep-up
             out_array(i)=sum(in_array(i-lw:i+up))
          end do
          do i=nstep-up+1,nstep
             out_array(i)=sum(in_array(i-lw:nstep))+&
                  sum(in_array(1:up-(nstep-i)))
          end do
          out_array=out_array/window
          return
       end if
    end if
    out_array=in_array

  end function running_avg
!###############################################################################


!###############################################################################
  function mean(in_array)
    !! Return the mean of a set of points.
    implicit none
    real(real32) :: mean
    !! Mean result.
    real(real32), dimension(:), intent(in) :: in_array
    !! Input array.

    mean=sum(in_array)/size(in_array)

  end function mean
!###############################################################################


!###############################################################################
  function median(in_array)
    !! Return the median of a set of points.
    implicit none
    integer :: i,loc
    !! Loop index and location of minimum element.
    real(real32) :: median,oddeven,rtmp1
    !! Median result, even/odd check value, and temporary variable.
    real(real32), allocatable, dimension(:) :: cp_array
    !! Copy of input array for sorting.
    real(real32), dimension(:), intent(in) :: in_array
    !! Input array.

    allocate(cp_array(size(in_array)))
    cp_array=in_array
    do i=1,size(cp_array)
       loc=minloc(cp_array(i:),dim=1)+i-1
       rtmp1=cp_array(i)
       cp_array(i)=cp_array(loc)
       cp_array(loc)=rtmp1
    end do

    oddeven=size(cp_array)/2.0
    if(abs(nint(oddeven)-oddeven).lt.tiny(0.0))then
       median=cp_array(nint(oddeven))
    else
       median=(cp_array(floor(oddeven)) + cp_array(ceiling(oddeven)))/2.0
    end if

  end function median
!###############################################################################


!###############################################################################
  function mode(in_array)
    !! Return the mode of a set of points.
    !!
    !! Currently only finds one mode, even if set is bimodal or multimodal.
    implicit none
    integer :: i,itmp1,maxcount
    !! Loop index, temporary count, and maximum count.
    real(real32) :: mode
    !! Mode result.
    real(real32), dimension(:), intent(in) :: in_array
    !! Input array.

    maxcount=0
    do i=1,size(in_array)
       itmp1=count(in_array.eq.in_array(i))
       itmp1=count(abs(in_array-in_array(i)).lt.1.E-8_real32)
       if(itmp1.gt.maxcount)then
          maxcount=itmp1
          mode=in_array(i)
       end if
    end do

  end function mode
!###############################################################################


!###############################################################################
  function range(in_array) result(output)
    !! Return the range of a set of points.
    implicit none
    real(real32) :: output
    !! Range result.
    real(real32), dimension(:), intent(in) :: in_array
    !! Input array.

    output=maxval(in_array)-minval(in_array)

  end function range
!###############################################################################


!###############################################################################
  function normalise(in_array) result(output)
    !! Return an array normalised to one.
    implicit none
    real(real32) :: sumval
    !! Sum of the input array elements.
    real(real32), dimension(:), intent(in) :: in_array
    !! Input array.
    real(real32), dimension(size(in_array)) :: output
    !! Normalised output array.
    
    sumval=sum(in_array)
    if(sumval.lt.1.E-8_real32)then
       output=in_array
    else
       output=in_array/sum(in_array)
    end if

  end function normalise
!###############################################################################


!###############################################################################
  function get_turn_points(invec,lperiodic,window) result(resvec)
    !! Find turning points, saving them in order of smallest to largest.
    !!
    !! Note: should check the turning point is sustained across the window.
    implicit none
    integer :: i,j,nturn,itmp1,itmp2
    !! Loop indices, number of turning points, and temporary variables.
    real(real32) :: l_grad,r_grad
    !! Left and right gradients.
    real(real32), dimension(:), intent(in) :: invec
    !! Input vector.
    integer, allocatable, dimension(:) :: tvec1,resvec
    !! Temporary turning point indices and result vector.
    integer, optional :: window
    !! Window size for reducing close turning points.
    logical, optional :: lperiodic
    !! Whether the input is periodic.


    nturn=0
    if(allocated(resvec)) deallocate(resvec)
    allocate(tvec1(size(invec)))
    l_grad=0._real32
    r_grad=invec(2)-invec(1)
    if(present(lperiodic))then
       if(lperiodic)then
          l_grad=invec(1)-invec(size(invec))
          if(sign(1._real32,l_grad).ne.sign(1._real32,r_grad).or.&
               (r_grad.eq.0._real32.and.l_grad.ne.r_grad))then
             nturn=nturn+1
             tvec1(nturn)=1
          end if
       end if
    end if


    do i=2,size(invec)-1
       l_grad=r_grad
       r_grad=invec(i+1)-invec(i)
       if(sign(1._real32,l_grad).ne.sign(1._real32,r_grad).or.&
            (r_grad.eq.0._real32.and.abs(l_grad-r_grad).gt.1.E-5_real32))then
          nturn=nturn+1
          tvec1(nturn)=i
       end if
    end do


    if(present(lperiodic))then
       if(lperiodic)then
          r_grad=invec(1)-invec(size(invec))
          if(sign(1._real32,l_grad).ne.sign(1._real32,r_grad).or.&
               (r_grad.eq.0._real32.and.l_grad.ne.r_grad))then
             nturn=nturn+1
             tvec1(nturn)=size(invec)
          end if
       end if
    end if

    if(present(window))then
       i=1
       reduceloop:do
          if(i.ge.nturn) exit reduceloop
          if(abs(tvec1(i)-tvec1(i+1)).lt.window)then
             itmp1=minloc((/invec(tvec1(i)),invec(tvec1(i+1))/),dim=1)
             tvec1(i+itmp1-1:nturn-1)=tvec1(i+itmp1:nturn)
             nturn=nturn-1
          else
             i=i+1
          end if
       end do reduceloop
    end if


    allocate(resvec(nturn))
    resvec(:nturn)=tvec1(:nturn)
    do i=1,nturn
       itmp1=minloc((/  (invec(resvec(j)),j=i,nturn)  /),dim=1)+i-1
       itmp2=resvec(i)
       resvec(i)=resvec(itmp1)
       resvec(itmp1)=itmp2
    end do

    
  end function get_turn_points
!###############################################################################


!###############################################################################
  function get_nth_plane(invec,nth,window,is_periodic) result(startend)
    !! Find the location of the nth plane as start and end coordinates.
    implicit none
    integer :: i,nstep,nplane,udef_window
    !! Loop index, number of steps, plane count, and user-defined window.
    real(real32) :: tol
    !! Tolerance for plane height variation.
    logical :: is_in_plane
    !! Flag indicating if currently within a plane.
    integer, dimension(2) :: startend
    !! Start and end indices of the nth plane.
    integer, allocatable, dimension(:,:) :: plane_loc
    !! Locations of identified planes.
    real(real32), dimension(:), intent(in) :: invec
    !! Input vector.
    integer, intent(in) :: nth
    !! Index of the plane to find.
    integer, optional, intent(in) :: window
    !! Window size for plane detection.
    logical, optional, intent(in) :: is_periodic
    !! Whether the input is periodic.


!-------------------------------------------------------------------------------
! Defines tolerance of plane height variation and initialises variables
!-------------------------------------------------------------------------------
    tol = 0.01_real32*(maxval(invec)-minval(invec))
    if(present(window))then
       udef_window=window
    else
       udef_window=10
    end if

    nplane=0
    nstep = size(invec,dim=1)
    allocate(plane_loc(nstep/udef_window,2))
    plane_loc=0
    is_in_plane=.false.


!-------------------------------------------------------------------------------
! Loops over points to identify planes
!-------------------------------------------------------------------------------
    i=0
    step_loop1: do while(i.le.nstep-udef_window)
       i = i + 1
       if(is_in_plane)then
          if(all(&
               abs(invec(plane_loc(nplane,1):plane_loc(nplane,2))-invec(i)).lt.&
               tol))then
             plane_loc(nplane,2)=i
             cycle step_loop1
          end if
       end if
       
       if(.not.is_in_plane)then
          if(all(abs(invec(i:i+udef_window-1)-invec(i)).lt.tol))then
             is_in_plane=.true.
             nplane = nplane + 1
             plane_loc(nplane,1) = i
             plane_loc(nplane,2) = i + udef_window -1
             i = i + udef_window
          end if
          cycle step_loop1
       end if

       is_in_plane=.false.
       
    end do step_loop1


!-------------------------------------------------------------------------------
! Handles the last few points depending on whether set is periodic
!-------------------------------------------------------------------------------
    if(present(is_periodic))then
       if(plane_loc(nplane,2).eq.nstep-udef_window)then
          step_loop2: do i=nstep-udef_window,nstep
             if(all(&
                  abs(invec(plane_loc(nplane,1):plane_loc(nplane,2))-&
                  invec(i)).lt.tol))then
                plane_loc(nplane,2) = i
             else
                exit step_loop2
             end if
          end do step_loop2
       end if
       if(is_periodic)then
          if(plane_loc(1,1).eq.1.and.&
               plane_loc(nplane,2).eq.nstep)then
             if(all(abs(invec(:plane_loc(1,2))-invec(nstep)).lt.tol))then
                plane_loc(1,1) = plane_loc(nplane,1)
                plane_loc(nplane,:) = 0
                nplane = nplane - 1
             end if
          else
             step_loop3: do i=nstep,nstep-udef_window,-1
                if(all(abs(invec(:plane_loc(1,2))-invec(i)).lt.tol))then
                   plane_loc(1,1)=i
                end if
             end do step_loop3
          end if
       end if

    end if


!-------------------------------------------------------------------------------
! Sets value of nth plane
!-------------------------------------------------------------------------------
    if(nplane.lt.nth)then
       startend=0
    else
       startend=plane_loc(nth,:)
    end if



  end function get_nth_plane
!###############################################################################


!###############################################################################
!###############################################################################
!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!###############################################################################
!###############################################################################


!###############################################################################
  function table_func(x,a) result(res)
    !! Compute a custom table function for a single point.
    !!
    !! Note: breaks when a = 1. Above this, res will always equal 1.
    !! Parameter a should be between -1 and 1.
    implicit none
    real(real32), intent(in) :: x,a
    !! Input value and shape parameter.
    real(real32) :: res
    !! Result of the table function.

    res=( ( cos(x) + a ) + abs( cos(x) - a ) - 2._real32 )/&
         ( 2._real32*a - 2._real32 )


  end function table_func
!###############################################################################


!###############################################################################
  function gauss_array(distance,in_array,sigma,tol,norm,mask) &
       result(gauss_func)
    !! Apply Gaussian distributions to a set of points in an array.
    implicit none
    integer :: i,n,init_step
    !! Loop indices and initial step position.
    real(real32) :: x,sigma,tol_,mult
    !! Squared distance, width, tolerance, and normalisation multiplier.
    real(real32), optional :: tol
    !! Optional tolerance for the Gaussian.
    logical, optional :: norm
    !! Optional flag to control normalisation.
    real(real32), dimension(:), intent(in) :: in_array,distance
    !! Input array of points and distance array.
    real(real32), dimension(size(distance)) :: gauss_func
    !! Resulting Gaussian function array.
    real(real32) :: pi = 4._real32*atan(1._real32)
    !! Value of pi.

    logical, dimension(size(distance)), optional, intent(in) :: mask
    !! Optional mask array for selective evaluation.


    tol_ = 38._real32
    if(present(tol)) tol_ = tol
    mult=(1._real32/(sqrt(pi*2._real32)*sigma))
    if(present(norm))then
       if(.not.norm) mult=1._real32
    end if
    
    gauss_func=0._real32
    do n=1,size(in_array)
       if(present(mask))then
          if(.not.mask(n)) cycle
       end if
       init_step=minloc(abs( distance(:) - in_array(n) ),dim=1)
       forward: do i=init_step,size(distance),1
          x=0.5_real32*(( distance(i) - in_array(n) )/sigma)**2._real32
          if(x.gt.tol_) exit forward
          gauss_func(i) = gauss_func(i) + exp(-x) * mult
       end do forward

       backward: do i=init_step-1,1,-1
          x=0.5_real32*(( distance(i) - in_array(n) )/sigma)**2._real32
          if(x.gt.tol_) exit backward
          gauss_func(i) = gauss_func(i) + exp(-x) * mult
       end do backward
    end do



  end function gauss_array
!###############################################################################


!###############################################################################
  function cauchy_array(distance,in_array,gamma,tol,norm) result(c_func)
    !! Apply Cauchy distributions to a set of points in an array.
    implicit none
    integer :: i,n,init_step
    !! Loop indices and initial step position.
    real(real32) :: x,gamma,tol_,mult
    !! Distance value, scale parameter, tolerance, and normalisation multiplier.
    real(real32), optional :: tol
    !! Optional tolerance for the Cauchy distribution.
    logical, optional :: norm
    !! Optional flag to control normalisation.
    real(real32), dimension(:), intent(in) :: in_array,distance
    !! Input array of points and distance array.
    real(real32), dimension(size(distance)) :: c_func
    !! Resulting Cauchy function array.
    real(real32) :: pi = 4._real32*atan(1._real32)
    !! Value of pi.


    tol_ = 1.E16_real32
    if(present(tol)) tol_=tol
    mult=(1._real32/(pi*gamma))
    if(present(norm))then
       if(.not.norm) mult=1._real32
    end if
    
    c_func=0._real32
    do n=1,size(in_array)
       init_step=minloc(abs( distance(:) - in_array(n) ),dim=1)
       forward: do i=init_step,size(distance),1
          x = 1._real32 + (( distance(i) - in_array(n) )/gamma)**2._real32
          if(x.gt.tol_) exit forward
          c_func(i) = c_func(i) + 1._real32/(x) * mult
       end do forward

       backward: do i=init_step-1,1,-1
          x = 1._real32 + (( distance(i) - in_array(n) )/gamma)**2._real32
          if(x.gt.tol_) exit backward
          c_func(i) = c_func(i) + 1._real32/x * mult
       end do backward
    end do



  end function cauchy_array
!###############################################################################


!###############################################################################
  function slater_array(distance,in_array,zeta,tol,norm) result(s_func)
    !! Apply Slater distributions to a set of points in an array.
    implicit none
    integer :: i,n,init_step
    !! Loop indices and initial step position.
    real(real32) :: x,zeta,tol_,mult
    !! Distance value, Slater exponent, tolerance, and normalisation multiplier.
    real(real32), optional :: tol
    !! Optional tolerance for the Slater distribution.
    logical, optional :: norm
    !! Optional flag to control normalisation.
    real(real32), dimension(:), intent(in) :: in_array,distance
    !! Input array of points and distance array.
    real(real32), dimension(size(distance)) :: s_func
    !! Resulting Slater function array.
    real(real32) :: pi = 4._real32*atan(1._real32)
    !! Value of pi.


    tol_ = 38._real32
    if(present(tol)) tol_=tol
    mult=((zeta**3._real32)/pi)**(0.5_real32)
    if(present(norm))then
       if(.not.norm) mult=1._real32
    end if
    
    s_func=0._real32
    do n=1,size(in_array)
       init_step=minloc(abs( distance(:) - in_array(n) ),dim=1)
       forward: do i=init_step,size(distance),1
          x = zeta*abs( distance(i) - in_array(n) )
          if(x.gt.tol_) exit forward
          s_func(i) = s_func(i) + exp(-x) * mult
       end do forward

       backward: do i=init_step-1,1,-1
          x = zeta*abs( distance(i) - in_array(n) )
          if(x.gt.tol_) exit backward
          s_func(i) = s_func(i) + exp(-x) * mult
       end do backward
    end do


  end function slater_array
!###############################################################################

end module artemis__misc_maths
