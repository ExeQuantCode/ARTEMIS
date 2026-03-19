module artemis__misc_linalg
  !! Module containing various linear algebra functions and subroutines.
  !!
  !! Includes vector operations (unit vector, projection, Gram-Schmidt,
  !! cross product matrix, outer product), matrix operations (determinant,
  !! inverse, trace, LU decomposition), geometric computations (angle, area,
  !! volume), equation solvers (simultaneous equations, transformation matrix),
  !! lattice reduction (LLL), vector rotation, and integer arithmetic utilities
  !! (GCD, LCM, fraction denominator, vector GCD reduction, group generation).
  use artemis__constants, only: real32
  use coreutils__linalg, only: cross, inverse_3x3, outer_product
  implicit none
  integer, parameter, private :: QuadInt_K = selected_int_kind (16)

  interface gcd
     procedure gcd_vec,gcd_num
  end interface gcd

  interface vec_mat_mul
     procedure ivec_dmat_mul,dvec_dmat_mul
  end interface vec_mat_mul

  interface det
     procedure idet,ddet,rec_det
  end interface det


contains
!###############################################################################
  function uvec(vec) result(output)
    !! Return the unit vector of an arbitrary-size vector.
    implicit none

    ! Arguments
    real(real32),dimension(:)::vec
    !! Input vector.
    real(real32),allocatable,dimension(:) :: output
    !! Unit vector result.

    allocate(output(size(vec)))
    output = vec/norm2(vec)
  end function uvec
!###############################################################################


!###############################################################################
  function proj(u,v) result(output)
    !! Return the projection of vector v onto vector u.
    implicit none

    ! Arguments
    real(real32), dimension(:) :: u,v
    !! Input vectors.
    real(real32), allocatable, dimension(:) :: output
    !! Projected vector result.

    allocate(output(size(u,dim=1)))
    output = u*dot_product(v,u)/dot_product(u,u)

  end function proj
!###############################################################################


!###############################################################################
  function GramSchmidt(basis,normalise,cmo) result(u)
    !! Evaluate the Gram-Schmidt orthogonal basis.
    !!
    !! Assumes basis(n,m) is a basis of n vectors, each of m dimensions.
    implicit none
    integer :: num,dim,i,j
    !! Number of vectors, dimension size, and loop counters.
    real(real32), allocatable, dimension(:) :: vtmp
    !! Temporary vector for accumulating projections.
    real(real32), dimension(:,:), intent(in) :: basis
    !! Input basis matrix (n vectors x m dimensions).
    real(real32), allocatable, dimension(:,:) :: u
    !! Output orthogonal basis.
    logical, optional, intent(in) :: cmo
    !! If true, use column major order (not yet implemented).
    logical, optional, intent(in) :: normalise
    !! If true, normalise the output basis.


    !! sets up array dimensions of Gram-Schmidt basis
    if(present(cmo))then
       if(cmo)then
          write(0,'("Column Major Order Gram-Schmidt &
               &not yet set up")')
          write(0,'("Stopping...")')
          stop
          num = size(basis(1,:),dim=1)
          dim = size(basis(:,1),dim=1)
          allocate(u(dim,num))
          goto 10
       end if
    end if
    num = size(basis(:,1),dim=1)
    dim = size(basis(1,:),dim=1)
    allocate(u(num,dim))
    
10  allocate(vtmp(dim))

    !! Evaluates the Gram-Schmidt basis
    u(1,:) = basis(1,:)
    do i=2,num
       vtmp = 0._real32
       do j=1,i-1,1
          vtmp(:) = vtmp(:) + proj(u(j,:),basis(i,:))
       end do
       u(i,:) = basis(i,:) - vtmp(:)
    end do


    !! Normalises new basis if required
    if(present(normalise))then
       if(normalise)then
          do i=1,num
             u(i,:) = u(i,:)/norm2(u(i,:))
          end do
       end if
    end if


  end function GramSchmidt
!###############################################################################


!###############################################################################
  function cross_matrix(a)
    !! Generate the cross product matrix of a 3D vector.
    !!
    !! For a = (a1,a2,a3), returns the skew-symmetric matrix [a]_x such that
    !! [a]_x * b = a x b.
    implicit none

    ! Arguments
    real(real32), dimension(3,3) :: cross_matrix
    !! The resulting 3x3 cross product matrix.
    real(real32), dimension(3), intent(in) :: a
    !! Input 3D vector.

    cross_matrix=0._real32

    cross_matrix(1,2) = -a(3)
    cross_matrix(1,3) =  a(2)
    cross_matrix(2,3) = -a(1)

    cross_matrix(2,1) =  a(3)
    cross_matrix(3,1) = -a(2)
    cross_matrix(3,2) =  a(1)

    return
  end function cross_matrix
!###############################################################################


!###############################################################################
  function ivec_dmat_mul(a,mat) result(vec)
    !! Multiply an integer vector with a real matrix.
    implicit none
    integer :: j
    !! Loop counter.
    integer, dimension(:) :: a
    !! Input integer vector.
    real(real32), dimension(:,:) :: mat
    !! Input real matrix.
    real(real32),allocatable,dimension(:) :: vec
    !! Output real vector.

    vec=0._real32
    allocate(vec(size(a)))
    do j=1,size(a)
       vec(:)=vec(:)+real(a(j),real32)*mat(j,:)
    end do

    return
  end function ivec_dmat_mul
!---------------------------------------------------------------------------
  function dvec_dmat_mul(a,mat) result(vec)
    !! Multiply a real vector with a real matrix.
    implicit none
    integer :: j
    !! Loop counter.
    real(real32), dimension(:) :: a
    !! Input real vector.
    real(real32), dimension(:,:) :: mat
    !! Input real matrix.
    real(real32),allocatable,dimension(:) :: vec
    !! Output real vector.

    vec=0._real32
    allocate(vec(size(a)))
    do j=1,size(a)
       vec(:)=vec(:)+a(j)*mat(j,:)
    end do

    return
  end function dvec_dmat_mul
!###############################################################################


!###############################################################################
  function get_vec_multiple(a,b) result(multi)
    !! Determine the scaling factor between two vectors.
    !!
    !! Returns the scalar multi such that b = multi * a, or 0 if no
    !! consistent scalar exists.
    implicit none
    integer :: i
    !! Loop counter.
    real(real32) :: multi
    !! Scaling factor result.
    real(real32), dimension(:) :: a,b
    !! Input vectors.
    
    multi=1._real32
    do i=1,size(a)
       if(abs(a(i)).lt.1.E-6_real32.or.abs(b(i)).lt.1.E-6_real32) cycle
       multi=b(i)/a(i)
       exit
    end do

    checkloop: do i=1,size(a)
       if(abs(a(i)).lt.1.E-6_real32.or.abs(b(i)).lt.1.E-6_real32) cycle
       if(abs(a(i)*multi-b(i)).gt.1.E-6_real32)then

          multi=0._real32
          exit checkloop
       end if
    end do checkloop

    return
  end function get_vec_multiple
!###############################################################################


!###############################################################################
!###############################################################################
!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!###############################################################################
!###############################################################################


!###############################################################################
  function get_angle(vec1,vec2) result(angle)
    !! Return the angle between two 3D vectors in radians.
    implicit none

    ! Arguments
    real(real32) :: angle
    !! Angle result.
    real(real32), dimension(3) :: vec1,vec2
    !! Input 3D vectors.

    angle = acos( dot_product(vec1,vec2)/&
         ( norm2(vec1) * norm2(vec2) ))
    if (isnan(angle)) angle = 0._real32

    return
  end function get_angle
!###############################################################################


!###############################################################################
  function get_area(a,b) result(area)
    !! Return the area of the parallelogram formed by two 3D vectors.
    implicit none
    real(real32) :: area
    !! Area result.
    real(real32), dimension(3) :: vec,a,b
    !! Cross product vector and input vectors.

    vec = cross(a,b)
    area = sqrt(dot_product(vec,vec))

    return
  end function get_area
!###############################################################################


!###############################################################################
  function get_vol(lat) result(vol)
    !! Return the volume of a 3x3 lattice matrix.
    implicit none
    integer :: n,i,j,k,l
    !! Loop counters and permutation indices.
    real(real32) :: vol,scale
    !! Volume result and sign scaling factor.
    real(real32), dimension(3,3) :: lat
    !! Input 3x3 lattice matrix.
    real(real32), dimension(3) :: a,b,c
    !! Lattice row vectors.

    a=lat(1,:)
    b=lat(2,:)
    c=lat(3,:)
    vol = 0._real32;scale = 1._real32
    i=1;j=2;k=3
1   do n=1,3
       vol = vol+scale*a(i)*b(j)*c(k)
       l=i;i=j;j=k;k=l
    end do
    i=2;j=1;k=3;scale=-scale
    if(scale<0._real32) goto 1

    return
  end function get_vol
!###############################################################################


!###############################################################################
  function trace(mat) result(output)
    !! Return the trace of an arbitrary-dimension square matrix.
    integer::j
    !! Loop counter.
    real(real32), dimension(:,:), intent(in) :: mat
    !! Input square matrix.
    real(real32) :: output
    !! Trace result.
    output = 0._real32
    do j = 1, size(mat,1)
      output = output + mat(j,j)
    end do
  end function trace
!###############################################################################


!###############################################################################
  function idet(mat) result(output)
    !! Return the determinant of a 3x3 integer matrix.
    integer :: output
    !! Determinant result.
    integer, dimension(3,3), intent(in) :: mat
    !! Input 3x3 integer matrix.

    output = mat(1,1)*mat(2,2)*mat(3,3)-mat(1,1)*mat(2,3)*mat(3,2)&
         - mat(1,2)*mat(2,1)*mat(3,3)+mat(1,2)*mat(2,3)*mat(3,1)&
         + mat(1,3)*mat(2,1)*mat(3,2)-mat(1,3)*mat(2,2)*mat(3,1)

  end function idet
!---------------------------------------------------------------------------
  function ddet(mat) result(output)
    !! Return the determinant of a 3x3 real matrix.
    real(real32) :: output
    !! Determinant result.
    real(real32), dimension(3,3), intent(in) :: mat
    !! Input 3x3 real matrix.

    output = mat(1,1)*mat(2,2)*mat(3,3)-mat(1,1)*mat(2,3)*mat(3,2)&
         - mat(1,2)*mat(2,1)*mat(3,3)+mat(1,2)*mat(2,3)*mat(3,1)&
         + mat(1,3)*mat(2,1)*mat(3,2)-mat(1,3)*mat(2,2)*mat(3,1)

  end function ddet
!###############################################################################


!###############################################################################
  pure function inverse(mat)
    !! Return the inverse of a 2x2 or 3x3 matrix.
    real(real32), dimension(:,:), intent(in) :: mat
    !! Input square matrix (2x2 or 3x3).
    real(real32), dimension(size(mat,dim=1),size(mat,dim=2)) :: inverse
    !! Inverse matrix result.

    select case(size(mat,dim=2))
    case(2)
       inverse = inverse_2x2(mat)
    case(3)
       inverse = inverse_3x3(mat)
    end select

  end function inverse
!###############################################################################


!###############################################################################
  pure function inverse_2x2(mat) result(output)
    !! Return the inverse of a 2x2 matrix.
    implicit none
    real(real32), dimension(2,2), intent(in) :: mat
    !! Input 2x2 matrix.
    real(real32), dimension(2,2) :: output
    !! Inverse matrix result.
    real(real32) :: inv_det
    !! Reciprocal of the determinant.

    associate(a => mat(1,1), b => mat(1,2), c => mat(2,1), d => mat(2,2))
       inv_det = 1._real32 / (a * d - b * c)

       output(1,1) =  d * inv_det
       output(1,2) = -b * inv_det
       output(2,1) = -c * inv_det
       output(2,2) =  a * inv_det
    end associate

  end function inverse_2x2
!###############################################################################


!###############################################################################
  recursive function rec_det(a,n) result(res)
    !! Return the determinant of an n x n matrix using cofactor expansion.
    integer :: i, sign
    !! Cofactor column index and sign toggle.
    real(real32) :: res
    !! Determinant result.
    integer, intent(in) :: n
    !! Matrix dimension.
    real(real32), dimension(n,n), intent(in) :: a
    !! Input n x n matrix.
    real(real32), dimension(n-1, n-1) :: tmp
    !! Submatrix for cofactor expansion.

    if(n.eq.1) then
       res = a(1,1)
    else
       res = 0._real32
       sign = 1
       do i=1, n
          tmp(:,:(i-1))=a(2:,:i-1)
          tmp(:,i:)=a(2:,i+1:)
          res=res+sign*a(1,i)*rec_det(tmp,n-1)
          sign=-1*sign
       end do
    end if

    return
  end function rec_det
!###############################################################################


!###############################################################################
  function LUdet(inmat)
    !! Return the determinant of an n x n matrix via LU decomposition.
    !!
    !! Computes LUdet = (-1)^N * prod(L(i,i)*U(i,i)).
    implicit none
    integer :: i,N
    !! Loop counter and matrix dimension.
    real(real32) :: LUdet
    !! Determinant result.
    real(real32), dimension(:,:) :: inmat
    !! Input n x n matrix.
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: L,U
    !! Lower and upper triangular matrices.

    L=0._real32
    U=0._real32
    N=size(inmat,1)
    call LUdecompose(inmat,L,U)

    LUdet=(-1._real32)**N
    do i=1,N
       LUdet=LUdet*L(i,i)*U(i,i)
    end do

    return
  end function LUdet
!###############################################################################


!###############################################################################
  function LUinv(inmat)
    !! Return the inverse of an n x n matrix using LU decomposition.
    !!
    !! Does not work if a diagonal element is zero.
    implicit none
    integer :: i,m,N
    !! Loop counters and matrix dimension.
    real(real32), dimension(:,:) :: inmat
    !! Input n x n matrix.
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: LUinv
    !! Inverse matrix result.
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: L,U
    !! Lower and upper triangular matrices.
    real(real32), dimension(size(inmat,1)) :: c,z,x
    !! Identity column vector, intermediate vector, and solution vector.

    L=0._real32
    U=0._real32
    N=size(inmat,1)
    call LUdecompose(inmat,L,U)

! Lz=c
! c are column vectors of the identity matrix
! uses forward substitution to solve
    do m=1,N
       c=0._real32
       c(m)=1._real32

       z(1)=c(1)
       do i=2,N
          z(i)=c(i)-dot_product(L(i,1:i-1),z(1:i-1))
       end do


! Ux=z
! x are the rows of the inversion matrix
! uses backwards substitution to solve
       x(N)=z(N)/U(N,N)
       do i=N-1,1,-1
          x(i)=z(i)-dot_product(U(i,i+1:N),x(i+1:N))
          x(i)= x(i)/U(i,i)
       end do

       LUinv(:,m)=x(:)
    end do

    return
  end function LUinv
!###############################################################################


!###############################################################################
  subroutine LUdecompose(inmat,L,U)
    !! Decompose a matrix into lower and upper triangular matrices (A = LU).
    !!
    !! Based on Doolittle LU factorization for Ax=b.
    !! Does not work if a diagonal element is zero.
    implicit none
    integer :: i,j,N
    !! Loop counters and matrix dimension.
    real(real32), dimension(:,:) :: inmat,L,U
    !! Input matrix, lower and upper triangular output matrices.
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: mat
    !! Working copy of input matrix.

    N=size(inmat,1)
    mat=inmat
    L=0._real32
    U=0._real32

    do j=1,N
       L(j,j)=1._real32
    end do
! Solves the lower matrix
    do j=1,N-1
       do i=j+1,N
          L(i,j)=mat(i,j)/mat(j,j)
          mat(i,j+1:N)=mat(i,j+1:N)-L(i,j)*mat(j,j+1:N)
       end do
    end do

! Equates upper half of remaining mat to upper matrix
    do j=1,N
       do i=1,j
          U(i,j)=mat(i,j)
       end do
    end do

    return
  end subroutine LUdecompose
!###############################################################################


!###############################################################################
!###############################################################################
!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!###############################################################################
!###############################################################################


!###############################################################################
  function find_tf(mat1,mat2) result(tf)
    !! Find the transformation matrix between two matrices.
    !!
    !! Computes T = A^-1 * B where A = mat1 and B = mat2.
    implicit none

    ! Arguments
    real(real32), dimension(:,:) :: mat1,mat2
    !! Input matrices A and B.
    real(real32), dimension(size(mat1,dim=1),size(mat1,dim=2)) :: tf
    !! Transformation matrix result.

    tf=matmul(inverse(mat1),mat2)

  end function find_tf
!---------------------------------------------------------------------------
  function find_tf_2x2(mat1,mat2) result(tf)
    !! Find the transformation matrix between two 2x2 matrices.
    implicit none

    ! Arguments
    real(real32), dimension(2,2) :: mat1,mat2
    !! Input 2x2 matrices.
    real(real32), dimension(2,2) :: tf
    !! Transformation matrix result.

    tf=matmul(inverse_2x2(mat1),mat2)

  end function find_tf_2x2
!###############################################################################


!###############################################################################
  function simeq(qX,qY)
    !! Solve simultaneous equations for n dimensions.
    !!
    !! Given x values qX and y values qY, returns the coefficients of the
    !! power series f(qX)=qY. Highest power coefficient is simeq(1).
    integer :: i,j,n,loc
    !! Loop counters, equation order, and pivot location.
    real(real32), dimension(:) :: qX,qY
    !! Input x values and y values.
    real(real32), dimension(size(qY)) :: funcY
    !! Working copy of y values.
    real(real32), dimension(size(qY)) :: simeq,tmpqY
    !! Coefficient result and temporary y values.
    real(real32), dimension(size(qY),size(qY)) :: P,invP,tmpP
    !! Power series matrix, its inverse, and temporary copy.


    n=size(qX)
    funcy=qY
    P=0._real32
    do i=1,n
       do j=1,n
          P(i,j)=(qX(i)**real(n-j,real32))
       end do
    end do
    !  P(1,1)=qX(1)**2 ;P(1,2)=qX(1)   ;P(1,3)=1.0;
    !  P(2,1)=qX(2)**2 ;P(2,2)=qX(2)   ;P(2,3)=1.0;
    !  P(3,2)=qX(3)**2 ;P(3,2)=qX(3)   ;P(3,3)=1.0;

    if(any(qX.lt.1.E-5_real32)) then
       loc=minloc(abs(qX),dim=1)
       tmpqY=funcY
       tmpP=P
       funcY(loc)=tmpqY(n)
       funcY(n)=tmpqY(loc)
       P(loc,:)=tmpP(n,:)
       P(n,:)=tmpP(loc,:)
    end if

    !  invP=inverse(P)
    invP=LUinv((P))
    !  invP=LUinv(real(P,real32))
    simeq=matmul(invP,funcY)

  end function simeq
!###############################################################################


!###############################################################################
  function LLL_reduce(basis,delta) result(obas)
    !! Perform Lenstra-Lenstra-Lovasz (LLL) lattice basis reduction.
    !!
    !! LLL algorithm based on Hoffstein, Pipher and Silverman 2008.
    implicit none
    integer :: num,dim,i,j,k,loc
    !! Number of vectors, dimension, loop counters, and location index.
    real(real32) :: d,dtmp
    !! Delta parameter and temporary scalar.
    real(real32), allocatable, dimension(:) :: vtmp,mag_bas
    !! Temporary vector and basis magnitudes.
    real(real32), allocatable, dimension(:,:) :: mu,GSbas,obas
    !! Gram-Schmidt coefficients, orthogonal basis, and output basis.

    real(real32), dimension(:,:), intent(in) :: basis
    !! Input basis to reduce.
    real(real32), optional, intent(in) :: delta
    !! LLL reduction parameter (default 0.75).


    !! set up the value for delta
    if(present(delta))then
       d = delta
    else
       d = 0.75_real32
    end if
    
    !! allocate and initialise arrays
    num = size(basis(:,1),dim=1)
    dim = size(basis(1,:),dim=1)
    allocate(vtmp(dim))
    allocate(mag_bas(num))
    allocate(obas(num,dim))
    obas = basis

    !! reduce the gcd of the vectors
    do i=1,num
       obas(i,:) = reduce_vec_gcd(obas(i,:))
       mag_bas(i) = norm2(obas(i,:))
    end do
    
    !! sort basis such that b1 is smallest
    do i=1,num-1,1
       loc = maxloc(mag_bas(i:num),dim=1) + i - 1
      if(loc.eq.i) cycle
       dtmp = mag_bas(i)
       mag_bas(i) = mag_bas(loc)
       mag_bas(loc) = dtmp
    
       vtmp = obas(i,:)
       obas(i,:) = obas(loc,:)
       obas(loc,:) = vtmp
    end do

    !! set up Gram-Schmidt process orthogonal basis
    allocate(GSbas(num,dim))
    GSbas = GramSchmidt(obas)

    !! set up the Gram-Schmidt coefficients
    allocate(mu(num,num))
    mu = get_mu(obas,GSbas)

    !! minimise the basis
    k = 2
    do while(k.le.num)

       jloop: do j=k-1,1!,-1
          if(abs(mu(k,j)).lt.0.5_real32)then
             obas(k,:) = obas(k,:) - &
                  nint(mu(k,j))*obas(j,:)
             !! only need to update GSbas(k:,:) and mu
             !GSbas = GramSchmidt(obas)
             !mu = get_mu(obas,GSbas)
             call update_GS_and_mu(GSbas,mu,obas,k)
          end if
       end do jloop

       if(dot_product(GSbas(k,:),GSbas(k,:)).ge.&
            (d - mu(k,k-1)**2._real32)*&
            dot_product(GSbas(k-1,:),GSbas(k-1,:)) )then
          k = k + 1
       else
          vtmp = obas(k,:)
          obas(k,:) = obas(k-1,:)
          obas(k-1,:) = vtmp
          !GSbas = GramSchmidt(obas)
          !mu = get_mu(obas,GSbas)
          if(k.eq.1)then
             call update_GS_and_mu(GSbas,mu,obas,k)
          else
             call update_GS_and_mu(GSbas,mu,obas,k-1)
          end if
          k = max(k-1,2)
       end if

    end do


! Separate functions for this to run efficiently
  contains
    function get_mu(bas1,bas2) result(mu)
      !! Return the Gram-Schmidt mu coefficient matrix.
      implicit none
      integer :: num1,num2
      real(real32), allocatable, dimension(:,:) :: mu,bas1,bas2
      num1 = size(bas1(:,1),dim=1)
      num2 = size(bas2(:,1),dim=1)

      allocate(mu(num1,num2))
      do i=1,num1
         do j=1,num2

            mu(i,j) = dot_product(bas1(i,:),bas2(j,:))/&
                 dot_product(bas2(j,:),bas2(j,:))

         end do
      end do

    end function get_mu


    subroutine update_GS_and_mu(GSbas,mu,basis,k)
      !! Update Gram-Schmidt vectors and mu values from index k onwards.
      implicit none
      integer :: num,dim,i,j
      real(real32), allocatable, dimension(:) :: vtmp

      integer, intent(in) :: k
      real(real32), allocatable, dimension(:,:) :: GSbas,basis,mu

      num = size(basis(:,1),dim=1)
      dim = size(basis(1,:),dim=1)

      allocate(vtmp(dim))
      
      !!update Gram-Schmidt vectors
      do i=k,num,1
         vtmp = 0._real32
         do j=1,i-1,1
            vtmp(:) = vtmp(:) + proj(GSbas(j,:),basis(i,:))
         end do
         GSbas(i,:) = basis(i,:) - vtmp(:)
      end do


      !!update mu values
      mu_loop1: do i=1,num,1
         mu_loop2: do j=1,num,1
      
            if(i.lt.k.and.j.lt.k) cycle mu_loop2
            
            mu(i,j) = dot_product(basis(i,:),GSbas(j,:))/&
                 dot_product(GSbas(j,:),GSbas(j,:))
            
      
         end do mu_loop2
      end do mu_loop1

    end subroutine update_GS_and_mu


  end function LLL_reduce
!###############################################################################


!###############################################################################
  function rotvec(a,theta,phi,psi,new_length)
    !! Rotate a 3D vector about the x, y, and z Cartesian axes.
    implicit none
    real(real32) :: magold,theta,phi,psi
    !! Old magnitude, and rotation angles about x, y, z axes.
    real(real32), dimension(3) :: a,rotvec
    !! Input vector and rotated result.
    real(real32), dimension(3,3) :: rotmat,rotmatx,rotmaty,rotmatz
    !! Combined and individual rotation matrices.
    real(real32), optional :: new_length
    !! If present, scale the rotated vector to this length.

    !  if(phi.ne.0._real32) phi=-phi

    rotmatx=reshape((/&
         1._real32,   0._real32,      0._real32,  &
         0._real32, cos(theta), -sin(theta),&
         0._real32, sin(theta),  cos(theta)/), shape(rotmatx))
    rotmaty=reshape((/&
         cos(phi), 0._real32, sin(phi),&
         0._real32,     1._real32,   0._real32,    &
         -sin(phi), 0._real32, cos(phi)/), shape(rotmaty))
    rotmatz=reshape((/&
         cos(psi), -sin(psi), 0._real32,&
         sin(psi), cos(psi), 0._real32,    &
         0._real32,        0._real32,       1._real32/), shape(rotmatz))


    rotmat=matmul(rotmaty,rotmatx)
    rotmat=matmul(rotmatz,rotmat)
    rotvec=matmul(a,transpose(rotmat))

    if(present(new_length))then
       magold=sqrt(dot_product(a,a))
       rotvec=rotvec*new_length/magold
    end if

    return
  end function rotvec
!###############################################################################


!###############################################################################
  function rot_arb_lat(a,lat,ang) result(vec)
    !! Rotate a 3D vector about arbitrary lattice axes.
    implicit none
    integer :: i
    !! Loop counter.
    real(real32), dimension(3) :: a,u,ang,vec
    !! Input vector, unit axis, rotation angles, and result.
    real(real32), dimension(3,3) :: rotmat,ident,lat
    !! Rotation matrix, identity matrix, and lattice matrix.


    ident=0._real32
    do i=1,3
       ident(i,i)=1._real32
    end do
   
    vec=a
    do i=1,3
       u=uvec(lat(i,:))
       rotmat=&
            (cos(ang(i))*ident)+&
            (sin(ang(i)))*cross_matrix(u)+&
            (1-cos(ang(i)))*outer_product(u,u)
       vec=matmul(vec,rotmat)
    end do


    return
  end function rot_arb_lat
!###############################################################################


!###############################################################################
!###############################################################################
!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!###############################################################################
!###############################################################################


!###############################################################################
  function gcd_num(numer,denom) result(gcd)
    !! Find the greatest common divisor of two integers.
    implicit none
    integer :: numer,denom
    !! Input numerator and denominator.
    integer :: a,b,c,gcd
    !! Working variables and GCD result.

    a=abs(numer)
    b=abs(denom)
    if(a.gt.b)then
       c=a
       a=b
       b=c
    end if

    if(a.eq.0)then
       gcd=b
       return
    end if

    do 
       c=mod(b,a)
       if(c.eq.0) exit
       b=a
       a=c
    end do
    gcd=a

    return
  end function gcd_num
!---------------------------------------------------------------------------
  function gcd_vec(vec) result(gcd)
    !! Find the greatest common divisor of an integer vector.
    implicit none
    integer :: i,a,b,c,dim,itmp1,loc
    !! Loop counter, working variables, dimension, temp, and location.
    integer :: gcd
    !! GCD result.
    integer, dimension(:),intent(in) :: vec
    !! Input integer vector.
    integer, allocatable, dimension(:) :: in_vec
    !! Working copy of input vector.


    dim=size(vec,dim=1)
    allocate(in_vec(dim))
    in_vec=abs(vec)
    do i=1,dim
       loc=maxloc(in_vec(i:dim),dim=1)+i-1
       itmp1=in_vec(i)
       in_vec(i)=in_vec(loc)
       in_vec(loc)=itmp1
    end do

    a=in_vec(2)
    do i=1,dim
       if(in_vec(i).eq.0) exit
       b=in_vec(i)
       do 
          c=mod(b,a)
          if(c.eq.0) exit
          b=a
          a=c
       end do
    end do
    gcd=a

    return
  end function gcd_vec
!###############################################################################


!###############################################################################
  function lcm(a,b)
    !! Find the lowest common multiple of two integers.
    implicit none
    integer :: a,b,lcm
    !! Input integers and LCM result.

    lcm=abs(a*b)/gcd(a,b)

    return
  end function lcm
!###############################################################################


!###############################################################################
  integer function get_frac_denom(val)
    !! Convert a decimal to a fraction and find the lowest denominator.
    implicit none
    integer :: i
    !! Iteration counter.
    real(real32) :: val
    !! Input decimal value.
    real(real32) :: a,b,c,tiny
    !! Working variables and tolerance.

    a=mod(val,1._real32)
    b=1._real32
    tiny = 1.E-6_real32
    i=0
    do 
       i=i+1
       if(abs(nint(1._real32/a)-(1._real32/a)).lt.tiny.and.&
            abs(nint(val*1._real32/a)-val*(1._real32/a)).lt.tiny) exit
       c=abs(b-a)
       b=a
       a=c
       if(i.ge.1000)then
          get_frac_denom=0
          return
       end if
    end do

    get_frac_denom=nint(1._real32/a)

    return
  end function get_frac_denom
!###############################################################################


!###############################################################################
  function reduce_vec_gcd(invec) result(vec)
    !! Reduce a real vector so that its GCD is 1.
    implicit none
    integer :: i,a
    !! Loop counter and integer GCD.
    real(real32) :: div,old_div,tol
    !! Divisor, previous divisor, and tolerance.
    real(real32), allocatable, dimension(:) :: vec,tvec
    !! Output vector and temporary vector.
    real(real32), dimension(:), intent(in) :: invec
    !! Input vector.


! MAKE IT DO SOMETHING IF IT CANNOT FULLY INTEGERISE

    tol=1.E-5_real32
    allocate(vec(size(invec)))
    vec=invec
    if(any(abs(vec(:)-nint(vec(:))).gt.tol))then
       div=abs(vec(1))
       do i=2,size(vec),1
          old_div=div
          if(min(abs(vec(i)),div).lt.tol)then
             div=max(abs(vec(i)),div)
             cycle
          end if
          div=abs(modulo(max(abs(vec(i)),div),min(abs(vec(i)),div)))
          if(abs(div).lt.tol) div=min(abs(vec(i)),old_div)
       end do
    else
       a=nint(vec(1))
       do i=2,size(vec)
          if(a.eq.0.and.int(vec(i)).eq.0) cycle
          a=gcd(a,int(vec(i)))
          if(abs(a).le.1)then
             a=1
             exit
          end if
       end do
       div=a
    end if

    if(div.eq.0._real32) return
    allocate(tvec(size(invec)))
    tvec=vec/div
    if(any(abs(tvec(:)-nint(tvec(:))).gt.tol)) return
    vec=tvec


  end function reduce_vec_gcd
!###############################################################################


!###############################################################################
  function gen_group(elem,mask,tol) result(group)
    !! Generate the entire group from a supplied subset of elements.
    implicit none
    integer :: i,j,k,nelem,ntot_elem,dim1,dim2,iter
    !! Loop counters, number of elements, total elements, dimensions, iteration.
    real(real32) :: tiny
    !! Tolerance for comparison.
    real(real32), allocatable, dimension(:,:) :: tmp_elem,cur_elem,apply_elem
    !! Temporary, current, and applied element matrices.
    real(real32), allocatable, dimension(:,:,:) :: tmp_group
    !! Temporary storage for group elements.

    real(real32), dimension(:,:,:), intent(in) :: elem
    !! Input subset of group elements.
    logical, dimension(:,:), optional, intent(in) :: mask
    !! Optional mask for wrapping elements.
    real(real32), allocatable, dimension(:,:,:) :: group
    !! Output full group.
    real(real32), optional, intent(in) :: tol
    !! Optional tolerance (default 1.E-5).


    if(present(tol))then
       tiny = tol
    else
       tiny = 1.E-5_real32
    end if
    nelem = size(elem(:,1,1))
    dim1 = size(elem(1,:,1))
    dim2 = size(elem(1,1,:))
    ! HARDCODED LIMIT OF A GROUP SIZE TO 10,000
    allocate(tmp_group(10000,dim1,dim2))
    allocate(tmp_elem(dim1,dim2))
    allocate(cur_elem(dim1,dim2))
    allocate(apply_elem(dim1,dim2))

    ntot_elem = 0
    elem_loop1: do i=1,nelem
       cur_elem(:,:) = elem(i,:,:)
       !write(0,*) "##########"
       !write(0,*)
       !write(0,*) i
       !write(0,'(2(2X,F9.6))') cur_elem(:,:)
       !write(0,*)
       if(present(mask))then
          where(mask.and.(cur_elem(:,:).lt.-tiny.or.cur_elem(:,:).ge.1._real32-tiny))
             cur_elem(:,:) = cur_elem(:,:) - floor(cur_elem(:,:)+tiny)
          end where
       end if
       do k=1,ntot_elem
          if(all(abs(tmp_group(k,:,:)-cur_elem(:,:)).lt.tiny)) cycle elem_loop1
       end do
       ntot_elem = ntot_elem + 1
       tmp_group(ntot_elem,:,:) = cur_elem(:,:)

       elem_loop2: do j=1,nelem
          tmp_elem(:,:) = cur_elem(:,:)
          apply_elem(:,:) = elem(j,:,:)
          iter = 0
          recursive_loop: do
             iter = iter + 1
             if(iter.ge.10)then
                write(0,'("ERROR: unending loop in mod_misc_linalg.f90")')
                write(0,'(2X,"subroutine gen_group in mod_misc_linalg.f90 encountered an unending loop")')
                write(0,'(2X,"Exiting...")')
                stop
             end if
             tmp_elem(:,:) = matmul((apply_elem(:,:)),tmp_elem(:,:))
             if(present(mask))then
                where(mask.and.(tmp_elem(:,:).lt.-tiny.or.tmp_elem(:,:).ge.1._real32-tiny))
                   tmp_elem(:,:) = tmp_elem(:,:) - floor(tmp_elem(:,:)+tiny)
                end where
             end if
             if(all(abs(cur_elem(:,:)-tmp_elem(:,:)).lt.tiny)) exit recursive_loop
             do k=1,ntot_elem
                if(all(abs(tmp_group(k,:,:)-tmp_elem(:,:)).lt.tiny)) cycle recursive_loop
             end do
             ntot_elem = ntot_elem + 1
             tmp_group(ntot_elem,:,:) = tmp_elem(:,:)
          end do recursive_loop
       end do elem_loop2
          
       
    end do elem_loop1

    allocate(group(ntot_elem,dim1,dim2))
    group(:,:,:) = tmp_group(:ntot_elem,:,:)
    return




  end function gen_group
!###############################################################################
  

end module artemis__misc_linalg
