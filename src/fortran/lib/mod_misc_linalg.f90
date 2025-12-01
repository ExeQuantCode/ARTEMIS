!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! module contains various linear algebra functions and subroutines.
!!! module includes the following functions and subroutines:
!!! uvec             (unit vector of vector of any size)
!!! modu             (magnitude of vector of any size)
!!! proj             (projection operator of one vector on another)
!!! GramSchmidt      (evaluates Gram-Schmidt orthogonal vectors)
!!! cross            (cross product of two vectors)
!!! cross_matrix     (generates cross product matrix of a vector)
!!! outer_product    (performs outer_product of two vectors)
!!! vec_mat_mul      (multiply a vector with a matrix)
!!! get_vec_multiple (determines the scaling factor between two vectors)
!!!##################
!!! get_angle        (get the angle between two vectors)
!!! get_area         (get the area made by two vectors)
!!! get_vol          (get the volume of a matrix)
!!! trace            (trace of a matrix of any size)
!!! det              (determinant of a 3x3 matrix)
!!! inverse          (inverse of a 3x3 matrix)
!!! rec_det          (determinant of a matrix of any size)
!!! LUdet            (determinant of a matrix of any size using LUdecomposition)
!!! LUinv            (inverse of a matrix of any size using LUdecomposition)
!!! LUdecompose      (decompose a matrix into upper and lower matrices. A=LU)
!!!##################
!!! find_tf          (transformation matrix to move between two matrices)
!!! simeq            (simultaneous equation solver)
!!! LLL_reduce       (performs LLL reduction on a basis)
!!! rotvec           (rotate vector in 3D space about x, y, z cartesian axes)
!!! rot_arb_lat      (rotate vector in 3D space about a, b, c arbitrary axes)
!!!##################
!!! gcd              (greatest common denominator (to reduce a fraction))
!!! lcm              (lowest common multiple)
!!! get_frac_denom   (convert decimal to fraction and finds lowest denominator)
!!! reduce_vec_gcd   (reduces the gcd of a vector to 1)
!!!##################
!!! gen_group        (generate group from a subset of elements)
!!!#############################################################################
module misc_linalg
  use artemis__constants, only: real32
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



!!!updated 2021/12/09


contains
!!!#####################################################
!!! finds unit vector of an arbitrary vector
!!!#####################################################
  function uvec(vec) result(output)
    implicit none
    real(real32),dimension(:)::vec
    real(real32),allocatable,dimension(:) :: output
    allocate(output(size(vec)))
    output = vec/modu(vec)
  end function uvec
!!!#####################################################


!!!#####################################################
!!! finds modulus of an arbitrary length vector
!!!#####################################################
  function modu(vec) result(output)
    implicit none
    real(real32),dimension(:)::vec
    real(real32)::output
    output = abs(sqrt(sum(vec(:)**2)))
  end function modu
!!!#####################################################


!!!#####################################################
!!! projection operator
!!!#####################################################
!!! projection of v on u
  function proj(u,v) result(output)
    implicit none
    real(real32), dimension(:) :: u,v
    real(real32), allocatable, dimension(:) :: output

    allocate(output(size(u,dim=1)))
    output = u*dot_product(v,u)/dot_product(u,u)

  end function proj
!!!#####################################################


!!!#####################################################
!!! Gram-Schmidt process
!!!#####################################################
!!! assumes basis(n,m) is a basis of n vectors, each ...
!!! ... of m-dimensions
!!! rmc = row major order
  function GramSchmidt(basis,normalise,cmo) result(u)
    implicit none
    integer :: num,dim,i,j
    real(real32), allocatable, dimension(:) :: vtmp
    real(real32), dimension(:,:), intent(in) :: basis
    real(real32), allocatable, dimension(:,:) :: u
    logical, optional, intent(in) :: cmo
    logical, optional, intent(in) :: normalise


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
             u(i,:) = u(i,:)/modu(u(i,:))
          end do
       end if
    end if


  end function GramSchmidt
!!!#####################################################


!!!#####################################################
!!! cross product
!!!#####################################################
  pure function cross(a,b) result(output)
    implicit none
    real(real32), dimension(3) :: output
    real(real32), dimension(3), intent(in) :: a,b

    output(1) = a(2)*b(3) - a(3)*b(2)
    output(2) = a(3)*b(1) - a(1)*b(3)
    output(3) = a(1)*b(2) - a(2)*b(1)

  end function cross
!!!#####################################################


!!!#####################################################
!!! cross product matrix
!!!#####################################################
!!! a = (a1,a2,a3)
!!! 
!!!         (  0  -a3  a2 )
!!! [a]_x = (  a3  0  -a1 )
!!!         ( -a2  a1  0  )
!!!#####################################################
  function cross_matrix(a)
    implicit none
    real(real32), dimension(3,3) :: cross_matrix
    real(real32), dimension(3), intent(in) :: a

    cross_matrix=0._real32

    cross_matrix(1,2) = -a(3)
    cross_matrix(1,3) =  a(2)
    cross_matrix(2,3) = -a(1)

    cross_matrix(2,1) =  a(3)
    cross_matrix(3,1) = -a(2)
    cross_matrix(3,2) =  a(1)

    return
  end function cross_matrix
!!!#####################################################


!!!#####################################################
!!! outer product
!!!#####################################################
  function outer_product(a,b)
    implicit none
    integer :: j
    real(real32), dimension(:) :: a,b
    real(real32),allocatable,dimension(:,:)::outer_product
   
    allocate(outer_product(size(a),size(b)))

    do j=1,size(b)
       outer_product(:,j)=a(:)*b(j)
    end do

    return
  end function outer_product
!!!#####################################################


!!!#####################################################
!!! function to multiply a vector and a matrix
!!!#####################################################
  function ivec_dmat_mul(a,mat) result(vec)
    implicit none
    integer :: j
    integer, dimension(:) :: a
    real(real32), dimension(:,:) :: mat
    real(real32),allocatable,dimension(:) :: vec

    vec=0._real32
    allocate(vec(size(a)))
    do j=1,size(a)
       vec(:)=vec(:)+real(a(j),real32)*mat(j,:)
    end do

    return
  end function ivec_dmat_mul
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  function dvec_dmat_mul(a,mat) result(vec)
    implicit none
    integer :: j
    real(real32), dimension(:) :: a
    real(real32), dimension(:,:) :: mat
    real(real32),allocatable,dimension(:) :: vec

    vec=0._real32
    allocate(vec(size(a)))
    do j=1,size(a)
       vec(:)=vec(:)+a(j)*mat(j,:)
    end do

    return
  end function dvec_dmat_mul
!!!#####################################################


!!!#####################################################
!!! get vec_multiple
!!!#####################################################
  function get_vec_multiple(a,b) result(multi)
    implicit none
    integer :: i
    real(real32) :: multi
    real(real32), dimension(:) :: a,b
    
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
!!!#####################################################


!!!#############################################################################
!!!#############################################################################
!!!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!!!#############################################################################
!!!#############################################################################


!!!#####################################################
!!! returns angle between two vectors
!!!#####################################################
  function get_angle(vec1,vec2) result(angle)
    implicit none
    real(real32) :: angle
    real(real32), dimension(3) :: vec1,vec2

    angle = acos( dot_product(vec1,vec2)/&
         ( modu(vec1) * modu(vec2) ))
    if (isnan(angle)) angle = 0._real32

    return
  end function get_angle
!!!#####################################################


!!!#####################################################
!!! returns area made by two vectors
!!!#####################################################
  function get_area(a,b) result(area)
    implicit none
    real(real32) :: area
    real(real32), dimension(3) :: vec,a,b

    vec = cross(a,b)
    area = sqrt(dot_product(vec,vec))

    return
  end function get_area
!!!#####################################################


!!!#####################################################
!!! returns volume of a lattice
!!!#####################################################
  function get_vol(lat) result(vol)
    implicit none
    integer :: n,i,j,k,l
    real(real32) :: vol,scale
    real(real32), dimension(3,3) :: lat
    real(real32), dimension(3) :: a,b,c

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
!!!#####################################################


!!!#####################################################
!!! finds trace of an arbitrary dimension square matrix
!!!#####################################################
  function trace(mat) result(output)
    integer::j
    real(real32), dimension(:,:), intent(in) :: mat
    real(real32) :: output
    output = 0._real32
    do j = 1, size(mat,1)
      output = output + mat(j,j)
    end do
  end function trace
!!!#####################################################


!!!#####################################################
!!! returns determinant of 3 x 3 matrix
!!!#####################################################
  function idet(mat) result(output)
    integer :: output
    integer, dimension(3,3), intent(in) :: mat

    output = mat(1,1)*mat(2,2)*mat(3,3)-mat(1,1)*mat(2,3)*mat(3,2)&
         - mat(1,2)*mat(2,1)*mat(3,3)+mat(1,2)*mat(2,3)*mat(3,1)&
         + mat(1,3)*mat(2,1)*mat(3,2)-mat(1,3)*mat(2,2)*mat(3,1)

  end function idet
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  function ddet(mat) result(output)
    real(real32) :: output
    real(real32), dimension(3,3), intent(in) :: mat

    output = mat(1,1)*mat(2,2)*mat(3,3)-mat(1,1)*mat(2,3)*mat(3,2)&
         - mat(1,2)*mat(2,1)*mat(3,3)+mat(1,2)*mat(2,3)*mat(3,1)&
         + mat(1,3)*mat(2,1)*mat(3,2)-mat(1,3)*mat(2,2)*mat(3,1)

  end function ddet
!!!#####################################################


!!!#####################################################
!!! returns inverse of 2x2 or 3x3 matrix
!!!#####################################################
  pure function inverse(mat)
    real(real32), dimension(:,:), intent(in) :: mat
    real(real32), dimension(size(mat,dim=1),size(mat,dim=2)) :: inverse

    select case(size(mat,dim=2))
    case(2)
       inverse = inverse_2x2(mat)
    case(3)
       inverse = inverse_3x3(mat)
    end select

  end function inverse
!!!#####################################################


!!!#####################################################
!!! returns inverse of 2 x 2 matrix
!!!#####################################################
  pure function inverse_2x2(mat) result(output)
    implicit none
    real(real32), dimension(2,2), intent(in) :: mat
    real(real32), dimension(2,2) :: output
    real(real32) :: inv_det

    associate(a => mat(1,1), b => mat(1,2), c => mat(2,1), d => mat(2,2))
       inv_det = 1._real32 / (a * d - b * c)

       output(1,1) =  d * inv_det
       output(1,2) = -b * inv_det
       output(2,1) = -c * inv_det
       output(2,2) =  a * inv_det
    end associate

  end function inverse_2x2
!!!#####################################################


!!!#####################################################
!!! returns inverse of 3 x 3 matrix
!!!#####################################################
  pure function inverse_3x3(mat) result(output)
  implicit none
  real(real32), dimension(3,3), intent(in) :: mat
  real(real32), dimension(3,3) :: output
  real(real32) :: inv_det
  real(real32) :: c00, c01, c02, c10, c11, c12, c20, c21, c22

  associate( &
    m11 => mat(1,1), m12 => mat(1,2), m13 => mat(1,3), &
    m21 => mat(2,1), m22 => mat(2,2), m23 => mat(2,3), &
    m31 => mat(3,1), m32 => mat(3,2), m33 => mat(3,3))

    ! Cofactors
    c00 =  m22 * m33 - m23 * m32
    c01 = -m21 * m33 + m23 * m31
    c02 =  m21 * m32 - m22 * m31

    c10 = -m12 * m33 + m13 * m32
    c11 =  m11 * m33 - m13 * m31
    c12 = -m11 * m32 + m12 * m31

    c20 =  m12 * m23 - m13 * m22
    c21 = -m11 * m23 + m13 * m21
    c22 =  m11 * m22 - m12 * m21

    inv_det = 1._real32 / (m11 * c00 + m12 * c01 + m13 * c02)

    ! Transpose cofactors into the inverse
    output(1,1) = c00 * inv_det
    output(2,1) = c01 * inv_det
    output(3,1) = c02 * inv_det

    output(1,2) = c10 * inv_det
    output(2,2) = c11 * inv_det
    output(3,2) = c12 * inv_det

    output(1,3) = c20 * inv_det
    output(2,3) = c21 * inv_det
    output(3,3) = c22 * inv_det

  end associate
end function inverse_3x3
!!!#####################################################


!!!#####################################################
!!! determinant function
!!!#####################################################
  recursive function rec_det(a,n) result(res)
    integer :: i, sign
    real(real32) :: res
    integer, intent(in) :: n
    real(real32), dimension(n,n), intent(in) :: a
    real(real32), dimension(n-1, n-1) :: tmp

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
!!!#####################################################


!!!#####################################################
!!! determinant of input matrix via LU decomposition
!!!#####################################################
!!! L = lower
!!! U = upper
!!! inmat = input nxn matrix
!!! LUdet = determinant of inmat
!!! LUdet = (-1)**N * prod(L(i,i)*U(i,i))
  function LUdet(inmat)
    implicit none
    integer :: i,N
    real(real32) :: LUdet
    real(real32), dimension(:,:) :: inmat
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: L,U

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
!!!#####################################################


!!!#####################################################
!!! inverse of n x n matrix
!!!#####################################################
!!! doesn't work if a diagonal element = 0
!!! L = lower
!!! U = upper
!!! inmat = input nxn matrix
!!! LUinv = output nxn inverse of matrix
!!! Lz=b
!!! Ux=z
!!! x=column vectors of the inverse matrix
  function LUinv(inmat)
    implicit none
    integer :: i,m,N
    real(real32), dimension(:,:) :: inmat
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: LUinv
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: L,U
    real(real32), dimension(size(inmat,1)) :: c,z,x

    L=0._real32
    U=0._real32
    N=size(inmat,1)
    call LUdecompose(inmat,L,U)

!!! Lz=c
!!! c are column vectors of the identity matrix
!!! uses forward substitution to solve
    do m=1,N
       c=0._real32
       c(m)=1._real32

       z(1)=c(1)
       do i=2,N
          z(i)=c(i)-dot_product(L(i,1:i-1),z(1:i-1))
       end do


!!! Ux=z
!!! x are the rows of the inversion matrix
!!! uses backwards substitution to solve
       x(N)=z(N)/U(N,N)
       do i=N-1,1,-1
          x(i)=z(i)-dot_product(U(i,i+1:N),x(i+1:N))
          x(i)= x(i)/U(i,i)
       end do

       LUinv(:,m)=x(:)
    end do

    return
  end function LUinv
!!!#####################################################


!!!#####################################################
!!! A=LU matrix decomposer
!!!#####################################################
!!! Method: Based on Doolittle LU factorization for Ax=b
!!! doesn't work if a diagonal element = 0
!!! L = lower
!!! U = upper
!!! inmat = input nxn matrix
  subroutine LUdecompose(inmat,L,U)
    implicit none
    integer :: i,j,N
    real(real32), dimension(:,:) :: inmat,L,U
    real(real32), dimension(size(inmat,1),size(inmat,1)) :: mat

    N=size(inmat,1)
    mat=inmat
    L=0._real32
    U=0._real32

    do j=1,N
       L(j,j)=1._real32
    end do
!!! Solves the lower matrix
    do j=1,N-1
       do i=j+1,N
          L(i,j)=mat(i,j)/mat(j,j)
          mat(i,j+1:N)=mat(i,j+1:N)-L(i,j)*mat(j,j+1:N)
       end do
    end do

!!! Equates upper half of remaining mat to upper matrix
    do j=1,N
       do i=1,j
          U(i,j)=mat(i,j)
       end do
    end do

    return
  end subroutine LUdecompose
!!!#####################################################


!!!#############################################################################
!!!#############################################################################
!!!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!!!#############################################################################
!!!#############################################################################


!!!#####################################################
!!! find transformation matrix between two matrices
!!! A=mat1; B=mat2; T=find_tf
!!! A T = B
!!! A^-1 A T = A^-1 B
!!! T = A^-1 B
!!!#####################################################
  function find_tf(mat1,mat2) result(tf)
    implicit none
    real(real32), dimension(:,:) :: mat1,mat2
    real(real32), dimension(size(mat1,dim=1),size(mat1,dim=2)) :: tf

    tf=matmul(inverse(mat1),mat2)

  end function find_tf
  function find_tf_2x2(mat1,mat2) result(tf)
    implicit none
    real(real32), dimension(2,2) :: mat1,mat2
    real(real32), dimension(2,2) :: tf

    tf=matmul(inverse_2x2(mat1),mat2)

  end function find_tf_2x2
!!!#####################################################


!!!#####################################################
!!! simultaneous equation solver for n dimensions
!!!#####################################################
!!! P     = power seiers equation in matrix
!!! invP  = inverse of the power seiers matrix
!!! qX    = the x values of the power seires with a ...
!!!         ... size equal to order
!!! qY    = the Y values for the n simult eqns
!!! simeq = the coefficients of the powers of ...
!!!                ... qX with highest power simeq(1)
!!! f(qX)=qY
!!! qA P(qX) = qY (in matrix form)
!!! hence, qA=qY P^-1
  function simeq(qX,qY)
    integer :: i,j,n,loc
    real(real32), dimension(:) :: qX,qY
    real(real32), dimension(size(qY)) :: funcY
    real(real32), dimension(size(qY)) :: simeq,tmpqY
    real(real32), dimension(size(qY),size(qY)) :: P,invP,tmpP


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
!!!#####################################################


!!!#####################################################
!!! Lenstra-Lenstra-Lovász reduction
!!!#####################################################
!!! LLL algorithm based on the one found on Wikipedia, ...
!!! ... which is based on Hoffstein, Pipher and Silverman 2008
!!! https://en.wikipedia.org/wiki/Lenstra–Lenstra–Lovász_lattice_basis_reduction_algorithm
  function LLL_reduce(basis,delta) result(obas)
    implicit none
    integer :: num,dim,i,j,k,loc
    real(real32) :: d,dtmp
    real(real32), allocatable, dimension(:) :: vtmp,mag_bas
    real(real32), allocatable, dimension(:,:) :: mu,GSbas,obas

    real(real32), dimension(:,:), intent(in) :: basis
    real(real32), optional, intent(in) :: delta


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
       mag_bas(i) = modu(obas(i,:))
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


!!! Separate functions for this to run efficiently
  contains
    !!function to get the mu values
    function get_mu(bas1,bas2) result(mu)
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


    !!subroutine to update Gram-Schmidt vectors and mu values
    subroutine update_GS_and_mu(GSbas,mu,basis,k)
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
!!!#####################################################


!!!#####################################################
!!! vector rotation
!!!#####################################################
  function rotvec(a,theta,phi,psi,new_length)
    implicit none
    real(real32) :: magold,theta,phi,psi
    real(real32), dimension(3) :: a,rotvec
    real(real32), dimension(3,3) :: rotmat,rotmatx,rotmaty,rotmatz
    real(real32), optional :: new_length

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
!!!#####################################################


!!!#####################################################
!!! vector rotation
!!!#####################################################
  function rot_arb_lat(a,lat,ang) result(vec)
    implicit none
    integer :: i
    real(real32), dimension(3) :: a,u,ang,vec
    real(real32), dimension(3,3) :: rotmat,ident,lat


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
!!!#####################################################


!!!#############################################################################
!!!#############################################################################
!!!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!!!#############################################################################
!!!#############################################################################


!!!#####################################################
!!! finds the greatest common denominator
!!!#####################################################
  function gcd_num(numer,denom) result(gcd)
    implicit none
    integer :: numer,denom
    integer :: a,b,c,gcd

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
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  function gcd_vec(vec) result(gcd)
    implicit none
    integer :: i,a,b,c,dim,itmp1,loc
    integer :: gcd
    integer, dimension(:),intent(in) :: vec
    integer, allocatable, dimension(:) :: in_vec


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
!!!#####################################################


!!!#####################################################
!!! finds the lowest common multiple
!!!#####################################################
  function lcm(a,b)
    implicit none
    integer :: a,b,lcm

    lcm=abs(a*b)/gcd(a,b)

    return
  end function lcm
!!!#####################################################


!!!#####################################################
!!! converts decimal into a fraction and finds the ...
!!! ... lowest denominator for it.
!!!#####################################################
  integer function get_frac_denom(val)
    implicit none
    integer :: i
    real(real32) :: val
    real(real32) :: a,b,c,tiny

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
!!!#####################################################


!!!#####################################################
!!! reduces the gcd of a vector to 1
!!!#####################################################
  function reduce_vec_gcd(invec) result(vec)
    implicit none
    integer :: i,a
    real(real32) :: div,old_div,tol
    real(real32), allocatable, dimension(:) :: vec,tvec
    real(real32), dimension(:), intent(in) :: invec


!!! MAKE IT DO SOMETHING IF IT CANNOT FULLY INTEGERISE

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
!!!#####################################################


!!!#####################################################
!!! generate entire group from supplied elements
!!!#####################################################
  function gen_group(elem,mask,tol) result(group)
    implicit none
    integer :: i,j,k,nelem,ntot_elem,dim1,dim2,iter
    real(real32) :: tiny
    real(real32), allocatable, dimension(:,:) :: tmp_elem,cur_elem,apply_elem
    real(real32), allocatable, dimension(:,:,:) :: tmp_group

    real(real32), dimension(:,:,:), intent(in) :: elem
    logical, dimension(:,:), optional, intent(in) :: mask
    real(real32), allocatable, dimension(:,:,:) :: group
    real(real32), optional, intent(in) :: tol


    if(present(tol))then
       tiny = tol
    else
       tiny = 1.E-5_real32
    end if
    nelem = size(elem(:,1,1))
    dim1 = size(elem(1,:,1))
    dim2 = size(elem(1,1,:))
    !!! HARDCODED LIMIT OF A GROUP SIZE TO 10,000
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
!!!#####################################################
  

end module misc_linalg
