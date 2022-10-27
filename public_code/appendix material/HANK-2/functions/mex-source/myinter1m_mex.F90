#include "fintrf.h"

!----------------------------------------------------------------------------
! Gateway routine
!----------------------------------------------------------------------------

subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    implicit none
    
    ! Arguments in mexFunction
      mwPointer  :: plhs(*), prhs(*)
      integer    :: nlhs, nrhs
    
    ! Function declarations
      mwPointer     mxCreateNumericArray ,mxGetPr
      mwPointer     mxGetDimensions
      integer(4)    mxClassIDFromClassName
      integer(4)    mxGetString
      mwSize        mxGetNumberOfDimensions
      
    ! Pointers to input/output mxArrays:
      ! Input
      mwPointer  :: x_pr, xi_pr, indicator_pr, y1_pr, y2_pr, y3_pr
       ! Output
      mwPointer  :: yi1_pr, yi2_pr, yi3_pr
      
 
    ! Array size information
      integer(4)   :: classid
      integer(4)   :: complexflag
      mwSize        :: nx                ! length of input array
	  mwSize        :: nxi
      mwSize        :: ndim_xi           ! number of dimensions of xi array
      mwSize        :: dims_x(2)         ! size array for x
      mwSize        :: dims_xi(2)       ! size array for xi
    
!--------------------------------------------------------------------------------    

    ! Retrieve size of x 
      call mxCopyPtrToPtrArray(mxGetDimensions(prhs(1)), dims_x, mxGetNumberOfDimensions(prhs(1)))
      nx = maxval(dims_x)
      
    ! Retrieve no of dimensions of xi
      ndim_xi = mxGetNumberOfDimensions(prhs(2))
    ! Retrieve size of xi
      call mxCopyPtrToPtrArray(mxGetDimensions(prhs(2)), dims_xi, mxGetNumberOfDimensions(prhs(2)))
	  nxi=maxval(dims_xi)
    
    ! Copy pointers to double precision input arrays
      x_pr = mxGetPr(prhs(1))
      xi_pr = mxGetPr(prhs(2))
      y1_pr = mxGetPr(prhs(3))
      if (nrhs .eq. 4) then
        y2_pr = mxGetPr(prhs(4))
      elseif (nrhs .eq. 5) then
        y2_pr = mxGetPr(prhs(4))
        y3_pr = mxGetPr(prhs(5))
      end if
    
            
      ! Create matrices for return arguments
      classid = mxClassIDFromClassName('double')
      complexflag = 0  
      
      plhs(1) = mxCreateNumericArray(ndim_xi, dims_xi, classid, complexflag)
      yi1_pr = mxGetPr(plhs(1))
      if (nrhs .eq. 4) then
        plhs(2) = mxCreateNumericArray(ndim_xi, dims_xi, classid, complexflag)
        yi2_pr = mxGetPr(plhs(2))
      elseif (nrhs .eq.  5) then
        plhs(2) = mxCreateNumericArray(ndim_xi, dims_xi, classid, complexflag)
        yi2_pr = mxGetPr(plhs(2))
        plhs(3) = mxCreateNumericArray(ndim_xi, dims_xi, classid, complexflag)
        yi3_pr = mxGetPr(plhs(3))
      end if
    
	if (nrhs .eq. 4) then
		call myinter3_2inp(nx,nxi,%val(x_pr),%val(xi_pr), &
		%val(y1_pr),%val(y2_pr),%val(yi1_pr),%val(yi2_pr))
	elseif (nrhs .eq. 3) then
		call myinter3_1inp(nx,nxi,%val(x_pr),%val(xi_pr),%val(y1_pr),%val(yi1_pr))
	elseif (nrhs .eq. 5) then  
		call myinter3_3inp(nx,nxi,%val(x_pr),%val(xi_pr), &
		%val(y1_pr),%val(y2_pr),%val(y3_pr),%val(yi1_pr),%val(yi2_pr),%val(yi3_pr))
	end if
   
      return
end subroutine

subroutine myinter3_1inp(nx,nxi,x,xi,fx1,fxi1)
    implicit none
    ! Input arguments
    mwSize, intent(in)      :: nx                             ! number of points in x-vector
    mwSize, intent(in)      :: nxi		                       ! size of xi-array
    real(8), intent(in)         :: x(nx)                          ! x-values
    real(8), intent(in)         :: xi(nxi)    ! interpolation points
    real(8), intent(in)         :: fx1(nx)   ! f1(x)-values
    ! Output arguments
    real(8), intent(out)        :: fxi1(nxi)  ! interpolations on f1
    
    ! Local variables
    integer                     :: i1
    integer						:: ind			                  ! index of nearest neighbor to left
	real(8)						:: weight		                  ! interpolation weight on left neighbor

  
	! Loop over xi and interpolate
	do i1 = 1,nxi
		!ind = max(min(int((xi(i1,i2,i3)-x(1))/incr)+1,nx-1),1)
		call locate(x,nx,xi(i1),ind)
		ind = max(min(ind,nx-1),1)
		weight = 1.d0-(xi(i1)-x(ind))/(x(ind+1)-x(ind))
		fxi1(i1) = weight * fx1(ind) + (1-weight) * fx1(ind+1)
	end do

    return
end subroutine

subroutine myinter3_3inp(nx,nxi,x,xi,fx1,fx2,fx3,fxi1,fxi2,fxi3)
    implicit none
    ! Input arguments
    mwSize, intent(in)      :: nx                             ! number of points in x-vector
    mwSize, intent(in)      :: nxi		                       ! size of xi-array
    real(8), intent(in)         :: x(nx)                          ! x-values
    real(8), intent(in)         :: xi(nxi)    ! interpolation points
    real(8), intent(in)         :: fx1(nx)   ! f1(x)-values
    real(8), intent(in)         :: fx2(nx)   ! f2(x)-values
    real(8), intent(in)         :: fx3(nx)   ! f2(x)-values
    ! Output arguments
    real(8), intent(out)        :: fxi1(nxi)  ! interpolations on f1
    real(8), intent(out)        :: fxi2(nxi)  ! interpolations on f2
    real(8), intent(out)        :: fxi3(nxi)  ! interpolations on f2
    ! Local variables
    integer                     :: i1
    integer						:: ind			                  ! index of nearest neighbor to left
	real(8)						:: weight		                  ! interpolation weight on left neighbor
    
  
	! Loop over xi and interpolate
	do i1 = 1,nxi
		!ind = max(min(int((xi(i1,i2,i3)-x(1))/incr)+1,nx-1),1)
		call locate(x,nx,xi(i1),ind)
		ind = max(min(ind,nx-1),1)
		weight = 1.d0-(xi(i1)-x(ind))/(x(ind+1)-x(ind))
		fxi1(i1) = weight * fx1(ind) + (1-weight) * fx1(ind+1)
		fxi2(i1) = weight * fx2(ind) + (1-weight) * fx2(ind+1)
		fxi3(i1) = weight * fx3(ind) + (1-weight) * fx3(ind+1)
	end do

    return
end subroutine

subroutine myinter3_2inp(nx,nxi,x,xi,fx1,fx2,fxi1,fxi2)
    implicit none
    ! Input arguments
    mwSize, intent(in)      :: nx                             ! number of points in x-vector
    mwSize, intent(in)      :: nxi		                       ! size of xi-array
    real(8), intent(in)         :: x(nx)                          ! x-values
    real(8), intent(in)         :: xi(nxi)    ! interpolation points
    real(8), intent(in)         :: fx1(nx)   ! f1(x)-values
    real(8), intent(in)         :: fx2(nx)   ! f2(x)-values
    ! Output arguments
    real(8), intent(out)        :: fxi1(nxi)  ! interpolations on f1
    real(8), intent(out)        :: fxi2(nxi)  ! interpolations on f2
    ! Local variables
    integer                     :: i1
    integer						:: ind			                  ! index of nearest neighbor to left
	real(8)						:: weight		                  ! interpolation weight on left neighbor
    

	! Loop over xi and interpolate
	do i1 = 1,nxi
				!ind = max(min(int((xi(i1,i2,i3)-x(1))/incr)+1,nx-1),1)
				call locate(x,nx,xi(i1),ind)
				ind = max(min(ind,nx-1),1)
				weight = 1.d0-(xi(i1)-x(ind))/(x(ind+1)-x(ind))
				fxi1(i1) = weight * fx1(ind) + (1-weight) * fx1(ind+1)
				fxi2(i1) = weight * fx2(ind) + (1-weight) * fx2(ind+1)
	end do

    return
end subroutine

subroutine locate(xx,N,x,j)
	implicit none
	real(8), intent(in), dimension(N)		:: xx  ! lookup table
	integer, intent(in)						:: N   ! no elements of lookup table
	real(8), intent(in)						:: x   ! Value whose neares neighbors are to be found
	integer, intent(out)					:: j   ! returns j if xx(j)<x<xx(j+1),
												   ! 0 if value is to the left of grid,
												   ! N if value is to the right of grid
	! locals
	integer :: jl, jm, ju
	
	jl = 0
	ju = N+1
	
	do
		if ((ju-jl)==1) exit
		jm = (ju+jl)/2
		if ((xx(N) > xx(1)) .eqv. (x > xx(jm))) then 
			jl = jm
		else
			ju = jm
		end if
	end do											   
	
	! Set output
	if (x == xx(1)) then
		j = 1
	elseif (x == xx(N)) then
		j = N
	else
		j = jl
	end if
end subroutine locate