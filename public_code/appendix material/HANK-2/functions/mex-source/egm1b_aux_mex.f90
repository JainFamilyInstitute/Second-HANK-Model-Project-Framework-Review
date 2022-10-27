#include "fintrf.h"

! [c,m]=EGM1b_AUX_MEX(grid_m,m_n,c_n)
!
! egm1b_aux_mex: performs the fortran analogue to
! for j=1:mpar.nk*mpar.ns*mpar.nh*mpar.nM*mpar.nH %loop over stochastic states, Parallel is not worth the overhea
!       [ c_update(:,j), m_update(:,j)] = myinter1m_mex(m_star_n(:,j),grid.m,c_n_aux(:,j),grid.m);
! end
!
!----------------------------------------------------------------------------
! Gateway routine
!----------------------------------------------------------------------------

subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    implicit none
    
    ! Arguments in mexFunction
      mwPointer  :: plhs(*), prhs(*)
      integer    :: nlhs, nrhs
    
    ! Function declarations
    mwPointer   ::  mxGetPr, mxCreateDoubleMatrix
    mwPointer     ::  mxGetN , mxGetM

    ! Pointers to input/output mxArrays:
    ! Input
      mwPointer  :: grid_m, m_n, c_n
    ! Output
      mwPointer  :: c, m
      
 
    ! Array size information
      mwSize    :: nrows,ncols            ! length of input array

    
!--------------------------------------------------------------------------------    


    grid_m = mxGetPr(prhs(1))
    m_n = mxGetPr(prhs(2))
    c_n = mxGetPr(prhs(3))



! Retrieve size of x
    nrows = mxGetM(prhs(2))
    ncols = mxGetN(prhs(2))

    plhs(1)=mxCreateDoubleMatrix(nrows, ncols, int(0,4)) ! Reserve memory
    plhs(2)=mxCreateDoubleMatrix(nrows, ncols, int(0,4))

    c = mxGetPr(plhs(1)) !copy pointers
    m = mxGetPr(plhs(2))

    call myinter_inp(nrows, ncols, %val(grid_m), %val(m_n), %val(c_n), %val(m), %val(c))
	
return
end subroutine



subroutine myinter_inp(nrows,ncols,grid_m,m_n,c_n,out_m,out_c)
    implicit none
    ! Input arguments
    mwSize, intent(in)      :: nrows                             ! number of points in x-vector
    mwSize, intent(in)      :: ncols		                      ! size of xi-array
    real(8), intent(in)         :: grid_m(nrows)                          ! x-values
    real(8), intent(in)         :: m_n(nrows,ncols)    ! interpolation points
    real(8), intent(in)         :: c_n(nrows,ncols)   ! f1(x)-values

    ! Output arguments
    real(8), intent(out)        :: out_m(nrows,ncols)  ! interpolations on f1
    real(8), intent(out)        :: out_c(nrows,ncols)  ! interpolations on f2
    ! Local variables
    integer                     :: i1,i2
    integer						:: ind			                  ! index of nearest neighbor to left
	real(8)						:: weight		                  ! interpolation weight on left neighbor



	! Loop over xi and interpolate
	do i1 = 1,ncols
        do i2 =1,nrows
				!ind = max(min(int((xi(i1,i2,i3)-x(1))/incr)+1,nx-1),1)
                call locate(m_n(:,i1),nrows,grid_m(i2),ind)
				ind = max(min(ind,nrows-1),1)
				weight = 1.d0-(grid_m(i2)-m_n(ind,i1))/(m_n(ind+1,i1)-m_n(ind,i1))
				out_m(i2,i1) = weight * grid_m(ind) + (1.d0-weight) * grid_m(ind+1)
				out_c(i2,i1) = weight * c_n(ind,i1) + (1.d0-weight) * c_n(ind+1,i1)
        end do
	end do

    return
end subroutine

subroutine locate(xx,N,x,j)
	implicit none
	real(8), intent(in), dimension(N)		:: xx  ! lookup table
	mwSize, intent(in)						:: N   ! no elements of lookup table
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

