!
!   myinterp3_mex.F90
!   
!
!   Created by Volker Tjaden on 23.08.2013.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
#include "fintrf.h"

!-----------------------------------------------------------------------
! Gateway routine
!-----------------------------------------------------------------------

subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    	implicit none
	! Arguments in mexFunction
    	mwPointer  :: plhs(*), prhs(*)
    	integer    :: nlhs, nrhs
	
    	! Function declarations
    	mwPointer   ::  mxGetPr, mxCreateNumericArray
    	mwPointer   ::  mxGetDimensions
	mwSize      ::  mxGetNumberOfDimensions
    	integer(4)  ::  mxClassIDFromClassName
    	integer(4)  ::  mxGetString
    	mwPointer   ::  mxGetN , mxGetM

    	! Pointers to input/output mxArrays:
    	! Input
    	mwPointer  :: x1,x2,x3,x4,fx,x1i,x2i,x3i,x4i
    	! Output
    	mwPointer  :: fint
	
    	! Array size information
    	mwSize	  :: ndim
    	mwSize	  :: dims(4)
    	mwSize	  :: nx1i,nx2i,nx3i,nx4i
    	integer(4)   :: classid
    	integer(4)   :: complexflag

!-----------------------------------------------------------------------

	! Copy pointers
	x1	= mxGetPr(prhs(1))
	x2	= mxGetPr(prhs(2))
	x3	= mxGetPr(prhs(3))
	x4	= mxGetPr(prhs(4))
	fx	= mxGetPr(prhs(5))
	x1i  	= mxGetPr(prhs(6))
	x2i  	= mxGetPr(prhs(7))
	x3i  	= mxGetPr(prhs(8))
	x4i  	= mxGetPr(prhs(9))
	
	! Retrieve size information
	ndim= mxGetNumberOfDimensions(prhs(5))
	call mxCopyPtrToPtrArray(mxGetDimensions(prhs(5)),&
	dims, ndim)
	nx1i=max(mxGetN(prhs(6)),mxGetM(prhs(6)))
	nx2i=max(mxGetN(prhs(7)),mxGetM(prhs(7)))
	nx3i=max(mxGetN(prhs(8)),mxGetM(prhs(8)))
	nx4i=max(mxGetN(prhs(9)),mxGetM(prhs(9)))
	
	! Create return arguments
	classid = mxClassIDFromClassName('double')
    	complexflag = 0  
      
    	plhs(1) = mxCreateNumericArray(ndim,(/nx1i,nx2i,nx3i,nx4i/), classid, complexflag)
    	fint = mxGetPr(plhs(1))

	! call subroutine
	call myinterp4(%val(fint),%val(fx),%val(x1),%val(x2),%val(x3),%val(x4),&
		%val(x1i),%val(x2i),%val(x3i),%val(x4i),dims(1),dims(2),dims(3),dims(4),&
		nx1i,nx2i,nx3i,nx4i)
	return
end subroutine

subroutine myinterp4(fint,fx,x1,x2,x3,x4,x1i,x2i,x3i,x4i,&
	nx1,nx2,nx3,nx4,nx1i,nx2i,nx3i,nx4i)
    implicit none
	double precision, intent(out) :: fint(nx1i,nx2i,nx3i,nx4i)
	double precision, intent(in)  :: fx(nx1,nx2,nx3,nx4)
	double precision, intent(in)  :: x1(nx1)
	double precision, intent(in)  :: x2(nx2)
	double precision, intent(in)  :: x3(nx3)
	double precision, intent(in)  :: x4(nx4)
	double precision, intent(in)  :: x1i(nx1i)
	double precision, intent(in)  :: x2i(nx2i)
	double precision, intent(in)  :: x3i(nx3i)
	double precision, intent(in)  :: x4i(nx4i)
	mwSize, intent(in)	      :: nx1,nx2,nx3,nx4,nx1i,nx2i,nx3i,nx4i
	! Local variables
	integer			      :: i1,i2,i3,i4,j1,j2,j3,j4
	integer			      :: idx1,idx2,idx3,idx4
	double precision	      :: wx1(2),wx2(2),wx3(2),wx4(2)
	
	! Loop over all dimensions to interpolate
	do i1=1,nx4i
		! find neigbor,
		call locate(x4,nx4,x4i(i1),idx4)
		idx4=max(min(idx4,nx4-1),1)
		! weight on right neigbor
		wx4(2)=(x4i(i1)-x4(idx4))/(x4(idx4+1)-x4(idx4))
		wx4(1)=1.d0-wx4(2)
		do i2=1,nx3i
			call locate(x3,nx3,x3i(i2),idx3)
			idx3=max(min(idx3,nx3-1),1)
			wx3(2)=(x3i(i2)-x3(idx3))/(x3(idx3+1)-x3(idx3))
			wx3(1)=1.d0-wx3(2)
			do i3=1,nx2i
				call locate(x2,nx2,x2i(i3),idx2)
				idx2=max(min(idx2,nx2-1),1)
				wx2(2)=(x2i(i3)-x2(idx2))/(x2(idx2+1)-x2(idx2))
				wx2(1)=1.d0-wx2(2)
!$OMP PARALLEL DO SHARED(fint, i1,i2,i3,wx2,wx3,wx4,fx,idx2,idx3,idx4,nx1,x1i) PRIVATE(i4,idx1,wx1) 
				do i4=1,nx1i
					call locate(x1,nx1,x1i(i4),idx1)
					idx1=max(min(idx1,nx1-1),1)
					wx1(2)=(x1i(i4)-x1(idx1))/(x1(idx1+1)-x1(idx1))
					wx1(1)=1.d0-wx1(2)
					
					! Interpolation step
                    fint(i4,i3,i2,i1)=fint(i4,i3,i2,i1)+&
                                      wx1(1)*wx2(1)*wx3(1)*wx4(1)*fx(idx1,idx2,idx3,idx4)+&
                                      wx1(1)*wx2(1)*wx3(1)*wx4(2)*fx(idx1,idx2,idx3,idx4+1)+&
                                      wx1(1)*wx2(1)*wx3(2)*wx4(1)*fx(idx1,idx2,idx3+1,idx4)+&
                                      wx1(1)*wx2(1)*wx3(2)*wx4(2)*fx(idx1,idx2,idx3+1,idx4+1)+&
                                      wx1(1)*wx2(2)*wx3(1)*wx4(1)*fx(idx1,idx2+1,idx3,idx4)+&
                                      wx1(1)*wx2(2)*wx3(1)*wx4(2)*fx(idx1,idx2+1,idx3,idx4+1)+&
                                      wx1(1)*wx2(2)*wx3(2)*wx4(1)*fx(idx1,idx2+1,idx3+1,idx4)+&
                                      wx1(1)*wx2(2)*wx3(2)*wx4(2)*fx(idx1,idx2+1,idx3+1,idx4+1)+&
                                      wx1(2)*wx2(1)*wx3(1)*wx4(1)*fx(idx1+1,idx2,idx3,idx4)+&
                                      wx1(2)*wx2(1)*wx3(1)*wx4(2)*fx(idx1+1,idx2,idx3,idx4+1)+&
                                      wx1(2)*wx2(1)*wx3(2)*wx4(1)*fx(idx1+1,idx2,idx3+1,idx4)+&
                                      wx1(2)*wx2(1)*wx3(2)*wx4(2)*fx(idx1+1,idx2,idx3+1,idx4+1)+&
                                      wx1(2)*wx2(2)*wx3(1)*wx4(1)*fx(idx1+1,idx2+1,idx3,idx4)+&
                                      wx1(2)*wx2(2)*wx3(1)*wx4(2)*fx(idx1+1,idx2+1,idx3,idx4+1)+&
                                      wx1(2)*wx2(2)*wx3(2)*wx4(1)*fx(idx1+1,idx2+1,idx3+1,idx4)+&
                                      wx1(2)*wx2(2)*wx3(2)*wx4(2)*fx(idx1+1,idx2+1,idx3+1,idx4+1)


					!do j1=1,2
					!do j2=1,2
					!do j3=1,2
					!do j4=1,2
					!fint(i4,i3,i2,i1)=fint(i4,i3,i2,i1)+&
					!	wx1(j4)*wx2(j3)*wx3(j2)*wx4(j1)*&
					!	fx(idx1+j4-1,idx2+j3-1,&
					!	idx3+j2-1,idx4+j1-1)
					!end do
					!end do
					!end do
					!end do

				end do
!$OMP END PARALLEL DO
			end do
		end do
	end do		
	
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
