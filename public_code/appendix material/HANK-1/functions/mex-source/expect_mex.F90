!
!   expect_mex.F90
!
!
!   Created by Volker Tjaden on 30.07.13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
#include "fintrf.h"
!#include "blas.h"
!function [Evalue] =  expect_mex(value,index,weight,beta)
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
    mwPointer     ::  mxGetN , mxGetM, mxGetData


! Pointers to input/output mxArrays:
! Input
    mwPointer  :: value,index,weight,beta
! Output
    mwPointer  :: Evalue

! Array size information
    mwSize    :: nrows,ncols,nind


!--------------------------------------------------------------------------------


    value   = mxGetPr(prhs(1))
    index = mxGetData(prhs(2))
    weight   = mxGetPr(prhs(3))
    beta = mxGetPr(prhs(4))


	! Retrieve size information
    nrows = mxGetM(prhs(1))
    ncols = mxGetN(prhs(1))
    nind  = mxGetN(prhs(2))


    plhs(1)=mxCreateDoubleMatrix(nrows, ncols, real(8)) ! Reserve memory
    Evalue = mxGetPr(plhs(1)) !copy pointers



	call interpE(nrows,nind,%val(Evalue),%val(value),%val(index),%val(weight),%val(beta))

end subroutine

subroutine interpE(nrows,nind,Evalue,value,index,weight,beta)

		implicit none
        ! Input Variables
		mwSize, intent(in)      :: nrows ! Size of Problem
		mwSize, intent(in)      :: nind ! Size of Problem
		real(8), intent(in)     :: beta

		integer(4),intent(in)   :: index(nrows,nind)

		real(8),intent(in)      :: value(nrows)
    real(8),intent(in)      :: weight(nrows,nind)

    !Output Variables
    real(8),intent(out)      :: Evalue(nrows)

    ! Local variables
		integer(8) :: i1,i2
    real(8)    :: aux


		do i1=1,nrows
        Evalue(i1)=dot_product(value(index(i1,:)), weight(i1,:))
    end do

    Evalue=beta*Evalue

	end subroutine
