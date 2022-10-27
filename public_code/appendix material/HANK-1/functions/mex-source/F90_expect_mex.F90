!
!   expect_mex.F90
!
!
!   Created by Volker Tjaden on 30.07.13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
!#include "fintrf.h"
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
    mwPointer     ::  mxGetN , mxGetM


! Pointers to input/output mxArrays:
! Input
    mwPointer  :: value,index,weight,beta
! Output
    mwPointer  :: Evalue

! Array size information
    mwSize    :: nrows,ncols,nind


!--------------------------------------------------------------------------------


    value   = mxGetPr(prhs(1))
    index = mxGetPr(prhs(2))
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
		!mwSignedIndex(8),intent(in)   :: index(nrows,nind)
		integer(8),intent(in)   :: index(nrows,nind)

		real(8),intent(in)      :: value(nrows)
    real(8),intent(in)      :: weight(nrows,nind)

    !Output Variables
    real(8),intent(out)      :: Evalue(nrows)

    ! Local variables
		integer(8) :: i1
    !real(8)    :: aux

		do i1=1,nrows
        Evalue(i1)=DOT_PRODUCT(value(index(i1,:)), weight(i1,:))
    !        aux=0
    !        do i2=1,nind

    !            aux=aux+value(index(i1,i2))*weight(i1,i2)

    !        end do
    !        Evalue(i1)=aux
    end do

    Evalue=beta*Evalue

	end subroutine
