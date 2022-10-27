!
!   expect_mex.F90
!
!
!   Created by Volker Tjaden on 30.07.13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
#include "fintrf.h"
!#include "blas.h"
!function [expected] =  expect_mex_new(P,P_KM,input)
!----------------------------------------------------------------------------
! Gateway routine
!----------------------------------------------------------------------------


subroutine mexFunction(nlhs,plhs,nrhs,prhs)

    implicit none

    mwPointer  :: plhs(*), prhs(*)
    integer    :: nlhs, nrhs

! Function declarations
    mwPointer   ::  mxGetPr, mxCreateNumericArray, mxGetDimensions
	  mwSize      ::  mxGetNumberOfDimensions
    integer(4)  ::  mxClassIDFromClassName,  mxGetString


! Pointers to input/output mxArrays:
! Input
    mwPointer  :: P, P_KM, input

! Output
    mwPointer  :: expected

! Array size information
  mwSize    :: ns,nh,nm,nk,ncapK,ncapM
	mwSize	  :: ndim
	mwSize	  :: dims(5), dims2(3)
	integer(4)   :: classid, complexflag


!--------------------------------------------------------------------------------


    P       = mxGetPr(prhs(1))
    P_KM    = mxGetPr(prhs(2))
    input   = mxGetPr(prhs(3))

  ! Retrieve size information
  ndim=mxGetNumberOfDimensions(prhs(3))
  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(3)), dims, ndim)
  nm=dims(1)
  nk=dims(2)
  nh=dims(3)
  ncapM=dims(4)
  ncapK=dims(5)
  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(2)), dims2, mxGetNumberOfDimensions(prhs(2)))
  ns=dims2(3)
  nh=nh/ns
  ! Prepare output
  classid 	= mxClassIDFromClassName('double')
  complexflag = 0
  plhs(1) 	= mxCreateNumericArray(ndim, dims, classid, complexflag)
  expected	= mxGetPr(plhs(1))



	call expectations(nm,nk,ns,nh,ncapK,ncapM,%val(expected),%val(input),%val(P),%val(P_KM))

end subroutine

subroutine expectations(nm,nk,ns,nh,ncapK,ncapM,output,input,P,P_KM)
  ! input/output variables
  implicit none
  mwSize, intent(in)  :: ns
  mwSize, intent(in)  :: nh
  mwSize, intent(in)  :: nm
  mwSize, intent(in)  :: nk
  mwSize, intent(in)  :: ncapK
  mwSize, intent(in)  :: ncapM
  real(8),intent(out)     :: output(nm,nk,ns*nh,ncapM,ncapK)
  real(8),intent(in)      :: input(nm,nk,ns*nh,ncapM,ncapK)
  real(8),intent(in)      :: P(ns*nh,ns*nh)
  real(8),intent(in)      :: P_KM(ncapM*ncapK,ncapM*ncapK,ns)

  ! local variables
  integer(4)              :: i1,i2,i3,i4,i5,i6,inext,Knext,Mnext
  integer(4)              :: shind !,ifail
  real(8)                 :: trprob,mprob,kprob, auxaux(nm*nk,ncapM*ncapK), auxaux2(nm*nk,ncapM*ncapK)
  real(8)                 :: temp(nm*nk,ns*nh,ncapM,ncapK), aux(nm*nk,ns*nh,ncapM,ncapK)


  temp=0.d0
  output=0.d0


  ! 1. Take expectations over sxh
  aux=reshape(input,(/nm*nk,ns*nh,ncapM,ncapK/))
  do i1=1,ncapK
    do i2=1,ncapM
        call DGEMM('N','T', nk*nm, ns*nh,ns*nh ,1.d0, aux(:,:,i2,i1),nk*nm, P, ns*nh,0.d0, temp(:,:,i2,i1), nk*nm)
    end do
  end do



  ! 2. Transitions in M and K
  aux=0.d0
  do i4=1,nh
    do i3=1,ns
      shind=ns*(i4-1)+i3
      auxaux=reshape(temp(:,shind,:,:),(/nm*nk,ncapM*ncapK/))
      !auxaux=MATMUL(auxaux,TRANSPOSE(P_KM(:,:,i3)))
      call DGEMM('N','T', nk*nm, ncapM*ncapK,ncapM*ncapK ,1.d0, auxaux,nk*nm, P_KM(:,:,i3), ncapM*ncapK,0.d0, auxaux2, nk*nm)
      do i1=1,ncapK
        do i2=1,ncapM
          aux(:,shind,i2,i1)=auxaux2(:,i2+(i1-1)*ncapM)
        end do
      end do
    end do
  end do
  output=reshape(aux,(/nm,nk,ns*nh,ncapM,ncapK/))

	end subroutine
