!
!   EGM_Step5b_mex.F90
!
!
!   Created by Volker Tjaden on 30.07.13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
#include "fintrf.h"

!function [ JD_next ] =  JDIteration(JD,idm_n,idm_a,idk_a,dist_m_n,dist_m_a,dist_k_a,AProb,P)
!----------------------------------------------------------------------------
! Gateway routine
!----------------------------------------------------------------------------


subroutine mexFunction(nlhs,plhs,nrhs,prhs)
  implicit none

  ! Arguments in mexFunction
  mwPointer  :: plhs(*), prhs(*)
  integer    :: nlhs, nrhs

  ! Function declarations
  mwPointer   ::  mxGetPr, mxCreateNumericArray
  mwPointer     mxGetDimensions
  mwSize        mxGetNumberOfDimensions
  integer(4)    mxClassIDFromClassName
  integer(4)    mxGetString


  ! Pointers to input/output mxArrays:
  ! Input
  mwPointer  :: JD,idm_n,idm_a,idk_a,AProb,P,dist_m_a,dist_k_a,dist_m_n

  ! Output
  mwPointer  :: JD_next

  ! Array size information
  mwSize    :: nh,nm,nk
  mwSize	  :: ndim
  mwSize	  :: dims(3)
  integer(4)   :: classid
  integer(4)   :: complexflag

  !--------------------------------------------------------------------------------


  JD    = mxGetPr(prhs(1))
  idm_n    = mxGetPr(prhs(2))
  idm_a   = mxGetPr(prhs(3))
  idk_a      = mxGetPr(prhs(4))
  dist_m_n   = mxGetPr(prhs(5))
  dist_m_a   = mxGetPr(prhs(6))
  dist_k_a   = mxGetPr(prhs(7))
  AProb       = mxGetPr(prhs(8))
  P           = mxGetPr(prhs(9))


  ! Retrieve size information
  ndim=mxGetNumberOfDimensions(prhs(1))
  call mxCopyPtrToPtrArray(mxGetDimensions(prhs(1)), dims, ndim)
  nm=dims(1)
  nk=dims(2)
  nh=dims(3)



  ! Create matrices for return arguments
  classid 	= mxClassIDFromClassName('double')
  complexflag = 0
  plhs(1) 	= mxCreateNumericArray(ndim, dims, classid, complexflag)
  JD_next	= mxGetPr(plhs(1))

  call JDIteration_aux(nm,nk,nh,%val(JD),&
  %val(idm_n),%val(idm_a),%val(idk_a),&
  %val(dist_m_n),%val(dist_m_a),%val(dist_k_a),&
  %val(Aprob),%val(P),%val(JD_next))

end subroutine

subroutine JDIteration_aux(nm,nk,nh,JD,idm_n,idm_a,idk_a,dist_m_n,dist_m_a,dist_k_a,Aprob,P,JD_next)
  implicit none
  mwSize, intent(in)      :: nm, nk, nh
  real(8), intent(in)     :: JD(nm,nk,nh)
  integer(4),intent(in)   :: idm_n(nm,nk,nh),idm_a(nm,nk,nh),idk_a(nm,nk,nh)
  real(8), intent(in)     :: dist_m_n(nm,nk,nh),dist_m_a(nm,nk,nh),dist_k_a(nm,nk,nh)
  real(8),intent(in)      :: AProb(nm,nk,nh), P(nh,nh)
  real(8), intent(out)     :: JD_next(nm,nk,nh)
  ! Local variables
  integer(4)              :: i1,i2,i3,i4,i5, i1_target,i2_target
  real(8)                 :: ap,aux

  ! intialize JD_next
  JD_next=0.d0
  ! 1. Calculate weight on right neighbour in interpolation
  do i4=1,nh
    do i3=1,nh
      do i2=1,nk
        do i1=1,nm
          ! no adjustment
          i1_target=idm_n(i1,i2,i3)
          i2_target=i2
          ap=1.0-Aprob(i1,i2,i3)
          aux=ap*P(i3,i4)*JD(i1,i2,i3)

          JD_next(i1_target,i2_target,i4)=JD_next(i1_target,i2_target,i4)+&
            aux*(1-dist_m_n(i1,i2,i3))

          i1_target=i1_target+1
          JD_next(i1_target,i2_target,i4)=JD_next(i1_target,i2_target,i4)+&
            aux*(dist_m_n(i1,i2,i3))

          ! Adjustment
          i1_target=idm_a(i1,i2,i3)
          i2_target=idk_a(i1,i2,i3)
          ap=Aprob(i1,i2,i3)
          aux=ap*P(i3,i4)*JD(i1,i2,i3)

          JD_next(i1_target,i2_target,i4)=JD_next(i1_target,i2_target,i4)+&
            aux*(1-dist_m_a(i1,i2,i3))*(1-dist_k_a(i1,i2,i3))

          i1_target=i1_target+1
          JD_next(i1_target,i2_target,i4)=JD_next(i1_target,i2_target,i4)+&
            aux*dist_m_a(i1,i2,i3)*(1-dist_k_a(i1,i2,i3))

          i2_target=i2_target+1
          JD_next(i1_target,i2_target,i4)=JD_next(i1_target,i2_target,i4)+&
            aux*dist_m_a(i1,i2,i3)*(dist_k_a(i1,i2,i3))

          i1_target=i1_target-1
          JD_next(i1_target,i2_target,i4)=JD_next(i1_target,i2_target,i4)+&
            aux*(1-dist_m_a(i1,i2,i3))*(dist_k_a(i1,i2,i3))

        end do
      end do
    end do
  end do



end subroutine
