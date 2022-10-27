!
!   egm5b_aux_mex.F90
!   
!
!   Created by Christian Bayer on 26.07.13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
#include "fintrf.h"

!function [ psi_new ] =  EGM_Step5b(E_rhs_mu,m_n_star,psi_guess,grid,idm_n,P,TT,beta,nu,REP,nk,nMxnH)
!----------------------------------------------------------------------------
! Gateway routine
!----------------------------------------------------------------------------
!
! TT=[index1 index2 weight] matrix of sparse-function entries (obtained by find(TT) in MATLAB.

subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    implicit none

    ! Arguments in mexFunction
    mwPointer  :: plhs(*), prhs(*)
    integer    :: nlhs, nrhs

! Function declarations
    mwPointer   ::  mxGetPr, mxCreateDoubleMatrix
    mwPointer   ::  mxGetN , mxGetM

! Pointers to input/output mxArrays:
! Input
    mwPointer  :: E_rhs_mu, m_n_star, psi_guess, gridm, idm, P, TT,TTw
    mwPointer  :: beta, nu, REP, nk, nAgg
! Output
    mwPointer  :: psi_new


! Array size information
    mwSize    :: nm,ncols,nTT,nsh            ! length of input array


!--------------------------------------------------------------------------------


    E_rhs_mu    = mxGetPr(prhs(1))
    m_n_star    = mxGetPr(prhs(2))
    psi_guess   = mxGetPr(prhs(3))
    gridm       = mxGetPr(prhs(4))
    idm         = mxGetPr(prhs(5))
    P           = mxGetPr(prhs(6))
    TT          = mxGetPr(prhs(7))
    TTw          = mxGetPr(prhs(8))
    beta        = mxGetPr(prhs(9))
    nu          = mxGetPr(prhs(10))
    REP         = mxGetPr(prhs(11))
    nk          = mxGetPr(prhs(12))
    nAgg        = mxGetPr(prhs(13))

    ! Retrieve size of x
    nm = mxGetM(prhs(2))
    ncols = mxGetN(prhs(2))
    nsh = mxGetM(prhs(6))
    nTT = mxGetM(prhs(7))


    plhs(1)=mxCreateDoubleMatrix(nm, ncols, int(0,4)) ! Reserve memory

    
    psi_new = mxGetPr(plhs(1))

    call Howards(nm, ncols, nsh, nTT, %val(nk),%val(nAgg),%val(E_rhs_mu), %val(m_n_star), %val(psi_guess),&
            %val(idm), %val(gridm), %val(P), %val(TT),%val(TTw), %val(beta), %val(nu), %val(REP), %val(psi_new))

    return
end subroutine


subroutine Howards(nm, ncols, nsh, nTT, nk, nAgg, E_rhs_mu, m_n_star, psi_guess, idm, gridm &
                ,P, TT, TTw,beta, nu, REP, psi_update)
    implicit none
! Input arguments
    mwSize, intent(in)      :: nm, ncols, nsh, nTT                            ! number of points in x-vector
    mwSize, intent(in)      :: nk, nAgg
    real(8), intent(in)                             :: gridm(nm)                          ! x-values
    real(8), dimension(nm,ncols), intent(in)        :: E_rhs_mu, m_n_star, psi_guess
    mwSize, dimension(nm,ncols), intent(in)        :: idm
    real(8), intent(in)                             :: P(nsh,nsh)    ! Transition Probabilities
    mwSize, intent(in)                              :: TT(nTT,2)   ! Transitions in Aggregate
    real(8), intent(in)                             :: TTw(nTT)   ! Transitions in Aggregate
    real(8), intent(in)                             :: beta, nu ! Parameters
    integer, intent(in)                             :: REP
    
! Output arguments
    real(8), intent(out)        :: psi_update(nm,ncols)  ! psi update

! Local variables
    integer                     :: i1,i2,i3,index
    real(8), dimension(nm,ncols) :: weight, mu_flow,ww, outaux	                  ! interpolation weight on left neighbor
    real(8)                     :: aux

    
    ww=0
    do i1 = 1,nm
        do i2 = 1,ncols
            index=idm(i1,i2)
            weight(i1,i2)= (m_n_star(i1,i2)-gridm(index))/(gridm(index+int(1,4))-gridm(index)) ! weights on higher gridpoint
            aux=E_rhs_mu(index+int(1,4),i2)-E_rhs_mu(index,i2)
            mu_flow(i1,i2) = E_rhs_mu(index,i2) + weight(i1,i2)*aux
        end do
    end do
    !call sub_mexWriteString('here')

    call Transitions(psi_guess,TT,TTw,P,nm, ncols, nsh, nk, nAgg,nTT, outaux)
    
    do i1 = 1,nm
        do i2 = 1,ncols
            index=idm(i1,i2)
            aux=outaux(index+int(1,4),i2)-outaux(index,i2)
            ww(i1,i2) = mu_flow(i1,i2)+ beta*(1-nu)*(outaux(index,i2) + weight(i1,i2)*aux)
        end do
    end do


    do i3=1,REP
        call Transitions(ww,TT,TTw,P,nm, ncols, nsh, nk, nAgg,nTT, outaux)
        !outaux=outaux
        do i1 = 1,nm
            do i2 = 1,ncols
                index=idm(i1,i2)
                aux=outaux(index+int(1,4),i2)-outaux(index,i2)
                ww(i1,i2) = mu_flow(i1,i2)+ beta*(1-nu)*(outaux(index,i2) + weight(i1,i2)*aux)
            end do
        end do
    end do
    
    psi_update=ww
return
end subroutine Howards

!----------------------------------------------------------------------------------------
subroutine Transitions(var2tran,TT,TTw,P,nm, ncols, nsh, nk, nAgg,nTT,varout)
!---------------------------------------------------------------------------------------
    implicit none
! Input arguments
    mwSize, intent(in)      :: nm, ncols, nsh,nk,nAgg,nTT                             ! number of points in x-vector

    real(8), dimension(nm,ncols) , intent(in)        :: var2tran
    real(8), intent(in)         :: P(nsh,nsh)    ! Transition Probabilities
    real(8), intent(in)         :: TTw(nTT)   ! Transitions in Aggregate
    mwSize, intent(in)          :: TT(nTT,2)   ! Transitions in Aggregate
    real(8), dimension(nm,ncols) , intent(out)        :: varout
! Local variables
    mwSize                     :: i1,i2,i3,i4,i5, index1,index2,index3,index4
    real(8), dimension(nm,ncols)        :: aux

    varout=0

! Stochastic transitions
    do i1 = 1,nm
        do i2=1,nk
            do i3=1,nsh
                do i4 = 1,nAgg
                    index1=i2 + (i3-1)*nk +(i4-1)*nk*nsh
                    do i5=1,nsh
                        index2= i2 + (i5-int(1,4))*nk +(i4-int(1,4))*nk*nsh
                        varout(i1,index1)= varout(i1,index1)+P(i3,i5)*var2tran(i1,index2)
                    end do
                end do
            end do
        end do
    end do



! Multiplication of sparse matrix and vector S*x
!    aux =0
!    do i1=1,nTT
!        index1 = TT(i1,1) - nm*int((TT(i1,1)-1)/nm)
!        index2 = int((TT(i1,1)-1)/nm)+1
!        index3 = TT(i1,2) - nm*aint(((TT(i1,2)-1)/nm))
!        index4 = int((TT(i1,2)-1)/nm)+1
!        aux(index1,index2)=aux(index3,index4) + TTw(i1)*varout(index3,index4)
!    end do
!    varout=aux
    return
end subroutine Transitions

subroutine sub_mexWriteString(String)
implicit none
character(*),intent(in):: String
call mexPrintf(achar(10))
call mexPrintf(trim(adjustl(String)))
end subroutine sub_mexWriteString
