!
!   EGM_Step5b_mex.F90
!
!
!   Created by Volker Tjaden on 30.07.13.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!
#include "fintrf.h"

!function [ psi_new ] =  EGM_Step5b(E_rhs_mu,m_n_star,psi_guess,grid,idm_n,P,P_M,P_K,beta,AProb,REP)
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
    mwPointer  :: E_rhs_mu, m_n_star, psi_guess, grid_m, idm_n, P
	mwPointer  :: P_M,P_K
    mwPointer  :: beta, AProb, rep
! Output
    mwPointer  :: psi_new

! Array size information
    mwSize    :: ns,nh,nm,nk,ncapK,ncapM
	mwSize	  :: ndim
	mwSize	  :: dims(5)
	mwSize	  :: dims2(3)
	integer(4)   :: classid
    integer(4)   :: complexflag

!--------------------------------------------------------------------------------


    E_rhs_mu    = mxGetPr(prhs(1))
    m_n_star    = mxGetPr(prhs(2))
    psi_guess   = mxGetPr(prhs(3))
    grid_m      = mxGetPr(prhs(4))
    idm_n       = mxGetPr(prhs(5))
    P           = mxGetPr(prhs(6))
    P_M         = mxGetPr(prhs(7))
    P_K         = mxGetPr(prhs(8))
    beta        = mxGetPr(prhs(9))
    AProb       = mxGetPr(prhs(10))
    REP         = mxGetPr(prhs(11))

	! Retrieve size information
	ndim=mxGetNumberOfDimensions(prhs(1))
	call mxCopyPtrToPtrArray(mxGetDimensions(prhs(1)), dims, ndim)
    nm=dims(1)
	nk=dims(2)
	nh=dims(3)
	ncapM=dims(4)
	ncapK=dims(5)
	call mxCopyPtrToPtrArray(mxGetDimensions(prhs(8)), dims2, mxGetNumberOfDimensions(prhs(8)))
	ns=dims2(3)
	nh=nh/ns

	! Create matrices for return arguments
    classid 	= mxClassIDFromClassName('double')
    complexflag = 0
	plhs(1) 	= mxCreateNumericArray(ndim, dims, classid, complexflag)
	psi_new		= mxGetPr(plhs(1))

	call egm_step5b_aux(nm,nk,ns,nh,ncapK,ncapM,&
		%val(psi_new),%val(e_rhs_mu),%val(m_n_star),%val(psi_guess),&
		%val(grid_m),%val(idm_n),%val(P),%val(P_M),%val(P_K),&
		%val(beta),%val(AProb),%val(rep))

end subroutine

	subroutine egm_step5b_aux(nm,nk,ns,nh,ncapK,ncapM,&
		psi_new,e_rhs_mu,m_n_star,psi_guess,grid_m,&
		idm_n,P,P_M,P_K,beta,AProb,rep)
		implicit none
		mwSize, intent(in)  :: ns
		mwSize, intent(in)  :: nh
		mwSize, intent(in)  :: nm
		mwSize, intent(in)  :: nk
		mwSize, intent(in)  :: ncapK
		mwSize, intent(in)  :: ncapM
		real(8), intent(in)     :: beta
		integer(4),intent(in)   :: rep
		real(8),intent(in)      :: e_rhs_mu(nm,nk,ns*nh,ncapM,ncapK)
		real(8),intent(in)      :: grid_m(nm)
		real(8),intent(in)      :: m_n_star(nm,nk,ns*nh,ncapM,ncapK)
        real(8),intent(in)      :: AProb(nm,nk,ns*nh,ncapM,ncapK)
		real(8),intent(in)      :: psi_guess(nm,nk,ns*nh,ncapM,ncapK)
		integer(4),intent(in)   :: idm_n(nm,nk,ns*nh,ncapM,ncapK)
		real(8),intent(in)      :: P(ns*nh,ns*nh)
		real(8),intent(in)      :: P_M(ncapM,ncapM,ns,ncapK)
		real(8),intent(in)      :: P_K(ncapK,ncapK,ns,ncapM)
		real(8),intent(out)      :: psi_new(nm,nk,ns*nh,ncapM,ncapK)


		 ! Local variables
		integer(4)              :: i1,i2,i3,i4,i5
		integer(4)              :: idx
		real(8),allocatable     :: dist_m_n(:,:,:,:,:)
		real(8),allocatable     :: mu_flow(:,:,:,:,:)
		real(8),allocatable     :: aux(:,:,:,:,:)
		real(8),allocatable     :: ww(:,:,:,:,:)
		real(8)					:: aux1, aux2
		! 0. Allocate local arrays
		allocate(dist_m_n(nm,nk,ns*nh,ncapM,ncapK))
		allocate(mu_flow(nm,nk,ns*nh,ncapM,ncapK))
		allocate(aux(nm,nk,ns*nh,ncapM,ncapK))
		allocate(ww(nm,nk,ns*nh,ncapM,ncapK))

		! 1. Calculate weight on right neighbour in interpolation

		do i1=1,ncapK
			do i2=1,ncapM
				do i3=1,ns*nh
					do i4=1,nk
						do i5=1,nm
							idx=idm_n(i5,i4,i3,i2,i1)
							aux1=grid_m(idx)
              aux2=grid_m(idx+1)
							dist_m_n(i5,i4,i3,i2,i1)=&
								(m_n_star(i5,i4,i3,i2,i1)-aux1)/(aux2-aux1)
						end do
					end do
				end do
			end do
		end do

		! 2. Interpolate e_rhs_mu to obtain mu_flow
		call interpolate_m_n(nm,nk,ns,nh,ncapK,ncapM,&
			mu_flow,e_rhs_mu,dist_m_n,idm_n)

		! 3. First Howard step / generate ww
		call expectations(nm,nk,ns,nh,ncapK,ncapM,&
				aux,psi_guess,P,P_M,P_K,beta,AProb)
		call interpolate_m_n(nm,nk,ns,nh,ncapK,ncapM,&
			ww,aux,dist_m_n,idm_n)

		ww=mu_flow+ww

		! 4. Howard's iteration steps

		do i1=1,rep
			call expectations(nm,nk,ns,nh,ncapK,ncapM,&
				aux,ww,P,P_M,P_K,beta,AProb)
			call interpolate_m_n(nm,nk,ns,nh,ncapK,ncapM,&
				ww,aux,dist_m_n,idm_n)
			ww=mu_flow+ww
		end do

		psi_new=ww

		deallocate(dist_m_n,mu_flow,aux,ww)
	end subroutine

	! subroutines: one for expectations/interpolation (first sh then PM)
     subroutine interpolate_m_n(nm,nk,ns,nh,ncapK,ncapM,&
        output,input,dist_m_n,idm_n)
        implicit none
        mwSize, intent(in)  :: ns
        mwSize, intent(in)  :: nh
        mwSize, intent(in)  :: nm
        mwSize, intent(in)  :: nk
        mwSize, intent(in)  :: ncapK
        mwSize, intent(in)  :: ncapM
        real(8),intent(out)     :: output(nm,nk,ns*nh,ncapM,ncapK)
        real(8),intent(in)      :: input(nm,nk,ns*nh,ncapM,ncapK)
        real(8),intent(in)      :: dist_m_n(nm,nk,ns*nh,ncapM,ncapK)
        integer(4),intent(in)   :: idm_n(nm,nk,ns*nh,ncapM,ncapK)
        ! local variables
        integer(4)              :: i1,i2,i3,i4,i5
        integer(4)              :: idx
        real(8)					:: aux1, aux2
        ! interpolation over policies
        !$OMP PARALLEL DO PRIVATE(i1, i2, i3, i4, i5, aux1, aux2) NUM_THREADS(3)
        do i1=1,ncapK
            do i2=1,ncapM
                do i3=1,ns*nh
                    do i4=1,nk
                        do i5=1,nm
                            idx=idm_n(i5,i4,i3,i2,i1)
                            aux1=input(idx,i4,i3,i2,i1)
                            aux2=input(idx+1,i4,i3,i2,i1)
                            output(i5,i4,i3,i2,i1)=aux1+&
                                dist_m_n(i5,i4,i3,i2,i1)*(aux2-aux1)
                        end do
                    end do
                end do
            end do
        end do
        !$ OMP END PARALLEL DO
    end subroutine

	subroutine expectations(nm,nk,ns,nh,ncapK,ncapM,&
		output,input,P,P_M,P_K,beta,AProb)
		!Use nag_library, Only: f01ctf
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
		real(8),intent(in)      :: P_M(ncapM,ncapM,ns,ncapK)
		real(8),intent(in)      :: P_K(ncapK,ncapK,ns,ncapM)
		real(8), intent(in)     :: beta
        real(8),intent(in)      :: AProb(nm,nk,ns*nh,ncapM,ncapK)

		! local variables
		integer(4)              :: i1,i2,i3,i4,i5,i6,inext,Knext,Mnext
		integer(4)              :: shind !,ifail
		real(8)                 :: trprob,mprob,kprob, auxaux(nm*nk), auxaux2(nm*nk)
		real(8)                 :: temp(nm*nk,ns*nh,ncapM,ncapK), aux(nm*nk,ns*nh,ncapM,ncapK)
        real(8)      :: CAProb(nm,nk,ns*nh,ncapM,ncapK)


		temp=0.d0
		output=0.d0
		CAProb=1.d0
        CAProb=CAProb-AProb

    ! 1. Take expectations over sxh
		aux=reshape(input,(/nm*nk,ns*nh,ncapM,ncapK/))
		do i1=1,ncapK
			do i2=1,ncapM
					call DGEMM('N','T', nk*nm, ns*nh,ns*nh ,1.d0, aux(:,:,i2,i1),nk*nm, P, ns*nh,0.d0, temp(:,:,i2,i1), nk*nm)
			end do
		end do

		! 2. Transitions in M and varH
		aux=0.d0
		do i1=1,ncapK
			do i2=1,ncapM
				do i4=1,nh
					do i3=1,ns
						shind=ns*(i4-1)+i3
						auxaux=0.d0
						do Knext=1,ncapK
							kprob=P_K(i1,Knext,i3,i2)
							do Mnext=1,ncapM
								mprob=P_M(i2,Mnext,i3,i1)
								auxaux2=temp(:,shind,Mnext,Knext)
								!do i5=1,nk
								!	do i6=1,nm
								!		call DAXPY(m, a, x, 1, yout, 1)
								!		call DAXPY(nt, trprob, auxaux2, 1, auxaux, 1)
										auxaux=kprob*mprob*auxaux2+auxaux
								!		aux(:,shind,i2,i1)=aux(:,shind,i2,i1)+&
								!			kprob*mprob*temp(:,shind,Mnext,Knext)
								!	end do
								!end do
							end do
						end do
						aux(:,shind,i2,i1)=auxaux
					end do
				end do

			end do
		end do
		output=reshape(aux,(/nm,nk,ns*nh,ncapM,ncapK/))
        !----CHECK THIS HERE----!
		output=beta*CAProb*output
        end subroutine
