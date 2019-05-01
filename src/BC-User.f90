!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-User.f90
!!!      AUTHOR: Paul Bartholomew
!!! DESCRIPTION: This module allows a user to setup their own flows.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module user_sim

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_user, boundary_conditions_user, postprocessing_user

contains

  subroutine init_user (rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    real(mytype) :: x
    integer :: k,j,i,is

    if (iin.eq.0) then !empty domain

       if (nrank==0) write(*,*) "Empty initial domain!"

       ux1=zero; uy1=zero; uz1=zero

    endif

    if (iin.eq.1) then !generation of a random noise

       !INIT FOR G AND U=MEAN FLOW + NOISE
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                x = real(i + xstart(1) - 2, mytype) * dx
                
                ux1(i,j,k) = one
                uy1(i,j,k) = zero
                uz1(i,j,k) = zero

                rho1(i,j,k,1) = one

                !! Set phi1 like level-set
                !! XXX level-set is a DISTANCE function
                if ((x.gt.(0.4_mytype*xlx)).and.(x.lt.(0.6_mytype*xlx))) then
                   phi1(i, j, k, :) = min(abs(x - 0.4_mytype * xlx), abs(x - 0.6_mytype))
                else
                   if (x.gt.half*xlx) then
                      phi1(i, j, k, :) = -min(abs(x - 0.6_mytype * xlx), abs((xlx - x) + 0.4_mytype * xlx))
                   else
                      phi1(i, j, k, :) = -min(abs(x - 0.4_mytype * xlx), abs(x + (xlx - 0.6_mytype * xlx)))
                   endif
                endif

                ! x = abs(x - half * xlx)
                ! phi1(i, j, k, :) = tanh((xlx / ten - x) / dx)

                !!------------------------------------------------------------
                !! XXX Do not touch this bit
                !!------------------------------------------------------------
                dux1(i,j,k,1)=ux1(i,j,k)
                duy1(i,j,k,1)=uy1(i,j,k)
                duz1(i,j,k,1)=uz1(i,j,k)

                drho1(i,j,k,1) = rho1(i,j,k,1)
                do is = 2, ntime
                   dux1(i,j,k,is)=dux1(i,j,k,is - 1)
                   duy1(i,j,k,is)=duy1(i,j,k,is - 1)
                   duz1(i,j,k,is)=duz1(i,j,k,is - 1)

                   drho1(i,j,k,is) = drho1(i,j,k,is - 1)
                enddo
                do is = 2, nrhotime
                   rho1(i,j,k,is) = rho1(i,j,k,is - 1)
                enddo
                !!------------------------------------------------------------
             enddo
          enddo
       enddo

       if (iscalar==1) then
          !! Reinitialise to smooth level-set
          call reinit_ls(phi1(:,:,:,1))

          dphi1(:,:,:,1,:) = phi1(:,:,:,:)
          do is = 2,ntime
             dphi1(:,:,:,is,:) = dphi1(:,:,:,is - 1,:)
          enddo
       endif

    endif

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_user

  subroutine boundary_conditions_user (ux,uy,uz,phi,ep)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    IF (nclx1.EQ.2) THEN
    ENDIF
    IF (nclxn.EQ.2) THEN
    ENDIF

    IF (ncly1.EQ.2) THEN
    ENDIF
    IF (nclyn.EQ.2) THEN
    ENDIF

    IF (nclz1.EQ.2) THEN
    ENDIF
    IF (nclzn.EQ.2) THEN
    ENDIF

  end subroutine boundary_conditions_user

  subroutine postprocessing_user(ux1,uy1,uz1,phi1,ep1)

    USE decomp_2d
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

  end subroutine postprocessing_user

end module user_sim
