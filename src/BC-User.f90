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
  PUBLIC :: init_user, boundary_conditions_user, postprocess_user

contains

  subroutine init_user (rho1,ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    real(mytype) :: x, y, z, r, x0, y0, z0, r0
    integer :: k,j,i,is

    if (iscalar==1) then
       phi1(:,:,:,:) = zero
    endif

    if (iin.eq.0) then !empty domain

       if (nrank==0) write(*,*) "Empty initial domain!"

       ux1=zero; uy1=zero; uz1=zero

    endif

    !! Set the circle centre
    x0 = zero
    y0 = zero !-yly / four
    z0 = zero

    !! Set the circle radius
    r0 = one / four
    
    if (iin.eq.1) then !generation of a random noise

       !INIT FOR G AND U=MEAN FLOW + NOISE
       do k=1,xsize(3)
          z = real(k + xstart(3) - 2, mytype) * dz - half * zlz
          do j=1,xsize(2)
             y = real(j + xstart(2) - 2, mytype) * dy - half * yly
             do i=1,xsize(1)
                x = real(i + xstart(1) - 2, mytype) * dx - half * xlx

                r = sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2) !! Sphere
                ! r = sqrt((x - x0)**2 + (y - y0)**2) !! Cylinder
                
                ux1(i,j,k) = u1
                uy1(i,j,k) = u2
                uz1(i,j,k) = zero

                rho1(i,j,k,1) = one

                !! Circle
                if (r.lt.r0) then
                   phi1(i, j, k, :) = (r0 - r)
                else if (r.gt.r0) then
                   phi1(i, j, k, :) = -(r - r0)
                else
                   phi1(i, j, k, :) = zero
                endif

                ! !! Square
                ! if ((abs(x).lt.0.5_mytype).and.(abs(y).lt.0.5_mytype)) then
                !    phi1(i, j, k, :) = one
                ! else
                !    phi1(i, j, k, :) = -one
                ! endif

                ! !! Free surface
                ! phi1(i, j, k, :) = y
             enddo
          enddo
       enddo

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
       bxx1(:,:) = zero
       bxy1(:,:) = zero
       bxz1(:,:) = zero
    ENDIF
    IF (nclxn.EQ.2) THEN
       bxxn(:,:) = zero
       bxyn(:,:) = zero
       bxzn(:,:) = zero
    ENDIF

    IF (ncly1.EQ.2) THEN
       byx1(:,:) = zero
       byy1(:,:) = u2
       byz1(:,:) = zero
    ENDIF
    IF (nclyn.EQ.2) THEN
       byxn(:,:) = zero
       byyn(:,:) = u2
       byzn(:,:) = zero
    ENDIF

    IF (nclz1.EQ.2) THEN
       bzx1(:,:) = zero
       bzy1(:,:) = zero
       bzz1(:,:) = zero
    ENDIF
    IF (nclzn.EQ.2) THEN
       bzxn(:,:) = zero
       bzyn(:,:) = zero
       bzzn(:,:) = zero
    ENDIF

  end subroutine boundary_conditions_user

  subroutine postprocess_user(ux1,uy1,uz1,phi1,ep1)

    USE decomp_2d
    USE MPI

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

  end subroutine postprocess_user

end module user_sim
