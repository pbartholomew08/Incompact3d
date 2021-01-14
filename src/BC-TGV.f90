!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module tgv

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes !só vai funcionar se a impressão for em relação ao lapis X!
  integer :: nprobes
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tgv, boundary_conditions_tgv, postprocess_tgv

contains

  subroutine init_tgv (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,fh,ierror,is,code
    integer (kind=MPI_OFFSET_KIND) :: disp
    integer, dimension (:), allocatable :: seed
    integer ::  isize

    if (iscalar==1) then

       phi1(:,:,:,:) = zero
    endif

    if (iin.eq.0) then !empty domain

       if (nrank==0) write(*,*) "Empty initial domain!"

       ux1=zero; uy1=zero; uz1=zero

    endif

    if (iin.eq.1) then !generation of a random noise

       if (nrank==0) write(*,*) "Filled initial domain!"

       ux1=zero; uy1=zero; uz1=zero

       do k=1,xsize(3)
          z=real((k+xstart(3)-1-1),mytype)*dz
          do j=1,xsize(2)
             y=real((j+xstart(2)-1-1),mytype)*dy
             do i=1,xsize(1)
                x=real(i-1,mytype)*dx

                if (.not.tgv_twod) then
                   ux1(i,j,k)=+sin(x)*cos(y)*cos(z)
                   uy1(i,j,k)=-cos(x)*sin(y)*cos(z)
                   if (iscalar==1) then
                      phi1(i,j,k,1:numscalar)=sin(x)*sin(y)*cos(z)
                   endif
                else
                   ux1(i,j,k)=+sin(x)*cos(y)
                   uy1(i,j,k)=-cos(x)*sin(y)
                   if (iscalar==1) then
                      phi1(i,j,k,1:numscalar)=sin(x)*sin(y)
                   endif
                endif
                uz1(i,j,k)=zero
             enddo
          enddo
       enddo

       call random_seed(size=isize)
       allocate (seed(isize))
       seed(:)=67
       call random_seed(put=seed)
       !     call random_number(ux1)
       !     call random_number(uy1)
       ! call random_number(uz1)

       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                !              ux1(i,j,k)=noise*(ux1(i,j,k)-half)
                !              uy1(i,j,k)=noise*(uy1(i,j,k)-half)
                ! uz1(i,j,k)=0.05*(uz1(i,j,k)-half)
             enddo
          enddo
       enddo

       !     !modulation of the random noise
       !     do k=1,xsize(3)
       !        do j=1,xsize(2)
       !           if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/two
       !           if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
       !           um=exp(-0.2*y*y)
       !           do i=1,xsize(1)
       !              ux1(i,j,k)=um*ux1(i,j,k)
       !              uy1(i,j,k)=um*uy1(i,j,k)
       !              uz1(i,j,k)=um*uz1(i,j,k)
       !           enddo
       !        enddo
       !     enddo

    endif

    !  bxx1(j,k)=zero
    !  bxy1(j,k)=zero
    !  bxz1(j,k)=zero

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)!+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)!+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)!+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_tgv
  !********************************************************************

  subroutine boundary_conditions_tgv (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
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

  end subroutine boundary_conditions_tgv

  !********************************************************************

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: dxdydz, dxdz, x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

    call alloc_x(vol1, opt_global=.true.)
    call alloc_x(volSimps1, opt_global=.true.)

    vol1 = zero
    volSimps1 = zero

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico
    dxdydz=dx*dy*dz
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=xstart(1),xend(1)
             vol1(i,j,k)=dxdydz
             if (i .eq. 1   .or. i .eq. nx) vol1(i,j,k) = vol1(i,j,k)/two
             if (j .eq. 1   .or. j .eq. ny)  vol1(i,j,k) = vol1(i,j,k)/two
             if (k .eq. 1   .or. k .eq. nz)  vol1(i,j,k) = vol1(i,j,k)/two
          end do
       end do
    end do

    !X PENCILS !Utilizar para integral volumétrica dentro do domínio físico (método de Simpson)
    do k=xstart(3),xend(3)
       do j=xstart(2),xend(2)
          do i=xstart(1),xend(1)
             volSimps1(i,j,k)=dxdydz
             if (i .eq. 1   .or. i .eq. nx) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (j .eq. 1   .or. j .eq. ny) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (k .eq. 1   .or. k .eq. nz) volSimps1(i,j,k) = volSimps1(i,j,k) * (five/twelve)
             if (i .eq. 2   .or. i .eq. nx-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
             if (j .eq. 2   .or. j .eq. ny-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
             if (k .eq. 2   .or. k .eq. nz-1) volSimps1(i,j,k) = volSimps1(i,j,k) * (thirteen/twelve)
          end do
       end do
    end do

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocess_tgv(ux1,uy1,uz1,phi1,ep1)

    USE decomp_2d
    USE decomp_2d_io
    USE MPI
    USE var, only: nut1, srt_smag 
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1,ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

!    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,ep1,diss1
!    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
!    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    real(mytype) :: mp(numscalar),mps(numscalar),vl,es,es1,ek,ek1,ds,ds1
    real(mytype) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9

    real(mytype) :: eek, enst, eps, eps2
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

    integer :: i,j,k,is,code,nvect1
    character(len=30) :: filename

    if ((ivisu.ne.0).and.(mod(itime, ioutput).eq.0)) then
    !! Write vorticity as an example of post processing
    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    !y-derivatives
    call transpose_x_to_y(ux1,td2)
    call transpose_x_to_y(uy1,te2)
    call transpose_x_to_y(uz1,tf2)
    call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    !!z-derivatives
    call transpose_y_to_z(td2,td3)
    call transpose_y_to_z(te2,te3)
    call transpose_y_to_z(tf2,tf3)
    call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1=0.
    di1(:,:,:)=-0.5*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2)-&
        td1(:,:,:)*tb1(:,:,:)-&
        tg1(:,:,:)*tc1(:,:,:)-&
        th1(:,:,:)*tf1(:,:,:)
    uvisu=0.
    call fine_to_coarseV(1,di1,uvisu)
    994 format('./data/critq',I5.5)
    write(filename, 994) itime/ioutput
    call decomp_2d_write_one(1,uvisu,filename,2)
    endif

  end subroutine postprocess_tgv
  !############################################################################
  subroutine suspended(phi1,vol1,mp1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: vol1
    real(mytype),intent(out) :: mp1(1:numscalar)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temp1
    real(mytype) :: mp(1:numscalar)
    integer :: is,code

    mp=zero; mp1=zero

    do is=1, numscalar
       temp1 = phi1(:,:,:,is)*vol1(:,:,:)
       mp(is)= sum(temp1)
    end do

    call MPI_REDUCE(mp,mp1,numscalar,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    return
  end subroutine suspended
  !############################################################################
  subroutine dissipation (ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1)

    USE param
    USE variables
    USE decomp_2d
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,diss1
    real(mytype),dimension(3,3,xsize(1),xsize(2),xsize(3)) :: A
    integer :: k,j,i,m,l

    !INSTANTANEOUS DISSIPATION RATE
    diss1=zero
    A(:,:,:,:,:)=zero
    A(1,1,:,:,:)=ta1(:,:,:) !du/dx=ta1
    A(2,1,:,:,:)=tb1(:,:,:) !dv/dx=tb1
    A(3,1,:,:,:)=tc1(:,:,:) !dw/dx=tc1
    A(1,2,:,:,:)=td1(:,:,:) !du/dy=td1
    A(2,2,:,:,:)=te1(:,:,:) !dv/dy=te1
    A(3,2,:,:,:)=tf1(:,:,:) !dw/dy=tf1
    A(1,3,:,:,:)=tg1(:,:,:) !du/dz=tg1
    A(2,3,:,:,:)=th1(:,:,:) !dv/dz=th1
    A(3,3,:,:,:)=ti1(:,:,:) !dw/dz=ti1

    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             do m=1,3
                do l=1,3
                   diss1(i,j,k)=diss1(i,j,k)+two*xnu*half*half*(A(l,m,i,j,k)+A(m,l,i,j,k))**two
                enddo
             enddo
          enddo
       enddo
    enddo

    return

  end subroutine dissipation
  !############################################################################
  subroutine error_tgv2D(ux1,uy1,phi1)

    USE decomp_2d
    USE MPI
    USE param, only : one, two, xnu, ifirst, itime
    USE variables, only : numscalar, sc

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    ! Analytical error
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: errx, erry
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errs
    ! Time-discrete error
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: errxd1, erryd1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errsd1
    ! Space-time discrete error
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: errxd2, erryd2
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errsd2
    ! Temporary storage array
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tmp

    integer :: i,j,k,l,it
    real(mytype) :: x, y, xl1, xl2, xlinf, yl1, yl2, ylinf
    real(mytype) :: solx0, soly0, sols0
    real(mytype) :: solxt, solyt, solst, solxd, solyd, solsd, solxdd, solydd, solsdd
    real(mytype) :: xdamping(3), ydamping(3), sdamping(3,numscalar)

    ! Compute the errors
    call compute_tgv2D_errors(xdamping, ydamping, sdamping)

    ! Compute solutions and errors
    do k = 1,xsize(3)
      do j = 1,xsize(2)
        y = (j+xstart(2)-1-1)*dy
        do i = 1,xsize(1)
          x = (i+xstart(1)-1-1)*dx
          ! Initial solution
          solx0 = sin(x) * cos(y)
          soly0 = - cos(x) * sin(y)
          ! Analytical solution
          solxt = solx0 * xdamping(1)
          solyt = soly0 * ydamping(1)
          ! Time discrete solution
          solxd = solx0 * xdamping(2)
          solyd = soly0 * ydamping(2)
          ! Space-time discrete solution
          solxdd = solx0 * xdamping(3)
          solydd = soly0 * ydamping(3)
          ! Errors
          errx(i,j,k) = ux1(i,j,k) - solxt
          erry(i,j,k) = uy1(i,j,k) - solyt
          errxd1(i,j,k) = ux1(i,j,k) - solxd
          erryd1(i,j,k) = uy1(i,j,k) - solyd
          errxd2(i,j,k) = ux1(i,j,k) - solxdd
          erryd2(i,j,k) = uy1(i,j,k) - solydd
        enddo
      enddo
    enddo
    if (iscalar==1) then
    do l = 1, numscalar
      do k = 1,xsize(3)
        do j = 1,xsize(2)
          y = (j+xstart(2)-1-1)*dy
          do i = 1,xsize(1)
            x = (i+xstart(1)-1-1)*dx
            ! Initial solution
            sols0 = sin(x) * sin(y)
            ! Analytical solution
            solst = sols0 * sdamping(1,l)
            ! Time discrete solution
            solsd = sols0 * sdamping(2,l)
            ! Space-time discrete solution
            solsdd = sols0 * sdamping(3,l)
            ! Errors
            errs(i,j,k,l) = phi1(i,j,k,l) - solst
            errsd1(i,j,k,l) = phi1(i,j,k,l) - solsd
            errsd2(i,j,k,l) = phi1(i,j,k,l) - solsdd
          enddo
        enddo
      enddo
    enddo
    endif

    !
    ! Compute and print synthetic errors indicator
    !
    ! Magnitude of the solution
    call error_L1_L2_Linf(ux1, xl1, xl2, xlinf)
    call error_L1_L2_Linf(uy1, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "2D Taylor-Green errors:"
      write(*,*) "  Solution amplitude, L1, sqrt(L2), Linf:"
      write(*,*) "            ", xl1, xl2, xlinf
      write(*,*) "            ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = phi1(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif
    ! Analytical solution
    call error_L1_L2_Linf(errx, xl1, xl2, xlinf)
    call error_L1_L2_Linf(erry, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "  1: Analytical, L1, sqrt(L2), Linf:"
      write(*,*) "     XERROR ", xl1, xl2, xlinf
      write(*,*) "     YERROR ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = errs(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif
    ! Time discrete solution
    call error_L1_L2_Linf(errxd1, xl1, xl2, xlinf)
    call error_L1_L2_Linf(erryd1, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "  2: Time-discrete, L1, sqrt(L2), Linf:"                                 
      write(*,*) "     XERROR ", xl1, xl2, xlinf
      write(*,*) "     YERROR ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = errsd1(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif
    ! Space-time discrete solution
    call error_L1_L2_Linf(errxd2, xl1, xl2, xlinf)
    call error_L1_L2_Linf(erryd2, yl1, yl2, ylinf)
    if (nrank==0) then
      write(*,*) "  3: Time-space-discrete, L1, sqrt(L2), Linf:"
      write(*,*) "     XERROR ", xl1, xl2, xlinf
      write(*,*) "     YERROR ", yl1, yl2, ylinf
    endif
    if (iscalar==1) then
    do l = 1, numscalar
      tmp(:,:,:) = errsd2(:,:,:,l)
      call error_L1_L2_Linf(tmp, xl1, xl2, xlinf)
      if (nrank==0) then
        write(*,*) "            ", xl1, xl2, xlinf
      endif
    enddo
    endif

  end subroutine error_tgv2D

  ! Compute the damping factors
  subroutine compute_tgv2D_errors(xdamping, ydamping, sdamping)

    use decomp_2d
    use param, only : one, two, xnu, ifirst, itime, itimescheme
    use variables, only : numscalar, sc

    implicit none

    real(mytype), intent(out) :: xdamping(3), ydamping(3), sdamping(3, numscalar)

    integer :: it, l
    real(mytype) :: ktgv, k2tgv

    ! Compute modified wavenumber
    ktgv = one
    call compute_k2(ktgv, k2tgv)

    ! Init
    xdamping(:) = one
    ydamping(:) = one
    sdamping(:,:) = one

    ! Compute analytical damping
    do it = ifirst, itime
      xdamping(1) = xdamping(1) * exp(-two*dt*xnu)
      ydamping(1) = ydamping(1) * exp(-two*dt*xnu)
      do l = 1, numscalar
        sdamping(1,l) = sdamping(1,l) * exp(-two*dt*xnu/sc(l))
      enddo
    enddo

    ! Compute time and space-time discrete damping
    !
    ! Explicit Euler
    if (itimescheme.eq.1) then

      ! Time discrete errors
      do it = ifirst, itime
        xdamping(2) = xdamping(2) * (one - two*dt*xnu)
        ydamping(2) = ydamping(2) * (one - two*dt*xnu)
        do l = 1, numscalar
          sdamping(2,l) = sdamping(2,l) * (one - two*dt*xnu/sc(l))
        enddo
      enddo

      ! Space-time discrete errors
      do it = ifirst, itime
        xdamping(3) = xdamping(3) * (one - two*k2tgv*dt*xnu)
        ydamping(3) = ydamping(3) * (one - two*k2tgv*dt*xnu)
        do l = 1, numscalar
          sdamping(3,l) = sdamping(3,l) * (one - two*k2tgv*dt*xnu/sc(l))
        enddo
      enddo

    else

      if (nrank==0) print *, "TGV2D: No discrete error implemented for this time scheme."
      xdamping(2:) = xdamping(1)
      ydamping(2:) = ydamping(1)
      do l = 1, numscalar
        sdamping(2:,l) = sdamping(1,l)
      enddo

    endif

  end subroutine compute_tgv2D_errors

  ! Compute the modified wavenumber for the second derivative
  ! Warning : we use the X momentum wavenumber for Y momentum and for the scalars
  subroutine compute_k2(kin, k2out)

  USE decomp_2d, only : mytype
  USE param
  USE derivX, only : alsaix, asix, bsix, csix, dsix

  implicit none

  real(mytype), intent(in) :: kin
  real(mytype), intent(out) :: k2out

  if (kin.lt.zero .or. kin.gt.pi/min(dx,dy)) then
    if (nrank==0) then
      write(*,*) "TGV2D: Warning, incorrect wavenumber provided."
    endif
  endif

  k2out = asix * two * (one - cos(kin*dx)) &
        + four * bsix * half * (one - cos(two*kin*dx)) &
        + nine * csix * (two / nine) * (one - cos(three*kin*dx)) &
        + sixteen * dsix * (one / eight) * (one - cos(four*kin*dx))
  k2out = k2out / (one + two * alsaix * cos(kin*dx))

  end subroutine compute_k2

  ! Compute L1, L2 and Linf norm of given 3D array
  subroutine error_L1_L2_Linf(err, l1, l2, linf)

    USE decomp_2d
    USE MPI
    
    implicit none
      
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: err
    real(mytype),intent(out) :: l1, l2, linf

    integer :: i,j,k,code
    real(mytype) :: ll1, ll2, llinf, ntot

    ll1 = zero
    ll2 = zero
    llinf = zero
    ntot = nx*ny*nz

    do k = 1,xsize(3)
      do j = 1,xsize(2)
        do i = 1,xsize(1)
          ll1 = ll1 + abs(err(i,j,k))
          ll2 = ll2 + err(i,j,k)*err(i,j,k)
          llinf = max(llinf, abs(err(i,j,k)))
        enddo
      enddo
    enddo

    ! Parallel
    call MPI_ALLREDUCE(ll1,l1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(ll2,l2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(llinf,linf,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)

    ! Rescaling
    l1 = l1 / ntot
    l2 = sqrt(l2 / ntot)

  end subroutine error_L1_L2_Linf

end module tgv
