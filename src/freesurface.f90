module freesurface

  use decomp_2d, only : mytype, xsize
  
  implicit none
  
  real(mytype), save, allocatable, dimension(:,:,:) :: S1

  private
  public :: reinit_ls, update_fluid_properties

contains

  !!------------------------------------------------------------------------
  !!  SUBROUTINE: reinit_ls
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Reinitialise the level-set distance function by computing
  !!              the steady-state of
  !!                d(ls)/dtau = S(ls_0) (1 - |grad(ls)|)
  !!              where ls_0 is the initially computed level-set, tau is a
  !!              pseudo timestep and S is a smoothing function
  !!                S(ls_0) = ls_0 / sqrt({ls_0}^2 + (|grad(ls)|dx)^2)
  !!------------------------------------------------------------------------
  subroutine reinit_ls(levelset1, neps, init)

    use MPI

    use decomp_2d, only : mytype, real_type, xsize, xstart, xend, ysize, zsize, nrank
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : alloc_x
    use param, only : zero, one, two, three, four, six, ten, sixteen
    use param, only : dt, dx, dy, dz
    use var, only : itime
    use var, only : xlx, yly, zlz
    use var, only : nclx1, ncly1, nclz1

    !! Derivatives
    use weno, only : weno5

    implicit none

    !! In
    integer, intent(in) :: neps
    logical, intent(in) :: init

    !! InOut
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: levelset1

    !! Local
    integer :: i, j, k

    integer :: iter, subiter
    logical :: converged
    integer :: ierr

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mag_grad_ls1, gradx_ls1, grady_ls1, &
         gradz_ls1, levelset1_old
    real(mytype), dimension(3) :: a, b, c

    real(mytype) :: dtau
    real(mytype) :: delta_ls, global_delta_ls, max_delta_ls, global_max_delta_ls, abs_delta_ls
    real(mytype) :: deltax
    integer :: ctr, global_ctr
    real(mytype) :: alpha

    real(mytype) :: eps
    integer :: niter

    real(mytype), parameter :: cfl = one / ten
    real(mytype) :: tol, global_tol

    if (init) then
       call alloc_x(S1)
       S1(:,:,:) = levelset1(:,:,:)
    endif

    deltax = min(dx, dy, dz)
    eps = deltax !! SUSSMAN1994
    alpha = sixteen / ten * eps

    dtau = cfl * deltax
    if (.not.init) then
       niter = int((float(neps) * alpha) / dtau)
       tol = alpha
    else
       niter = int(max(xlx, yly, zlz) / dtau)
       niter = int(one / dtau)
       tol = max(maxval(levelset1), abs(minval(levelset1)))
       call MPI_ALLREDUCE(tol,global_tol,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
       tol = global_tol
    endif

    a(1) = one
    a(2) = three / four
    a(3) = one / three

    b(1) = zero
    b(2) = one / four
    b(3) = two / three

    c(1) = one
    c(2) = one / four
    c(3) = two / three

    !! Step to steady state
    converged = .false.
    iter = 0
    dtau = dt
    do while (.not.converged)

       !! Update level-set
       levelset1_old(:,:,:) = levelset1(:,:,:)

       do subiter = 1, 3
          call reinit_ls_grad(gradx_ls1, grady_ls1, gradz_ls1, levelset1, S1)
          mag_grad_ls1 = sqrt(gradx_ls1**2 + grady_ls1**2 + gradz_ls1**2)
          call compute_reinit_smooth(S1, levelset1, mag_grad_ls1, alpha, eps)

          levelset1(:,:,:) = a(subiter) * levelset1_old(:,:,:) + b(subiter) * levelset1(:,:,:) &
               + c(subiter) * dtau * S1(:,:,:) * (one - mag_grad_ls1(:,:,:))

          !! Bounding at the boundaries
          if (nclx1.ne.0) then
             levelset1(1,:,:) = levelset1(2,:,:)
             levelset1(xsize(1),:,:) = levelset1(xsize(1)-1,:,:)
          endif
          if (ncly1.ne.0) then
             if (xstart(2).eq.1) then
                levelset1(:,1,:) = levelset1(:,2,:)
             endif
             if (xend(2).eq.ysize(2)) then
                levelset1(:,xsize(2),:) = levelset1(:,xsize(2)-1,:)
             endif
          endif
          if (nclz1.ne.0) then
             if (xstart(3).eq.1) then
                levelset1(:,:,1) = levelset1(:,:,2)
             endif
             if (xend(3).eq.zsize(3)) then
                levelset1(:,:,xsize(3)) = levelset1(:,:,xsize(3)-1)
             endif
          endif
       enddo

       !! Test convergence (SUSSMAN1994)
       delta_ls = zero
       max_delta_ls = zero
       ctr = 0
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                abs_delta_ls = abs(levelset1(i, j, k) - levelset1_old(i, j, k))
                max_delta_ls = max(max_delta_ls, abs_delta_ls)
                if (abs(levelset1_old(i, j, k)).le.tol) then
                   delta_ls = delta_ls + abs_delta_ls
                   ctr = ctr + 1
                endif
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(delta_ls,global_delta_ls,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(max_delta_ls,global_max_delta_ls,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(ctr,global_ctr,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

       if ((iter.ge.1).or.(itime.gt.1)) then 
          if ((global_ctr==0).or.(global_delta_ls < dtau * (eps**2) * float(global_ctr))) then
             !! SUSSMAN TEST
             converged = .true.
          else if (global_max_delta_ls.lt.0.0005) then
             !! HYDRO3D TEST
             converged = .true.
          endif
       endif

       iter = iter + 1

       if (iter==niter) then
          converged = .true.
       endif

       if (nrank.eq.0) then
          print *, "Level-set reinitialisation: ", iter, " / ", niter, " : ", global_delta_ls, &
               " / ", dtau * (eps**3) * float(global_ctr), dtau, eps**2
       endif
    enddo

  end subroutine reinit_ls

  subroutine reinit_ls_grad (gradx_ls1, grady_ls1, gradz_ls1, levelset1, S1)

    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

    use weno, only : weno5

    use var, only : zero, one
    use var, only : wx1 => ta1, &
         gradx_ls1m => tb1, gradx_ls1p => tc1, grady_ls1m => td1, grady_ls1p => te1, &
         gradz_ls1m => tf1, gradz_ls1p => tg1
    use var, only : levelset2 => ta2, wy2 => tb2, grady_ls2 => tc2, gradz_ls2 => td2
    use var, only : levelset3 => ta3, wz3 => tb3, gradz_ls3 => tc3
    use var, only : nclx1, ncly1, nclz1
    
    use variables, only : nz

    implicit none

    !! In
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: S1

    !! Out
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: gradx_ls1, grady_ls1, gradz_ls1

    !! Locals
    integer :: i, j, k
    integer :: fsbc

    call transpose_x_to_y(levelset1, levelset2)
    call transpose_y_to_z(levelset2, levelset3)

    ! Z
    if (nz.gt.1) then
       if (nclz1.eq.0) then
          fsbc = 0
       else
          fsbc = 2
       endif
          
       wz3(:,:,:) = one
       call weno5 (gradz_ls3, levelset3, wz3, 3, fsbc, fsbc, zsize(1), zsize(2), zsize(3), 1)
       call transpose_z_to_y(gradz_ls3, gradz_ls2)
       call transpose_y_to_x(gradz_ls2, gradz_ls1m)

       wz3(:,:,:) = -one
       call weno5 (gradz_ls3, levelset3, wz3, 3, fsbc, fsbc, zsize(1), zsize(2), zsize(3), 1)
       call transpose_z_to_y(gradz_ls3, gradz_ls2)
       call transpose_y_to_x(gradz_ls2, gradz_ls1p)
    else
       gradz_ls1m(:,:,:) = zero
       gradz_ls1p(:,:,:) = zero
    endif

    ! Y
    if (ncly1.eq.0) then
       fsbc = 0
    else
       fsbc = 2
    endif
    
    wy2(:,:,:) = one
    call weno5 (grady_ls2, levelset2, wy2, 2, fsbc, fsbc, ysize(1), ysize(2), ysize(3), 1)
    call transpose_y_to_x(grady_ls2, grady_ls1m)

    wy2(:,:,:) = -one
    call weno5 (grady_ls2, levelset2, wy2, 2, fsbc, fsbc, ysize(1), ysize(2), ysize(3), 1)
    call transpose_y_to_x(grady_ls2, grady_ls1p)

    ! X
    if (nclx1.eq.0) then
       fsbc = 0
    else
       fsbc = 2
    endif
    
    wx1 = one
    call weno5 (gradx_ls1m, levelset1, wx1, 1, fsbc, fsbc, xsize(1), xsize(2), xsize(3), 1)
    wx1 = -one
    call weno5 (gradx_ls1p, levelset1, wx1, 1, fsbc, fsbc, xsize(1), xsize(2), xsize(3), 1)
    
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             !! X
             if ((S1(i,j,k) * gradx_ls1m(i,j,k) > zero) &
                  .and. (S1(i,j,k) * gradx_ls1p(i,j,k) > -S1(i,j,k) * gradx_ls1m(i,j,k))) then
                gradx_ls1(i,j,k) = gradx_ls1m(i,j,k)
             elseif ((S1(i,j,k) * gradx_ls1p(i,j,k) < zero) &
                  .and. (S1(i,j,k) * gradx_ls1m(i,j,k) < -S1(i,j,k) * gradx_ls1p(i,j,k))) then
                gradx_ls1(i,j,k) = gradx_ls1p(i,j,k)
             else
                gradx_ls1(i,j,k) = zero
             endif

             !! Y
             if ((S1(i,j,k) * grady_ls1m(i,j,k) > zero) &
                  .and. (S1(i,j,k) * grady_ls1p(i,j,k) > -S1(i,j,k) * grady_ls1m(i,j,k))) then
                grady_ls1(i,j,k) = grady_ls1m(i,j,k)
             elseif ((S1(i,j,k) * grady_ls1p(i,j,k) < zero) &
                  .and. (S1(i,j,k) * grady_ls1m(i,j,k) < -S1(i,j,k) * grady_ls1p(i,j,k))) then
                grady_ls1(i,j,k) = grady_ls1p(i,j,k)
             else
                grady_ls1(i,j,k) = zero
             endif

             !! Z
             if ((S1(i,j,k) * gradz_ls1m(i,j,k) > zero) &
                  .and. (S1(i,j,k) * gradz_ls1p(i,j,k) > -S1(i,j,k) * gradz_ls1m(i,j,k))) then
                gradz_ls1(i,j,k) = gradz_ls1m(i,j,k)
             elseif ((S1(i,j,k) * gradz_ls1p(i,j,k) < zero) &
                  .and. (S1(i,j,k) * gradz_ls1m(i,j,k) < -S1(i,j,k) * gradz_ls1p(i,j,k))) then
                gradz_ls1(i,j,k) = gradz_ls1p(i,j,k)
             else
                gradz_ls1(i,j,k) = zero
             endif
          enddo
       enddo
    enddo

  endsubroutine reinit_ls_grad

  subroutine compute_reinit_smooth (S1, levelset1, mag_grad_ls1, alpha, eps)

    use decomp_2d, only : mytype, xsize

    use var, only : zero, one, two

    implicit none

    !! In
    real(mytype), intent(in) :: eps, alpha
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1, mag_grad_ls1

    !! Out
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: S1

    !! Local
    integer :: i, j, k

    ! !! SUSSMAN1994
    ! S1(:,:,:) = levelset1(:,:,:) / sqrt(levelset1(:,:,:)**2 + eps**2)

    !! MCSHERRY2017
    S1(:,:,:) = levelset1(:,:,:) &
         / (sqrt(levelset1(:,:,:)**2 + ((mag_grad_ls1(:,:,:))**2 * eps**2)))

    ! do k = 1, xsize(3)
    !    do j = 1, xsize(2)
    !       do i = 1, xsize(1)
    !          if (abs(levelset1(i, j, k)).gt.alpha) then
    !             S1(i, j, k) = zero
    !          endif
    !       enddo
    !    enddo
    ! enddo

  endsubroutine compute_reinit_smooth

  subroutine update_fluid_properties(rho1, mu1, phi1)

    use decomp_2d, only : mytype, xsize
    use var, only : half, one, two, three, pi, ten, sixteen
    use var, only : dx, dy, dz
    use var, only : numscalar
    use param, only : nrhotime, ilevelset
    use param, only : dens1, dens2, visc1, visc2

    implicit none

    !! Input
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1

    !! InOut
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    !! Local
    integer :: i, j, k
    real(mytype) :: eps

    eps = (sixteen / ten) * (dx * dy * dz)**(one / three)

    !!---------------------------------
    if (ilevelset.gt.0) then
       !! Use levelset to set density field
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                if (phi1(i, j, k, ilevelset).gt.eps) then
                   !! Fluid 1
                   rho1(i, j, k, 1) = dens1
                   mu1(i, j, k) = visc1
                else if (phi1(i, j, k, ilevelset).lt.-eps) then
                   !! Fluid 2
                   rho1(i, j, k, 1) = dens2
                   mu1(i, j, k) = visc2
                else
                   !! Interface: smooth properties
                   ! rho1(i, j, k, 1) = ((dens1 + dens2) &
                   !      + (dens1 - dens2) * sin(pi * phi1(i, j, k, ilevelset) / (two * eps))) &
                   !      / (two * dens1)
                   ! mu1(i, j, k) = ((visc1 + visc2) &
                   !      + (visc1 - visc2) * sin(pi * phi1(i, j, k, ilevelset) / (two * eps))) &
                   !      / (two * visc1)
                   rho1(i, j, k, 1) = dens2 + (dens1 - dens2) * half * (one + phi1(i, j, k, ilevelset) / eps &
                        + (sin(pi * phi1(i, j, k, ilevelset) / eps)) / pi)
                   mu1(i, j, k) = visc2 + (visc1 - visc2) * half * (one + phi1(i, j, k, ilevelset) / eps &
                        + (sin(pi * phi1(i, j, k, ilevelset) / eps)) / pi)
                endif
             enddo
          enddo
       enddo
    endif

  endsubroutine update_fluid_properties


  subroutine surface_tension_force(dux1, duy1, duz1, rho1, levelset1)

    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

    use variables, only : cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : sx, sy, sz

    use var, only : zero, one, two, three, pi
    use var, only : stfx1 => ta1, stfy1 => tb1, stfz1 => tc1, &
         gradx_ls1 => td1, grady_ls1 => te1, gradz_ls1 => tf1, mag_grad_ls1 => tg1, &
         kdelta => th1
    use var, only : levelset2 => ta2, grady_ls2 => tb2, gradz_ls2 => tc2
    use var, only : levelset3 => ta3, gradz_ls3 => tb3
    use var, only : ux1, uy2, uz3
    use var, only : ep1, drho1, divu3
    use var, only : dv3, ppi3, dip3
    use var, only : pp2, ppi2, dip2
    use var, only : pp1, ppi1 => ti1, di1
    use var, only : ph2, ph3
    use var, only : nxmsize, nymsize, nzmsize

    use div, only : divergence

    use param, only : dx, dy, dz
    use param, only : nclx1, nclxn, ncly1, nclyn, nclz1, nclzn
    use param, only : nrhotime
    use param, only : interface_thickness

    use weno, only : weno5
    
    implicit none

    !! IN
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1

    !! INOUT
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: dux1, duy1, duz1

    !! LOCAL
    real(mytype), parameter :: macheps = epsilon(1._mytype)
    integer :: nlock
    integer :: i, j, k
    real(mytype) :: alpha

    alpha = interface_thickness * (dx * dy * dz)**(one / three)

    !! Fist, compute grad(ls) and normalise
    call transpose_x_to_y(levelset1, levelset2)
    call transpose_y_to_z(levelset2, levelset3)

    call weno5 (gradz_ls3, levelset3, uz3, 3, nclz1, nclzn, zsize(1), zsize(2), zsize(3), 1)
    call weno5 (grady_ls2, levelset2, uy2, 2, ncly1, nclyn, ysize(1), ysize(2), ysize(3), 1)
    call weno5 (gradx_ls1, levelset1, ux1, 1, nclx1, nclxn, xsize(1), xsize(2), xsize(3), 1)

    call transpose_z_to_y(gradz_ls3, gradz_ls2)
    call transpose_y_to_x(grady_ls2, grady_ls1)
    call transpose_y_to_x(gradz_ls2, gradz_ls1)
    
    mag_grad_ls1(:,:,:) = sqrt(gradx_ls1**2 + grady_ls1**2 + gradz_ls1**2) + macheps
    gradx_ls1(:,:,:) = gradx_ls1(:,:,:) / mag_grad_ls1(:,:,:)
    grady_ls1(:,:,:) = grady_ls1(:,:,:) / mag_grad_ls1(:,:,:)
    gradz_ls1(:,:,:) = gradz_ls1(:,:,:) / mag_grad_ls1(:,:,:)

    !! Compute div(grad(ls) / |grad(ls)|)
    nlock = 0
    call divergence(dv3,rho1,gradx_ls1,grady_ls1,gradz_ls1,ep1,drho1,divu3,nlock)
    call interzpv(ppi3,dv3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    call transpose_z_to_y(ppi3,pp2,ph3)
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call transpose_y_to_x(ppi2,pp1,ph2)
    call interxpv(ppi1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)

    !! Recompute grad(ls) - memory will have been stomped on :(
    call transpose_x_to_y(levelset1, levelset2)
    call transpose_y_to_z(levelset2, levelset3)

    call weno5 (gradz_ls3, levelset3, uz3, 3, nclz1, nclzn, zsize(1), zsize(2), zsize(3), 1)
    call weno5 (grady_ls2, levelset2, uy2, 2, ncly1, nclyn, ysize(1), ysize(2), ysize(3), 1)
    call weno5 (gradx_ls1, levelset1, ux1, 1, nclx1, nclxn, xsize(1), xsize(2), xsize(3), 1)

    call transpose_z_to_y(gradz_ls3, gradz_ls2)
    call transpose_y_to_x(grady_ls2, grady_ls1)
    call transpose_y_to_x(gradz_ls2, gradz_ls1)

    !! Compute forcing
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             if (abs(levelset1(i, j, k)).lt.alpha) then
                kdelta(i, j, k) = (one / (two * alpha)) * (one + cos(pi * levelset1(i, j, k) / alpha))
             else
                kdelta(i, j, k) = zero
             endif
          enddo
       enddo
    enddo
    stfx1(:,:,:) = ppi1(:,:,:) * kdelta(:,:,:) * gradx_ls1(:,:,:)
    stfy1(:,:,:) = ppi1(:,:,:) * kdelta(:,:,:) * grady_ls1(:,:,:)
    stfz1(:,:,:) = ppi1(:,:,:) * kdelta(:,:,:) * gradz_ls1(:,:,:)
    
    !! Add contribution to forcing terms
    dux1(:,:,:) = dux1(:,:,:) + stfx1(:,:,:) / rho1(:,:,:,1)
    duy1(:,:,:) = duy1(:,:,:) + stfy1(:,:,:) / rho1(:,:,:,1)
    duz1(:,:,:) = duz1(:,:,:) + stfz1(:,:,:) / rho1(:,:,:,1)
    
  end subroutine surface_tension_force

endmodule freesurface