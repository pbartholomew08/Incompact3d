module freesurface

  use decomp_2d, only : mytype, xsize
  
  implicit none
  
  real(mytype), save, allocatable, dimension(:,:,:) :: nx1, ny1, nz1
  real(mytype), save, allocatable, dimension(:,:,:) :: S1
  real(mytype), save, allocatable, dimension(:,:,:) :: wx1, wy1, wz1

  private
  public :: reinit_ls, update_fluid_properties, &
       surface_tension_force

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
    use param, only : zero, one, two, three, four, ten, sixteen
    use param, only : dt, dx, dy, dz
    use var, only : nclx1, ncly1, nclz1
    use var, only : interface_thickness

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
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mask
    real(mytype), dimension(3) :: a, b, c

    real(mytype) :: dtau
    real(mytype) :: delta_ls, global_delta_ls, max_delta_ls, global_max_delta_ls, abs_delta_ls
    real(mytype) :: deltax
    integer :: ctr, global_ctr
    real(mytype) :: alpha

    real(mytype) :: eps
    integer :: maxiter, miniter

    real(mytype), parameter :: cfl = one / ten
    real(mytype) :: maxdist, global_maxdist

    if (init) then
       call alloc_x(S1)
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                if (levelset1(i, j, k).lt.zero) then
                   S1(i, j, k) = -one
                else if (levelset1(i, j, k).gt.zero) then
                   S1(i, j, k) = one
                else
                   S1(i, j, k) = zero
                endif
             end do
          end do
       end do

       call init_wls(levelset1, S1)
    endif

    deltax = min(dx, dy, dz)
    eps = deltax !! SUSSMAN1994
    alpha = interface_thickness * deltax

    dtau = cfl * deltax
    if (.not.init) then
       maxdist = alpha
       maxiter = int((float(neps) * maxdist) / dtau)
       miniter = 0
    else
       maxdist = zero
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                if (levelset1(i, j, k).gt.zero) then
                   maxdist = maxdist + dx * dy * dz
                endif
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(maxdist,global_maxdist,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       global_maxdist = global_maxdist**(one / three)
       maxdist = max(global_maxdist, alpha, maxval(levelset1), abs(minval(levelset1)))
       miniter = int(float(neps) * alpha / dtau)
       maxiter = int(float(neps) * maxdist / dtau)
    endif

    !! Set the mask
    !  XXX It is only necessary to solve reinitialisation in the vicinity of the interface.
    mask(:,:,:) = one
    if (.not.init) then
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                if (abs(levelset1(i, j, k)).gt.maxdist) then
                   !! Point is far enough away from interface to ignore
                   mask(i, j, k) = zero
                endif
             enddo
          enddo
       enddo
    endif

    a(1) = one ! one
    a(2) = one ! three / four
    a(3) = zero ! one / three

    b(1) = zero ! zero
    b(2) = zero ! one / four
    b(3) = one ! two / three

    c(1) = three / four ! one
    c(2) = one ! one / four
    c(3) = zero ! two / three

    !! Step to steady state
    converged = .false.
    iter = 0
    dtau = dt
    gradx_ls1(:,:,:) = zero; grady_ls1(:,:,:) = zero; gradz_ls1(:,:,:) = zero
    do while (.not.converged)

       !! Update level-set
       levelset1_old(:,:,:) = levelset1(:,:,:)

       do subiter = 1, 2
          call reinit_ls_grad(gradx_ls1, grady_ls1, gradz_ls1, levelset1)
          mag_grad_ls1(:,:,:) = sqrt((gradx_ls1(:,:,:))**2 &
               + (grady_ls1(:,:,:))**2 &
               + (gradz_ls1(:,:,:))**2)
          call compute_reinit_smooth(levelset1, mag_grad_ls1, eps)

          levelset1(:,:,:) = a(subiter) * levelset1_old(:,:,:) &
               + b(subiter) * levelset1(:,:,:) &
               + mask(:,:,:) * (c(subiter) * dtau * S1(:,:,:) * (one - mag_grad_ls1(:,:,:)))

          !! Update pseudo-velocities
          wx1(:,:,:) = S1(:,:,:) * gradx_ls1(:,:,:) / (mag_grad_ls1(:,:,:) + epsilon(eps))
          wy1(:,:,:) = S1(:,:,:) * grady_ls1(:,:,:) / (mag_grad_ls1(:,:,:) + epsilon(eps))
          wz1(:,:,:) = S1(:,:,:) * gradz_ls1(:,:,:) / (mag_grad_ls1(:,:,:) + epsilon(eps))

          ! !! Bounding at the boundaries
          ! if (nclx1.ne.0) then
          !    i = 1
          !    levelset1(i,:,:) = levelset1(i + 1,:,:)
          !    i = xsize(1)
          !    levelset1(i,:,:) = levelset1(i - 1,:,:)
          ! endif
          ! if (ncly1.ne.0) then
          !    if (xstart(2).eq.1) then
          !       j = 1
          !       levelset1(:, j,:) = levelset1(:, j + 1,:)
          !    endif
          !    if (xend(2).eq.ysize(2)) then
          !       j = xsize(2)
          !       levelset1(:, j,:) = levelset1(:, j - 1,:)
          !    endif
          ! endif
          ! if (nclz1.ne.0) then
          !    if (xstart(3).eq.1) then
          !       k = 1
          !       levelset1(:,:, k) = levelset1(:,:, k + 1)
          !    endif
          !    if (xend(3).eq.zsize(3)) then
          !       levelset1(:,:, k) = levelset1(:,:, k - 1)
          !    endif
          ! endif
       enddo !! End time sub-iterations

       !! Test convergence (SUSSMAN1994)
       delta_ls = zero
       max_delta_ls = zero
       ctr = 0
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                abs_delta_ls = abs(levelset1(i, j, k) - levelset1_old(i, j, k))
                max_delta_ls = max(max_delta_ls, abs_delta_ls)
                if (abs(levelset1_old(i, j, k)).le.maxdist) then
                   delta_ls = delta_ls + abs_delta_ls
                   ctr = ctr + 1
                endif
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(delta_ls,global_delta_ls,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(max_delta_ls,global_max_delta_ls,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(ctr,global_ctr,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

       if (iter.ge.miniter) then 
          if ((global_ctr==0).or.(global_delta_ls < dtau * (eps**2) * float(global_ctr))) then
             !! SUSSMAN TEST
             converged = .true.
          else if (global_max_delta_ls.lt.0.0005) then
             !! HYDRO3D TEST
             converged = .true.
          endif
       endif

       iter = iter + 1

       if (iter==maxiter) then
          converged = .true.
       elseif (init) then
          converged = .false.
       elseif (iter.lt.miniter) then
          converged = .false.
       endif

       if (nrank.eq.0) then
          print *, "Level-set reinitialisation: ", iter, " / ", maxiter, "(", miniter, ")", " ; ", &
               global_delta_ls, " / ", dtau * (eps**2) * float(global_ctr), " ; ", &
               global_max_delta_ls, " / 0.0005"
       endif
    enddo

  end subroutine reinit_ls

  subroutine reinit_ls_grad (gradx_ls1, grady_ls1, gradz_ls1, levelset1)

    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

    use weno, only : weno5

    use var, only : zero, one, two
    use var, only : gradx_ls1m => tb1, gradx_ls1p => tc1, grady_ls1m => td1, grady_ls1p => te1, &
         gradz_ls1m => tf1, gradz_ls1p => tg1
    use var, only : levelset2 => ta2, wy2 => tb2, grady_ls2 => tc2, gradz_ls2 => td2, wz2 => te2
    use var, only : levelset3 => ta3, wz3 => tb3, gradz_ls3 => tc3
    use var, only : nclx1, ncly1, nclz1
    
    use variables, only : nz
    use var, only : dx, dy

    implicit none

    !! In
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1

    !! Out
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: gradx_ls1, grady_ls1, gradz_ls1

    !! Locals
    integer :: i, j, k
    integer :: fsbc

    call transpose_x_to_y(levelset1, levelset2)
    call transpose_y_to_z(levelset2, levelset3)

    call transpose_x_to_y(wy1, wy2)
    call transpose_x_to_y(wz1, wz2)
    call transpose_y_to_z(wz2, wz3)

    ! Z
    if (nz.gt.1) then
       if (nclz1.eq.0) then
          fsbc = 0
       else
          fsbc = 2
       endif
          
       ! wz3(:,:,:) = one
       ! call weno5 (gradz_ls3, levelset3, wz3, 3, fsbc, fsbc, zsize(1), zsize(2), zsize(3), 1)
       ! call transpose_z_to_y(gradz_ls3, gradz_ls2)
       ! call transpose_y_to_x(gradz_ls2, gradz_ls1m)

       ! wz3(:,:,:) = -one
       ! call weno5 (gradz_ls3, levelset3, wz3, 3, fsbc, fsbc, zsize(1), zsize(2), zsize(3), 1)
       ! call transpose_z_to_y(gradz_ls3, gradz_ls2)
       ! call transpose_y_to_x(gradz_ls2, gradz_ls1p)

       call weno5 (gradz_ls3, levelset3, wz3, 3, fsbc, fsbc, zsize(1), zsize(2), zsize(3), 1)
       call transpose_z_to_y(gradz_ls3, gradz_ls2)
       call transpose_y_to_x(gradz_ls2, gradz_ls1)
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

    ! wy2(:,:,:) = one
    ! call weno5 (grady_ls2, levelset2, wy2, 2, fsbc, fsbc, ysize(1), ysize(2), ysize(3), 1)
    ! call transpose_y_to_x(grady_ls2, grady_ls1m)
       
    ! wy2(:,:,:) = -one
    ! call weno5 (grady_ls2, levelset2, wy2, 2, fsbc, fsbc, ysize(1), ysize(2), ysize(3), 1)
    ! call transpose_y_to_x(grady_ls2, grady_ls1p)

    call weno5 (grady_ls2, levelset2, wy2, 2, fsbc, fsbc, ysize(1), ysize(2), ysize(3), 1)
    call transpose_y_to_x(grady_ls2, grady_ls1)

    ! X
    if (nclx1.eq.0) then
       fsbc = 0
    else
       fsbc = 2
    endif

    ! wx1 = one
    ! call weno5 (gradx_ls1m, levelset1, wx1, 1, fsbc, fsbc, xsize(1), xsize(2), xsize(3), 1)
    ! wx1 = -one
    ! call weno5 (gradx_ls1p, levelset1, wx1, 1, fsbc, fsbc, xsize(1), xsize(2), xsize(3), 1)

    call weno5 (gradx_ls1, levelset1, wx1, 1, fsbc, fsbc, xsize(1), xsize(2), xsize(3), 1)

    ! do k = 1, xsize(3)
    !    do j = 1, xsize(2)
    !       do i = 1, xsize(1)
    !          !! X
    !          if ((S1(i,j,k) * gradx_ls1m(i,j,k) > zero) &
    !               .and. (S1(i,j,k) * gradx_ls1p(i,j,k) > -S1(i,j,k) * gradx_ls1m(i,j,k))) then
    !             gradx_ls1(i,j,k) = gradx_ls1m(i,j,k)
    !          else if ((S1(i,j,k) * gradx_ls1p(i,j,k) < zero) &
    !               .and. (S1(i,j,k) * gradx_ls1m(i,j,k) < -S1(i,j,k) * gradx_ls1p(i,j,k))) then
    !             gradx_ls1(i,j,k) = gradx_ls1p(i,j,k)
    !          else
    !             gradx_ls1(i,j,k) = zero
    !          endif

    !          !! Y
    !          if ((S1(i,j,k) * grady_ls1m(i,j,k) > zero) &
    !               .and. (S1(i,j,k) * grady_ls1p(i,j,k) > -S1(i,j,k) * grady_ls1m(i,j,k))) then
    !             grady_ls1(i,j,k) = grady_ls1m(i,j,k)
    !          else if ((S1(i,j,k) * grady_ls1p(i,j,k) < zero) &
    !               .and. (S1(i,j,k) * grady_ls1m(i,j,k) < -S1(i,j,k) * grady_ls1p(i,j,k))) then
    !             grady_ls1(i,j,k) = grady_ls1p(i,j,k)
    !          else
    !             grady_ls1(i,j,k) = zero
    !          endif

    !          !! Z
    !          if ((S1(i,j,k) * gradz_ls1m(i,j,k) > zero) &
    !               .and. (S1(i,j,k) * gradz_ls1p(i,j,k) > -S1(i,j,k) * gradz_ls1m(i,j,k))) then
    !             gradz_ls1(i,j,k) = gradz_ls1m(i,j,k)
    !          elseif ((S1(i,j,k) * gradz_ls1p(i,j,k) < zero) &
    !               .and. (S1(i,j,k) * gradz_ls1m(i,j,k) < -S1(i,j,k) * gradz_ls1p(i,j,k))) then
    !             gradz_ls1(i,j,k) = gradz_ls1p(i,j,k)
    !          else
    !             gradz_ls1(i,j,k) = zero
    !          endif
    !       enddo
    !    enddo
    ! enddo

  endsubroutine reinit_ls_grad

  subroutine compute_reinit_smooth (levelset1, mag_grad_ls1, eps)

    use decomp_2d, only : mytype, xsize

    implicit none

    !! In
    real(mytype), intent(in) :: eps
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1, mag_grad_ls1

    !! Local
    real(mytype), parameter :: macheps = epsilon(eps)

    ! !! SUSSMAN1994
    ! S1(:,:,:) = levelset1(:,:,:) / sqrt(levelset1(:,:,:)**2 + eps**2)

    !! MCSHERRY2017
    S1(:,:,:) = levelset1(:,:,:) &
         / (sqrt(levelset1(:,:,:)**2 + (mag_grad_ls1(:,:,:) * eps)**2) + macheps)

  endsubroutine compute_reinit_smooth

  subroutine update_fluid_properties(rho1, mu1, levelset1)

    use decomp_2d, only : mytype, xsize, ysize, zsize
    use var, only : half, one, two, three, pi, ten, sixteen
    use var, only : dx, dy, dz
    use var, only : numscalar
    use var, only : interface_thickness
    use param, only : nrhotime, ilevelset
    use param, only : dens1, dens2, visc1, visc2
    use var, only : nclx1, ncly1, nclz1, nclxn, nclyn, nclzn

    use var, only : lsp1 => pgy1, lspv21 => pp1
    use var, only : lsp12 => uyp2, lsp2 => upi2, lspv32 => pp2, lspv2 => ppi2
    use var, only : lsp23 => uzp3, lsp3 => po3, lspv3 => ppi3, prop3 => dv3

    use var, only : di1, sx, cifxp6, cisxp6, ciwxp6, nxmsize
    use var, only : dipp2, sy, cifyp6, cisyp6, ciwyp6, nymsize
    use var, only : dipp3, sz, cifzp6, ciszp6, ciwzp6, nzmsize
    use var, only : dip3, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use var, only : dip2, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use var, only : cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x

    use var, only : ph1, ph2, ph3, ph4

    implicit none

    !! Input
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1

    !! InOut
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    !! Local
    integer :: i, j, k
    real(mytype) :: deltax, alpha
    real(mytype) :: hk, k99

    deltax = min(dx, dy, dz)
    alpha = interface_thickness * deltax

    k99 = -(one / (two * alpha)) * log(one / 0.99_mytype - one)

    !! Interpolate levelset to pressure nodes
    if ((nclx1.ne.0).and.(nclxn.ne.0)) then
       do i = 1, nxmsize
          lsp1(i,:,:) = 0.5 * (levelset1(i + 1,:,:) + levelset1(i,:,:))
       enddo
    else
       do i = 1, nxmsize - 1
          lsp1(i,:,:) = 0.5 * (levelset1(i + 1,:,:)) + levelset1(i,:,:)
       enddo
       lsp1(nxmsize,:,:) = 0.5 * (levelset1(1,:,:) + levelset1(nxmsize,:,:))
    endif
    call transpose_x_to_y(lsp1,lsp12,ph4)

    if ((ncly1.ne.0).and.(nclyn.ne.0)) then
       do j = 1, nymsize
          lsp2(:,j,:) = 0.5 * (lsp12(:,j+1,:) + lsp12(:,j,:))
       enddo
    else
       do j = 1, nymsize - 1
          lsp2(:,j,:) = 0.5 * (lsp12(:,j+1,:) + lsp12(:,j,:))
       enddo
       lsp2(:,nymsize,:) = 0.5 * (lsp12(:,1,:) + lsp12(:,nymsize,:))
    endif
    call transpose_y_to_z(lsp2,lsp23,ph3)

    if ((nclz1.ne.0).and.(nclzn.ne.0)) then
       do k = 1, nzmsize
          lsp3(:,:,k) = 0.5 * (lsp23(:,:,k+1) + lsp23(:,:,k))
       enddo
    else
       do k = 1, nzmsize - 1
          lsp3(:,:,k) = 0.5 * (lsp23(:,:,k+1) + lsp23(:,:,k))
       enddo
       lsp3(:,:,nzmsize) = 0.5 * (lsp23(:,:,1) + lsp23(:,:,nzmsize))
    endif

    !! Use levelset to set fluid properties (on pressure nodes)
    do k = 1, nzmsize
       do j = ph1%zst(2), ph1%zen(2)
          do i = ph1%zst(1), ph1%zen(1)
             if (lsp3(i, j, k).gt.alpha) then
                prop3(i, j, k) = dens1
             else if (lsp3(i, j, k).lt.-alpha) then
                prop3(i, j, k) = dens2
             else
                hk = one / (one + exp(-two * k99 * lsp3(i, j, k)))
                prop3(i, j, k) = dens2 + hk * (dens1 - dens2)
             endif
          enddo
       enddo
    enddo
    
    !! Interpolate fluid properties back to velocity mesh
    if ((nclz1.ne.0).and.(nclzn.ne.0)) then
       do k = 2, nzmsize
          lspv3(:,:,k) = 0.5 * (prop3(:,:,k) + prop3(:,:,k - 1))
       enddo

       k = 1
       if (nclz1.eq.1) then
          !! Zero gradient
          lspv3(:,:,k) = prop3(:,:,k)
       else
          lspv3(:,:,k) = 2.0 * prop3(:,:,k) - lspv3(:,:,k + 1)
       endif

       k = nzmsize + 1
       if (nclzn.eq.1) then
          !! Zero gradient
          lspv3(:,:,k) = prop3(:,:,k - 1)
       else
          lspv3(:,:,k) = 2.0 * prop3(:,:,k - 1) - lspv3(:,:,k - 1)
       endif
    else
       do k = 2, nzmsize
          lspv3(:,:,k) = 0.5 * (prop3(:,:,k - 1) + prop3(:,:,k))
       enddo

       k = 1
       lspv3(:,:,k) = 0.5 * (prop3(:,:,k) + prop3(:,:,nzmsize))
    endif
    call transpose_z_to_y(lspv3,lspv32,ph3)

    if ((ncly1.ne.0).and.(nclyn.ne.0)) then
       do j = 2, nymsize
          lspv2(:,j,:) = 0.5 * (lspv32(:,j - 1,:) + lspv32(:,j,:))
       enddo

       j = 1
       if (ncly1.eq.1) then
          lspv2(:,j,:) = lspv32(:,j,:)
       else
          lspv2(:,j,:) = 2.0 * lspv32(:,j,:) - lspv32(:,j + 1,:)
       endif

       j = nymsize + 1
       if (nclyn.eq.1) then
          lspv2(:,j,:) = lspv32(:,j - 1,:)
       else
          lspv2(:,j,:) = 2.0 * lspv32(:,j - 1,:) - lspv2(:,j - 1,:)
       endif
    else
       do j = 2, nymsize
          lspv2(:,j,:) = 0.5 * (lspv32(:,j,:) + lspv32(:,j - 1,:))
       enddo

       j = 1
       lspv2(:,j,:) = 0.5 * (lspv32(:,j,:) + lspv32(:,nymsize,:))
    endif
    call transpose_y_to_x(lspv2,lspv21,ph2) !nxm ny nz

    if ((nclx1.ne.0).and.(nclxn.ne.0)) then
       do i = 2, nxmsize
          rho1(i,:,:,1) = 0.5 * (lspv21(i - 1,:,:) + lspv21(i,:,:))
       enddo

       i = 1
       if (nclx1.eq.1) then
          rho1(i,:,:,1) = lspv21(i,:,:)
       else
          rho1(i,:,:,1) = 2.0 * lspv21(i,:,:) - rho1(i + 1,:,:,1)
       endif

       i = nxmsize + 1
       if (nclxn.eq.1) then
          rho1(i,:,:,1) = lspv21(i - 1,:,:)
       else
          rho1(i,:,:,1) = 2.0 * lspv21(i - 1,:,:) - rho1(i - 1,:,:,1)
       endif
    else
       do i = 2, nxmsize
          rho1(i,:,:,1) = 0.5 * (lspv21(i - 1,:,:) + lspv21(i,:,:))
       enddo

       i = 1
       rho1(i,:,:,1) = 0.5 * (lspv21(i,:,:) + lspv21(nxmsize,:,:))
    endif

    !! Use levelset to set fluid properties (on pressure nodes)
    do k = 1, nzmsize
       do j = ph1%zst(2), ph1%zen(2)
          do i = ph1%zst(1), ph1%zen(1)
             if (lsp3(i, j, k).gt.alpha) then
                prop3(i, j, k) = visc1
             else if (lsp3(i, j, k).lt.-alpha) then
                prop3(i, j, k) = visc2
             else
                hk = one / (one + exp(-two * k99 * lsp3(i, j, k)))
                prop3(i, j, k) = visc2 + hk * (visc1 - visc2)
             endif
          enddo
       enddo
    enddo
    
    !! Interpolate fluid properties back to velocity mesh
    if ((nclz1.ne.0).and.(nclzn.ne.0)) then
       do k = 2, nzmsize
          lspv3(:,:,k) = 0.5 * (prop3(:,:,k) + prop3(:,:,k - 1))
       enddo

       k = 1
       if (nclz1.eq.1) then
          !! Zero gradient
          lspv3(:,:,k) = prop3(:,:,k)
       else
          lspv3(:,:,k) = 2.0 * prop3(:,:,k) - lspv3(:,:,k + 1)
       endif

       k = nzmsize + 1
       if (nclzn.eq.1) then
          !! Zero gradient
          lspv3(:,:,k) = prop3(:,:,k - 1)
       else
          lspv3(:,:,k) = 2.0 * prop3(:,:,k - 1) - lspv3(:,:,k - 1)
       endif
    else
       do k = 2, nzmsize
          lspv3(:,:,k) = 0.5 * (prop3(:,:,k - 1) + prop3(:,:,k))
       enddo

       k = 1
       lspv3(:,:,k) = 0.5 * (prop3(:,:,k) + prop3(:,:,nzmsize))
    endif
    call transpose_z_to_y(lspv3,lspv32,ph3)

    if ((ncly1.ne.0).and.(nclyn.ne.0)) then
       do j = 2, nymsize
          lspv2(:,j,:) = 0.5 * (lspv32(:,j - 1,:) + lspv32(:,j,:))
       enddo

       j = 1
       if (ncly1.eq.1) then
          lspv2(:,j,:) = lspv32(:,j,:)
       else
          lspv2(:,j,:) = 2.0 * lspv32(:,j,:) - lspv32(:,j + 1,:)
       endif

       j = nymsize + 1
       if (nclyn.eq.1) then
          lspv2(:,j,:) = lspv32(:,j - 1,:)
       else
          lspv2(:,j,:) = 2.0 * lspv32(:,j - 1,:) - lspv2(:,j - 1,:)
       endif
    else
       do j = 2, nymsize
          lspv2(:,j,:) = 0.5 * (lspv32(:,j,:) + lspv32(:,j - 1,:))
       enddo

       j = 1
       lspv2(:,j,:) = 0.5 * (lspv32(:,j,:) + lspv32(:,nymsize,:))
    endif
    call transpose_y_to_x(lspv2,lspv21,ph2) !nxm ny nz

    if ((nclx1.ne.0).and.(nclxn.ne.0)) then
       do i = 2, nxmsize
          mu1(i,:,:) = 0.5 * (lspv21(i - 1,:,:) + lspv21(i,:,:))
       enddo

       i = 1
       if (nclx1.eq.1) then
          mu1(i,:,:) = lspv21(i,:,:)
       else
          mu1(i,:,:) = 2.0 * lspv21(i,:,:) - mu1(i + 1,:,:)
       endif

       i = nxmsize + 1
       if (nclxn.eq.1) then
          mu1(i,:,:) = lspv21(i - 1,:,:)
       else
          mu1(i,:,:) = 2.0 * lspv21(i - 1,:,:) - mu1(i - 1,:,:)
       endif
    else
       do i = 2, nxmsize
          mu1(i,:,:) = 0.5 * (lspv21(i - 1,:,:) + lspv21(i,:,:))
       enddo

       i = 1
       mu1(i,:,:) = 0.5 * (lspv21(i,:,:) + lspv21(nxmsize,:,:))
    endif

  endsubroutine update_fluid_properties

  subroutine init_wls(levelset1, S1)

    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : alloc_x

    use var, only : gradx_ls1 => ta1, grady_ls1 => tb1, gradz_ls1 => tc1, mag_grad_ls1 => td1
    use var, only : levelset2 => ta2, grady_ls2 => tb2, gradz_ls2 => tc2
    use var, only : levelset3 => ta3, gradz_ls3 => tb3

    use var, only : two
    use var, only : dx, dy, dz

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1, S1

    integer :: i, j, k
    
    call alloc_x(wx1)
    call alloc_x(wy1)
    call alloc_x(wz1)

    call transpose_x_to_y(levelset1, levelset2)
    call transpose_y_to_z(levelset2, levelset3)

    do k = 2, zsize(3)
       gradz_ls3(:, :, k) = (levelset3(:, :, k + 1) - levelset3(:, :, k - 1)) / (two * dz)
    enddo
    gradz_ls3(:, :, 1) = (levelset3(:, :, 2) - levelset3(:, :, zsize(3))) / (two * dz)
    gradz_ls3(:, :, zsize(3)) = (levelset3(:, :, 1) - levelset3(:, :, zsize(3) - 1)) / (two * dz)

    call transpose_z_to_y(gradz_ls3, gradz_ls2)

    do j = 2, ysize(2)
       grady_ls2(:, j, :) = (levelset2(:, j + 1, :) - levelset2(:, j - 1, :)) / (two * dy)
    enddo
    grady_ls2(:, 1, :) = (levelset2(:, 2, :) - levelset2(:, ysize(2), :)) / (two * dy)
    grady_ls2(:, ysize(2), :) = (levelset2(:, 1, :) - levelset2(:, ysize(2) - 1, :)) / (two * dy)

    call transpose_y_to_x(grady_ls2, grady_ls1)
    call transpose_y_to_x(gradz_ls2, gradz_ls1)

    do i = 2, xsize(1)
       gradx_ls1(i, :, :) = (levelset1(i + 1, :, :) - levelset1(i - 1, :, :)) / (two * dx)
    enddo
    gradx_ls1(1, :, :) = (levelset1(2, :, :) - levelset1(xsize(1), :, :)) / (two * dx)
    gradx_ls1(xsize(1), :, :) = (levelset1(1, :, :) - levelset1(xsize(1) - 1, :, :)) / (two * dx)

    mag_grad_ls1(:,:,:) = sqrt(gradx_ls1**2 + grady_ls1**2 + gradz_ls1**2) + epsilon(two)

    wx1(:,:,:) = S1(:,:,:) * gradx_ls1(:,:,:) / mag_grad_ls1(:,:,:)
    wy1(:,:,:) = S1(:,:,:) * grady_ls1(:,:,:) / mag_grad_ls1(:,:,:)
    wz1(:,:,:) = S1(:,:,:) * gradz_ls1(:,:,:) / mag_grad_ls1(:,:,:)
    
  endsubroutine init_wls

  subroutine surface_tension_force(sx1, sy1, sz1, rho1, levelset1)

    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : xstart, xend
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : alloc_x

    use div, only : divergence
    use navier, only : gradp

    use variables, only : derx, dery, derz
    use variables, only : sx, sy, sz

    use var, only : ep1, divu3, drho1
    
    use var, only : zero, half, one, two, three, pi
    use var, only : stfx1 => ta1, stfy1 => tb1, stfz1 => tc1, curvature1 => te1
    use var, only : di1
    use var, only : ny, nz

    use var, only : ls2 => ta2, ny2 => tb2, nz2 => tc2
    use var, only : ls3 => ta3, nz3 => tb3

    use var, only : curvature2 => td2
    use var, only : curvature3 => tc3

    use variables, only : ffz, fsz, fwz, ffy, fsy, fwy, ffx, fsx, fwx

    use var, only : di2
    use var, only : di3
    use variables, only : ffxp, fsxp, fwxp
    use variables, only : ffyp, fsyp, fwyp, ppy
    use variables, only : ffzp, fszp, fwzp

    use param, only : nclx1, ncly1, nclz1
    use param, only : dx, dy, dz, xlx, yly, zlz
    use param, only : nrhotime
    use param, only : interface_thickness
    use param, only : dens1, dens2
    use param, only : itime, ifirst

    use var, only : lsp1 => pgy1, lspv21 => pp1
    use var, only : lsp12 => uyp2, lsp2 => upi2, lspv32 => pp2, lspv2 => ppi2
    use var, only : lsp23 => uzp3, lsp3 => qq3, curvaturep3 => dv3
    use var, only : curvaturepi3 => ppi3
    use var, only : curvaturepi32 => pp2, curvaturepi2 => ppi2
    use var, only : curvaturepi21 => pp1

    use var, only : dip3, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use var, only : dip2, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use var, only : cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6

    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ph1, ph2, ph3, ph4
    use var, only : nclx1, ncly1, nclz1, nclxn, nclyn, nclzn

    use weno, only : weno5
    
    implicit none

    !! IN
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1

    !! INOUT
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: sx1, sy1, sz1

    !! LOCAL
    real(mytype) :: normmag
    integer :: nlock
    integer :: i, j, k
    real(mytype) :: alpha
    real(mytype) :: sigma ! Surface tension coefficient
    real(mytype) :: rhomean
    real(mytype) :: x, y, z
    real(mytype) :: ddelta

    alpha = interface_thickness * (dx * dy * dz)**(one / three)
    sigma = one / 121._mytype

    rhomean = (dens1 + dens2) / two !! Recommended by Dodd&Ferrante2014 for smoothing

    if (itime.eq.ifirst) then
       call alloc_x(nx1); call alloc_x(ny1); call alloc_x(nz1)
    endif

    !! Get levelset to pressure mesh
    if ((nclx1.ne.0).and.(nclxn.ne.0)) then
       do i = 1, nxmsize
          lsp1(i,:,:) = 0.5 * (levelset1(i + 1,:,:) + levelset1(i,:,:))
       enddo
    else
       do i = 1, nxmsize - 1
          lsp1(i,:,:) = 0.5 * (levelset1(i + 1,:,:)) + levelset1(i,:,:)
       enddo
       lsp1(nxmsize,:,:) = 0.5 * (levelset1(1,:,:) + levelset1(nxmsize,:,:))
    endif
    call transpose_x_to_y(lsp1,lsp12,ph4)

    if ((ncly1.ne.0).and.(nclyn.ne.0)) then
       do j = 1, nymsize
          lsp2(:,j,:) = 0.5 * (lsp12(:,j+1,:) + lsp12(:,j,:))
       enddo
    else
       do j = 1, nymsize - 1
          lsp2(:,j,:) = 0.5 * (lsp12(:,j+1,:) + lsp12(:,j,:))
       enddo
       lsp2(:,nymsize,:) = 0.5 * (lsp12(:,1,:) + lsp12(:,nymsize,:))
    endif
    call transpose_y_to_z(lsp2,lsp23,ph3)

    if ((nclz1.ne.0).and.(nclzn.ne.0)) then
       do k = 1, nzmsize
          lsp3(:,:,k) = 0.5 * (lsp23(:,:,k+1) + lsp23(:,:,k))
       enddo
    else
       do k = 1, nzmsize - 1
          lsp3(:,:,k) = 0.5 * (lsp23(:,:,k+1) + lsp23(:,:,k))
       enddo
       lsp3(:,:,nzmsize) = 0.5 * (lsp23(:,:,1) + lsp23(:,:,nzmsize))
    endif

    !! Compute the normals
    call gradp(nx1, ny1, nz1, lsp3)
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             normmag = sqrt(nx1(i,j,k)**2 + ny1(i,j,k)**2 + nz1(i,j,k)**2)
             if (normmag.gt.zero) then
                nx1(i,j,k) = nx1(i,j,k) / normmag
                ny1(i,j,k) = ny1(i,j,k) / normmag
                nz1(i,j,k) = nz1(i,j,k) / normmag
             else
                nx1(i,j,k) = zero
                ny1(i,j,k) = zero
                nz1(i,j,k) = zero
             endif
          enddo
       enddo
    enddo
    
    !! Compute curvature = div(n)
    nlock = -1
    CALL divergence(curvaturep3,rho1,nx1,ny1,nz1,ep1,drho1,divu3,nlock)

    !! Evaluate Dirac delta on pressure mesh and apply to curvature
    do k = 1, size(lsp3, 3)
       do j = 1, size(lsp3, 2)
          do i = 1, size(lsp3, 1)
             ! !! COS Dirac delta
             ! if (abs(levelset1(i, j, k)).lt.alpha) then
             !    ddelta1(i, j, k) = half * (one + cos(pi * levelset1(i, j, k) / alpha)) / alpha
             ! else
             !    ddelta1(i, j, k) = zero
             ! endif

             !! GAUSS Dirac delta
             ddelta = exp(-((lsp3(i, j, k) / alpha)**2)) / (alpha * sqrt(pi))
             
             curvaturep3(i, j, k) = ddelta * curvaturep3(i, j, k)
          enddo
       enddo
    enddo

    !! Get curvature to velocity mesh
    call interypv(curvaturepi2,curvaturepi32,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call transpose_y_to_x(curvaturepi2,curvaturepi21,ph2) !nxm ny nz
    call interxpv(curvature1,curvaturepi21,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    do k = 2, nzmsize
       do j = 1, size(curvaturep3, 2)
          do i = 1, size(curvaturep3, 1)
             if (curvaturep3(i,j,k).eq.zero) then
                curvaturepi3(i,j,k) = curvaturep3(i,j,k - 1)
             elseif (curvaturep3(i,j,k-1).eq.zero) then
                curvaturepi3(i,j,k) = curvaturep3(i,j,k)
             else
                curvaturepi3(i,j,k) = 0.5 * (curvaturep3(i,j,k) + curvaturep3(i,j,k - 1))
             endif
          enddo
       enddo
    enddo
    if ((nclz1.ne.0).and.(nclzn.ne.0)) then
       k = 1
       if (nclz1.eq.1) then
          !! Zero gradient
          curvaturepi3(:,:,k) = curvaturep3(:,:,k)
       else
          curvaturepi3(:,:,k) = 2.0 * curvaturep3(:,:,k) - curvaturepi3(:,:,k + 1)
       endif

       k = nzmsize + 1
       if (nclzn.eq.1) then
          !! Zero gradient
          curvaturepi3(:,:,k) = curvaturep3(:,:,k - 1)
       else
          curvaturepi3(:,:,k) = 2.0 * curvaturep3(:,:,k - 1) - curvaturepi3(:,:,k - 1)
       endif
    else
       k = 1
       do j = 1, size(curvaturep3, 2)
          do i = 1, size(curvaturep3, 1)
             if (curvaturep3(i,j,k).eq.zero) then
                curvaturepi3(i,j,k) = curvaturep3(i,j,nzmsize)
             elseif (curvaturep3(i,j,nzmsize).eq.zero) then
                curvaturepi3(i,j,k) = curvaturep3(i,j,k)
             else
                curvaturepi3(i,j,k) = 0.5 * (curvaturep3(i,j,k) + curvaturep3(i,j,nzmsize))
             endif
          enddo
       enddo
    endif
    call transpose_z_to_y(curvaturepi3,curvaturepi32,ph3)

    do k = 1, size(curvaturepi32, 3)
       do j = 2, nymsize
          do i = 1, size(curvaturepi32, 1)
             if (curvaturepi32(i,j,k).eq.zero) then
                curvaturepi2(i,j,k) = curvaturepi32(i,j-1,k)
             elseif (curvaturepi32(i,j-1,k).eq.zero) then
                curvaturepi2(i,j,k) = curvaturepi32(i,j,k)
             else
                curvaturepi2(:,j,:) = 0.5 * (curvaturepi32(:,j - 1,:) + curvaturepi32(:,j,:))
             endif
          enddo
       enddo
    enddo
    
    if ((ncly1.ne.0).and.(nclyn.ne.0)) then
       j = 1
       if (ncly1.eq.1) then
          curvaturepi2(:,j,:) = curvaturepi32(:,j,:)
       else
          curvaturepi2(:,j,:) = 2.0 * curvaturepi32(:,j,:) - curvaturepi32(:,j + 1,:)
       endif

       j = nymsize + 1
       if (nclyn.eq.1) then
          curvaturepi2(:,j,:) = curvaturepi32(:,j - 1,:)
       else
          curvaturepi2(:,j,:) = 2.0 * curvaturepi32(:,j - 1,:) - curvaturepi2(:,j - 1,:)
       endif
    else
       j = 1
       do k = 1, size(curvaturepi32, 3)
          do i = 1, size(curvaturepi32, 1)
             if (curvaturepi32(i,j,k).eq.zero) then
                curvaturepi2(i,j,k) = curvaturepi32(i,nymsize,k)
             elseif (curvaturepi32(i,nymsize,k).eq.zero) then
                curvaturepi2(i,j,k) = curvaturepi32(i,j,k)
             else
                curvaturepi2(:,j,:) = 0.5 * (curvaturepi32(:,nymsize,:) + curvaturepi32(:,j,:))
             endif
          enddo
       enddo
    endif
    call transpose_y_to_x(curvaturepi2,curvaturepi21,ph2) !nxm ny nz

    do k = 1, size(curvaturepi21, 3)
       do j = 1, size(curvaturepi21, 2)
          do i = 2, nxmsize
             if (curvaturepi21(i,j,k).eq.zero) then
                curvature1(i,j,k) = curvaturepi21(i-1,j,k)
             elseif (curvaturepi21(i-1,j,k).eq.zero) then
                curvature1(i,j,k) = curvaturepi21(i,j,k)
             else
                curvature1(i,:,:) = 0.5 * (curvaturepi21(i - 1,:,:) + curvaturepi21(i,:,:))
             endif
          enddo
       enddo
    enddo
    if ((nclx1.ne.0).and.(nclxn.ne.0)) then
       i = 1
       if (nclx1.eq.1) then
          curvature1(i,:,:) = curvaturepi21(i,:,:)
       else
          curvature1(i,:,:) = 2.0 * curvaturepi21(i,:,:) - curvature1(i + 1,:,:)
       endif

       i = nxmsize + 1
       if (nclxn.eq.1) then
          curvature1(i,:,:) = curvaturepi21(i - 1,:,:)
       else
          curvature1(i,:,:) = 2.0 * curvaturepi21(i - 1,:,:) - curvature1(i - 1,:,:)
       endif
    else
       i = 1
       do k = 1, size(curvaturepi21, 3)
          do j = 1, size(curvaturepi21, 2)
             do i = 2, nxmsize
                if (curvaturepi21(i,j,k).eq.zero) then
                   curvature1(i,j,k) = curvaturepi21(nxmsize,j,k)
                elseif (curvaturepi21(nxmsize,j,k).eq.zero) then
                   curvature1(i,j,k) = curvaturepi21(i,j,k)
                else
                   curvature1(i,:,:) = 0.5 * (curvaturepi21(nxmsize,:,:) + curvaturepi21(i,:,:))
                endif
             enddo
          enddo
       enddo
    endif

    !! Evaluate forcing on velocity mesh
    stfx1(:,:,:) = sigma * curvature1(:,:,:) * nx1(:,:,:)
    stfy1(:,:,:) = sigma * curvature1(:,:,:) * ny1(:,:,:)
    stfz1(:,:,:) = sigma * curvature1(:,:,:) * nz1(:,:,:)

    !! Add contribution to forcing terms
    sx1(:,:,:) = sx1(:,:,:) - stfx1(:,:,:) * rho1(:,:,:,1) / rhomean
    sy1(:,:,:) = sy1(:,:,:) - stfy1(:,:,:) * rho1(:,:,:,1) / rhomean
    sz1(:,:,:) = sz1(:,:,:) - stfz1(:,:,:) * rho1(:,:,:,1) / rhomean

  end subroutine surface_tension_force

endmodule freesurface
