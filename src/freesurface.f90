module freesurface

  implicit none

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
  subroutine reinit_ls(levelset1, neps)

    use MPI

    use decomp_2d, only : mytype, real_type, xsize, nrank
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use param, only : zero, one, two, three, four, six, ten
    use param, only : dt, dx, dy, dz

    !! Derivatives
    use weno, only : weno5

    implicit none

    !! In
    integer, intent(in) :: neps

    !! InOut
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(inout) :: levelset1

    !! Local
    integer :: i, j, k

    integer :: iter
    logical :: converged
    integer :: ierr

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: S1, &
         mag_grad_ls1, gradx_ls1, grady_ls1, gradz_ls1, &
         levelset1_old, &
         rhs_old, rhs_1, rhs_2

    real(mytype) :: dtau
    real(mytype) :: delta_ls, global_delta_ls
    real(mytype) :: deltax
    integer :: ctr, global_ctr
    real(mytype) :: alpha

    real(mytype) :: eps
    integer :: niter

    real(mytype), parameter :: cfl = one / ten

    deltax = (dx * dy * dz)**(one / three)
    eps = deltax !! SUSSMAN1994
    alpha = two * deltax

    dtau = cfl * deltax
    niter = int((float(neps) * deltax) / dtau)

    !! Step to steady state
    converged = .false.
    iter = 0
    dtau = dt
    do while (.not.converged)
       if (nrank.eq.0) then
          print *, "Level-set reinitialisation: ", iter
       endif

       if (iter == 0) then
          S1(:,:,:) = levelset1(:,:,:)
       else
          call compute_reinit_smooth(S1, levelset1, mag_grad_ls1, eps)
       endif
       call reinit_ls_grad(gradx_ls1, grady_ls1, gradz_ls1, levelset1, S1)
       mag_grad_ls1 = sqrt(gradx_ls1**2 + grady_ls1**2 + gradz_ls1**2)

       !! Update level-set
       levelset1_old(:,:,:) = levelset1(:,:,:)

       ! Solve an ODE
       rhs_old(:,:,:) = S1(:,:,:) * (one - mag_grad_ls1(:,:,:))
       levelset1(:,:,:) = levelset1_old(:,:,:) + dtau * rhs_old

       call compute_reinit_smooth(S1, levelset1, mag_grad_ls1, eps)
       call reinit_ls_grad(gradx_ls1, grady_ls1, gradz_ls1, levelset1, S1)
       mag_grad_ls1 = sqrt(gradx_ls1**2 + grady_ls1**2 + gradz_ls1**2)
       rhs_1(:,:,:) = S1(:,:,:) * (one - mag_grad_ls1(:,:,:))
       levelset1(:,:,:) = levelset1_old(:,:,:) + (dtau / four) * (rhs_old(:,:,:) + rhs_1(:,:,:))

       call compute_reinit_smooth(S1, levelset1, mag_grad_ls1, eps)
       call reinit_ls_grad(gradx_ls1, grady_ls1, gradz_ls1, levelset1, S1)
       mag_grad_ls1 = sqrt(gradx_ls1**2 + grady_ls1**2 + gradz_ls1**2)
       rhs_2(:,:,:) = S1(:,:,:) * (one - mag_grad_ls1(:,:,:))
       levelset1(:,:,:) = levelset1_old(:,:,:) &
            + (dtau / six) * (rhs_old(:,:,:) + four * rhs_2(:,:,:) + rhs_1(:,:,:))

       !! Test convergence (SUSSMAN1994)
       delta_ls = zero
       ctr = 0
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                delta_ls = delta_ls + abs(levelset1(i, j, k) - levelset1_old(i, j, k))
                ctr = ctr + 1
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(delta_ls,global_delta_ls,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       global_ctr = ctr !call MPI_ALLREDUCE(ctr,global_ctr,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

       if ((global_ctr==0).or.(global_delta_ls < 0.0005_mytype)) then
          converged = .true.
       endif

       iter = iter + 1

       if (iter==niter) then
          converged = .true.
       endif

       if (converged.and.(nrank.eq.0)) then
          print *, "Level-set delta: ", global_delta_ls
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

    use variables, only : nz

    implicit none

    !! In
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: S1

    !! Out
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: gradx_ls1, grady_ls1, gradz_ls1

    !! Locals
    integer :: i, j, k
    integer, parameter :: fixedbc = 2

    call transpose_x_to_y(levelset1, levelset2)
    call transpose_y_to_z(levelset2, levelset3)

    ! Z
    if (nz.gt.1) then
       wz3(:,:,:) = one
       call weno5 (gradz_ls3, levelset3, wz3, 3, fixedbc, fixedbc, zsize(1), zsize(2), zsize(3), 1)
       call transpose_z_to_y(gradz_ls3, gradz_ls2)
       call transpose_y_to_x(gradz_ls2, gradz_ls1m)

       wz3(:,:,:) = -one
       call weno5 (gradz_ls3, levelset3, wz3, 3, fixedbc, fixedbc, zsize(1), zsize(2), zsize(3), 1)
       call transpose_z_to_y(gradz_ls3, gradz_ls2)
       call transpose_y_to_x(gradz_ls2, gradz_ls1p)
    else
       gradz_ls1m(:,:,:) = zero
       gradz_ls1p(:,:,:) = zero
    endif

    ! Y
    wy2(:,:,:) = one
    call weno5 (grady_ls2, levelset2, wy2, 2, fixedbc, fixedbc, ysize(1), ysize(2), ysize(3), 1)
    call transpose_y_to_x(grady_ls2, grady_ls1m)

    wy2(:,:,:) = -one
    call weno5 (grady_ls2, levelset2, wy2, 2, fixedbc, fixedbc, ysize(1), ysize(2), ysize(3), 1)
    call transpose_y_to_x(grady_ls2, grady_ls1p)

    ! X
    wx1 = one
    call weno5 (gradx_ls1m, levelset1, wx1, 1, fixedbc, fixedbc, xsize(1), xsize(2), xsize(3), 1)
    wx1 = -one
    call weno5 (gradx_ls1p, levelset1, wx1, 1, fixedbc, fixedbc, xsize(1), xsize(2), xsize(3), 1)
    
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

  subroutine compute_reinit_smooth (S1, levelset1, mag_grad_ls1, eps)

    use decomp_2d, only : mytype, xsize

    use var, only : one

    implicit none

    !! In
    real(mytype), intent(in) :: eps
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: levelset1, mag_grad_ls1

    !! Out
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: S1

    !! Local
    real(mytype), parameter :: macheps = epsilon(one)

    ! !! SUSSMAN1994
    ! S1(:,:,:) = levelset1(:,:,:) / sqrt(levelset1(:,:,:)**2 + eps**2)

    !! MCSHERRY2017
    S1(:,:,:) = levelset1(:,:,:) &
         / (sqrt(levelset1(:,:,:)**2 + (mag_grad_ls1(:,:,:) * eps)**2) + macheps)

  endsubroutine compute_reinit_smooth

  subroutine update_fluid_properties(rho1, mu1, phi1)

    use decomp_2d, only : mytype, xsize
    use var, only : one, two, three, pi, five
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
    real(mytype) :: alpha

    alpha = five * (dx * dy * dz)**(one / three)

    !!---------------------------------
    if (ilevelset.gt.0) then
       !! Use levelset to set density field
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                if (phi1(i, j, k, ilevelset).gt.alpha) then
                   !! Fluid 1
                   rho1(i, j, k, 1) = dens1
                   mu1(i, j, k) = visc1
                else if (phi1(i, j, k, ilevelset).lt.-alpha) then
                   !! Fluid 2
                   rho1(i, j, k, 1) = dens2
                   mu1(i, j, k) = visc2
                else
                   !! Interface: smooth properties
                   rho1(i, j, k, 1) = ((dens1 + dens2) &
                        + (dens1 - dens2) * sin(pi * phi1(i, j, k, ilevelset) / (two * alpha))) &
                        / (two * dens1)
                   mu1(i, j, k) = ((visc1 + visc2) &
                        + (visc1 - visc2) * sin(pi * phi1(i, j, k, ilevelset) / (two * alpha))) &
                        / (two * visc1)
                endif
             enddo
          enddo
       enddo
    endif

  endsubroutine update_fluid_properties

endmodule freesurface
