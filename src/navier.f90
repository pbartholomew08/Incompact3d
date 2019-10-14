module navier

  implicit none

  private

  public :: solve_poisson, calc_divu_constraint
  public :: pre_correc, cor_vel
  public :: momentum_to_velocity, velocity_to_momentum
  public :: gradp

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  SUBROUTINE: solve_poisson
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Takes the intermediate momentum field as input,
  !!              computes div and solves pressure-Poisson equation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE solve_poisson(pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3)

    USE decomp_2d, ONLY : mytype, xsize, zsize, ph1
    USE decomp_2d_poisson, ONLY : poisson
    USE var, ONLY : itr, itime
    USE var, ONLY : nzmsize
    USE var, ONLY : dv3
    USE var, ONLY : two
    USE param, ONLY : ntime, nrhotime, npress
    USE param, ONLY : ilmn, ivarcoeff, ifreesurface

    USE div, ONLY : divergence

    IMPLICIT NONE

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime), INTENT(IN) :: drho1
    REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3

    !! Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    !! Locals
    INTEGER :: nlock, poissiter, piter
    LOGICAL :: converged
    REAL(mytype) :: atol, rtol

    nlock = 1 !! Corresponds to computing div(u*)
    converged = .FALSE.
    poissiter = 0

    atol = 1.0e-12_mytype !! Absolute tolerance for Poisson solver
    rtol = 1.0e-14_mytype !! Relative tolerance for Poisson solver

    IF (ilmn.AND.ivarcoeff) THEN
       !! Variable-coefficient Poisson solver works on div(u), not div(rho u)
       !! rho u -> u
       CALL momentum_to_velocity(rho1, ux1, uy1, uz1)
    ENDIF

    CALL divergence(pp3(:,:,:,1),rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)
    IF ((ilmn.AND.ivarcoeff).OR.ifreesurface) THEN
       dv3(:,:,:) = pp3(:,:,:,1)

       IF ((itime.GT.1).AND.(itr.EQ.1)) THEN
          !! Extrapolate pressure
          pp3(:,:,:,1) = two * pp3(:,:,:,2) - pp3(:,:,:,3)
          CALL gradp(px1,py1,pz1,pp3(:,:,:,1))

          !! Store previous pressure field
          pp3(:,:,:,3) = pp3(:,:,:,2)
       ENDIF
    ENDIF

    DO WHILE(.NOT.converged)
       IF (ivarcoeff.OR.ifreesurface) THEN

          !! Test convergence
          CALL test_varcoeff(converged, pp3, dv3, atol, rtol, poissiter)

          IF (.NOT.converged) THEN
             !! Store current pressure field
             pp3(:,:,:,2) = pp3(:,:,:,1)
             
             !! Evaluate additional RHS terms
             CALL calc_varcoeff_rhs(pp3(:,:,:,1), rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, &
                  poissiter)
          ENDIF
       ENDIF

       IF (.NOT.converged) THEN
          CALL poisson(pp3(:,:,:,1))

          !! Need to update pressure gradient here for varcoeff
          CALL gradp(px1,py1,pz1,pp3(:,:,:,1))

          IF (((.NOT.ilmn).AND.(.NOT.ivarcoeff)).AND.(.NOT.ifreesurface)) THEN
             !! Once-through solver
             !! - Incompressible flow
             !! - LMN - constant-coefficient solver
             converged = .TRUE.
          ENDIF
       ENDIF

       poissiter = poissiter + 1
    ENDDO

    IF (ilmn.AND.ivarcoeff) THEN
       !! Variable-coefficient Poisson solver works on div(u), not div(rho u)
       !! u -> rho u
       CALL velocity_to_momentum(rho1, ux1, uy1, uz1)
    ENDIF

  END SUBROUTINE solve_poisson

  !********************************************************************
  !subroutine COR_VEL
  !Correction of u* by the pressure gradient to get a divergence free
  !field
  ! input : px,py,pz
  ! output : ux,uy,uz
  !written by SL 2018
  !********************************************************************
  subroutine cor_vel (rho,ux,uy,uz,px,py,pz,pp3)

    USE MPI
    
    USE decomp_2d
    USE variables
    USE param

    USE var, only : nzmsize

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(inout) :: px,py,pz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),intent(in) :: rho
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    integer :: ierr
    real(mytype) :: rhomin, rho0

    if (.not.ifreesurface) then
       ux(:,:,:)=ux(:,:,:)-px(:,:,:)
       uy(:,:,:)=uy(:,:,:)-py(:,:,:)
       uz(:,:,:)=uz(:,:,:)-pz(:,:,:)
    else
       !! Compute rho0
       rhomin = minval(rho)

       call MPI_ALLREDUCE(rhomin,rho0,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)

       !! First correct according to new pressure gradient
       ux(:,:,:) = ux(:,:,:) - px(:,:,:) / rho0
       uy(:,:,:) = uy(:,:,:) - py(:,:,:) / rho0
       uz(:,:,:) = uz(:,:,:) - pz(:,:,:) / rho0

       !! Get previous iter pressure gradient
       CALL gradp(px,py,pz,pp3(:,:,:,2))

       !! Not correct with the previous pressure gradient
       ux(:,:,:)=ux(:,:,:)-(one / rho(:,:,:,1) - one / rho0) * px(:,:,:)
       uy(:,:,:)=uy(:,:,:)-(one / rho(:,:,:,1) - one / rho0) * py(:,:,:)
       uz(:,:,:)=uz(:,:,:)-(one / rho(:,:,:,1) - one / rho0) * pz(:,:,:)

       !! Reset pressure gradient
       CALL gradp(px,py,pz,pp3(:,:,:,1))
    endif

    return
  end subroutine cor_vel


  !********************************************************************
  !subroutine GRADP
  !Computation of the pressure gradient from the pressure mesh to the
  !velocity mesh
  !Saving pressure gradients on boundaries for correct imposition of
  !BCs on u* via the fractional step methodi (it is not possible to
  !impose BC after correction by pressure gradient otherwise lost of
  !incompressibility--> BCs are imposed on u*
  !
  ! input: pp3 - pressure field (on pressure mesh)
  ! output: px1, py1, pz1 - pressure gradients (on velocity mesh)
  !written by SL 2018
  !********************************************************************
  subroutine gradp(px1,py1,pz1,pp3)

    USE param
    USE decomp_2d
    USE variables
    USE MPI
    USE var, only: pp1,pgy1,pgz1,di1,pp2,ppi2,pgy2,pgz2,pgzi2,dip2,&
         pgz3,ppi3,dip3,nxmsize,nymsize,nzmsize

    USE forces, only : ppi1

    implicit none

    integer :: i,j,k

    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1

    !WORK Z-PENCILS
    call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    call derzpv(pgz3,pp3,dip3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

    !WORK Y-PENCILS
    call transpose_z_to_y(pgz3,pgz2,ph3) !nxm nym nz
    call transpose_z_to_y(ppi3,pp2,ph3)

    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call derypv(pgy2,pp2,dip2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call interypv(pgzi2,pgz2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

    !WORK X-PENCILS

    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call transpose_y_to_x(pgy2,pgy1,ph2)
    call transpose_y_to_x(pgzi2,pgz1,ph2)

    call derxpv(px1,pp1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(py1,pgy1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(pz1,pgz1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)

    if (iforces) then
       call interxpv(ppi1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)
    endif

    !we are in X pencils:
    if (nclx1.eq.2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyx1(j,k)=py1(1,j,k)/gdt(itr)
             dpdzx1(j,k)=pz1(1,j,k)/gdt(itr)
          enddo
       enddo
    endif
    if (nclxn.eq.2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyxn(j,k)=py1(nx,j,k)/gdt(itr)
             dpdzxn(j,k)=pz1(nx,j,k)/gdt(itr)
          enddo
       enddo
    endif

    if (ncly1.eq.2) then
       if (xsize(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxy1(i,k)=px1(i,1,k)/gdt(itr)
                dpdzy1(i,k)=pz1(i,1,k)/gdt(itr)
             enddo
          enddo
       endif
    endif
    if (nclyn.eq.2) then
       if (xsize(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=px1(i,ny,k)/gdt(itr)
                dpdzyn(i,k)=pz1(i,ny,k)/gdt(itr)
             enddo
          enddo
       endif
    endif

    if (nclz1.eq.2) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxz1(i,j)=py1(i,j,1)/gdt(itr)
                dpdyz1(i,j)=pz1(i,j,1)/gdt(itr)
             enddo
          enddo
       endif
    endif
    if (nclzn.eq.2) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxzn(i,j)=py1(i,j,xsize(3))/gdt(itr)
                dpdyzn(i,j)=pz1(i,j,xsize(3))/gdt(itr)
             enddo
          enddo
       endif
    endif

    return
  end subroutine gradp
  !*******************************************************************
  subroutine pre_correc(ux,uy,uz,ep)

    USE decomp_2d
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    integer :: i,j,k,is
    real(mytype) :: ut,ut1,utt,ut11

    integer :: code
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, code)

    !********NCLX==2*************************************
    !we are in X pencils:
    if ((itype.eq.itype_channel).and.(nclx1==2.and.nclxn==2)) then

       !Computatation of the flow rate Inflow/Outflow
       ut1=zero
       do k=1,xsize(3)
          do j=1,xsize(2)
             ut1=ut1+bxx1(j,k)
          enddo
       enddo
       call MPI_ALLREDUCE(ut1,ut11,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       ut11=ut11/(real(ny*nz,mytype))
       ut=zero
       do k=1,xsize(3)
          do j=1,xsize(2)
             ut=ut+bxxn(j,k)
          enddo
       enddo
       call MPI_ALLREDUCE(ut,utt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       utt=utt/(real(ny*nz,mytype))
       if (nrank==0) print *,'Flow rate x I/O/O-I',real(ut11,4),real(utt,4),real(utt-ut11,4)
       do k=1,xsize(3)
          do j=1,xsize(2)
             bxxn(j,k)=bxxn(j,k)-utt+ut11
          enddo
       enddo

    endif

    if (nclx1==2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
             dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
          enddo
       enddo
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=bxx1(j,k)
             uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
             uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
          enddo
       enddo
    endif
    if (nclxn==2) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
             dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
          enddo
       enddo
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=bxxn(j,k)
             uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
             uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
          enddo
       enddo
    endif


    !********NCLX==1*************************************
    if (nclx1==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(1 ,j,k)=zero
          enddo
       enddo
    endif
    if (nclxn==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             ux(nx,j,k)=zero
          enddo
       enddo
    endif

    !********NCLY==2*************************************
    if (ncly1==2) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
                dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
             enddo
          enddo
          do k=1,xsize(3)
             do i=1,xsize(1)
                ux(i,1,k)=byx1(i,k)+dpdxy1(i,k)
                uy(i,1,k)=byy1(i,k)
                uz(i,1,k)=byz1(i,k)+dpdzy1(i,k)
             enddo
          enddo
       endif
    endif

    if (nclyn==2) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
                dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
             enddo
          enddo
       endif
       if (dims(1)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                ux(i,xsize(2),k)=byxn(i,k)+dpdxyn(i,k)
                uy(i,xsize(2),k)=byyn(i,k)
                uz(i,xsize(2),k)=byzn(i,k)+dpdzyn(i,k)
             enddo
          enddo
       elseif (ny - (nym / dims(1)) == xstart(2)) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                ux(i,xsize(2),k)=byxn(i,k)+dpdxyn(i,k)
                uy(i,xsize(2),k)=byyn(i,k)
                uz(i,xsize(2),k)=byzn(i,k)+dpdzyn(i,k)
             enddo
          enddo
       endif
    endif

    !********NCLY==1*************************************
    if (ncly1==1) then
       if (xstart(2)==1) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,1,k)=zero
             enddo
          enddo
       endif
    endif

    if (nclyn==1) then
       if (xend(2)==ny) then
          do k=1,xsize(3)
             do i=1,xsize(1)
                uy(i,xsize(2),k)=zero
             enddo
          enddo
       endif
    endif


    !********NCLZ==2*************************************
    if (nclz1==2) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
                dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
             enddo
          enddo
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux(i,j,1)=bzx1(i,j)+dpdxz1(i,j)
                uy(i,j,1)=bzy1(i,j)+dpdyz1(i,j)
                uz(i,j,1)=bzz1(i,j)
             enddo
          enddo
       endif
    endif

    if (nclzn==2) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
                dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
             enddo
          enddo
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux(i,j,xsize(3))=bzxn(i,j)+dpdxzn(i,j)
                uy(i,j,xsize(3))=bzyn(i,j)+dpdyzn(i,j)
                uz(i,j,xsize(3))=bzzn(i,j)
             enddo
          enddo
       endif
    endif
    !********NCLZ==1************************************* !just to reforce free-slip condition
    if (nclz1==1) then
       if (xstart(3)==1) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,1)=zero
             enddo
          enddo
       endif
    endif

    if (nclzn==1) then
       if (xend(3)==nz) then
          do j=1,xsize(2)
             do i=1,xsize(1)
                uz(i,j,xsize(3))=zero
             enddo
          enddo
       endif
    endif

    if (iibm==1) then !solid body old school
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
       call body(ux1,uy1,uz1,ep1,1)
       call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
    endif

    return
  end subroutine pre_correc
  !*******************************************************************

  !! Convert to/from conserved/primary variables
  SUBROUTINE primary_to_conserved(rho1, var1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = rho1(:,:,:,1) * var1(:,:,:)

  ENDSUBROUTINE primary_to_conserved
  SUBROUTINE velocity_to_momentum (rho1, ux1, uy1, uz1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime
    USE var, ONLY : ilmn

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    IF (.NOT.ilmn) THEN
       RETURN
    ENDIF

    CALL primary_to_conserved(rho1, ux1)
    CALL primary_to_conserved(rho1, uy1)
    CALL primary_to_conserved(rho1, uz1)

  ENDSUBROUTINE velocity_to_momentum
  SUBROUTINE conserved_to_primary(rho1, var1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: var1

    var1(:,:,:) = var1(:,:,:) / rho1(:,:,:,1)

  ENDSUBROUTINE conserved_to_primary
  SUBROUTINE momentum_to_velocity (rho1, ux1, uy1, uz1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : nrhotime
    USE var, ONLY : ilmn

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    IF (.NOT.ilmn) THEN
       RETURN
    ENDIF

    CALL conserved_to_primary(rho1, ux1)
    CALL conserved_to_primary(rho1, uy1)
    CALL conserved_to_primary(rho1, uz1)

  ENDSUBROUTINE momentum_to_velocity

  !! Calculate velocity-divergence constraint
  SUBROUTINE calc_divu_constraint(divu3, rho1, phi1)

    USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
    USE decomp_2d, ONLY : transpose_x_to_y, transpose_y_to_z
    USE param, ONLY : nrhotime, zero, ilmn, pressure0, imultispecies, massfrac, mol_weight
    USE param, ONLY : ibirman_eos
    USE param, ONLY : xnu, prandtl
    USE param, ONLY : one
    USE variables

    USE var, ONLY : ta1, tb1, tc1, td1, di1
    USE var, ONLY : phi2, ta2, tb2, tc2, td2, te2, di2
    USE var, ONLY : phi3, ta3, tb3, tc3, td3, rho3, di3

    IMPLICIT NONE

    INTEGER :: is

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    REAL(mytype), INTENT(OUT), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3

    IF (ilmn.and.(.not.ibirman_eos)) THEN
       !!------------------------------------------------------------------------------
       !! X-pencil

       !! We need temperature
       CALL calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

       CALL derxx (tb1, ta1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)
       IF (imultispecies) THEN
          tb1(:,:,:) = (xnu / prandtl) * tb1(:,:,:) / ta1(:,:,:)

          !! Calc mean molecular weight
          td1(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                td1(:,:,:) = td1(:,:,:) + phi1(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          td1(:,:,:) = one / td1(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL derxx (tc1, phi1(:,:,:,is), di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)
                tb1(:,:,:) = tb1(:,:,:) + (xnu / sc(is)) * (td1(:,:,:) / mol_weight(is)) * tc1(:,:,:)
             ENDIF
          ENDDO
       ENDIF

       CALL transpose_x_to_y(ta1, ta2)        !! Temperature
       CALL transpose_x_to_y(tb1, tb2)        !! d2Tdx2
       IF (imultispecies) THEN
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))
             ENDIF
          ENDDO
       ENDIF

       !!------------------------------------------------------------------------------
       !! Y-pencil
       CALL deryy (tc2, ta2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
       IF (imultispecies) THEN
          tc2(:,:,:) = (xnu / prandtl) * tc2(:,:,:) / ta2(:,:,:)

          !! Calc mean molecular weight
          te2(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                te2(:,:,:) = te2(:,:,:) + phi2(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          te2(:,:,:) = one / te2(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL deryy (td2, phi2(:,:,:,is), di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
                tc2(:,:,:) = tc2(:,:,:) + (xnu / sc(is)) * (te2(:,:,:) / mol_weight(is)) * td2(:,:,:)
             ENDIF
          ENDDO
       ENDIF
       tb2(:,:,:) = tb2(:,:,:) + tc2(:,:,:)

       CALL transpose_y_to_z(ta2, ta3)        !! Temperature
       CALL transpose_y_to_z(tb2, tb3)        !! d2Tdx2 + d2Tdy2
       IF (imultispecies) THEN
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL transpose_y_to_z(phi2(:,:,:,is), phi3(:,:,:,is))
             ENDIF
          ENDDO
       ENDIF

       !!------------------------------------------------------------------------------
       !! Z-pencil
       CALL derzz (divu3, ta3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
       IF (imultispecies) THEN
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / ta3(:,:,:)

          !! Calc mean molecular weight
          td3(:,:,:) = zero
          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                td3(:,:,:) = td3(:,:,:) + phi3(:,:,:,is) / mol_weight(is)
             ENDIF
          ENDDO
          td3(:,:,:) = one / td3(:,:,:)

          DO is = 1, numscalar
             IF (massfrac(is)) THEN
                CALL derzz (tc3, phi3(:,:,:,is), di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
                divu3(:,:,:) = divu3(:,:,:) + (xnu / sc(is)) * (td3(:,:,:) / mol_weight(is)) * tc3(:,:,:)
             ENDIF
          ENDDO
       ENDIF
       divu3(:,:,:) = divu3(:,:,:) + tb3(:,:,:)

       IF (imultispecies) THEN
          !! Thus far we have computed rho * divu, want divu
          CALL calc_rho_eos(rho3, ta3, phi3, tb3, zsize(1), zsize(2), zsize(3))
          divu3(:,:,:) = divu3(:,:,:) / rho3(:,:,:)
       ELSE
          divu3(:,:,:) = (xnu / prandtl) * divu3(:,:,:) / pressure0
       ENDIF
    ELSE
       divu3(:,:,:) = zero
    ENDIF

  ENDSUBROUTINE calc_divu_constraint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: test_varcoeff
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Tests convergence of the variable-coefficient Poisson solver
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE test_varcoeff(converged, pp3, dv3, atol, rtol, poissiter)

    USE MPI
    USE decomp_2d, ONLY: mytype, ph1, real_type, nrank
    USE var, ONLY : nzmsize
    USE param, ONLY : npress
    USE variables, ONLY : nxm, nym, nzm

    IMPLICIT NONE

    !! INPUTS
    REAL(mytype), INTENT(INOUT), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype), INTENT(IN), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3
    REAL(mytype), INTENT(IN) :: atol, rtol
    INTEGER, INTENT(IN) :: poissiter

    !! OUTPUTS
    LOGICAL, INTENT(OUT) :: converged

    !! LOCALS
    INTEGER :: ierr
    REAL(mytype) :: errloc, errglob
    REAL(mytype), SAVE :: divup3norm

    IF (poissiter.EQ.0) THEN
       errloc = SUM(dv3**2)
       CALL MPI_ALLREDUCE(errloc,divup3norm,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       divup3norm = SQRT(divup3norm / nxm / nym / nzm)

       IF (nrank.EQ.0) THEN
          PRINT *, "Solving variable-coefficient Poisson equation:"
          PRINT *, "+ RMS div(u*) - div(u): ", divup3norm
       ENDIF
    ELSE
       !! Compute RMS change
       errloc = SUM((pp3(:,:,:,1) - pp3(:,:,:,2))**2)
       CALL MPI_ALLREDUCE(errloc,errglob,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       errglob = SQRT(errglob / nxm / nym / nzm)

       IF (nrank.EQ.0) THEN
          PRINT *, "+ RMS change in pressure: ", errglob
       ENDIF

       IF (errglob.LE.atol) THEN
          converged = .TRUE.
          IF (nrank.EQ.0) THEN
             PRINT *, "- Converged: atol"
          ENDIF
       ENDIF

       !! Compare RMS change to size of |div(u*) - div(u)|
       IF (errglob.LT.(rtol * divup3norm)) THEN
          converged = .TRUE.
          IF (nrank.EQ.0) THEN
             PRINT *, "- Converged: rtol"
          ENDIF
       ENDIF
    ENDIF

  ENDSUBROUTINE test_varcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: calc_varcoeff_rhs
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Computes RHS of the variable-coefficient Poisson solver
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calc_varcoeff_rhs(pp3, rho1, px1, py1, pz1, dv3, drho1, ep1, divu3, poissiter)

    USE MPI

    USE decomp_2d

    USE param, ONLY : nrhotime, ntime
    USE param, ONLY : one

    USE var, ONLY : td1, te1, tf1
    USE var, ONLY : nzmsize

    USE div, ONLY : divergence

    IMPLICIT NONE

    !! INPUTS
    INTEGER, INTENT(IN) :: poissiter
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ep1
    REAL(mytype), INTENT(IN), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3
    REAL(mytype), INTENT(IN), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: dv3

    !! OUTPUTS
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: pp3

    !! LOCALS
    INTEGER :: nlock, ierr
    REAL(mytype) :: rhomin
    REAL(mytype), SAVE :: rho0

    IF (poissiter.EQ.0) THEN
       !! Compute rho0
       rhomin = MINVAL(rho1)

       CALL MPI_ALLREDUCE(rhomin,rho0,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    ENDIF

    td1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * px1(:,:,:)
    te1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * py1(:,:,:)
    tf1(:,:,:) = (one - rho0 / rho1(:,:,:,1)) * pz1(:,:,:)

    nlock = -1 !! Don't do any funny business with LMN
    CALL divergence(pp3,rho1,td1,te1,tf1,ep1,drho1,divu3,nlock)

    !! lapl(p) = div((1 - rho0/rho) grad(p)) + rho0(div(u*) - div(u))
    !! dv3 contains div(u*) - div(u)
    pp3(:,:,:) = pp3(:,:,:) + rho0 * dv3(:,:,:)

  ENDSUBROUTINE calc_varcoeff_rhs

endmodule navier
