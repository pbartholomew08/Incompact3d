module div

  implicit none
  
  private
  public :: divergence

contains
  
  !********************************************************************
  !subroutine DIVERGENCe
  !Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
  ! input : ux1,uy1,uz1,ep1 (on velocity mesh)
  ! output : pp3 (on pressure mesh)
  !written by SL 2018
  !********************************************************************
  subroutine divergence (pp3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,nlock)

    USE param
    USE decomp_2d
    USE variables
    USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
         duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
    USE MPI

    implicit none

    !  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

    !X PENCILS NX NY NZ  -->NXM NY NZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime),intent(in) :: drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),intent(in) :: rho1
    !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)),intent(in) :: divu3
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

    integer :: nvect3,i,j,k,nlock
    integer :: code
    real(mytype) :: tmax,tmoy,tmax1,tmoy1

    nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

    if (iibm.eq.0) then
       ta1(:,:,:) = ux1(:,:,:)
       tb1(:,:,:) = uy1(:,:,:)
       tc1(:,:,:) = uz1(:,:,:)
    else
       ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)
       tb1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)
       tc1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)
    endif

    !WORK X-PENCILS

    call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

    if (ilmn.and.(nlock.gt.0)) then
       if ((nlock.eq.1).and.(.not.ivarcoeff)) then
          !! Approximate -div(rho u) using ddt(rho)
          call extrapol_drhodt(ta1, rho1, drho1)
       elseif ((nlock.eq.2).or.ivarcoeff) then
          !! Need to check our error against divu constraint
          !! Or else we are solving the variable-coefficient Poisson equation
          call transpose_z_to_y(-divu3, ta2)
          call transpose_y_to_x(ta2, ta1)
       endif
       call interxvp(pgy1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
       pp1(:,:,:) = pp1(:,:,:) + pgy1(:,:,:)
    endif

    call interxvp(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    call interxvp(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

    call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
    call transpose_x_to_y(pgy1,uyp2,ph4)
    call transpose_x_to_y(pgz1,uzp2,ph4)

    !WORK Y-PENCILS
    call interyvp(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
    call deryvp(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

    !! Compute sum dudx + dvdy
    duydypi2(:,:,:) = duydypi2(:,:,:) + upi2(:,:,:)

    call interyvp(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

    call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
    call transpose_y_to_z(upi2,uzp3,ph3)

    !WORK Z-PENCILS
    call interzvp(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
    call derzvp(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

    !! Compute sum dudx + dvdy + dwdz
    pp3(:,:,:) = pp3(:,:,:) + po3(:,:,:)

    if (nlock==2) then
       pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
    endif

    tmax=-1609._mytype
    tmoy=0._mytype
    do k=1,nzmsize
       do j=ph1%zst(2),ph1%zen(2)
          do i=ph1%zst(1),ph1%zen(1)
             if (abs(pp3(i,j,k)).gt.tmax) tmax=abs(pp3(i,j,k))
             tmoy=tmoy+abs(pp3(i,j,k))
          enddo
       enddo
    enddo
    tmoy=tmoy/nvect3

    call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    if ((nrank==0).and.(nlock.gt.0)) then
       if (nlock==2) then
          print *,'DIV U  max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
       else
          print *,'DIV U* max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
       endif
    endif

    return
  end subroutine divergence

  SUBROUTINE extrapol_drhodt(drhodt1_next, rho1, drho1)

    USE decomp_2d, ONLY : mytype, xsize, nrank
    USE param, ONLY : ntime, nrhotime, itime, itimescheme, itr, dt, gdt, irestart
    USE param, ONLY : half, three, four
    USE param, ONLY : ibirman_eos

    IMPLICIT NONE

    INTEGER :: subitr

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(OUT), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    IF (itimescheme.EQ.1) THEN
       !! EULER
       drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
    ELSEIF (itimescheme.EQ.2) THEN
       !! AB2
       IF ((itime.EQ.1).AND.(irestart.EQ.0)) THEN
          drhodt1_next(:,:,:) = (rho1(:,:,:,1) - rho1(:,:,:,2)) / dt
       ELSE
          drhodt1_next(:,:,:) = three * rho1(:,:,:,1) - four * rho1(:,:,:,2) + rho1(:,:,:,3)
          drhodt1_next(:,:,:) = half * drhodt1_next(:,:,:) / dt
       ENDIF
       ! ELSEIF (itimescheme.EQ.3) THEN
       !    !! AB3
       ! ELSEIF (itimescheme.EQ.4) THEN
       !    !! AB4
    ELSEIF (itimescheme.EQ.5) THEN
       !! RK3
       IF (itime.GT.1) THEN
          drhodt1_next(:,:,:) = rho1(:,:,:,2)
          DO subitr = 1, itr
             drhodt1_next(:,:,:) = drhodt1_next(:,:,:) + (gdt(subitr) / dt) &
                  * (rho1(:,:,:,2) - rho1(:,:,:,3))
          ENDDO
       ELSE
          drhodt1_next(:,:,:) = drho1(:,:,:,1)
       ENDIF
    ELSE
       IF (nrank.EQ.0) THEN
          PRINT *, "Extrapolating drhodt not implemented for timescheme:", itimescheme
          STOP
       ENDIF
    ENDIF

    IF (ibirman_eos) THEN
       CALL birman_drhodt_corr(drhodt1_next, rho1)
    ENDIF

  ENDSUBROUTINE extrapol_drhodt

  SUBROUTINE birman_drhodt_corr(drhodt1_next, rho1)

    USE decomp_2d, ONLY : mytype, xsize, ysize, zsize
    USE decomp_2d, ONLY : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    USE variables, ONLY : derxx, deryy, derzz
    USE param, ONLY : nrhotime
    USE param, ONLY : xnu, prandtl

    USE var, ONLY : td1, te1, di1, sx, sfxp, ssxp, swxp
    USE var, ONLY : rho2, ta2, tb2, di2, sy, sfyp, ssyp, swyp
    USE var, ONLY : rho3, ta3, di3, sz, sfzp, sszp, swzp

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1_next

    REAL(mytype) :: invpe

    invpe = xnu / prandtl

    CALL transpose_x_to_y(rho1(:,:,:,1), rho2)
    CALL transpose_y_to_z(rho2, rho3)

    !! Diffusion term
    CALL derzz (ta3,rho3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    CALL transpose_z_to_y(ta3, tb2)

    CALL deryy (ta2,rho2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)
    CALL transpose_y_to_x(ta2, te1)

    CALL derxx (td1,rho1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
    td1(:,:,:) = td1(:,:,:) + te1(:,:,:)

    drhodt1_next(:,:,:) = drhodt1_next(:,:,:) - invpe * td1(:,:,:)

  ENDSUBROUTINE birman_drhodt_corr
  
endmodule div
