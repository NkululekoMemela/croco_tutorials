      subroutine lmd_skpp_tile (Istr,Iend,Jstr,Jend, Kv,Kt,Ks,
     &                                       my_hbl, Bo,Bosol,
     &                     Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS,
     &                          Bfsfc_bl,Cr,FC,wrk1,wrk2,wrk3,
     &                                       swr_frac, my_kbl)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=40,   MMm0=80,   N=20)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=NPP)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Msrc
      parameter (Msrc=2)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=LLm, Mm=MMm)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=1)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real h(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real hinv(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real f(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real fomn(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real xp(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real xr(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real yp(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real yr(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /grid_xr/xr /grid_xp/xp /grid_yp/yp /grid_yr/yr
      real pm(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pn(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real om_r(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real on_r(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real om_u(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real on_u(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real om_v(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real on_v(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real om_p(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real on_p(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pn_u(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pm_v(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pm_u(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pn_v(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real pmon_p(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pmon_r(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pmon_u(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pnom_p(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pnom_r(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pnom_v(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real grdscl(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pmask(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real umask(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real vmask(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real pmask2(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real zob(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /Z0B_VAR/zob
      real u(0:Lm+1+padd_X,-1:Mm+2+padd_E,N,3)
      real v(0:Lm+1+padd_X,-1:Mm+2+padd_E,N,3)
      real t(0:Lm+1+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      real Hz_bak(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      real z_r(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      real z_w(0:Lm+1+padd_X,-1:Mm+2+padd_E,0:N)
      real Huon(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      real Hvom(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(0:Lm+1+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real rho1(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      real rho(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(0:Lm+1+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real sustr(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real svstr(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real bustr(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real bvstr(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(0:Lm+1+padd_X,-1:Mm+2+padd_E,2)
      real bvstrg(0:Lm+1+padd_X,-1:Mm+2+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(0:Lm+1+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_stflx/stflx
      real btflx(0:Lm+1+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_btflx/btflx
      real srflx(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /forces_srflx/srflx
      real Akv(0:Lm+1+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt(0:Lm+1+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(0:Lm+1+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real ustar(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_ustar/ustar
      integer*4 kbl(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      integer*4 kbbl(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real hbbl(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_kbl/ kbl
      common /lmd_kpp_hbbl/ hbbl
      common /lmd_kpp_kbbl/ kbbl
      real hbls(0:Lm+1+padd_X,-1:Mm+2+padd_E,2)
      common /lmd_kpp_hbl/ hbls
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
      logical ldefhis
      logical got_tini(NT)
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , weight
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , got_tini
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, khbl
      integer*4 imin,imax,jmin,jmax
      integer*4 my_kbl(Istr-2:Iend+2     )
      real    Kv    (Istr-2:Iend+2,Jstr-2:Jend+2,0:N ),
     &        Kt    (Istr-2:Iend+2,Jstr-2:Jend+2,0:N ),
     &        Ks    (Istr-2:Iend+2,Jstr-2:Jend+2,0:N ),
     &      swr_frac(Istr-2:Iend+2,Jstr-2:Jend+2,0:N ),
     &        my_hbl(Istr-2:Iend+2,Jstr-2:Jend+2     ),
     &        Bo    (Istr-2:Iend+2,Jstr-2:Jend+2     ),
     &        Bosol (Istr-2:Iend+2,Jstr-2:Jend+2     ),
     &           Gm1(Istr-2:Iend+2     ),
     &        dGm1dS(Istr-2:Iend+2     ),
     &           Gt1(Istr-2:Iend+2     ),
     &        dGt1dS(Istr-2:Iend+2     ),
     &           Gs1(Istr-2:Iend+2     ),
     &        dGs1dS(Istr-2:Iend+2     ),
     &      Bfsfc_bl(Istr-2:Iend+2     ),
     &        Cr    (Istr-2:Iend+2,0:N ),
     &        FC    (Istr-2:Iend+2,0:N ),
     &          wrk1(Istr-2:Iend+2,0:N ),
     &          wrk2(Istr-2:Iend+2,0:N ),
     &          wrk3(Istr-2:Iend+2,Jstr-2:Jend+2     )
      real eps
      parameter (eps=1.D-20)
      real ustar3,Bfsfc,zscale,zetahat, ws,wm,z_bl
      real Av_bl,dAv_bl, f1,At_bl,dAt_bl,As_bl,dAs_bl
      real Kern, Vtc, Vtsq,  sigma, cff,cff1, cff_up,cff_dn
      real zeta_m, a_m, c_m, zeta_s, a_s,c_s,r2,r3,r4
      real a1,a2,a3
      real nu0c, Cv, Ric, Ri_inv,
     &     betaT,epssfc,Cstar,Cg,C_Ek
      real ustarb
      real su_r, sv_r, ustokes
      parameter (
     &   Cv=1.8D0,
     &   Ric=0.45D0,
     &   Ri_inv=1.D0/Ric,
     &   betaT=-0.2D0,
     &   epssfc=0.1D0,
     &   Cstar=10.D0,
     &   nu0c=0.1D0,
     &   C_Ek=258.D0,
     &   zeta_m=-0.2D0,
     &   a_m=1.257D0,
     &   c_m=8.360D0,
     &   zeta_s=-1.0D0,
     &   a_s=-28.86D0,
     &   c_s=98.96D0,
     &   r2=0.5D0, r3=1.D0/3.D0, r4=0.25D0
     &                                                           )
      if (istr.eq.1) then
        imin=Istr
      else
        imin=Istr-1
      endif
      if (iend.eq.Lm) then
        imax=Iend
      else
        imax=Iend+1
      endif
      Cg=Cstar*vonKar*(c_s*vonKar*epssfc)**r3
      Vtc= Cv * sqrt(-betaT/(c_s*epssfc)) / (Ric*vonKar**2)
      call alfabeta_tile (Istr,Iend,Jstr,Jend, Bosol,Bo)
      do j=Jstr-1,Jend+1
        do i=imin,imax
          Bo(i,j)= 0
     &                 +g*( Bosol(i,j)*(stflx(i,j,itemp)-srflx(i,j)))
     &                 -g*Bo(i,j)*stflx(i,j,isalt)
          Bosol(i,j)=g*Bosol(i,j)*srflx(i,j)
          su_r=0.5D0*(sustr(i,j)+sustr(i+1,j))
          sv_r=0.5D0*(svstr(i,j)+svstr(i,j+1))
          ustar(i,j)=sqrt(sqrt(su_r**2+sv_r**2))
        enddo
      enddo
      do k=0,N-1
        do j=Jstr-1,Jend+1
          do i=imin,imax
            wrk3(i,j)=z_w(i,j,k)-z_w(i,j,N)
          enddo
        enddo
        call lmd_swfrac_ER_tile (Istr,Iend,Jstr,Jend,1.D0,wrk3,my_hbl)
        do j=Jstr-1,Jend+1
          do i=imin,imax
            swr_frac(i,j,k)=my_hbl(i,j)
          enddo
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=imin,imax
          swr_frac(i,j,N)=1.D0
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=imin,imax
          my_hbl(i,j)  = hbls(i,j,nstp)
          my_kbl(i)    = 0
          FC(i,N)      = 0.D0
          Cr(i,N)      = 0.D0
          Cr(i,0)      = 0.D0
        enddo
        do k=1,N-1
         do i=imin,imax
          cff=1.D0/(Hz(i,j,k)+Hz(i,j,k+1))
          wrk1(i,k)=cff*( u(i,j,k+1,nstp)+u(i+1,j,k+1,nstp)
     &                 -u(i,j,k  ,nstp)-u(i+1,j,k  ,nstp))
          wrk2(i,k)=cff*( v(i,j,k+1,nstp)+v(i,j+1,k+1,nstp)
     &                 -v(i,j,k  ,nstp)-v(i,j+1,k  ,nstp))
         enddo
        enddo
        do i=imin,imax
          wrk1(i,N)=wrk1(i,N-1)
          wrk2(i,N)=wrk2(i,N-1)
          wrk1(i,0)=wrk1(i,  1)
          wrk2(i,0)=wrk2(i,  1)
        enddo
        do k=N,1,-1
          do i=imin,imax
            zscale = z_w(i,j,N)-z_w(i,j,k-1)
            Kern   = zscale/(zscale+epssfc*my_hbl(i,j))
            Bfsfc  = Bo(i,j) +Bosol(i,j)*(1.D0-swr_frac(i,j,k-1))
          if (Bfsfc .lt. 0.D0) zscale=min(zscale,
     &                                  my_hbl(i,j)*epssfc)
          zscale=zscale*rmask(i,j)
          zetahat=vonKar*zscale*Bfsfc
          ustar3=ustar(i,j)**3
          if (zetahat .ge. 0.D0) then
            ws=vonKar*ustar(i,j)*ustar3/max(ustar3+5.D0*zetahat,
     &                                                 1.D-20)
          elseif (zetahat .gt. zeta_s*ustar3) then
            ws=vonKar*( (ustar3-16.D0*zetahat)/ustar(i,j) )**r2
          else
            ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
          endif
            cff=bvf(i,j,k)*bvf(i,j,k-1)
            if (cff.gt.0.D0) then
              cff=cff/(bvf(i,j,k)+bvf(i,j,k-1))
            else
              cff=0.D0
            endif
            FC(i,k-1)=FC(i,k) + Kern*Hz(i,j,k)*
     &       ( 0.375D0*( wrk1(i,k)**2 + wrk1(i,k-1)**2
     &               + wrk2(i,k)**2 + wrk2(i,k-1)**2 )
     &       + 0.25D0 *( wrk1(i,k-1)*wrk1(i,k)
     &               + wrk2(i,k-1)*wrk2(i,k) )
     &       - Ri_inv*( cff +
     &              0.25D0*(bvf(i,j,k)+bvf(i,j,k-1)) )
     &       - C_Ek*f(i,j)*f(i,j) )
            Vtsq=Vtc*ws*sqrt(max(0.D0, bvf(i,j,k-1)))
            Cr(i,k-1)=FC(i,k-1) +Vtsq
            if (my_kbl(i).eq.0 .and. Cr(i,k-1).lt.0.D0)
     &                          my_kbl(i)=k
          enddo
        enddo
        do i=imin,imax
          if (my_kbl(i).gt.0) then
            k=my_kbl(i)
            my_hbl(i,j)=z_w(i,j,N)-( z_w(i,j,k-1)*Cr(i,k)
     &                              -z_w(i,j,k)*Cr(i,k-1)
     &                              )/(Cr(i,k)-Cr(i,k-1))
          else
            my_hbl(i,j)=z_w(i,j,N)-z_w(i,j,0)
          endif
        enddo
      enddo
      if (istr.eq.1) then
        do j=Jstr-1,Jend+1
          my_hbl(Istr-1,j)=my_hbl(Istr,j)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=Jstr-1,Jend+1
          my_hbl(Iend+1,j)=my_hbl(Iend,j)
        enddo
      endif
      do j=Jstr,Jend+1
        do i=Istr,Iend+1
          wrk3(i,j)=0.25D0*(my_hbl(i,j)  +my_hbl(i-1,j)
     &                  +my_hbl(i,j-1)+my_hbl(i-1,j-1))
     &                   *pmask2(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          cff=0.25D0*(pmask2(i,j)   +pmask2(i+1,j)
     &             +pmask2(i,j+1) +pmask2(i+1,j+1))
          my_hbl(i,j)=(1.D0-cff)*my_hbl(i,j)+
     &              0.25D0*(wrk3(i,j)  +wrk3(i+1,j)
     &                   +wrk3(i,j+1)+wrk3(i+1,j+1))
          my_hbl(i,j)=my_hbl(i,j)*rmask(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=istr,iend
          kbl(i,j)=N
        enddo
        do k=N-1,1,-1
          do i=istr,iend
            my_hbl(i,j)=min(my_hbl(i,j),z_w(i,j,N)-z_w(i,j,0))
     &                                             *rmask(i,j)
            if (z_w(i,j,k) .gt. z_w(i,j,N)-my_hbl(i,j)) kbl(i,j)=k
          enddo
        enddo
        do i=istr,iend
          k=kbl(i,j)
          z_bl=z_w(i,j,N)-my_hbl(i,j)
          zscale=my_hbl(i,j)
          if (swr_frac(i,j,k-1).gt. 0.D0) then
            Bfsfc=Bo(i,j) +Bosol(i,j)*( 1.D0 -swr_frac(i,j,k-1)
     &              *swr_frac(i,j,k)*(z_w(i,j,k)-z_w(i,j,k-1))
     &               /( swr_frac(i,j,k  )*(z_w(i,j,k)   -z_bl)
     &                 +swr_frac(i,j,k-1)*(z_bl -z_w(i,j,k-1))
     &                                                      ))
          else
            Bfsfc=Bo(i,j)+Bosol(i,j)
          endif
            if (Bfsfc.lt.0.D0) zscale=min(zscale,
     &                       my_hbl(i,j)*epssfc)
            zscale=zscale*rmask(i,j)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
            if (zetahat.ge.0.D0) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5.D0*zetahat,
     &                                                    1.D-20)
              ws=wm
            else
              if (zetahat .gt. zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16.D0*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat .gt. zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16.D0*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
          f1=5.0D0 * max(0.D0, Bfsfc) * vonKar/(ustar(i,j)**4+eps)
          cff=1.D0/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)
         if(k.eq.1) then
            ustarb=SQRT(SQRT((0.5D0*(bustr(i,j)+bustr(i+1,j)))**2+
     &                       (0.5D0*(bvstr(i,j)+bvstr(i,j+1)))**2))
            dAv_bl=vonKar*ustarb
            Av_bl=dAv_bl*(z_bl-z_w(i,j,0))
            dAt_bl=vonKar*ustarb
            At_bl=dAt_bl*(z_bl-z_w(i,j,0))
            dAs_bl=vonKar*ustarb
            As_bl=dAs_bl*(z_bl-z_w(i,j,0))
          else
            Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
            dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
            At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
            dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
            As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
            dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          endif
          Gm1(i)=Av_bl/(my_hbl(i,j)*wm+eps)
          dGm1dS(i)=min(0.D0, Av_bl*f1-dAv_bl/(wm+eps))
          Gt1(i)=At_bl/(my_hbl(i,j)*ws+eps)
          dGt1dS(i)=min(0.D0, At_bl*f1-dAt_bl/(ws+eps))
          Gs1(i)=As_bl/(my_hbl(i,j)*ws+eps)
          dGs1dS(i)=min(0.D0, As_bl*f1-dAs_bl/(ws+eps))
          Bfsfc_bl(i)=Bfsfc
        enddo
        do i=istr,iend
          khbl=kbl(i,j)
          do k=N-1,khbl,-1
            Bfsfc=Bfsfc_bl(i)
            zscale=z_w(i,j,N)-z_w(i,j,k)
            if (Bfsfc.lt.0.D0) zscale=min(zscale,
     &                       my_hbl(i,j)*epssfc)
            zscale=zscale*rmask(i,j)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
            if (zetahat.ge.0.D0) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5.D0*zetahat,
     &                                                    1.D-20)
              ws=wm
            else
              if (zetahat .gt. zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16.D0*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat .gt. zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16.D0*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
            sigma=(z_w(i,j,N)-z_w(i,j,k))/max(my_hbl(i,j),eps)
            a1=sigma-2.D0
            a2=3.D0-2.D0*sigma
            a3=sigma-1.D0
            if (sigma.lt.0.07D0) then
              cff=0.5D0*(sigma-0.07D0)**2/0.07D0
            else
              cff=0.D0
            endif
            Kv(i,j,k)=wm*my_hbl(i,j)*( cff + sigma*( 1.D0+sigma*(
     &               a1+a2*Gm1(i)+a3*dGm1dS(i) )))
            Kt(i,j,k)=ws*my_hbl(i,j)*( cff + sigma*( 1.D0+sigma*(
     &               a1+a2*Gt1(i)+a3*dGt1dS(i) )))
            Ks(i,j,k)=ws*my_hbl(i,j)*( cff + sigma*( 1.D0+sigma*(
     &               a1+a2*Gs1(i)+a3*dGs1dS(i) )))
          enddo
          do k=khbl-1,1,-1
          enddo
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          hbls(i,j,3-nstp)=my_hbl(i,j)
        enddo
      enddo
      if (istr.eq.1) then
        do j=jstr,jend
          hbls(istr-1,j,3-nstp)=hbls(istr,j,3-nstp)
          kbl (istr-1,j)       =kbl(istr,j)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=jstr,jend
          hbls(iend+1,j,3-nstp)=hbls(iend,j,3-nstp)
          kbl(iend+1,j)=kbl(iend,j)
        enddo
      endif
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                     hbls(0,-1,3-nstp))
      return
      end
