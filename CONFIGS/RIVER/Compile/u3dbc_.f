      subroutine u3dbc_tile (Istr,Iend,Jstr,Jend,grad)
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
      real rhoA(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real rhoS(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real rvfrc(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real rufrc_bak(0:Lm+1+padd_X,-1:Mm+2+padd_E,2)
      real rvfrc_bak(0:Lm+1+padd_X,-1:Mm+2+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real DU_avg1(0:Lm+1+padd_X,-1:Mm+2+padd_E,5)
      real DV_avg1(0:Lm+1+padd_X,-1:Mm+2+padd_E,5)
      real DU_avg2(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      real DV_avg2(0:Lm+1+padd_X,-1:Mm+2+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
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
      real bry_time(2)
      common /bry_indices_array/ bry_time
      real bry_cycle
      common /bry_indices_real/ bry_cycle
      integer*4 bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      common /bry_indices_integer/ bry_id, bry_time_id, bry_ncycle,
     &                             bry_rec, itbry, ntbry
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
      real    grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real    cff,eps,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.D-20)
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      if (istr.eq.1) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (iend.eq.Lm) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (jstr.eq.1) then
        JstrR=Jstr-2
      else
        JstrR=Jstr
      endif
      if (jend.eq.Mm) then
        JendR=Jend+2
      else
        JendR=Jend
      endif
      if (istr.eq.1) then
        do k=1,N
          do j=jstr,jend
            u(istr,j,k,nnew)=0.D0
          enddo
        enddo
      endif
      if (iend.eq.Lm) then
        do k=1,N
          do j=jstr,jend
            u(iend+1,j,k,nnew)=0.D0
          enddo
        enddo
      endif
      return
      end
