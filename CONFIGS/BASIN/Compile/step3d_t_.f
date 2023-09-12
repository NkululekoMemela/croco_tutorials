      subroutine step3d_t (tile)
      implicit none
      integer*4 tile
      integer*4 trd, omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=60,   MMm0=50,   N=10)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=NPP)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
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
      parameter (ntrc_salt=0)
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
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,7,0:NPP-1)
      common/private_scratch/ A2d,A3d
      real u(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3)
      real v(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3)
      real t(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real Hz_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real z_r(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real z_w(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real Huon(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real Hvom(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real rho1(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real rho(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lm+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lm)/2
      chunk_size_E=(Mm+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mm)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lm)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mm)
      trd=omp_get_thread_num()
      call step3d_t_tile (Istr,Iend,Jstr,Jend,
     &                    A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                    A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd),
     &                    A2d(1,7,trd), A2d(1,8,trd), A2d(1,9,trd),
     &                                                A3d(1,1,trd))
      return
      end
      subroutine step3d_t_tile (Istr,Iend,Jstr,Jend,
     &                          FX,FE, WORK, FC,CF,BC,DC,EC,GC, swdk)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=60,   MMm0=50,   N=10)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=NPP)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
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
      parameter (ntrc_salt=0)
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
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real h(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real hinv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real f(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real fomn(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real xp(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real xr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real yp(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real yr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_xr/xr /grid_xp/xp /grid_yp/yp /grid_yr/yr
      real pm(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pm_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pm_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real pmon_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmon_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmon_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real grdscl(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real zob(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /Z0B_VAR/zob
      real u(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3)
      real v(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3)
      real t(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real Hz_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real z_r(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real z_w(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real Huon(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real Hvom(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real rho1(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real rho(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real visc2_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real visc2_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real visc2_sponge_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real visc2_sponge_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real diff2(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real Akv(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real Akt(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
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
      real R0,T0,S0, Tcoef, Scoef
      real weight(6,0:NWEIGHT)
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , levbfrc,   levsfrc
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
     &                      , R0,T0,S0,  Tcoef,   Scoef
     &                      , weight
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , levbfrc,   levsfrc
     &                      , got_tini
     &                      , ldefhis
      real Akv_bak
      common /scalars_akv/ Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akt_bak
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
     &        , avg_vol, avg_rho
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
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
      real sustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real svstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real bustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real bvstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real bvstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /forces_stflx/stflx
      real stflxg(0:Lm+1+padd_X,0:Mm+1+padd_E,2,NT)
      common /stfdat_stflxg/stflxg
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer*4 itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer*4 lstfgrd(NT), stf_tid(NT), stf_id(NT)
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat3/  stf_tid, stf_id
      real btflx(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /forces_btflx/btflx
      real srflx(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_srflx/srflx
      real srflxg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /srfdat_srflxg/srflxg
      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer*4 itsrf, srf_ncycle, srf_rec
      integer*4 lsrfgrd, srf_tid, srf_id
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ itsrf,srf_ncycle,srf_rec,lsrfgrd,srf_tid,srf_id
      integer*4 Istr,Iend,Jstr,Jend, itrc, i,j,k, indx, kmld
     &       ,imin,imax,jmin,jmax,iAkt,nadv
      real FX(Istr-2:Iend+2,Jstr-2:Jend+2),
     &     FE(Istr-2:Iend+2,Jstr-2:Jend+2),   cff,
     &     WORK(Istr-2:Iend+2,Jstr-2:Jend+2), epsil,
     &     FC(Istr-2:Iend+2,0:N),
     &     CF(Istr-2:Iend+2,0:N),
     &     BC(Istr-2:Iend+2,0:N),
     &     DC(Istr-2:Iend+2,0:N),
     &     EC(Istr-2:Iend+2,0:N),
     &     GC(Istr-2:Iend+2,0:N),
     &     swdk(Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
      real cff1,cff2,gama,dRz,hbltmp,sig,dXmax,dEmax,
     &     dpth,smax,amax,amaxx
      parameter (epsil=1.D-16)
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
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
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (jend.eq.Mm) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      nadv = 3
      do itrc=1,NT
        do k=1,N
          do j=Jstr,Jend
            do i=max(Istr-1,1),min(Iend+2,Lm+1)
              FX(i,j)=(t(i,j,k,nadv,itrc)-t(i-1,j,k,nadv,itrc))
            enddo
          enddo
          if (istr.eq.1) then
            do j=Jstr,Jend
              FX(0,j)=FX(1,j)
            enddo
          endif
          if (iend.eq.Lm) then
             do j=Jstr,Jend
              FX(Lm+2,j)=FX(Lm+1,j)
            enddo
          endif
          do j=Jstr,Jend
            do i=Istr-1,Iend+1
              WORK(i,j)=FX(i+1,j)-FX(i,j)
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend+1
              if (Huon(i,j,k) .gt. 0.D0) then
                cff=WORK(i-1,j)
              else
                cff=WORK(i,j)
              endif
              FX(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i-1,j,k,nadv,itrc)
     &                           -0.333333333333D0*cff )*Huon(i,j,k)
            enddo
          enddo
          do j=max(Jstr-1,1),min(Jend+2,Mm+1)
            do i=Istr,Iend
              FE(i,j)=(t(i,j,k,nadv,itrc)-t(i,j-1,k,nadv,itrc))
            enddo
          enddo
          if (jstr.eq.1) then
            do i=Istr,Iend
              FE(i,0)=FE(i,1)
            enddo
          endif
          if (jend.eq.Mm) then
            do i=Istr,Iend
              FE(i,Mm+2)=FE(i,Mm+1)
            enddo
          endif
          do j=Jstr-1,Jend+1
            do i=Istr,Iend
              WORK(i,j)=FE(i,j+1)-FE(i,j)
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              if (Hvom(i,j,k) .gt. 0.D0) then
                cff=WORK(i,j-1)
              else
                cff=WORK(i,j)
              endif
              FE(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i,j-1,k,nadv,itrc)
     &                          -0.333333333333D0*cff )*Hvom(i,j,k)
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend
              t(i,j,k,nnew,itrc)=Hz_bak(i,j,k)*t(i,j,k,nstp,itrc)
     &                    -dt*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                         +FE(i,j+1)-FE(i,j))
            enddo
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do itrc=1,NT
          do k=1,N-1
            do i=istr,iend
              FC(i,k)=t(i,j,k+1,nadv,itrc)-t(i,j,k,nadv,itrc)
            enddo
          enddo
          do i=istr,iend
            FC(i,0)=FC(i,1)
            FC(i,N)=FC(i,N-1)
          enddo
          do k=1,N
            do i=istr,iend
              cff=2.D0*FC(i,k)*FC(i,k-1)
              if (cff.gt.epsil) then
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              else
                CF(i,k)=0.D0
              endif
            enddo
          enddo
          do k=1,N-1
            do i=istr,iend
              FC(i,k)=0.5D0*(   t(i,j,k  ,nadv,itrc)
     &                +       t(i,j,k+1,nadv,itrc)
     &                -0.333333333333D0*(CF(i,k+1)-CF(i,k))
     &                                                  )*We(i,j,k)
            enddo
          enddo
          do i=istr,iend
            FC(i,0)=0.D0
            FC(i,N)=0.D0
            CF(i,0)=dt*pm(i,j)*pn(i,j)
          enddo
            do k=1,N
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &                             -CF(i,0)*(FC(i,k)-FC(i,k-1))
              enddo
            enddo
          do i=Istr,Iend
            FC(i,N)=dt*stflx(i,j,itrc)
            FC(i,0)=-dt*btflx(i,j,itrc)
          enddo
          if (itrc.eq.itemp) then
            do k=1,N-1
              do i=Istr,Iend
                FC(i,k)=0.D0
              enddo
            enddo
          endif
          if (itrc.eq.itemp) then
            do k=1,N
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+FC(i,k )
     &                                               -FC(i,k-1)
              enddo
            enddo
          endif
          indx=min(itrc,itemp)
          do i=istr,iend
            FC(i,1)=dt* Akt(i,j,1,indx)
     &                               /( z_r(i,j,2)-z_r(i,j,1) )
            cff=1.D0/(Hz(i,j,1)+FC(i,1))
            CF(i,1)= cff*FC(i,1)
            DC(i,1)= cff* t(i,j,1,nnew,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=dt* Akt(i,j,k,indx)
     &                              /( z_r(i,j,k+1)-z_r(i,j,k) )
              cff=1.D0/(Hz(i,j,k)+FC(i,k)+FC(i,k-1)*(1.D0-CF(i,k-1)))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(t(i,j,k,nnew,itrc)+FC(i,k-1)*DC(i,k-1)
     &                                                          )
            enddo
          enddo
          do i=istr,iend
             t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc)
     &                           +FC(i,N-1)*DC(i,N-1) )
     &                        /(Hz(i,j,N)+FC(i,N-1)*(1.D0-CF(i,N-1)))
        enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo
        enddo
      enddo
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
      enddo
      do itrc=1,NT
      enddo
      return
      end
