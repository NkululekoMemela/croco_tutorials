      subroutine step3d_uv2 (tile)
      implicit none
      integer*4 tile, trd, omp_get_thread_num
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
      call step3d_uv2_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd),
     &                                  A2d(1,2,trd), A2d(1,3,trd),
     &                                                A2d(1,4,trd)
     &                                                            )
      return
      end
      subroutine step3d_uv2_tile (Istr,Iend,Jstr,Jend, BC,CF,FC,DC
     &                                                            )
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
      real zeta(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      real ubar(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      real vbar(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
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
      real rhoA(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rhoS(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rvfrc(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rufrc_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real rvfrc_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real DU_avg1(0:Lm+1+padd_X,0:Mm+1+padd_E,5)
      real DV_avg1(0:Lm+1+padd_X,0:Mm+1+padd_E,5)
      real DU_avg2(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real DV_avg2(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
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
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
      real BC(Istr-2:Iend+2,0:N),
     &     CF(Istr-2:Iend+2,0:N),
     &     FC(Istr-2:Iend+2,0:N), cff,
     &     DC(Istr-2:Iend+2,0:N)
      real grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real dpth,aa,cc
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
      do j=Jstr,Jend
        do i=IstrU,Iend
          FC(i,0)=0.D0
          FC(i,N)=0.D0
        enddo
        do k=1,N-1
          do i=IstrU,Iend
            FC(i,k)=-dt*(Akv(i,j,k)+Akv(i-1,j,k))
     &                 /( z_r(i,j,k+1)+z_r(i-1,j,k+1)
     &                   -z_r(i,j,k  )-z_r(i-1,j,k  ))
          enddo
        enddo
        do k=1,N
          do i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))-FC(i,k)-FC(i,k-1)
          enddo
        enddo
        do i=IstrU,Iend
          cff=1.D0/BC(i,1)
          CF(i,1)=cff*( FC(i,1)
     &                                 )
          DC(i,1)=cff*DC(i,1)
        enddo
        do k=2,N-1
          do i=IstrU,Iend
              aa = FC(i,k-1)
              cc = FC(i,k  )
            cff=1.D0/(BC(i,k)-aa*CF(i,k-1))
            CF(i,k)=cff*cc
            DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
          enddo
        enddo
        do i=IstrU,Iend
              aa = FC(i,N-1)
          DC(i,N)=(DC(i,N)-aa*DC(i,N-1))
     &           /(BC(i,N)-aa*CF(i,N-1))
          CF(i,0)=0.5D0*(Hz(i,j,N)+Hz(i-1,j,N))
          DC(i,0)=CF(i,0)*DC(i,N)
        enddo
        do k=N-1,1,-1
          do i=IstrU,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            cff=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))
            CF(i,0)=CF(i,0)+cff
            DC(i,0)=DC(i,0)+cff*DC(i,k)
          enddo
        enddo
        do i=IstrU,Iend
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j,nnew))
     &                             /(CF(i,0)*on_u(i,j))
        enddo
        do k=1,N
          do i=IstrU,Iend
            u(i,j,k,nnew)=(DC(i,k)-DC(i,0))
          enddo
        enddo
        do k=1,N
          do i=IstrU,Iend
          enddo
        enddo
        if (j.ge.JstrV) then
          do i=Istr,Iend
            FC(i,0)=0.D0
            FC(i,N)=0.D0
          enddo
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=-dt*(Akv(i,j,k)+Akv(i,j-1,k))
     &                   /( z_r(i,j,k+1)+z_r(i,j-1,k+1)
     &                     -z_r(i,j,k  )-z_r(i,j-1,k  ))
            enddo
          enddo
          do k=1,N
            do i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))-FC(i,k)-FC(i,k-1)
            enddo
          enddo
          do i=Istr,Iend
            cff=1.D0/BC(i,1)
            CF(i,1)=cff*( FC(i,1)
     &                                 )
            DC(i,1)=cff*DC(i,1)
          enddo
          do k=2,N-1
            do i=Istr,Iend
              aa = FC(i,k-1)
              cc = FC(i,k  )
              cff=1.D0/(BC(i,k)-aa*CF(i,k-1))
              CF(i,k)=cff*cc
              DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
            enddo
          enddo
          do i=Istr,Iend
              aa = FC(i,N-1)
            DC(i,N)=(DC(i,N)-aa*DC(i,N-1))
     &             /(BC(i,N)-aa*CF(i,N-1))
            CF(i,0)=0.5D0*(Hz(i,j,N)+Hz(i,j-1,N))
            DC(i,0)=CF(i,0)*DC(i,N)
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              cff=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))
              CF(i,0)=CF(i,0)+cff
              DC(i,0)=DC(i,0)+cff*DC(i,k)
            enddo
          enddo
          do i=Istr,Iend
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j,nnew))
     &                               /(CF(i,0)*om_v(i,j))
          enddo
          do k=1,N
            do i=Istr,Iend
              v(i,j,k,nnew)=(DC(i,k)-DC(i,0))
          enddo
        enddo
          do k=1,N
            do i=Istr,Iend
          enddo
        enddo
        endif
      enddo
      call v3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call u3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      do j=JstrR,JendR
        do i=Istr,IendR
          DC(i,0)=0.D0
          CF(i,0)=0.D0
          FC(i,0)=0.D0
        enddo
        do k=1,N,+1
          do i=Istr,IendR
            DC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo
        do i=Istr,IendR
          DC(i,0)=1.D0/DC(i,0)
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j,nnew))
          ubar(i,j,knew)=DC(i,0)*DU_avg1(i,j,nnew)
        enddo
        do k=N,1,-1
          do i=Istr,IendR
            u(i,j,k,nnew)=(u(i,j,k,nnew)-CF(i,0))
            FC(i,k)=0.75D0*Huon(i,j,k) +
     &              0.125D0*DC(i,k)*(u(i,j,k,nstp)+u(i,j,k,nnew))
            FC(i,0)=FC(i,0)+FC(i,k)
          enddo
        enddo
        do i=Istr,IendR
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        enddo
        do k=1,N,+1
          do i=Istr,IendR
            Huon(i,j,k)=FC(i,k)-DC(i,k)*FC(i,0)
          enddo
        enddo
        if (j.ge.Jstr) then
          do i=IstrR,IendR
            DC(i,0)=0.D0
            CF(i,0)=0.D0
            FC(i,0)=0.D0
          enddo
          do k=1,N,+1
            do i=IstrR,IendR
              DC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
          do i=IstrR,IendR
            DC(i,0)=1.D0/DC(i,0)
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j,nnew))
            vbar(i,j,knew)=DC(i,0)*DV_avg1(i,j,nnew)
          enddo
          do k=N,1,-1
            do i=IstrR,IendR
              v(i,j,k,nnew)=(v(i,j,k,nnew)-CF(i,0))
              FC(i,k)=0.75D0*Hvom(i,j,k)+
     &                0.125D0*DC(i,k)*(v(i,j,k,nstp)+v(i,j,k,nnew))
              FC(i,0)=FC(i,0)+FC(i,k)
            enddo
          enddo
          do i=IstrR,IendR
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          enddo
          do k=1,N,+1
            do i=IstrR,IendR
              Hvom(i,j,k)=FC(i,k)-DC(i,k)*FC(i,0)
            enddo
          enddo
        endif
      enddo
      return
      end
