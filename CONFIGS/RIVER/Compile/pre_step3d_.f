      subroutine pre_step3d (tile)
      implicit none
      integer*4 tile,  trd,omp_get_thread_num
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
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,7,0:NPP-1)
      integer*4 B2d(N2d,0:NPP-1)
      common/private_scratch/ A2d,A3d
      common/private_scratch_bis/ B2d
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
      call pre_step3d_tile (Istr,Iend,Jstr,Jend,
     &                A3d(1,1,trd), A3d(1,2,trd), A3d(1,3,trd),
     &                A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd)
     &                                          , A3d(1,4,trd)
     &                                                        )
      return
      end
      subroutine pre_step3d_tile (Istr,Iend,Jstr,Jend, ru,rv,rw,
     &                                                 FC,CF,DC,
     &                                               FX,FE,WORK
     &                                                 ,Hz_half
     &                                                         )
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
      integer*4 Istr,Iend,Jstr,Jend, itrc, i,j,k, indx
     &       ,imin,imax,jmin,jmax,nadv,iAkt
     &       ,is
      real   ru(Istr-2:Iend+2,Jstr-2:Jend+2,N),    cff,
     &       rv(Istr-2:Iend+2,Jstr-2:Jend+2,N),    cff1,
     &       rw(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &       FC(Istr-2:Iend+2,0:N),  cff2,
     &       CF(Istr-2:Iend+2,0:N),
     &       DC(Istr-2:Iend+2,0:N),  gamma,
     &       FX(Istr-2:Iend+2,Jstr-2:Jend+2),      epsil,
     &       FE(Istr-2:Iend+2,Jstr-2:Jend+2),        cdt,
     &     WORK(Istr-2:Iend+2,Jstr-2:Jend+2)
      real Hz_half(Istr-2:Iend+2,Jstr-2:Jend+2,N)
      parameter (gamma=1.D0/6.D0, epsil=1.D-16)
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
      real zeta(0:Lm+1+padd_X,-1:Mm+2+padd_E,4)
      real ubar(0:Lm+1+padd_X,-1:Mm+2+padd_E,4)
      real vbar(0:Lm+1+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
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
      real Qbar0(Msrc)
      common /sources_Qbar0/ Qbar0
      real Qbar(Msrc)
      common /sources_Qbar/ Qbar
      real Qsrc(Msrc,N)
      common /source_Qsrc/ Qsrc
      real Qshape(Msrc,N)
      common /source_Qshape/ Qshape
      real Tsrc(Msrc,N,NT)
      common /source_Tsrc/ Tsrc
      real Tsrc0(Msrc,NT)
      common /source_Tsrc0/ Tsrc0
      real lasrc(Msrc)
      common /source_lasrc/ lasrc
      real losrc(Msrc)
      common /source_losrc/ losrc
      integer*4 Nsrc
      common /source_Nsrc/ Nsrc
      integer*4 Dsrc(Msrc)
      common /source_Dsrc/ Dsrc
      integer*4 Isrc(Msrc)
      common /source_Isrc/ Isrc
      integer*4 Jsrc(Msrc)
      common /source_Jsrc/ Jsrc
      logical Lsrc(Msrc,NT)
      common /source_Lsrc/ Lsrc
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
      indx=3-nstp
      nadv=  nstp
      if (iic.eq.ntstart) then
        cff=0.5D0*dt
        cff1=1.D0
        cff2=0.D0
      else
        cff=(1.D0-gamma)*dt
        cff1=0.5D0+gamma
        cff2=0.5D0-gamma
      endif
      do k=1,N
        do j=Jstr-1,Jend
          do i=IstrU-1,Iend
            Hz_half(i,j,k)=cff1*Hz(i,j,k)+cff2*Hz_bak(i,j,k)
     &        -cff*pm(i,j)*pn(i,j)*( Huon(i+1,j,k)-Huon(i,j,k)
     &                              +Hvom(i,j+1,k)-Hvom(i,j,k)
     &                                  +We(i,j,k)-We(i,j,k-1)
     &                                                       )
          enddo
        enddo
      enddo
      if (istr.eq.1) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (iend.eq.Lm) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
      jmin=Jstr-2
      jmax=Jend+2
      do itrc=1,NT
        do k=1,N
          do j=Jstr,Jend
            do i=max(Istr-1,1),min(Iend+2,Lm+1)
              FX(i,j)=(t(i,j,k,nadv,itrc)-t(i-1,j,k,nadv,itrc))
     &                                               *umask(i,j)
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
              WORK(i,j)=0.5D0*(FX(i+1,j)+FX(i,j))
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend+1
              FX(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i-1,j,k,nadv,itrc)
     &                     -0.333333333333D0*(WORK(i,j)-WORK(i-1,j))
     &                                                )*Huon(i,j,k)
            enddo
          enddo
          do j=Jstr-1,Jend+2
            do i=Istr,Iend
              FE(i,j)=(t(i,j,k,nadv,itrc)-t(i,j-1,k,nadv,itrc))
     &                                               *vmask(i,j)
            enddo
          enddo
          do j=Jstr-1,Jend+1
            do i=Istr,Iend
              WORK(i,j)=0.5D0*(FE(i,j+1)+FE(i,j))
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              FE(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i,j-1,k,nadv,itrc)
     &                     -0.333333333333D0*(WORK(i,j)-WORK(i,j-1))
     &                                               )*Hvom(i,j,k)
            enddo
          enddo
          do is=1,Nsrc
           i=Isrc(is)
           j=Jsrc(is)
            if (Istr.le.i .and. i.le.Iend+1
     &                   .and. Jstr.le.j .and. j.le.Jend+1) then
              if (Dsrc(is).eq.0) then
                if (Lsrc(is,itrc)) then
                  FX(i,j)=Huon(i,j,k)*Tsrc(is,k,itrc)
                else
                  if (rmask(i,j).eq.0 .and. rmask(i-1,j).eq.1) then
                    FX(i,j)=Huon(i,j,k)*t(i-1,j,k,nstp,itrc)
                  elseif (rmask(i,j).eq.1.D0 .and. rmask(i-1,j).eq.0) 
     &                                                              then
                    FX(i,j)=Huon(i,j,k)*t(i  ,j,k,nstp,itrc)
                  endif
                endif
              elseif(Dsrc(is).eq.1) then
                if (Lsrc(is,itrc)) then
                  FE(i,j)=Hvom(i,j,k)*Tsrc(is,k,itrc)
                else
                  if (rmask(i,j).eq.0 .and. rmask(i,j-1).eq.1) then
                    FE(i,j)=Hvom(i,j,k)*t(i,j-1,k,nstp,itrc)
                  elseif (rmask(i,j).eq.1 .and. rmask(i,j-1).eq.0) then
                    FE(i,j)=Hvom(i,j,k)*t(i,j  ,k,nstp,itrc)
                  endif
                endif
              endif
            endif
          enddo
          if (iic.eq.ntstart) then
            cff=0.5D0*dt
            do j=Jstr,Jend
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=Hz(i,j,k)*t(i,j,k,nstp,itrc)
     &                      -cff*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                            +FE(i,j+1)-FE(i,j))
              enddo
            enddo
          else
            cff=(1.D0-gamma)*dt
            cff1=0.5D0+gamma
            cff2=0.5D0-gamma
            do j=Jstr,Jend
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=cff1*Hz(i,j,k)*t(i,j,k,nstp,itrc)
     &                            +cff2*Hz_bak(i,j,k)*t(i,j,k,indx,itrc)
     &                         -cff*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                               +FE(i,j+1)-FE(i,j))
              enddo
            enddo
          endif
        enddo
      enddo
      if (iic.eq.ntstart) then
            cdt=0.5D0*dt
      else
            cdt=(1.D0-gamma)*dt
      endif
      do j=Jstr,Jend
        do i=Istr,Iend
           do k=1,N
             DC(i,k)=1.D0/Hz_half(i,j,k)
           enddo
           DC(i,0)=cdt*pn(i,j)*pm(i,j)
        enddo
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
              t(i,j,k,nnew,itrc)=DC(i,k)*( t(i,j,k,nnew,itrc)
     &               -DC(i,0)*(FC(i,k)-FC(i,k-1)))
            enddo
          enddo
        enddo
        do i=IstrU,Iend
          DC(i,0)=pm_u(i,j)*pn_u(i,j)
        enddo
        if (iic.eq.ntstart) then
          do k=1,N
            do i=IstrU,Iend
              cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i-1,j,k))
              u(i,j,k,nnew)=(u(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))
     &                                           +cdt*DC(i,0)*ru(i,j,k)
     &                      )*cff
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+
     &                                         Hz(i-1,j,k))
            enddo
          enddo
        else
          cff1=0.5D0+gamma
          cff2=0.5D0-gamma
          do k=1,N
            do i=IstrU,Iend
              cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i-1,j,k))
              u(i,j,k,nnew)=( cff1*u(i,j,k,nstp)*0.5D0*(Hz(i  ,j,k)+
     &                                                Hz(i-1,j,k))
     &                       +cff2*u(i,j,k,indx)*0.5D0*(Hz_bak(i  ,j,k)+
     &                                                Hz_bak(i-1,j,k))
     &                       +cdt*DC(i,0)*ru(i,j,k)
     &                                                     )*cff
              u(i,j,k,indx)=u(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+
     &                                         Hz(i-1,j,k))
            enddo
          enddo
         endif
        if (j.ge.Jstr) then
          do i=Istr,Iend
            DC(i,0)=pm_v(i,j)*pn_v(i,j)
          enddo
          if (iic.eq.ntstart) then
            do k=1,N
              do i=Istr,Iend
                cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i,j-1,k))
                v(i,j,k,nnew)=(v(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+Hz(i,j-1,
     &                                                               k))
     &                                           +cdt*DC(i,0)*rv(i,j,k)
     &                           )*cff
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5D0*(Hz(i,j  ,k)+
     &                                           Hz(i,j-1,k))
              enddo
            enddo
          else
            cff1=0.5D0+gamma
            cff2=0.5D0-gamma
            do k=1,N
              do i=Istr,Iend
                cff = 2.D0/(Hz_half(i,j,k)+Hz_half(i,j-1,k))
                v(i,j,k,nnew)=( cff1*v(i,j,k,nstp)*0.5D0*(Hz(i,j  ,k)+
     &                                                  Hz(i,j-1,k))
     &                         +cff2*v(i,j,k,indx)*0.5D0*(Hz_bak(i,j  
     &                                                              ,k)+
     &                                                  Hz_bak(i,j-1,k))
     &                         +cdt*DC(i,0)*rv(i,j,k)
     &                                                           )*cff
                v(i,j,k,indx)=v(i,j,k,nstp)*0.5D0*(Hz(i,j,k)+
     &                                           Hz(i,j-1,k))
              enddo
            enddo
          endif
        endif
      enddo
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(0,-1,1,nnew,itrc))
      enddo
      call u3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
      call v3dbc_tile (Istr,Iend,Jstr,Jend, WORK)
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(0,-1,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(0,-1,1,nnew))
      do j=Jstr,Jend
        do i=IstrR,IendR
          zeta(i,j,knew)=Zt_avg1(i,j)
        enddo
      enddo
       call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          zeta(0,-1,knew))
      return
      end
