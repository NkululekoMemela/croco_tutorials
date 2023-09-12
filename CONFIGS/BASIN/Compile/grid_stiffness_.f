      subroutine grid_stiffness (tile)
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
      integer*4 tile, trd
C$    integer*4 omp_get_thread_num
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
      call grid_stiffness_tile (Istr,Iend,Jstr,Jend)
      return
      end
      subroutine grid_stiffness_tile (Istr,Iend,Jstr,Jend)
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, NSUB
      real my_rx0, my_rx1
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
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
      my_rx0=0.D0
      my_rx1=0.D0
      do j=Jstr,Jend
        do i=IstrU,Iend
            my_rx0=max(my_rx0, abs( (z_w(i,j,0)-z_w(i-1,j,0))
     &                                   /(z_w(i,j,0)+z_w(i-1,j,0))
     &                                                            ))
            do k=1,N
              my_rx1=max(my_rx1, abs(
     &         (z_w(i,j,k)-z_w(i-1,j,k)+z_w(i,j,k-1)-z_w(i-1,j,k-1))
     &        /(z_w(i,j,k)+z_w(i-1,j,k)-z_w(i,j,k-1)-z_w(i-1,j,k-1))
     &                                                            ))
            enddo
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
            my_rx0=max(my_rx0, abs( (z_w(i,j,0)-z_w(i,j-1,0))
     &                                   /(z_w(i,j,0)+z_w(i,j-1,0))
     &                                                            ))
            do k=1,N
              my_rx1=max(my_rx1, abs(
     &         (z_w(i,j,k)-z_w(i,j-1,k)+z_w(i,j,k-1)-z_w(i,j-1,k-1))
     &        /(z_w(i,j,k)+z_w(i,j-1,k)-z_w(i,j,k-1)-z_w(i,j-1,k-1))
     &                                                            ))
            enddo
        enddo
      enddo
      if (Iend-Istr+Jend-Jstr.eq.Lm+Mm-2) then
        NSUB=1
      else
        NSUB=NSUB_X*NSUB_E
      endif
C$OMP CRITICAL (grd_stff_cr_rgn)
        if (tile_count.eq.0) then
          rx0=my_rx0
          rx1=my_rx1
        else
          rx0=max(rx0, my_rx0)
          rx1=max(rx1, my_rx1)
        endif
        tile_count=tile_count+1
        if (tile_count.eq.NSUB) then
          tile_count=0
           write(stdout,*)
     &     'Maximum grid stiffness ratios:   rx0 =',rx0, 'rx1 =',rx1
        endif
C$OMP END CRITICAL (grd_stff_cr_rgn)
      return
      end
