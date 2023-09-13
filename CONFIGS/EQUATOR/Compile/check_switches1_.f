      subroutine check_switches1 (ierr)
      implicit none
      integer*4 ierr, is,ie, iexample
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=40,   MMm0=32,   N=32)
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
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
       write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
      iexample=iexample+1
       write(stdout,'(10x,A)') 'EQUATOR'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EQUATOR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'SOLVE3D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLVE3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_COR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_COR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'QCORRECTION'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QCORRECTION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_VIS2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_VADV_SPLINES'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_DIF2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_VADV_AKIMA'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_AKIMA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'LMD_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'LMD_RIMIX'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_RIMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'LMD_CONVEC'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_CONVEC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_BSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_BTFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BTFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_INITIAL'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_INITIAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_SSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_STFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_STFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_GRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_GRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_SMFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SMFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'NO_FRCFILE'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NO_FRCFILE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_SRFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_SST'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'MPI_COMM_WORLD'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI_COMM_WORLD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'M2FILTER_POWER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_POWER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TRACERS'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TRACERS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TEMPERATURE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TEMPERATURE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'HZR'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HZR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'VAR_RHO_2D'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_RHO_2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'RESET_RHO0'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RESET_RHO0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'PGF_FLAT_BOTTOM'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PGF_FLAT_BOTTOM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'NTRA_T3DMIX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NTRA_T3DMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'NF_CLOBBER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_CLOBBER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(/)')
      if (iexample.eq.0) then
         write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
         write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99   write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
