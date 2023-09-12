      subroutine setup_kwds (ierr)
      implicit none
      integer*4 ierr, is,ie
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
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      is=1
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='title'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='time_stepping'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='S-coord'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='initial'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='restart'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='history'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='primary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='auxiliary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='rho0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='lateral_visc'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bottom_drag'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='gamma2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='tracer_diff2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='vertical_mixing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bodyforce'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='lin_EOS_cff'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      return
  99   write(stdout,'(/1x,A,A/14x,A)')
     &  'SETUP_KWDS ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
