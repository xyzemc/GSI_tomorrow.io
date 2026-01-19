subroutine read_tms(mype,val_tovs,ithin,isfcalc,&
     rmesh,jsatid,gstime,infile,lunout,obstype,&
     nread,ndata,nodata,twind,sis, &
     mype_root,mype_sub,npe_sub,mpi_comm_sub,nobs, &
     nrec_start,dval_use,radmod)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_tms                  read tms 1b data from tomorrow.io
!   prgmmr: eliu          org: np23                date: 2024-08-26
!
! abstract:  This routine reads BUFR format TMS radiance 
!            (brightness temperature) files.
!
!            ptionally, the data 
!            are thinned to a specified resolution using simple 
!            quality control checks.
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!  2024-08-26  eliu
!
!   input argument list:
!     mype     - mpi task id
!     val_tovs - weighting factor applied to super obs
!     ithin    - flag to thin data
!     isfcalc  - method to calculate surface fields within FOV
!                when one, calculate accounting for size/shape of FOV.
!                otherwise, use bilinear interpolation.
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window(hours)
!     sis      - sensor/instrument/satellite indicator
!     mype_root - "root" task for sub-communicator
!     mype_sub - mpi task id within sub-communicator
!     npe_sub  - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!     nrec_start - first subset with useful information
!
!   output argument list:
!     nread    - number of BUFR TMS 1b observations read
!     ndata    - number of BUFR TMS 1b profiles retained for further processing
!     nodata   - number of BUFR TMS 1b observations retained for further processing
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,destroygrids,checkob, &
      finalcheck,map2tgrid,score_crit
  use satthin, only: radthin_time_info,tdiff2crit
  use obsmod,  only: time_window_max, ta2tb
  use radinfo, only: iuse_rad,newchn,cbias,nusis,jpch_rad,air_rad,ang_rad, &
      use_edges,radedge1,radedge2,nusis,radstart,radstep,newpc4pred,maxscan
  use radinfo, only: adp_anglebc
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,tll2xy,txy2ll,rlats,rlons
  use constants, only: deg2rad,zero,one,two,three,rad2deg,r60inv,r100,rearth_equator
  use crtm_module, only : max_sensor_zenith_angle
  use calc_fov_crosstrk, only : instrument_init, fov_cleanup, fov_check
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen
  use deter_sfc_mod, only: deter_sfc_fov,deter_sfc
  use atms_spatial_average_mod, only : atms_spatial_average
  use gsi_nstcouplermod, only: nst_gsi,nstinfo
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth,gsi_nstcoupler_deter
  use mpimod, only: npe
  use radiance_mod, only: rad_obs_type

  implicit none

! Declare passed variables
  character(len=*),intent(in   ) :: infile,obstype,jsatid
  character(len=20),intent(in  ) :: sis
  integer(i_kind) ,intent(in   ) :: mype,lunout,ithin
  integer(i_kind) ,intent(in   ) :: nrec_start
  integer(i_kind) ,intent(inout) :: isfcalc
  integer(i_kind) ,intent(inout) :: nread
  integer(i_kind),dimension(npe) ,intent(inout) :: nobs
  integer(i_kind) ,intent(  out) :: ndata,nodata
  real(r_kind)    ,intent(in   ) :: rmesh,gstime,twind
  real(r_kind)    ,intent(inout) :: val_tovs
  integer(i_kind) ,intent(in   ) :: mype_root
  integer(i_kind) ,intent(in   ) :: mype_sub
  integer(i_kind) ,intent(in   ) :: npe_sub
  integer(i_kind) ,intent(in   ) :: mpi_comm_sub
  logical         ,intent(in   ) :: dval_use
  type(rad_obs_type),intent(in ) :: radmod

! Declare local parameters

  character(8),parameter:: fov_flag="crosstrk"
  integer(i_kind),parameter:: n1bhdr=14
  integer(i_kind),parameter:: n2bhdr=5
!emily integer(i_kind),parameter:: maxobs = 800000
  integer(i_kind),parameter:: maxobs = 5000000
  integer(i_kind),parameter:: max_chanl = 12 
  real(r_kind),parameter:: r360=360.0_r_kind
  real(r_kind),parameter:: tbmin=50.0_r_kind
  real(r_kind),parameter:: tbmax=550.0_r_kind
  ! The next two are one minute in hours
  real(r_kind),parameter:: one_minute=0.01666667_r_kind
  real(r_kind),parameter:: minus_one_minute=-0.01666667_r_kind

! Declare local variables
  logical outside,iuse,assim,valid

  character(8) subset
  character(80) hdr1b,hdr2b

  integer(i_kind) ireadsb,ireadmg,irec
  integer(i_kind) i,j,k,ntest,iob,llll
  integer(i_kind) iret,idate,nchanl,n,idomsfc(1)
  integer(i_kind) ich1,ich2,ich3,ich4,ich5,ich6
  integer(i_kind) ich7,ich8,ich9,ich10,ich11,ich12
  integer(i_kind) kidsat,kidsatsub,maxinfo  !emily
  integer(i_kind) nmind,itx,nreal,nele,itt,num_obs
  integer(i_kind) iskip,scanline
  integer(i_kind) lnbufr,ksatid,ksatidsub,isflg !emily
  integer(i_kind) ilat,ilon, nadir
  integer(i_kind),dimension(5):: idate5
  integer(i_kind) instr,ichan
  integer(i_kind):: ierr
  integer(i_kind):: radedge_min, radedge_max
  integer(i_kind), POINTER :: ifov
  integer(i_kind), TARGET :: ifov_save(maxobs)
  integer(i_kind), ALLOCATABLE :: IScan(:)

  real(r_kind) cosza,sfcr
  real(r_kind) ch1,ch2,ch3,d0,d1,d2,ch16,qval
  real(r_kind) expansion
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr
  real(r_kind) :: satellite_height, rato

  real(r_kind) pred
  real(r_kind) dlat,scanang,dlon,tdiff
  real(r_kind) dlon_earth_deg,dlat_earth_deg,r01
  real(r_kind) step,start,dist1
  real(r_kind) tt,lzaest
  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),allocatable :: Relative_Time_In_Seconds(:)
  real(r_kind),allocatable,dimension(:,:):: data_all
  real(r_kind), POINTER :: bt_in(:), crit1,rsat, t4dv, solzen, solazi
  real(r_kind), POINTER :: dlon_earth, dlat_earth, satazi, lza, panglr

  integer(i_kind), ALLOCATABLE, TARGET :: it_mesh_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: rsat_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: t4dv_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: dlon_earth_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: dlat_earth_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: crit1_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: lza_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: panglr_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: satazi_save(:)
  real(r_kind), ALLOCATABLE, TARGET :: solzen_save(:) 
  real(r_kind), ALLOCATABLE, TARGET :: solazi_save(:) 
  real(r_kind), ALLOCATABLE, TARGET :: bt_save(:,:)

  integer(i_kind),allocatable,dimension(:):: nrec
  real(r_double),allocatable,dimension(:):: data1b8
  real(r_double),dimension(n1bhdr):: bfr1bhdr
  real(r_double),dimension(n2bhdr):: bfr2bhdr

  real(r_kind) cdist,disterr,disterrmax,dlon00,dlat00

  logical :: critical_channels_missing
  real(r_kind)    :: ptime,timeinflat,crit0
  integer(i_kind) :: ithin_time,n_tbin
  integer(i_kind),pointer :: it_mesh => null()

  real(r_double),dimension(12):: flags !xzhang
  integer(i_kind), ALLOCATABLE, TARGET :: qc_flag(:,:) !xzhang
  real(r_double), parameter    :: Missing_Value=1.e11_r_double !xzhang
  logical :: tms_qcflag=.true. !xzhang
  integer(i_kind),parameter:: mxib=100
  integer(i_kind) ibit(mxib),nib,n_bad


!**************o***********************************************************
! Initialize variables

  maxinfo=32
  lnbufr = 15
  disterrmax=zero
  ntest=0
  ndata  = 0
  nodata  = 0
  nread  = 0

  ilon=3
  ilat=4

  if(nst_gsi>0) then
     call gsi_nstcoupler_skindepth(obstype,zob)
  endif

  call radthin_time_info(obstype, jsatid, sis, ptime, ithin_time)
  if( ptime > 0.0_r_kind) then
     n_tbin=nint(2*time_window_max/ptime)
  else
     n_tbin=1
  endif
! Make thinning grids
  call makegrids(rmesh,ithin,n_tbin=n_tbin)

! Set nadir position based on value of maxscan
  nadir=41

! Set various variables depending on type of data to be read

  if (obstype /= 'tms') then
     write(6,*) 'READ_TMS called for obstype '//obstype//': RETURNING'
     return
  end if

!  instrument specific variables
  d1 =  0.754_r_kind
  d2 = -2.265_r_kind 
  r01 = 0.01_r_kind
  ich1   = 1   !1
  ich2   = 2   !2
  ich3   = 3   !3
  ich4   = 4   !4
  ich5   = 5   !5
  ich6   = 6   !6
  ich7   = 7   !7
  ich8   = 8   !8
  ich9   = 9   !9
  ich10  = 10  !10
  ich11  = 11  !11
  ich12  = 12  !12

  ! Set array index for surface-sensing channels

  if(jsatid == 'tomorrow-s01_v4') then
     kidsat    = 769
     kidsatsub = 101
  else if (jsatid == 'tomorrow-s02_v4') then
     kidsat    = 769
     kidsatsub = 102
  else if (jsatid == 'tomorrow-s03_v4') then
     kidsat    = 769
     kidsatsub = 103
  else if (jsatid == 'tomorrow-s04_v4') then
     kidsat    = 769
     kidsatsub = 104
  else if (jsatid == 'tomorrow-s05_v4') then
     kidsat    = 769
     kidsatsub = 105
  else if (jsatid == 'tomorrow-s06_v4') then
     kidsat    = 769
     kidsatsub = 106
  else if (jsatid == 'tomorrow-s07_v4') then
     kidsat    = 769
     kidsatsub = 107
  else 
     write(6,*) 'READ_TMS: Unrecognized value for jsatid '//jsatid//': RETURNING'
     return
  end if

  radedge_min = 0
  radedge_max = 1000
  do i=1,jpch_rad
     if (trim(nusis(i))==trim(sis)) then
        step  = radstep(i)
        start = radstart(i)
        if (radedge1(i)/=-1 .or. radedge2(i)/=-1) then
           radedge_min=radedge1(i)
           radedge_max=radedge2(i)
        end if
        exit 
     endif
  end do 

! Allocate arrays to hold all data for given satellite
  nchanl=12
  if(dval_use) maxinfo = maxinfo+2
  nreal = maxinfo + nstinfo
  if(tms_qcflag) nreal = maxinfo+nstinfo+nchanl
  nele  = nreal   + nchanl
  allocate(data_all(nele,itxmax),nrec(itxmax))
  nrec=999999

! IFSCALC setup
  if (isfcalc==1) then
     instr=20                    
     ichan=16                    ! pick a surface sens. channel
     expansion=2.9_r_kind        ! use almost three for microwave sensors.
  endif
! Set rlndsea for types we would prefer selecting
  rlndsea(0) = zero
  rlndsea(1) = 15._r_kind
  rlndsea(2) = 10._r_kind
  rlndsea(3) = 15._r_kind
  rlndsea(4) = 100._r_kind
     
! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_tovs=zero

! Initialize variables for use by FOV-based surface code.
  if (isfcalc == 1) then
     call instrument_init(instr,jsatid,expansion,valid)
     if (.not. valid) then
       if (assim) then
         write(6,*)'READ_TMS:  ***ERROR*** IN SETUP OF FOV-SFC CODE. STOP'
         call stop2(71)
       else
         call fov_cleanup
         isfcalc = 0
         write(6,*)'READ_TMS:  ***ERROR*** IN SETUP OF FOV-SFC CODE'
       endif
     endif
  endif

! This eventually needs to be spit between MHS and AMSU-A like channels
  if (isfcalc==1) then
!   if (amsub.or.mhs)then
!     rlndsea(4) = max(rlndsea(0),rlndsea(1),rlndsea(2),rlndsea(3))
!   else
      rlndsea=0
!   endif
  endif

! Allocate arrays for BUFR I/O
  ALLOCATE(data1b8(nchanl))
  ALLOCATE(rsat_save(maxobs))
  ALLOCATE(t4dv_save(maxobs))
  ALLOCATE(dlon_earth_save(maxobs))
  ALLOCATE(dlat_earth_save(maxobs))
  ALLOCATE(crit1_save(maxobs))
  ALLOCATE(it_mesh_save(maxobs))
  ALLOCATE(lza_save(maxobs))
  ALLOCATE(panglr_save(maxobs))
  ALLOCATE(satazi_save(maxobs))
  ALLOCATE(solzen_save(maxobs)) 
  ALLOCATE(solazi_save(maxobs)) 
  ALLOCATE(bt_save(max_chanl,maxobs))
  ALLOCATE(qc_flag(max_chanl,maxobs))

  qc_flag=0
! Read in data from bufr into arrays first      
! Open unit to satellite bufr file
  iob=1
  n_bad=0
  open(lnbufr,file=trim(infile),form='unformatted',status = 'old', iostat = ierr)
  call openbf(lnbufr,'IN',lnbufr)
  call datelen(10)

  write(6,*)'READ_TMS: emily checking - reading ', trim(infile)

  hdr1b ='SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON CLATH CLONH HMSL SASBID'
  hdr2b ='SAZA SOZA BEARAZ SOLAZI SANGST'
! hdr2b ='SAZA SOZA BEARAZ SOLAZI'
   
! Loop to read bufr file
  irec=0
  read_subset: do while(ireadmg(lnbufr,subset,idate)>=0 .AND. iob < maxobs)
     irec = irec + 1
     if (irec < nrec_start) cycle read_subset
     read_loop: do while (ireadsb(lnbufr)==0 .and. iob < maxobs)

        rsat       => rsat_save(iob)
        t4dv       => t4dv_save(iob)
        dlon_earth => dlon_earth_save(iob)
        dlat_earth => dlat_earth_save(iob)
        crit1      => crit1_save(iob)
        it_mesh    => it_mesh_save(iob)
        ifov       => ifov_save(iob)
        lza        => lza_save(iob)
        panglr     => panglr_save(iob)
        satazi     => satazi_save(iob)
        solzen     => solzen_save(iob)
        solazi     => solazi_save(iob)

!       initialize inflate selection value
        crit0 = 0.01_r_kind

        call ufbint(lnbufr,bfr1bhdr,n1bhdr,1,iret,hdr1b)

!       Extract satellite id.  If not the one we want, read next record
        rsat=bfr1bhdr(1) 
        ksatid=nint(bfr1bhdr(1))
        ksatidsub=nint(bfr1bhdr(14))
!        write(6,*)'emily checking ksatid    = ', ksatid
!        write(6,*)'emily checking kidsat    = ', kidsat
!        write(6,*)'emily checking ksatidsub = ', ksatidsub
!        write(6,*)'emily checking kidsatsub = ', kidsatsub 

        if(ksatid /= kidsat) cycle read_subset
        if(ksatidsub /= kidsatsub) cycle read_subset

!       Extract observation location and other required information
        if(abs(bfr1bhdr(11)) <= 90._r_kind .and. abs(bfr1bhdr(12)) <= r360)then
           dlat_earth = bfr1bhdr(11)
           dlon_earth = bfr1bhdr(12)
        elseif(abs(bfr1bhdr(9)) <= 90._r_kind .and. abs(bfr1bhdr(10)) <= r360)then
           dlat_earth = bfr1bhdr(9)
           dlon_earth = bfr1bhdr(10)
        else
           cycle read_loop
        end if
        if(dlon_earth<zero)  dlon_earth = dlon_earth+r360
        if(dlon_earth>=r360) dlon_earth = dlon_earth-r360

!       Extract date information.  If time outside window, skip this obs
        idate5(1) = bfr1bhdr(3) !year
        idate5(2) = bfr1bhdr(4) !month
        idate5(3) = bfr1bhdr(5) !day
        idate5(4) = bfr1bhdr(6) !hour
        idate5(5) = bfr1bhdr(7) !minute
        call w3fs21(idate5,nmind)
        t4dv= (real((nmind-iwinbgn),r_kind) + bfr1bhdr(8)*r60inv)*r60inv    ! add in seconds
        tdiff=t4dv+(iwinbgn-gstime)*r60inv

        if (l4dvar.or.l4densvar) then
           if (t4dv<minus_one_minute .OR. t4dv>winlen+one_minute) &
               cycle read_loop
        else
           if(abs(tdiff) > twind+one_minute) cycle read_loop
        endif

        timeinflat=two
        call tdiff2crit(tdiff,ptime,ithin_time,timeinflat,crit0,crit1,it_mesh)


!       Get various angles 
        call ufbint(lnbufr,bfr2bhdr,n2bhdr,1,iret,hdr2b)

        satazi=bfr2bhdr(3)
        if (abs(satazi) > r360) then
           satazi=zero
        endif
        scanang=bfr2bhdr(5)            ! sensor scan angle
        ifov = nint(bfr1bhdr(2))       ! field of view number
        lza = bfr2bhdr(1)*deg2rad      ! local zenith angle
      ! if(ifov <= nadir)  lza=-lza    ! emily check here

        panglr=scanang*deg2rad        
        write(6,'(a28, 2x, i6, 2x, 2(f12.5,2x))') &
             'READ_TMS ifov lza panglr ', ifov, abs(lza)*rad2deg, panglr*rad2deg

        if(abs(lza)*rad2deg > MAX_SENSOR_ZENITH_ANGLE) then  ! MAX_SENSOR_ZENITH_ANGLE = 80
           write(6,'(a28, 2x, i6, 2x, 2(f12.5,2x))') &
                'READ_TMS WARNING lza error ',ifov, abs(lza)*rad2deg, panglr*rad2deg
           cycle read_loop
        end if

        solzen_save(iob)=bfr2bhdr(2) 
        solazi_save(iob)=bfr2bhdr(4) 
!       Read TMSTBR flags for all channels
        call ufbrep(lnbufr, flags, 1, nchanl, iret, 'TMSF')
       
        do i = 1,12
           call upftbv(lnbufr,'TMSF',flags(i),mxib,ibit,nib)
           if (nib > 0 )then
             do j=1,nib
               if (ibit(j) == 23) then !v1 before 12/3/2025
               !if (ibit(j) == 9) then 
                 qc_flag(i,iob) = 1
                 n_bad=n_bad+1
               end if
             end do
           end if
        end do

!       Read data record.  Increment data counter
        call ufbrep(lnbufr,data1b8,1,nchanl,iret,'TMBR')


        bt_save(1:nchanl,iob) = data1b8(1:nchanl)


        iob=iob+1

     end do read_loop
  end do read_subset
  call closbf(lnbufr)
  close(lnbufr)
  deallocate(data1b8)

  num_obs = iob-1
  write(6,*) 'READ_TMS: emily checking num_obs = ', num_obs

  if (num_obs <= 0) then
     write(6,*) 'READ_TMS: No TMS Data were read in'
     return
  end if

! Call filtering code 

  ALLOCATE(Relative_Time_In_Seconds(Num_Obs))
  ALLOCATE(IScan(Num_Obs))
  Relative_Time_In_Seconds = 3600.0_r_kind*T4DV_Save(1:Num_Obs)
! write(6,*) 'Calling ATMS_Spatial_Average'
!  CALL ATMS_Spatial_Average(Num_Obs, NChanl, IFOV_Save(1:Num_Obs), &
!       Relative_Time_In_Seconds, BT_Save(1:nchanl,1:Num_Obs), IScan, IRet)
! write(6,*) 'ATMS_Spatial_Average Called with IRet=',IRet
  DEALLOCATE(Relative_Time_In_Seconds)
  
!  IF (IRet /= 0) THEN
!     write(6,*) 'Error Calling ATMS_Spatial_Average from READ_TMS'
!     RETURN
!  END IF

! Complete read_tms thinning and QC steps

  ObsLoop: do iob = 1, num_obs  

     rsat       => rsat_save(iob)
     t4dv       => t4dv_save(iob)
     dlon_earth => dlon_earth_save(iob)
     dlat_earth => dlat_earth_save(iob)
     crit1      => crit1_save(iob)
     it_mesh    => it_mesh_save(iob)
     ifov       => ifov_save(iob)
     lza        => lza_save(iob)
     panglr     => panglr_save(iob)
     satazi     => satazi_save(iob)
     solzen     => solzen_save(iob)
     solazi     => solazi_save(iob)
     bt_in      => bt_save(1:nchanl,iob)
     
     dlat_earth_deg = dlat_earth
     dlon_earth_deg = dlon_earth
     dlat_earth = dlat_earth*deg2rad
     dlon_earth = dlon_earth*deg2rad   

!    Regional case
     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if(diagnostic_reg) then
           call txy2ll(dlon,dlat,dlon00,dlat00)
           ntest=ntest+1
           cdist=sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00))
           cdist=max(-one,min(cdist,one))
           disterr=acos(cdist)*rad2deg
           disterrmax=max(disterrmax,disterr)
        end if
           
!       Check to see if in domain
        if(outside) cycle ObsLoop
           
!    Global case
     else
        dlat=dlat_earth
        dlon=dlon_earth
        call grdcrd1(dlat,rlats,nlat,1)
        call grdcrd1(dlon,rlons,nlon,1)
     endif

! Check time window
     if (l4dvar.or.l4densvar) then
        if (t4dv<zero .OR. t4dv>winlen) cycle ObsLoop
     else
        tdiff=t4dv+(iwinbgn-gstime)*r60inv
        if(abs(tdiff) > twind) cycle ObsLoop
     endif
 
!    Map obs to thinning grid
     call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis,it_mesh=it_mesh)
     if(.not. iuse)cycle ObsLoop

!
!    Check FOV and scan-edge usage
     if (.not. use_edges .and. (ifov < radedge_min .OR. ifov > radedge_max )) &
          cycle ObsLoop

     nread=nread+nchanl
     
!    Transfer observed brightness temperature to work array.  If any
!    temperature exceeds limits, reset observation to "bad" value
     iskip=0
     critical_channels_missing = .false.
     do j=1,nchanl
        if (bt_in(j) < tbmin .or. bt_in(j) > tbmax) then
           iskip = iskip + 1
           
!          Flag profiles where key channels are bad  
           if(j == ich1 .or. j == ich12)  critical_channels_missing = .true.
        endif
     end do
!emily     if (iskip >= nchanl) cycle ObsLoop

!    Determine surface properties based on 
!    sst and land/sea/ice mask   
!
!    isflg    - surface flag
!               0 sea
!               1 land
!               2 sea ice
!               3 snow
!               4 mixed                       

!    FOV-based surface code requires fov number.  if out-of-range, then
!    skip this ob.

     if (isfcalc == 1) then
        call fov_check(ifov,instr,ichan,valid)
        if (.not. valid) cycle ObsLoop

!    When isfcalc is one, calculate surface fields based on size/shape of fov.
!    Otherwise, use bilinear method.

        call deter_sfc_fov(fov_flag,ifov,instr,ichan,satazi,dlat_earth_deg,&
             dlon_earth_deg,expansion,t4dv,isflg,idomsfc(1), &
             sfcpct,vfr,sty,vty,stp,sm,ff10,sfcr,zz,sn,ts,tsavg)
     else
        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg, &
             idomsfc(1),sfcpct,ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
     endif

     crit1 = crit1 + rlndsea(isflg) + 10._r_kind*real(iskip,r_kind) + 0.01_r_kind * abs(zz)
     call checkob(dist1,crit1,itx,iuse)
     if(.not. iuse)cycle ObsLoop

     if (critical_channels_missing) then
         pred=1.0e8_r_kind
     else
         pred=zero
     endif

!    Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"
     crit1 = crit1+pred 
     call finalcheck(dist1,crit1,itx,iuse)
     if(.not. iuse)cycle ObsLoop
     
!    interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
     if(nst_gsi>0) then
        tref  = ts(0)
        dtw   = zero
        dtc   = zero
        tz_tr = one
        if(sfcpct(0)>zero) then
           call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
        endif
     endif

!    Re-calculate look angle
!    panglr=(start+real(ifov-1,r_kind)*step)*deg2rad

!    Load selected observation into data array
              
     data_all(1 ,itx)= rsat                      ! satellite ID
     data_all(2 ,itx)= t4dv                      ! time
     data_all(3 ,itx)= dlon                      ! grid relative longitude
     data_all(4 ,itx)= dlat                      ! grid relative latitude
     data_all(5 ,itx)= lza                       ! local zenith angle
     data_all(6 ,itx)= satazi                    ! local azimuth angle
     data_all(7 ,itx)= panglr                    ! look angle
     data_all(8 ,itx)= ifov                      ! scan position
     data_all(9 ,itx)= solzen                    ! solar zenith angle
     data_all(10,itx)= solazi                    ! solar azimuth angle
     data_all(11,itx) = sfcpct(0)                ! sea percentage of
     data_all(12,itx) = sfcpct(1)                ! land percentage
     data_all(13,itx) = sfcpct(2)                ! sea ice percentage
     data_all(14,itx) = sfcpct(3)                ! snow percentage
     data_all(15,itx)= ts(0)                     ! ocean skin temperature
     data_all(16,itx)= ts(1)                     ! land skin temperature
     data_all(17,itx)= ts(2)                     ! ice skin temperature
     data_all(18,itx)= ts(3)                     ! snow skin temperature
     data_all(19,itx)= tsavg                     ! average skin temperature
     data_all(20,itx)= vty                       ! vegetation type
     data_all(21,itx)= vfr                       ! vegetation fraction
     data_all(22,itx)= sty                       ! soil type
     data_all(23,itx)= stp                       ! soil temperature
     data_all(24,itx)= sm                        ! soil moisture
     data_all(25,itx)= sn                        ! snow depth
     data_all(26,itx)= zz                        ! surface height
     data_all(27,itx)= idomsfc(1) + 0.001_r_kind ! dominate surface type
     data_all(28,itx)= sfcr                      ! surface roughness
     data_all(29,itx)= ff10                      ! ten meter wind factor
     data_all(30,itx) = dlon_earth_deg           ! earth relative longitude (deg)
     data_all(31,itx) = dlat_earth_deg           ! earth relative latitude (deg)
     data_all(32,itx) = scanline                 ! scan line number
     
     if(dval_use) then
        data_all(33,itx)= val_tovs
        data_all(34,itx)= itt
     end if

     
     if(nst_gsi>0) then
        data_all(maxinfo+1,itx) = tref            ! foundation temperature
        data_all(maxinfo+2,itx) = dtw             ! dt_warm at zob
        data_all(maxinfo+3,itx) = dtc             ! dt_cool at zob
        data_all(maxinfo+4,itx) = tz_tr           ! d(Tz)/d(Tr)
     endif

     if (tms_qcflag) then
       do i=1,nchanl
          data_all(nreal-nchanl+i,itx)=qc_flag(i,iob)
       end do
     end if  

     do i=1,nchanl
        data_all(i+nreal,itx)=bt_in(i)
     end do
     nrec(itx)=iob

  end do ObsLoop
  print*,'iob, n_bad = ', iob, n_bad


  DEALLOCATE(iscan)
! DEAllocate I/O arrays
  DEALLOCATE(rsat_save)
  DEALLOCATE(t4dv_save)
  DEALLOCATE(dlon_earth_save)
  DEALLOCATE(dlat_earth_save)
  DEALLOCATE(crit1_save)
  DEALLOCATE(it_mesh_save)
  DEALLOCATE(lza_save)
  DEALLOCATE(panglr_save)
  DEALLOCATE(satazi_save)
  DEALLOCATE(solzen_save) 
  DEALLOCATE(solazi_save) 
  DEALLOCATE(bt_save)
  DEALLOCATE(qc_flag)

  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
       nele,itxmax,nread,ndata,data_all,score_crit,nrec)

! 
  if(mype_sub==mype_root)then
     do n=1,ndata
        do i=1,nchanl
           if(data_all(i+nreal,n) > tbmin .and. &
                data_all(i+nreal,n) < tbmax)nodata=nodata+1
        end do
     end do

     if(dval_use .and. assim)then
        do n=1,ndata
           itt=nint(data_all(33,n))
           super_val(itt)=super_val(itt)+val_tovs
        end do
     end if
     
!    Write final set of "best" observations to output file
     call count_obs(ndata,nele,ilat,ilon,data_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  end if
     
! Deallocate local arrays
  deallocate(data_all,nrec)

! Deallocate satthin arrays
  call destroygrids

! Deallocate FOV surface code arrays and nullify pointers.
  if (isfcalc == 1) call fov_cleanup

  if(diagnostic_reg.and.ntest>0) write(6,*)'READ_TMS:  ',&
     'mype,ntest,disterrmax=',mype,ntest,disterrmax

  write(6,*)'READ_TMS: emily checking - reading DONE ', trim(infile)
! End of routine
  return

end subroutine read_tms
