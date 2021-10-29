!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSpireObs
! \label{readSpireObs}
!
! !INTERFACE: 
subroutine readSpireObs(source)
! 
! !USES:   
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, & 
       LVT_releaseUnitNumber
  use LVT_timeMgrMod,   only : LVT_get_julss
  use Spire_dataMod 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for the Spire Soil 
! Moisture dataset 

! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  ! Local variables
  integer, parameter :: = ntiles = 200000		!max samples
  logical           :: file_exists

  character*200     :: filename
  real              :: lat(ntiles)
  real							:: lon(ntiles)
  real							:: lh(ntiles)
  real							:: lh2d(LVT_rc%lnc,LVT_rc%lnr)
  integer						:: nlh2d(LVT_rc%lnc,LVT_rc%lnr)
  
  real              :: sfsm(ntiles)
  real							:: sfsm2d(LVT_rc%lnc,LVT_rc%lnr)
  integer						:: nsfsm2d(LVT_rc%lnc,LVT_rc%lnr)

  integer						:: latid,lonid
  integer						:: lhid, sfsmid
  real							:: col,row
  integer           :: c,r
  integer           :: t, iret, ftn

	integer           :: yr, mo, dy, hr, mn, ss
  

  ! Initialize variables
  lh2d		= 0
  nlhd2d	= 0
  
  sfsm2d	= 0
  nsfsm2d	= 0

  
  ! Define datetime 
  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  dy = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  mn = 0
  ss = 0
  
  ! Create filename for reading
  call create_Spire_filename(spiredata(source)%odir, &
                             yr,mo,dy,hr,filename)
  inquire(file=trim(filename),exist=file_exists)

	! Read Spire data
  if(file_exists) then 
  	write(LVT_logunit,*) '[INFO] Reading Spire data file ',trim(filename)
  	
		#if(defined USE_NETCDF3 || defined USE_NETCDF4)
			iret = nf90_open(path=trim(filename), mode=NF90_NOWRITE, ncid=ftn)
											
			if (iret.eq.0) then
		
				call LVT_verify(nf90_inq_varid(ftn,"sp_lat",latid),&
							'nf90_inq_varid failed for lat')								! variable sp_lat
				call LVT_verify(nf90_inq_varid(ftn,"sp_lon",lonid),&
							'nf90_inq_varid failed for lon')								! variable sp_lon
				call LVT_verify(nf90_inq_varid(ftn,"cssm",sfsmid),&
							'nf90_inq_varid failed for SFSM')								! variable cssm

				call LVT_verify(nf90_get_var(ftn,latid,lat),&
							'nf90_inq_varid failed for lat')								! variable sp_lat
				call LVT_verify(nf90_inq_varid(ftn,lonid,lon),&
							'nf90_inq_varid failed for lon')								! variable sp_lon
				call LVT_verify(nf90_inq_varid(ftn,sfsmid,sfsm),&
							'nf90_inq_varid failed for SFSM')								! variable cssm
			
				do t=1,ntiles
					call latlon_to_ij(LVT_domain%lvtproj,&
														lat(t),lon(t),col,row)
					c=nint(col)
					r=nint(row)
					if(c.ge.1 .AND. c.le.LVT_rc%lnc .AND. &
					 r.ge.1 .AND. r.le.LVT_rc%lnr) then
					 lh2d(c,r) = lh2d(c,r) + lh(t)
					 nlh2d(c,r) = nlh2d(c,r) + 1
				 
					 sfsm2d(c,r) = sfsm2d(c,r) + sfsm(t)
					 nsfsm2d(c,r) = nsfsm2d(c,r) + 1
				 
					endif 
				enddo
			
			endif
			iret = nf90_close(ftn)
		
		#endif
	else
		write(LVT_logunit,*) '[WARN] Warning: Spire data file ', trim(filename), 'does not exist' 
		lh = -9999.0
	endif
			
	do r=1,LVT_rc%lnr
		do c=1,LVT_rc%lnc
			if(nlh2d(c,r).gt.0) then
				lh2d(c,r) = lh2d(c,r)/nlh2d(c,r)
				sfsm2d(c,r) = sfsm2d(c,r)/nsfsm2d(c,r)
			else
				lh2d(c,r) = LVT_rc%udef
				sfsm2d(c,r) = LVT_rc%udef
			endif
		enddo
	enddo
	 
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source, &
       sfsm2d,vlevel=1,units="cm3/cm3")

end subroutine readSpireObs


!------------------------------------------------------------------------------

!BOP
! !ROUTINE: create_Spire_filename
! \label{create_Spire_filename}
! 
! !INTERFACE: 
subroutine create_Spire_filename(odir, &
                             yr,mo,dy,hr,filename)
! !USES:   
  use Spire_dataMod
  use LVT_logMod
  
  
  implicit none
  
  
	! Arguments
	character(len=*), intent(in) :: odir 
	integer, intent(in) :: yr, mo, dy, hr
	character(len=*), intent(out) :: filename

	! Local variables
	integer :: uyr, umo, udy, uhr, uyradd, umoadd,udyadd, uhradd 
	character*4   :: cyr,cyradd
	character*2   :: cmo, cdy, chr, cmoadd,cdyadd, chradd 
  character*100 :: forb, fstem 
	
	uyr = yr
	umo = mo
	udy = dy
  uhr = 1*(hr/1)  !hour needs to be a multiple of 1 hour
	
	! Calculate the end datetime for Spire's hourly increment
	
	if (uhr == 23) then				!end of day
		if (udy == 31) then			!end of month
			if (umo == 12) then		!end of year
				uhradd = 0
				udyadd = 1
				umoadd = 1
				uyradd = uyr + 1
			else 		!Months 1,3,5,7,8,10, of 31 days
				uhradd = 0
				udyadd = 1
				umoadd = umo + 1 
				uyradd = uyr
			endif
		elseif (udy == 30) then			 
			if (umo == 4 .OR. umo == 6 .OR. umo == 9 .OR. umo == 11) then	!Months of 30 days
				uhradd = 0
				udyadd = 1
				umoadd = umo + 1
				uyradd = uyr
			else 		!Months 1,3,5,7,8,10, of 31 days
				uhradd = 0
				udyadd = udy + 1
				umoadd = umo 
				uyradd = uyr
			endif
		elseif ((udy == 28 .OR. udy == 29) .AND. umo == 2) then  !Feb-end 
			uhradd = 0
			udyadd = 1
			umoadd = umo + 1
			uyradd = uyr
		else			!not the end of month
			uhradd = 0
			udyadd = udy + 1
			umoadd = umo 
			uyradd = uyr
		endif

	else				!hr < 23
		uhradd = uhr + 1
		udyadd = udy 
		umoadd = umo 
		uyradd = uyr
	endif

	

	write(cyr, '(I4.4)') uyr
	write(cmo, '(I2.2)') umo
	write(cdy, '(I2.2)') udy
	write(chr, '(I2.2)') uhr
	write(cyradd, '(I4.4)') uyradd 
	write(cmoadd, '(I2.2)') umoadd
	write(cdyadd, '(I2.2)') udyadd
	write(chradd, '(I2.2)') uhradd


	! Create filename:
	! odir = /css/csda-spire/gnss-r-0.3.14_evaluation/L2-land
	! orbit,date = FM109/2021-03-31/
	! file = Spire_GNSSR_L2A-LAND_FM109_C0_s2021-03-31T17-00-00Z_e2021-03-31T18-00-00Z_v0.5.0.nc

	if(spiredata(1)%spireorbit == 109) then
		forb = 'FM109'
		fstem = 'Spire_GNSSR_L2A-LAND_FM109_C0_s'
	elseif(spiredata(1)%spireorbit == 110) then 
		forb = 'FM110'
		fstem = 'Spire_GNSSR_L2A-LAND_FM110_C0_s'
	else
		write(LVT_logunit,*) "[ERR] Invalid Spire orbit was chosen."
		write(LVT_logunit,*) "[ERR] Please choose orbit as either '109' or '110'."
		call LVT_endrun()
	endif
  
  
  filename = trim(odir)//'/'//&
       			 trim(forb)//'/'//&
       			 cyr//'-'//cmo//'-'//cdy//'/'//&
       			 trim(fstem)//cyr//'-'//cmo//'-'//cdy//'T'//chr//'-00-00Z_e'//&
       			 cyradd//'-'//cmoadd//'-'//cdyadd//'T'//chradd//'-00-00Z_v0.5.0.nc'
       

end subroutine create_Spire_filename
