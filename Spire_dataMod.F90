!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
module Spire_dataMod

   ! Defaults
   implicit none
   private

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the Spire soil moisture datasets
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
   ! Public routines
   public :: Spire_datainit

   ! Public types
   public :: Spiredata

   type, public :: spiredatadec
   
      character*100 							:: odir
      integer  										:: spireorbit

 			
			
   end type spiredatadec
   
   type(spiredatadec), allocatable :: spiredata(:)

contains
   !---------------------------------------------------------------------------
   subroutine Spire_datainit(i)
      ! Imports
      use ESMF
      use LVT_coreMod
      use LVT_logMod
      use LVT_histDataMod
      use LVT_timeMgrMod

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: i

      ! Local variables
      integer              :: status
      real                 :: gridDesci(50)
      integer              :: yr,mo,dy,hr,mn,ss

      if(.not.allocated(Spiredata)) then 
         allocate(spiredata(LVT_rc%nDataStreams))
      endif
      ! Get top level Spire data directory
      call ESMF_ConfigGetAttribute(LVT_Config, spiredata(i)%odir, &
           label='Spire data directory:', rc=status)
      call LVT_verify(status, 'Spire data directory: not defined')

      call ESMF_ConfigGetAttribute(LVT_config, spiredata(i)%spireorbit, &
      		 label='Spire orbit (109 or 110):', rc=status)
      call LVT_verify(status, 'Spire orbit: not defined')

    	call LVT_update_timestep(LVT_rc, 3600)  	! Spire data is hourly

   
		 
 
   end subroutine Spire_datainit
end module Spire_dataMod
