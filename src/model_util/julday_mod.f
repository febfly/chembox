! $id: julday_mod.f,v 1.1 2009/11/20 21:43:04 bmy exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: julday_mod.f
!
! !DESCRIPTION: Module JULDAY\_MOD contains routines used to convert from 
!  month/day/year to Astronomical Julian Date and back again.
!\\
!\\
! !INTERFACE: 
!
      module julday_mod
! 
! !USES:
!
      implicit none
      private
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  :: julday
      public  :: caldate
!
! !PRIVATE MEMBER FUNCTIONS:
!
      private :: mint
!
! !REVISION HISTORY:
!  (1 ) Moved JULDAY, MINT, CALDATE here from "bpch2_mod.f" (bmy, 11/20/01)
!  (2 ) Bug fix: now compute NHMS correctly.  Also use REAL*4 variables to
!        avoid roundoff errors. (bmy, 11/26/01)
!  (3 ) Updated comments (bmy, 5/28/02)
!  (4 ) Renamed arguments for clarity (bmy, 6/26/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX Headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      contains
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: julday  
!
! !DESCRIPTION: Function JULDAY returns the astronomical Julian day. 
!\\
!\\
! !INTERFACE:
!
      function julday( yyyy, mm, dd ) result( julianday )
!
! !INPUT PARAMETERS: 
!
      integer, intent(in) :: yyyy        ! year (must be in 4-digit format!)
      integer, intent(in) :: mm          ! month (1-12)
      real*8,  intent(in) :: dd          ! day of month (may be fractional!)
!
! !RETURN VALUE:
! 
      real*8              :: julianday   ! astronomical julian date
!
! !REMARKS:
!  (1) Algorithm taken from "Practical Astronomy With Your Calculator",
!       Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
!  (2) Requires the external function MINT.F. 
!  (3) JulDay will compute the correct Julian day for any BC or AD date.
!  (4) For BC dates, subtract 1 from the year and append a minus sign.  
!       For example, 1 BC is 0, 2 BC is -1, etc.  This is necessary for 
!       the algorithm.  
!
! !REVISION HISTORY: 
!  26 Nov 2001 - R. Yantosca - Initial version
!  Changed YEAR to YYYY, MONTH to MM, and DAY to DD for documentation
!   purposes. (bmy, 6/26/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
      integer             :: year1, month1
      real*8              :: x1, a, b, c, d
      logical             :: isgregorian
   
      !==================================================================
      ! JULDAY begins here!
      !
      ! Follow algorithm from Peter Duffett-Smith (1992)
      !==================================================================
   
      ! Compute YEAR and MONTH1
      if ( ( mm == 1 ) .or. ( mm == 2 ) ) then
         year1  = yyyy  - 1
         month1 = mm    + 12 
      else
         year1  = yyyy
         month1 = mm
      endif
   
      ! Compute the "A" term. 
      x1 = dble( yyyy ) / 100.0d0
      a  = mint( x1 )
   
      ! The Gregorian calendar begins on 10 October 1582
      ! Any dates prior to this will be in the Julian calendar
      if ( yyyy > 1582 ) then
         isgregorian = .true.
      else
         if ( ( yyyy   == 1582 )  .and. 
     &        ( month1 >= 10   )  .and. 
     &        ( dd     >= 15.0 ) ) then 
            isgregorian = .true.
         else
            isgregorian = .false.
         endif
      endif
            
      ! Compute the "B" term according to Gregorian or Julian calendar
      if ( isgregorian ) then
         b = 2.0d0 - a + mint( a / 4.0d0 )
      else
         b = 0.0d0
      endif
   
      ! Compute the "C" term for BC dates (YEAR1 <= 0 ) 
      ! or AD dates (YEAR1 > 0)
      if ( year1 < 0 ) then
         x1 = ( 365.25d0 * year1 ) - 0.75d0
         c  = mint( x1 )
      else
         x1 = 365.25d0 * year1
         c  = mint( x1 ) 
      endif
   
      ! Compute the "D" term    
      x1 = 30.6001d0 * dble( month1 + 1 )
      d  = mint( x1 )
   
      ! Add the terms to get the Julian Day number 
      julianday = b + c + d + dd + 1720994.5d0
   
      end function julday
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mint
!
! !DESCRIPTION: Function MINT is the modified integer function.
!\\
!\\
! !INTERFACE:
!
      function mint( x ) result ( value )
!
! !INPUT PARAMETERS: 
!
      real*8, intent(in) :: x
!
! !RETURN VALUE:
!
      real*8             :: value
!
! !REMARKS:
!  The modified integer function is defined as follows:
!
!            { -INT( ABS( X ) )   for X < 0
!     MINT = {
!            {  INT( ABS( X ) )   for X >= 0
!
! !REVISION HISTORY: 
!  20 Nov 2001 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      if ( x < 0d0 ) then 
         value = -int( abs( x ) )        
      else
         value =  int( abs( x ) )        
      endif
   
      end function mint
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: caldate
!
! !DESCRIPTION: Subroutine CALDATE converts an astronomical Julian day to 
!  the YYYYMMDD and HHMMSS format.
!\\
!\\
! !INTERFACE:
!
      subroutine caldate( julianday, yyyymmdd, hhmmss )
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      real*8,  intent(in)  :: julianday  ! astronomical julian date 
!
! !OUTPUT PARAMETERS: 
!
      integer, intent(out) :: yyyymmdd   ! date in yyyy/mm/dd format
      integer, intent(out) :: hhmmss     ! time in hh:mm:ss format
!
! !REMARKS:
!   Algorithm taken from "Practical Astronomy With Your Calculator",
!   Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
!
! !REVISION HISTORY: 
!  (1 ) Now compute HHMMSS correctly.  Also use REAL*4 variables HH, MM, SS
!        to avoid roundoff errors. (bmy, 11/21/01)
!  (2 ) Renamed NYMD to YYYYMMDD and NHMS to HHMMSS for documentation
!        purposes (bmy, 6/26/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
      real*4               :: hh, mm, ss
      real*8               :: a, b, c, d, day, e, f 
      real*8               :: fday, g, i, j, jd, m, y

      !=================================================================
      ! CALDATE begins here!
      ! See "Practical astronomy with your calculator", Peter Duffett-
      ! Smith 1992, for an explanation of the following algorithm.
      !=================================================================
      jd = julianday + 0.5d0
      i  = int( jd )
      f  = jd - int( i )
    
      if ( i > 2299160d0 ) then
         a = int( ( i - 1867216.25d0 ) / 36524.25d0 )
         b = i + 1 + a - int( a / 4 )
      else
         b = i
      endif
    
      c = b + 1524d0
      
      d = int( ( c - 122.1d0 ) / 365.25d0 )
    
      e = int( 365.25d0 * d )
    
      g = int( ( c - e ) / 30.6001d0 )
   
      ! DAY is the day number
      day  = c - e + f - int( 30.6001d0 * g ) 
    
      ! FDAY is the fractional day number
      fday = day - int( day )
      
      ! M is the month number
      if ( g < 13.5d0 ) then
         m = g - 1d0
      else
         m = g - 13d0
      endif
    
      ! Y is the year number
      if ( m > 2.5d0 ) then
         y = d - 4716d0
      else
         y = d - 4715d0
      endif
    
      ! Year-month-day value
      yyyymmdd = ( int( y ) * 10000 ) + ( int( m ) * 100 ) + int( day )
      
      ! Hour-minute-second value
      ! NOTE: HH, MM, SS are REAL*4 to avoid numerical roundoff errors
      hh     = fday * 24d0 
      mm     = ( hh - int( hh ) ) * 60d0
      ss     = ( mm - int( mm ) ) * 60d0   
      hhmmss = ( int( hh ) * 10000 ) + ( int( mm ) * 100 ) + int( ss )
      
      end subroutine caldate
!EOC
      end module julday_mod
