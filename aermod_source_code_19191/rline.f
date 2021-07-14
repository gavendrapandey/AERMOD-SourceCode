      SUBROUTINE RLCALC
C***********************************************************************
C                 RLCALC Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Introduce RLINE source. Based on RLINEv1.2 released 
C                     November 2013. 
C
C        PROGRAMMER:  M. Snyder, R. Cleary, J. Paumier, Wood
C
C        DATE:        July 20, 2018
C
C        INPUTS:      None
C
C        OUTPUTS:    
C
C          MODIFIED:  Added Urban option for RLINE sources. 
C                     Wood, 03/04/2019
C
C        CALLED FROM: MAIN
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: IREC, AZFLAG, HRVAL, NUMTYP,
     &                 NUMGRP, IGRP, ISRC, IREC, NUMREC, L_HRLYSIG,
     &                 ITYP, ARM2, CHI, EVONLY, AQS, AHS, ASZINI, QFLAG,
     &                 URBSRC, NUMURB, STABLE, L_MorningTrans,
     &                 ZI, ZIURB, ZIRUR, ZIMECH, URBAN, SRCTYP,
     &                 USTAR, URBUSTR, RURUSTR, URBSTAB,
     &                 IURB, OBULEN, URBOBULEN, RUROBULEN, IURBGRP,
     &                 WSTAR, URBWSTR, EMIFAC
      USE RLINE_DATA
      IMPLICIT NONE
      INTEGER  :: I
      DOUBLE PRECISION  :: ERROR
      DOUBLE PRECISION  :: CONCD
      DOUBLE PRECISION  :: ZB
      DOUBLE PRECISION  :: CONCENTRATION
C     ERROR         = error in numerical integration
C     CONCD         = dummy concentration at receptor
C     ZB            = release height for barriers
C     CONCENTRATION = concentration
      
C     Variable Initializations
      CONCD = 0.0D0
      CONCENTRATION = 0.0D0
     
C     Initialize __VAL arrays
      HRVAL = 0.0D0

C     Allow for HOUREMIS input to override the input QEMIS on the 
C     SO SRCPARAM card. 
      IF(QFLAG(ISRC) .EQ. 'HOURLY') THEN
C        Set hourly variable source height and sigmaz-initial (optional)
         IF(L_HRLYSIG(ISRC)) THEN
           RLSOURCE(ISRC)%ZSB = AHS(ISRC)
           RLSOURCE(ISRC)%ZSE = AHS(ISRC)
           RLSOURCE(ISRC)%INIT_SIGMAZ = ASZINI(ISRC) 
         END IF
C        Set hourly variable emission rate (required)	  
         IF(SRCTYP(ISRC) .EQ. 'RLINEXT') THEN
           RLSOURCE(ISRC)%QEMIS = AQS(ISRC)
         ELSE !RLINE source with Lnemis inputs
           RLSOURCE(ISRC)%QEMIS = AQS(ISRC)*RLSOURCE(ISRC)%WIDTH
         END IF
      END IF
C     Perform first hour calculations
      IF(RLFIRSTHR) THEN
C        Create exponential tables                                           --- CALL CREATE_EXP_TABLE
         CALL CREATE_EXP_TABLE
C        Perform MOVES to RLINE unit conversion                              --- CALL RLEMCONV
         CALL RLEMCONV
         RLFIRSTHR = .FALSE.
      END IF

C        Set Mixing Height and Adjust L & USTAR for Urban Option if Needed
         IF (URBSRC(ISRC) .EQ. 'Y') THEN
C           Find Urban Area Index for This Source
            DO I = 1, NUMURB
               IF (IURBGRP(ISRC,I) .EQ. 1) THEN
                  IURB = I
                  EXIT
               END IF
            END DO
            IF (STABLE .OR. L_MorningTrans(IURB)) THEN
               URBSTAB = .TRUE.
               ZI = MAX( ZIURB(IURB), ZIMECH )
               OBULEN = DABS( URBOBULEN(IURB) )
               USTAR  = URBUSTR(IURB)
               RLWSTAR = URBWSTR(IURB)
            ELSE
               URBSTAB = .FALSE.
               ZI = ZIRUR
               OBULEN = RUROBULEN
               USTAR  = RURUSTR
               RLWSTAR = WSTAR
            END IF
         ELSE IF (URBAN .AND. URBSRC(ISRC) .EQ. 'N') THEN
            URBSTAB = .FALSE.
            ZI = ZIRUR
            OBULEN = RUROBULEN
            USTAR  = RURUSTR
            RLWSTAR = WSTAR
         ELSE
C ---       Rural
            URBSTAB = .FALSE.
            RLWSTAR = WSTAR
         END IF
  
C        Calculate SIGMAV: Compute SIGMAV, UEFF, and THETAW                  --- CALL COMPUTE_MET
         CALL COMPUTE_MET 

      IF(.NOT. RLPROCESSED) THEN
C        Translate and rotate the line source to align with wind direction   --- CALL TRANSLTE_ROTATE
         CALL TRANSLATE_ROTATE
      END IF

C     Begin loop over receptors
      RECEPTOR_LOOP: DO IREC = 1, NUMREC
C        Rotate X, Y receptors.  Z receptor not rotated.
         XR_ROT = XRCP_ROT(IREC)
         YR_ROT = YRCP_ROT(IREC)
         ZRECEP = AZFLAG(IREC)

C        Get translated an rotated source information
         SIGMAZ0 = RLSOURCE(ISRC)%INIT_SIGMAZ 
         XSBEGIN = XSB_ROT(ISRC)
         YSBEGIN = YSB_ROT(ISRC)
         XSEND   = XSE_ROT(ISRC)  
         YSEND   = YSE_ROT(ISRC)

C        Calculate a vertical displacement of the source to  
C        account for the effect of a roadside barrier                        --- CALL BARRIER_DISPLACEMENT
         IF(RLSOURCE(ISRC)%HTWALL .GT. 0.0) THEN
            CALL BARRIER_DISPLACEMENT(ZB)
            ZSBEGIN = ZB
            ZSEND   = ZB
         ELSE 
            ZSBEGIN = RLSOURCE(ISRC)%ZSB
            ZSEND = RLSOURCE(ISRC)%ZSE
         END IF

C        Calculate the contribution to concentration at a 
C        receptor due to a line source using Romberg integration             --- CALL NUMERICAL_LINE_SOURCE
         SIGMAY0 = 0.0D0
         SIGMAY0 = DABS(0.5D0 * (RLSOURCE(ISRC)%WIDTH)
     &      * DSIN(THETAW + SM_NUM))
         CALL NUMERICAL_LINE_SOURCE(CONCD, ERROR)

C        Convert Qemis from MOVES units (g/hr/link) to RLINE units (g/m/s)
C        EMIFAC(1) is for concentration
         CONCENTRATION = CONCD * EMIFAC(1) * RLSOURCE(ISRC)%QEMIS *
     &      RLEMISCONV(ISRC) 
         HRVAL(:) = CONCENTRATION
                     
C        For the initial integration of R-LINE v1.2 into AERMOD, 
C        only including ARM2 chemistry options.  OLM and PVMRM not included.
         IF (ARM2) THEN
C           Store conc by source and receptor for ARM2 options
            DO ITYP = 1, NUMTYP
               CHI(IREC,ISRC,ITYP) = HRVAL(ITYP)
            END DO

C           Initialize __VAL arrays (1:NUMTYP)
            HRVAL   = 0.0D0

         END IF
                
C        Sum HRVAL to AVEVAL and ANNVAL Arrays                               --- CALL SUMVAL
         IF (EVONLY) THEN
            CALL EV_SUMVAL
         ELSE
            DO IGRP = 1, NUMGRP
               CALL SUMVAL
            END DO
         END IF

C        Initialize __VAL arrays (1:NUMTYP)
         HRVAL(:) = 0.0D0

C        Reset concentration value
         CONCENTRATION = 0.0D0

      END DO RECEPTOR_LOOP

      END SUBROUTINE RLCALC

      SUBROUTINE BARRIER_DISPLACEMENT(HEFF)
C***********************************************************************
C       BARRIER_DISPLACEMENT Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate a vertical displacement of the source to
C                     account for the effect of a roadside barrier. 
C
C        PROGRAMMER:  M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      RLINE source, receptor and barrier parameters
C
C        OUTPUTS:    
C
C        CALLED FROM: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: ISRC
      USE RLINE_DATA, ONLY: RLSOURCE, XSBEGIN, XSEND, YSBEGIN, YSEND, 
     &                      YR_ROT, XR_ROT, DISPHT 
      IMPLICIT NONE

      DOUBLE PRECISION  :: HEFF
      DOUBLE PRECISION  :: XD, THETA_LINE, SBDIST
      DOUBLE PRECISION  :: HTWALL, FUP, DW, M, B, ZAVE
C     HEFF        = effective height
C     XD          = perpendicular distance between source and receptor
C     THETA_LINE  = angle between wind direction and source
C     SBDIST      = barrier downwind distance from line
C     HTWALL      = barrier height
C     FUP         = barrier height factor defining vertical wake zone behind barrier
C     DW          = barrier distance from source
C     M           = y-intercept of extended source
C     B           = slope of source
C     ZAVE        = average height
      
      FUP        = 1.5D0
      HTWALL     = RLSOURCE(ISRC)%HTWALL
      DW         = (RLSOURCE(ISRC)%DCLWALL - RLSOURCE(ISRC)%DCL);
      ZAVE       = (RLSOURCE(ISRC)%ZSB + RLSOURCE(ISRC)%ZSE) / 2.0D0;
      THETA_LINE = DATAN2(YSEND - YSBEGIN, XSEND - XSBEGIN)
      SBDIST     = DW*DSIN(THETA_LINE) 

      IF(SBDIST .GT. 0.0) THEN
         IF((XSEND-XSBEGIN) .EQ. 0.0) THEN
            XD = DABS(XR_ROT - XSBEGIN)
         ELSE
            M  = (YSEND - YSBEGIN) / (XSEND - XSBEGIN)
            B  = (YSBEGIN * XSEND - YSEND * XSBEGIN) / (XSEND - XSBEGIN)
            XD = DABS(YR_ROT-M * XR_ROT - B) / DSQRT(M * M + 1.0D0)
         END IF
         HEFF = MAX(ZAVE, FUP * HTWALL - DABS(FUP * HTWALL - DISPHT) /
     &          (2.0D0 * DW)* XD)
      ELSE
         HEFF = ZAVE
      END IF
 
      END SUBROUTINE BARRIER_DISPLACEMENT

      SUBROUTINE COMPUTE_MET
C***********************************************************************
C        COMPUTE_MET Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate SIGMAV using USTAR and WSTAR. 
C
C        PROGRAMMER:  A. Venkatram, M. Snyder
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      Meteorological variables
C
C        OUTPUTS:     SIGMAV
C
C        CALLED FROM: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: USTAR, WDREF, SFCZ0
      USE RLINE_DATA, ONLY: SIGMAV, THETAW, PI, FAC_DISPHT, DISPHT, 
     &                 RLWSTAR
      IMPLICIT NONE

      DOUBLE PRECISION  :: ANGLE, SIGMAV_CALC, WSTAR_LOC
C     ANGLE       = 270.0 - wind direction
C     SIGMAV_CALC = SIGMAV calculated from WSTAR and USTAR
C     WSTAR_LOC   = check for WSTAR variable


C     Check for WSTAR from metext.f
      WSTAR_LOC = MAX(RLWSTAR, 0.0D0)
    
C     Calculate SIGMAV from WSTAR and USTAR variables
C     Calculate standard deviation of wind direction
      SIGMAV_CALC = DSQRT((0.6D0 * WSTAR_LOC)**2.0D0 + 
     &              (1.9D0 * USTAR)**2.0D0)
      SIGMAV      = MAX(SIGMAV_CALC, 0.2D0) 

      ANGLE = 270.0D0 - WDREF
      IF (ANGLE > 180.0D0) THEN            
         ANGLE = ANGLE - 360.0D0
      END IF
      THETAW = ANGLE*PI / 180.0D0
      
C     Calculate displacement height 
      DISPHT = FAC_DISPHT * SFCZ0

      END SUBROUTINE COMPUTE_MET

      SUBROUTINE CREATE_EXP_TABLE
C***********************************************************************
C        CREATE_EXP_TABLE Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Create a lookup table based on arguments of the
C                     built-in exponential function to improve 
C                     computation time.
C
C        PROGRAMMER:  A. Venkatram
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      
C
C        OUTPUTS:     
C
C        CALLED FROM: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE RLINE_DATA, ONLY: XEXP, AEXP, BEXP, DELEXP
      IMPLICIT NONE

      INTEGER  :: IND
      DOUBLE PRECISION, DIMENSION(1000)  :: EXT
C     IND         = local source index
C     EXT         = exponent

    
      DELEXP  = 20.0D0 / 999.0D0
      XEXP(1) = -20.0D0
      EXT(1)  = DEXP(XEXP(1))

      DO IND=2, 1000
         XEXP(IND) = XEXP(IND-1) + DELEXP
         EXT(IND)  = DEXP(XEXP(IND))
      END DO

      DO IND=1,999
         BEXP(IND) = (EXT(IND+1) - EXT(IND)) / (XEXP(IND+1) - XEXP(IND))
         AEXP(IND) = EXT(IND) - BEXP(IND) * XEXP(IND)
      END DO

      BEXP(1000) = BEXP(999)
      AEXP(1000) = AEXP(999)

      END SUBROUTINE CREATE_EXP_TABLE

      DOUBLE PRECISION FUNCTION DEPRESSED_DISPLACEMENT(THETA_LINE,IND)
C***********************************************************************
C        DEPRESSED_DISPLACEMENT Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Computes transformation for a depressed roadway,
C                     shifting the roadway upwind, and compresses or expands 
C                     the roadway based on the fractional width and distance
C                     from the centerline. 
C
C        PROGRAMMER:  M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:  THETA_LINE, IND
C
C         RETURNED:     
C
C        CALLING ROUTINES: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE RLINE_DATA, ONLY: RLSOURCE, THETAW
      IMPLICIT NONE

      INTEGER  :: IND
      DOUBLE PRECISION  :: THETA_LINE
      DOUBLE PRECISION  :: DEPTH, WTOP, WBOTTOM, DCL 
      DOUBLE PRECISION  :: EFFD, RELD, FRACW, EFFW, THETA_REL, DCRIT, F
C     IND         = local source index
C     THETA_LINE  = angle between wind direction and source
C     DEPTH       = depth of depression
C     WTOP        = width of top of depression
C     WBOTTOM     = width of bottom of depression
C     DCL         = offset distance from center line
C     EFFD        = effective depth
C     RELD        = relative roadway depth
C     FRACW       = fractional width
C     EFFW        = effective width
C     THETA_REL   = relative angle between roadway and wind direction
C     DCRIT       = critical depth
C     F           = effective wind fraction

      DEPTH     = RLSOURCE(IND)%DEPTH
      WBOTTOM   = RLSOURCE(IND)%WBOTTOM
      WTOP      = RLSOURCE(IND)%WTOP     
      DCL       = RLSOURCE(IND)%DCL

      THETA_REL = THETA_LINE - THETAW 

      EFFD      = (WBOTTOM * DABS(DEPTH) + ((WTOP - WBOTTOM) /
     &            2D0 * DABS(DEPTH))) / WTOP
      RELD      = EFFD / WBOTTOM
 
      IF (RELD .GE. 0.0483) THEN
         FRACW = 1.0D0
      ELSE
         FRACW = -0.0776D0 + DSQRT(1.506D0 - 7.143D0 * RELD)
      END IF
 
      EFFW    = FRACW**(1.0D0 - (DCOS(DABS(THETA_REL)) * 
     &          DCOS(DABS(THETA_REL)))) * WBOTTOM
      DCRIT   = 0.2064D0 * WTOP * WBOTTOM / (0.5D0 * (WTOP + WBOTTOM))
      F       = MIN(1.0D0, WBOTTOM / WTOP * 
     &          (1.0D0 + DABS(DEPTH) / DCRIT))
      
      DEPRESSED_DISPLACEMENT = ((WTOP * F - EFFW) / 2.0D0)*
     &                         ((DSIN(THETA_REL))**2) * 
     &                         DSIGN(1.0D0, DSIN(THETA_REL)) -
     &                         (EFFW / WBOTTOM * DCL) * 
     &                         DSIGN(1.0D0, SIN(THETA_LINE))

      END FUNCTION

      SUBROUTINE EFFECTIVE_WIND(XD,HEFF)
C***********************************************************************
C        EFFECTIVE_WIND Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate the effective wind speed at mean
C                     plume height.

C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      Meteorological variables
C
C        OUTPUTS:     
C
C        CALLED FROM: MEANDER
C                     POINT_CONC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: UREFHT, USTAR, OBULEN, UREF, SFCZ0
      USE RLINE_DATA, ONLY: UEFF, DISPHT, SIGMAV
      IMPLICIT NONE

      INTEGER  :: ITER
      DOUBLE PRECISION  :: ERF
      DOUBLE PRECISION  :: SZ, SZ_NEW, ERR, UREFCALC, ZBAR
      DOUBLE PRECISION, INTENT(IN)  :: XD, HEFF
C     ITER        = iteration
C     ERF         = error function
C     SZ          = effective SIGMAZ
C     SZ_NEW      = intermediate vertical dispersion
C     ERR         = error in each successive calculation
C     UREFCALC    = calculated wind speed at reference height
C     ZBAR        = mean plume height
C     XD          = perpendicular distance between source and receptor
C     HEFF        = effective height

C     Computation Parameters
      DOUBLE PRECISION, PARAMETER  :: SQ2PI = 0.797885D0

C     External Functions:
      DOUBLE PRECISION, EXTERNAL  :: SIGMAZ, MOST_WIND

C     Computation parameters
      SZ   = SIGMAZ(XD)
      ERR  = 10
      ITER = 1

C     MOST_WIND function      
      UREFCALC = MOST_WIND(UREFHT, SFCZ0, DISPHT, USTAR, OBULEN)
      UEFF     = MOST_WIND(HEFF, SFCZ0, DISPHT, USTAR, OBULEN)
     & * UREF / UREFCALC
      
      DO WHILE ((ERR > 1.0E-02) .AND. (ITER < 20)) 
         ZBAR   = SQ2PI * SZ * DEXP(-0.5 * (HEFF / SZ)**2.0D0) +
     &            HEFF * ERF(HEFF / (DSQRT(2.0D0) * SZ))
C        MOST_WIND function 
         UEFF   = MOST_WIND(MAX(ZBAR, HEFF), SFCZ0, DISPHT, USTAR,
     &            OBULEN) * UREF / UREFCALC
         UEFF   = DSQRT(2.0D0 * SIGMAV**2.0D0  + UEFF**2.0D0)
         SZ_NEW = SIGMAZ(XD)
         ERR    = DABS((SZ_NEW - SZ) / SZ)
         SZ     = SZ_NEW
         ITER   = ITER + 1
      END DO
  
      END SUBROUTINE EFFECTIVE_WIND

      DOUBLE PRECISION FUNCTION EXPX(XP)
C***********************************************************************
C        EXPX Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Computes the exponential function using a table. 
C
C        PROGRAMMER:  A. Venkatram
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:  XP
C
C         RETURNED:     
C
C        CALLING ROUTINES: MEANDER
C                          POINT_CONC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE RLINE_DATA, ONLY: DELEXP, AEXP, BEXP, XEXP
      IMPLICIT NONE

      INTEGER  :: P
      DOUBLE PRECISION  :: XPD
      DOUBLE PRECISION, INTENT(IN) :: XP 
C     P           = exponential table index
C     XPD         = closest precalculated exponent value 
C     XP          = input exponent value


      XPD  = XP
      XPD  = MAX(XEXP(1), XPD)
      P    = FLOOR((XPD + 20.0D0) / DELEXP) + 1.0D0
      EXPX = AEXP(P) + BEXP(P) * XPD
      
      END FUNCTION

      DOUBLE PRECISION FUNCTION MEANDER(X, Y, Z)
C***********************************************************************
C        MEANDER Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate the contribution of a point source 
C                     at (X,Y,Z) to a receptor at (Xr_rot,Yr_rot,Zrcp),
C                     assuming that the material spreads out radially
C                     in all directions. 
C
C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:   X, Y, Z
C
C         RETURNED:     
C
C        CALLING ROUTINES: POINT_CONC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE RLINE_DATA, ONLY: XR_ROT, YR_ROT, XD_MIN, ZRECEP, UEFF,
     &                      PI, XD_MIN
      IMPLICIT NONE

      DOUBLE PRECISION  :: R, VERT, HORZ, SZ, HEFF
      DOUBLE PRECISION, INTENT(IN)  :: X, Y, Z
C     R           = radial distance to receptor
C     VERT        = vertical component of concentration
C     HORZ        = horizontal component of concentration
C     SZ          = effective SIGMAZ
C     HEFF        = effective height

C     External functions:
      DOUBLE PRECISION :: SIGMAZ, EXPX


      R       = DSQRT((XR_ROT - X)**2.0D0 + (YR_ROT - Y)**2.0D0)   
      R       = MAX(R, XD_MIN)
      HEFF    = Z

C     Calculate effective wind speed                                         --- CALL EFFECTIVE_WIND
      CALL EFFECTIVE_WIND(R,HEFF)

C     Account for source height
      SZ      = SIGMAZ(R)
      VERT    = DSQRT(2.0D0 / PI) * (EXPX(-0.5D0 * ((HEFF - ZRECEP)
     &          / SZ)**2.0D0) + EXPX(-0.5D0 * 
     &          ((HEFF + ZRECEP) / SZ)**2.0D0))/(2.0D0*SZ*UEFF)

      HORZ    = 1.0D0 / (2.0D0 * PI * R)
      MEANDER = VERT * HORZ

      END FUNCTION

      DOUBLE PRECISION FUNCTION MOST_WIND(Z, Z0, DH, UST, L)
C***********************************************************************
C        MOST_WIND Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate the wind speed from similarity theory. 
C
C        PROGRAMMER:  A. Venkatram
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:   Z, Z0, DH, UST, L
C
C         RETURNED:     
C
C        CALLING ROUTINES: EFFECTIVE_WIND
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: STABLE
      USE RLINE_DATA, ONLY: PI, SIGMAV
      IMPLICIT NONE

      DOUBLE PRECISION  :: X1, X2, PSI1, PSI2, L_LOC
      DOUBLE PRECISION, INTENT(IN)  :: Z, Z0, DH, UST, L
C     X1          = computation within PSI1
C     X2          = computation within PSI2
C     PSI1        = stability function
C     PSI2        = stability function
C     L_LOC       = check for L variable
C     Z           = height
C     Z0          = surface roughness length
C     DH          = displacement height
C     UST         = surface friction velocity
C     L           = Monin-Obukhov length

C     Computation Parameters
      DOUBLE PRECISION, PARAMETER :: KAPPA=0.4D0


      IF (DH .GE. Z) THEN
         MOST_WIND = DSQRT(2.0D0) * SIGMAV
      ELSE
         IF (STABLE) THEN
            L_LOC = DABS(L)
            PSI1 = -17.0D0 * (1.0D0 - DEXP(-0.29D0 * (Z - DH) / L_LOC))  
            PSI2 = -17.0D0 * (1.0D0 - DEXP(-0.29D0 * Z0 / L_LOC))
         ELSE
            L_LOC = -1.0D0*DABS(L)
            X1   = (1.0D0 - 16.0D0 * (Z - DH) / L_LOC)**0.25D0
            X2   = (1.0D0 - 16.0D0 * Z0 / L_LOC)**0.25D0     
            PSI1 = 2.0D0 * DLOG((1.0D0 + X1) / 2.0D0) + 
     &             DLOG((1.0 + X1 * X1) / 2.0D0) -
     &             2.0D0 * DATAN(X1) + PI / 2.0D0
            PSI2 = 2.0D0 * DLOG((1.0D0 + X2) / 2.0D0) +
     &             DLOG((1.0 + X2 * X2) / 2.0D0) -
     &             2.0D0 * DATAN(X2) + PI / 2.0D0
         END IF
         MOST_WIND = UST*(DLOG((Z - DH) / Z0) - PSI1 + PSI2) / KAPPA
      END IF

      END FUNCTION

      SUBROUTINE NUMERICAL_LINE_SOURCE(CONC_NUM,ERR)
C***********************************************************************
C     NUMERICAL_LINE_SOURCE Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate the contribution to concentration 
C                     at a receptor due to a line source using Romberg 
C                     integration wind speed from similarity theory. 
C
C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      
C
C        OUTPUTS:     
C
C        CALLED FROM: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE RLINE_DATA, ONLY: XSBEGIN, YSBEGIN, ZSBEGIN, XSEND, YSEND,
     &                      ZSEND, XR_ROT, YR_ROT, SM_NUM,
     &                      ERROR_LIMIT, XD_MIN
      IMPLICIT NONE

      INTEGER  :: NO_POINTS, J, IS, MINJ, IT_LIM
      INTEGER  :: ALLOCATESTATUS, ALLOCERROR
      INTEGER  :: ST, FI
C     NO_POINTS       = number of points added each step 
C     J               = index
C     IS              = index for integration of point sources
C     MINJ            = minimum number of iterations
C     IT_LIM          = maximum number of iterations
C     ALLOCATESTATUS  = flag for intermediate array allocation error
C     ALLOCERROR      = status flag for error of allocation
C     ST              = starting indice
C     FI              = finishing indice
      
      DOUBLE PRECISION  :: XDIF, YDIF, ZDIF
      DOUBLE PRECISION  :: DISP
      DOUBLE PRECISION  :: HDUM(3), CONCDUM(3) 
      DOUBLE PRECISION  :: CONC_INT
      DOUBLE PRECISION  :: X, Y, Z, DELT, THETA, PHI
      DOUBLE PRECISION  :: T, TMAX, COST, SINT, SINP, COSP, XREC
      DOUBLE PRECISION  :: DR
      DOUBLE PRECISION  :: XTEMP, YTEMP, ZTEMP
C     XDIF           = x integral limit 
C     YDIF           = y integral limit
C     ZDIF           = z integral limit
C     DISP           = dummy variable used to store integral
C     HDUM           = step size
C     CONCDUM        = successive concentration approximations
C     CONC_INT       = numerical integrals
C     X              = x-coordinate of receptor
C     Y              = y-coordinate of receptor
C     Z              = z-coordinate of receptor
C     DELT           = distance between points used to estimate line source
C     THETA          = horizontal wind direction
C     PHI            = angle of line elevation
C     T              = half of DELT
C     TMAX           = 3-dimensional 
C     COST           = cosine of theta
C     SINT           = sin of theta
C     SINP           = sin of phi
C     COSP           = cosine of phi
C     XREC           = x-coordinate of point on source directly upwind of receptor
C     DR             = distance beetween source and receptor in wind direction
C     XTEMP          = temporary x-coordinate of receptor for line orientation
C     YTEMP          = temporary y-coordinate of receptor for line orientation
C     ZTEMP          = temporary z-coordinate of receptor for line orientation
      
      DOUBLE PRECISION, INTENT(OUT)  :: CONC_NUM
      DOUBLE PRECISION, INTENT(OUT)  :: ERR
C     CONC_NUM       = numerical routine output concentration 
C     ERR            = integration is set to an arbitrarily large value before it is reduced
      
      DOUBLE PRECISION, ALLOCATABLE  :: H(:), CONC(:)
C     H              =  step size
C     CONC           =  successive concentration approximations

C     2K is the order of Romberg integration.  Nmax needs to be greater than K for Romberg integration to be used, 
C     otherwise trapezoidal integration is used
      INTEGER, PARAMETER  :: K=3 

C     Computation Parameters
      DOUBLE PRECISION, PARAMETER :: XINTERP = 0.0D0

C     External unctions:
      DOUBLE PRECISION  :: POINT_CONC

C     Variable Initializations
      CONC_NUM = 0.0

C     Orient the end points so the begining has a lower Y value
      IF (YSEND < YSBEGIN) THEN
         XTEMP   = XSEND
         YTEMP   = YSEND
         ZTEMP   = ZSEND
         XSEND   = XSBEGIN
         YSEND   = YSBEGIN
         ZSEND   = ZSBEGIN
         XSBEGIN = XTEMP
         YSBEGIN = YTEMP
         ZSBEGIN = ZTEMP 
      END IF

      XDIF  = XSEND - XSBEGIN  
      YDIF  = YSEND - YSBEGIN
      ZDIF  = ZSEND - ZSBEGIN
      THETA = DATAN2(YDIF, XDIF)
      TMAX  = DSQRT(XDIF * XDIF + YDIF * YDIF + ZDIF * ZDIF)    
      PHI   = DASIN(ZDIF / (TMAX + SM_NUM))  
      SINP  = DSIN(PHI)
      COSP  = DCOS(PHI)
      COST  = DCOS(THETA)
      SINT  = DSIN(THETA)

C     Find x-coordinate of point on source line directly upwind of receptor
      XREC = (XSEND - (YSEND - YR_ROT) * (XSEND - XSBEGIN) / 
     &       (YSEND -YSBEGIN))
      DR   = DABS(XR_ROT - XREC)
C     Prevent user from placing receptor on source
      DR   = MAX(XD_MIN, DR)

C     Convergence Criteria: Minimum Iterations    
      IF ((YR_ROT > YSBEGIN - TMAX / 2.0D0 * DABS(COST)) .AND. 
     &   (YR_ROT < YSEND + TMAX / 2.0D0 * DABS(COST))) THEN
         MINJ = CEILING(DLOG(2.0D0 * TMAX / 
     &          (MAX(XD_MIN, DR * DABS(SINT)))
     &          - 2.0D0) / DLOG(2.0D0)) + 2.0D0
      ELSE
C        Set MINJ = 0 so if receptor is upwind, the conc will converge quickly  
         MINJ = 0
      END IF

C     If receptor is upwind of the line 
      IF ((XR_ROT < XSBEGIN) .AND. (XR_ROT < XSEND)) MINJ = 0
         IT_LIM = MAX(10, 2 * MINJ)

       ALLOCATE(H(IT_LIM), STAT = ALLOCATESTATUS)
       ALLOCATE(CONC(IT_LIM), STAT = ALLOCATESTATUS)

C     Compute concentration at receptor
C     Initialize concentrations         
      CONC(:)    = 0.0D0
      H(:)    = 0.0D0
      
      DISP    = (POINT_CONC(XSBEGIN, YSBEGIN, ZSBEGIN) + 
     &          POINT_CONC(XSEND, YSEND, ZSEND)) * 0.5D0
C     Calculate first approximation of the integration.  Set relative size of
C     integration interval
      CONC(1) = DISP * TMAX
      H(1)    = 1.0D0
C     Trapezoidal integration 
      DO J = 2, IT_LIM 
         NO_POINTS = 2.0D0**(J - 2)    
         DELT      = TMAX / NO_POINTS    
         T         = DELT / 2.0D0
         DISP      = 0.0D0
         DO IS = 1, NO_POINTS
            X      = T * COST * COSP + XSBEGIN
            Y      = T * SINT * COSP + YSBEGIN
            Z      = T * SINP + ZSBEGIN
            DISP   = DISP + POINT_CONC(X, Y, Z)
            T      = T + DELT
         END DO
         CONC(J) = (DISP * DELT + CONC(J - 1)) / 2.0D0
C        See page 134 in "Numerical Receipes" for an explanation
         H(J)    = 0.25D0 * H(J - 1)

C        Romberg integration is invoked if (J >= K)
         IF (J >= K) THEN
            ST       = J - K + 1
            FI       = ST + K - 1
            HDUM     = H(ST:FI)
            CONCDUM  = CONC(ST:FI)
C           Extrapolate to H=0.0 to compute integral                         --- CALL POLYINTERP
            CALL POLYINTERP(CONC_INT, ERR, HDUM, CONCDUM, XINTERP, K)
            CONC_NUM = DABS(Conc_Int)

C           Check convergence criteria          
            IF ((DABS(ERR) < ERROR_LIMIT) .AND. (J > MINJ)) THEN
               DEALLOCATE(H, CONC, STAT=ALLOCERROR)
               RETURN
            END IF

         END IF
         
        CONC_NUM = DABS(CONC(J))
      
      END DO
      
      END SUBROUTINE NUMERICAL_LINE_SOURCE

      DOUBLE PRECISION FUNCTION POINT_CONC(X, Y, Z)
C***********************************************************************
C        POINT_CONC Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate the direct plume contribution 
C                     of a point source using Gaussian dispersion and
C                     combine the direct plume and meander contributions 
C                     to determine the total concentration from a point
C                     to the receptor. 
C
C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:   X, Y, Z
C
C         RETURNED:     
C
C        CALLING ROUTINES: NUMERICAL_LINE_SOURCE
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE RLINE_DATA, ONLY: XR_ROT, YR_ROT, ZRECEP, UEFF, PI,
     &                      XD_MIN, SIGMAV
      IMPLICIT NONE

      DOUBLE PRECISION  :: CONC_M, CONC_P
      DOUBLE PRECISION  :: SY, SZ
      DOUBLE PRECISION  :: FRAN
      DOUBLE PRECISION  :: VERT, HORZ, XD, YD 
      DOUBLE PRECISION  :: HEFF
      DOUBLE PRECISION, INTENT(IN)  :: X, Y, Z
C     CONC_M      = meander concentration
C     CONC_P      = direct plume concentration
C     SY          = horizontal dispersion
C     SZ          = vertical dispersion
C     FRAN        = fraction between meander and wind following plume
C     VERT        = vertical component of concentration
C     HORZ        = horizontal component of concentration
C     XD          = distance between source and receptor in direction parallel to wind
C     YD          = distance between source and receptor in direction perpendicular to wind
C     HEFF        = effective height
C     X           = x-coordinate of receptor
C     Y           = y-coordinate of receptor
C     Z           = z-coordinate of receptor

C     External functions:
      DOUBLE PRECISION :: MEANDER, EXPX
      DOUBLE PRECISION :: SIGMAY, SIGMAZ

C     Variable Initializations
      CONC_P = 0.0D0
      CONC_M = 0.0D0
      FRAN   = 0.0D0
      SY     = 0.0D0
      XD     = XR_ROT - X
      HEFF   = Z  

C     Calculate effective wind speed                                         --- CALL EFFECTIVE_WIND     
      CALL EFFECTIVE_WIND(XD, HEFF)
      
C     Calculate direct plume concentration
      IF (XD < 0.0001D0) THEN     
         CONC_P = 0.0D0   
      ELSE
         XD     = MAX(XD, XD_MIN)
         YD     = YR_ROT - Y
         SZ     = SIGMAZ(XD) 
         SY     = SIGMAY(XD)
         VERT   = DSQRT(2.0D0 / PI) *
     &            (EXPX(-0.5D0 * ((HEFF - ZRECEP) / SZ)**2.0D0) +
     &            EXPX(-0.5D0 * ((HEFF + ZRECEP) / SZ)**2.0D0)) /
     &            (2.0D0 * SZ * UEFF)
         HORZ   = 1.0D0 / (DSQRT(2.0D0 * PI) * SY) *
     &            EXPX(-0.5D0*(YD/SY)**2.0D0)
         CONC_P = VERT * HORZ 
      END IF
      
C     Calculate meander concentration
      FRAN   = 2.0D0 * SIGMAV * SIGMAV / (UEFF * UEFF)
      CONC_M = MEANDER(X, Y, Z)

C     Combine direct plume and meander contributions
      POINT_CONC = CONC_P * (1.0D0 - FRAN) + CONC_M * FRAN

      END FUNCTION

      SUBROUTINE POLYINTERP(Y,ERR,XA,YA,X,N)
C***********************************************************************
C        POLYINTERP Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Given vectors xa and ya, y is interpolated value 
C                     for x; err is the error in interpolation.   
C                     Uses polynomial interpolation from "Numerical 
C                     Recipes in Fortran" pages 103-104. 
C
C        PROGRAMMER:  A. Venkatram
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      
C
C        OUTPUTS:     
C
C        CALLED FROM: NUMERICAL_LINE_SOURCE
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C        
C                    "Numerical Recipes in Fortran", Press et al., (1992)
C***********************************************************************

C     Variable Declarations
      IMPLICIT NONE

      INTEGER  :: NS
      INTEGER  :: I, M
C     NS          = position of x in array 
C     I           = counting index
C     M           = counting index

      DOUBLE PRECISION  :: DIF, DIFT
      DOUBLE PRECISION  :: HO, HP, W, DEN
      DOUBLE PRECISION  :: DELTAY
      DOUBLE PRECISION, DIMENSION(N)  :: C, D
C     DIF         = difference used in calculations
C     DIFT        = difference used in calculations
C     HO          = polyinterpolation point
C     HP          = polyinterpolation point
C     W           = weighting between polynomial interpolation points
C     DEN         = difference between consecutive polynomial interpolation points
C     DELTAY      = polynomial interpolation error

      DOUBLE PRECISION, INTENT(IN)  :: X
      DOUBLE PRECISION, INTENT(IN)  :: XA(N), YA(N)
      DOUBLE PRECISION, INTENT(OUT)  :: ERR
      DOUBLE PRECISION, INTENT(OUT)  :: Y
C     XA(N)       = table of XA values used in interpolation
C     YA(N)       = table of YA values used in interpolation
C     ERR         = error in interpolation
C     Y           = interpolated value at X
      
      INTEGER, INTENT(IN)  :: N
C     N           = length of XA and YA arrays

C     Computation Parameters
      DOUBLE PRECISION  :: EPS=1.0E-10
      

C     Variable Initializations
      DELTAY = 0.0

      NS  = 1   
      DIF = DABS(X - XA(1))
      DO I = 1, N
         DIFT = DABS(X - XA(I))
         IF (DIFT < DIF) THEN
            NS  = I
            DIF = DIFT
         END IF
         C(I) = YA(I)
         D(I) = YA(I)
      END DO
      Y  = YA(NS)
      NS = NS - 1    
      DO M = 1, N-1
         DO I = 1, N-M
            HO   = XA(I) - X
            HP   = XA(I + M) - X
            W    = C( I + 1) - D(I)
            DEN  = HO - HP
            D(I) = W * HP / DEN
            C(I) = W * HO / DEN
         END DO
         IF (2 * NS < (N - M)) THEN
            DELTAY = C(NS + 1)
         ELSE
            DELTAY = D(NS)
            NS = NS - 1
         END IF
         Y = Y + DELTAY
      END DO

      ERR = DABS(DELTAY / (Y + EPS))
   
      END SUBROUTINE POLYINTERP

      SUBROUTINE RLEMCONV
C***********************************************************************
C        RLEMCONV Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Fills RLMOVESCONV array with all entries equal
C                     to 1 if FALSE; or a conversion from MOVES 
C                     (g/hr/link) to RLINE native units of (g/m/s),
C                     based on length, if TRUE. 
C
C        PROGRAMMER:  M. Snyder, Wood
C
C        DATE:        May 24, 2018
C    
C        MODIFIED:    
C
C        INPUTS:      
C
C        OUTPUTS:     
C
C        CALLED FROM: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: ISRC, NUMSRC, SRCTYP
      USE RLINE_DATA, ONLY: RLEMISCONV, RLSOURCE, RLMOVESCONV
      IMPLICIT NONE

      INTEGER  :: INDEX
      DOUBLE PRECISION  :: LENGTH
C     INDEX       = index 
C     LENGTH      = length of RLINE source
     
      RLEMISCONV(:) = 1.0d0

C     Perform conversion of MOVES units (g/hr/link) to RLINE units (g/m/s)
C     only for RLINE sources.
      IF(RLMOVESCONV) THEN
         DO INDEX = ISRC, NUMSRC
            IF (SRCTYP(INDEX) .EQ. 'RLINE') THEN
               LENGTH = SQRT((RLSOURCE(INDEX)%XSB - 
     &                  RLSOURCE(INDEX)%XSE)**2 +
     &                  (RLSOURCE(INDEX)%YSB - 
     &                  RLSOURCE(INDEX)%YSE)**2)
               RLEMISCONV(INDEX) = 1.0d0 / LENGTH /3600d0
            END IF
         END DO
      END IF

      END SUBROUTINE RLEMCONV

      DOUBLE PRECISION FUNCTION SIGMAY(XD)
C***********************************************************************
C        SIGMAY Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate horizontal plume spread. 
C
C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:   XD
C
C         RETURNED:     
C
C        CALLING ROUTINES: POINT_CONC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: OBULEN, USTAR, STABLE
      USE RLINE_DATA, ONLY: SIGMAV, SIGMAY0
      IMPLICIT NONE

      DOUBLE PRECISION  :: SZ
      DOUBLE PRECISION, INTENT(IN)  :: XD
C     SZ          = vertical dispersion
C     XD          = distance between source and receptor in direction parallel to wind

C     External functions:
      DOUBLE PRECISION ::  SIGMAZ


      SZ = SIGMAZ(XD)

C     Calculate SIGMAY based on stability
      IF (STABLE) THEN
         SIGMAY = 1.6D0 * SIGMAV / USTAR * SZ *
     &           (1.0D0 + 2.5D0 * SZ / DABS(OBULEN))  
      ELSE
         SIGMAY = 1.6D0 * SIGMAV / USTAR * SZ *
     &           (1.0D0 + 1.0D0 * SZ / DABS(OBULEN))**(-1.0D0 / 2.0D0) 
      END IF

      SIGMAY = DSQRT(SIGMAY**2.0D0 + SIGMAY0**2.0D0)
      
      END FUNCTION

      DOUBLE PRECISION FUNCTION SIGMAZ(XD) 
C***********************************************************************
C        SIGMAZ Function of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate vertical plume spread,
C                     including source configuration effects. 
C
C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        ARGUMENTS:
C           PASSED:   XD
C
C         RETURNED:     
C
C        CALLING ROUTINES: EFFECTIVE_WIND
C                          MEANDER
C                          POINT_CONC 
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: OBULEN, USTAR, ZI, ISRC, STABLE
      USE RLINE_DATA, ONLY: RLSOURCE, UEFF, SIGMAZ0
      IMPLICIT NONE

      DOUBLE PRECISION  :: SIGMAZ_MAX, XBAR, UTEMP, XDABS
      DOUBLE PRECISION  :: SIGZ
      DOUBLE PRECISION  :: SIGMAZD
      DOUBLE PRECISION  :: SIGMAZB, HB, DW
      DOUBLE PRECISION, INTENT(IN)  :: XD
C     SIGMAZ_MAX  = maximum vertical dispersion
C     XBAR        = absolute value of XD/L
C     UTEMP       = temporary wind speed
C     XDABS       = absolute value of XD
C     SIGZ        = vertical dispersion initial calculation
C     SIGMAZD     = vertical dispersion due to depression
C     SIGMAZB     = vertical dispersion due to barrier
C     HB          = height of barrier
C     DW          = distance between barrier and source

C     Computation Parameters
      DOUBLE PRECISION, PARAMETER :: SQ2PI=0.797885D0

      XBAR       = DABS(XD/OBULEN)
      XDABS      = DABS(XD)
      SIGMAZ_MAX = SQ2PI*ZI
      SIGMAZD    = 0.0D0
      SIGMAZB    = 0.0D0
      UTEMP      = UEFF
      HB         = RLSOURCE(ISRC)%HTWALL
      DW         = (RLSOURCE(ISRC)%DCLWALL - RLSOURCE(ISRC)%DCL)

C     Calculate vertical dispersion curve based on stability
      IF (STABLE) THEN
         SIGZ = 0.57D0 * (USTAR * XDABS / UTEMP) /
     &          (1.0D0 + 3.0D0 * (USTAR / UTEMP) * 
     &          (XBAR)**(2.0D0 / 3.0D0))
      ELSE 
         SIGZ = 0.57D0 * (USTAR * XDABS / UTEMP) *
     &          (1.0D0 + 1.5D0 * (USTAR / UTEMP * XBAR))
      END IF

C     Adjust for depressed roadway, if used
      IF (RLSOURCE(ISRC)%DEPTH < 0.0) THEN
         SIGMAZD = DSQRT((RLSOURCE(ISRC)%DEPTH / 2.15D0)**2)
      END IF
        
C     Adjust for barrier, if used
      IF ((HB .GT. 0.0) .AND. (DABS(DW) .LT. 11.0D0 * HB)) THEN
         SIGMAZB = 0.5D0 * SQ2PI * (HB - 0.125D0 *
     &             (DABS(DW) - 3.0D0 * HB))
      END IF

      SIGZ = DSQRT(SIGMAZD * SIGMAZD + SIGZ *
     &       SIGZ + SIGMAZ0 * SIGMAZ0) + SIGMAZB
      SIGMAZ = MIN(SIGZ, SIGMAZ_MAX)
      
      END FUNCTION

      SUBROUTINE  TRANSLATE_ROTATE
C***********************************************************************
C        TRANSLATE_ROTATE Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Translate and rotate the coordinates so that x-axis
C                     lies along the wind. In addition, this subroutine 
C                     allows the user to specify sources based on a 
C                     centerline and the distance from the 
C                     centerline (DCL).. ie an offset.
C
C        PROGRAMMER:  M. Snyder, Wood
C
C        DATE:        May 24, 2018
C    
C        MODIFIED:    
C
C        INPUTS:      
C
C        OUTPUTS:     
C
C        CALLED FROM: RLCALC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: AXR, AYR, NUMREC, NUMSRC, ISRC, SRCTYP
      USE RLINE_DATA, ONLY: X0, Y0, XSB_ROT, YSB_ROT, XSE_ROT, YSE_ROT,
     &                      THETAW, XRCP_ROT, YRCP_ROT, RLSOURCE
      IMPLICIT NONE

      INTEGER  :: INDEX
      DOUBLE PRECISION :: XR_TRAN, YR_TRAN, THETA_LINE
      DOUBLE PRECISION :: XSB_TRAN, YSB_TRAN, XSE_TRAN, YSE_TRAN
      DOUBLE PRECISION :: DCL
C     INDEX       = index 
C     XR_TRAN     = x-coordinate for translated receptor
C     YR_TRAN     = y-coordinate for translated receptor      
C     THETA_LINE  = angle between wind direction and source
C     XSB_TRAN    = beginning x-coordinate for translated source
C     YSB_TRAN    = beginning y-coordinate for translated source      
C     XSE_TRAN    = end x-coordinate for translated source
C     YSE_TRAN    = end y-coordinate for translated source
C     DCL         = offset distance from center line

C     External functions
      DOUBLE PRECISION, EXTERNAL   :: DEPRESSED_DISPLACEMENT
                 

C     Translate line source origin
      X0 = RLSOURCE(ISRC)%XSB    
      Y0 = RLSOURCE(ISRC)%YSB

C     Translate source and receptor coordinates and then rotate them along wind direction
      DO INDEX = ISRC, NUMSRC
         IF (SRCTYP(INDEX) .EQ. 'RLINE' .OR. 
     &       SRCTYP(INDEX) .EQ. 'RLINEXT') THEN
C           Initialize variables used
            DCL = RLSOURCE(INDEX)%DCL 
  
C           Move initial user coordinate system so the origin is at the beginning of first source             
            XSB_TRAN = RLSOURCE(INDEX)%XSB - X0
            YSB_TRAN = RLSOURCE(INDEX)%YSB - Y0
            XSE_TRAN = RLSOURCE(INDEX)%XSE - X0
            YSE_TRAN = RLSOURCE(INDEX)%YSE - Y0
            THETA_LINE = DATAN2(YSE_TRAN - YSB_TRAN, XSE_TRAN - 
     &                   XSB_TRAN)

C           Account for due east and north lines 
            IF (DSIN(THETA_LINE) .EQ. 0) THEN
               DCL = -DCL
            END IF

C           Determine location of the line that is not within a depression,
C           but is specified in source file with the centerline and distance
C           from the centerline
            IF (DCL .NE. 0.0 .AND. RLSOURCE(INDEX)%DEPTH .EQ. 0.0) THEN 
               XSE_TRAN = XSE_TRAN + DCL * DSIN(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
               YSE_TRAN = YSE_TRAN - DCL * COS(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
               XSB_TRAN = XSB_TRAN + DCL * SIN(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
               YSB_TRAN = YSB_TRAN - DCL * COS(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
            END IF
    
C           Adjustments for near source configurations (depression)
            IF (RLSOURCE(INDEX)%DEPTH < 0.0) THEN 
               XSE_TRAN = XSE_TRAN -
     &                    DEPRESSED_DISPLACEMENT(THETA_LINE, INDEX) *
     &                    DSIN(THETA_LINE)
               YSE_TRAN = YSE_TRAN +
     &                    DEPRESSED_DISPLACEMENT(THETA_LINE, INDEX) *
     &                    DCOS(THETA_LINE)
               XSB_TRAN = XSB_TRAN -
     &                    DEPRESSED_DISPLACEMENT(THETA_LINE, INDEX) *
     &                    DSIN(THETA_LINE)
               YSB_TRAN = YSB_TRAN +
     &                    DEPRESSED_DISPLACEMENT(THETA_LINE, INDEX) *
     &                    DCOS(THETA_LINE)
            END IF

            XSB_ROT(INDEX) = XSB_TRAN * DCOS(THETAW) +
     &                       YSB_TRAN * DSIN(THETAW)
            YSB_ROT(INDEX) = -XSB_TRAN * DSIN(THETAW) +
     &                       YSB_TRAN * DCOS(THETAW)
            XSE_ROT(INDEX) = XSE_TRAN * DCOS(THETAW) +
     &                       YSE_TRAN * DSIN(THETAW)
            YSE_ROT(INDEX) = -XSE_TRAN * DSIN(THETAW) +
     &                       YSE_TRAN * DCOS(THETAW)

         END IF 
      END DO

      DO INDEX = 1, NUMREC
         XR_TRAN = AXR(INDEX) - X0
         YR_TRAN = AYR(INDEX) - Y0
         XRCP_ROT(INDEX) = XR_TRAN * DCOS(THETAW) + 
     &                     YR_TRAN * DSIN(THETAW)
         YRCP_ROT(INDEX) = -XR_TRAN * DSIN(THETAW) + 
     &                     YR_TRAN * DCOS(THETAW)
      END DO
      
      END SUBROUTINE TRANSLATE_ROTATE
