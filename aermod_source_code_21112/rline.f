      SUBROUTINE RLCALC
C***********************************************************************
C                 RLCALC Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Introduce RLINE source. Based on RLINEv1.2 released 
C                     November 2013. 
C
C        PROGRAMMER:  M. Snyder, R. Cleary, J. Paumier, R. Wagoner, Wood
C
C        DATE:        July 20, 2018
C
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
C
C        MODIFIED:    Added Urban option for RLINE sources. 
C                     Wood, 03/04/2019
C                     Corrected processing of EMISFACT for RLINE sources.
C                     D42, Wood, 06/19/2020
C
C        INPUTS:      None
C
C        OUTPUTS:    
C
C        CALLED FROM: MAIN
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1, ONLY: IREC, AZFLAG, HRVAL, NUMTYP,
     &                 NUMGRP, IGRP, ISRC, IREC, NUMREC, L_HRLYSIG,
     &                 ITYP, ARM2, CHI, EVONLY, AQS, AHS, ASZINI, QFLAG,
     &                 URBSRC, NUMURB, STABLE, L_MorningTrans,
     &                 ZI, ZIURB, ZIRUR, ZIMECH, URBAN, SRCTYP,
     &                 USTAR, URBUSTR, RURUSTR, URBSTAB,
     &                 IURB, OBULEN, URBOBULEN, RUROBULEN, IURBGRP,
     &                 WSTAR, URBWSTR, EMIFAC, QTK
C     rline_emisfact_d42_Wood
C     Added import of QTK from MAIN1 above, and new local QEMIS below
      USE RLINE_DATA
      IMPLICIT NONE
      
      INTEGER  :: I
      DOUBLE PRECISION  :: ERROR
      DOUBLE PRECISION  :: CONCD
      DOUBLE PRECISION  :: CONCENTRATION
	    DOUBLE PRECISION  :: QEMIS
C     ERROR         = error in numerical integration
C     CONCD         = dummy concentration at receptor
C     CONCENTRATION = concentration
C     QEMIS         = emission rate
      
      DOUBLE PRECISION  :: XTEMP, YTEMP, ZTEMP
C     XTEMP          = temporary x-coordinate of receptor for line orientation
C     YTEMP          = temporary y-coordinate of receptor for line orientation
C     ZTEMP          = temporary z-coordinate of receptor for line orientation

C     Variable Initializations:
      CONCD = 0.0D0
      CONCENTRATION = 0.0D0
     
C     Initialize __VAL arrays
      HRVAL = 0.0D0

C     rline_emisfact_d42_Wood
C     Local QEMIS will hold the hourly source specific emission to be
C     used in the calculation of concentration below, replacing the use
C     of RLSOURCE(ISRC)%QEMIS in those expressions.
      QEMIS = RLSOURCE(ISRC)%QEMIS
      IF(QFLAG(ISRC) .EQ. 'HOURLY') THEN
C        Set hourly variable source height and sigmaz-initial (optional)
         IF(L_HRLYSIG(ISRC)) THEN
           RLSOURCE(ISRC)%ZSB = AHS(ISRC)
           RLSOURCE(ISRC)%ZSE = AHS(ISRC)
           RLSOURCE(ISRC)%INIT_SIGMAZ = ASZINI(ISRC) 
         END IF
C        Set hourly variable emission rate (required)	  
         IF(SRCTYP(ISRC) .EQ. 'RLINEXT') THEN
           QEMIS = AQS(ISRC)
         ELSE ! RLINE source with Lnemis inputs
           QEMIS = AQS(ISRC)*RLSOURCE(ISRC)%WIDTH
         END IF
C     rline_emisfact_d42_Wood
C     The following block is added to process the EMISFACT keyword with
C     RLINE sources, set QTK equal to appropriate emission factor, and
C     multiply by emission rate.
      ELSE IF ((QFLAG(ISRC) .EQ. 'MONTH') .OR.
     &         (QFLAG(ISRC) .EQ. 'HROFDY') .OR.
     &         (QFLAG(ISRC) .EQ. 'WSPEED') .OR.
     &         (QFLAG(ISRC) .EQ. 'SEASON') .OR.
     &         (QFLAG(ISRC) .EQ. 'SEASHR') .OR.
     &         (QFLAG(ISRC) .EQ. 'HRDOW') .OR.
     &         (QFLAG(ISRC) .EQ. 'HRDOW7') .OR.
     &         (QFLAG(ISRC) .EQ. 'SHRDOW') .OR.
     &         (QFLAG(ISRC) .EQ. 'SHRDOW7') .OR.
     &         (QFLAG(ISRC) .EQ. 'MHRDOW') .OR.
     &         (QFLAG(ISRC) .EQ. 'MHRDOW7')) THEN
         CALL EMFACT(1.0D0)
         QEMIS = RLSOURCE(ISRC)%QEMIS*QTK
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
  
      IF(.NOT. RLPROCESSED) THEN
C        Translate and rotate the line source to align with wind direction   --- CALL TRANSLATE_ROTATE
         CALL TRANSLATE_ROTATE
      END IF

C        Get translated an rotated source information
         SIGMAZ0 = RLSOURCE(ISRC)%INIT_SIGMAZ 
         XSBEGIN = XSB_ROT(ISRC)
         YSBEGIN = YSB_ROT(ISRC)
         XSEND   = XSE_ROT(ISRC)  
         YSEND   = YSE_ROT(ISRC)
         ZSBEGIN = RLSOURCE(ISRC)%ZSB
         ZSEND   = RLSOURCE(ISRC)%ZSE

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
      THETA_LINE = DATAN2(YSEND - YSBEGIN, XSEND - XSBEGIN)

C        Calculate a vertical displacement of the source to  
C        account for the effect of a roadside barrier                        --- CALL BARRIER_DISPLACEMENT
      IF((RLSOURCE(ISRC)%HTWALL  .GT. 0.0D0) .OR. 
     &   (RLSOURCE(ISRC)%HTWALL2 .GT. 0.0D0)) THEN
          CALL BARRIER_DISPLACEMENT
      END IF

C     Calculate initial sigma-y from the road width
      SIGMAY0 = 0.0D0
      SIGMAY0 = DABS(0.5D0 * (RLSOURCE(ISRC)%WIDTH) *
     &          DCOS(THETA_LINE))
     
C     Calculate met parameters                                               --- CALL COMPUTE_MET
      CALL COMPUTE_MET 

C     Begin loop over receptors
      RECEPTOR_LOOP: DO IREC = 1, NUMREC
C        Rotate X, Y receptors.  Z receptor not rotated.
         XR_ROT = XRCP_ROT(IREC)
         YR_ROT = YRCP_ROT(IREC)
         ZRECEP = AZFLAG(IREC)

C        Calculate the contribution to concentration at a 
C        receptor due to a line source using Romberg integration             --- CALL NUMERICAL_LINE_SOURCE

         CALL NUMERICAL_LINE_SOURCE(CONCD, ERROR)

C        Convert Qemis from MOVES units (g/hr/link) to RLINE units (g/m/s)
C        EMIFAC(1) is for concentration
C        rline_emisfact_d42_Wood, using local QEMIS
         CONCENTRATION = CONCD * EMIFAC(1) * QEMIS * RLEMISCONV(ISRC)
         HRVAL(:) = CONCENTRATION
                     
C        For the initial integration of R-LINE v1.2 into AERMOD, 
C        only including ARM2 chemistry options.  OLM, PVMRM and GRSM not included.
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

      SUBROUTINE BARRIER_DISPLACEMENT
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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
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

C     Variable Declarations:
CCRT 3/24/2021, D058 2 Barriers - add variables for error handling
      USE MAIN1, ONLY: ISRC, SFCZ0, SRCID, PATH
      USE RLINE_DATA, ONLY: RLSOURCE, ZSBEGIN, ZSEND, 
     &                      THETA_LINE, SM_NUM, NBARR, BDW_FLAG, 
     &                      DWU, DW_PERU, HBU, DWD, DW_PERD, HBD, 
     &                      XSHIFT, SHIFT_FLAG, ALPHA_U, ALPHA_D

      IMPLICIT NONE

C     Local Variables:
      DOUBLE PRECISION      :: HB(2), DW(2), DW_PER(2)
      DOUBLE PRECISION      :: RECIRC_LENGTH, HBMAX, Z0
      INTEGER               :: I   
C     HB            = height of barrier; for up to 2 barriers
C     DW            = along wind distance between source and barrier; for up to 2 barriers
C     DW_PER        = perpendicular distance between source and barrier; for up to 2 barriers
C     RECIRC_LENGTH = length of recirculation zone; either 0, 4H, or 6.5H
C     HBMAX         = maximum HB; needed for 2 barrier case to pick the highest barrier height
C     Z0            = local surface roughness for calculation of UH
C     I             = index to indicate if barrier 1 or barrier 2
      
C     Initialize Variables:
      I        = 2     ! assume barrier 2 is present, correct if necessary
      HBD      = 0.0D0 ! height of downwind barrier
      DW_PERD  = 0.0D0 ! perpendicular distance between source and downwind barrier
      DWD      = 0.0D0 ! along wind distance between source and downwind barrier
      
      HBU      = 0.0D0 ! height of upwind barrier
      DW_PERU  = 0.0D0 ! perpendicular distance between source and upwind barrier
      DWU      = 0.0D0 ! along wind distance between source and upwind barrier
      
      NBARR    = 0     ! number of barriers present
      
      SHIFT_FLAG    = .FALSE.
      XSHIFT        = 0.0D0
      RECIRC_LENGTH = 0.0D0 
      
      ALPHA_U = 1.0D0
      ALPHA_D = 1.0D0
      Z0      = SFCZ0  ! initialize local zrough for UH calculations
      HBMAX   = 0.0D0  ! needed for 2 barrier case
      
C     Check for existence of barriers in user input
      IF((RLSOURCE(ISRC)%HTWALL  .GT. 0.0D0) .OR.
     &    (RLSOURCE(ISRC)%HTWALL2 .GT. 0.0D0)) THEN 
       
C     Calculate barriers distances and source height
          ZSBEGIN     = 0.5D0 * (ZSBEGIN + ZSEND)
          ZSEND       = ZSBEGIN

          HB(1)       = RLSOURCE(ISRC)%HTWALL
          DW_PER(1)   = DABS(RLSOURCE(ISRC)%DCLWALL - 
     &                      RLSOURCE(ISRC)%DCL)
          DW(1)       = DW_PER(1) / (DABS(DSIN(THETA_LINE)) + SM_NUM)     ! distance between source and barrier along wind direction
          
          HB(2)       = RLSOURCE(ISRC)%HTWALL2
          DW_PER(2)   = DABS(RLSOURCE(ISRC)%DCLWALL2 - 
     &                      RLSOURCE(ISRC)%DCL)
          DW(2)       = DW_PER(2) / (DABS(DSIN(THETA_LINE)) + SM_NUM)
  
C     Determine number of barriers (0, 1, or 2)
          IF ((HB(1) > 0.0D0) .AND. (HB(2) > 0.0D0)) THEN           ! TWO BARRIERS
            NBARR = 2
          ELSE IF ((HB(1) == 0.0D0) .AND. (HB(2) == 0.0D0)) THEN    ! NO BARRIERS
            NBARR = 0
          ELSE IF ((HB(1) > 0.0D0) .OR. (HB(2) == 0.0D0)) THEN      ! ONE BARRIER
            NBARR = 1
          ELSE
            NBARR = -1
          END IF
          
C     Assigning barriers 1 or 2 to upwind or downwind; setting barrier displacement variables for each case 
        SELECT CASE(NBARR) 
          CASE(0) ! NO BARRIERS
            DW_PERU   = 0.0D0
            DWU       = 0.0D0
            DW_PERD   = 0.0D0
            DWD       = 0.0D0
            
          CASE(1) ! ONE BARRIER
            IF (HB(1) > 0.0D0) I = 1
            IF(BDW_FLAG(ISRC,I) .EQ. 1) THEN                            ! located on downwind side
              HBD     = HB(I)
              DWD     = DW(I)
              DW_PERD = DW_PER(I)
              HBU     = 0.0D0
              DWU     = 0.0D0
              DW_PERU = 0.0D0
              Z0      = MAX(HBD / 9.0D0, SFCZ0)                         ! adjusting surface reference for presence of barrier
              ALPHA_D = MAX(1.0D0, EXP(0.14D0 * LOG(Z0 / SFCZ0)))       ! Venkatram & Schulte (2018) alpha; Note: exp(n*log(x)) = x**n
            ELSE                                                        ! located on upwind side
              HBD     = 0.0D0
              DWD     = 0.0D0
              DW_PERD = 0.0D0
              HBU     = HB(I)
              DWU     = DW(I)
              DW_PERU = DW_PER(I) 
              Z0      = MAX(HBU / 9.0D0, SFCZ0)                         ! adjusting surface reference for presence of barrier
              ALPHA_U = MAX(1.0D0, EXP(0.14D0 * LOG(Z0 / SFCZ0)))       ! Venkatram & Schulte (2018) alpha; Note: exp(n*log(x)) = x**n          
              RECIRC_LENGTH = 6.5D0 * HBU                               ! recirculation zone is 6.5h when one barrier present
            END IF
            
          CASE(2) ! TWO BARRIERS
            IF(BDW_FLAG(ISRC, 1) .EQ. BDW_FLAG(ISRC, 2)) THEN           ! both are on same side of source
              PRINT *, "WARNING: Both barriers associated with source ",
     &               SRCID(ISRC)," are on the same side of the source."
              PRINT *, "Barrier closer to source will be used."
CCRT 3/24/2021, D058 2 Barriers - write warning to .out/.err file
              CALL ERRHDL(PATH,'RLINE','W','620',SRCID(ISRC))
              NBARR     = 1
              I         = 2                                             ! assume barrier 2 is close; change in next line if not
              IF(DABS(DW_PER(1)) < DABS(DW_PER(2))) I = 1               ! dw_per(i) is closer to source
              
              IF (BDW_FLAG(ISRC, I) .EQ. 1) THEN                        ! barriers are downwind of source but only choosing closest barrier
                HBD     = HB(I)
                DWD     = DW(I)
                DW_PERD = DW_PER(I)
                HBU     = 0.0D0
                DWU     = 0.0D0
                DW_PERU = 0.0D0
                Z0      = MAX(HBD / 9.0D0, SFCZ0)                       ! adjusting surface reference for presence of barrier
                ALPHA_D = MAX(1.0D0, EXP(0.14D0 * LOG(Z0 / SFCZ0)))     ! Venkatram & Schulte (2018) alpha
              ELSE                                                      ! barriers are upwind of the source but only choosing closest barrier
                HBD     = 0.0D0
                DWD     = 0.0D0
                DW_PERD = 0.0D0
                HBU     = HB(I)
                DWU     = DW(I)
                DW_PERU = DW_PER(I)
                Z0      = MAX(HBU / 9.0D0, SFCZ0)                       ! adjusting surface reference for presence of barrier
                ALPHA_U = MAX(1.0D0, EXP(0.14D0 * LOG(Z0 / SFCZ0)))     ! Venkatram & Schulte (2018) alpha          
                RECIRC_LENGTH = 6.5D0 * HBU                             ! assuming one upwind barrier so recirculation zone is 6.5H
              END IF
            
            ELSE                                                        ! barriers are on opposite sides of source
              IF (BDW_FLAG(ISRC,1) .EQ. 1) THEN                         ! barrier 1 is downwind, barrier 2 is upwind
                HBD     = HB(1)
                DW_PERD = DW_PER(1)
                DWD     = DW(1)
                HBU     = HB(2)
                DW_PERU = DW_PER(2)
                DWU     = DW(2)
                RECIRC_LENGTH = 4.0D0 * HBU                             ! recirculation zone is 4H when two barriers present
              ELSE                                                      ! barrier 2 is downwind, barrier 1 is upwind
                HBD     = HB(2)
                DW_PERD = DW_PER(2)
                DWD     = DW(2)
                HBU     = HB(1)
                DW_PERU = DW_PER(1)
                DWU     = DW(1)
                RECIRC_LENGTH = 4.0D0 * HBU                             ! recirculation zone is 4H when two barriers present
              END IF
              Z0      = MAX(HBU / 9.0D0, SFCZ0)                         ! adjusting surface reference for presence of barrier
              ALPHA_U = MAX(1.0D0, EXP(0.14D0 * LOG(Z0 / SFCZ0)))       ! Venkatram & Schulte (2018) alpha
              Z0      = MAX(HBD / 9.0D0, SFCZ0)                         ! adjusting surface reference for presence of barrier
              ALPHA_D = MAX(1.0D0, EXP(0.14D0 * LOG(Z0 / SFCZ0)))       ! Venkatram & Schulte (2018) alpha
            END IF
            
          CASE DEFAULT 
          PRINT *, "WARNING: Barriers input problem with source ", ISRC
            
        END SELECT
  
C     Set Shift_flag for cases when there is an upwind barrier present 
        IF((HBU > 0.0D0) .AND. (DW_PERU < RECIRC_LENGTH)) THEN
           SHIFT_FLAG = .TRUE.
           XSHIFT     = DWU
           DWD        = DWD + XSHIFT                                    ! adjusting distance to downwind barrier due to shifting source to upwind barrier
        ELSE
          SHIFT_FLAG  = .FALSE.
          XSHIFT      = 0.0D0
        END IF
      
      END IF
 
C      Possible future debugging print statement
C      PRINT *, "Xshift = ", Xshift, "Shift_flag = ", Shift_flag 

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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
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

C     Variable Declarations:
      USE MAIN1, ONLY: USTAR, SFCZ0, OBULEN, UREFHT, UREF
      USE RLINE_DATA, ONLY: SIGMAV, FAC_DISPHT, DISPHT, RLWSTAR, 
     &                      WSPD_ADJ, UH, HBD, HBU, I_ALPHA,
     &                      ALPHA_D, ALPHA_U,
     &                      Z0_A, DH_A, UST_A, LMO_A
      IMPLICIT NONE

C     External Functions:
      DOUBLE PRECISION, EXTERNAL  :: MOST_WIND

C     Local Variables:
      DOUBLE PRECISION  :: SIGMAV_CALC, WSTAR_LOC, UREFCALC
      DOUBLE PRECISION  :: ALPHAS(3)
C     SIGMAV_CALC = SIGMAV calculated from WSTAR and USTAR
C     WSTAR_LOC   = check for WSTAR variable
C     UREFCALC    = theoretical value of wind speed at z = UREFHT
C     ALPHAS      = enhancement of ustar (UST) due to barrier presence 

C     Variable Initialization:
      UH        = 0.0D0   ! wind speed at barrier height (H)

C     Check for WSTAR from metext.f
      WSTAR_LOC = MAX(RLWSTAR, 0.0D0)
    
C     Calculate SIGMAV from WSTAR and USTAR variables
C     Calculate standard deviation of wind direction
      SIGMAV_CALC = DSQRT((0.6D0 * WSTAR_LOC)**2 + 
     &              (1.9D0 * USTAR)**2)
      SIGMAV      = MAX(SIGMAV_CALC, 0.2D0) 

C     Calculate z0 and dh for three cases - no barrier and two possible barrier heights
      Z0_A(1)  = SFCZ0
      Z0_A(2)  = MAX(HBD / 9.0D0, SFCZ0)                                  ! adjusting surface reference for presence of barrier
      Z0_A(3)  = MAX(HBU / 9.0D0, SFCZ0)                                  ! adjusting surface reference for presence of barrier
      DH_A(:)  = FAC_DISPHT * Z0_A(:)
      DISPHT = DH_A(1)
      ALPHAS(1) = 1.0D0
      ALPHAS(2) = ALPHA_D
      ALPHAS(3) = ALPHA_U
      UST_A(:)    = USTAR * ALPHAS(:)
      LMO_A(:)    = OBULEN * ALPHAS(:)**3
C                                                                            --- CALL MOST_WIND
      I_ALPHA  = 1
      UREFCALC = MOST_WIND(UREFHT,Z0_A(I_ALPHA),DH_A(I_ALPHA),
     &           UST_A(I_ALPHA), LMO_A(I_ALPHA))
      WSPD_ADJ = UREF / UREFCALC
      
      IF (HBD .GT. 0.0D0) THEN        
C                                                                            --- CALL MOST_WIND
          I_ALPHA = 2
          UH      = MOST_WIND(UREFHT,Z0_A(I_ALPHA),DH_A(I_ALPHA), 
     &              UST_A(I_ALPHA), LMO_A(I_ALPHA)) * WSPD_ADJ
      END IF

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

C     Variable Declarations:
      USE RLINE_DATA, ONLY: XEXP, AEXP, BEXP, DELEXP
      IMPLICIT NONE

C     Local Variables:
      INTEGER  :: IND
      DOUBLE PRECISION, DIMENSION(1000)  :: EXT
C     IND         = local source index
C     EXT         = exponent

C     Create look-up table
      DELEXP  = 20.0D0 / 999.0D0
      XEXP(1) = -20.0D0
      EXT(1)  = DEXP(XEXP(1))

      DO IND = 2, 1000
         XEXP(IND) = XEXP(IND - 1) + DELEXP
         EXT(IND)  = DEXP(XEXP(IND))
      END DO

      DO IND = 1, 999
         BEXP(IND) = (EXT(IND + 1) - EXT(IND)) / 
     &               (XEXP(IND + 1) - XEXP(IND))
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
C        INPUTS:      THETA_LINE, IND
C
C        OUTPUTS:      
C                    
C        CALLED FROM: TRANSLATE_ROTATE
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE RLINE_DATA, ONLY: RLSOURCE, THETAW
      IMPLICIT NONE

C     Local Variables:
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
 
      IF (RELD .GE. 0.0483D0) THEN
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
     &                         DSIGN(1.0D0, DSIN(THETA_LINE))

      END FUNCTION

      SUBROUTINE EFFECTIVE_WIND(XD,HEFF,HSHIFT)
C***********************************************************************
C        EFFECTIVE_WIND Module of the AMS/EPA Regulatory Model - AERMOD
C
C        PURPOSE:     Calculate the effective wind speed at mean
C                     plume height.
C
C        PROGRAMMER:  A. Venkatram, M. Snyder, D. Heist
C
C        DATE:        November, 2013
C    
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
C
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      XD, HEFF, HSHIFT
C
C        OUTPUTS:     
C
C        CALLED FROM: MEANDER
C                     POINT_CONC
C
C        CALLING
C        ROUTINES:    EXPX
C                     MOST_WIND
C                     SIGMAZ
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1, ONLY: RTOF2, RT2BYPI
      USE RLINE_DATA, ONLY: UEFF, SIGMAV, WSPD_ADJ, DISPHT, 
     &    I_ALPHA, Z0_A, DH_A, UST_A, LMO_A 

      IMPLICIT NONE

C     External Functions:
      DOUBLE PRECISION, EXTERNAL  :: SIGMAZ, MOST_WIND, EXPX
      
C     Local Variables:
      INTEGER  :: ITER
      DOUBLE PRECISION  :: ERF
      DOUBLE PRECISION  :: SZ, SZ_NEW, ERR, ZBAR
      DOUBLE PRECISION, INTENT(IN)  :: XD, HEFF, HSHIFT
C     ITER        = iteration
C     ERF         = error function
C     SZ          = effective SIGMAZ
C     SZ_NEW      = intermediate vertical dispersion
C     ERR         = error in each successive calculation
C     ZBAR        = mean plume height
C     XD          = perpendicular distance between source and receptor
C     HEFF        = effective source height
C     HSHIFT      = vertical shift in source height for barriers

C     Initialize variables:
      ZBAR = 0.0D0
      ERR  = 10.0D0
      ITER = 1 
      SZ   = SIGMAZ(XD)
C                                                                            --- CALL MOST_WIND    
      UEFF = MOST_WIND(HEFF,Z0_A(I_ALPHA),DH_A(I_ALPHA), 
     &       UST_A(I_ALPHA), LMO_A(I_ALPHA)) * WSPD_ADJ
C                                                                            --- CALL EXPX      
      DO WHILE ((ERR > 1.0E-02) .AND. (ITER < 20)) 
         ZBAR   = RT2BYPI * SZ * EXPX(-0.5D0 * (HEFF / SZ)**2) +
     &            HEFF * ERF(HEFF / (RTOF2 * SZ)) + HSHIFT                   ! Venkatram et al. (2013)
C                                                                            --- CALL MOST_WIND
         UEFF   = MOST_WIND(MAX(ZBAR, HEFF), Z0_A(I_ALPHA),
     &            DH_A(I_ALPHA), UST_A(I_ALPHA), LMO_A(I_ALPHA)) *
     &            WSPD_ADJ
         UEFF   = DSQRT(2.0D0 * SIGMAV**2  + UEFF**2)
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
C        INPUTS:      XP
C
C        OUTPUTS:     
C
C        CALLED FROM: EFFECTIVE_WIND
C                     MEANDER
C                     POINT_CONC
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE RLINE_DATA, ONLY: DELEXP, AEXP, BEXP, XEXP
      IMPLICIT NONE

C     Local Variables:
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

      DOUBLE PRECISION FUNCTION MEANDER(X, Y, Z, HSHIFT)
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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      X, Y, Z, HSHIFT
C
C        OUTPUTS:
C
C        CALLED FROM: POINT_CONC     
C
C        CALLING 
C        ROUTINES:    EXPX
C                     EFFECTIVE_WIND
C                     SIGMAZ
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1,      ONLY: PI, RT2BYPI
      USE RLINE_DATA, ONLY: XR_ROT, YR_ROT, XD_MIN, ZRECEP, UEFF,
     &                      XD_MIN, HBD, HBU
      IMPLICIT NONE

C     External functions:
      DOUBLE PRECISION, EXTERNAL :: SIGMAZ, EXPX

C     Local Variables:
      DOUBLE PRECISION  :: R, VERT, HORZ, SZ, HEFF
      DOUBLE PRECISION, INTENT(IN)  :: X, Y, Z, HSHIFT
C     R           = radial distance to receptor
C     VERT        = vertical component of concentration
C     HORZ        = horizontal component of concentration
C     SZ          = effective SIGMAZ
C     HEFF        = effective height
C     X           = x-coordinate of source location
C     Y           = y-coordinate of source location
C     Z           = z-coordinate of source location
C     HSHIFT      = vertical shift in source height for barriers

      R       = DSQRT((XR_ROT - X)**2 + (YR_ROT - Y)**2)                     ! radial distance to the receptor
      R       = MAX(R, XD_MIN)
      HEFF    = MAX(Z, 0.75D0 * MAX(HBU, HBD))                               ! if no barrier, heff = Z; if barrier, heff = 0.75*H using largest H 

C     Calculate effective wind speed                                         --- CALL EFFECTIVE_WIND
      CALL EFFECTIVE_WIND(R,HEFF,HSHIFT)

C     Account for source height                                              --- CALL EXPX
      SZ      = SIGMAZ(R)
      VERT    = RT2BYPI * (EXPX(-0.5D0 * ((HEFF - ZRECEP)
     &          / SZ)**2) + EXPX(-0.5D0 * 
     &          ((HEFF + ZRECEP) / SZ)**2))/(2.0D0*SZ*UEFF)

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
C        INPUTS:      Z, Z0, DH, UST, L
C
C        OUTPUTS:     
C
C        CALLED FROM: COMPUTE_MET
C                     EFFECTIVE_WIND
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations
      USE MAIN1, ONLY: STABLE, PI, RTOF2, VONKAR
      USE RLINE_DATA, ONLY:  SIGMAV
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

      IF (DH .GE. Z) THEN
         MOST_WIND = RTOF2 * SIGMAV
      ELSE
         IF (STABLE) THEN
            L_LOC = DABS(L)
            PSI1 = -17.0D0 * (1.0D0 - DEXP(-0.29D0 * (Z - DH) / L_LOC))
            PSI2 = -17.0D0 * (1.0D0 - DEXP(-0.29D0 * Z0 / L_LOC))
         ELSE
            L_LOC = -1.0D0*DABS(L)
            X1   = SQRT(SQRT(1.0d0 - 16.0d0 * (Z - DH) / L_LOC))
            X2   = SQRT(SQRT(1.0D0 - 16.0D0 * Z0 / L_LOC))
            PSI1 = 2.0D0 * DLOG((1.0D0 + X1) / 2.0D0) + 
     &             DLOG((1.0 + X1 * X1) / 2.0D0) -
     &             2.0D0 * DATAN(X1) + PI / 2.0D0
            PSI2 = 2.0D0 * DLOG((1.0D0 + X2) / 2.0D0) +
     &             DLOG((1.0 + X2 * X2) / 2.0D0) -
     &             2.0D0 * DATAN(X2) + PI / 2.0D0
         END IF
         MOST_WIND = UST*(DLOG((Z - DH) / Z0) - PSI1 + PSI2) / VONKAR
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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
C
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      
C
C        OUTPUTS:     CONC_NUM, ERR 
C
C        CALLED FROM: RLCALC
C
C        CALLING 
C        ROUTINES:    POINT_CONC
C                     POLYINTERP  
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE RLINE_DATA, ONLY: XSBEGIN, YSBEGIN, ZSBEGIN, XSEND, YSEND,
     &                      ZSEND, XR_ROT, YR_ROT, SM_NUM, XPER,
     &                      ERROR_LIMIT, XD_MIN, THETA_LINE
      IMPLICIT NONE

C     External Functions:
      DOUBLE PRECISION, EXTERNAL  :: POINT_CONC

C     Local Variables:
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
      DOUBLE PRECISION  :: X, Y, Z, DELT, PHI
      DOUBLE PRECISION  :: T, TMAX, COST, SINT, SINP, COSP, XREC
      DOUBLE PRECISION  :: DR
C     XDIF           = x integral limit 
C     YDIF           = y integral limit
C     ZDIF           = z integral limit
C     DISP           = dummy variable used to store integral
C     HDUM           = step size
C     CONCDUM        = successive concentration approximations
C     CONC_INT       = numerical integrals
C     X              = x-coordinate of source location
C     Y              = y-coordinate of source location
C     Z              = z-coordinate of source location
C     DELT           = distance between points used to estimate line source
C     PHI            = angle of line elevation
C     T              = half of DELT
C     TMAX           = 3-dimensional 
C     COST           = cosine of theta_line
C     SINT           = sin of theta_line
C     SINP           = sin of phi
C     COSP           = cosine of phi
C     XREC           = x-coordinate of point on source directly upwind of receptor
C     DR             = distance beetween source and receptor in wind direction
      
      DOUBLE PRECISION, INTENT(OUT)  :: CONC_NUM
      DOUBLE PRECISION, INTENT(OUT)  :: ERR
C     CONC_NUM       = numerical routine output concentration 
C     ERR            = integration is set to an arbitrarily large value before it is reduced
      
      DOUBLE PRECISION, ALLOCATABLE  :: H(:), CONC(:)
C     H              =  step size
C     CONC           =  successive concentration approximations

C     2K is the order of Romberg integration.  Nmax needs to be greater than K for Romberg integration to be used, 
C     otherwise trapezoidal integration is used
      INTEGER, PARAMETER  :: K = 3 

C     Computation Parameters
      DOUBLE PRECISION, PARAMETER :: XINTERP = 0.0D0, A = 1.0d0

C     Variable Initializations
      CONC_NUM = 0.0D0

      XDIF  = XSEND - XSBEGIN  
      YDIF  = YSEND - YSBEGIN
      ZDIF  = ZSEND - ZSBEGIN

      TMAX  = DSQRT(XDIF * XDIF + YDIF * YDIF + ZDIF * ZDIF)    
      PHI   = DASIN(ZDIF / (TMAX + SM_NUM))  
      SINP  = DSIN(PHI)
      COSP  = DCOS(PHI)
      COST  = DCOS(THETA_LINE)
      SINT  = DSIN(THETA_LINE)

C     Find x-coordinate of point on source line directly upwind of receptor
      XREC = (XSEND - (YSEND - YR_ROT) * (XSEND - XSBEGIN) / 
     &       (YSEND - YSBEGIN))
      DR   = DABS(XR_ROT - XREC)
      
C     Prevent user from placing receptor on source
      DR   = MAX(XD_MIN, DR)
      XPER = DSIGN(A, XR_ROT - XREC) * DABS(-1.0D0 * YDIF * 
     &       (XR_ROT - XSBEGIN) + XDIF * (YR_ROT - YSBEGIN)) /
     &       DSQRT(XDIF**2 + YDIF**2)

C     Convergence Criteria: Minimum Iterations    
      IF ((YR_ROT > YSBEGIN - TMAX / 2.0D0 * DABS(COST)) .AND. 
     &   (YR_ROT < YSEND + TMAX / 2.0D0 * DABS(COST))) THEN
         MINJ = CEILING(DLOG(2.0D0 * TMAX / 
     &          (MAX(XD_MIN, DR * DABS(SINT))) - 
     &          2.0D0) / DLOG(2.0D0)) + 2.0D0
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
      CONC(:) = 0.0D0
      H(:)    = 0.0D0
C                                                                            --- CALL POINT_CONC
      DISP    = (POINT_CONC(XSBEGIN, YSBEGIN, ZSBEGIN) + 
     &          POINT_CONC(XSEND, YSEND, ZSEND)) * 0.5D0
C     Calculate first approximation of the integration.  Set relative size of
C     integration interval
      CONC(1) = DISP * TMAX
      H(1)    = 1.0D0
C     Trapezoidal integration 
      DO J = 2, IT_LIM 
         NO_POINTS = 2**(J - 2)    
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
            CONC_NUM = DABS(CONC_INT)

C           Check convergence criteria          
            IF ((DABS(ERR) < ERROR_LIMIT) .AND. (J > MINJ)) THEN
               DEALLOCATE(H, CONC, STAT = ALLOCERROR)
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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      X, Y, Z
C
C        OUTPUTS: 
C       
C        CALLED FROM: NUMERICAL_LINE_SOURCE    
C
C        CALLING 
C        ROUTINES:   EFFECTIVE_WIND 
C                    EXPX
C                    MEANDER
C                    SIGMAY
C                    SIGMAZ
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1,      ONLY: RT2BYPI, SRT2PI
      USE RLINE_DATA, ONLY: XR_ROT, YR_ROT, ZRECEP, UEFF, 
     &                      XSBEGIN, XSEND, XR_ROT, YR_ROT, 
     &                      ZRECEP, HBD, DWD, DW_PERD, HBU, 
     &                      SIGMAZ0, SIGMAV, SZB, UEFF, XPER, 
     &                      NBARR, UH, XSHIFT, SHIFT_FLAG, XD_MIN,
     &                      I_ALPHA, ALPHA, ALPHA_U, ALPHA_D
      IMPLICIT NONE

C     External Functions:
      DOUBLE PRECISION, EXTERNAL  :: MEANDER, EXPX
      DOUBLE PRECISION, EXTERNAL  :: SIGMAY, SIGMAZ

C     Local Variables:
      DOUBLE PRECISION  :: CONC_M, CONC_P, VERT, HORZ
      DOUBLE PRECISION  :: SY, SZ, SIGMAZ0_ORIG
      DOUBLE PRECISION  :: FRAN
      DOUBLE PRECISION  :: XD, YD, XMAX, DWDT
      DOUBLE PRECISION  :: HEFF, HSHIFT, HBMAX, HBDEFF
      DOUBLE PRECISION  :: CQ, FQ, UEFF_BU
      DOUBLE PRECISION  :: F_UEFF, F_UH, HMAX, A
      INTEGER           :: IA_MAX 
      DOUBLE PRECISION, INTENT(IN)  :: X, Y, Z
C     CONC_M      = meander concentration
C     CONC_P      = direct plume concentration
C     VERT        = vertical component of concentration
C     HORZ        = horizontal component of concentration
C     SY            = horizontal dispersion coefficient
C     SZ            = vertical dispersion coefficient
C     SIGMAZ0_ORIG  = reset sigmaz0 to the value read from the source input
C     FRAN          = fraction between meander and direct plume
C     XD            = distance between source and receptor in direction parallel to wind
C     YD            = distance between source and receptor in direction perpendicular to wind
C     XMAX          = distance between downwind barrier (at the most upwind point) and the receptor
C     DWDT          = distance between the source and downwind barrier
C     HEFF          = effective source height
C     HSHIFT        = vertical shift in source height for barriers
C     HBMAX         = taller barrier height, either HBU or HBD
C     HBDEFF        = effective downwind barrier height in the mixed-wake algorithm
C     CQ            = variable used to compute concentration for mixed-wake algorithm; Venkatram et al. (2021)
C     FQ            = scaling factor for mixed-wake algorithm; Venkatram et al. (2021)
C     UEFF_BU       = factor for reducing ueff immediately downwind of the upwind barrier
C     F_UEFF        = factor for reducing ueff in upwind barrier shift cases
C     F_UH          = factor for UH
C     HMAX          = tallest barrier height for scaling of UH
C     A             = variable used in UH
C     IA_MAX        = index for alpha for the barrier with greater height (i.e., upwind or downwind barrier)
C     X             = x-coordinate of source location
C     Y             = y-coordinate of source location
C     Z             = z-coordinate of source location

C     Declare flags:
      logical ::  Direct_flag     ! TRUE = direct plume; FALSE = meander only plume
      logical ::  Gaussian_flag   ! TRUE = gaussian mode; FALSE = mixed-wake mode

C     ----------------------------------------------------------------------------------
C     Explaining All Case Combinations:
C     1 --  Flat (no barrier), receptor upwind of source; meander only
C     2 --  Flat (no barrier), receptor downwind of source; direct, no shift, gaussian
C     3 --  1 downwind barrier, receptor upwind of source; meander only
C     4 --  1 downwind barrier, receptor between source and barrier; direct, no shift, gaussian
C     5 --  1 downwind barrier, receptor downwind of barrier; direct, no shift, mixed-wake
C     6 --  1 upwind barrier, receptor upwind of barrier; meander only
C     7 --  1 upwind barrier, receptor downwind of source, source within recirc zone; direct, shift, gaussian
C     8 --  1 upwind barrier, receptor upwind of source; meander only
C     9 --  1 upwind barrier, receptor downwind of source, source outside of recirc zone; direct, no shift, gaussian
C     10 -- 2 barriers, receptor upwind of upwind_barrier; meander only
C     11 -- 2 barriers, receptor between source and downwind barrier, source within recirc zone; direct, shift, gaussian
C     12 -- 2 barriers, receptor upwind of source; meander only
C     13 -- 2 barriers, receptor between source and downwind barrier, source outside of recirc zone; direct, no shift, gaussian
C     14 -- 2 barriers, receptor downwind of downwind_barrier, source within recirc zone; direct, shift, mixed-wake
C     15 -- 2 barriers, receptor downwind downwind_barrier, source outside of recirc zone; direct, no shift, mixed-wake
C 
C     NOTE: Recirculation zone only occurs with an upwind barrier, and only if the source is
C           located between upwind barrier and 6.5H downwind (for 1 barrier case) or 4H (for
C           2 barrier case). So, there will only be a shift if all of this is true.
C     ----------------------------------------------------------------------------------

C     Initialize local variables
      DIRECT_FLAG   = .FALSE.
      GAUSSIAN_FLAG = .FALSE.
      CONC_P  = 0.0D0
      CONC_M  = 0.0D0
      FRAN    = 0.0D0
      SZB     = 0.0D0
      SY      = 0.0D0
      SZ      = 0.0D0
      XMAX    = 0.0D0
      VERT    = 0.0D0
      HORZ    = 0.0D0
      DWDT    = 0.0D0
      ALPHA   = 1.0D0
      I_ALPHA = 1
      IA_MAX  = 2                                                ! assume HBD > HBU (will correct this below if it's not)
      HEFF    = Z  
      HBDEFF  = 0.0D0
      HSHIFT  = 0.0D0    
      SIGMAZ0_ORIG = SIGMAZ0                                     ! store intial sigmaz value from input file      
      F_UEFF  = 1.0D0                                            ! factor for lowering ueff and then slowly increasing back to original ueff

C     Compute distance from source to barrier, adjusting for upwind barrier when appropriate      
      XD     = XR_ROT - (X - XSHIFT)                             ! shift in x is applied here
      YD     = YR_ROT - Y
      
C     Setting up options for barrier combinations
      IF(NBARR .GT. 0) THEN
        HBMAX   = MAX(HBU, HBD)
        HBDEFF  = HBD * 1.25D0                                   ! effective barrier height in the mixed-wake algorithm

C       For nearly parallel winds, prevent dwDt from exceeding the distance to the end of the line source      
        XMAX    = MAX(DABS(XR_ROT - (XSBEGIN - XSHIFT)), 
     &                DABS(XR_ROT - (XSEND - XSHIFT)))
        IF(XMAX .LT. DW_PERD) DWDT = DWD
        IF(XMAX .GE. DW_PERD) DWDT = MIN(XMAX, DWD)
      END IF


C     Calculate direct plume concentration ----------------------------
     
C     Set flags for direct/meander and gaussian/mixed-wake
      IF (XD .LT. 0.0001D0) THEN     
C       NO DIRECT PLUME
        DIRECT_FLAG = .FALSE.   
      ELSE
C       DIRECT PLUME
        XD            = MAX(XD, XD_MIN)                                      ! shift in x has already occurred if needed
        DIRECT_FLAG   = .TRUE.
        GAUSSIAN_FLAG = .TRUE.
        
C       If there is a downwind barrier and receptor is downwind of it then mixed-wake mode
        IF((HBD .GT. 0.0D0) .AND. (XD .GT. DWDT)) THEN
          GAUSSIAN_FLAG = .FALSE.
        END IF

      END IF
      
C     Setting heff, szB, i_alpha, and alpha (and F_UEFF for upwind barriers)
      IF(DIRECT_FLAG) THEN
        IF(SHIFT_FLAG) THEN
C         Corresponding to case numbers 7, 11, 14 (shift cases)
          HEFF    = MAX(Z, 0.75D0 * HBU)                                     ! source height is adjusted upwards due to the upwind barrier
          SZB     = 0.0D0 
          I_ALPHA = 3                                                        ! set index for ALPHA_U
          ALPHA   = ALPHA_U                                                  ! set alpha
C                                                                            --- CALL EXPX
          F_UEFF  = 1.0D0 - 1.0D0 * EXPX(-1.0D0 * XD / (8.0D0 * HBU))        ! factor for reducing ueff in upwind barrier shift cases
        ELSE
C         Corresponding to case numbers 2, 4, 5, 9, 13, 15 (no shift cases)
          HEFF    = Z
          SZB     = 0.0D0
          I_ALPHA = 1                                                        ! set index for ALPHA = 1.0D0
          ALPHA   = 1.0D0                                                    ! set alpha
        END IF
      ELSE
C       Corresponding to case numbers 1, 3, 6, 8, 10, 12 (meander only cases)
C       Note that heff is set in MEANDER subroutine
        IF(NBARR .GT. 0) THEN
          SZB     = 0.25D0 * HBMAX
          IF(HBU .GT. HBD) IA_MAX = 3                                        ! If HBU > HBD, use ALPHA_U
          I_ALPHA = IA_MAX                                                   ! set index for MAX(ALPHA_D, ALPHA_U)
          ALPHA   = MAX(ALPHA_D, ALPHA_U)                                    ! set alpha
        ELSE
          I_ALPHA = 1
          ALPHA   = 1.0D0
        END IF
      END IF


C     Calculate vertical plume
      IF((DIRECT_FLAG) .AND. 
     &  (GAUSSIAN_FLAG)) THEN
C       Calculate with gaussian mode equations
C       Corresponding to case numbers 2, 4, 7, 9, 11, 13 (gaussian cases)      
        HSHIFT = 0.0D0
C                                                                            --- CALL EFFECTIVE_WIND
        CALL EFFECTIVE_WIND(XD, Z, HSHIFT)                                
         SZ     = SIGMAZ(XD) 
         SY     = SIGMAY(XD)
        UEFF_BU = F_UEFF * UEFF                                              ! reduce ueff immediately downwind of the barrier and then slowly increase back to ueff

C       Calculate vertical plume for gaussian mode                           --- CALL EXPX
        VERT   = RT2BYPI * (EXPX(-0.5D0 * ((HEFF - ZRECEP) / SZ)**2)
     &           + EXPX(-0.5D0 * ((HEFF + ZRECEP) / SZ)**2)) /
     &           (2.0D0 * SZ * UEFF_BU)
C      Calculate horizontal plume                                            --- CALL EXPX
        HORZ   = 1.0D0 / (SRT2PI * SY) * EXPX(-0.5D0 * (YD / SY)**2)

      ELSE IF((DIRECT_FLAG) .AND. (.NOT. GAUSSIAN_FLAG)) THEN
C       Corresponding to case numbers 5, 14, 15 (mixed-wake cases)
C       Calculate with mixed-wake mode equations in 2 steps
        
C       Step 1: calculate plume spread between source and downwind barrier   --- CALL EFFECTIVE_WIND
        HSHIFT  = HBD             
        CALL EFFECTIVE_WIND(DWDT, HEFF, HSHIFT)
        SZB     = SIGMAZ(DWDT) + 0.1D0 * HBD                                 ! assigns vertical spread of plume at dwDt to szB 
        SIGMAZ0 = 0.0D0                                                      ! set to zero to not double count intial sigmaz
        
C       Step 2: calculate plume spread beyond the downwind barrier           --- CALL EFFECTIVE_WIND
        I_ALPHA = 2                                                          ! set index for ALPHA_D      
        ALPHA   = ALPHA_D                                                    ! set alpha                                
        CALL EFFECTIVE_WIND(XD - DWDT, HEFF, HSHIFT)
        SZ      = SIGMAZ(XD - DWDT) 
        SY      = SIGMAY(XD - DWDT)
        SIGMAZ0 = SIGMAZ0_ORIG                                               ! reset sigmaz0 for next time thru subroutine

C       Calculate vertical plume for mixed-wake mode                         --- CALL EXPX
        CQ      = (1.0D0 / (SRT2PI * SZ * UEFF)) * 2.0D0 *
     &             EXPX(-0.5D0 * (HEFF / SZ)**2)                             ! Venkatram et al. (2021)
        
        HMAX    = MAX(HBD,12.0D0)                                            ! tallest barrier height affected by downwind barrier effect
        A       = (1.0D0 - (0.75D0 * HBD)/HMAX) 
C                                                                            --- CALL EXPX        
        F_UH    = 1.0D0 - A * EXPX(-(XD - DWDT) / (9.0D0 * HBD))             ! UH factor; 9H is the scaling factor controlling rise back to original ueff
          
        FQ      = 1.0D0 / (1.0D0 + F_UH * 0.5D0 * UH * CQ * HBDEFF)          ! Venkatram et al. (2021); assuming unit emissions        

C                                                                            --- CALL EXPX
        IF(ZRECEP .GT. HBDEFF) THEN
C         Note: reflection term is about Z = HBDEFF, so the reflective source is at z = HBDEFF-HEFF
          VERT  = FQ * (1.0D0 / (SRT2PI * SZ * UEFF)) * 
     &           (EXPX(-0.5D0 * ((ZRECEP - HBDEFF - HEFF) / SZ)**2) + 
     &            EXPX(-0.5D0 * ((ZRECEP - HBDEFF + HEFF) / SZ)**2))         ! Venkatram et al. (2021)
        ELSE
          VERT  = FQ * CQ                                                    ! Venkatram et al. (2021)
        END IF
C      Calculate horizontal plume                                            --- CALL EXPX
        HORZ   = 1.0D0 / (SRT2PI * SY) * EXPX(-0.5D0 * (YD / SY)**2)
        
      ELSE
C       Corresponding to case numbers 1, 3, 6, 8, 10, 12 (meander only cases)          
        HSHIFT = 0.0D0
        CALL EFFECTIVE_WIND(XD, HEFF, HSHIFT)
        VERT = 0.0D0                                                         ! no direct plume, so CONC_P should be zero
        HORZ = 0.0D0

      END IF
      
C     Calculate total direct plume
      CONC_P = VERT * HORZ 
      
C     Calculate meander concentration                                        --- CALL MEANDER
      FRAN   = 2.0D0 * SIGMAV * SIGMAV / (UEFF * UEFF)
      CONC_M = MEANDER(X, Y, Z, HSHIFT)

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
C        INPUTS:      XA, YA, X, N
C
C        OUTPUTS:     Y, ERR
C
C        CALLED FROM: NUMERICAL_LINE_SOURCE
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C        
C                    "Numerical Recipes in Fortran", Press et al., (1992)
C***********************************************************************

C     Variable Declarations:
      IMPLICIT NONE

C     Local Variables:
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

      DOUBLE PRECISION, INTENT(IN)   :: X
      DOUBLE PRECISION, INTENT(IN)   :: XA(N), YA(N)
      DOUBLE PRECISION, INTENT(OUT)  :: ERR
      DOUBLE PRECISION, INTENT(OUT)  :: Y
C     XA(N)       = table of XA values used in interpolation
C     YA(N)       = table of YA values used in interpolation
C     ERR         = error in interpolation
C     Y           = interpolated value at X
      
      INTEGER, INTENT(IN)  :: N
C     N           = length of XA and YA arrays

C     Computation Parameters:
      DOUBLE PRECISION  :: EPS = 1.0E-10
      
C     Variable Initializations:
      DELTAY = 0.0D0


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
      DO M = 1, N - 1
         DO I = 1, N - M
            HO   = XA(I) - X
            HP   = XA(I + M) - X
            W    = C(I + 1) - D(I)
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

C     Variable Declarations:
      USE MAIN1, ONLY: ISRC, NUMSRC, SRCTYP
      USE RLINE_DATA, ONLY: RLEMISCONV, RLSOURCE, RLMOVESCONV
      IMPLICIT NONE

C     Local Variables:
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
               LENGTH = DSQRT((RLSOURCE(INDEX)%XSB - 
     &                  RLSOURCE(INDEX)%XSE)**2 +
     &                  (RLSOURCE(INDEX)%YSB - 
     &                  RLSOURCE(INDEX)%YSE)**2)
               RLEMISCONV(INDEX) = 1.0d0 / LENGTH / 3600d0
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
C        INPUTS:      XD
C
C        OUTPUTS: 
C
C        CALLED FROM: POINT_CONC 
C
C        CALLING
C        ROUTINES:   
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1, ONLY: OBULEN, USTAR, STABLE
      USE RLINE_DATA, ONLY: SIGMAV, SIGMAY0, SIGZ_Y
      IMPLICIT NONE

C     Local Variables:
      DOUBLE PRECISION  :: SZ
      DOUBLE PRECISION, INTENT(IN)  :: XD
C     SZ          = vertical dispersion
C     XD          = distance between source and receptor in direction parallel to wind

C     Set sigmaz to SIGZ from SIGMAZ function before minimum taken
      SZ = SIGZ_Y

C     Calculate SIGMAY based on stability
      IF (STABLE) THEN
         SIGMAY = 1.6D0 * SIGMAV / USTAR * SZ *
     &           (1.0D0 + 2.5D0 * SZ / DABS(OBULEN))  
      ELSE
         SIGMAY = 1.6D0 * SIGMAV / USTAR * SZ /
     &            DSQRT(1.0D0 + 1.0D0 * SZ / DABS(OBULEN)) 
      END IF

      SIGMAY = DSQRT(SIGMAY**2 + SIGMAY0**2)
      
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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021
C    
C        MODIFIED:    Code integrated into AERMOD for RLINE source.
C                     Wood, 07/20/2018
C
C        INPUTS:      XD
C
C        OUTPUTS:    
C
C        CALLED FROM: EFFECTIVE_WIND
C                     MEANDER
C                     POINT_CONC
C                     SIGMAY
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1, ONLY: OBULEN, USTAR, ZI, ISRC, STABLE,
     &                 RT2BYPI, TWOTHIRDS
      USE RLINE_DATA, ONLY: RLSOURCE, UEFF, SIGMAZ0, SZB, ALPHA, 
     &                      NBARR, SIGZ_Y
      IMPLICIT NONE

C     Local Variables:
      DOUBLE PRECISION  :: SIGMAZ_MAX, XBAR, XDABS, URATIO
      DOUBLE PRECISION  :: SIGZ
      DOUBLE PRECISION  :: SIGMAZD
      DOUBLE PRECISION  :: SIGMAZB
      DOUBLE PRECISION, INTENT(IN)  :: XD
C     SIGMAZ_MAX  = maximum vertical dispersion
C     XBAR        = absolute value of XD/L
C     XDABS       = absolute value of XD
C     URATIO      = USTAR / UEFF
C     SIGZ        = vertical dispersion initial calculation
C     SIGMAZD     = vertical dispersion due to depression
C     SIGMAZB     = vertical dispersion due to barrier
C     XD          = distance between source and receptor in direction parallel to wind

      XBAR       = DABS(XD / OBULEN)
      XDABS      = DABS(XD)
      SIGMAZ_MAX = RT2BYPI * ZI
      SIGMAZD    = 0.0D0
      SIGMAZB    = 0.0D0
      URATIO     = USTAR / UEFF

C     Calculate vertical dispersion curve based on stability
      IF (STABLE) THEN
         SIGZ = ALPHA * 0.57D0 * (URATIO * XDABS) /
     &          (1.0D0 + 3.0D0 * (URATIO / ALPHA) *
     &          (EXP(TWOTHIRDS * LOG(XBAR))))
      ELSE 
         SIGZ = ALPHA * 0.57D0 * (URATIO * XDABS) *
     &          (1.0D0 + 1.5D0 * ((URATIO / ALPHA**2) * XBAR))
      END IF

C     Adjust for depressed roadway, if used
      IF (RLSOURCE(ISRC)%DEPTH < 0.0D0) THEN
         SIGMAZD = -1.0D0 * RLSOURCE(ISRC)%DEPTH / 2.15D0
      END IF
        
C     Adjust for barrier, if barrier is present
      IF(NBARR > 0) THEN
           SIGMAZB = SZB
      END IF

      SIGZ   = DSQRT(SIGZ * SIGZ + SIGMAZ0 * SIGMAZ0 +
     &         SIGMAZD * SIGMAZD + SIGMAZB * SIGMAZB) 
      SIGZ_Y = SIGZ
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
C        MODIFIED:    Incorporated RLINE2 barrier algorithms.
C                     Dianna Francisco & David Heist, US EPA ORD, 02/12/2021 
C
C        INPUTS:      
C
C        OUTPUTS:     
C
C        CALLED FROM: RLCALC
C
C        CALLING 
C        ROUTINES:    DEPRESSED_DISPLACEMENT
C
C        References: "RLINE: A Line Source Dispersion Model for Near-Surface
C                    Releases", Snyder et. al (2013)
C***********************************************************************

C     Variable Declarations:
      USE MAIN1, ONLY: AXR, AYR, NUMREC, NUMSRC, ISRC, SRCTYP, WDREF, PI
      USE RLINE_DATA, ONLY: X0, Y0, XSB_ROT, YSB_ROT, XSE_ROT, YSE_ROT,
     &                      THETAW, XRCP_ROT, YRCP_ROT, RLSOURCE,
     &                      BDW_FLAG, THETA_LINE
      IMPLICIT NONE

C     External Functions:
      DOUBLE PRECISION, EXTERNAL :: DEPRESSED_DISPLACEMENT

C     Local Variables:
      INTEGER  :: INDEX, I
      DOUBLE PRECISION :: XR_TRAN, YR_TRAN, ANGLE
      DOUBLE PRECISION :: XSB_TRAN, YSB_TRAN, XSE_TRAN, YSE_TRAN
      DOUBLE PRECISION :: DCL_LOC, DCLWALL_LOC(2)
      DOUBLE PRECISION :: XBB, YBB, XBB_ROT
C     INDEX       = index 
C     I           = index
C     XR_TRAN     = x-coordinate for translated receptor
C     YR_TRAN     = y-coordinate for translated receptor      
C     ANGLE       = 270.0 - wind direction
C     XSB_TRAN    = beginning x-coordinate for translated source
C     YSB_TRAN    = beginning y-coordinate for translated source      
C     XSE_TRAN    = end x-coordinate for translated source
C     YSE_TRAN    = end y-coordinate for translated source
C     DCL_LOC        = offset distance from center line
C     DCLWALL_LOC(2) = barrier distance from center line of road; DCLWALL(1) and DCLWALL(2)
C     XBB            = beginning coordinates of the barrier; after rotation
C     YBB            = beginning coordinates of the barrier; after rotation
C     XBB_ROT        = rotated beginning coordinates of the barrier 

C     Initialize flag for 'barrier is downwind of source' (1 if true, 0 if false)
      BDW_FLAG(:,:) = 0
                     
      ANGLE = 270.0D0 - WDREF
      IF (ANGLE > 180.0D0) THEN            
         ANGLE = ANGLE - 360.0D0
      END IF
      THETAW = ANGLE*PI / 180.0D0

C     Translate line source origin
      X0 = RLSOURCE(ISRC)%XSB    
      Y0 = RLSOURCE(ISRC)%YSB

C     Translate source and receptor coordinates and then rotate them along wind direction
      DO INDEX = ISRC, NUMSRC
         IF (SRCTYP(INDEX) .EQ. 'RLINE' .OR. 
     &       SRCTYP(INDEX) .EQ. 'RLINEXT') THEN
C           Initialize variables used
            DCL_LOC        = RLSOURCE(INDEX)%DCL 
            DCLWALL_LOC(1) = RLSOURCE(INDEX)%DCLWALL
            DCLWALL_LOC(2) = RLSOURCE(INDEX)%DCLWALL2      
  
C           Move initial user coordinate system so the origin is at the beginning of first source             
            XSB_TRAN = RLSOURCE(INDEX)%XSB - X0
            YSB_TRAN = RLSOURCE(INDEX)%YSB - Y0
            XSE_TRAN = RLSOURCE(INDEX)%XSE - X0
            YSE_TRAN = RLSOURCE(INDEX)%YSE - Y0
            THETA_LINE = DATAN2(YSE_TRAN - YSB_TRAN, XSE_TRAN - 
     &                   XSB_TRAN)

C           Account for due east and north source lines 
            IF (DSIN(THETA_LINE) .EQ. 0.0D0) THEN
               DCL_LOC        = -1.0D0 * DCL_LOC              ! needed for lines running West-East; + is North
               DCLWALL_LOC(1) = -1.0D0 * DCLWALL_LOC(1)      
               DCLWALL_LOC(2) = -1.0D0 * DCLWALL_LOC(2)      
            END IF

C           Determine location of the line that is not within a depression,
C           but is specified in source file with the centerline and distance
C           from the centerline
            IF (DCL_LOC .NE. 0.0D0 .AND. 
     &          RLSOURCE(INDEX)%DEPTH .EQ. 0.0D0) THEN
               XSE_TRAN = XSE_TRAN + DCL_LOC * DSIN(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
               YSE_TRAN = YSE_TRAN - DCL_LOC * DCOS(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
               XSB_TRAN = XSB_TRAN + DCL_LOC * DSIN(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
               YSB_TRAN = YSB_TRAN - DCL_LOC * DCOS(THETA_LINE) *
     &                    DSIGN(1.0D0, DSIN(THETA_LINE))
            END IF
    
C           Adjustments for near source configurations (depression)
            IF (RLSOURCE(INDEX)%DEPTH < 0.0D0) THEN 
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
            YSB_ROT(INDEX) = -1.0D0 * XSB_TRAN * DSIN(THETAW) +
     &                       YSB_TRAN * DCOS(THETAW)
            XSE_ROT(INDEX) = XSE_TRAN * DCOS(THETAW) +
     &                       YSE_TRAN * DSIN(THETAW)
            YSE_ROT(INDEX) = -1.0D0 * XSE_TRAN * DSIN(THETAW) +
     &                       YSE_TRAN * DCOS(THETAW)

C           For barrier algorithm: Calculate the beginning coordinates of the barrier.
C           After rotation, determine if the barrier is downwind of the roadway. If so, set a flag.
            DO I = 1, 2
               XBB     = XSB_TRAN + DCLWALL_LOC(I) * DSIN(THETA_LINE) *
     &                   DSIGN(1.0D0, DSIN(THETA_LINE))
               YBB     = YSB_TRAN - DCLWALL_LOC(I) * DCOS(THETA_LINE) *
     &                   DSIGN(1.0D0, DSIN(THETA_LINE))
               XBB_ROT = XBB * DCOS(THETAW) + YBB * DSIN(THETAW)  
               IF(XBB_ROT > XSB_ROT(INDEX)) BDW_FLAG(INDEX, I) = 1
            END DO

         END IF 
      END DO

      DO INDEX = 1, NUMREC
         XR_TRAN = AXR(INDEX) - X0
         YR_TRAN = AYR(INDEX) - Y0
         XRCP_ROT(INDEX) = XR_TRAN * DCOS(THETAW) + 
     &                     YR_TRAN * DSIN(THETAW)
         YRCP_ROT(INDEX) = -1.0D0 * XR_TRAN * DSIN(THETAW) + 
     &                     YR_TRAN * DCOS(THETAW)
      END DO
      
      END SUBROUTINE TRANSLATE_ROTATE
