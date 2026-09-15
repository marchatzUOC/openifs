MODULE TEGEN_DUST_SCHEME_MOD

CONTAINS

SUBROUTINE TEGEN_DUST_SCHEME( YDEPHY, YDEAERMAP, YDEAERSRC,                 &
                         & KIDIA, KFDIA, KLON, KLEV, KTILES, KSW,        &
                         & PLSM , PWIND, PSNS, PZ0M,                     &
                         & SP, PTL, PSOIL_TYPE,                          &
                         & PFRTI, PCVL, PCVH, KTVL, KTVH,                &
                         & EMIS_MASS, EMIS_NUMBER ,PAERFLX,PGLON, PGLAT, &
                         & PRWPWP,PRWSAT,PAERMAP,PALB,PALBD,PWS1,PHSDFOR,&
                         & IMM,ISOILPH1, ISOILPH2, ISOILPH3, ISOILPH4, ISOILPH5, &
                         & IZ0AM, IPOTSRC, ISOILTYPE, IAREA, ICULT,IZ0M, IFPAR, GPGAW,&
                         & ILAI_MAX,ILAI_AVG)

! --- IFS/OpenIFS modules ------------------------------------------------------

USE TYPE_MODEL,ONLY : MODEL
USE YOMLUN,    ONLY : NULOUT
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST,    ONLY : RPI

! -- M7 modules ----------------------------------------------------------------
USE TM5M7_DATA,      ONLY: NMOD, MODE_ACI, MODE_COI, sigma, sigma_lognormal,   &
                         & iacci,icoai
!                       
USE TM5M7_EMIS_DATA, ONLY: MODAL_EMISSIONS,                           &
                         & nsoilph, nfpar,    &
                         & vkarman!,         &

USE YOEPHY   , ONLY : TEPHY
USE YOEAERMAP, ONLY : TEAERMAP
USE YOEAERSRC, ONLY : TEAERSRC

!------------------------------------------------------------------------------!
!             0.6 ARGUMENTS TEGEN                                              !
!                                                                              !
!------------------------------------------------------------------------------!

! parameters for online dust calculations
INTEGER, PARAMETER            :: ntraced=8                     ! number of coarse-grained bins
                                                                 ! in the original emission model
INTEGER, PARAMETER            :: nbin=24                       ! number of discretization points per bin
INTEGER, PARAMETER            :: nclass=ntraced*nbin           ! total number of discretization points
INTEGER, PARAMETER            :: nats=12                       ! number of soil types
INTEGER, PARAMETER            :: nmode=4                       ! number of particle size distributions in soils,
                                                                 ! which distinguishes between clay, silt,
                                                                 ! medium/fine sand, and coarse sand
INTEGER, PARAMETER            :: nspe=nmode*3+2                ! for explanation, see below
REAL(KIND=JPRB), PARAMETER    :: xmair=28.94 ! mass of air, g/mol
REAL(KIND=JPRB), PARAMETER    :: xmdust=xmair
! Constants used in the parameterization of the efficient friction velocity ratio,
! see Eqs. (17-20) in MB95:
REAL(KIND=JPRB), PARAMETER    :: aeff=0.35
REAL(KIND=JPRB), PARAMETER    :: xeff=10.
REAL(KIND=JPRB), PARAMETER    :: u1fac=0.6    ! 0.7 in EC-Earth 3.2.3
REAL(KIND=JPRB), PARAMETER    :: ddcal=0.1   

REAL(KIND=JPRB), PARAMETER    :: cd=1.2507E-06                 ! flux dimensioning parameter [g s^2/cm^4]
REAL(KIND=JPRB), PARAMETER    :: z0_min=1.e-2
REAL(KIND=JPRB), PARAMETER    :: lai_lim=0.25
REAL(KIND=JPRB), PARAMETER    :: lai_lim2=0.5
REAL(KIND=JPRB), PARAMETER    :: d_thrsld=2.31e-6           ! threshold value
REAL(KIND=JPRB), PARAMETER    :: Dmin=2.0210403762e-5          ! diameter (cm) at first discretization point
REAL(KIND=JPRB), PARAMETER    :: Dmax=0.126667434757           ! diameter (cm) at last discretization point
REAL(KIND=JPRB), PARAMETER    :: Dstep=0.04577551202           ! diameter increment in log-space
REAL(KIND=JPRB), PARAMETER    :: grav =  9.80665               ! m/s2
! Constants in the parameterization of the Reynolds number,
! see Eq. (5) in MB95:
REAL(KIND=JPRB), PARAMETER    :: a_rnolds=1331.647             ! Reynolds constant
REAL(KIND=JPRB), PARAMETER    :: b_rnolds=0.38194              ! Reynolds constant
REAL(KIND=JPRB), PARAMETER    :: x_rnolds=1.561228             ! Reynolds constant
REAL(KIND=JPRB), PARAMETER    :: roa=0.001227                  ! reference air density (g/cm^3)

REAL(KIND=JPRB)               :: rho_air                       ! variable air density (g/cm^3)
REAL(KIND=JPRB), PARAMETER    :: rgas =8.3144
REAL(KIND=JPRB), PARAMETER    :: airfac=1./rgas*xmair*1.e-6    ! factor for rho_air
REAL(KIND=JPRB)               :: airdens_ratio, airdens_ratio2
REAL(KIND=JPRB), PARAMETER    :: umin=13.75                    ! minimum threshold friction velocity (cm/s)
REAL(KIND=JPRB), PARAMETER    :: ZZ=1000.                      ! wind measurement height (cm)
REAL(KIND=JPRB), PARAMETER    :: ddust   = 2.650              ! Density          du     [g cm-3]
REAL(KIND=JPRB), PARAMETER    :: dust_density = ddust * 1.e3

INTEGER(KIND=JPIM), PARAMETER :: min_ai=1
INTEGER(KIND=JPIM), PARAMETER :: max_ai=1
! Boundaries for Coa. mode
INTEGER(KIND=JPIM), PARAMETER :: min_ci=2
INTEGER(KIND=JPIM), PARAMETER :: max_ci=4
REAL(KIND=JPRB), PARAMETER    :: mf_acc_r1 = 0.313758
REAL(KIND=JPRB), PARAMETER    :: mf_acc_r2 = 0.684043
REAL(KIND=JPRB), PARAMETER    :: mf_coa_r1 = 0.00518309
REAL(KIND=JPRB), PARAMETER    :: mf_coa_r2 = 0.980634

REAL(KIND=JPRB), PARAMETER    :: ratio_coa = mf_coa_r1/mf_coa_r2
REAL(KIND=JPRB), PARAMETER    :: ratio_acc = mf_acc_r2/mf_acc_r1
REAL(KIND=JPRB), PARAMETER    :: denom_acc_inv = 1./(mf_acc_r1-ratio_coa*mf_acc_r2)
REAL(KIND=JPRB), PARAMETER    :: denom_coa_inv = 1./(mf_coa_r2-ratio_acc*mf_coa_r1)
REAL(KIND=JPRB), PARAMETER    :: mf_acc_r12_inv = 1./(mf_acc_r1+mf_acc_r2)
REAL(KIND=JPRB), PARAMETER    :: mf_coa_r12_inv = 1./(mf_coa_r1+mf_coa_r2)

REAL(KIND=JPRB), PARAMETER    :: mmr_ai=0.37E-4  ! cm
REAL(KIND=JPRB), PARAMETER    :: mmr_ci=1.75E-4

!----------------------------------------------------------------

!-----------------------------------------------------------------------
!*     0.1   ARGUMENTS
!            ---------
INTEGER(KIND=JPIM),     INTENT(IN)    :: IMM  ! not used
TYPE(TEPHY),           INTENT(IN)    :: YDEPHY
TYPE(TEAERMAP),        INTENT(INOUT) :: YDEAERMAP
TYPE(TEAERSRC),        INTENT(IN)    :: YDEAERSRC

INTEGER(KIND=JPIM),    INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),    INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),    INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),    INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),    INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),    INTENT(IN)    :: KSW

REAL(KIND=JPRB),       INTENT(IN)    :: GPGAW(KLON)
REAL(KIND=JPRB),       INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB),       INTENT(IN)    :: PWIND(KLON)        ! 10m wind speed, see tm5m7_src.F90
REAL(KIND=JPRB),       INTENT(IN)    :: PSNS(KLON)         ! Snow depth
REAL(KIND=JPRB),       INTENT(IN)    :: PZ0M(KLON)         ! Roughness length [m]
REAL(KIND=JPRB),       INTENT(IN)    :: SP(KLON)           ! Surface pressure
REAL(KIND=JPRB),       INTENT(IN)    :: PTL(KLON)          ! surface temperature
REAL(KIND=JPRB),       INTENT(IN)    :: PSOIL_TYPE(KLON)
REAL(KIND=JPRB),       INTENT(IN)    :: PFRTI(KLON,KTILES) ! Tile fraction (0-1)
!  1 : Water                      5 : Snow on low-veg + bare-soil 
!  2 : Ice                        6 : Dry snow-free high veg
!  3 : Wet skin                   7 : snow under high-veg
!  4 : Dry snow-free low-veg      8 : bare soil
REAL(KIND=JPRB),       INTENT(IN)    :: PCVL(KLON), PCVH(KLON) ! Low/High vegetation cover
INTEGER(KIND=JPIM),    INTENT(IN)    :: KTVL(KLON), KTVH(KLON) ! Low/High vegetation type
! M7 
TYPE(MODAL_EMISSIONS), INTENT(INOUT) :: emis_mass(NMOD)
TYPE(MODAL_EMISSIONS), INTENT(INOUT) :: emis_number(NMOD)
REAL(KIND=JPRB),       INTENT(INOUT) :: PAERFLX(KLON,12,9) !diagnostic array/not used.
REAL(KIND=JPRB),       INTENT(IN)    :: PGLON(KLON),PGLAT(KLON)
REAL(KIND=JPRB),       INTENT(INOUT) :: PRWPWP, PRWSAT, PAERMAP(KLON,5)
REAL(KIND=JPRB),       INTENT(IN)    :: PALB(KLON), PALBD(KLON,KSW)
REAL(KIND=JPRB),       INTENT(IN)    :: PWS1(KLON),PHSDFOR(KLON)

REAL(KIND=JPRB),         INTENT(IN) :: ISOILPH1(KLON), ISOILPH2(KLON), ISOILPH3(KLON), ISOILPH4(KLON), ISOILPH5(KLON), &
                                       & IZ0AM(KLON), IPOTSRC(KLON), IAREA(KLON), ICULT(KLON)
REAL(KIND=JPRB),        INTENT(IN) :: IZ0M(KLON), IFPAR(KLON)
REAL(KIND=JPRB),        INTENT(IN) :: ILAI_MAX(KLON) ,ILAI_AVG(KLON) 
REAL(KIND=JPRB),        INTENT(IN) :: ISOILTYPE(KLON)

!*    0.5   LOCAL VARIABLES
!           ---------------
REAL(KIND=JPRB)               :: exp_Dstep, sqrt_2pi
INTEGER(KIND=JPIM), PARAMETER ::  KBINDD=3 
INTEGER(KIND=JPIM) :: JL, ID, JAER, INBAER

REAL(KIND=JPRB)    :: FLUX_AI(KLON), FLUX_CI(KLON),FNUM_AI(KLON),FNUM_CI(KLON)
REAL(KIND=JPRB)    :: FLUXTOT(NTRACED),FDUST(NTRACED) 
REAL(KIND=JPRB)    :: FLUXTYP(NCLASS)
REAL(KIND=JPRB)    :: ZDEPTILE
REAL(KIND=JPRB)    :: TV_DAT(20) ! Local grid box fractions (0-1) for each of 
                                 ! presumeably 20 IFS vegetation types
! RCHG -> Here it i simportant to explain what are 9 , 12  
!         => PROBABLY related to PAERFLUX dimensions 
REAL(KIND=JPRB)    :: ZFLX_SDUST(KLON,9,12)
REAL(KIND=JPRB)    :: ZSCC2(KLON), ZDEP2(KLON) 
REAL(KIND=JPRB)    :: ZLTS2(KLON), ZLTSMIN(KLON), ZLTSMAX(KLON)
REAL(KIND=JPRB)    :: ZWND3(KLON) 
REAL(KIND=JPRB)    :: ZDUEMPOT(KLON,3)
REAL(KIND=JPRB)    :: ZDEGRAD, ZFSWET, ZSWETN
REAL(KIND=JPRB)    :: ZRWPWP, ZRWSAT 
REAL(KIND=JPRB)    :: ZEPSSNO, ZEPSARE
REAL(KIND=JPRB)    :: ZREFSPD, ZRADREF, ZREFRAD
REAL(KIND=JPRB)    :: ZAERDUB
REAL(KIND=JPRB)    :: RDDUSRC(9)
LOGICAL            :: LLDUST(KLON,12), LLPDUSTS(KLON)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
LOGICAL            :: TEGEN
CHARACTER(LEN=45)  :: CLAERWND(0:3)
!----------------------------------------------------------------
! SOIL CARACTERISTICS:
! ZOBLER texture classes:
!----------------------------------------------------------------
!! nats =12
!! nspe =nmode*3+2  = 14 
!! nmode=4
INTEGER :: jp 

!!!!>>>>>
REAL(KIND=JPRB), DIMENSION(nats,nspe) :: solspe
!--     soil type 1 : Coarse
DATA (solspe(1,jp),jp=1,nspe)/  &
     0.0707, 2.,  0.43 ,      &
     0.0158, 2.,  0.4 ,       &
     0.0015, 2.,  0.17 ,      &
     0.0002 ,2.,  0. ,        &
     2.1E-06,   0.2/
!--     soil type 2 : Medium
DATA (solspe(2,jp),jp=1,nspe)/  &
     0.0707, 2.,  0. ,            &
     0.0158, 2.,  0.37 ,          &
     0.0015, 2.,  0.33 ,          &
     0.0002, 2.,  0.3 ,           &
     4.0e-6,    0.25/
!--     soil type 3 : Fine
DATA (solspe(3,jp),jp=1,nspe)/  &
     0.0707, 2.,  0. ,            &
     0.0158, 2.,  0. ,            &
     0.0015, 2.,  0.33 ,          &
     0.0002, 2.,  0.67 ,          &
     !>>> TvN
     ! 33% x 1e-5 + 67% x 1e-7 = 3.367e-6 cm^-1
     !1.E-07,   0.5/
     3.4e-6,   0.5/
     !<<< TvN
!--     soil type 4 : Coarse Medium
DATA (solspe(4,jp),jp=1,nspe)/  &
     0.0707, 2.,  0.1 ,           &
     0.0158, 2.,  0.5 ,           &
     0.0015, 2.,  0.2 ,           &
     0.0002, 2.,  0.2 ,           &
     2.7E-06,   0.23/
!--     soil type 5 : Coarse Fine
DATA (solspe(5,jp),jp=1,nspe)/  &
     0.0707, 2.,  0. ,            &
     0.0158, 2.,  0.5 ,           &
     0.0015, 2.,  0.12 ,          &
     0.0002, 2.,  0.38 ,          &
     !>>> TvN
     ! 50% x 1e-6 + 12% x 1e-5 + 38% x 1e-6 = 2.08e-6 cm^-1
     !2.8E-06,   0.25/
     2.1e-6,   0.25/
     !<<< TvN
!--     soil type 6 : Medium Fine
DATA (solspe(6,jp),jp=1,nspe)/  &
     0.0707, 2.,  0.   ,          &
     0.0158, 2.,  0.27 ,          &
     0.0015, 2.,  0.25 ,          &
     0.0002, 2.,  0.48 ,          &
     !>>> TvN
     ! 27% x 1e-6 + 25% x 1e-5 + 48% x 1e-7 = 2.818e-6 cm^-1
     !1e-07,   0.36/
     2.8e-6,   0.36/
     !<<< TvN
!--     soil type 7 : Coarse, Medium, Fine
DATA (solspe(7,jp),jp=1,nspe)/  &
     0.0707, 2.,  0.23 ,          &
     0.0158, 2.,  0.23 ,          &
     0.0015, 2.,  0.19 ,          &
     0.0002, 2.,  0.35 ,          &
     2.5E-06,  0.25/
!--     soil type 8 : Organic
DATA (solspe(8,jp),jp=1,nspe)/  &
     0.0707, 2.,  0.25 ,          &
     0.0158, 2.,  0.25 ,          &
     0.0015, 2.,  0.25 ,          &
     0.0002, 2.,  0.25 ,          &
     0.,   0.5/
!--     soil type 9 : Ice
DATA (solspe(9,jp),jp=1,nspe)/  &
     0.0707,  2.,  0.25 ,         &
     0.0158,  2.,  0.25 ,         &
     0.0015,  2.,  0.25 ,         &
     0.0002,  2.,  0.25 ,         &
     0.,       0.5/
!--     soil type 10 : Potential Lakes (additional)
!       GENERAL CASE
DATA (solspe(10,jp),jp=1,nspe)/  &
     0.0707,  2.,  0. ,            &
     0.0158,  2.,  0. ,            &
     0.0015,  2.,  1. ,            &
     0.0002,  2.,  0. ,            &
     1.E-05,  0.25/
!--     soil type 11 : Potential Lakes (clay)
!       GENERAL CASE
DATA (solspe(11,jp),jp=1,nspe)/  &
     0.0707,  2.,  0. ,            &
     0.0158,  2.,  0. ,            &
     0.0015,  2.,  0. ,            &
     0.0002,  2.,  1. ,            &
     1.E-05,  0.25/
!--     soil type 12 : Potential Lakes Australia
DATA (solspe(12,jp),jp=1,nspe)/  &
     0.0707,  2.,  0. ,            &
     0.0158,  2.,  0. ,            &
     0.0027,  2.,  1. ,            &
     0.0002,  2.,  0. ,            &
     1.E-05,  0.25/

!!!!!!<<<<<
!------------CRITICAL ARRAYS-------------
REAL(KIND=JPRB)    :: SOIL_TYPE(KLON)
REAL(KIND=JPRB)    :: POT_SOURCE(KLON)  ! Local potencial sources are calculated 
REAL(KIND=JPRB)    :: CULT(KLON)        ! Local copy of cultivation 
REAL(KIND=JPRB)    :: Z0(KLON)          ! Local copy of roughness lengthi
REAL(KIND=JPRB)    :: FPAR(KLON)        ! Local copy of fraction photochem/radiation
REAL(KIND=JPRB)    :: SOILPH(KLON)      ! Local copy of  [THIS SHOULD BE 5 different types] 

REAL(KIND=JPRB) ::    UTH  (     NCLASS)
REAL(KIND=JPRB) ::    SREL (NATS,NCLASS)
REAL(KIND=JPRB) ::    SRELV(NATS,NCLASS)
REAL(KIND=JPRB) ::    SU_SRELV(NATS,NCLASS)

REAL(KIND=JPRB)    :: SNOWCOVER(KLON), DESERT(KLON)
REAL(KIND=JPRB)    :: LAI_EFF(KLON),UMIN2(KLON), ALPHA(KLON), C_EFF(KLON)
REAL(KIND=JPRB)    :: AREA(KLON)

INTEGER(KIND=JPIM) :: NN, ND, NS, KK, NM, NSI, NP
REAL(KIND=JPRB)    :: DP, STOTAL,STOTALV
REAL(KIND=JPRB)    :: su_class(nclass), su_classv(nclass), utest(nats)

REAL(KIND=JPRB)    :: VEGET, LAI_MAX, LAI_AVG, LAI_CUR, Z0S, DPD, FLUX_DIAM, CULTFAC1, DLAST
REAL(KIND=JPRB)    :: AAA, BB, CCC, FF, FEFF, DBSTART, UTHP, WIND10M, USTAR
REAL(KIND=JPRB)    :: XK, DDD, EE, FDP1, FDP2,temp_val
REAL(KIND=JPRB)    :: SU, SUV, SU_LOC, SU_LOCV, XL, XM, XN, XNV
REAL(KIND=JPRB)    :: FLUX_R1, FLUX_R2

REAL(KIND=JPRB) :: log_dp, log_mmd, log_stdv
REAL(KIND=JPRB), PARAMETER :: small_number = 1.0E-10

INTEGER(KIND=JPIM) :: ISTAT, REGION
INTEGER(KIND=JPIM) :: I, J, I_S1, I_S11, I_S111, IDUST, LAI_FLAG, MONTH, IVEG
INTEGER(KIND=JPIM) :: KKK, KFIRST, KKMIN
INTEGER(KIND=JPIM) :: I01, J01, I02, J02
INTEGER(KIND=JPIM) :: I1, J1, I2, J2, ACCESS_MODE
! saving the status of being called
LOGICAL, SAVE :: initial = .TRUE.
#include "abor1.intfb.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (LHOOK) CALL DR_HOOK('TEGEN_DUST_SCHEME',0,ZHOOK_HANDLE)

ASSOCIATE( NDUSRCP       => YDEAERMAP%NDUSRCP, RDDUAER => YDEAERMAP%RDDUAER,   &
         & RDUSRCP       => YDEAERMAP%RDUSRCP, NDDUST  => YDEAERSRC%NDDUST,    &
         & NALBEDOSCHEME => YDEPHY%NALBEDOSCHEME, DCAL => YDEAERSRC%DCAL,     &
         & NAERWND => YDEAERSRC%NAERWND) ! LE4ALB to NALBEDOSCHEME

CLAERWND(0) = '10-M WIND AS PREDICTOR FOR SS AND DU         '
CLAERWND(1) = 'PREDICTORS: WIND GUST FOR SS, 10M-WIND FOR DU'
CLAERWND(2) = 'PREDICTORS: WIND GUST FOR DU, 10M-WIND FOR SS'
CLAERWND(3) = 'WIND GUST AS PREDICTORS FOR SS AND DU        '

! =========================== INIT
!IF( initial ) THEN

       !---------------------------------------------------------------------------------------
       !        initializations : This should be done idealy once per day.
       !---------------------------------------------------------------------------------------
       uth      = 0.
       srel     = 0.          ! fraction of the grid area correspondent to each soil population
       srelV    = 0.          ! fraction of volume
       su_srelV = 0.
       utest    = 0.

       !---------------------------------------------------------------------------------------
       !       Uth calculation
       !       Threshold friction velocity dependent on the particle diameter
       !       following Eqs. (3-5) in MB95.
       !---------------------------------------------------------------------------------------
       nn = 0
       dp = Dmin
   !dp = MAX(Dmin, 1.0E-6)   ! avoid zero or negative dp (more robust than 1e-10)
      
   DO WHILE (dp <= Dmax + 1.0E-5)
   
       nn = nn + 1
   
       ! -----------------------------------------------------------------------------------
       ! Protect BB calculation: avoid dp ** x_rnolds if dp <= 0
       ! -----------------------------------------------------------------------------------
       IF (dp > 0.0) THEN
           BB = a_rnolds * (dp ** x_rnolds) + b_rnolds
       ELSE
           !BB = b_rnolds
           CALL ABOR1('ABORT: dp is negative in BB caclulation')
       END IF
       ! remove the check of negative  
       ! -----------------------------------------------------------------------------------
       ! XK computation - safe sqrt
       ! -----------------------------------------------------------------------------------
       IF (roa > 0.0 .AND. ddust > 0.0 .AND. grav > 0.0) THEN
           XK = SQRT(MAX(0.0, ddust * grav * 100.0 * dp / roa))
       ELSE
           CALL ABOR1('ABORT: roa or ddust or grav is Negative in XK calculation')
       END IF
   
       ! -----------------------------------------------------------------------------------
       ! CCC computation - safe sqrt
       ! -----------------------------------------------------------------------------------
       IF (dp > 0.0) THEN
           CCC = SQRT(MAX(0.0, 1.0 + d_thrsld / (dp ** 2.5)))
       ELSE
          CALL ABOR1('ABORT: dp is negative in CCC caclulation')
       END IF
   
       ! -----------------------------------------------------------------------------------
       ! Uth computation
       ! -----------------------------------------------------------------------------------
       IF (BB < 10.0) THEN
           ! Safe DDD calculation
           temp_val = 1.928 * (BB ** 0.092) - 1.0
           IF (temp_val > 0.0) THEN
               DDD = SQRT(temp_val)
               IF (DDD /= 0.0) THEN
                   Uth(nn) = 0.129 * XK * CCC / DDD
               ELSE
                   CALL ABOR1('ABORT:DDD is 0 negative Uth(nn)')
                   !Uth(nn) = 0.0 ! PRINT UTH is zeros ABORT! IDIALY NO IF ! OR SET IT TO UMIN =! ... 
               END IF
           ELSE
               CALL ABOR1('ABORT: negative BB is 0 negative Uth(nn)')
               !Uth(nn) = 0.0 !PRINT IS ZEROS ? ! ABORT 
           END IF
       ELSE
           EE = -0.0617 * (BB - 10.0)
           FF = 1.0 - 0.0858 * EXP(EE)
           Uth(nn) = 0.12 * XK * CCC * FF
       END IF
   
       ! -----------------------------------------------------------------------------------
       ! Advance dp
       ! -----------------------------------------------------------------------------------
       dp = dp * EXP(Dstep)
   
        !IF (dp > Dmax+1e.-05) CALL ABOR1("[TM5M7_SRC_DUST_INIT] NCLASS inconsistent with [Dmin,Dmax]")
   END DO 
       !THERE is BUG is the loop is not consistance with nn Uth has unidentified  variables 
       !---------------------------------------------------------------------------------------
       !       surface calculation - calculation of the soil size distribution
       !       Through all soil particle diameter the calculation of the relative contribution
       !       in surface and volume of the soil population independently of the grid
       !---------------------------------------------------------------------------------------
     ! ============================
! Soil Types Loop
! ============================
       DO ns = 1, nats ! soil types

          Stotal    = 0.
          StotalV   = 0.
          su_class  = 0.
          su_classV = 0.

          kk = 0
          dp = Dmin
          !dp = MAX(Dmin, 1.0E-10)
          DO WHILE( dp <= Dmax + 1.0E-5 )
             kk  = kk + 1
             su  = 0.
             suV = 0.
             DO nm = 1, Nmode            ! particle size populations in soils
                nd  = ((nm - 1) *3 ) + 1 ! index to mass median diameter
                nsi = nd + 1             ! index to standard deviation
                np  = nd + 2             ! index to relative contribution
                !
                !   based on soil type and contribution of population of the soil type the soil size
                !   distribution population is calculated
                !

                !>>> TvN
                ! Bug in the original code: nd should be np
                ! Since solspe(ns,nd) is never zero
                ! and the final result is proportional to solspe(ns,np),
                ! the bug has no impact on the results.
                !IF (solspe(ns,nd).EQ.0.) THEN
                !IF (solspe(ns,np).EQ.0.) THEN
                IF (solspe(ns,np).EQ.0. .or. solspe(ns,nsi).EQ.0. .or. solspe(ns,nd).EQ.0.) THEN
                 !<<< TvN
                   su_loc = 0.
                   su_locV=0.
                ELSE
                   xk      = solspe(ns,np)/(SQRT(2.* RPI)*LOG(solspe(ns,nsi)))
                   xl      = ( (LOG(dp) - LOG( solspe(ns,nd ) ))**2 ) / &
                        (2.*(LOG( solspe(ns,nsi) ))**2 )
                   xm      =  xk * EXP(-xl)         ! value of the lognormal mass size distribution
                                                    ! dM/dln(dp) in Eq. (29) in MB95
                                                    ! (Aerosol Sci. Technol., 1994)
                   xn      =  ddust*(2./3.)*(dp/2.) ! surface
                                                    ! cf. the denominator in Eq. (30) in MB95
                                                    ! The factor 2 difference is irrelevant,
                                                    ! since only relative contributions are used.
                   xnV     =  1. !volume
                   su_loc  = (xm*Dstep/xn)          ! Eq. (30) in MB95
                   su_locV = (xm*Dstep/xnV)
                END IF !
                su  = su  + su_loc
                suV = suV + su_locV
             END DO !Nmode

             su_class(kk)   = su
             su_classV(kk)  = suV
             Stotal         = Stotal + su
             StotalV        = StotalV + suV
             dp             = dp * EXP(Dstep)
          END DO !dp

          DO nn = 1,Nclass
             IF (Stotal.EQ.0.)THEN
                srel (ns,nn) = 0.
                srelV(ns,nn) = 0.
             ELSE
                srel    (ns,nn) = su_class(nn)/Stotal
                srelV   (ns,nn) = su_classV(nn)/StotalV
                utest   (ns   ) = utest(ns)+srelV(ns,nn)
                su_srelV(ns,nn) = utest(ns)
             END IF
          END DO !j=1,nclass
       END DO !ns (soil type)

    initial = .FALSE.
!END IF ! =========================== INIT

! ifs vegetation                        
!  
!1)  L ! Crops, Mixed Farming           
!2)  L ! Short Grass                    
!3)  H ! Evergreen Needleleaf Trees     
!4)  H ! Deciduous Needleleaf Trees     
!5)  H ! Deciduous Broadleaf Trees      
!6)  H ! Evergreen Broadleaf Trees      
!7)  L ! Tall Grass                     
!8)    ! Desert                         
!9)  L ! Tundra                         
!10) L ! Irrigated Crops                
!11) L ! Semidesert                     
!12)   ! Ice Caps and Glaciers
!13) L ! Bogs and Marshes               
!14)   ! Inland Water
!15)   ! Ocean
!16) L ! Evergreen Shrubs               
!17) L ! Deciduous Shrubs               
!18) H ! Mixed Forest/woodland          
!19) H ! Interrupted Forest             
!20) L ! Water and Land Mixtures        
!PAERFLX(KIDIA:KFDIA,1:12,1:9)=0._JPRB

ZFLX_SDUST(KIDIA:KFDIA,1:9,1:12)=0._JPRB
    
DUSTOPT: If (NDDUST==8) then
  
  ! Make local copy:
  uthp = 0._JPRB
  SOIL_TYPE(KIDIA:KFDIA)= ISOILTYPE(KIDIA:KFDIA)
  POT_SOURCE(KIDIA:KFDIA)= IPOTSRC(KIDIA:KFDIA)
  FPAR(KIDIA:KFDIA) = IFPAR(KIDIA:KFDIA)
  Z0(KIDIA:KFDIA) = IZ0M(KIDIA:KFDIA)
       
 ! calculation of snow cover from snow dept
 ! Tegen et al. fraction (0-1)
  snowcover(KIDIA:KFDIA) = PSNS(KIDIA:KFDIA) / 0.015
  WHERE( snowcover(KIDIA:KFDIA) > 1. ) snowcover(KIDIA:KFDIA) = 1.


  !!Maria G Checking the decimal values ? 
  !DO JL = KIDIA, KFDIA

  !  IF (SOIL_TYPE(JL) > 1.0 .AND. SOIL_TYPE(JL) < 2.0) THEN
  !      WRITE(NULOUT,*) 'Decimal soil type between 1 and 2:', SOIL_TYPE(JL), &
  !                      'at lat r =', PGLAT(JL), 'lon r =', PGLON(JL)

  !  ELSE IF (SOIL_TYPE(JL) > 2.0 .AND. SOIL_TYPE(JL) < 3.0) THEN
  !      WRITE(NULOUT,*) 'Decimal soil type between 2 and 3:', SOIL_TYPE(JL), &
  !                      'at lat r =', PGLAT(JL), 'lon r =', PGLON(JL)

  !  ELSE IF (SOIL_TYPE(JL) > 4.0 .AND. SOIL_TYPE(JL) < 5.0) THEN
  !      WRITE(NULOUT,*) 'Decimal soil type between 4 and 5:', SOIL_TYPE(JL), &
  !                      'at lat r =', PGLAT(JL), 'lon r =', PGLON(JL)

  !  END IF

  !END DO
  !
  !
  !---------------------------------------------------------------------------------------
  !       Prepare the flux calculation
  !---------------------------------------------------------------------------------------
  !
  !       Calculations done on monthly fields

  ! default: no dust source due to 
  !          - vegetation
  !          - not a desert pixel or 
  !          - no pure land grid cell
  lai_eff(KIDIA:KFDIA) = 0. 

  ! per grid box
    DO JL=KIDIA,KFDIA
         TV_DAT(:)=0. ! Fraction IFS land type in grid cell, between 0-1
         ! VH identify dominant ifs land use type.
         DO ID=1,KTILES
           ZDEPTILE=PFRTI(JL,ID)
           IF (ZDEPTILE < 0.01) CYCLE !skip if not contributing
           SELECT CASE(ID)
            CASE(1) ! Water
               TV_DAT(15)=TV_DAT(15)+ZDEPTILE
               ! TV_DAT(14)=ZDEPTILE (alternative: inland water?)
            CASE(2) ! ICE
               TV_DAT(12)=TV_DAT(12)+ZDEPTILE
            CASE(3) ! wet skin
              IF (PCVL(JL) + PCVH(JL) < 0.5) THEN
                TV_DAT(8)=TV_DAT(8)+ZDEPTILE
              ELSE 
                TV_DAT(KTVL(JL))=TV_DAT(KTVL(JL))+PCVL(JL)
                TV_DAT(KTVH(JL))=TV_DAT(KTVH(JL))+PCVH(JL)
              ENDIF
            CASE(4,5) ! Low veg, with/without snow
              TV_DAT(KTVL(JL))=TV_DAT(KTVL(JL))+ZDEPTILE ! make sure to filter out snow-events below
            CASE(6,7) ! high veg, with/without snow
              TV_DAT(KTVH(JL))=TV_DAT(KTVL(JL))+ZDEPTILE ! make sure to filter out snow-events below
            CASE(8) ! Bare soil
              TV_DAT(8)= TV_DAT(8)+ZDEPTILE
            END SELECT
         ENDDO


         !---------------------------------------------------------------------------------------
         !       Selection of potential dust sources areas
         !---------------------------------------------------------------------------------------
         !      Preferential Sources = Potential lakes
 
         !>>> TvN
         ! If monthly surface roughness is not available
         ! use the annual mean value, if available.
         ! Since the annual mean is calculated
         ! based on all available months,
         ! it has a much better spatial coverage 
         ! than the individual months.
         IF( Z0(JL) <= 0. .AND. IZ0AM(JL) > 0. ) THEN
               Z0(JL) = IZ0AM(JL) 
         !ELSE IF ( Z0(JL) <= 0. .AND. IZ0AM(JL) < 0. ) THEN 
         !      Z0(JL) = z0_min 
         ENDIF

         !<<< TvN
         IF( pot_source(JL) > 0.5 ) THEN 
            ! if the potential lake area is > 50%, it is a pot. lake grid
            SOIL_TYPE(JL) = 10.             
         !>>> TvN
         ! Use minimum value for roughness length.
         ! Since there are only few potential source areas
         ! where the annual mean is not available,
         ! this will only have a limited impact.
         !IF( z0(JL,idate(2)) <= 0. ) z0(JL,idate(2)) = 0.001 !! if z0 is not valid or missing (cm), PhD thesis Marticorena p.85
            IF( Z0(JL) <= 0. ) Z0(JL) = z0_min
            !<<< TvN
         END IF
         !---------------------------------------------------------------------------------------
         !       Calculation of the ratio: horizontal/vertical flux (alpha)
         !---------------------------------------------------------------------------------------
         !---------------------------------------------------------------------------------------
         !       Test on the vegetation type
         !---------------------------------------------------------------------------------------
         !  When cult=0, the cultivation field info is not used. Otherwise: cult(JL)=3
!!$         cult(JL)   = 0.

         desert(JL) = isoilph3(JL) + isoilph4(JL)
         !desert(JL)=TV_DAT(8)+TV_DAT(11)
         veget=0.
         veget = veget + PFRTI(JL,4)+PFRTI(JL,6)+PFRTI(JL,7) ! dry low veg + dry high veg + snow under high veg

         ! default: no dust emissions
         idust = 0 
         ! dust emissions only when 
         ! 1) there is only land (almost)
         ! 2) 'desert' is positive or vegetation active
         IF( PLSM(JL) >= 0.99 .AND. (desert(JL) > 0.001 .OR. veget > TINY(veget)) ) &
              idust = 1

         ! here is dust uptake possible
         IF( idust == 1 ) THEN
            !---------------------------------------------------------------------------------------
            !--  Calculate effective surface for fpar < lai_lim (as proxy for
            !--  veg. cover), shrubby vegetation is determined by max
            !--  annual fpar, grassy by monthly fpar (Tegen et al.2002)
            !---------------------------------------------------------------------------------------

            ! so we start with no vegetation --> full area available
            lai_eff(JL) = 1. 

            !--    get max/mean fpar of the full year --> needed for shrub land
            !lai_max = MAXVAL(ifpar(JL,1:12))
            !lai_avg =    SUM(ifpar(JL,1:12)) / 12. 
            lai_max = ILAI_MAX(JL)
            lai_avg = ILAI_AVG(JL)
            lai_cur = IFPAR(JL)


            ! ---------------------------------------------
            ! 3 classes: grass, shrub, mixed{grass,shrub}
            ! ---------------------------------------------
            ! HERE 
            ! first: grass dominated (tv(2) and tv(7))
            !        current fpar determines available area
            !VH IF( (tv_dat(iglbsfc,2)%data(JL,1) + tv_dat(iglbsfc,7)%data(JL,1)) > 50 ) THEN 
            !VH: over 50% tile fraction is low veg, with dominant veg type being agricultural land or range land: 
            IF ((TV_DAT(2) + TV_DAT(7)) > 0.5 ) THEN 

              lai_eff  (JL) = 1. - lai_cur / lai_lim

              ! second: shrub dominated (tv(16) and tv(17))
              !         if max(fpar) > 0.25 --> no dust 
              !         else max(fpar) determines area
            ELSEIF( (tv_dat(16) + tv_dat(17)) > 0.5 ) THEN 

              ! lai_eff is zero for lai_max > lai_min and 
              ! [0,1] for lai_max < lai_lim
              lai_eff  (JL) = 1. - lai_max / lai_lim

              ! third: mixtures of grass and shrub land
              !        if mean(fpar) > 0.5 --> shrub dominated --> use max(fpar) for scaling
              !        else grass dominated --> use current(fpar) for scaling
            ELSE

              IF( lai_avg > lai_lim2 ) THEN 
                lai_eff  (JL) = 1. - lai_max / lai_lim
              ELSE
                lai_eff  (JL) = 1. - lai_cur / lai_lim
              END IF

            END IF

            ! limit to valid range [0,1]
            lai_eff(JL) = MAX( 0., MIN( 1., lai_eff(JL) ) )

         END IF    ! if idust=1

         !---------------------------------------------------------------------------------------
         !     Lowering the threshold friction velocity depending on the presence of cultivations
         !---------------------------------------------------------------------------------------
         !       Factors according to dsf increase seen in data **
         !---------------------------------------------------------------------------------------
         umin2(JL) = umin
         ! 
         !---------------------------------------------------------------------------------------
         IF( icult(JL) <= 0.5 .AND. icult(JL) > 0.08 ) THEN
            IF( desert(JL) > 0. .OR. tv_dat(16) > 0.5 .OR. tv_dat(17) > 0.5 ) & 
                 umin2(JL) = umin * 0.93
            ! 
            !---------------------------------------------------------------------------------------
            IF( tv_dat(2) > 0.5 .OR. tv_dat(7) > 0.5 ) & 
                 umin2(JL) = umin * 0.99
         END IF !cult=2

         !  
         !---------------------------------------------------------------------------------------
         IF( icult(JL) > 0.5 ) THEN
            IF( ( desert(JL) > 0 ) .OR. ( tv_dat(16) > 0.5 ) .OR. ( tv_dat(17) > 0.5 ) ) &
                 umin2(JL) = umin * 0.73                 
         END IF !cult=1
         !---------------------------------------------------------------------------------------
         !       Daily z0 and efficient fraction feff
         !---------------------------------------------------------------------------------------

         i_s1 = INT( SOIL_TYPE(JL) )         ! soil type index for the calcl. of horiz. dust flux
         IF( i_s1 == 0 ) i_s1 = 9            ! set it the same as ice if the soil type is not defined
         !PAERFLX(JL,3,2)=i_s1
         ! Roughness length [cm] of the surface without obstacles, i.e. of the smooth surface:
         Z0S  = 0.001 !! en cm, these Marticorena p.85    ! optimum value for the calculation of energy loss

          
         ! Soil-type dependent saltation efficiency,
         ! i.e. the ratio between vertical and horizontal fluxes,
         ! (see  Eq. (42) in MB95; Eq. (3) in Heinold et al.):
        
         alpha(JL) = solspe(i_s1,nmode*3+1)
         !PAERFLX(JL,3,4)=alpha(JL) !=2 on land

         ! for now moist is not included but when it is done then:

         !---------------------------------------------------------------------------------------
         !       Calculation of the threshold soil moisture (w')  [Fecan, F. et al., 1999] 
         !---------------------------------------------------------------------------------------
         !          when moist is included   !!!!!!!!!!!!!!!!!!
         !          w_str(j,i,1) = 0.0014*(solspe(i_s1,nmode*3)*100)**2 + 0.17*(solspe(i_s1,nmode*3)*100)
         !          W0   = 0.99           ! used by Bernd solspe(i_s1,nmode*3+2)
         feff = 0.
         !          * partition of energy between the surface and the elements of rugosity *
         !           these pp 111-112

         IF( Z0(JL) <= 0. ) THEN     ! if there are no info on z0 and no potential sources
            Z0(JL) = 1.             ! then z0 is set to 1 and no dust can be produced
            feff = 0.
         ELSE
            !>>> TvN
            ! Use minimum value for roughness length.
	    ! VH convert PZ0M from [m] to [cm]
            !z0(JL) = z0_min !max(z0_min,PZ0M(JL)*100._JPRB )
            Z0(JL) = max(z0_min,Z0(JL))
            !write(3000,*)z0(JL),z0_min
            !<<< TvN
            ! Eq. (20) in MB95:
            AAA = LOG( z0(JL) / Z0S )
            BB  = LOG( aeff * (xeff / Z0S)**0.8)
            !write(5547,*)aeff,xeff,z0s
            CCC = 1. - AAA/BB
            !          * partition between Z01 and Z02 * which are z0 of larger stone which cannot be mobilized
            FF = 1.    ! we do not separate roughness length between soil which
                       ! gives dust and solid material which is not mobilised
            ! total efficient friction velocity ratio:
            feff = FF * CCC
            !!PAERFLX(JL,6,1) = feff
            !!PAERFLX(JL,6,2) = AAA
            !!PAERFLX(JL,6,3) = BB 
            !!PAERFLX(JL,6,4) = CCC
            ! restrict to [0,1]
            feff = MIN( 1., feff )
            feff = MAX( 0., feff )
         END IF

         c_eff(JL) = feff  ! scaling parameter for the threshold friction velocity

         ! due to energy loss
         !---------------------------------------------------------------------------------------
       END DO     ! JL
       !---------------------------------------------------------------------------------------
       !      End of daily base calculations

!VH    END IF ! newday 

    ! reset flux masses 
    flux_ai(KIDIA:KFDIA) = 0. 
    flux_ci(KIDIA:KFDIA) = 0. 

    DO JL = KIDIA,KFDIA

      !-- initialisation of the fields
      !   size: ntraced
      fluxtot = 0. 
      fdust   = 0.


      !----- --------------------------------------------------------------------------
      !     Calculation of dust emission flux
      !     dependent on the 3 hourly wind fields
      !----------------------------------------------------------------------
      IF( c_eff(JL) > 0. ) THEN 

         ! Calculation of ustar

         ! AS: initialise ustar (for those cases where if statement(s) are not fulfilled)
         ustar = 0. 

         IF( PLSM(JL) > 0. ) THEN 
            ! wind10m = SQRT(u10m_dat(iglbsfc)%data(JL,1)**2 + &
            !                v10m_dat(iglbsfc)%data(JL,1)**2) * 100. ! cm/s
            ustar = (vKarman * PWIND(JL)*100._JPRB) / ( log( ZZ / Z0(JL) ) ) ! cm/s
         ENDIF

         IF( Ustar > 0 .AND. (Ustar > umin2(JL) / c_eff(JL)) ) THEN

            !>>> TvN 
            rho_air = SP(JL)/PTL(JL)*airfac ! g/cm3
            airdens_ratio  = rho_air/roa
            airdens_ratio2 = sqrt(roa/rho_air)
            !<<< TvN

            !-- initialisation of the fields
            !   size: ntraced
            !dbmin   = 0. 
            !dbmax   = 0. 
            !    size: nclass
            fluxtyp = 0.


            ! soil type index for the calcl. of horiz. dust flux
            i_s1 = INT( SOIL_TYPE(JL) )            
            ! set it the same as ice
            IF( i_s1 == 0 ) i_s1 = 9            
            ! to separate from now on between saltation and mobilisation
            i_s11 = i_s1                  
            ! to separate between mobilisation and saltation and dust particles
            IF( i_s1 == 10 .OR. i_s1 == 12 ) i_s11 = 11 
            kk = 0
            dp = Dmin
            DO WHILE( dp <= Dmax+1E-5)
               kk    = kk+1
               uthp  = uth(kk) * umin2(JL) / umin * u1fac !reduce saltation threshold for cultivated soils
               !>>> TvN
               ! Include correction factor for variable air density
               uthp = uthp * airdens_ratio2
               !<<< TvN

               ! See Eq. (28) in MB95; Eq. (6) in Tegen et al.; Eq. (2) in Heinold et al.
               ! Note that (1+R)^2 * (1-R) = (1+R) * (1-R^2)
               fdp1 = (1.-(Uthp/(c_eff(JL) * Ustar)))           ! component of the horiz. flux
               fdp2 = (1.+(Uthp/(c_eff(JL) * Ustar)))**2.      !    

               IF( fdp1 > 0 .AND. fdp2 > 0) THEN

                  ! vertical flux dust weighted by the surface area relative to each soil type
                  flux_diam = srel(i_s1,kk) * fdp1 * fdp2 * cd * Ustar**3 * alpha(JL)
                  !>>> TvN
                  ! Include correction factor for variable air density
                  flux_diam = flux_diam * airdens_ratio
                  !<<< TvN

                  !----------------------------------------------------------------------
                  !   all particles even the small ones can be mobilised by saltation
                  !----------------------------------------------------------------------
                  dbstart = dmin

                  IF( dbstart >= dp ) THEN 
                     fluxtyp(kk) = fluxtyp(kk) + flux_diam
                  ELSE
                     !----------------------------------------------------------------------
                     !  loop over dislocated dust particle sizes
                     !----------------------------------------------------------------------
                     dpd    = dmin
                     kkk    = 0
                     kfirst = 0
                     DO WHILE( dpd <= dp+1e-5 )
                        kkk = kkk + 1
                        IF( dpd >= dbstart ) THEN                      ! the particles produced by saltation are put
                           IF( kfirst == 0 ) kkmin = kkk               ! in finer bins
                           kfirst = 1
                           !----------------------------------------------------------------------
                           !  scaling with relative contribution of dust size  fraction
                           !  we take into account the volume contribution of the particle types:
                           !  all the particles from soil type 10 are put into the 11 soil type when
                           !  we are in the production region
                           !----------------------------------------------------------------------
                           IF( kk > kkmin ) THEN
                             ! remember: i_s11 puts the mobilised
                             fluxtyp(kkk) = fluxtyp(kkk) + flux_diam * srelV(i_s11,kkk) / &
                             (su_srelV(i_s11,kk) - su_srelV(i_s11,kkmin) )
                             ! particles in smaller bins
                           END IF !kk.gt.kmin
                        END IF !dpd.gt.dbstart
                        dpd = dpd * EXP(dstep)
                     END DO !dpd
                     !----------------------------------------------------------------------
                     !  end of saltation loop
                     !----------------------------------------------------------------------
                  END IF !dbstart.lt.dp

               END IF !fdp1

               dp = dp * EXP(Dstep)

            END DO !dp   
            !----------------------------------------------------------------------
            !  assign fluxes to bins: flux is in g cm-2 s-1 for each bin
            !  192 sub-bins are put into 8 bins
            !----------------------------------------------------------------------
            dp    = dmin   
            dlast = dmin
            nn    = 1
            kk    = 0
            DO WHILE( dp <= dmax+1e-5 )  
               kk = kk+1
               ! add to total
               IF( nn <= ntraced ) fluxtot(nn) = fluxtot(nn) + fluxtyp(kk) 

               IF( MOD(kk,nbin) == 0 ) THEN
                   !dbmax(nn) = dp * 10000. * 0.5  !radius in um
                   !dbmin(nn) = dlast * 10000. * 0.5
                   !dpk(nn)   = SQRT( dbmax(nn) * dbmin(nn) )
                   nn        = nn+1
                   dlast     = dp
               END IF

               dp = dp * EXP(Dstep)

            END DO !dp      

         END IF   !ustar
      END IF   !c_eff 

      ! Masking the area covered by snow, vegetation and [...?...]
      cultfac1 = 1.

      DO nn = 1, ntraced
         !        fluxtot: g/cm2/sec 
         !    MASK: Effective area determined by cultfac1/snow
         fdust(nn) = fluxtot(nn) * cultfac1 * (1.-snowcover(JL))

         !    MASK: Effective area determined by fpar:

         fdust(nn) = fdust(nn) * lai_eff(JL) ! turn off vegetation limitation here!
         ! TvN: an alternative approach based on surface roughness
         ! is applied by Laurent et al. (JGR, 2006).


         !    MASK: Soil moisture threshold, using w0
         !        when moisture is included    !!!!!!!!!!!!!!!!!!
         !           IF(qrsur(JL).GE.w0) THEN
         !         fdust(JL,nn)=0.
         !           END IF
       !soil size distribution
       !       Through all soil particle diameter the calculation of the relative contribution
       !       in surface and volume of the soil population independently of the grid
      END DO
      ! ------------------------------------------------------------------------------
      ! Grouping into 2 modes: 1sec accumulation
      !
      !>>> TvN
      !   Accumulation
      flux_r1 = 0.
      DO nn = min_ai, max_ai
       !flux_ai(JL) = flux_ai(JL) + fdust(nn)
       flux_r1 = flux_r1 + fdust(nn)
      END DO

      !   Coarse
      flux_r2 = 0.
      DO nn = min_ci, max_ci
         !flux_ci(JL) = flux_ci(JL) + fdust(nn)
         flux_r2 = flux_r2 + fdust(nn)
      END DO

      ! The solution of the system of linear equations
      ! (see comments above).
      ! For special conditions, 
      ! the solution can give a negative mass flux 
      ! in either the accumulation or coarse mode.
      ! In those case, all mass is put into
      ! the other mode.
     
      !units : 
      !
     
      flux_ai(JL) = flux_r1 - ratio_coa * flux_r2
      flux_ci(JL) = flux_r2 - ratio_acc * flux_r1
      IF (flux_ai(JL) .gt. 0. .AND. flux_ci(JL) .gt. 0.) THEN
        flux_ai(JL) = flux_ai(JL) * denom_acc_inv
        flux_ci(JL) = flux_ci(JL) * denom_coa_inv
      ELSEIF (flux_ai(JL) .lt. 0.) THEN
        flux_ai(JL) = 0.
        flux_ci(JL) = (flux_r1 + flux_r2) * mf_coa_r12_inv
      ELSEIF (flux_ci(JL) .lt. 0.) THEN
        flux_ai(JL) = (flux_r1 + flux_r2) * mf_acc_r12_inv
        flux_ci(JL) = 0.
      ENDIF
      !<<< TvN

      ! now scale the emissions
      ! convert from g/cm2/s to  g/m2/s to kg/m2/s
      flux_ai(JL) = flux_ai(JL) *  1.E01  *DCAL
      flux_ci(JL) = flux_ci(JL) *  1.E01  *DCAL

      !----------------------------------------------------------
      ! Ensure fluxes are non-negative
      IF (flux_ai(JL) < 0.0_JPRB) THEN
        flux_ai(JL) = 0.0_JPRB
      END IF
      IF (flux_ci(JL) < 0.0_JPRB) THEN
        flux_ci(JL) = 0.0_JPRB
      END IF


      !-------------------------------------------------------------------------------
      !  Calculating number flux (#/m2/sec)
      !  (kg/m2/s) / g/cm3 *cm3 = (kg/s)/g *1e3 = #/m2/s =   
      !   Accumulation
      fnum_ai(JL) =  1e3*flux_ai(JL) * 3. / (4.*RPI*ddust*mmr_ai**3) * EXP(4.5*LOG(sigma(iacci))**2)
      !   Coarse
      fnum_ci(JL) = 1e3*flux_ci(JL) * 3. / (4.*RPI*ddust*mmr_ci**3) * EXP(4.5*LOG(sigma(icoai))**2)
      
      IF (fnum_ai(JL) < 0.0_JPRB) THEN
        fnum_ai(JL)= 0.0_JPRB
      END IF
      IF (fnum_ci(JL) < 0.0_JPRB) THEN
        fnum_ci(JL)= 0.0_JPRB
      END IF
      ! ------------------------------
      ! accumulation mode
      ! number
      emis_number(mode_aci)%d3(JL,KLEV,1)   =   emis_number(mode_aci)%d3(JL,KLEV,1) +fnum_ai(JL) 

      ! mass
      emis_mass(mode_aci)%d3(JL,KLEV,1)   = emis_mass(mode_aci)%d3(JL,KLEV,1) +flux_ai(JL)
      ! ------------------------------
      ! coarse mode
      ! number
      emis_number(mode_coi)%d3(JL,KLEV,1)   = emis_number(mode_coi)%d3(JL,KLEV,1) + fnum_ci(JL)
      ! mass
      emis_mass(mode_coi)%d3(JL,KLEV,1)   = emis_mass(mode_coi)%d3(JL,KLEV,1)+flux_ci(JL)
    ENDDO

END IF DUSTOPT
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TEGEN_DUST_SCHEME',1,ZHOOK_HANDLE)
END SUBROUTINE TEGEN_DUST_SCHEME

END MODULE TEGEN_DUST_SCHEME_MOD
