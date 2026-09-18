SUBROUTINE ECMWF_DUST_SCHEME( YDEPHY, YDEAERMAP, YDEAERSRC,                 &
                         & KIDIA, KFDIA, KLON, KLEV, KTILES, KSW,        &
                         & PLSM , PWIND, PSNS, PZ0M,                     &
                         & SP, PTL, PSOIL_TYPE,                          &
                         & PFRTI, PCVL, PCVH, KTVL, KTVH,                &
                         & EMIS_MASS, EMIS_NUMBER ,PAERFLX,PGLON, PGLAT, &
                         & PRWPWP,PRWSAT,PAERMAP,PALB,PALBD,PWS1,PHSDFOR)

! RCHG -> Here a dependence is KLEV => it is that ok?

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                      (updated 04-Jun-2024) │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │  *tm5m7_src_dust* - SOURCE TERMS FOR MINERAL DUST AEROSOLS                 │
! │                                                                            │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *tm5m7_src_dust* is called from tm5m7_src                                │
! │                                                                            │
! │                                                                            │
! │ Input :                                                                    │
! │ -----                                                                      │
! │                                                                            │
! │                                                                            │
! │ Output :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │                                                                            │
! │ Externals :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Method :                                                                   │
! │ ------                                                                     │
! │  Online dust emissions based on Tegen/Vignati/Strunk                       │
! │                                                                            │
! │  Please read the section above for background information about the        │
! │  underlying approach. An improved and modified online implementation has   │
! │  been accomplished from which. It can be activated by setting              │
! │                                                                            │
! │    input.emis.dust : ONLINE                                                │
! │                                                                            │
! │  in the rc-file. An additional netcdf file is needed for some input        │
! │  parameters. The path to which needs to be defined in the key              │
! │                                                                            │
! │    input.emis.dust.dir :                                                   │
! │    /ms_perm/TM/TM5/emissions/other/Dust_online/onlinedust.nc               │
! │                                                                            │
! │  For every time step there will be particles emitted, scaled to monthly    │
! │  amounts (both mass and numbers) in order to keep compliance with          │
! │  assumption sabout the aerosol emissions in sedimentation.F90.             │
! │                                                                            │
! │ Reference :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Author :                                                                   │
! │ -------                                                                    │
! │     Orginal version: T. van Noije et al. (KNMI)                            │ 
! │     Nov 2011 - Achim Strunk - v0                                           │
! │     Vincent Huijen (KNMI) adapted to OpenIFS                               │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │     Jun.  2024 - R. Checa-Garcia: revision for CY48r1 and refactory        │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯


! --- IFS/OpenIFS modules ------------------------------------------------------

USE TYPE_MODEL,ONLY : MODEL
USE YOMLUN,    ONLY : NULOUT
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST, ONLY : RPI

! -- M7 modules ----------------------------------------------------------------
USE TM5M7_DATA,      ONLY: NMOD, MODE_ACI, MODE_COI, sigma, sigma_lognormal,   &
                         & ddust,iacci,icoai
USE TM5M7_EMIS_DATA, ONLY: MODAL_EMISSIONS, NTRACED,                           &
                         & nclass, nmode,nbin, nats,solspe, nsoilph, nfpar,    &
                         & umin, z0_min, lai_lim, lai_lim2, ROA,               &
                         & U1FAC, cd,  min_ai, max_ai, min_ci, max_ci,         &
                         & vkarman, zz, airfac, DMIN, DMAX, DSTEP,             &
                         & UTH, SREL, SRELV, SU_SRELV,                         &
                         & ratio_coa, ratio_acc, denom_acc_inv, denom_coa_inv, &
                         & mf_coa_r12_inv, mf_acc_r12_inv,                     &
                         & mmr_ai, mmr_ci
                       
USE YOEPHY   , ONLY : TEPHY
USE YOEAERMAP, ONLY : TEAERMAP
USE YOEAERSRC, ONLY : TEAERSRC


!-----------------------------------------------------------------------
!*     0.1   ARGUMENTS
!            ---------

TYPE(TEPHY),           INTENT(IN)    :: YDEPHY
TYPE(TEAERMAP),        INTENT(INOUT) :: YDEAERMAP
TYPE(TEAERSRC),        INTENT(IN)    :: YDEAERSRC

INTEGER(KIND=JPIM),    INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),    INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),    INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),    INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),    INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),    INTENT(IN)    :: KSW

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
REAL(KIND=JPRB),       INTENT(INOUT) :: PAERFLX(KLON,12,9)
REAL(KIND=JPRB),       INTENT(IN)    :: PGLON(KLON),PGLAT(KLON)
REAL(KIND=JPRB),       INTENT(INOUT) :: PRWPWP, PRWSAT, PAERMAP(KLON,5)
REAL(KIND=JPRB),       INTENT(IN)    :: PALB(KLON), PALBD(KLON,KSW)
REAL(KIND=JPRB),       INTENT(IN)    :: PWS1(KLON),PHSDFOR(KLON)

!*    0.5   LOCAL VARIABLES
!           ---------------

INTEGER(KIND=JPIM), PARAMETER ::  KBINDD=3 
INTEGER(KIND=JPIM) :: JL, I_S1,I_S11, ID, JAER, INBAER
INTEGER(KIND=JPIM) :: IDUST,KK,KKK,KFIRST,KKMIN,NN
REAL(KIND=JPRB)    :: FDP1,FDP2
REAL(KIND=JPRB)    :: VEGET, z0s, dpd, uthp, flux_diam, cultfac1
REAL(KIND=JPRB)    :: AAA,BB, CCC, AEFF,FF, XEFF, feff, DBSTART, USTAR, rho_air
REAL(KIND=JPRB)    :: DLAST,DP
REAL(KIND=JPRB)    :: AIRDENS_RATIO,AIRDENS_RATIO2
REAL(KIND=JPRB)    :: FLUX_R1,FLUX_R2
REAL(KIND=JPRB)    :: SNOWCOVER(KLON), DESERT(KLON)
REAL(KIND=JPRB)    :: LAI_EFF(KLON),UMIN2(KLON), ALPHA(KLON), C_EFF(KLON)
REAL(KIND=JPRB)    :: Z0(KLON)        ! Local copy of roughness length
REAL(KIND=JPRB)    :: SOIL_TYPE(KLON) ! Local copy of soil type

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

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ECMWF_DUST_SCHEME',0,ZHOOK_HANDLE)

ASSOCIATE( NDUSRCP       => YDEAERMAP%NDUSRCP, RDDUAER => YDEAERMAP%RDDUAER,   &
         & RDUSRCP       => YDEAERMAP%RDUSRCP, NDDUST  => YDEAERSRC%NDDUST,    &
         & NALBEDOSCHEME => YDEPHY%NALBEDOSCHEME) ! LE4ALB to NALBEDOSCHEME

ZFLX_SDUST(KIDIA:KFDIA,1:9,1:12)=0._JPRB
IF (NDDUST == 3 ) THEN ! case ECMWF formulation
    !end if
    
!!!!! AER emissions
!!$CALL SURF_INQ(YSURF,PRWPWP=ZRWPWP)
!!$CALL SURF_INQ(YSURF,PRWSAT=ZRWSAT)

!!$ZETAH(KIDIA:KFDIA,0)=0._JPRB
!!$DO JK=1,KLEV
!!$  DO JL=KIDIA,KFDIA
!!$    ZETA(JL,JK) =PAP(JL,JK) /PAPH(JL,KLEV)
!!$    ZETAH(JL,JK)=PAPH(JL,JK)/PAPH(JL,KLEV)
!!$  ENDDO
!!$ENDDO
 
!ITEST=0
!ZDEGRAD= 180._JPRB/RPI
!ZDLAT  = 180._JPRB / NDGLG      ! distance in degrees between latitude lines
!ZGRDLAT= RPI / NDGLG                          ! distance in radians between latitude lines
!ZGRDLAT2=ZGRDLAT*0.55_JPRB

!ZDDUAER(:) = 1.00_JPRB
!KBINDD=3

RDDUAER(:) = 0.0_JPRB
RDDUSRC(:)= 0.0_JPRB
NDUSRCP(:) = 1
RDUSRCP(:,:) = 0.0_JPRB

!* Default values for RDDUAER
RDDUAER(1)=1.0_JPRB
RDDUAER(2)=1.0_JPRB
RDDUAER(3)=0.5_JPRB
RDDUAER(4)=0.6_JPRB
RDDUAER(5)=0.6_JPRB
RDDUAER(6)=1.0_JPRB
RDDUAER(7)=1.0_JPRB
RDDUAER(8)=1.0_JPRB
RDDUAER(9)=1.0_JPRB
RDDUAER(10)=1.0_JPRB
RDDUAER(11)=1.0_JPRB
RDDUAER(12)=1.0_JPRB
RDDUAER(13)=0.0_JPRB
RDDUAER(14)=0.2_JPRB
RDDUAER(15)=0.5_JPRB
RDDUAER(16)=1.0_JPRB
RDDUAER(17)=1.0_JPRB
RDDUAER(18)=0.5_JPRB
RDDUAER(19)=0.5_JPRB
RDDUAER(20)=0.8_JPRB
RDDUAER(21)=1.2_JPRB
RDDUAER(22)=1.5_JPRB
RDDUAER(23)=1.5_JPRB
RDDUAER(24)=0.5_JPRB
RDDUAER(25)=1.0_JPRB
RDDUAER(26)=0.5_JPRB
RDDUAER(27)=1.0_JPRB
RDDUAER(28)=0.5_JPRB
RDDUAER(29)=1.0_JPRB
RDDUAER(30)=1.0_JPRB
RDDUAER(31)=0.1_JPRB
RDDUAER(32)=0.1_JPRB
RDDUAER(33)=0.3_JPRB
RDDUAER(34)=0.8_JPRB
RDDUAER(35)=0.5_JPRB
RDDUAER(36)=0.4_JPRB
RDDUAER(37)=0.6_JPRB

!* default reference values for threshold speed and reference particle radius
! -- 1 N & S America, Europe
RDUSRCP(1,1) = 6.0_JPRB
RDUSRCP(1,2) = 5.0_JPRB
! -- 2 Russia, Urals
RDUSRCP(2,1) = 6.0_JPRB
RDUSRCP(2,2) = 5.0_JPRB
! -- 3  Africa, Sahara, S. Africa
RDUSRCP(3,1) = 6.0_JPRB
RDUSRCP(3,2) = 5.0_JPRB
! -- 4 Australasia
RDUSRCP(4,1) = 4.0_JPRB
RDUSRCP(4,2) = 5.0_JPRB
! -- 5 Asian deserts
RDUSRCP(5,1) = 3.5_JPRB
RDUSRCP(5,2) = 5.0_JPRB
! -- 6 dry lands of S.America
RDUSRCP(6,1) = 4.0_JPRB
RDUSRCP(6,2) = 5.0_JPRB
! -- 7 the rest (Japan, Greenland, Antarctica)
RDUSRCP(7,1) = 4.0_JPRB
RDUSRCP(7,2) = 5.0_JPRB
RDDUSRC(1)=0.3_JPRB
RDDUSRC(2)=0.8_JPRB
RDDUSRC(3)=5.5_JPRB
!-- default values are for use of 10-m wind as predictor for SS and DU
!RFCTDU     = 1.0_JPRB
!RFCTSS     = 1.0_JPRB  
!* New defaults taken from previous namelist; stj - 27-10-2010
!RFCTDUR    = 1.0_JPRB
!RFCTSSR    = 1.0_JPRB
!!$CLAERWND(0) = '10-M WIND AS PREDICTOR FOR SS AND DU         '
!!$CLAERWND(1) = 'PREDICTORS: WIND GUST FOR SS, 10M-WIND FOR DU'
!!$CLAERWND(2) = 'PREDICTORS: WIND GUST FOR DU, 10M-WIND FOR SS'
!!$CLAERWND(3) = 'WIND GUST AS PREDICTORS FOR SS AND DU        '


    !*       0.6   EMPIRICAL EFFICIENCY FACTORS FOR SOURCES
!              ----------------------------------------
! N.B.: Security parameters
ZEPSISS=1.E-09_JPRB
ZEPSIDD=1.E-12_JPRB
ZEPSIRA=1.E-06_JPRB
ZEPSISS=0.E+00_JPRB
ZEPSIDD=0.E+00_JPRB
ZEPSIRA=0.E+00_JPRB
ZEPSSNO=1.E-03_JPRB
ZEPSARE=1.E-03_JPRB

!PAERMAP(KIDIA:KFDIA,1:5) = 0._JPRB
DO JL=KIDIA,KFDIA
  ZLAT=PGLAT(JL)
  ZLON=PGLON(JL) 
  ZBNDA= 30._JPRB+(36._JPRB -ZLAT)*14._JPRB/24._JPRB
  ZBNDB= 30._JPRB+(36._JPRB -ZLAT)*40._JPRB/16._JPRB
  ZBNDC= 38._JPRB+(ZLON-124._JPRB)*12._JPRB/29._JPRB
  ZBNDD= 32._JPRB-(ZLON-243._JPRB)* 6._JPRB/21._JPRB

!-- Eastern border Canada/USA
  ZBNDE= 49._JPRB
  IF (ZLON > 268._JPRB .AND. ZLON < 277._JPRB) THEN
    ZBNDE= 49._JPRB-(ZLON-268._JPRB)*7._JPRB/9._JPRB
  ELSEIF (ZLON >= 277._JPRB .AND. ZLON < 285._JPRB) THEN
    ZBNDE= 42._JPRB+(ZLON-277._JPRB)*2._JPRB/8._JPRB
  ELSEIF (ZLON >= 285._JPRB .AND. ZLON < 310._JPRB) THEN
    ZBNDE= 44._JPRB+(ZLON-285._JPRB)*3._JPRB/25._JPRB
  ENDIF

!-- limits Britain
  ZLONGB=-9999._JPRB
  IF (ZLON > 354._JPRB .AND. ZLON < 360._JPRB) THEN
    ZLONGB=ZLON
  ELSEIF (ZLON >= 0._JPRB .AND. ZLON < 3._JPRB) THEN
    ZLONGB=ZLON+360._JPRB
  ENDIF
  ZBNDF= 47._JPRB+(ZLONGB-349._JPRB)*4.5_JPRB/14._JPRB

!-- limits Ireland
  ZBNDG= 61._JPRB-(ZLON-349._JPRB)*7._JPRB/6._JPRB
  ZBNDH= 45._JPRB+(ZLON-349._JPRB)*9._JPRB/6._JPRB

!-- Western border Brazil
  IF (ZLAT <= 4._JPRB .AND. ZLAT > 2._JPRB) THEN
    ZBNDI= 296._JPRB
  ELSEIF (ZLAT <= 2._JPRB .AND. ZLAT > -4._JPRB) THEN
    ZBNDI= 290._JPRB
  ELSEIF (ZLAT <= -4._JPRB .AND. ZLAT > -7._JPRB) THEN
    ZBNDI= 290._JPRB-(-4._JPRB-ZLAT)*4._JPRB/3._JPRB
  ELSEIF (ZLAT <= -7._JPRB .AND. ZLAT > -11._JPRB) THEN
    ZBNDI= 286._JPRB+(-7._JPRB-ZLAT)*4._JPRB/4._JPRB
  ENDIF

  IF (ZLAT <= -11._JPRB .AND. ZLAT > -18._JPRB) THEN
    ZBNDJ= 294._JPRB+(-11._JPRB-ZLAT)*8._JPRB/7._JPRB
  ELSEIF (ZLAT <= -18._JPRB .AND. ZLAT > -27._JPRB) THEN
    ZBNDJ= 302._JPRB+(-18._JPRB-ZLAT)*4._JPRB/9._JPRB
  ELSEIF (ZLAT <= -27._JPRB .AND. ZLAT > -30._JPRB) THEN
    ZBNDJ= 306._JPRB-(-27._JPRB-ZLAT)*3._JPRB/3._JPRB
  ELSEIF (ZLAT <= -30._JPRB .AND. ZLAT >= -34._JPRB) THEN
    ZBNDJ= 303._JPRB+(-30._JPRB-ZLAT)*4._JPRB/4._JPRB
  ENDIF

!-- Northern border India
  IF (ZLON > 70._JPRB .AND. ZLON <= 90._JPRB) THEN
    ZBNDK= 35._JPRB-(ZLON-70._JPRB)*0.5_JPRB
  ENDIF

!-- South border of Asian deserts
  IF (ZLON > 90._JPRB .AND. ZLON <= 135._JPRB) THEN
    ZBNDL= 25._JPRB+(ZLON-90._JPRB)*15._JPRB/45._JPRB
  ENDIF

!-- North limit of the Argentinian pampas
  IF (ZLON > 285._JPRB .AND. ZLON <= 297._JPRB) THEN
    ZBNDM= -42._JPRB+(ZLON-285._JPRB)*6._JPRB/12._JPRB 
  ENDIF

  IFF=0
  ITYPDU=0
 
!-- North America
!  ITYPDU=1
!----- Canada
  IF ( ZLAT >= ZBNDE .AND.&
    &      (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=1
!----- USA
  ELSEIF ( (ZLAT >= ZBNDD .AND. ZLAT < ZBNDE )&
    & .AND. (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=3
  ENDIF
!-- Alaska
  IF ( (ZLAT < 72._JPRB .AND. ZLAT > 52._JPRB)&
    & .AND. (ZLON > 190._JPRB .AND. ZLON <= 219._JPRB) ) THEN
    IFF=2
  ENDIF

!-- Central America
  IF (ZLAT < ZBNDD .AND.&
    &     (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=4
  ENDIF

!-- South America
  IF ( ZLAT < 12._JPRB .AND.&
    &      (ZLON > 190._JPRB .AND. ZLON < 330._JPRB) ) THEN
    IFF=5
  ENDIF
!-- Brazil
  IF ( (ZLAT <= 4._JPRB .AND. ZLAT > 2._JPRB)&
    &      .AND. (ZLON >= 296._JPRB .AND. ZLON <= 300._JPRB) ) THEN
    IFF=6 
  ENDIF
  IF (ZLAT <= 2._JPRB .AND. ZLAT > -11._JPRB) THEN
    IF (ZLON >= ZBNDI .AND. ZLON < 330._JPRB) THEN
      IFF=6
    ENDIF
  ENDIF
  IF (ZLAT <= -11._JPRB .AND. ZLAT >= -34._JPRB) THEN
    IF (ZLON >= ZBNDJ .AND. ZLON < 330._JPRB) THEN
      IFF=6
    ENDIF
  ENDIF

!-- Western Europe
  IF ( ZLAT > 36._JPRB .AND. ( ZLON >= 330._JPRB .OR. ZLON <= 30._JPRB) ) THEN
    IFF=10
  ENDIF

!----- Iceland
  IF ( (ZLAT < 67._JPRB .AND. ZLAT > 63._JPRB)&
    &      .AND. ( ZLON > 335._JPRB .AND. ZLON < 353._JPRB) ) THEN
    IFF=7
  ENDIF
!----- Britain  
  IF ( (ZLAT < 63._JPRB .AND. ZLAT > ZBNDF)&
    &      .AND. ( ZLON > 354._JPRB .OR. ZLON < 3._JPRB) ) THEN
    IFF=9
  ENDIF
!----- Ireland
  IF ( (ZLAT < ZBNDG .AND. ZLAT > ZBNDH)&
    &      .AND. ( ZLON > 349._JPRB .AND. ZLON < 355._JPRB) ) THEN
    IFF=8
  ENDIF
  ITYPDU=1
  IF ( (IFF >= 1 .AND. IFF <= 10) .OR. IFF == 16) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF        

! ITYPDU=2 
!-- Russia to Urals
  IF ( ZLON > 30._JPRB .AND. ZLON <= 70._JPRB ) THEN
    IF ( ZLAT > 51._JPRB ) THEN
      IFF=11
    ELSEIF ( ZLAT > 36._JPRB ) THEN
      IFF=12
    ENDIF
  ENDIF
  ITYPDU=2
  IF ( IFF >= 11 .AND. IFF <= 12 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF        

! ITYPDU=3 
!-- Northern Sahara
!  if ( ( zlat <= 36._JPRB .and. zlat >= 21._JPRB) &
!    & .and. ( zlon >= 330._JPRB .or. zlon <= zbnda) ) then
!    iff=13
!-- Northern Sahara (West)
  IF ( ( ZLAT <= 36._JPRB .AND. ZLAT >= 21._JPRB)&
    & .AND. ( ZLON >= 330._JPRB .OR. ZLON < 2._JPRB) ) THEN
    IFF=36
!-- Northern Sahara (East)
  ELSEIF ( ( ZLAT <= 36._JPRB .AND. ZLAT >= 21._JPRB)&
    & .AND. ( ZLON >= 2._JPRB .OR. ZLON <= ZBNDA) ) THEN
    IFF=37
!-- Southern Sahara (West)
  ELSEIF ( ( ZLAT < 21._JPRB .AND. ZLAT >= 12._JPRB)&
    & .AND. ( ZLON >= 330._JPRB .OR. ZLON < 8._JPRB) ) THEN
    IFF=34
!-- Southern Sahara (East)
  ELSEIF ( ( ZLAT < 21._JPRB .AND. ZLAT >= 12._JPRB)&
    & .AND. ( ZLON >= 8._JPRB .AND. ZLON <= ZBNDA) ) THEN
    IFF=35
!-- Central Africa
  ELSEIF ( ( ZLAT < 12._JPRB .AND. ZLAT >= -12._JPRB)&
    &      .AND. ( ZLON >= 330._JPRB .OR. ZLON <= 60._JPRB) ) THEN
    IFF=14
!-- Southern Africa
  ELSEIF ( ZLAT < -12._JPRB .AND. ZLAT >= -60._JPRB&
    &      .AND. ( ZLON >= 330._JPRB .OR. ZLON <= 60._JPRB) ) THEN
    IFF=15
  ENDIF
  ITYPDU=3
  IF ( (IFF >= 13 .AND. IFF <= 15) .OR. (IFF >= 34 .AND. IFF <= 37) ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF        

!  ITYPDU=4
!-- Australasia
  IF (ZLON > 70._JPRB .AND. ZLON <= 190._JPRB) THEN
    IFF=26

!-- Siberia
    IF (ZLAT <= 90._JPRB .AND. ZLAT > 51._JPRB) THEN
      IFF=16

!-- South Australasia
!---- Tropical Pacific Islands
    ELSEIF ( ZLAT > -10.5_JPRB) THEN
      IFF=27
!---- Australia
    ELSEIF ( ZLAT <= -10.5_JPRB .AND. ZLAT >= -60._JPRB) THEN
      IFF=28
    ENDIF
  ENDIF
  ITYPDU=4
  IF ( IFF >= 26 .AND. IFF <= 28 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF    

! ITYPDU=5
!-- Asian deserts
  IF (ZLON > 90._JPRB .AND. ZLON <= 135._JPRB) THEN
    IF (ZLAT <= 51._JPRB .AND. ZLAT > ZBNDL) THEN
      IFF=17
    ENDIF
  ENDIF

!-- Saudi Arabia
  IF ((ZLAT <= 36._JPRB .AND. ZLAT >= 12._JPRB)&
    & .AND.( ZLON > ZBNDA .AND. ZLON < ZBNDB) ) THEN
    IFF=18
  ENDIF
!-- Irak, Iran, Pakistan
  IF ((ZLAT <= 36._JPRB   .AND. ZLAT >= 20._JPRB)&
    & .AND.( ZLON > ZBNDB .AND. ZLON < 70._JPRB) ) THEN
    IFF=19
  ENDIF

!-- Central Asia and India
  IF ( ZLON > 70._JPRB .AND. ZLON <= 90._JPRB) THEN
!----- Central Asia: Taklamakan
    IF (ZLAT <= 43._JPRB .AND. ZLAT >= ZBNDK) THEN
      IFF=20
!----- India
    ELSEIF (ZLAT <= ZBNDK .AND. ZLAT > 7._JPRB) THEN
      IFF=21
    ENDIF
  ENDIF
!-- other Gobi(s) in South Mongolia and Central China
  IF ( ZLAT <= 49._JPRB .AND. ZLAT > 35._JPRB) THEN
    IF (ZLON > 90._JPRB .AND. ZLON <= 110._JPRB)  THEN
      IFF=22
    ELSEIF (ZLON > 110._JPRB .AND. ZLON <= 125._JPRB) THEN
      IFF=23
    ENDIF
  ENDIF

!-- South China  
  IF (ZLON > 90._JPRB .AND. ZLON <= 135._JPRB) THEN
    IF (ZLAT <= ZBNDL .AND. ZLAT > 7._JPRB) THEN
      IFF=24
    ENDIF
  ENDIF
  ITYPDU=5
  IF ( IFF >= 17 .AND. IFF <= 24 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF    


! ITYPDU=7
!-- Japan and S.Korea
  IF ( (ZLON > 124._JPRB .AND. ZLON < 153._JPRB)&
    & .AND. (ZLAT > 24._JPRB .AND. ZLAT < ZBNDC) ) THEN
    IFF=25
  ENDIF

!-- Greenland
  IF (ZLAT > 50._JPRB) THEN
    ZINCLAT=(90._JPRB-ZLAT)/40._JPRB*45._JPRB
    ZLONW=270._JPRB +ZINCLAT
    ZLONE=360._JPRB -ZINCLAT
    IF ( ZLON > ZLONW .AND. ZLON <  ZLONE ) THEN
      IFF=29
    ENDIF
  ENDIF

!-- Antarctica
  IF (ZLAT < -60._JPRB) THEN
    IFF=30
  ENDIF
  ITYPDU=7
  IF ( (IFF >= 29 .AND. IFF <= 30) .OR. IFF == 25 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF

!-- awaiting a proper recoding, new areas are set between iff=31 and 35
! ITYPDU=6 
  IF ( ZLON > 285._JPRB .AND. ZLON < 295._JPRB) THEN
!- Atacama desert and Salar de Uyuni
    IF ( ZLAT < -16._JPRB  .AND. ZLAT > -28._JPRB) THEN
      IFF=31
    ENDIF
!- Salar de Pipanaco and other small ones
    IF ( ZLAT <= -28._JPRB .AND. ZLAT > ZBNDM) THEN
      IFF=32
    ENDIF
  ENDIF
!- Argentinian pampas
  IF (ZLON > 285._JPRB .AND. ZLON < 297._JPRB) THEN
    IF ( ZLAT <= ZBNDM ) THEN
      IFF=33
    ENDIF
  ENDIF
  ITYPDU=6
  IF ( IFF >= 31 .AND. IFF <= 33 ) THEN
    NDUSRCP(IFF)=ITYPDU
  ENDIF
       
  IF (IFF /= 0) THEN
    PAERMAP(JL,1)=IFF*PLSM(JL)                                   ! area index
    PAERMAP(JL,3)=RDUSRCP(NDUSRCP(IFF),1)                        ! reference speed
    PAERMAP(JL,4)=RDUSRCP(NDUSRCP(IFF),2)                        ! reference particule radius
    DO JAER=1,KBINDD
       !ZDUEMPOT(JL,JAER)=RDDUAER(IFF)*RDDUSRC(IFF,JAER)*PLSM(JL)  ! dust emission potential factor (including land-sea mask)
       ZDUEMPOT(JL,JAER)=RDDUAER(IFF)*RDDUSRC(JAER)*PLSM(JL)  ! dust emission potential factor (including land-sea mask)
    ENDDO
    PAERMAP(JL,2)=ZDUEMPOT(JL,1)                                 ! for diagnostics only
  ELSE
    WRITE(NULOUT,FMT='(''aer_src: Unassigned grid for Lat,Lon='',2F8.2)') ZLAT,ZLON
    PAERMAP(JL,:)=0._JPRB
  ENDIF
ENDDO


!-----------------------------------------------------------------------

!*       0.3   SURFACE WIND VARIABLE RELEVANT FOR SS AND DU EMISSIONS
!              ------------------------------------------------------

!!$IF (NAERWND == 0) THEN
!!$!-- no gust accounted for
!!$  ZWNDDU(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
!!$  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
!!$ELSEIF (NAERWND == 1) THEN
!!$!-- gust only for SS, 10-m wind for DU
!!$  ZWNDDU(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
!!$  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
!!$ELSEIF (NAERWND == 2) THEN
!!$!-- gust only for DU, 10-m wind for SS
!!$  ZWNDDU(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
!!$  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
!!$ELSEIF (NAERWND == 3) THEN
!!$!-- gust for both SS and DU
!!$  ZWNDDU(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
!!$  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
!!$ENDIF

! correction to account for the decrease of mean wind and gusts with decreasing
! time step
!!$IF (PTSPHY < 1000) THEN
!!$  ZWNDDU(KIDIA:KFDIA)=1.06_JPRB*ZWNDDU(KIDIA:KFDIA)
!!$  ZWNDSS(KIDIA:KFDIA)=1.08_JPRB*ZWNDSS(KIDIA:KFDIA)
!!$ENDIF
ZRWPWP=PRWPWP
ZRWSAT=PRWSAT
!*       2.0   DESERT DUST
!              -----------

!- Simplistic lifting from surface based on 10-m wind and surface albedo
ZHDD=MAX(1.0_JPRB,8434._JPRB/1000._JPRB)

!PAERLIF(KIDIA:KFDIA,1:9)=0._JPRB
PAERFLX(KIDIA:KFDIA,1:12,1:9)=0._JPRB
ZFLX_SDUST(KIDIA:KFDIA,1:9,1:12)=0._JPRB

 !-----------------------------------------------


INBAER=0
RAERDUB=1.E-11_JPRB

 !- ECMWF dust emission fluxes come in either 3- or 10-size bins
 ! 0.03 - 0.55 - 0.9 - 20.
 ! 0.03 - 0.06 - 0.12 - 0.24 - 0.48 - 0.96 - 1.92 - 3.84 - 7.68 - 15.36 - 30.72

 !!-- for potential dust sources, select land points, snow-free, and zero ice, no wet skin cover
 !!   with fraction of bare soil > 10%, no high vegetation, possible low vegetation < 50% but 
 !!   with soil moisture below moisture corresponding to twice the wilting point (0.171), and 
 !!   a flatish surface (st.dev.orog < 50) with total albedo < 50%

 !-- for potential dust sources, select land points, snow-free, and zero ice, 
 !   no wet skin cover, with some fraction of bare soil, with test on soil 
 !   moisture, and a flatish surface (st.dev.orog < 50) with total albedo < 52%

  DO JL=KIDIA,KFDIA    
!-- default values for non-land points
    LLPDUSTS(JL)=.FALSE.
    PAERMAP(JL,5)=0.0_JPRB
    ZSCC2(JL)=0._JPRB
    ZDEP2(JL)=0._JPRB
    ZLTS2(JL)=0._JPRB  
    IF (PLSM(JL) >= 0.99_JPRB) THEN
      ZREFSPD = PAERMAP(JL,3)
      ZREFRAD = PAERMAP(JL,4)
      ZRADREF = ZREFSPD * ZREFRAD**0.25_JPRB
!-- default min and max of LTS correspond to PWS1 = ZRWPWP and PSW1 = ZRWSAT
      ZLTSMIN(JL) = 0.6_JPRB * ZRADREF       ! ZFSWET = 0.6 
      ZLTSMAX(JL) = 1.2_JPRB * ZRADREF       ! ZFSWET = 1.2
      ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB, (PWS1(JL)-ZRWPWP)/(ZRWSAT-ZRWPWP) ) )
      ZFSWET = 1.2_JPRB+0.2_JPRB*LOG10(ZSWETN)
!-- background lifting threshold speed (defined for all land points)
      ZLTS2(JL) = MIN( ZLTSMAX(JL), MAX( ZLTSMIN(JL), ZFSWET * ZRADREF ))
!-- replace  by simpler test on:
!     absence of snow
!     flatish surface
!     total albedo < 0.52 (no permanent ice)
!     type 8 fraction bare soil > 10%
!     type 4 cover by dry snow-free low vegetated < 50%
!     all other cover types < 0.1%
      IF (PSNS(JL) < ZEPSSNO .AND. PHSDFOR(JL) <= 50._JPRB .AND. PALB(JL) < 0.52_JPRB .AND.&
        & PFRTI(JL,2) < ZEPSARE .AND. PFRTI(JL,3) < ZEPSARE .AND.&! no ice, no wet skin
        & PFRTI(JL,5) < ZEPSARE .AND. PFRTI(JL,6) < ZEPSARE .AND.&! no snow under bare soil-low veg, no dry high veg
        & PFRTI(JL,7) < ZEPSARE .AND.&! no snow under high veg
        & PFRTI(JL,8) > 0.1_JPRB .AND. PFRTI(JL,4) < 0.5_JPRB ) THEN

         LLPDUSTS(JL)=.TRUE.

         PAERMAP(JL,5)=RAERDUB * ZDUEMPOT(JL,1)                         ! for diagnostics only
      ENDIF
    ENDIF
  ENDDO   

! ZFLX_SDUST is positive in kg m-2 s-1
! but ECMWF conventions have PCFLX as a negative upward flux
! input parameters from climatology are:
!  -- soil clay content        (%)
!     dust emission potential  (kg s2 m-5)
!     lifting thereshold speed (m s-1)
!     
  DO JAER=1,KBINDD

!-- surface source of dust is assumed if LLPTDUSTS is true, and/or UVis albedo > 0.11
!                                        10m wind > threshold = f(soil wetness, mean particle radius)
     DO JL=KIDIA,KFDIA
        LLDUST(JL,:)=.FALSE.
        ZFLX_SDUST(JL,JAER,1:12)=0._JPRB
        
!---------------------------------------------------------------------
!-- ECMWF formulation
        
        IF (LLPDUSTS(JL)) THEN
           ZDEP2(JL)= RAERDUB * ZDUEMPOT(JL,JAER)
           ZSCC2(JL)= 20._JPRB

!-- Present formulation in MACC (June'11, still kept June'13)
!      use a formula of threshold wind velocity modified from Ginoux et al., 2001
!      based on 1st layer soil wetness and an averaged particle radius
!--    All of the above limes computed above
        !PAERLIF(JL,JAER)=ZLTS2(JL)    ! for diagnostics only

!- compare surface 10-m wind with threshold wind velocity

           ZWND3(JL) = MAX(0._JPRB, (PWIND(JL)-ZLTS2(JL)) *PWIND(JL)*PWIND(JL) )

!- preferred approach: flux is based on MODIS-derived UVis_Alb (0.3-0.7 um)
!        IF (LE4ALB) THEN
        IF (NALBEDOSCHEME>0) THEN
!          IF (PALBD(JL,1) >= 0.20_JPRB .AND. PALBD(JL,1) < 0.55_JPRB ) THEN
          IF (PALBD(JL,1) >= 0.08_JPRB .AND. PALBD(JL,1) < 0.55_JPRB ) THEN
            ZFLX_SDUST(JL,JAER,3)= ZDEP2(JL) * PALBD(JL,1) * ZWND3(JL)
          ENDIF
!-- alternate approach, if MODIS-derived albedo not available, use total albedo
        ELSE 
           ZFLX_SDUST(JL,JAER,3)= ZDEP2(JL) * PALB(JL) * ZWND3(JL)
           
        ENDIF

           LLDUST(JL,3)=.TRUE.
        PAERFLX(JL,3,JAER) = ZFLX_SDUST(JL,JAER,3)
        !write(9504,*)ZFLX_SDUST(JL,JAER,3)
        !if (NSTEP==2.or.NSTEP==3 )write(9501,*)JAER,ZDEP2(JL), ZDUEMPOT(JL,JAER),PALB(JL)
        !write(9502,*)PALB(JL) 
        !write(9503,*)ZWND3(JL)
     ENDIF
  ENDDO
ENDDO
ENDIF  ! case NDDUST == 3

!-- PCFLX in kg m-2 s-1

DO JAER=1,KBINDD
   INBAER=INBAER+1
   DO JL=KIDIA,KFDIA
      IF (LLDUST(JL,NDDUST) .AND. ZFLX_SDUST(JL,JAER,NDDUST) > 0._JPRB) THEN
         ZFLX_SDUST(JL,JAER,NDDUST)=ZFLX_SDUST(JL,JAER,NDDUST)
      ENDIF
      !PCFLX(JL,KAERO(INBAER))=-ZFLX_SDUST(JL,JAER,NDDUST) * 1.E+00_JPRB
      if (JAER<2) then
         !----
         ! accumulation mode
         ! number
         emis_number(mode_aci)%d3(JL,KLEV,1)   = emis_number(mode_aci)%d3(JL,KLEV,1) +ZFLX_SDUST(JL,JAER,NDDUST)* 3./(4.*RPI*ddust*mmr_ai**3) * EXP(4.5*LOG(sigma(iacci))**2)*1.E+3
         ! mass
         emis_mass(mode_aci)%d3(JL,KLEV,1)   = emis_mass(mode_aci)%d3(JL,KLEV,1)+ZFLX_SDUST(JL,JAER,NDDUST)!flux_ai(KIDIA:KFDIA)
      else if(JAER>=2 )then
         
         ! ------------------------------
         ! coarse mode
         ! number
         emis_number(mode_coi)%d3(JL,KLEV,1)   = emis_number(mode_coi)%d3(JL,KLEV,1) +ZFLX_SDUST(JL,JAER,NDDUST)* 3./(4.*RPI*ddust*mmr_ci**3) * EXP(4.5*LOG(sigma(icoai))**2)*1.E+3
         ! mass
         emis_mass(mode_coi)%d3(JL,KLEV,1)   = emis_mass(mode_coi)%d3(JL,KLEV,1) +ZFLX_SDUST(JL,JAER,NDDUST)
      end if
   ENDDO
END DO
!-- if no vertical diffusion, distribute the flux in layers with scale height
!-- between half-levels IHTST-1 and KLEV
!!$    IF (.NOT.LVDFTRAC) THEN
!!$      DO JK=IHTST,KLEV
!!$        DO JL=KIDIA,KFDIA
!!$          ZDETAH(JL,JK) = ZETAH(JL,JK)**ZHDD - ZETAH(JL,JK-1)**ZHDD
!!$          PTENC(JL,JK,KAERO(INBAER)) = PTENC(JL,JK,KAERO(INBAER))+PCFLX(JL,KAERO(INBAER))*ZDETAH(JL,JK)
!!$        ENDDO
!!$      ENDDO
!!$    ENDIF
!!$  ENDDO


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ECMWF_DUST_SCHEME',1,ZHOOK_HANDLE)
END SUBROUTINE ECMWF_DUST_SCHEME
