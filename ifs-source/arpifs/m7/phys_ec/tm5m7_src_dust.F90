SUBROUTINE TM5M7_SRC_DUST( YDEPHY, YDEAERMAP, YDEAERSRC,                 &
                         & KIDIA, KFDIA, KLON, KLEV, KTILES, KSW,        &
                         & PLSM , PWIND, PSNS, PZ0M,                     &
                         & SP, PTL, PSOIL_TYPE,                          &
                         & PFRTI, PCVL, PCVH, KTVL, KTVH,                &
                         & EMIS_MASS, EMIS_NUMBER ,PAERFLX,PGLON, PGLAT, &
                         & PRWPWP,PRWSAT,PAERMAP,PALB,PALBD,PWS1,PHSDFOR,&
                         & IMM,ISOILPH1, ISOILPH2, ISOILPH3, ISOILPH4, ISOILPH5, &
                         & IZ0AM, IPOTSRC, ISOILTYPE, IAREA, ICULT,IZ0M, IFPAR, GPGAW,&
                         & ILAI_MAX,ILAI_AVG)

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                      (updated 15-Sep-2026) │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │  *tm5m7_src_dust* - dispatch dust source terms to the configured scheme    │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *tm5m7_src_dust* is called from tm5m7_src                                │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │     Sep. 2026 - split Tegen and ECMWF dust schemes into separate modules   │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯

USE PARKIND1,        ONLY : JPIM, JPRB
USE YOMHOOK,         ONLY : LHOOK, DR_HOOK, JPHOOK
USE TM5M7_DATA,      ONLY : NMOD
USE TM5M7_EMIS_DATA, ONLY : MODAL_EMISSIONS
USE YOEPHY,          ONLY : TEPHY
USE YOEAERMAP,       ONLY : TEAERMAP
USE YOEAERSRC,       ONLY : TEAERSRC
USE TEGEN_DUST_SCHEME_MOD, ONLY : TEGEN_DUST_SCHEME
USE ECMWF_DUST_SCHEME_MOD, ONLY : ECMWF_DUST_SCHEME

INTEGER(KIND=JPIM),     INTENT(IN)    :: IMM
TYPE(TEPHY),            INTENT(IN)    :: YDEPHY
TYPE(TEAERMAP),         INTENT(INOUT) :: YDEAERMAP
TYPE(TEAERSRC),         INTENT(IN)    :: YDEAERSRC
INTEGER(KIND=JPIM),     INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),     INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),     INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),     INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),     INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),     INTENT(IN)    :: KSW
REAL(KIND=JPRB),        INTENT(IN)    :: GPGAW(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PWIND(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PSNS(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PZ0M(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: SP(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PTL(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PSOIL_TYPE(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: PFRTI(KLON,KTILES)
REAL(KIND=JPRB),        INTENT(IN)    :: PCVL(KLON), PCVH(KLON)
INTEGER(KIND=JPIM),     INTENT(IN)    :: KTVL(KLON), KTVH(KLON)
TYPE(MODAL_EMISSIONS),  INTENT(INOUT) :: EMIS_MASS(NMOD)
TYPE(MODAL_EMISSIONS),  INTENT(INOUT) :: EMIS_NUMBER(NMOD)
REAL(KIND=JPRB),        INTENT(INOUT) :: PAERFLX(KLON,12,9)
REAL(KIND=JPRB),        INTENT(IN)    :: PGLON(KLON), PGLAT(KLON)
REAL(KIND=JPRB),        INTENT(INOUT) :: PRWPWP, PRWSAT, PAERMAP(KLON,5)
REAL(KIND=JPRB),        INTENT(IN)    :: PALB(KLON), PALBD(KLON,KSW)
REAL(KIND=JPRB),        INTENT(IN)    :: PWS1(KLON), PHSDFOR(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: ISOILPH1(KLON), ISOILPH2(KLON), ISOILPH3(KLON), ISOILPH4(KLON), ISOILPH5(KLON), &
                                        & IZ0AM(KLON), IPOTSRC(KLON), IAREA(KLON), ICULT(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: IZ0M(KLON), IFPAR(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: ILAI_MAX(KLON), ILAI_AVG(KLON)
REAL(KIND=JPRB),        INTENT(IN)    :: ISOILTYPE(KLON)
REAL(KIND=JPHOOK)                     :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TM5M7_SRC_DUST',0,ZHOOK_HANDLE)

SELECT CASE (YDEAERSRC%NDDUST)
CASE (8)
  CALL TEGEN_DUST_SCHEME( YDEAERSRC,                                   &
                        & KIDIA, KFDIA, KLON, KLEV, KTILES,            &
                        & PLSM , PWIND, PSNS,                          &
                        & SP, PTL,                                     &
                        & PFRTI, PCVL, PCVH, KTVL, KTVH,                &
                        & EMIS_MASS, EMIS_NUMBER,                       &
                        & ISOILPH3, ISOILPH4,                           &
                        & IZ0AM, IPOTSRC, ISOILTYPE, ICULT, IZ0M, IFPAR,&
                        & ILAI_MAX,ILAI_AVG)
CASE (3)
  CALL ECMWF_DUST_SCHEME( YDEPHY, YDEAERMAP,                            &
                        & KIDIA, KFDIA, KLON, KLEV, KTILES, KSW,        &
                        & PLSM , PWIND, PSNS,                           &
                        & PFRTI,                                        &
                        & EMIS_MASS, EMIS_NUMBER, PGLON, PGLAT,         &
                        & PRWPWP,PRWSAT,PAERMAP,PALB,PALBD,PWS1,PHSDFOR)
END SELECT

IF (LHOOK) CALL DR_HOOK('TM5M7_SRC_DUST',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_SRC_DUST
