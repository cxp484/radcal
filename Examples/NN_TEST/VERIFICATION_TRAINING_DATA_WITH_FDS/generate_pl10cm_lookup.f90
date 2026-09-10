!> Compare the FDS absorption-coefficient lookup approach against direct RadCal.
!>
!> This program builds five single-species RadCal tables, then evaluates every
!> state in an nn_database.f90 output file with the FDS-style table query.
!> Command-line arguments are the database header, valid-gas-index stream, and
!> kappa stream, in that order.  A compact per-temperature error summary and a
!> bounded scatter-data CSV are written in the directory from which the
!> executable is launched.
!>
!> The implementation intentionally reproduces the relevant radi.f90 behavior:
!>   * 50 logarithmically spaced composition entries in slots 1:NX;
!>   * an unpopulated, zero-valued slot 0 for below-range composition values;
!>   * 45 temperature slices, selected by floor-and-clamp (not interpolation);
!>   * linear interpolation in composition only; and
!>   * AMEAN during table construction, matching nn_database.f90's direct K10.
!>
!> nn_database.f90 stores cm^-1 values and uses gas mole fractions plus soot
!> volume fraction as the RadCal inputs.  Those conventions are retained here
!> so the direct database is a like-for-like comparison.  radi.f90 converts its
!> completed table to m^-1; this utility deliberately leaves values in cm^-1.
PROGRAM GENERATE_PL10CM_LOOKUP
USE ISO_FORTRAN_ENV, ONLY: INT32,INT64
USE PRECISION_PARAMETERS
USE RADCAL_CALC
USE RADCAL_VAR
IMPLICIT NONE
INTEGER,PARAMETER :: NX=50,NTAB=44,NS=5
INTEGER,PARAMETER :: SID(NS)=(/I_CO2,I_H2O,I_CO,I_C2H4,I_FV/)
REAL(EB),PARAMETER :: XMIN=1.E-5_EB,XMAX=1._EB,CMIN=2.E-7_EB,CMAX=.2_EB,TMIN=270._EB,TMAX=2470._EB,P=1._EB
REAL(EB) :: TAB(NS,0:NX,0:NTAB),XFAC,CFAC,LXFAC,LCFAC,TFAC
CHARACTER(512) :: HF,GF,KF,CSVF,SCATTERF,INTERPOLATED_SCATTERF
INTEGER :: NT,NG,NC,UIN,UCSV,USCATTER,UINTERPOLATED,IT,IG,IC,IOS
INTEGER(INT64) :: NV,NSAMP,NRECORD,RECORD,STRIDE
INTEGER(INT32),ALLOCATABLE :: GIDX(:,:)
REAL(EB),ALLOCATABLE :: TG(:),GA(:),CA(:)
REAL(FB) :: KP,KD
REAL(EB) :: X(5),KL,KL_INTERPOLATED,ERR,ERR_INTERPOLATED,ESUM,EMAX,T,TESUM,TEMAX
INTEGER(INT64) :: TNSAMP

IF (COMMAND_ARGUMENT_COUNT()/=3) ERROR STOP 'Usage: generate_pl10cm_lookup header.txt valid_gas_indices.bin kappa.bin'
CALL GET_COMMAND_ARGUMENT(1,HF); CALL GET_COMMAND_ARGUMENT(2,GF); CALL GET_COMMAND_ARGUMENT(3,KF)
CALL READ_HEADER(TRIM(HF),NT,NG,NC,NV,TG,GA,CA)
ALLOCATE(GIDX(4,NV))
OPEN(NEWUNIT=UIN,FILE=TRIM(GF),STATUS='OLD',ACTION='READ',ACCESS='STREAM',FORM='UNFORMATTED'); READ(UIN) GIDX; CLOSE(UIN)
CALL INIT_RC; CALL BUILD_FDS_TABLE
! Write results to the caller's working directory, independent of the input
! database location. This makes repeated comparisons easy to collect or rename.
CSVF='fds_lookup_comparison.csv'
SCATTERF='fds_lookup_scatter.csv'
INTERPOLATED_SCATTERF='fds_lookup_scatter_temperature_interpolated.csv'
OPEN(NEWUNIT=UIN,FILE=TRIM(KF),STATUS='OLD',ACTION='READ',ACCESS='STREAM',FORM='UNFORMATTED')
OPEN(NEWUNIT=UCSV,FILE=TRIM(CSVF),STATUS='REPLACE',ACTION='WRITE')
OPEN(NEWUNIT=USCATTER,FILE=TRIM(SCATTERF),STATUS='REPLACE',ACTION='WRITE')
OPEN(NEWUNIT=UINTERPOLATED,FILE=TRIM(INTERPOLATED_SCATTERF),STATUS='REPLACE',ACTION='WRITE')
WRITE(UCSV,'(A)') 'temperature_k,samples,mean_absolute_percent_error,maximum_absolute_percent_error'
! Keep the sampled state variables so the scatter plot can be traced back to
! the temperature, soot fraction, and gas mixture responsible for an outlier.
WRITE(USCATTER,'(A)') 'temperature_k,x_co2,x_h2o,x_co,x_c2h4,soot_volume_fraction,'// &
                       'direct_cm-1,fds_lookup_cm-1,abs_percent_error'
WRITE(UINTERPOLATED,'(A)') 'temperature_k,x_co2,x_h2o,x_co,x_c2h4,soot_volume_fraction,'// &
                            'direct_cm-1,fds_lookup_cm-1,abs_percent_error'
ESUM=0._EB; EMAX=0._EB; NSAMP=0_INT64
NRECORD=INT(NT,INT64)*NV*INT(NC,INT64)
STRIDE=MAX(1_INT64,NRECORD/100000_INT64)
RECORD=0_INT64
DO IT=1,NT
 T=TG(IT)
 TESUM=0._EB; TEMAX=0._EB; TNSAMP=0_INT64
 DO IG=1,NV
  X=0._EB; X(1)=GA(GIDX(1,IG)+1); X(2)=GA(GIDX(2,IG)+1); X(3)=GA(GIDX(3,IG)+1); X(5)=GA(GIDX(4,IG)+1)
  DO IC=1,NC
   X(4)=CA(IC); READ(UIN,IOSTAT=IOS) KP,KD
   IF (IOS/=0) ERROR STOP 'Unexpected end of kappa database'
   RECORD=RECORD+1_INT64
   KL=GET_KAPPA_FDS(X,T)
   KL_INTERPOLATED=GET_KAPPA_TEMPERATURE_INTERPOLATED(X,T)
   IF (KD>0._FB) THEN
    ERR=100._EB*ABS(KL-REAL(KD,EB))/REAL(KD,EB)
    ESUM=ESUM+ERR; EMAX=MAX(EMAX,ERR); NSAMP=NSAMP+1_INT64
    TESUM=TESUM+ERR; TEMAX=MAX(TEMAX,ERR); TNSAMP=TNSAMP+1_INT64
   ELSE; ERR=0._EB
   ENDIF
   IF (KD>0._FB) THEN
      ERR_INTERPOLATED=100._EB*ABS(KL_INTERPOLATED-REAL(KD,EB))/REAL(KD,EB)
   ELSE
      ERR_INTERPOLATED=0._EB
   ENDIF
   IF (MOD(RECORD-1_INT64,STRIDE)==0_INT64) &
      WRITE(USCATTER,'(9(ES24.16E3,:,","))') T,X(1),X(2),X(3),X(5),X(4),REAL(KD,EB),KL,ERR
   IF (MOD(RECORD-1_INT64,STRIDE)==0_INT64) &
      WRITE(UINTERPOLATED,'(9(ES24.16E3,:,","))') T,X(1),X(2),X(3),X(5),X(4),REAL(KD,EB),KL_INTERPOLATED,ERR_INTERPOLATED
  ENDDO
 ENDDO
 WRITE(UCSV,'(ES24.16E3,",",I0,",",ES24.16E3,",",ES24.16E3)') T,TNSAMP,TESUM/REAL(TNSAMP,EB),TEMAX
ENDDO
CLOSE(UIN); CLOSE(UCSV); CLOSE(USCATTER); CLOSE(UINTERPOLATED)
WRITE(*,'(A)') 'FDS-style lookup comparison complete.'
WRITE(*,'(A,I0)') 'Samples: ',NSAMP
WRITE(*,'(A,ES12.4E3,A)') 'Mean absolute percent error: ',ESUM/REAL(NSAMP,EB),' %'
WRITE(*,'(A,ES12.4E3,A)') 'Maximum absolute percent error: ',EMAX,' %'
WRITE(*,'(A)') 'Per-temperature summary CSV: '//TRIM(CSVF)
WRITE(*,'(A)') 'Bounded scatter-data CSV: '//TRIM(SCATTERF)
WRITE(*,'(A)') 'Temperature-interpolated scatter CSV: '//TRIM(INTERPOLATED_SCATTERF)
CALL CLOSE_RC
CONTAINS
!> Allocate and initialize the one-point RadCal calculation used for table construction.
SUBROUTINE INIT_RC
 NPT=1; OMMIN=50._EB; OMMAX=10000._EB; TWALL=1173.15_EB
 ALLOCATE(PARTIAL_PRESSURES_ATM(N_RADCAL_SPECIES,1),TEMP_GAS(1),SEGMENT_LENGTH_M(1),TOTAL_PRESSURE_ATM(1))
 SEGMENT_LENGTH_M=.1_EB; CALL RCALLOC; CALL INIT_RADCAL
END SUBROUTINE INIT_RC
SUBROUTINE CLOSE_RC
 CALL RCDEALLOC2; CALL RCDEALLOC
 DEALLOCATE(PARTIAL_PRESSURES_ATM,TEMP_GAS,SEGMENT_LENGTH_M,TOTAL_PRESSURE_ATM)
END SUBROUTINE CLOSE_RC
!> Populate the lookup tables with pure absorber plus nitrogen background states.
!> Each table value is AMEAN in cm^-1, matching the direct K10 database.
!> Unlike radi.f90, this comparison utility does not apply MIN(AMEAN,AP0),
!> because that would compare different RadCal quantities.
SUBROUTINE BUILD_FDS_TABLE
 INTEGER :: I,J,S
 REAL(EB) :: XX,YY,TT,AMEAN,AP0,RAD,TRAN
 ! Exact radi.f90 spacing: 50 populated locations, in slots 1:50; slot 0 is zero.
 XFAC=(XMAX/XMIN)**(1._EB/REAL(NX-1,EB)); LXFAC=LOG(XFAC)
 CFAC=(CMAX/CMIN)**(1._EB/REAL(NX-1,EB)); LCFAC=LOG(CFAC)
 TFAC=REAL(NTAB,EB)/(TMAX-TMIN); TAB=0._EB
 DO J=0,NTAB
  TT=TMIN+REAL(J,EB)*(TMAX-TMIN)/REAL(NTAB,EB)
  DO I=0,NX-1
   XX=MIN(1._EB,XMIN*XFAC**I); YY=MIN(1._EB,CMIN*CFAC**I)
   DO S=1,NS
    PARTIAL_PRESSURES_ATM=0._EB
    IF (SID(S)==I_FV) THEN
     PARTIAL_PRESSURES_ATM(I_FV,1)=YY; PARTIAL_PRESSURES_ATM(I_N2,1)=P
    ELSE
     PARTIAL_PRESSURES_ATM(SID(S),1)=XX*P; PARTIAL_PRESSURES_ATM(I_N2,1)=(1._EB-XX)*P
    ENDIF
    TEMP_GAS=TT; TOTAL_PRESSURE_ATM=P; CALL SUB_RADCAL(AMEAN,AP0,RAD,TRAN)
    TAB(S,I+1,J)=AMEAN
   ENDDO
  ENDDO
 ENDDO
END SUBROUTINE BUILD_FDS_TABLE
!> Evaluate the lookup table using the GET_KAPPA interpolation/indexing rules in radi.f90.
!> XIN order is CO2, H2O, CO, soot volume fraction, C2H4, matching nn_database.f90.
REAL(EB) FUNCTION GET_KAPPA_FDS(XIN,TMP)
 REAL(EB),INTENT(IN)::XIN(5),TMP
 INTEGER::TY
 TY=MAX(0,MIN(NTAB,INT((TMP-TMIN)*TFAC)))
 GET_KAPPA_FDS=GET_KAPPA_AT_TINDEX(XIN,TY)
END FUNCTION GET_KAPPA_FDS

!> Evaluate the same composition interpolation, but linearly interpolate the
!> two adjacent temperature slices. This is intentionally not FDS GET_KAPPA;
!> it isolates the accuracy benefit of temperature interpolation.
REAL(EB) FUNCTION GET_KAPPA_TEMPERATURE_INTERPOLATED(XIN,TMP)
 REAL(EB),INTENT(IN)::XIN(5),TMP
 INTEGER::TY_LOW,TY_HIGH
 REAL(EB)::TFRAC
 IF (TMP<=TMIN) THEN
  GET_KAPPA_TEMPERATURE_INTERPOLATED=GET_KAPPA_AT_TINDEX(XIN,0)
 ELSEIF (TMP>=TMAX) THEN
  GET_KAPPA_TEMPERATURE_INTERPOLATED=GET_KAPPA_AT_TINDEX(XIN,NTAB)
 ELSE
  TY_LOW=MAX(0,MIN(NTAB-1,INT((TMP-TMIN)*TFAC)))
  TY_HIGH=TY_LOW+1
  TFRAC=(TMP-TMIN)*TFAC-REAL(TY_LOW,EB)
  GET_KAPPA_TEMPERATURE_INTERPOLATED=GET_KAPPA_AT_TINDEX(XIN,TY_LOW)+TFRAC* &
    (GET_KAPPA_AT_TINDEX(XIN,TY_HIGH)-GET_KAPPA_AT_TINDEX(XIN,TY_LOW))
 ENDIF
END FUNCTION GET_KAPPA_TEMPERATURE_INTERPOLATED

!> Apply FDS's composition interpolation at one specified temperature index.
REAL(EB) FUNCTION GET_KAPPA_AT_TINDEX(XIN,TY)
 REAL(EB),INTENT(IN)::XIN(5)
 INTEGER,INTENT(IN)::TY
 INTEGER::S,L,U
 REAL(EB)::F,XV
 GET_KAPPA_AT_TINDEX=0._EB
 DO S=1,NS
  SELECT CASE(SID(S))
  CASE(I_CO2);  XV=XIN(1)
  CASE(I_H2O);  XV=XIN(2)
  CASE(I_CO);   XV=XIN(3)
  CASE(I_C2H4); XV=XIN(5)
  CASE(I_FV);   XV=XIN(4)
  END SELECT
  IF (XV<=0._EB) CYCLE
  IF (SID(S)==I_FV) THEN; F=LOG(XV/CMIN)/LCFAC+1._EB
  ELSE; F=LOG(XV/XMIN)/LXFAC+1._EB
  ENDIF
  F=MAX(0._EB,MIN(REAL(NX,EB),F)); L=INT(F); F=F-L; L=MIN(L,NX); U=MIN(L+1,NX)
  GET_KAPPA_AT_TINDEX=GET_KAPPA_AT_TINDEX+TAB(S,L,TY)+F*(TAB(S,U,TY)-TAB(S,L,TY))
 ENDDO
END FUNCTION GET_KAPPA_AT_TINDEX
!> Read the text metadata and axes written by nn_database.f90.
SUBROUTINE READ_HEADER(FN,NT,NG,NC,NV,TG,GA,CA)
 CHARACTER(*),INTENT(IN)::FN; INTEGER,INTENT(OUT)::NT,NG,NC; INTEGER(INT64),INTENT(OUT)::NV
 REAL(EB),ALLOCATABLE,INTENT(OUT)::TG(:),GA(:),CA(:)
 CHARACTER(1024)::L; INTEGER::U,I,IDX; REAL(EB)::V
 OPEN(NEWUNIT=U,FILE=FN,STATUS='OLD',ACTION='READ')
 DO
  READ(U,'(A)') L
  IF(INDEX(L,'N_TEMPERATURE=')>0) READ(L(INDEX(L,'=')+1:),*) NT
  IF(INDEX(L,'N_GAS_AXIS=')>0) READ(L(INDEX(L,'=')+1:),*) NG
  IF(INDEX(L,'N_SOOT_AXIS=')>0) READ(L(INDEX(L,'=')+1:),*) NC
  IF(INDEX(L,'N_VALID_GAS=')>0) READ(L(INDEX(L,'=')+1:),*) NV
  IF(INDEX(L,'[TEMPERATURE_K]')>0) EXIT
 ENDDO
 ALLOCATE(TG(NT),GA(NG),CA(NC)); DO I=1,NT; READ(U,*)IDX,V; TG(I)=V; ENDDO
 DO; READ(U,'(A)')L; IF(INDEX(L,'[GAS_MOLE_FRACTION]')>0)EXIT; ENDDO
 DO I=1,NG; READ(U,*)IDX,V; GA(I)=V; ENDDO
 DO; READ(U,'(A)')L; IF(INDEX(L,'[SOOT_VOLUME_FRACTION]')>0)EXIT; ENDDO
 DO I=1,NC; READ(U,*)IDX,V; CA(I)=V; ENDDO; CLOSE(U)
END SUBROUTINE READ_HEADER
END PROGRAM GENERATE_PL10CM_LOOKUP
