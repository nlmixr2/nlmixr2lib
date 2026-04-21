$PROB NONMEM code of final model
$SIZES LTH=60 LVR=60 LNP4=10000 
$INPUT ID TIME CONC LNCONC=DV AMT EVID CMT MDV BLQ TYPE PERIOD
;TYPE: TYPE=0, parent drug agomelatine ; TYPE=1, metabolite 3-hydroxy-agomelatine; TYPE=2, metabolite 7-desmethyl-agomelatine
$DATA  …..
IGNORE=@
$SUBROUTINES ADVAN13 TOL=6
$MODEL
COMP (DEPOT1)
COMP (DEPOT2)
COMP (LIVER)
              COMP (CENTPRNT,DEFOBS) ;Central CMT of agomelatine
              COMP (3OHMETB)  ;Central CMT of 3-hydroxy-agomelatine
COMP (7DMMETB) ;Central CMT of 7-desmethyl-agomelatine
COMP (7DMPERI)   ;Peripheral CMT of 7-desmethyl-agomelatine
COMP (ALMTPERI) ;Peripheral CMT of agomelatine
$PK
;----------------------------Set up dosing times----------------------------
IF(AMT.GT.0) TDOS = TIME ; If AMT > 0, set TDOS to TIME
TAD = TIME-TDOS ; Set time after dose
;----------------------------OCCASIONS---------------------------------------
OOC1               = 0
 IF(PERIOD.EQ.1) OOC1 = 1
OOC2               = 0
 IF(PERIOD.EQ.2) OOC2 = 1
OOC3               = 0
 IF(PERIOD.EQ.3) OOC3 = 1
OOC4               = 0
 IF(PERIOD.EQ.4) OOC4 = 1
;-------------------------------IOV-----------------------------------------------
IOV1=OOC1*ETA(16)+OOC2*ETA(17)+OOC3*ETA(18)+OOC4*ETA(19)
IOV2=OOC1*ETA(20)+OOC2*ETA(21)+OOC3*ETA(22)+OOC4*ETA(23)
IOV3=OOC1*ETA(24)+OOC2*ETA(25)+OOC3*ETA(26)+OOC4*ETA(27)
IOV4=OOC1*ETA(28)+OOC2*ETA(29)+OOC3*ETA(30)+OOC4*ETA(31)
IOV5=OOC1*ETA(32)+OOC2*ETA(33)+OOC3*ETA(34)+OOC4*ETA(35)
;-----------------------PK parameters for agomelatine-------------------
TVK13 =THETA(1)
MU_1=LOG(TVK13)
K13   =EXP(MU_1+ETA(1))*EXP(IOV1)
TVV4 = THETA(2)
MU_2=LOG(TVV4)
V4   =EXP(MU_2+ETA(2))  ;Central volume of distribution of agomelatine
TVCLint = THETA(3)
MU_3=LOG(TVCLint)
CLint  = EXP(MU_3+ETA(3))*EXP(IOV4) ;Intrinsic clearance
TVALAG1 = THETA(8)
MU_8=LOG(TVALAG1)
ALAG1 = EXP(MU_8+ETA(8))
TVK23 = THETA(9)
MU_9=LOG(TVK23)
K23 = EXP(MU_9+ETA(9))*EXP(IOV3)
TVALAG2 = THETA(10)
MU_10=LOG(TVALAG2)
ALAG2 = EXP(MU_10+ETA(10))*EXP(IOV2) 
MU_11=LOG(THETA(11)/(1-THETA(11))) ; logit of population F1
EXPP=MU_11+ETA(11) ; add random effect 
F1 = EXP(EXPP+IOV5)/(1+EXP(EXPP+IOV5)) ; individual F1 (Fraction of the dose absorbed by the DEPOT1)
F2 = 1-F1 ; fraction of the dose absorbed through the DEPOT2
;--------------------------------------Liver volume and blood flow------------------------------------------
LV  = 0.05012*WT**0.78        ; Liver volume (LV) in liters, Body weight (WT) in kg
V3=LV
PBR = 1/0.69                             ; plasma to blood ratio: CP/CB=1/0.69     
QH  = 50.4*LV                           ; Liver blood flow (male)
;------------------------------------Well-stirred model for LIVER CMT-------------------------------------
FU=0.05  ;unbound fraction of agomelatine in human plasma 
EH  = FU*PBR*CLint/(QH +FU*PBR*CLint)  ; Liver extraction ratio (ER)
FH=1-EH    ;fraction of absorbed drug escaping hepatic first-pass extraction (FH)
CLH = QH*EH/PBR  ;hepatic plasma clearance
CL  = CLH  ;CL  = CLH + CLR, CLR (renal CL) near to zero since only 0.01% of the dose excreted from urine
;------------------------------------PK parameters for metabolites----------------------------------------
MU_4=THETA(4)
BA3=MU_4+ETA(4)
MU_5=THETA(5)
BA4=MU_5+ETA(5)
FM3OH=EXP(BA3)/(1+EXP(BA3)+EXP(BA4)) ;Fractional formation of 3-hydroxy-agomelatine
FM7DM=EXP(BA4)/(1+EXP(BA3)+EXP(BA4)) ;Fractional formation of 7-desmethyl-agomelatine
TVCL3OH=THETA(6)
MU_6=LOG(TVCL3OH)
CL3OH=EXP(MU_6+ETA(6))    ;Clearance of 3-hydroxy-agomelatine
TVCL7DM=THETA(7)
MU_7=LOG(TVCL7DM)
CL7DM=EXP(MU_7+ETA(7))    ;Clearance of 7-desmethyl-agomelatine
TVQ7DM=THETA(12)
MU_12=LOG(TVQ7DM)
Q7DM=EXP(MU_12+ETA(12))    ;Compartmental clearance of 7-desmethyl-agomelatine
TVV7DM=THETA(13)
MU_13=LOG(TVV7DM)
V7DM=EXP(MU_13+ETA(13))    ;Peripheral volume of distribution of 7-desmethyl-agomelatine
TVQALMT=THETA(14)
MU_14=LOG(TVQALMT)
QALMT=EXP(MU_14+ETA(14))    ;Compartmental clearance of agomelatine
TVVALMT=THETA(15)
MU_15=LOG(TVVALMT)
VALMT=EXP(MU_15+ETA(15))    ;Peripheral volume of distribution of agomelatine
;-------------------------metabolite-to-parent ratio of molecular weights-----------------------------
MPR1=259/243 ;3-hydroxy-agomelatine to parent drug
MPR2=229/243 ;7-desmethyl-agomelatine to parent drug
V5=V4  
V6=V4  
S4=V4
S5=V5
S6=V6
$DES
DADT(1) = -K13*A(1)
DADT(2) = -K23*A(2)
DADT(3) = K13*A(1)+K23*A(2)-QH*FH*A(3)/V3+QH*A(4)/V4-CLH*A(3)/V3 
DADT(4) = QH*FH*A(3)/V3-QH*A(4)/V4-A(4)*QALMT/V4+A(8)*QALMT/VALMT DADT(5)=FM3OH*CLH*A(3)/V3*MPR1-A(5)*CL3OH/V5  
DADT(6)=FM7DM*CLH*A(3)/V3*MPR2-A(6)*CL7DM/V6-A(6)*Q7DM/V6+A(7)*Q7DM/V7DM DADT(7)=A(6)*Q7DM/V6-A(7)*Q7DM/V7DM  
DADT(8)=A(4)*QALMT/V4-A(8)*QALMT/VALMT
$ERROR (OBSERVATION ONLY)
CALMT=A(4)/V4
C3OH=A(5)/V5
C7DM=A(6)/V6
IPRED = -6
IF(F.GT.0) IPRED = LOG(F)
SDALMT=THETA(16)
SD3OH=THETA(17)
SD7DM=THETA(18)
LLOQALMT=LOG(0.0457)  ;natural log of LLOQ of ALMT
LLOQ7DM=LOG(0.1372)    ;natural log of LLOQ of 7DM
CUMDALMT=PHI((LLOQALMT-LOG(CALMT))/SDALMT)
CUMD7DM=PHI((LLOQ7DM-LOG(C7DM))/SD7DM)
;------------------- Agomelatine plasma records: above LLOQ (observations) ------------------
IF(TYPE.EQ.0.AND.BLQ.EQ.0) THEN
F_FLAG=0 ;regular likelihood for measured concentrations
Y = LOG(CALMT) + SDALMT*EPS(1)
ENDIF
;------------------- Agomelatine plasma records: below LLOQ -------------------------------------
IF(TYPE.EQ.0.AND.BLQ.EQ.1) THEN
F_FLAG=1 ;probability that F is less than LLOQ for missing concentrations
Y=CUMDALMT
MDVRES=1
ENDIF
;---------------------------------- 3-hydroxy-agomelatine ----------------------------------------------
IF (TYPE.EQ.1) THEN
Y = LOG(C3OH) + SD3OH*EPS(2)
ENDIF
;-----------7-desmethyl-agomelatine plasma records: above LLOQ (observations) --------
IF(TYPE.EQ.2.AND.BLQ.EQ.0) THEN
F_FLAG=0 ;regular likelihood for measured concentrations
Y = LOG(C7DM)+SD7DM*EPS(3) 
ENDIF
;-------------------7-desmethyl-agomelatine plasma records: below LLOQ-------------------
IF(TYPE.EQ.2.AND.BLQ.EQ.1) THEN
F_FLAG=1 ;probability that F is less than LLOQ for missing concentrations
Y=CUMD7DM
MDVRES=1
ENDIF
$THETA ;final estimates
(0, 4.54) ; K13 (1/h)
(0, 64.6) ; V4 (L)
(0, 111000) ;CLint (L/h)
(-1.78) ; BA3
(-3.95) ; BA4
(0, 44.9) ;CL3OH (L/h)
(0, 52.9) ;CL7DM (L/h)
(0, 0.185) ;ALAG1 (h)
(0, 4.23) ; K23 (1/h)
(0, 0.305) ;ALAG2 (h)
(0.681) ;F1
(0, 28.1) ;Q7DM (L/h)
(0, 536) ;V7DM (L)
(0, 4.36) ;QALMT (L/h)
(0, 157) ;VALMT (L)
(0.39) ;SDALMT
(0.228) ;SD3OH
(0.297) ;SD7DM
$OMEGA
 0.101 ;IIV_K13
 0.0374 ;IIV_V4
 0.998 ;IIV_CLint
$OMEGA BLOCK(2)
 0.175  ;IIV_BA3
 0.153  0.231 ;IIV_BA4
 $OMEGA BLOCK(2)
 0.0414  ;IIV_CL3OH
 0.0492  0.0774 ;IIV_CL7DM
 $OMEGA
 0.0628 ;IIV_ALAG1
 3.78 ;IIV_K23
 3.3 ;IIV_ALAG2
 0.526 ;IIV_F1
 0.001 FIX  ;IIV_Q7DM
 0.001 FIX  ;IIV_V7DM
 1.46 ;IIV_QALMT
 0.001 FIX  ;IIV_VALMT
$OMEGA  BLOCK(1) 1.52 ;IOV1 on K13
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) 4.32 ;IOV2 on ALAG2
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) 5.01 ;IOV3 on K23
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) 0.0779 ;IOV4 on CLint
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) 2.32 ;IOV5 on F1
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$SIGMA
 1 FIX  ; RUV_AMLT
 1 FIX  ; RUV_3OH
 1 FIX  ; RUV_7DM
$EST METHOD=SAEM LAPLACIAN INTERACTION NUMERICAL SLOW NBURN=2000 ISAMPLE=2 CTYPE=3 GRD=TS(16-18) NITER=1000 PRINT=50
$EST METHOD=IMPMAP LAPLACIAN INTERACTION EONLY=1 ISAMPLE=3000 NITER=8 GRD=TS(16-18) PRINT=1
$COV UNCONDITIONAL PRINT=E
