Wu_2025_paracetamol <- function() {
  description <- "Parent-and-three-metabolites population PK model for intravenous and rectal paracetamol (PCM) and its glucuronide (PCM-GLU), sulfate (PCM-SULF), and combined oxidative metabolites (PCM-cysteine + PCM-mercapturate, PCM-OXI, denoted with the canonical cysmer suffix) from preterm and term neonates through infants, children, and adults (Wu 2025). Two-compartment plasma disposition for parent PCM with three parallel formation clearances and parallel renal elimination of unchanged parent; one-compartment plasma disposition for each metabolite with renal elimination expressed as a fraction of glomerular filtration rate (GFR). The preterm-and-term-neonate-to-adult (PTNA) maturation equation (Wu 2024) is applied to each formation clearance and to a separate PCM-SULF renal-secretion clearance; an additional adult-only correction factor scales the renal clearance of PCM-GLU in subjects >= 18 years. Rectal absorption parameters (Ka, Tlag, F) are fixed from Wang 2014."
  reference <- paste(
    "Wu Y, Voller S, Goulooze SC, Allegaert K, Sherwin CMT,",
    "van Rongen A, Roofthooft DWE, Simons SHP, Tibboel D, Flint RB,",
    "van den Anker JN, Knibbe CAJ (2025).",
    "A Novel Maturation Equation for Hepatic Clearance Across Preterm,",
    "Term Neonates, Children, and Adults: Application to Paracetamol",
    "and Its Metabolite. J Clin Pharmacol 65(12):1829-1843.",
    "doi:10.1002/jcph.70080.",
    "GFR-maturation backbone reproduced from Wu Y et al.",
    "Pharm Res 2024;41:637-649 (doi:10.1007/s11095-024-03677-3).",
    "Rectal-absorption parameters (Ka, Tlag, F) reproduced from",
    "Wang C et al. J Clin Pharmacol 2014;54(6):619-629",
    "(doi:10.1002/jcph.259).",
    sep = " "
  )
  vignette <- "Wu_2025_paracetamol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight (time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling on PCM disposition (V_P, V_PP, Q) with reference 1.08 kg (paper's BWc reference 1080 g, Table 3 'Paracetamol' block), and on PTNA-equation CL_max terms (formation CL of PCM-GLU, PCM-SULF, PCM-OXI and PCM-SULF renal secretion) with reference 1.75 kg (Table 3). The same WT enters the Wu 2024 GFR maturation equation with reference 1.75 kg (Equation 3). Source column 'CW' on input.",
      source_name        = "CW"
    ),
    WT_BIRTH = list(
      description        = "Body weight at birth (time-fixed per subject)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling on the PTNA-equation CL_birth terms for the three formation clearances (Table 3) and on the Wu 2024 GFR maturation equation (Equation 3), both with reference 1.75 kg. Source column 'Bwb' on input.",
      source_name        = "Bwb"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling of the PTNA maturation PNA50 term with reference 34 weeks. Per paper Methods, GA was set to 40 weeks for datasets 6, 7, and 8 (infants, children, adults) where the true GA was unavailable.",
      source_name        = "GA"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The canonical PNA carries months; the paper expresses PNA in days throughout. Inside model() PNA is converted to days (PNA_days = PNA * 30.4375) for use in the PTNA sigmoidal-Emax maturation terms (paper's PNA50 estimates in days are kept verbatim). PNA in years (PNA / 12) is also used to flag the 'adult' (>= 18 years) subset to which the f_GLU,adult correction factor on renal PCM-GLU CL is applied (Table 3, paper Results: Model Refinement).",
      source_name        = "PNA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 298L,
    n_studies      = 8L,
    age_range      = "0 days to 50 years (postnatal age); 23-41 weeks gestational age in the neonatal subset (datasets 1-5, n = 235); GA fixed at 40 weeks for the infant/child/adult subset (datasets 6-8, n = 63).",
    age_median     = "Pooled PNA median 5 days for the 235 neonates and 3 years for the 63 infants/children/adults (Tables 1-2).",
    weight_range   = "460 g (preterm neonate, dataset 2) to 91,700 g (adult, dataset 8).",
    weight_median  = "1,190 g for neonates (n = 235); 15,000 g for infants/children/adults (n = 63).",
    sex_female_pct = 55.7,
    race_ethnicity = NA_character_,
    disease_state  = "Mixed cohort: NICU preterm and term neonates with clinical indication for IV PCM or IV propacetamol analgesia (datasets 1-5, n = 235); infants after craniofacial surgery receiving rectal PCM or IV propacetamol (dataset 6, n = 26); children after adenotonsillectomy receiving rectal PCM (dataset 7, n = 29); healthy adults undergoing oral/maxillofacial surgery receiving IV PCM (dataset 8, n = 8).",
    dose_range     = "Single or multiple IV PCM doses 10, 15, or 20 mg/kg (neonates); IV propacetamol 20-40 mg/kg (10-20 mg PCM equivalent /kg); rectal PCM 20-40 mg/kg (infants/children); IV PCM 2000 mg followed by 1000 mg every 6 h (adults). One propacetamol unit = 0.5 units of PCM for dosing input.",
    regions        = "Netherlands (datasets 1, 2, 6, 8), United States (dataset 3), Belgium (datasets 4, 5), and one international children cohort (dataset 7).",
    notes          = "Pooled retrospective + prospective IV/rectal PCM/propacetamol PK data; plasma observations for all datasets, urine observations also collected for datasets 1, 3, and 5. 6428 total observations. Demographics from Tables 1-2 of Wu 2025."
  )

  ini({
    # =====================================================================
    # PCM disposition (Table 3, "Paracetamol" block).
    # Reference current weight 1.08 kg (paper's BWc/1080 g normalisation).
    # =====================================================================
    lvp        <- log(1.01)    ; label("PCM central volume V_P at WT = 1.08 kg (L)")            # Table 3 TVvp = 1.01 L (RSE 2.8%)
    e_wt_vp    <- 0.9781       ; label("Allometric exponent of WT on PCM V_P (unitless)")       # Table 3 theta_BWc(V_P) = 0.9781 (RSE 2%)
    lvpp       <- log(0.2707)  ; label("PCM peripheral volume V_PP at WT = 1.08 kg (L)")        # Table 3 TVvpp = 0.2707 L (RSE 9.8%)
    e_wt_vpp   <- fixed(0)     ; label("Allometric exponent of WT on PCM V_PP (fixed)")         # Table 3 theta_BWc(V_PP) = 0 FIXED
    lq         <- log(0.08921) ; label("PCM inter-compartmental CL Q at WT = 1.08 kg (L/h)")    # Table 3 TVQ = 0.08921 L/h (RSE 23%)
    e_wt_q     <- 2.212        ; label("Allometric exponent of WT on PCM Q (unitless)")         # Table 3 theta_BWc(Q) = 2.212 (RSE 26%) -- estimated >> 0.75; combined with fixed V_PP this collapses peripheral kinetics for larger subjects

    # Rectal absorption (datasets 6 and 7). Ka, Tlag, F all fixed from
    # Wang 2014 (J Clin Pharmacol 54:619-629) per Methods page 1833.
    lka        <- fixed(log(0.275))   ; label("Rectal first-order absorption rate Ka (1/h)")    # Table 3 Ka = 0.275 /h FIXED a (Wang 2014)
    ltlag      <- fixed(log(0.0103))  ; label("Rectal absorption lag time Tlag (h)")            # Table 3 Tlag = 0.0103 h FIXED a (Wang 2014)
    lfdepot    <- fixed(log(0.96))    ; label("Rectal bioavailability F (unitless)")            # Table 3 F = 0.96 FIXED a (Wang 2014)

    # =====================================================================
    # Formation CL of PCM-GLU via the PTNA equation.
    # IIV is applied multiplicatively to the whole derived CL_form,GLU
    # expression (Table 3); the lcl_gluc placeholder anchors etalcl_gluc.
    # Reference birthweight and current weight both 1.75 kg.
    # =====================================================================
    lcl_gluc            <- fixed(log(1))  ; label("Multiplicative anchor for CL_form(GLU) IIV (fixed at 1)") # placeholder so etalcl_gluc has a paired log-parameter (unit-multiplier in model())
    lclbirth_gluc       <- log(0.007961)  ; label("PCM-GLU formation CL at birth at WT_BIRTH = 1.75 kg (L/h)") # Table 3 TVCLbirth(GLU) = 0.007961 L/h (RSE 18%)
    e_bwbirth_cl_gluc   <- 1.515          ; label("WT_BIRTH exponent on CL_birth(GLU) (unitless)")             # Table 3 theta_BWb(GLU) = 1.515 (RSE 16%)
    lclmax_gluc         <- log(0.3378)    ; label("PCM-GLU formation CL_max at WT = 1.75 kg (L/h)")             # Table 3 TVCLmax(GLU) = 0.3378 L/h (RSE 15%)
    e_wt_cl_gluc        <- fixed(0.738)   ; label("WT exponent on CL_max(GLU) (fixed)")                         # Table 3 theta_BWc(GLU) = 0.738 FIXED (paper Methods: PTNA equation)
    lpna50_gluc         <- log(62.67)     ; label("PCM-GLU PNA50 (days)")                                       # Table 3 TVPNA50(GLU) = 62.67 days (RSE 23%)
    e_ga_pna50_gluc     <- -4.625         ; label("GA exponent on PNA50(GLU) (unitless)")                       # Table 3 theta_GAPNA50(GLU) = -4.625 (RSE 14%)
    lhill_gluc          <- log(1.553)     ; label("PCM-GLU PTNA Hill coefficient (unitless)")                   # Table 3 Hill(GLU) = 1.553 (RSE 8.9%)

    # =====================================================================
    # Formation CL of PCM-SULF via the PTNA equation (no-GA variant:
    # GA effect on PNA50 dropped per supplement Table S1 model 42.1).
    # =====================================================================
    lcl_sulf            <- fixed(log(1))  ; label("Multiplicative anchor for CL_form(SULF) IIV (fixed at 1)") # placeholder so etalcl_sulf has a paired log-parameter
    lclbirth_sulf       <- log(0.1854)    ; label("PCM-SULF formation CL at birth at WT_BIRTH = 1.75 kg (L/h)") # Table 3 TVCLbirth(SULF) = 0.1854 L/h (RSE 4.2%)
    e_bwbirth_cl_sulf   <- 1.173          ; label("WT_BIRTH exponent on CL_birth(SULF) (unitless)")             # Table 3 theta_BWb(SULF) = 1.173 (RSE 5%)
    lclmax_sulf         <- log(0.3787)    ; label("PCM-SULF formation CL_max at WT = 1.75 kg (L/h)")             # Table 3 TVCLmax(SULF) = 0.3787 L/h (RSE 11%)
    e_wt_cl_sulf        <- fixed(0.738)   ; label("WT exponent on CL_max(SULF) (fixed)")                         # Table 3 theta_BWc(SULF) = 0.738 FIXED
    lpna50_sulf         <- log(25.92)     ; label("PCM-SULF PNA50 (days)")                                       # Table 3 TVPNA50(SULF) = 25.92 days (RSE 24%)
    e_ga_pna50_sulf     <- fixed(0)       ; label("GA exponent on PNA50(SULF) (fixed)")                          # Table 3 theta_GAPNA50(SULF) = 0 FIXED (supplement Table S1 model 42.1 -- PTNAnoGA)
    lhill_sulf          <- log(1.955)     ; label("PCM-SULF PTNA Hill coefficient (unitless)")                   # Table 3 Hill(SULF) = 1.955 (RSE 23%)

    # =====================================================================
    # Formation CL of PCM-OXI (PCM-cysmer = cysteine + mercapturate) via the
    # PTNA equation (no-GA variant). The canonical metabolite suffix
    # 'cysmer' (per vanRongen 2016 precedent and rxode2 registered
    # metabolites list) carries the paper's combined PCM-OXI species.
    # =====================================================================
    lcl_cysmer            <- fixed(log(1)) ; label("Multiplicative anchor for CL_form(cysmer/OXI) IIV (fixed at 1)") # placeholder so etalcl_cysmer has a paired log-parameter
    lclbirth_cysmer       <- log(0.02047)  ; label("PCM-cysmer formation CL at birth at WT_BIRTH = 1.75 kg (L/h)")   # Table 3 TVCLbirth(OXI) = 0.02047 L/h (RSE 13%)
    e_bwbirth_cl_cysmer   <- 0.978         ; label("WT_BIRTH exponent on CL_birth(cysmer/OXI) (unitless)")           # Table 3 theta_BWb(OXI) = 0.978 (RSE 17%)
    lclmax_cysmer         <- log(0.06377)  ; label("PCM-cysmer formation CL_max at WT = 1.75 kg (L/h)")              # Table 3 TVCLmax(OXI) = 0.06377 L/h (RSE 12%)
    e_wt_cl_cysmer        <- fixed(0.738)  ; label("WT exponent on CL_max(cysmer/OXI) (fixed)")                      # Table 3 theta_BWc(OXI) = 0.738 FIXED
    lpna50_cysmer         <- log(10.47)    ; label("PCM-cysmer PNA50 (days)")                                         # Table 3 TVPNA50(OXI) = 10.47 days (RSE 12%)
    e_ga_pna50_cysmer     <- fixed(0)      ; label("GA exponent on PNA50(cysmer/OXI) (fixed)")                        # Table 3 theta_GAPNA50(OXI) = 0 FIXED (supplement Table S1 model 53 -- OXI_PTNAnoGA)
    lhill_cysmer          <- log(2.972)    ; label("PCM-cysmer PTNA Hill coefficient (unitless)")                     # Table 3 Hill(OXI) = 2.972 (RSE 17%)

    # =====================================================================
    # Renal elimination as fractions of GFR (Wu 2024 PTNA equation).
    # PCM-GLU additionally carries an adult-only multiplicative correction
    # factor; PCM-SULF additionally carries a secretion CL term that
    # matures via a simplified PTNA form.
    # =====================================================================
    lfrn_pcm        <- log(0.08486)   ; label("Fraction of GFR for renal CL of unchanged PCM (unitless)")        # Table 3 f_PCM = 0.08486 (RSE 9.7%)
    lfrn_gluc       <- log(0.4475)    ; label("Fraction of GFR for renal CL of PCM-GLU (unitless)")              # Table 3 f_GLU = 0.4475 (RSE 6.5%)
    lf_gluc_adult   <- log(1.663)     ; label("Adult-only correction factor on PCM-GLU renal CL (unitless)")     # Table 3 f_GLU,adult = 1.663 (RSE 18%) -- applied for subjects >= 18 years
    lfrn_sulf       <- log(0.3059)    ; label("Fraction of GFR for renal CL of PCM-SULF (unitless)")             # Table 3 f_SULF = 0.3059 (RSE 12%)
    lfrn_cysmer     <- log(0.7393)    ; label("Fraction of GFR for renal CL of PCM-cysmer/OXI (unitless)")       # Table 3 f_OXI = 0.7393 (RSE 6.8%)

    # PCM-SULF renal secretion CL (PTNA-like; CL_birth and theta_BWb both
    # fixed to 0, Hill fixed to 1 -- paper Results "Improvement of the
    # Description of Renal Elimination Clearance of PCM-SULF"). The
    # secretion equation re-uses the GFR PTNA structure, so TVCLmax is
    # also reported in mL/min (same unit convention as Wu 2024 GFR)
    # despite the Table 3 column header listing L/h. Apply the same
    # mL/min -> L/h conversion (x0.06) so that f_SULF * GFR + secretion
    # in L/h reproduces the paper's "31% -> 163% of GFR" maturation
    # narrative (11.92 mL/min * 0.06 * (70/1.75)^0.738 ~= 10.9 L/h at
    # adult; 0.3059 * GFR_adult ~= 2.5 L/h gives total ~13.4 L/h ~=
    # 163% of GFR_adult ~= 8.2 L/h).
    lclmax_secr_sulf      <- log(11.92 * 0.06) ; label("PCM-SULF renal secretion CL_max at WT = 1.75 kg (L/h)")  # Table 3 TVCLmax (secretion) = 11.92 (RSE 20%); paper's Table 3 header says L/h but the secretion re-uses the GFR PTNA structure whose original publication (Wu 2024 Pharm Res 41:637-649) confirms mL/min, so we apply the same x0.06 conversion here so renal SULF = f_SULF * GFR + secretion is in L/h
    e_wt_cl_secr_sulf     <- fixed(0.738)      ; label("WT exponent on PCM-SULF secretion CL_max (fixed)")        # Table 3 theta_BWc (secretion) = 0.738 FIXED
    lpna50_secr_sulf      <- log(80.03)        ; label("PCM-SULF secretion PNA50 (days)")                          # Table 3 TVPNA50 (secretion) = 80.03 days (RSE 28%)
    e_ga_pna50_secr_sulf  <- -5.849            ; label("GA exponent on PNA50 for PCM-SULF secretion (unitless)")   # Table 3 theta_GAPNA50 (secretion) = -5.849 (RSE 11%)

    # =====================================================================
    # Metabolite volumes of distribution -- fractions of PCM central V_P.
    # =====================================================================
    lf_vol_gluc    <- log(0.3299)  ; label("Fraction of PCM V_P for PCM-GLU volume V_G (unitless)")           # Table 3 f_VG = 0.3299 (RSE 6.8%)
    lf_vol_sulf    <- log(0.3476)  ; label("Fraction of PCM V_P for PCM-SULF volume V_S (unitless)")          # Table 3 f_VS = 0.3476 (RSE 3.8%)
    lf_vol_cysmer  <- log(0.7048)  ; label("Fraction of PCM V_P for PCM-cysmer/OXI volume V_O (unitless)")    # Table 3 f_VO = 0.7048 (RSE 6.9%)

    # =====================================================================
    # GFR maturation (Wu 2024 PTNA -- all parameters fixed by paper). The
    # original publication confirms output is in mL/min despite the Wu 2025
    # equation 3 label of "L/h"; see Wu et al. Pharm Res 2024;41:637-649.
    # CL_birth = 1.26 mL/min and CL_max = 8.98 mL/min at 1.75 kg
    # birthweight / current weight. We pre-convert mL/min -> L/h via
    # x60/1000 = 0.06 and carry the result in log-space:
    #   1.26 * 0.06 = 0.0756 L/h
    #   8.98 * 0.06 = 0.5388 L/h
    # =====================================================================
    lclbirth_gfr        <- fixed(log(1.26 * 0.06))  ; label("GFR contribution at birth at WT_BIRTH = 1.75 kg (L/h)")  # Wu 2024 eq. 1; 1.26 mL/min converted to L/h
    lclmax_gfr          <- fixed(log(8.98 * 0.06))  ; label("GFR contribution at maturation at WT = 1.75 kg (L/h)")    # Wu 2024 eq. 1; 8.98 mL/min converted to L/h
    e_wt_gfr            <- fixed(0.738)             ; label("WT exponent on GFR_max (fixed)")                          # Wu 2024 eq. 1; (CW/1.75)^0.738
    lpna50_gfr          <- fixed(log(34))           ; label("GFR PNA50 (days)")                                         # Wu 2024 eq. 1; PNA50 = 34 days
    e_ga_pna50_gfr      <- fixed(-3.61)             ; label("GA exponent on PNA50 for GFR (unitless)")                  # Wu 2024 eq. 1; (GA/34)^-3.61
    lhill_gfr           <- fixed(log(1.03))         ; label("GFR PTNA Hill coefficient (unitless)")                      # Wu 2024 eq. 1; PNA^1.03

    # =====================================================================
    # Inter-individual variability (Table 3 reports IIV as omega% (CV%) +
    # shrinkage; omega^2 = log(1 + CV^2) for the lognormal parameterisation).
    # =====================================================================
    etalvp           ~ 0.07227   # IIV(V_P)            = 27.38%  -> log(1 + 0.2738^2) = 0.07227 (shrinkage 19.7%)
    etalka           ~ 0.34636   # IIV(Ka, rectal)     = 64.43%  -> log(1 + 0.6443^2) = 0.34636 (shrinkage 64.6%)
    etalcl_gluc      ~ 0.27084   # IIV(CL_form,GLU)    = 55.99%  -> log(1 + 0.5599^2) = 0.27084 (shrinkage 24.1%)
    etalcl_sulf      ~ 0.12150   # IIV(CL_form,SULF)   = 36.10%  -> log(1 + 0.3610^2) = 0.12150 (shrinkage 18.8%)
    etalcl_cysmer    ~ 0.39646   # IIV(CL_form,OXI)    = 70.49%  -> log(1 + 0.7049^2) = 0.39646 (shrinkage 25.6%)
    etalfrn_pcm      ~ 0.45914   # IIV(f_PCM)          = 78.20%  -> log(1 + 0.7820^2) = 0.45914 (shrinkage 50.5%)
    etalfrn_gluc     ~ 0.15469   # IIV(f_GLU)          = 40.89%  -> log(1 + 0.4089^2) = 0.15469 (shrinkage 48.2%)
    etalfrn_sulf     ~ 0.15967   # IIV(f_SULF + secr)  = 41.50%  -> log(1 + 0.4150^2) = 0.15967 (shrinkage 34.8%)
    etalfrn_cysmer   ~ 0.16438   # IIV(f_OXI)          = 42.38%  -> log(1 + 0.4238^2) = 0.16438 (shrinkage 48.6%)

    # =====================================================================
    # Residual error -- exponential model on log-transformed concentrations
    # (Methods page 1834 "an additive error model for log-transformed
    # concentrations was used"). Table 3 reports variance on the log-scale,
    # so SD = sqrt(variance); we use the linear-space proportional
    # equivalent inside model().
    # =====================================================================
    propSd          <- 0.2657   ; label("Plasma PCM proportional residual SD (fraction)")              # Table 3 plasma PCM variance = 0.0706 -> SD = sqrt(0.0706) = 0.2657
    propSd_gluc     <- 0.5101   ; label("Plasma PCM-GLU proportional residual SD (fraction)")          # Table 3 plasma GLU variance = 0.2602 -> SD = 0.5101
    propSd_sulf     <- 0.2866   ; label("Plasma PCM-SULF proportional residual SD (fraction)")         # Table 3 plasma SULF variance = 0.08214 -> SD = 0.2866
    propSd_cysmer   <- 0.3947   ; label("Plasma PCM-cysmer/OXI proportional residual SD (fraction)")   # Table 3 plasma OXI variance = 0.1558 -> SD = 0.3947
  })

  model({
    # -------------------------------------------------------------------
    # Covariate unit conversions -- the paper's PTNA equations are written
    # with PNA in days but the canonical PNA carries months.
    # -------------------------------------------------------------------
    PNA_days  <- PNA * 30.4375
    PNA_years <- PNA / 12
    is_adult  <- PNA_years >= 18

    # -------------------------------------------------------------------
    # GFR maturation (Wu 2024 PTNA equation), output in L/h after the
    # mL/min -> L/h pre-conversion absorbed into lclbirth_gfr / lclmax_gfr.
    # -------------------------------------------------------------------
    gfr_birth_term <- exp(lclbirth_gfr) * (WT_BIRTH / 1.75)
    gfr_max_term   <- exp(lclmax_gfr)   * (WT       / 1.75)^e_wt_gfr
    pna50_gfr      <- (GA / 34)^e_ga_pna50_gfr * exp(lpna50_gfr)
    hill_gfr       <- exp(lhill_gfr)
    matur_gfr      <- PNA_days^hill_gfr / (pna50_gfr^hill_gfr + PNA_days^hill_gfr)
    gfr            <- gfr_birth_term + (gfr_max_term - gfr_birth_term) * matur_gfr

    # -------------------------------------------------------------------
    # Formation CL of PCM-GLU (PTNA), with multiplicative IIV on the
    # entire formation-CL expression. The lcl_gluc anchor is fixed at
    # log(1) so exp(lcl_gluc + etalcl_gluc) = exp(etalcl_gluc).
    # -------------------------------------------------------------------
    clbirth_gluc <- exp(lclbirth_gluc) * (WT_BIRTH / 1.75)^e_bwbirth_cl_gluc
    clmax_gluc   <- exp(lclmax_gluc)   * (WT       / 1.75)^e_wt_cl_gluc
    pna50_gluc   <- (GA / 34)^e_ga_pna50_gluc * exp(lpna50_gluc)
    hill_gluc    <- exp(lhill_gluc)
    matur_gluc   <- PNA_days^hill_gluc / (pna50_gluc^hill_gluc + PNA_days^hill_gluc)
    cl_form_gluc <- exp(lcl_gluc + etalcl_gluc) * (clbirth_gluc + (clmax_gluc - clbirth_gluc) * matur_gluc)

    # -------------------------------------------------------------------
    # Formation CL of PCM-SULF (PTNA no-GA: GA effect = 0 FIXED so
    # (GA/34)^0 = 1 and pna50_sulf collapses to the estimated TVPNA50).
    # -------------------------------------------------------------------
    clbirth_sulf <- exp(lclbirth_sulf) * (WT_BIRTH / 1.75)^e_bwbirth_cl_sulf
    clmax_sulf   <- exp(lclmax_sulf)   * (WT       / 1.75)^e_wt_cl_sulf
    pna50_sulf   <- (GA / 34)^e_ga_pna50_sulf * exp(lpna50_sulf)
    hill_sulf    <- exp(lhill_sulf)
    matur_sulf   <- PNA_days^hill_sulf / (pna50_sulf^hill_sulf + PNA_days^hill_sulf)
    cl_form_sulf <- exp(lcl_sulf + etalcl_sulf) * (clbirth_sulf + (clmax_sulf - clbirth_sulf) * matur_sulf)

    # -------------------------------------------------------------------
    # Formation CL of PCM-cysmer/OXI (PTNA no-GA).
    # -------------------------------------------------------------------
    clbirth_cysmer <- exp(lclbirth_cysmer) * (WT_BIRTH / 1.75)^e_bwbirth_cl_cysmer
    clmax_cysmer   <- exp(lclmax_cysmer)   * (WT       / 1.75)^e_wt_cl_cysmer
    pna50_cysmer   <- (GA / 34)^e_ga_pna50_cysmer * exp(lpna50_cysmer)
    hill_cysmer    <- exp(lhill_cysmer)
    matur_cysmer   <- PNA_days^hill_cysmer / (pna50_cysmer^hill_cysmer + PNA_days^hill_cysmer)
    cl_form_cysmer <- exp(lcl_cysmer + etalcl_cysmer) * (clbirth_cysmer + (clmax_cysmer - clbirth_cysmer) * matur_cysmer)

    # -------------------------------------------------------------------
    # PCM-SULF additional renal secretion CL (PTNA-like, CL_birth = 0 and
    # Hill = 1 by simplification per paper Results).
    # -------------------------------------------------------------------
    clmax_secr_sulf <- exp(lclmax_secr_sulf) * (WT / 1.75)^e_wt_cl_secr_sulf
    pna50_secr_sulf <- (GA / 34)^e_ga_pna50_secr_sulf * exp(lpna50_secr_sulf)
    cl_secr_sulf    <- clmax_secr_sulf * PNA_days / (pna50_secr_sulf + PNA_days)

    # -------------------------------------------------------------------
    # Renal elimination CL of unchanged PCM and the three metabolites.
    # The adult-only correction multiplies f_GLU by f_GLU,adult for
    # subjects >= 18 years.
    # -------------------------------------------------------------------
    frn_pcm      <- exp(lfrn_pcm    + etalfrn_pcm)
    frn_gluc     <- exp(lfrn_gluc   + etalfrn_gluc) * (1 + (exp(lf_gluc_adult) - 1) * is_adult)
    frn_sulf     <- exp(lfrn_sulf   + etalfrn_sulf)
    frn_cysmer   <- exp(lfrn_cysmer + etalfrn_cysmer)

    cl_renal_pcm    <- frn_pcm    * gfr
    cl_renal_gluc   <- frn_gluc   * gfr
    cl_renal_sulf   <- frn_sulf   * gfr + cl_secr_sulf
    cl_renal_cysmer <- frn_cysmer * gfr

    # -------------------------------------------------------------------
    # PCM disposition individual parameters.
    # -------------------------------------------------------------------
    vc_pcm <- exp(lvp  + etalvp) * (WT / 1.08)^e_wt_vp
    vpp    <- exp(lvpp)          * (WT / 1.08)^e_wt_vpp
    q_pcm  <- exp(lq)            * (WT / 1.08)^e_wt_q

    # Metabolite volumes.
    vc_gluc   <- exp(lf_vol_gluc)   * vc_pcm
    vc_sulf   <- exp(lf_vol_sulf)   * vc_pcm
    vc_cysmer <- exp(lf_vol_cysmer) * vc_pcm

    # Rectal absorption.
    ka   <- exp(lka + etalka)
    F1   <- exp(lfdepot)
    Tlag <- exp(ltlag)

    # -------------------------------------------------------------------
    # Micro-constants for the ODE system.
    # -------------------------------------------------------------------
    kfg <- cl_form_gluc   / vc_pcm     # PCM -> central_gluc
    kfs <- cl_form_sulf   / vc_pcm     # PCM -> central_sulf
    kfc <- cl_form_cysmer / vc_pcm     # PCM -> central_cysmer
    krp <- cl_renal_pcm   / vc_pcm     # PCM -> renal
    k12 <- q_pcm / vc_pcm               # PCM central -> peripheral
    k21 <- q_pcm / vpp                  # PCM peripheral -> central
    keg <- cl_renal_gluc   / vc_gluc    # PCM-GLU    -> urine
    kes <- cl_renal_sulf   / vc_sulf    # PCM-SULF   -> urine
    kec <- cl_renal_cysmer / vc_cysmer  # PCM-cysmer -> urine

    # -------------------------------------------------------------------
    # ODE system. The rectal absorption depot feeds the PCM central
    # compartment; IV doses go directly to central. Each metabolite
    # accumulates via its formation rate and is eliminated renally.
    # -------------------------------------------------------------------
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - (kfg + kfs + kfc + krp + k12) * central + k21 * peripheral1
    d/dt(peripheral1)    <-  k12 * central - k21 * peripheral1
    d/dt(central_gluc)   <-  kfg * central - keg * central_gluc
    d/dt(central_sulf)   <-  kfs * central - kes * central_sulf
    d/dt(central_cysmer) <-  kfc * central - kec * central_cysmer

    # Rectal bioavailability and absorption lag (IV doses go to central
    # and are unaffected by these depot adjustments).
    f(depot)    <- F1
    alag(depot) <- Tlag

    # -------------------------------------------------------------------
    # Observation variables (mg/L PCM-equivalent for the three
    # metabolites, per Methods page 1833: "All metabolite concentrations
    # were expressed in PCM equivalents (mg/L) via conversion based on
    # molecular weights").
    # -------------------------------------------------------------------
    Cc        <- central        / vc_pcm
    Cc_gluc   <- central_gluc   / vc_gluc
    Cc_sulf   <- central_sulf   / vc_sulf
    Cc_cysmer <- central_cysmer / vc_cysmer

    Cc        ~ prop(propSd)
    Cc_gluc   ~ prop(propSd_gluc)
    Cc_sulf   ~ prop(propSd_sulf)
    Cc_cysmer ~ prop(propSd_cysmer)
  })
}
