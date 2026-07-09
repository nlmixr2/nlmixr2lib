Han_2015_bevacizumab <- function() {
  description <- "Two-compartment population PK model for IV bevacizumab in paediatric cancer patients (0.5-21 years) with fixed allometric body-weight scaling on all four disposition parameters, plus additional covariate effects of sex, baseline serum albumin, and primary-CNS-tumour-vs-sarcoma indication on CL and V1 (Han 2015, Table 2 final model)"
  reference <- "Han K, Peyret T, Quartino A, Gosselin NH, Gururangan S, Casanova M, Merks JHM, Massimino M, Grill J, Daw NC, Navid F, Jin J, Allison DE. Bevacizumab dosing strategy in paediatric cancer patients based on population pharmacokinetic analysis with external validation. Br J Clin Pharmacol. 2016;81(1):148-160. doi:10.1111/bcp.12778"
  vignette <- "Han_2015_bevacizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight; used for allometric power scaling on all four disposition parameters with reference weight 44 kg (paediatric cohort median, Han 2015 Figure 3 caption).",
      source_name        = "BWT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline serum albumin; power scaling on CL with reference 39 g/L (paediatric cohort median, Han 2015 Figure 3 caption). Paper reports ALBU range 24-52 g/L in the model-building cohort (Table 1).",
      source_name        = "ALBU"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Han 2015 reports the effect using a male-indicator convention: CL is multiplied by 1.11 and V1 by 1.14 when the patient is male. The paper's typical values (CL = 9.90 mL/h, V1 = 2850 mL) correspond to a female patient (paper's Figure 3 base = 44-kg female paediatric patient with albumin 39 g/L and sarcomas). To keep those typicals visible verbatim, the model file applies the male-vs-female contrast on the (1 - SEXF) axis via 1.11^(1 - SEXF) on CL and 1.14^(1 - SEXF) on V1.",
      source_name        = "SEX (male indicator; derive SEXF = 1 - SEX_MALE)"
    ),
    TUMTP_CNS_PRIM = list(
      description        = "Primary CNS tumour vs non-CNS-primary tumour indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-primary-CNS tumour; in Han 2015 = sarcomas)",
      notes              = "Han 2015 Table 1 splits the paediatric cohort into primary CNS tumours (study AVF3842s; n = 76) and sarcomas (studies AVF2771s + AVF4117s + BO20924; n = 76). The primary CNS tumour umbrella pools glioma / medulloblastoma / ependymoma / brainstem-glioma / atypical-teratoid-rhabdoid histologies (specific breakdown not enumerated in the paper). Paper's typical CL and V1 are anchored to the sarcoma reference; multiplicative factors 0.725 on CL and 0.854 on V1 apply when TUMTP_CNS_PRIM = 1.",
      source_name        = "Primary CNS tumour indicator (paper narrative)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 152L,
    n_studies      = 4L,
    age_range      = "0.5-21 years",
    age_median     = "10.8 years",
    weight_range   = "5.94-125 kg",
    weight_median  = "43.8 kg",
    sex_female_pct = 46.1,
    race_ethnicity = "Not reported (multi-region study; race distribution not tabulated in Han 2015)",
    disease_state  = "Paediatric cancer patients: refractory sarcomas (AVF2771s), primary CNS tumours (AVF3842s), newly diagnosed osteosarcoma (AVF4117s), and metastatic soft-tissue sarcoma (BO20924 interim data). Cohort split evenly between primary CNS tumours (n = 76) and sarcomas (n = 76).",
    dose_range     = "5-15 mg/kg IV infusion over 30-90 minutes; regimens Q2W (5, 10, 15 mg/kg) or Q3W (7.5, 15 mg/kg).",
    regions        = "Multi-region: US, Europe (Italy, Netherlands, France)",
    n_observations = "1427 quantifiable bevacizumab serum concentrations (of 1592 collected in the model-building population).",
    co_medication  = "Bevacizumab given as single agent in AVF2771s and in combination with chemotherapy in all other studies.",
    notes          = "Table 1 model-building cohort. Nine infants (0-2 years) included; 12 patients under 3 years across the full n = 232 cohort. Serum LLOQ = 78 ng/mL (ELISA); values below LLOQ omitted from the fit. Median half-life 19.6 days (range 9-78 days) reported in Discussion. Model-building N = 152; external validation N = 80 (studies BO25041 primary CNS + BO20924 remaining sarcoma data).",
    external_validation = "N = 80 paediatric patients (544 concentrations) from BO25041 (primary CNS tumours) and remaining BO20924 (metastatic soft-tissue sarcoma). Mean prediction error for concentrations 3.54%, for CL -1.84%, for V1 -0.06% (Han 2015 Results, External validation)."
  )

  ini({
    # Structural parameters -- typical values for the paediatric cohort median
    # patient (44 kg female paediatric patient with albumin 39 g/L and sarcomas;
    # Han 2015 Figure 3 caption). Reported in mL/h and mL by the paper Table 2;
    # converted here to L/day and L to match the day/mg/ug/mL time-and-unit
    # convention used across the other bevacizumab entries in nlmixr2lib.
    #   CL = 9.90 mL/h  * 24 / 1000 = 0.2376 L/day
    #   V1 = 2850 mL          / 1000 = 2.85   L
    #   Q  = 28.0  mL/h * 24 / 1000 = 0.6720 L/day
    #   V2 = 2564 mL          / 1000 = 2.564  L
    lcl <- log(0.2376); label("Clearance for median 44 kg female paediatric sarcoma patient (CL, L/day)")  # Table 2 row CL = 9.90 mL/h (RSE 4.1%)
    lvc <- log(2.85);   label("Central volume of distribution for median patient (V1, L)")                 # Table 2 row V1 = 2850 mL (RSE 3.0%)
    lq  <- log(0.6720); label("Inter-compartmental clearance for median patient (Q, L/day)")               # Table 2 row Q  = 28.0 mL/h (RSE 10.4%)
    lvp <- log(2.564);  label("Peripheral volume of distribution for median patient (V2, L)")              # Table 2 row V2 = 2564 mL (RSE 5.8%)

    # Fixed allometric exponents on body weight (BWT, reference 44 kg = cohort
    # median). Han 2015 Table 2: BWT-on-CL and BWT-on-Q are fixed at the
    # theoretical 0.75; BWT-on-V1 and BWT-on-V2 were estimated in the base
    # model (0.701 and 0.766) and then held at those values in the final model
    # (Table 2 asterisk footnote: 'value estimated for the base model';
    # Bootstrapping-median column reads 'Fixed').
    e_wt_cl <- fixed(0.75);  label("Allometric exponent on log(BWT/44) for CL (unitless)")   # Table 2 row 'BWT on CL' = 0.75 Fixed
    e_wt_vc <- fixed(0.701); label("Allometric exponent on log(BWT/44) for V1 (unitless)")   # Table 2 row 'BWT on V1' = 0.701 Fixed (base-model estimate)
    e_wt_q  <- fixed(0.75);  label("Allometric exponent on log(BWT/44) for Q  (unitless)")   # Table 2 row 'BWT on Q'  = 0.75 Fixed
    e_wt_vp <- fixed(0.766); label("Allometric exponent on log(BWT/44) for V2 (unitless)")   # Table 2 row 'BWT on V2' = 0.766 Fixed (base-model estimate)

    # Covariate effects on CL (paper Table 2 rows below the allometric block).
    # ALB enters as a power scaling with reference 39 g/L; SEX and primary-CNS
    # indicator enter as multiplicative factors on the (1 - SEXF) and
    # TUMTP_CNS_PRIM axes so the paper's typical values stay verbatim.
    e_alb_cl      <- -0.300; label("Power exponent on log(ALB/39) for CL (unitless)")                     # Table 2 row 'ALBU on CL'                    = -0.300 (RSE 50.6%)
    e_sex_cl      <-  1.11;  label("Multiplicative male-vs-female factor on CL (unitless)")               # Table 2 row 'Male on CL'                    =  1.11  (RSE 4.1%)
    e_cns_prim_cl <-  0.725; label("Multiplicative primary-CNS-vs-sarcoma factor on CL (unitless)")       # Table 2 row 'Primary CNS tumour on CL'      =  0.725 (RSE 4.3%)

    # Covariate effects on V1
    e_sex_vc      <-  1.14;  label("Multiplicative male-vs-female factor on V1 (unitless)")               # Table 2 row 'Male on V1'                    =  1.14  (RSE 3.5%)
    e_cns_prim_vc <-  0.854; label("Multiplicative primary-CNS-vs-sarcoma factor on V1 (unitless)")       # Table 2 row 'Primary CNS tumour on V1'      =  0.854 (RSE 3.7%)

    # Full-block IIV on CL, V1, and V2 (Han 2015 Table 2 lower block).
    # Diagonal entries are omega^2 taken as CV^2 (NONMEM approximation matching
    # how the paper reports IIV(%) = 100 * sqrt(omega^2)); the raw off-diagonal
    # covariances are copied verbatim from the paper's Var. IIVs rows.
    #   omega^2 CL = 0.214^2 = 0.045796   Table 2 IIV CL 21.4 percent (RSE 7.4)
    #   omega^2 V1 = 0.176^2 = 0.030976   Table 2 IIV V1 17.6 percent (RSE 9.2)
    #   omega^2 V2 = 0.580^2 = 0.336400   Table 2 IIV V2 58.0 percent (RSE 13.3)
    #   cov(CL,V1) = 0.0248               Table 2 Var.IIVs CL and V1 (RSE 17.6)
    #   cov(CL,V2) = 0.0228               Table 2 Var.IIVs CL and V2 (RSE 75.2)
    #   cov(V1,V2) = 0.0248               Table 2 Var.IIVs V1 and V2 (RSE 51.0)
    # Implied correlations approximately r(CL,V1)=0.66, r(CL,V2)=0.18, r(V1,V2)=0.24.
    etalcl + etalvc + etalvp ~ c(
      0.045796,
      0.0248,   0.030976,
      0.0228,   0.0248,   0.336400
    )

    # Combined additive + proportional residual error (Han 2015 Table 2 bottom).
    # Additive term reported in ug/mL = mg/L (matches the units$concentration
    # declaration).
    propSd <- 0.139; label("Proportional residual error (fraction)")           # Table 2 row 'Prop. error (%)'          = 13.9  (RSE  6.9%)
    addSd  <- 3.06;  label("Additive residual error (ug/mL)")                  # Table 2 row 'Add. error (ug/mL)'      = 3.06  (RSE 25.1%)
  })
  model({
    # Individual PK parameters. Allometric power scaling with reference 44 kg
    # (paediatric cohort median). Covariate effects on CL and V1 exactly
    # reproduce the equations in Han 2015 text (paragraph following Table 2):
    #   CL_i = 9.90 * (BWT/44)^0.75 * (ALBU/39)^(-0.300) * 1.11^(male) * 0.725^(CNS)
    #   V1_i = 2850 * (BWT/44)^0.701                       * 1.14^(male) * 0.854^(CNS)
    #   Q_i  = 28.0 * (BWT/44)^0.75
    #   V2_i = 2564 * (BWT/44)^0.766
    # SEXF (canonical) = 1 for female / 0 for male; the paper's SEX column
    # is a male indicator, so (1 - SEXF) recovers the paper's convention.
    cl <- exp(lcl + etalcl) * (WT / 44)^e_wt_cl *
          (ALB / 39)^e_alb_cl *
          e_sex_cl^(1 - SEXF) *
          e_cns_prim_cl^TUMTP_CNS_PRIM
    vc <- exp(lvc + etalvc) * (WT / 44)^e_wt_vc *
          e_sex_vc^(1 - SEXF) *
          e_cns_prim_vc^TUMTP_CNS_PRIM
    q  <- exp(lq)           * (WT / 44)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 44)^e_wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV infusion enters `central` directly; no depot compartment. The
    # infusion rate / duration is supplied in the event dataset.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
