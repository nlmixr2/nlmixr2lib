Solms_2019_rivaroxaban <- function() {
  description <- "Population PK/PT model for oral rivaroxaban with a direct linear prothrombin-time response (Solms 2019; 1,779 patients with time-matched PK and PT samples from four phase II VTE-prevention / VTE-treatment trials and the ROCKET AF phase III PK substudy). The PK layer is the integrated popPK model of Willmann 2018a with every PK parameter fixed to its final estimate: one-compartment disposition with first-order absorption and elimination, dose-dependent relative bioavailability, Tietz-truncated creatinine clearance and body weight on CL/F, five multiplicative comedication factors and a study / indication factor on CL/F, and body weight, age and sex on V/F. Prothrombin time (PT, Neoplastine reagent) rises linearly with the individual predicted concentration, PT = BASE + SLOPE * Cc, with exponential interindividual variability on BASE and SLOPE and a proportional residual error on PT."
  reference <- paste(
    "Solms A, Frede M, Berkowitz SD, Hermanowski-Vosatka A, Kubitza D,",
    "Mueck W, Spiro TE, Willmann S, Yan X, Zhang L, Garmann D. Enhancing the Quality of Rivaroxaban Exposure Estimates Using",
    "Prothrombin Time in the Absence of Pharmacokinetic Sampling.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(11):805-814.",
    "doi:10.1002/psp4.12444"
  )
  vignette <- "Solms_2019_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # What each ODE state holds, in what amount units, in what biological matrix.
  compartmentData <- list(
    depot = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation, then Tietz-truncated inside model()",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw (untruncated) Cockcroft-Gault creatinine clearance in mL/min, not",
        "BSA-normalized. model() caps it at 140 * BSA / 1.73 mL/min, exactly as",
        "the popPK/PT control stream does ('CRCL_MAX = 140*BSA2/1.73'), and it",
        "then enters CL/F as (CRCL_truncated / 93)^0.405575. The PK layer is",
        "the Willmann 2018a integrated popPK model with all PK parameters",
        "fixed (Solms 2019 Supplementary Model Code, popPK/PT control stream)."
      ),
      source_name = "CRCL"
    ),
    BSA = list(
      description = "Body surface area, used only to set the Tietz truncation ceiling on creatinine clearance",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not a covariate on any PK parameter; it appears only in the Tietz",
        "creatinine-clearance ceiling 140 * BSA / 1.73 mL/min. The control",
        "stream imputes a missing BSA as 1.92 m^2 ('IF (BSA.LT.0) BSA2 =",
        "1.92'); supply a real BSA where available."
      ),
      source_name = "BSA"
    ),
    WT = list(
      description = "Total body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as (WT / 81)^-0.277974 and V/F as (WT / 81)^0.21584",
        "(popPK/PT control stream THETA(13) and THETA(5), both FIX). The",
        "negative exponent on CL/F refines the weight dependence already",
        "carried by the Cockcroft-Gault creatinine clearance (Willmann 2018a)."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters V/F only, as (AGE / 61)^-0.189009 (popPK/PT control stream THETA(6) FIX).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male; V/F multiplier 1 by definition).",
      notes = paste(
        "Multiplicative shift on V/F of 0.889284^SEXF (popPK/PT control stream",
        "THETA(7) FIX). The control stream writes THETA(7)**(SEX-1) on a",
        "1 = male / 2 = female source coding, so SEXF = SEX - 1."
      ),
      source_name = "SEX (1 = male, 2 = female)"
    ),
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant P-glycoprotein inhibitor).",
      notes = paste(
        "Multiplicative CL/F factor 0.965478^CONMED_PGP_INH (THETA(8) FIX).",
        "The five comedication factors multiply (control stream COMM).",
        "Solms 2019 Methods: comedication factors were assumed constant over",
        "the study period."
      ),
      source_name = "PGPINH"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A4 inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant strong CYP3A4 inhibitor).",
      notes = paste(
        "Multiplicative CL/F factor 0.977924^CONMED_CYP3A4_INH_STRONG",
        "(THETA(9) FIX). Willmann 2018a cautions that this coefficient rests on",
        "6 patients and must not be used to predict a dedicated strong-inhibitor",
        "interaction."
      ),
      source_name = "SCYP3A4INH"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Concomitant moderate ('medium') CYP3A4 inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant moderate CYP3A4 inhibitor).",
      notes = "Multiplicative CL/F factor 0.863422^CONMED_CYP3A4_INH_MOD (THETA(10) FIX; 'mCYPi' in the control stream).",
      source_name = "MCYP3A4INH"
    ),
    CONMED_CYP3A4_INH_WEAK = list(
      description = "Concomitant weak CYP3A4 inhibitor at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant weak CYP3A4 inhibitor).",
      notes = "Multiplicative CL/F factor 0.939123^CONMED_CYP3A4_INH_WEAK (THETA(11) FIX).",
      source_name = "WCYP3A4INH"
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer at the record",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant CYP3A4 inducer).",
      notes = "Multiplicative CL/F factor 1.30299^CONMED_CYP3A4_IND (THETA(12) FIX), a 30% increase in apparent clearance.",
      source_name = "CYP3A4IND"
    ),
    DIS_AF = list(
      description = "Nonvalvular atrial fibrillation cohort indicator (ROCKET AF, study 3001)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the VTE-treatment reference cohort, ODIXa-DVT and EINSTEIN-DVT, CL/F factor 1).",
      notes = paste(
        "Multiplicative CL/F factor 0.84858^DIS_AF (THETA(14) FIX; control",
        "stream 'IF (STU.EQ.3001) INDEFF = THETA(14)'). Exactly one of DIS_AF,",
        "DIS_ACS and DIS_VTE_P may be 1; all three 0 selects VTE treatment."
      ),
      source_name = "STU = 3001"
    ),
    DIS_ACS = list(
      description = "Acute coronary syndrome cohort indicator (study 2001)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the VTE-treatment reference cohort, CL/F factor 1).",
      notes = paste(
        "Multiplicative CL/F factor 1.13897^DIS_ACS (THETA(15) FIX). Carried",
        "from the integrated popPK model; no ACS patients contributed to the",
        "Solms 2019 PK-PT data set (PT was not measured in ATLAS ACS 2-TIMI 51)."
      ),
      source_name = "STU = 2001"
    ),
    DIS_VTE_P = list(
      description = "VTE-prevention (post elective hip or knee replacement thromboprophylaxis) cohort indicator (studies 10944, 10945, 11527)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the VTE-treatment reference cohort, CL/F factor 1).",
      notes = paste(
        "Time-varying CL/F factor: 1.03711 at or before 72 h after the first",
        "dose and 1.2848 thereafter (THETA(16) and THETA(17) FIX; control",
        "stream 'IF (TIME.LE.72) INDEFF = THETA(16)'). model() keys the split",
        "on solver time t, so the first dose must be placed at t = 0."
      ),
      source_name = "STU in (10944, 10945, 11527)"
    ),
    DOSE_RIVAROXABAN_MG = list(
      description = "Administered rivaroxaban dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The dose of a single administration (not the daily total), supplied",
        "as a data column so the dose-dependent relative bioavailability can be",
        "evaluated inside model(), as the control stream does with its DOSE",
        "column. Set it equal to amt on each dose record."
      ),
      source_name = "DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1779,
    n_studies = 5,
    age_range = "adults (demographics of the PK-PT data set are not tabulated in Solms 2019)",
    race_ethnicity = "not reported",
    disease_state = paste(
      "VTE prevention after elective hip or knee replacement (ODIXa-HIP2,",
      "ODIXa-KNEE; 831 patients), VTE treatment of acute symptomatic DVT",
      "(ODIXa-DVT, EINSTEIN-DVT dose-ranging; 788 patients) and nonvalvular",
      "atrial fibrillation (ROCKET AF PK substudy; 161 patients)."
    ),
    dose_range = "2.5-30 mg twice daily and 20-40 mg once daily (phase II dose ranging); 20 mg (15 mg at CrCL 30-49 mL/min) once daily in ROCKET AF",
    regions = "global phase II / III program",
    notes = paste(
      "Patients with valid time-matched PK and PT samples, Table 1: ODIXa-HIP2",
      "508, ODIXa-KNEE 322, ODIXa-DVT 437, EINSTEIN-DVT 351, ROCKET AF 161",
      "(sum 1,779; the Results text gives 831 + 788 + 161 = 1,780). PT was",
      "measured centrally with a rivaroxaban-sensitive thromboplastin (STA",
      "Neoplastine or STA Neoplastine CI Plus). PT values > 50 s, PT more than",
      "30 h after the last dose, PT within 72 h of surgery (VTE prevention) and",
      "PT in the first 3 weeks of VTE treatment were excluded. Estimation in",
      "NONMEM 7.3, FOCE with interaction, PK parameters fixed."
    )
  )

  ini({
    # PK layer -- the Willmann 2018a integrated popPK model. Every value is the
    # FIX value of the popPK/PT control stream in Solms 2019 Supplementary
    # Model Code ('In the following model file the popPK parameters are fixed
    # to parameter estimates of final popPK model'). Typical values apply at
    # CRCL 93 mL/min, WT 81 kg, AGE 61 years, male, no comedication, VTE
    # treatment.
    lka <- fixed(log(0.820743))
    label("Typical first-order absorption rate constant (1/h)") # popPK/PT $THETA 1 TKA 0.820743 FIX
    lcl <- fixed(log(6.58214))
    label("Typical apparent clearance CL/F (L/h)") # popPK/PT $THETA 2 TCL 6.58214 FIX
    lvc <- fixed(log(62.5172))
    label("Typical apparent central volume V/F (L)") # popPK/PT $THETA 3 TV 62.5172 FIX

    # Dose-dependent relative bioavailability,
    # F1 = F1MIN + (F1MAX - F1MIN) * exp(-log(2) / D50 * DOSE).
    fdepot_min <- fixed(0.58955)
    label("Asymptotic minimum relative bioavailability at high dose (fraction)") # popPK/PT $THETA 19 F1MIN 0.58955 FIX
    fdepot_max <- fixed(1.25)
    label("Asymptotic maximum relative bioavailability as dose approaches zero (fraction)") # popPK/PT $THETA 18 F1MAX 1.25 FIX
    led50 <- fixed(log(14.3671))
    label("Dose at which relative bioavailability has fallen halfway from Fmax to Fmin, D50 (mg)") # popPK/PT $THETA 20 D50 14.3671 FIX

    # Covariate effects on CL/F.
    e_crcl_cl <- fixed(0.405575)
    label("Power exponent on (Tietz-truncated CRCL / 93 mL/min) for CL/F (unitless)") # popPK/PT $THETA 4 CRCL_CL 0.405575 FIX
    e_wt_cl <- fixed(-0.277974)
    label("Power exponent on (WT / 81 kg) for CL/F (unitless)") # popPK/PT $THETA 13 CL~WT -0.277974 FIX
    e_pgp_inh_cl <- fixed(0.965478)
    label("Multiplicative CL/F factor for concomitant P-glycoprotein inhibitor (fraction)") # popPK/PT $THETA 8 CL~PGP 0.965478 FIX
    e_cyp3a4_inh_strong_cl <- fixed(0.977924)
    label("Multiplicative CL/F factor for concomitant strong CYP3A4 inhibitor (fraction)") # popPK/PT $THETA 9 CL~sCYPi 0.977924 FIX
    e_cyp3a4_inh_mod_cl <- fixed(0.863422)
    label("Multiplicative CL/F factor for concomitant moderate CYP3A4 inhibitor (fraction)") # popPK/PT $THETA 10 CL~mCYPi 0.863422 FIX
    e_cyp3a4_inh_weak_cl <- fixed(0.939123)
    label("Multiplicative CL/F factor for concomitant weak CYP3A4 inhibitor (fraction)") # popPK/PT $THETA 11 CL~wCYPi 0.939123 FIX
    e_cyp3a4_ind_cl <- fixed(1.30299)
    label("Multiplicative CL/F factor for concomitant CYP3A4 inducer (fraction)") # popPK/PT $THETA 12 CL~CYPind_CL 1.30299 FIX
    e_af_cl <- fixed(0.84858)
    label("Multiplicative CL/F factor for the atrial fibrillation cohort vs VTE treatment (fraction)") # popPK/PT $THETA 14 CL~AF 0.84858 FIX
    e_acs_cl <- fixed(1.13897)
    label("Multiplicative CL/F factor for the acute coronary syndrome cohort vs VTE treatment (fraction)") # popPK/PT $THETA 15 CL~ACS 1.13897 FIX
    e_vte_p_le72_cl <- fixed(1.03711)
    label("Multiplicative CL/F factor for VTE prevention at or before 72 h after first dose, vs VTE treatment (fraction)") # popPK/PT $THETA 16 VTE_LE3days 1.03711 FIX
    e_vte_p_gt72_cl <- fixed(1.2848)
    label("Multiplicative CL/F factor for VTE prevention after 72 h from first dose, vs VTE treatment (fraction)") # popPK/PT $THETA 17 VTE_GT3days 1.2848 FIX

    # Covariate effects on V/F.
    e_wt_vc <- fixed(0.21584)
    label("Power exponent on (WT / 81 kg) for V/F (unitless)") # popPK/PT $THETA 5 V~WT 0.21584 FIX
    e_age_vc <- fixed(-0.189009)
    label("Power exponent on (AGE / 61 years) for V/F (unitless)") # popPK/PT $THETA 6 V~AGE -0.189009 FIX
    e_sexf_vc <- fixed(0.889284)
    label("Multiplicative V/F factor for female sex vs male (fraction)") # popPK/PT $THETA 7 V~SEX 0.889284 FIX

    # Linear PK/PT relationship, PT = BASE + SLOPE * Cc. Final estimates from
    # the Results ('The estimated BASE PT value was 13.1 seconds ... and the
    # SLOPE value was 0.031 seconds/ug/L') and Supplementary Figure S2
    # (RSE 0.27% and 1.11%).
    lrbase <- log(13.1)
    label("Baseline prothrombin time before rivaroxaban (s)") # Results / Figure S2: BASE 13.1 s (RSE 0.27%)
    lslope <- log(0.031)
    label("Linear slope of prothrombin time on rivaroxaban concentration (s per ug/L)") # Results / Figure S2: SLOPE 0.031 s/ug/L (RSE 1.11%)

    # Interindividual variability. The PK etas are FIX in the popPK/PT
    # control stream; a FIX anywhere in a NONMEM $OMEGA BLOCK fixes the whole
    # block, so the CL/F-V/F covariance and V/F variance are held too. The
    # values equal the Willmann 2018a final estimates.
    etalka ~ fixed(0.628188) # popPK/PT $OMEGA 1 IIV_KA 0.628188 FIX
    etalcl + etalvc ~ fixed(
      0.166888,
      0.067385, 0.039115
    ) # popPK/PT $OMEGA BLOCK(2) IIV_CL 0.166888 FIX, covariance 0.067385, IIV_V 0.039115
    etalrbase ~ 0.006 # Results: IIV on BASE 0.006 (SD 0.0757, CV 7.76%); Figure S2 RSE 13%
    etalslope ~ 0.136 # Results: IIV on SLOPE 0.136 (SD 0.369, CV 38.2%); Figure S2 RSE 9.5%

    # Residual error. PK: proportional ('Y = IPRED + W*EPS(1)', W = IPRED). The
    # PK $SIGMA carries no FIX flag, so it was re-estimated with the PT data,
    # but its final value is not reported; the control-stream value 0.203389
    # (the Willmann 2018a final variance, sqrt = 0.450555) is used.
    # PT: proportional, 'PT(t) = (PT0 + PT1 Conc(t)) * (1 + eps)' (Figure S2
    # legend), final CV 11.4%.
    propSd <- 0.450555
    label("Proportional residual SD on rivaroxaban concentration (fraction)") # popPK/PT $SIGMA 1 RUV_PK_prop 0.203389, sqrt = 0.450555
    propSd_PT <- 0.114
    label("Proportional residual SD on prothrombin time (fraction)") # Figure S2: residual unexplained error CV 11.4% (RSE 12.8%)
  })

  model({
    # Tietz truncation of creatinine clearance (control stream CRCL_MAX).
    crclCap <- 140 * BSA / 1.73
    crclTrunc <- CRCL * (CRCL <= crclCap) + crclCap * (CRCL > crclCap)

    # Comedication factor: the five indicators multiply (control stream COMM);
    # each term is 1 when its indicator is 0.
    comed <- e_pgp_inh_cl^CONMED_PGP_INH *
      e_cyp3a4_inh_strong_cl^CONMED_CYP3A4_INH_STRONG *
      e_cyp3a4_inh_mod_cl^CONMED_CYP3A4_INH_MOD *
      e_cyp3a4_inh_weak_cl^CONMED_CYP3A4_INH_WEAK *
      e_cyp3a4_ind_cl^CONMED_CYP3A4_IND

    # Study / indication factor (control stream INDEFF). VTE treatment is the
    # reference; the VTE-prevention factor switches at 72 h after the first
    # dose (NONMEM TIME), so t must be referenced to the first dose.
    vteFactor <- e_vte_p_le72_cl +
      (e_vte_p_gt72_cl - e_vte_p_le72_cl) * (t > 72)
    study <- 1 +
      DIS_AF * (e_af_cl - 1) +
      DIS_ACS * (e_acs_cl - 1) +
      DIS_VTE_P * (vteFactor - 1)

    ka <- exp(lka + etalka)
    ed50 <- exp(led50)
    cl <- exp(lcl + etalcl) *
      (crclTrunc / 93)^e_crcl_cl *
      (WT / 81)^e_wt_cl *
      comed *
      study
    vc <- exp(lvc + etalvc) *
      (WT / 81)^e_wt_vc *
      (AGE / 61)^e_age_vc *
      e_sexf_vc^SEXF

    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # Dose-dependent relative bioavailability (control stream F1).
    f(depot) <- fdepot_min +
      (fdepot_max - fdepot_min) *
        exp(-log(2) / ed50 * DOSE_RIVAROXABAN_MG)

    # mg / L x 1000 = ug/L; the control stream folds the same factor into F1.
    Cc <- central / vc * 1000

    # Direct linear PK/PT relationship ('IPRED = BASE+SLOPE*F' for DVID 2),
    # driven by the error-free individual concentration.
    rbase <- exp(lrbase + etalrbase)
    slope <- exp(lslope + etalslope)
    PT <- rbase + slope * Cc

    Cc ~ prop(propSd)
    PT ~ prop(propSd_PT)
  })
}
