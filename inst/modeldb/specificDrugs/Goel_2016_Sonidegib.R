Goel_2016_Sonidegib <- function() {
  description <- "Two-compartment population PK model for sonidegib (LDE225) in healthy subjects and patients with advanced solid tumors with first-order absorption, lag time, linear elimination, and dose-dependent bioavailability (Goel 2016)"
  reference   <- "Goel V, Hurh E, Stein A, Nedelman J, Zhou J, Chiparus O, Huang P-H, Gogov S, Sellami D. Population pharmacokinetics of sonidegib (LDE225), an oral inhibitor of hedgehog pathway signaling, in healthy subjects and in patients with advanced solid tumors. Cancer Chemother Pharmacol. 2016;77(4):745-755. doi:10.1007/s00280-016-2982-1"
  vignette    <- "Goel_2016_Sonidegib"
  units       <- list(time = "hr", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F with reference 58 years (median pooled cohort age, Goel 2016 Table 1).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F (exponent -0.093) and on Vc/F (exponent 0.731); the Vp/F effect is fixed equal to the Vc/F effect per the paper. Reference 73 kg (median pooled cohort weight, Goel 2016 Table 1).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effect on CL/F: 0.846^SEXF. Goel 2016 Table 2 reports a 15% lower CL/F in females.",
      source_name        = "Female"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F with reference 93.3 mL/min (median pooled cohort CrCL, Goel 2016 Table 1). Goel 2016 reports CrCL in mL/min (not BSA-normalized to mL/min/1.73 m^2). The exponent of -0.043 is small; renal function did not reach clinical / statistical relevance in the pooled analysis.",
      source_name        = "CRCL"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F (exponent -1.44), Vc/F (exponent -2.33), and Vp/F (exponent -0.087). Reference 43 g/L (median pooled cohort albumin, Goel 2016 Table 1).",
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Baseline alanine aminotransferase, expressed as a ratio to the upper limit of normal (ULN)",
      units              = "fraction of ULN",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Goel 2016 Methods: 'The values of total bilirubin and ALT levels were normalized to the upper limit of normal (ULN).' The pooled-cohort median normalized ALT (column ALTN in the source NONMEM dataset) was 0.42 (Table 1), used here as the reference value. Values supplied in this canonical column should likewise be ratios (raw ALT in U/L divided by the assay's ULN), not raw U/L. Power effect on CL/F (exponent 0.120).",
      source_name        = "ALTN"
    ),
    TBILI = list(
      description        = "Baseline total bilirubin, expressed as a ratio to the upper limit of normal (ULN)",
      units              = "fraction of ULN",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Goel 2016 Methods: 'The values of total bilirubin and ALT levels were normalized to the upper limit of normal (ULN).' The pooled-cohort median normalized bilirubin (column BILN in the source NONMEM dataset) was 0.38 (Table 1), used here as the reference value. Values supplied in this canonical column should be ratios (raw bilirubin in mg/dL or umol/L divided by the assay's ULN), not raw concentrations. Power effect on CL/F (exponent 0.064).",
      source_name        = "BILN"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-vs-Western enrollment indicator (Goel 2016 ethnicity covariate)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subjects enrolled in Western countries)",
      notes              = "Multiplicative effect on CL/F: 0.905^RACE_JAPANESE. Goel 2016 collapses 'Japanese vs Western enrollment' into a single ethnicity binary; the canonical RACE_JAPANESE column is reused with that paper-specific operational definition.",
      source_name        = "JPN"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patient with advanced solid tumor or basal cell carcinoma)",
      notes              = "Multiplicative effect on CL/F: 2.96^DIS_HEALTHY. Healthy volunteers had ~3-fold higher apparent CL/F than cancer patients in Goel 2016. The reference complement is the pooled cancer-patient cohort (X2101 + X1101 + A2201).",
      source_name        = "HV"
    ),
    CONMED_PPI = list(
      description        = "Significant proton-pump-inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no significant CONMED_PPI use)",
      notes              = "Multiplicative effect on F: 0.696^CONMED_PPI (~30% lower bioavailability). Goel 2016 Methods defines 'significant' use as duration of CONMED_PPI use exceeding 80% of the PK assessment phase (i.e., at least 80% of the time between the first and last PK sample for the subject).",
      source_name        = "CONMED_PPI"
    ),
    CONMED_H2RA = list(
      description        = "Significant H2-receptor-antagonist coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no significant CONMED_H2RA use)",
      notes              = "Multiplicative effect on F: 0.996^CONMED_H2RA (no clinically meaningful effect). Operational definition: CONMED_H2RA use covering >= 80% of the PK assessment phase (Goel 2016 Methods).",
      source_name        = "H2"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat meal at dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other meal condition; in Goel 2016 the reference is 2 h post light meal for cancer-patient doses or overnight fast for the HV-fasted arm doses)",
      notes              = "Per-dose-record indicator. Multiplicative effect on F: 5.74^FED_HIGHFAT (~5.7-fold higher bioavailability after a high-fat meal vs the 2-h-post-light-meal reference) and on Ka: 1.01^FED_HIGHFAT (no meaningful effect). The same dataset column FATM = 1 / Fatmeal = 1 carries the indicator for both Ka and F covariate effects in Goel 2016 (Table 2 theta18 on Ka, theta23 on F).",
      source_name        = "Fatmeal / FATM"
    ),
    FED = list(
      description        = "Fed-state-at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted at dosing). In Goel 2016 the healthy-volunteer fasted-arm effect on F is derived from DIS_HEALTHY = 1 AND FED = 0 (overnight fast in studies A1102 / A2114); patient records typically have FED = 1 (2 h post light meal); the high-fat-meal arm uses FED = 1 AND FED_HIGHFAT = 1.",
      notes              = "Per-dose-record indicator. The Goel 2016 e_healthy_fast_f effect on F is applied via (DIS_HEALTHY * (1 - FED)) rather than the retired composite HV_FAST column; HV_FAST was removed from the canonical register on 2026-05-11.",
      source_name        = "Fed / HV.Fasting"
    ),
    DOSE = list(
      description        = "Per-dose-event sonidegib dose level (mg)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on F with reference 100 mg: F multiplier = (DOSE/100)^(-0.342). Sonidegib bioavailability decreases super-linearly with dose because of solubility-limited absorption. Goel 2016 fits the dose effect across the 100-3000 mg single-dose and 200-1500 mg QD multi-dose ranges. Supply DOSE per dose record at the same value as the rxode2 amt column.",
      source_name        = "Dose"
    ),
    MULTI_DOSE_PT = list(
      description        = "Multiple-dose-phase-in-cancer-patients indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single-dose / run-in / healthy-volunteer dose)",
      notes              = "Per-dose-record indicator (Goel 2016 dataset column FMDD). Multiplicative effect on F: 1.16^MULTI_DOSE_PT (~16% higher apparent F during the multiple-dose phase relative to the first dose, attributed by the paper to occasional non-fasting compliance during long-term once-daily dosing). Set MULTI_DOSE_PT = 1 for every dose event from the start of the multiple-dose phase onward in cancer-patient studies (X2101 / X1101 / A2201); set MULTI_DOSE_PT = 0 for run-in single doses, healthy-volunteer doses, and the first patient dose.",
      source_name        = "FMDD"
    )
  )

  population <- list(
    n_subjects     = 436L,                                         # Goel 2016 Results / Table 1, pooled across the five studies
    n_studies      = 5L,                                           # A1102, A2114, X2101, X1101, A2201 (Goel 2016 Methods)
    age_range      = "20-93 years (median 58)",                    # Goel 2016 Table 1, pooled
    age_median     = "58 years",
    weight_range   = "42-181 kg (median 73)",                      # Goel 2016 Table 1, pooled
    weight_median  = "73 kg",
    sex_female_pct = 32.6,                                         # Goel 2016 Table 1: 142 female / 436 total
    race_ethnicity = "Approximately 94% White, 6% Japanese; race effects beyond the Japanese-vs-Western contrast were not tested (Goel 2016 Results).",
    disease_state  = "Pooled healthy subjects (n = 85; studies A1102 + A2114) and patients with advanced solid tumors or locally-advanced/metastatic basal cell carcinoma (n = 351; studies X2101 + X1101 + A2201).",
    dose_range     = "Single doses 100, 200, 250, 400, 600, 800, 1200, 1500, 3000 mg; QD multiple doses 200, 400, 600, 800, 1000, 1500, 3000 mg; BID multiple doses 250, 400, 450, 750 mg. Approved indication is 200 mg QD (Goel 2016 Methods + label).",
    regions        = "Multi-national; 94% Western enrollment, 6% Japanese (Goel 2016 Results).",
    n_observations = 6510L,                                        # 6510 of 6680 reported PK samples were above LLOQ and used (Goel 2016 Results)
    cohort_split   = "Healthy subjects 85 / cancer patients 351 (Goel 2016 Methods).",
    studies        = c("A1102", "A2114", "X2101", "X1101", "A2201"),
    notes          = "Pooled population PK analysis. Demographics summarized from Goel 2016 Table 1; the median (min, max) figures above are for the 'Pooled' column. Healthy subjects in A1102 (Japanese, fasted) and A2114 (non-Japanese, fasted or high-fat-meal arms) received single doses; cancer patients dosed daily 2 h after a light meal in X2101 / X1101 / A2201 with run-in single dose followed by daily continuous dosing."
  )

  ini({
    # Structural PK parameters at the reference covariate set:
    #   AGE = 58 yr, WT = 73 kg, SEXF = 0 (male), CRCL = 93.3 mL/min, ALB = 43 g/L,
    #   ALT = 0.42 (ULN-normalized), TBILI = 0.38 (ULN-normalized), RACE_JAPANESE = 0
    #   (Western), DIS_HEALTHY = 0 (cancer patient), CONMED_PPI = CONMED_H2RA = FED_HIGHFAT = 0, FED = 1,
    #   MULTI_DOSE_PT = 0 (single-dose / first-dose), DOSE = 100 mg.
    lka  <- log(0.219); label("Absorption rate at the reference covariate set (Ka, 1/h)")                     # Goel 2016 Table 2 full-model theta17
    lcl  <- log(10.2);  label("Apparent oral clearance at the reference covariate set (CL/F, L/h)")           # Goel 2016 Table 2 full-model theta1
    lvc  <- log(145);   label("Apparent central volume of distribution at the reference covariate set (Vc/F, L)") # Goel 2016 Table 2 full-model theta11
    lq   <- log(215);   label("Apparent inter-compartmental clearance at the reference covariate set (Q/F, L/h)") # Goel 2016 Table 2 full-model theta14
    lvp  <- log(8414);  label("Apparent peripheral volume of distribution at the reference covariate set (Vp/F, L)") # Goel 2016 Table 2 full-model theta15
    ltlag <- log(0.474); label("Absorption lag time (Tlag, h)")                                              # Goel 2016 Table 2 full-model theta19

    # Covariate effects on CL/F (Goel 2016 Table 2 full model thetas 2-10)
    e_age_cl       <- -0.375; label("Power exponent of AGE on CL/F (unitless; reference 58 y)")              # Goel 2016 Table 2 theta2
    e_wt_cl        <- -0.093; label("Power exponent of WT on CL/F (unitless; reference 73 kg)")              # Goel 2016 Table 2 theta3
    e_sexf_cl      <-  0.846; label("Multiplicative female-vs-male CL/F ratio (applied as ratio^SEXF)")      # Goel 2016 Table 2 theta4
    e_crcl_cl      <- -0.043; label("Power exponent of CRCL on CL/F (unitless; reference 93.3 mL/min)")      # Goel 2016 Table 2 theta5
    e_alb_cl       <- -1.44;  label("Power exponent of ALB on CL/F (unitless; reference 43 g/L)")            # Goel 2016 Table 2 theta6
    e_alt_cl       <-  0.120; label("Power exponent of ULN-normalized ALT on CL/F (unitless; reference 0.42)") # Goel 2016 Table 2 theta7
    e_tbili_cl     <-  0.064; label("Power exponent of ULN-normalized TBILI on CL/F (unitless; reference 0.38)") # Goel 2016 Table 2 theta8
    e_race_japanese_cl <- 0.905; label("Multiplicative Japanese-vs-Western CL/F ratio (applied as ratio^RACE_JAPANESE)") # Goel 2016 Table 2 theta9
    e_healthy_cl    <-  2.96;  label("Multiplicative HV-vs-patient CL/F ratio (applied as ratio^DIS_HEALTHY)")     # Goel 2016 Table 2 theta10

    # Covariate effects on Vc/F (Goel 2016 Table 2 full model thetas 12-13)
    e_wt_vc        <-  0.731; label("Power exponent of WT on Vc/F (unitless; reference 73 kg)")              # Goel 2016 Table 2 theta12
    e_alb_vc       <- -2.33;  label("Power exponent of ALB on Vc/F (unitless; reference 43 g/L)")            # Goel 2016 Table 2 theta13

    # Covariate effects on Vp/F (Goel 2016 Table 2 full model theta16; the WT effect on Vp/F
    # is fixed equal to the WT effect on Vc/F per the paper, so no separate parameter)
    e_alb_vp       <- -0.087; label("Power exponent of ALB on Vp/F (unitless; reference 43 g/L)")            # Goel 2016 Table 2 theta16

    # Covariate effects on Ka (Goel 2016 Table 2 full model theta18)
    e_fed_highfat_ka <-  1.01;  label("Multiplicative high-fat-meal Ka ratio (applied as ratio^FED_HIGHFAT)")    # Goel 2016 Table 2 theta18

    # Covariate effects on F (Goel 2016 Table 2 full model thetas 20-25; the structural F is
    # fixed at 1 at the reference covariate set, with covariate-driven multiplicative shifts)
    e_conmed_h2ra_f         <-  0.996; label("Multiplicative significant-CONMED_H2RA F ratio (applied as ratio^CONMED_H2RA)")            # Goel 2016 Table 2 theta20
    e_conmed_ppi_f          <-  0.696; label("Multiplicative significant-CONMED_PPI F ratio (applied as ratio^CONMED_PPI)")              # Goel 2016 Table 2 theta21
    e_healthy_fast_f      <-  0.855; label("Multiplicative healthy-fasted F ratio (applied as ratio^(DIS_HEALTHY * (1 - FED)))")  # Goel 2016 Table 2 theta22
    e_fed_highfat_f    <-  5.74;  label("Multiplicative high-fat-meal F ratio (applied as ratio^FED_HIGHFAT)")          # Goel 2016 Table 2 theta23
    e_dose_f         <- -0.342; label("Power exponent of DOSE on F (unitless; reference 100 mg)")                   # Goel 2016 Table 2 theta24
    e_multi_dose_pt_f <- 1.16;  label("Multiplicative multi-dose-phase F ratio (applied as ratio^MULTI_DOSE_PT)")   # Goel 2016 Table 2 theta25

    # Inter-individual variability. Goel 2016 Table 2 reports IIV as %CV; the stored
    # log-normal variance follows omega^2 = log(CV^2 + 1).
    # CL/F 64.8% CV -> 0.3505; Vc/F 213% CV -> 1.7115; Q/F 106% CV -> 0.7531;
    # Vp/F 78.9% CV -> 0.4839; Ka 44.2% CV -> 0.1786.
    # Q/F-Vp/F correlation 0.692 -> covariance 0.692 * sqrt(0.7531 * 0.4839) = 0.4178.
    etalka ~ 0.1786                                              # Goel 2016 Table 2 IIV Ka 44.2% CV
    etalcl ~ 0.3505                                              # Goel 2016 Table 2 IIV CL/F 64.8% CV
    etalvc ~ 1.7115                                              # Goel 2016 Table 2 IIV Vc/F 213% CV
    etalq + etalvp ~ c(0.7531, 0.4178, 0.4839)                   # Goel 2016 Table 2 IIV Q/F 106% CV, Vp/F 78.9% CV, CorrVp/F-Q/F = 0.692

    # Residual error. Goel 2016 Table 2 reports the additive component as 1.1e-4 ng/mL with
    # large RSE and footnotes it as 'negligible'; only the proportional component (31.8% CV
    # on the linear concentration scale) is retained here.
    propSd <- 0.318;  label("Proportional residual error (fraction)")                                        # Goel 2016 Table 2 sigma_mult 31.8% CV
  })

  model({
    # Individual PK parameters with covariate effects (Goel 2016 Eq. 1 power-form covariate
    # model: continuous covariates as (cov / cov_ref)^theta; dichotomous covariates as
    # theta^cov_di). The reference covariate values appear inside this block as numeric
    # constants matching the ini() reference set.
    ka <- exp(lka + etalka) *
      e_fed_highfat_ka^FED_HIGHFAT
    cl <- exp(lcl + etalcl) *
      (AGE / 58)^e_age_cl *
      (WT / 73)^e_wt_cl *
      e_sexf_cl^SEXF *
      (CRCL / 93.3)^e_crcl_cl *
      (ALB / 43)^e_alb_cl *
      (ALT / 0.42)^e_alt_cl *
      (TBILI / 0.38)^e_tbili_cl *
      e_race_japanese_cl^RACE_JAPANESE *
      e_healthy_cl^DIS_HEALTHY
    vc <- exp(lvc + etalvc) *
      (WT / 73)^e_wt_vc *
      (ALB / 43)^e_alb_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) *
      (WT / 73)^e_wt_vc *
      (ALB / 43)^e_alb_vp
    tlag_depot <- exp(ltlag)
    fdepot <- e_conmed_h2ra_f^CONMED_H2RA *
      e_conmed_ppi_f^CONMED_PPI *
      e_healthy_fast_f^(DIS_HEALTHY * (1 - FED)) *
      e_fed_highfat_f^FED_HIGHFAT *
      (DOSE / 100)^e_dose_f *
      e_multi_dose_pt_f^MULTI_DOSE_PT

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: depot -> central -> peripheral1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability and absorption-lag time apply to the depot compartment.
    f(depot)    <- fdepot
    alag(depot) <- tlag_depot

    # Concentration: dose in mg, volume in L -> central / vc in mg/L.
    # Multiply by 1000 to express Cc in ng/mL to match Goel 2016 Tables 3-4.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
