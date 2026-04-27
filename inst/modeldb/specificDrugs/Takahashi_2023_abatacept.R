Takahashi_2023_abatacept <- function() {
  description <- "Two-compartment IV population PK model for abatacept (CTLA4-Ig Fc-fusion) pooled across 685 adult/pediatric patients with rheumatoid arthritis or polyarticular juvenile idiopathic arthritis and adult/pediatric patients receiving allogeneic hematopoietic cell transplantation in the ABA2 trial (Takahashi 2023). Linear elimination, allometric weight scaling on CL/Vc/Vp/Q with estimated exponents, and a three-level cohort categorical (RA/JIA reference, ABA2 HLA 7/8, ABA2 HLA 8/8) on CL and Vc."
  reference <- "Takahashi T, Al-Kofahi M, Jaber M, Bratrude B, Betz K, Suessmuth Y, Yu A, Neuberg DS, Choi SW, Davis J, Duncan C, Giller R, Grimley M, Harris AC, Jacobsohn D, Lalefar N, Farhadfar N, Pulsipher MA, Shenoy S, Petrovic A, Schultz KR, Yanik GA, Blazar BR, Horan JT, Watkins B, Langston A, Qayed M, Kean LS. Higher abatacept exposure after transplant decreases acute GVHD risk without increasing adverse events. Blood. 2023 Aug 24;142(8):700-710. doi:10.1182/blood.2023020035"
  vignette <- "Takahashi_2023_abatacept"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg (Takahashi 2023 Supplemental Table 4: CL70kg, Vc70kg, Vp70kg, Q70kg are typical values for a 70 kg subject). Power scaling on CL (estimated exponent 0.65), Vc (0.70), Vp (1.02), and Q (0.63) per Supplemental Table 4 footnote equations.",
      source_name        = "WT"
    ),
    STUDY_ABA2_HLA78 = list(
      description        = "ABA2 trial HLA 7/8 (one-allele-mismatched donor) hematopoietic-cell-transplant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is in the RA/JIA pooled reference cohort or in the ABA2 HLA 8/8 cohort; the latter is captured separately by STUDY_ABA2_HLA88)",
      notes              = "Takahashi 2023 Supplemental Table 4: multiplicative Ratio = 0.70 on CL (95% CI 0.65 - 0.77) and 0.99 on Vc (95% CI 0.91 - 1.08) vs the RA/JIA reference (Ratio 1 fixed). RA/JIA reference is reproduced by STUDY_ABA2_HLA78 = STUDY_ABA2_HLA88 = 0.",
      source_name        = "Cohort (subset == 'ABA2 7/8')"
    ),
    STUDY_ABA2_HLA88 = list(
      description        = "ABA2 trial HLA 8/8 (allele-matched donor) hematopoietic-cell-transplant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is in the RA/JIA pooled reference cohort or in the ABA2 HLA 7/8 cohort; the latter is captured separately by STUDY_ABA2_HLA78)",
      notes              = "Takahashi 2023 Supplemental Table 4: multiplicative Ratio = 0.91 on CL (95% CI 0.86 - 0.97) and 1.32 on Vc (95% CI 1.23 - 1.42) vs the RA/JIA reference (Ratio 1 fixed). RA/JIA reference is reproduced by STUDY_ABA2_HLA78 = STUDY_ABA2_HLA88 = 0.",
      source_name        = "Cohort (subset == 'ABA2 8/8')"
    )
  )

  population <- list(
    n_subjects        = 685L,
    n_observations    = 4872L,
    n_studies         = 8L,
    age_range         = "6 - 84 years (median 45)",
    weight_range      = "14.4 - 186.8 kg (median 67.9)",
    sex_female_pct    = 67,
    race_ethnicity    = c(White = 85),
    disease_state     = "Pooled adult rheumatoid arthritis and pediatric polyarticular juvenile idiopathic arthritis (n = 570 across studies IM103002, IM101100, IM101101, IM101102, IM101029, IM101031, IM101033) plus adult and pediatric (>= 6 years) hematopoietic-cell-transplant recipients in the ABA2 trial (IM101311 / NCT01743131): HLA 7/8 cohort n = 42, HLA 8/8 cohort n = 73.",
    dose_range        = "IV abatacept 0.5 - 10 mg/kg infusions over 0.5 - 1 h across the seven RA/JIA studies (multiple regimens per study; weight-tier dosing approximating 10 mg/kg with 500 / 750 / 1000 mg vial-rounding for adult Q4W maintenance). ABA2 trial: weight-tiered IV abatacept 10 mg/kg (maximum 1000 mg) over 1 h on days -1, +5, +14, and +28 relative to graft.",
    regions           = "Multi-regional pooled analysis (RA/JIA studies global; ABA2 multi-center US/Canada).",
    reference_subject = "70 kg subject in the RA/JIA reference cohort (STUDY_ABA2_HLA78 = STUDY_ABA2_HLA88 = 0). Takahashi 2023 Supplemental Table 4 footnote: CL70kg / Vc70kg / Vp70kg / Q70kg are 'for a 70 kg subject with RA' (RA/JIA cohort serves as the reference category).",
    notes             = "Pooled population PK analysis of 4872 abatacept serum concentrations from 685 patients across 8 studies (Takahashi 2023 Supplemental Tables 1 and 2). 4.9% of measurements (237 / 4872) were below the assay's accuracy/precision thresholds and retained at the reported value (Byon 2008 approach). Estimation by FOCE-I in NONMEM 7.5; parameter uncertainty quantified by sampling-importance-resampling. Forward p < 0.01, backward p < 0.001 stepwise covariate modeling identified weight (estimated exponents) and a three-level cohort categorical (RA/JIA, ABA2 7/8, ABA2 8/8) as the retained PK covariates after dropping clinically-not-meaningful sex, eGFR, and albumin effects (Supplemental Table 3 'Final' row)."
  )

  ini({
    # Structural PK parameters - Takahashi 2023 Supplemental Table 4 'Final
    # Model Results / Estimate' column. CL and Q reported in L/h; converted to
    # L/day (x 24) because this model keeps time in days. Vc and Vp unchanged.
    lcl     <- log(0.023 * 24); label("Typical clearance CL70kg for a 70 kg RA/JIA reference subject (L/day)")            # Takahashi 2023 Suppl. Tab. 4: CL70 kg = 0.023 L/h
    lvc     <- log(3.06);       label("Typical central volume Vc70kg for a 70 kg RA/JIA reference subject (L)")           # Takahashi 2023 Suppl. Tab. 4: Vc70 kg = 3.06 L
    lvp     <- log(5.34);       label("Typical peripheral volume Vp70kg for a 70 kg subject (L)")                          # Takahashi 2023 Suppl. Tab. 4: Vp70 kg = 5.34 L
    lq      <- log(0.031 * 24); label("Typical inter-compartmental clearance Q70kg for a 70 kg subject (L/day)")           # Takahashi 2023 Suppl. Tab. 4: Q70 kg = 0.031 L/h

    # Allometric weight exponents (estimated, not fixed; Supplemental Table 4
    # 'Exponent of weight effect on ...' rows).
    e_wt_cl <- 0.65; label("Power exponent of (WT/70 kg) on CL (unitless)")                                                # Takahashi 2023 Suppl. Tab. 4: theta_WT_CL = 0.65
    e_wt_vc <- 0.70; label("Power exponent of (WT/70 kg) on Vc (unitless)")                                                # Takahashi 2023 Suppl. Tab. 4: theta_WT_Vc = 0.70
    e_wt_vp <- 1.02; label("Power exponent of (WT/70 kg) on Vp (unitless)")                                                # Takahashi 2023 Suppl. Tab. 4: theta_WT_Vp = 1.02
    e_wt_q  <- 0.63; label("Power exponent of (WT/70 kg) on Q (unitless)")                                                 # Takahashi 2023 Suppl. Tab. 4: theta_WT_Q  = 0.63

    # Cohort multiplicative ratios (RA/JIA reference; Ratio 1 fixed for both
    # parameters). Stored as ratio values (matching Supplemental Table 4
    # 'Ratio' column) and applied multiplicatively in model() so the parameter
    # values match the paper directly. STUDY_ABA2_HLA78 and STUDY_ABA2_HLA88
    # are mutually exclusive (at most one is 1 per subject); both 0 reproduces
    # the RA/JIA reference factor 1.
    r_hct78_cl <- 0.70; label("Multiplicative ratio on CL for ABA2 HLA 7/8 cohort (unitless)")                             # Takahashi 2023 Suppl. Tab. 4: theta_Cohort_CL (HLA 7/8) = 0.70
    r_hct88_cl <- 0.91; label("Multiplicative ratio on CL for ABA2 HLA 8/8 cohort (unitless)")                             # Takahashi 2023 Suppl. Tab. 4: theta_Cohort_CL (HLA 8/8) = 0.91
    r_hct78_vc <- 0.99; label("Multiplicative ratio on Vc for ABA2 HLA 7/8 cohort (unitless)")                             # Takahashi 2023 Suppl. Tab. 4: theta_Cohort_Vc (HLA 7/8) = 0.99
    r_hct88_vc <- 1.32; label("Multiplicative ratio on Vc for ABA2 HLA 8/8 cohort (unitless)")                             # Takahashi 2023 Suppl. Tab. 4: theta_Cohort_Vc (HLA 8/8) = 1.32

    # Inter-individual variability - Takahashi 2023 Supplemental Table 4 reports
    # IIV as CV(%); converted to log-normal omega^2 = log(CV^2 + 1). No IIV on Q
    # in the final model. IIVs are independent (no block correlation reported in
    # Supplemental Table 4).
    etalcl ~ 0.067858; label("IIV variance on log-CL (Takahashi 2023 Suppl. Tab. 4: 26.5%% CV; log(0.265^2 + 1))")
    etalvc ~ 0.036946; label("IIV variance on log-Vc (Takahashi 2023 Suppl. Tab. 4: 19.4%% CV; log(0.194^2 + 1))")
    etalvp ~ 0.167530; label("IIV variance on log-Vp (Takahashi 2023 Suppl. Tab. 4: 42.7%% CV; log(0.427^2 + 1))")

    # Residual error - Takahashi 2023 Supplemental Table 4 'Residual variability'
    # block: combined additive (0.049 ug/mL) + proportional (25.0% CV).
    addSd  <- 0.049; label("Additive residual error (ug/mL)")                                                              # Takahashi 2023 Suppl. Tab. 4: Additional = 0.049 ug/mL
    propSd <- 0.250; label("Proportional residual error (fraction)")                                                       # Takahashi 2023 Suppl. Tab. 4: Proportional CV = 25.0%
  })

  model({
    # Cohort multiplicative factors. STUDY_ABA2_HLA78 and STUDY_ABA2_HLA88 are
    # mutually exclusive binary indicators; both 0 reproduces the RA/JIA reference
    # (factor 1.0). The exponentiation form r^indicator yields r when the indicator
    # is 1 and 1 when the indicator is 0, exactly matching Supplemental Table 4's
    # multiplicative Ratio column.
    cohort_cl <- (r_hct78_cl ^ STUDY_ABA2_HLA78) * (r_hct88_cl ^ STUDY_ABA2_HLA88)
    cohort_vc <- (r_hct78_vc ^ STUDY_ABA2_HLA78) * (r_hct88_vc ^ STUDY_ABA2_HLA88)

    # Individual PK parameters (Takahashi 2023 Supplemental Table 4 footnote
    # equations: CLpop = CL70kg * theta_Cohort_CL * (WT/70)^theta_WT_CL, etc.).
    cl <- exp(lcl + etalcl) * cohort_cl * (WT / 70) ^ e_wt_cl
    vc <- exp(lvc + etalvc) * cohort_vc * (WT / 70) ^ e_wt_vc
    vp <- exp(lvp + etalvp)             * (WT / 70) ^ e_wt_vp
    q  <- exp(lq)                       * (WT / 70) ^ e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc gives mg/L = ug/mL,
    # matching Takahashi 2023 reporting units (39 ug/mL Ctrough_1 threshold).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
