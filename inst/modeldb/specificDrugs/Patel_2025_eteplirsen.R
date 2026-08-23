Patel_2025_eteplirsen <- function() {
  description <- "Three-compartment IV-infusion population PK model for eteplirsen, an exon-51-skipping phosphorodiamidate morpholino oligomer (PMO), in 157 male patients with Duchenne muscular dystrophy aged 6 months to 16.4 years pooled from six clinical studies (Patel 2025). Clearance uses an age-cutoff structure at 4 years: a separate typical clearance is estimated in each age stratum (6.98 L/h for age > 4 years, 4.97 L/h for age <= 4 years), both normalized to a 37 kg reference weight and an eGFR of 145 mL/min/1.73 m^2. Body weight enters every disposition parameter allometrically with exponents fixed at 0.75 on the three clearance terms and 1 on the three volumes; cystatin-C-based CKD-EPI eGFR enters CL as a power effect with an estimated exponent of 1.60. Interindividual variability is a full 4x4 block on CL, V1, V2 and Q2 (none on V3 or Q3), and residual error is additive on the log scale. Both age strata come from a single joint NONMEM fit."
  reference <- "Patel Y, Orogun L, Yocum N, Rodino-Klapac LR, East L. A population pharmacokinetic model to inform extension of the eteplirsen dosing regimen across the broad DMD population. CPT Pharmacometrics Syst Pharmacol. 2025;14(5):891-901. doi:10.1002/psp4.70001"
  vignette <- "Patel_2025_eteplirsen"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "eteplirsen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "eteplirsen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "eteplirsen", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Table 1 (All studies): median 27.1 kg, range 6.80-68.9. The reference weight of 37 kg is NOT the cohort median; Results section 3.4 states 'A body weight of 37 kg was selected as a reference male DMD patient (> 4 years of age) with an eGFR of 145 mL/min/1.73 m^2 to normalize the younger pediatric population (<= 4 years old) to an older reference population where the majority of PK data were collected from.' Enters all six disposition parameters as (WT/37)^k with k fixed at 0.75 for CL, Q2 and Q3 and at 1 for V1, V2 and V3 (Results section 3.3: 'The body weight effects on volume and clearance parameters using fixed allometric exponents of 0.75 for clearance parameters and 1 for volume parameters described data well'; Eq. 1-6).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Table 1 (All studies): median 8.37 years, range 0.55-16.4. Used ONLY as the 4-year cutoff selecting between the two typical clearances, exactly as written in Eq. 1: 'ifelse(Age_i > 4 yrs, theta1, theta2)'. It does not enter as a continuous term: Results section 3.4 states that because body weight and age were correlated (r = 0.781, Figure S1), 'inclusion of a continuous covariate effect of age in addition to body weight resulted in poor characterization of covariate effects', and that a <= 3 vs > 3 year cutoff was rejected for having too few subjects (about 5.7% of the population) while the <= 4 vs > 4 year cutoff 'described the data well and was selected in the final model'.",
      source_name        = "Age"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from serum cystatin C by the CKD-EPI equation (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3: 'Because creatinine is a byproduct of the muscle breakdown occurring in patients with DMD, creatinine clearance was not an appropriate marker of renal function for this analysis. Instead, continuous serum cystatin C was used to estimate the eGFR using the CKD-EPI equation.' Table 1 (All studies): median 145 mL/min/1.73 m^2, range 85.5-180. The reference value 145 in Eq. 1 is that cohort median. eGFR was not measured in Study 4658-28 and was imputed there using the median eGFR observed in patients of the same age (Methods section 2.3). Enters CL as (CRCL/145)^1.60. Stored under canonical CRCL, which admits BSA-normalized eGFR variants; the estimating equation is documented here per the register's per-model requirement. Note this is the cystatin-C-based CKD-EPI variant, not the creatinine-based one, and it is deliberately supranormal for the DMD population (no subject had eGFR below 85.5 mL/min/1.73 m^2).",
      source_name        = "eGFR"
    )
  )

  # Covariates the paper screened (Methods section 2.3) but did not retain in
  # the final model. Documentation only: Patel 2025 reports no point estimate
  # for any of them, so none can be encoded. The paper also screened baseline
  # ambulatory status, race, ethnicity, DMD genotype (exon-deletion pattern),
  # single-vs-multiple steroid use and steroid strength, and kidney injury
  # molecule 1; those have no clean canonical column here and are recorded in
  # the vignette instead.
  covariatesDataExcluded <- list(
    LBM = list(
      description        = "Lean body mass. Screened as an alternative body-size descriptor and rejected.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3 lists lean body mass among the body-size covariates evaluated. Results section 3.3: 'Alternative models to characterize body size effect on the model parameters resulted in non-identifiability, unreasonable estimates in some of the parameters, higher objective function values, higher model condition number, increased residual error, and/or poor model diagnostics.' Total body weight was retained instead. No point estimate published."
    ),
    BMI = list(
      description        = "Body mass index. Screened as an alternative body-size descriptor and rejected.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3 and Results section 3.3, same rejection rationale as LBM. No point estimate published."
    ),
    BSA = list(
      description        = "Body surface area. Screened as an alternative body-size descriptor and rejected.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3 and Results section 3.3, same rejection rationale as LBM. Discussion notes a further reason it is hard to separate: 'eGFR is normalized to body surface area, and given that body weight was included in the CL analysis using allometry, it is possible that eGFR does not have a significant effect considering the co-linearity between body surface area and body weight.' No point estimate published."
    ),
    CYSC = list(
      description        = "Serum cystatin C. Measured, but entered the model only through the derived CKD-EPI eGFR rather than as a covariate in its own right.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3 lists 'renal biomarkers if available (e.g., cystatin C or kidney injury molecule 1)' among the screened covariates. The final model uses the CKD-EPI eGFR derived from cystatin C (stored as CRCL), not the raw concentration. No point estimate published for a direct CYSC effect."
    ),
    ALB = list(
      description        = "Serum albumin. Screened as a liver-function covariate and rejected.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3 lists albumin among the liver-function covariates evaluated. Results section 3.4: 'No specific trends observed in the plots of interindividual random effect versus demographic, body size, and laboratory parameter were observed in the final model, which suggested that no additional covariate effects remained to account for in the model (Figures S2 and S3).' No point estimate published."
    ),
    ALT = list(
      description        = "Alanine aminotransferase. Screened as a liver-function covariate and rejected.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3; rejected per Results section 3.4 (Figures S2 and S3). Note that in DMD, ALT and AST are largely of skeletal-muscle rather than hepatic origin. No point estimate published."
    ),
    AST = list(
      description        = "Aspartate aminotransferase. Screened as a liver-function covariate and rejected.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3; rejected per Results section 3.4 (Figures S2 and S3). As with ALT, elevated values in DMD reflect muscle breakdown. No point estimate published."
    ),
    TBILI = list(
      description        = "Total bilirubin. Screened as a liver-function covariate and rejected.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patel 2025 Methods section 2.3 lists bilirubin among the liver-function covariates evaluated; rejected per Results section 3.4 (Figures S2 and S3). The paper does not state whether total or direct bilirubin was used, nor the reporting units. No point estimate published."
    ),
    CONMED_STEROID = list(
      description        = "Concomitant corticosteroid use. Screened as a concomitant-medication covariate and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant corticosteroid)",
      notes              = "Patel 2025 Methods section 2.3 screened 'concomitant medications (steroid strength, single or multiple steroids)'. Corticosteroids are standard of care in DMD, so the untreated group is small. Rejected per Results section 3.4 (Figure S3, categorical-covariate eta plots). No point estimate published."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 157,
    n_studies      = 6,
    n_observations = 3258,
    age_range      = "0.55-16.4 years (6 months to 16.4 years)",
    age_median     = "8.37 years",
    weight_range   = "6.80-68.9 kg",
    weight_median  = "27.1 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 84, Black = 2, Asian = 6, `Pacific Islander` = 1, Other = 3, Missing = 4),
    disease_state  = "Duchenne muscular dystrophy with a confirmed mutation amenable to exon 51 skipping",
    renal_function = "eGFR (cystatin-C CKD-EPI) median 145 mL/min/1.73 m^2, range 85.5-180; no subject below 85.5, so the model is not informed about renal impairment",
    dose_range     = "0.5-50 mg/kg weekly as an approximately 1-h IV infusion; 30 mg/kg/week is the approved regimen",
    regions        = "Not reported by region; the six studies are NCT03218995 (4658-102), NCT01396239 (4658-201), NCT01540409 (4658-202), NCT02420379 (4658-203), NCT00844597 (4658-28) and NCT02255552 (4658-301)",
    notes          = "Patel 2025 Table 1 (study design, N, age, weight and eGFR by study) and Table S1 (ethnicity, race and exon-deletion mutations by study). All patients are male because DMD is X-linked recessive. 37 samples (1.1%) below the limit of quantification were excluded. Ethnicity: 8% Hispanic or Latino, 87% not Hispanic or Latino, 5% unknown."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters. Patel 2025 Table 2 ("Estimate" column) with
    # the covariate structure of Eq. 1-6 (p. 896):
    #
    #   CL_i  = exp(ifelse(Age_i > 4 yrs, theta1, theta2) + eta1_i)
    #             * (WT_i/37)^0.75 * (eGFR_i/145)^theta8          (1)
    #   V1_i  = exp(theta3 + eta2_i) * (WT_i/37)^1                (2)
    #   V2_i  = exp(theta4 + eta3_i) * (WT_i/37)^1                (3)
    #   Q2_i  = exp(theta5 + eta4_i) * (WT_i/37)^0.75             (4)
    #   V3_i  = exp(theta6)          * (WT_i/37)^1                (5)
    #   Q3_i  = exp(theta7)          * (WT_i/37)^0.75             (6)
    #
    # Methods section 2.5.1 states the parameters were "modeled in the
    # log-domain" with P_i = exp(P_hat + eta_Pi), so the thetas inside
    # exp() are log-scale. Table 2 reports the BACK-TRANSFORMED values
    # (its Units column reads L/h and L), which is why each entry below
    # is written log(<Table 2 estimate>).
    #
    # NONMEM V1/V2/V3 and Q2/Q3 map onto the nlmixr2lib canonical names
    # as V1 -> vc, V2 -> vp, V3 -> vp2, Q2 -> q, Q3 -> q2.
    #
    # This is ONE jointly fit model, not two: only the typical clearance
    # is age-stratum-specific, so it alone carries the symmetric
    # `_agele4` / `_agegt4` suffix pair. Every volume, every
    # inter-compartmental clearance, the eGFR exponent, the whole IIV
    # block and the residual error are shared across strata.
    # -----------------------------------------------------------------

    lcl_agegt4 <- log(6.98)  ; label("Typical clearance for age > 4 years at WT = 37 kg, CRCL = 145 mL/min/1.73 m^2 (L/h)")   # Patel 2025 Table 2: theta1 = 6.98 L/h (95% CI 6.68-7.31); Eq. 1
    lcl_agele4 <- log(4.97)  ; label("Typical clearance for age <= 4 years at WT = 37 kg, CRCL = 145 mL/min/1.73 m^2 (L/h)")  # Patel 2025 Table 2: theta2 = 4.97 L/h (95% CI 4.48-5.47); Eq. 1
    lvc        <- log(2.21)  ; label("Central volume of distribution at WT = 37 kg (V1, L)")                                  # Patel 2025 Table 2: theta3 = 2.21 L (95% CI 1.93-2.55); Eq. 2
    lvp        <- log(2.23)  ; label("First peripheral volume of distribution at WT = 37 kg (V2, L)")                         # Patel 2025 Table 2: theta4 = 2.23 L (95% CI 1.90-2.64); Eq. 3
    lq         <- log(0.334) ; label("Inter-compartmental clearance to peripheral1 at WT = 37 kg (Q2, L/h)")                  # Patel 2025 Table 2: theta5 = 0.334 L/h (95% CI 0.262-0.422); Eq. 4
    lvp2       <- log(4.12)  ; label("Second peripheral volume of distribution at WT = 37 kg (V3, L)")                        # Patel 2025 Table 2: theta6 = 4.12 L (95% CI 3.81-4.43); Eq. 5
    lq2        <- log(3.49)  ; label("Inter-compartmental clearance to peripheral2 at WT = 37 kg (Q3, L/h)")                  # Patel 2025 Table 2: theta7 = 3.49 L/h (95% CI 3.29-3.70); Eq. 6

    # -----------------------------------------------------------------
    # Covariate effects.
    #
    # The eGFR exponent is the only ESTIMATED covariate coefficient
    # (Table 2 gives it a 95% CI). The six body-weight exponents were
    # held constant at the canonical allometric values and carry no
    # uncertainty in Table 2, so they are fixed().
    # -----------------------------------------------------------------

    e_crcl_cl <- 1.60      ; label("Power exponent on (CRCL/145) for CL (unitless)")  # Patel 2025 Table 2: theta8 (eGFR on CL) = 1.60 (95% CI 1.33-1.92); Eq. 1

    e_wt_cl   <- fixed(0.75) ; label("Power exponent on (WT/37) for CL (unitless)")   # Patel 2025 Eq. 1 exponent 0.75; Results section 3.3 "fixed allometric exponents of 0.75 for clearance parameters"
    e_wt_q    <- fixed(0.75) ; label("Power exponent on (WT/37) for Q2 (unitless)")   # Patel 2025 Eq. 4 exponent 0.75; Results section 3.3
    e_wt_q2   <- fixed(0.75) ; label("Power exponent on (WT/37) for Q3 (unitless)")   # Patel 2025 Eq. 6 exponent 0.75; Results section 3.3
    e_wt_vc   <- fixed(1.0)  ; label("Power exponent on (WT/37) for V1 (unitless)")   # Patel 2025 Eq. 2 exponent 1; Results section 3.3 "and 1 for volume parameters"
    e_wt_vp   <- fixed(1.0)  ; label("Power exponent on (WT/37) for V2 (unitless)")   # Patel 2025 Eq. 3 exponent 1; Results section 3.3
    e_wt_vp2  <- fixed(1.0)  ; label("Power exponent on (WT/37) for V3 (unitless)")   # Patel 2025 Eq. 5 exponent 1; Results section 3.3

    # -----------------------------------------------------------------
    # Interindividual variability: a full 4x4 block on CL, V1, V2 and Q2.
    # Results section 3.3: "IIV on all PK parameters except on V3 and Q3"
    # and "Interindividual random effect distributions were modeled using
    # an additive model on log domain on CL, V1, V2, and Q2, with
    # covariance terms between all effects."
    #
    # Table 2 reports the four DIAGONAL terms as CV% and the six
    # off-diagonal terms as covariances on the log (variance) scale.
    # Methods section 2.5.1 prints the conversion the authors used,
    #     %CV = 100 * sqrt(exp(omega^2) - 1),
    # so each diagonal is recovered as omega^2 = log(1 + CV^2):
    #
    #   CL : log(1 + 0.235^2) = 0.05375412
    #   V1 : log(1 + 0.498^2) = 0.22154497
    #   V2 : log(1 + 0.705^2) = 0.40347978
    #   Q2 : log(1 + 1.030^2) = 0.72314292
    #
    # The resulting 4x4 matrix is positive definite (smallest eigenvalue
    # 0.0119) with correlations 0.50-0.87, so it simulates directly.
    #
    # nlmixr2 block order is lower-triangle, row-wise:
    #   var(cl),
    #   cov(cl,vc), var(vc),
    #   cov(cl,vp), cov(vc,vp), var(vp),
    #   cov(cl,q),  cov(vc,q),  cov(vp,q), var(q)
    #
    # Source of each of the ten entries, in that order (all Patel 2025
    # Table 2; 95% CI in parentheses). Comments cannot be interleaved
    # inside the c(...) call itself -- rxode2's ini() reparser turns a
    # trailing comment into a label and fails on a multi-line block --
    # so the trace lives here:
    #
    #   var(cl)      0.05375412  omega_CL = 23.5 CV% (20.2-27.4), eta-shrinkage 11.0%
    #   cov(cl,vc)   0.0541      Cov(CL,V1) = 0.0541 (0.0264-0.0924)
    #   var(vc)      0.22154497  omega_V1 = 49.8 CV% (40.1-61.7), eta-shrinkage 22.1%
    #   cov(cl,vp)   0.110       Cov(CL,V2) = 0.110 (0.0719-0.159)
    #   cov(vc,vp)   0.193       Cov(V1,V2) = 0.193 (0.103-0.315)
    #   var(vp)      0.40347978  omega_V2 = 70.5 CV% (55.6-91.5), eta-shrinkage 24.9%
    #   cov(cl,q)    0.167       Cov(CL,Q2) = 0.167 (0.116-0.234)
    #   cov(vc,q)    0.295       Cov(V1,Q2) = 0.295 (0.178-0.447)
    #   cov(vp,q)    0.470       Cov(V2,Q2) = 0.470 (0.310-0.695)
    #   var(q)       0.72314292  omega_Q2 = 103 CV% (79.8-137), eta-shrinkage 20.4%
    # -----------------------------------------------------------------

    etalcl + etalvc + etalvp + etalq ~ c(
      0.05375412,
      0.0541, 0.22154497,
      0.110, 0.193, 0.40347978,
      0.167, 0.295, 0.470, 0.72314292
    )

    # -----------------------------------------------------------------
    # Residual error. Methods section 2.5.1: "Plasma concentrations were
    # log-transformed for model-fitting purposes" and "The residual error
    # model was described by an additive error model on the log domain",
    # i.e. log(C_ij) = log(C_hat_ij) + eps_ij. That is exactly nlmixr2's
    # lnorm() error structure, whose parameter is the log-scale SD.
    #
    # Table 2 quotes sigma_RV as 35.1 CV%, presented under the same
    # sqrt(exp(sigma^2) - 1) transform as the IIVs, so
    #   sigma = sqrt(log(1 + 0.351^2)) = 0.3408558.
    # -----------------------------------------------------------------

    expSd <- 0.3408558 ; label("Log-normal (log-scale additive) residual SD")  # Patel 2025 Table 2: sigma_RV = 35.1 CV% (95% CI 34.1-36.1) back-transformed via sqrt(log(1 + CV^2)); epsilon-shrinkage 67.7%
  })

  model({
    # ---- 1. Age-stratum indicator (Eq. 1, cutoff at 4 years) ----
    is_agele4 <- (AGE <= 4.0)

    # ---- 2. Individual PK parameters ----
    # Only the typical clearance switches on the age stratum; the weight
    # allometry, the eGFR power effect and the IIV are shared.
    cl_typ <- exp(lcl_agele4) * is_agele4 +
      exp(lcl_agegt4) * (1.0 - is_agele4)
    cl  <- cl_typ * exp(etalcl) * (WT / 37.0)^e_wt_cl * (CRCL / 145.0)^e_crcl_cl

    vc  <- exp(lvc  + etalvc) * (WT / 37.0)^e_wt_vc
    vp  <- exp(lvp  + etalvp) * (WT / 37.0)^e_wt_vp
    q   <- exp(lq   + etalq)  * (WT / 37.0)^e_wt_q
    vp2 <- exp(lvp2)          * (WT / 37.0)^e_wt_vp2
    q2  <- exp(lq2)           * (WT / 37.0)^e_wt_q2

    # ---- 3. Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 4. ODE system. Three compartments with linear elimination from
    #         the central compartment. Eteplirsen is given as an
    #         approximately 1-h IV infusion (Methods section 2.2), which
    #         is the "zero-order absorption" of Results section 3.3: it
    #         is a property of the dosing record (rate / dur), not an
    #         absorption compartment, so doses go straight into central.
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ---- 5. Observation and residual error ----
    # Dose in mg, vc in L -> central/vc has units mg/L, i.e. ug/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
