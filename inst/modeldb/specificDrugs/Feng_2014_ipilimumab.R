Feng_2014_ipilimumab <- function() {
  description <- "Two-compartment population PK model for intravenous ipilimumab (anti-CTLA-4 IgG1) in patients with unresectable stage III or IV melanoma (Feng 2014)"
  reference <- "Feng Y, Masson E, Dai D, Parker SM, Berman D, Roy A. Model-based clinical pharmacology profiling of ipilimumab in patients with advanced melanoma. Br J Clin Pharmacol. 2014;78(1):106-117. doi:10.1111/bcp.12323"
  vignette <- "Feng_2014_ipilimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL (exponent CL_BW) and on Vc (exponent V_cBW) with reference weight 80 kg (Feng 2014 Table 2 footnote; reference subject defined in Figure 1 caption).",
      source_name        = "BW"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Effect on CL via the literal Feng 2014 form (log(LDH)/log(206))^CL_LDH -- a power of a ratio of logs, not the conventional (LDH/ref)^exponent. Feng 2014 Results state 'The value for LDH was log-transformed due to its right-skewed distribution'; reference LDH = 206 IU/L (Table 2 footnote). The same unusual literal form was carried forward by the same Bristol-Myers Squibb modeling group in the later ipilimumab popPK (Sanghavi 2020), giving a small (<20%) effect across the observed LDH range -- the conventional power-of-ratio form would inflate the effect well above the paper's narrative.",
      source_name        = "LDH"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 499L,
    n_studies      = 4L,
    age_range      = "26-86 years",
    age_median     = "57.75 years (SD 12.91) in the index dataset",
    weight_range   = "mean 80.11 kg (SD 16.87) in the index dataset",
    weight_median  = "not tabulated; mean 80.11 kg (SD 16.87) in the index dataset",
    sex_female_pct = 37.4,
    disease_state  = "Unresectable stage III or IV melanoma",
    dose_range     = "0.3, 3, or 10 mg/kg IV as a 90-minute infusion every 3 weeks for up to 4 induction doses, followed by maintenance every 12 weeks beginning at week 24 in eligible patients",
    regions        = "Multinational phase II studies CA184-007, CA184-008, CA184-022, CA184-004",
    ecog_distribution = "ECOG 0 65.0%, ECOG 1 34.5%, ECOG 2 0.5% (index dataset)",
    renal_function = "Baseline MDRD eGFR mean 86.66 (SD 25.78) mL/min/1.73 m^2 (index dataset)",
    hepatic_function = "Baseline total bilirubin mean 0.48 mg/dL; baseline direct bilirubin mean 0.16 mg/dL; ALT mean 23.65 IU/L (index dataset)",
    immunogenicity = "4.29% of index-dataset patients were ADA-positive at baseline or post-ipilimumab",
    co_medication = "Concomitant budesonide in 13.81% of index-dataset patients (CA184-007 prophylactic-toxicity sub-study)",
    notes          = "Pooled index dataset of 1,767 ipilimumab serum concentrations from 420 patients in three phase II studies (CA184-007, CA184-008, CA184-022) plus an external validation dataset of 328 concentrations from 79 patients in CA184-004. Final parameter estimates were obtained from the combined analysis (index + external) dataset. Baseline LDH mean 326.74 (SD 375.19) IU/L; reference LDH 206 IU/L (approximate median; Table 2 footnote). Demographics from Feng 2014 Table 1; values are for the index dataset (n = 420) used for model development."
  )

  ini({
    # Structural parameters at the reference patient: 80 kg BW and 206 IU/L
    # LDH (Feng 2014 Table 2 footnote: reference values selected to
    # approximate the population medians).
    # CL and Q are reported in the source as L/h and converted to L/day
    # (L/h * 24 = L/day) so this model keeps time in days, matching the
    # convention of the later Bristol-Myers Squibb anti-CTLA-4 / anti-PD-1
    # extractions (Bajaj 2017 nivolumab, Sanghavi 2020 ipilimumab).
    lcl <- log(0.0150 * 24); label("Reference clearance CL_REF at WT = 80 kg, LDH = 206 IU/L (L/day)") # Feng 2014 Table 2: CL_REF = 0.0150 L/h
    lvc <- log(4.15);        label("Reference central volume Vc_REF at WT = 80 kg (L)")               # Feng 2014 Table 2: Vc_REF = 4.15 L
    lq  <- log(0.0411 * 24); label("Reference intercompartmental clearance Q_REF (L/day)")            # Feng 2014 Table 2: Q_REF = 0.0411 L/h
    lvp <- log(3.11);        label("Reference peripheral volume Vp_REF (L)")                          # Feng 2014 Table 2: Vp_REF = 3.11 L

    # Allometric power exponents on baseline body weight (Feng 2014 Table 2
    # and Results "PPK model development" final-model equations).
    e_wt_cl  <- 0.642;  label("Power exponent of WT on CL (unitless)") # Feng 2014 Table 2: CL_BW  = 0.642
    e_wt_vc  <- 0.708;  label("Power exponent of WT on Vc (unitless)") # Feng 2014 Table 2: V_cBW  = 0.708

    # Power exponent on the ratio of logs: (log(LDH)/log(206))^e_logldh_cl.
    # Feng 2014 Results state the LDH value was log-transformed because of
    # its right-skewed distribution; the same modelling group's later
    # ipilimumab popPK (Sanghavi 2020 doi:10.1002/psp4.12477) writes the
    # equivalent term explicitly as (log(BLDH)/log(217))^CL_logBLDH. This
    # literal "power of log ratio" form is consistent with the paper's
    # narrative that the LDH-driven AUCss change is small (Figure 5A) and
    # not clinically meaningful.
    e_logldh_cl <- 1.13; label("Power exponent on log(LDH)/log(206) for CL (unitless)") # Feng 2014 Table 2: CL_LDH = 1.13

    # IIV: log-normal on CL and Vc with a single CL:Vc covariance.
    # Variance / covariance entered in lower-triangular row-major order:
    #   row 1: omega^2_CL                 = 0.125
    #   row 2: cov(CL,Vc), omega^2_Vc     = 0.0254, 0.0223
    # Feng 2014 Table 2 reports the correlation as 0.452, and
    #   0.0254 / sqrt(0.125 * 0.0223) = 0.481 ~= 0.452 (within rounding).
    # IIV on Q and Vp were fixed to zero in the source because the
    # variances could not be reliably estimated (Feng 2014 Results).
    etalcl + etalvc ~ c(0.125,
                        0.0254, 0.0223)                                                          # Feng 2014 Table 2: omega^2_CL = 0.125, cov_CL:Vc = 0.0254, omega^2_Vc = 0.0223

    # Combined proportional + additive residual error model
    # (Feng 2014 Results "PPK model development" and Table 2).
    propSd <- 0.157;  label("Proportional residual error (fraction)") # Feng 2014 Table 2: proportional error = 15.7%
    addSd  <- 0.244;  label("Additive residual error (ug/mL)")        # Feng 2014 Table 2: additive error = 0.244 ug/mL
  })
  model({
    # Individual CL and Vc with covariate adjustments (Feng 2014 Results
    # "PPK model development" final-model equations):
    #   CL_i = CL_REF * (BW/80)^CL_BW * (log(LDH)/log(206))^CL_LDH * exp(eta_CL)
    #   Vc_i = Vc_REF * (BW/80)^Vc_BW * exp(eta_Vc)
    cl <- exp(lcl + etalcl) *
      (WT / 80)^e_wt_cl *
      (log(LDH) / log(206))^e_logldh_cl

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc

    # Q and Vp have no covariate effects and no IIV (Feng 2014 Results;
    # IIV in Q and Vp were fixed to zero because they could not be
    # reliably estimated).
    q  <- exp(lq)
    vp <- exp(lvp)

    # Two-compartment micro-constants and ODEs (linear, time-invariant;
    # Feng 2014 Conclusion: "linear, two compartment, zero order i.v.
    # infusion model").
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, V in L -> central/Vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
