Nath_2010_melphalan_unbound <- function() {
  description <- "Two-compartment IV-infusion population PK model for unbound (ultrafiltrate) plasma melphalan in adults with multiple myeloma undergoing high-dose therapy and autologous stem-cell transplant (Nath 2010); additive non-renal + renal CL with hematocrit, fat-free mass, and creatinine-clearance covariates."
  reference <- "Nath CE, Shaw PJ, Trotman J, Zeng L, Duffull SB, Hegarty G, McLachlan AJ, Gurney H, Kerridge I, Kwan YL, Presgrave P, Tiley C, Joshua D, Earl J. Population pharmacokinetics of melphalan in patients with multiple myeloma undergoing high dose therapy. Br J Clin Pharmacol. 2010;69(5):484-497. doi:10.1111/j.1365-2125.2010.03638.x"
  vignette <- "Nath_2010_melphalan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on non-renal CL (fixed exponent 0.75) and V1 (fixed exponent 1) with reference 50 kg (Nath 2010 final equations p. 489). Note: the screening grid in Table 4 used FFM/53 (the cohort median per Table 2); the final-model equations on p. 489 use FFM/50, which we follow verbatim. Computed per Janmahasatian et al. Clin Pharmacokinet 2005;44:1051-1065 from total body weight, height, and sex.",
      source_name        = "FFM"
    ),
    HCT = list(
      description        = "Hematocrit on the day of melphalan administration (or closest prior value if unavailable)",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on non-renal CL with reference 34 % (cohort median per Nath 2010 Table 2). Population range 20-45 %.",
      source_name        = "HCT"
    ),
    CRCL = list(
      description        = "Estimated creatinine clearance by Cockcroft and Gault, normalized to a standard weight of 70 kg",
      units              = "mL/min/70 kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper applied Cockcroft-Gault using actual body weight, then divided by total body weight and multiplied by 70 to express the result as mL/min/70 kg (Nath 2010 Methods, page 486). Linear effect on renal CL with reference 88 mL/min/70 kg (cohort median per Table 2).",
      source_name        = "CLcr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 1L,
    n_observations = 691L,
    age_range      = "36-73 years",
    age_median     = "57 years",
    weight_range   = "42-132 kg",
    weight_median  = "78 kg",
    sex_female_pct = 41,
    race_ethnicity = "Not reported",
    disease_state  = "Multiple myeloma undergoing high-dose melphalan and autologous stem cell transplantation",
    dose_range     = "150-450 mg IV (median 368 mg; 115-216 mg/m^2, median 192 mg/m^2); infusion over a median 35 min (range 15-95 min)",
    regions        = "Australia (six participating hospitals; ACTRN0126000231549)",
    ffm_range      = "34.4-80.5 kg (median 53.3)",
    bsa_range      = "1.3-2.6 m^2 (median 1.9)",
    crcl_range     = "26-205 mL/min/70 kg (median 88); 29-234 mL/min (median 97)",
    hct_range      = "20-45 % (median 34); 35 % of patients had HCT < 33 %",
    notes          = "Nath 2010 Table 2 baseline characteristics. Unbound melphalan concentrations were measured in 5-6 samples per patient (n = 691 total observations), timed according to the optimal-design windows in Table 1. Total plasma melphalan concentrations from the same cohort are modelled separately (Nath_2010_melphalan_total). Linearity of melphalan plasma-protein binding over the observed concentration range was verified (Table 8 one-way ANOVA, NS)."
  )

  ini({
    # Structural parameters from Nath 2010 Table 5 (final covariate model with
    # CL_NR + CL_R parameterisation and CRCL in mL/min/70 kg). Final-model
    # equations for unbound melphalan are on p. 489:
    #   CL_NR = 79.7 * (HCT/34)^0.679 * (FFM/50)^0.75
    #   CL_R  = 50.7 * (CRCL/88)
    #   V1    = 63.8 * (FFM/50)
    #   Q     = 152
    #   V2    = 71.6
    lcl_nonren <- log(79.7);  label("Non-renal melphalan clearance at reference HCT 34 % and FFM 50 kg (L/h)")  # Table 5: theta1 = 79.7 (95% CI 64.8-97.3)
    lcl_renal  <- log(50.7);  label("Renal melphalan clearance at reference CRCL 88 mL/min/70 kg (L/h)")        # Table 5: CLR = 50.7 (95% CI 34.8-66.1)
    lvc        <- log(63.8);  label("Central volume of distribution at reference FFM 50 kg (L)")               # Table 5: V1 = 63.8
    lq         <- log(152);   label("Intercompartmental clearance (L/h)")                                       # Table 5: Q = 152
    lvp        <- log(71.6);  label("Peripheral volume of distribution (L)")                                    # Table 5: V2 = 71.6

    # Covariate-effect parameters (final model, p. 489)
    e_hct_cl_nonren  <- 0.679;        label("HCT power exponent on non-renal CL (reference 34 %)")              # Table 5: theta2 = 0.679 (95% CI 0.284-1.070)
    e_ffm_cl_nonren  <- fixed(0.75);  label("FFM power exponent on non-renal CL (fixed allometric value)")      # Methods p. 486: allometric exponent fixed at 0.75
    e_ffm_vc         <- fixed(1.0);   label("FFM power exponent on V1 (fixed allometric value)")                # Methods p. 486: allometric exponent fixed at 1.0

    # Inter-individual variability. CV% reported in Table 5; omega^2 = log(CV^2 + 1).
    # The exponential random effect applies to the typical-value PK parameter
    # (Methods p. 486). For CL, the single eta multiplies the total typical CL
    # (CL_NR + CL_R), as Table 5 reports a single omega_CL for unbound melphalan.
    # CL: CV 29.8% -> omega^2 = log(1 + 0.298^2) = 0.08509
    # V1: CV 38.7% -> omega^2 = log(1 + 0.387^2) = 0.13957
    # Q : CV 49.6% -> omega^2 = log(1 + 0.496^2) = 0.22000
    # V2: CV 35.4% -> omega^2 = log(1 + 0.354^2) = 0.11811
    etalcl ~ 0.08509  # Table 5: omega_CL CV 29.8%
    etalvc ~ 0.13957  # Table 5: omega_V1 CV 38.7%
    etalq  ~ 0.22000  # Table 5: omega_Q  CV 49.6%
    etalvp ~ 0.11811  # Table 5: omega_V2 CV 35.4%

    # Residual error: combined proportional + additive (Methods p. 486:
    # Y = Yhat * (1 + epsilon1) + epsilon2). SDs in Table 5.
    propSd <- 0.138;  label("Proportional residual error (fraction)")    # Table 5: sigma1 = 0.138
    addSd  <- 0.027;  label("Additive residual error (mg/L)")           # Table 5: sigma2 = 0.027
  })
  model({
    # Covariate-driven typical-value clearance components (Nath 2010 p. 489).
    cl_nonren <- exp(lcl_nonren) * (HCT / 34)^e_hct_cl_nonren * (FFM / 50)^e_ffm_cl_nonren
    cl_renal  <- exp(lcl_renal)  * (CRCL / 88)

    # Individual PK: single eta on total CL, separate etas on V1/Q/V2.
    cl <- (cl_nonren + cl_renal) * exp(etalcl)
    vc <- exp(lvc + etalvc) * (FFM / 50)^e_ffm_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Linear two-compartment IV model with first-order elimination from central
    # (Nath 2010 base model definition, p. 486). Unbound (ultrafiltrate)
    # melphalan concentrations were measured by HPLC after ultrafiltration of
    # plasma samples (Nath 2010 Methods, p. 486).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Unbound plasma melphalan concentration: dose mg / volume L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
