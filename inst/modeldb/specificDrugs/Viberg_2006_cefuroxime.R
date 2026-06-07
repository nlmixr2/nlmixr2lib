Viberg_2006_cefuroxime <- function() {
  description <- "Two-compartment population PK model for intravenous cefuroxime in adult patients with bacterial infections and a wide range of renal function (Viberg 2006); reciprocal serum cystatin C (1/CYSC) and body weight enter as centred-linear covariates on clearance, and body weight enters as a centred-linear covariate on the central volume of distribution."
  reference <- paste(
    "Viberg A, Lannergard A, Larsson A, Cars O, Karlsson MO, Sandstrom M.",
    "A population pharmacokinetic model for cefuroxime using cystatin C",
    "as a marker of renal function.",
    "Br J Clin Pharmacol. 2006;62(3):297-303.",
    "doi:10.1111/j.1365-2125.2006.02652.x"
  )
  vignette <- "Viberg_2006_cefuroxime"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CYSC = list(
      description        = "Serum cystatin C concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model via the reciprocal transform 1/CYSC, centred at 0.758 (mg/L)^-1 (Viberg 2006 Table 4 footnote). The reciprocal value 0.758 corresponds to CYSC = 1.32 mg/L (close to the population median of 1.18-1.32 across the four renal-function strata in Table 1). Time-fixed at baseline -- Viberg 2006 found that the first measurement alone was sufficient and that adding the second measurement did not improve the fit (Results paragraph 3).",
      source_name        = "CysC"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred-linear effect on CL and V1 with reference weight 74 kg (Viberg 2006 Table 4 footnote). Time-fixed at baseline. For the three patients with no recorded weight, Viberg 2006 imputed the population median.",
      source_name        = "WT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 97L,
    n_studies        = 1L,
    n_observations   = 409L,
    age_range        = "approximately 24-95 years (range across the four CLcr strata in Table 1; medians 56, 74, 82, 78 years for CLcr > 80, 41-80, 21-40, < 20 mL/min)",
    weight_range     = "approximately 35-137 kg (range across the four CLcr strata in Table 1; medians 85, 74, 70, 68 kg)",
    sex_female_pct   = 43.3,
    disease_state    = "Hospitalised adult patients with symptoms and signs indicating bacterial infection believed to be susceptible to cefuroxime, deliberately recruited to span a wide range of renal capacity (CLcr 6-115 mL/min by Cockcroft-Gault).",
    renal_function   = "CLcr range 6-115 mL/min. Stratified into four cohorts (CLcr > 80, 41-80, 21-40, < 20 mL/min) with corresponding dose adjustments per Viberg 2006 Table 1; the < 20 mL/min stratum (n = 10) carried the highest cystatin C values (median 4.51 mg/L).",
    cohort_strata    = c(`CLcr_gt_80_n` = 20L, `CLcr_41_80_n` = 40L, `CLcr_21_40_n` = 27L, `CLcr_lt_20_n` = 10L),
    dose_range       = "Intravenous cefuroxime 750-1500 mg per dose, administered as a 5-15 min intravenous injection. Frequency individualised by renal function: 1500 mg x 3/day or 750 mg x 3/day in CLcr > 80 patients; 750 mg x 3/day in CLcr 41-80; 750 mg x 2/day in CLcr 21-40; 750 mg x 1/day in CLcr < 20 (Viberg 2006 Table 1).",
    regions          = "Sweden -- Uppsala University Hospital (Departments of Infectious Diseases and Nephrology) and Karlstad Central Hospital (Department of Nephrology).",
    exclusions       = "Haemodialysis patients, patients with chronic inflammatory diseases, and patients who had received cefuroxime in the previous two weeks were excluded.",
    notes            = "Baseline demographics per Viberg 2006 Table 1. Final-model parameter estimates per Viberg 2006 Table 4. Modelling performed in NONMEM VI beta with FOCE on log-transformed data (Methods 'Data analysis')."
  )

  ini({
    # Structural parameters -- typical values for a 74 kg patient with the
    # population-typical 1/CYSC = 0.758 (mg/L)^-1 (CYSC ~ 1.32 mg/L).
    # Two-compartment ADVAN3 TRANS4 parameterisation: CL, V1 (= Vc),
    # Q (= intercompartmental CL), V2 (= Vp). Intravenous administration
    # (5-15 min injection) -- no absorption phase.
    lcl <- log(6.00); label("Clearance for the 74 kg, 1/CYSC = 0.758 reference (CL, L/h)")          # Viberg 2006 Table 4 (final estimate)
    lvc <- log(11.4); label("Central volume of distribution for the 74 kg reference (V1, L)")       # Viberg 2006 Table 4 (final estimate)
    lvp <- log(5.11); label("Peripheral volume of distribution (V2, L)")                            # Viberg 2006 Table 4 (final estimate)
    lq  <- log(3.65); label("Intercompartmental clearance (Q, L/h)")                                # Viberg 2006 Table 4 (final estimate)

    # Covariate effects -- centred-linear form per Viberg 2006 Table 4 footnote:
    #   CL = 6.00 * (1 + 1.43  * (1/CYSC - 0.758)) * (1 + 0.0108 * (WT - 74))
    #   V1 = 11.4 * (1 + 0.0097 * (WT - 74))
    # The cystatin C effect is parameterised on the reciprocal 1/CYSC scale
    # because Viberg 2006 found that 1/CYSC produced the largest drop in
    # objective function value among all renal-function covariate transforms
    # tested (Results Table 3: -154.0 OFV for 1/CYSC vs -131.4 for CLcr and
    # -75.4 for 1/SCr).
    e_cysc_cl <-  1.43;   label("Centred-linear coefficient for 1/CYSC on CL (per (mg/L)^-1)")     # Viberg 2006 Table 4 (1/cystatin C column, CL row)
    e_wt_cl   <-  0.0108; label("Centred-linear coefficient for WT on CL (per kg)")                # Viberg 2006 Table 4 (Body weight column, CL row)
    e_wt_vc   <-  0.0097; label("Centred-linear coefficient for WT on V1 (per kg)")                # Viberg 2006 Table 4 (Body weight column, V1 row)

    # Inter-individual variability (IIV) -- diagonal omega (Viberg 2006
    # Results paragraph 4: "Different variance/covariance structures of IIV
    # were assessed, but the model did not benefit from any block.").
    # Conversion of CV% to log-scale variance: omega^2 = log(1 + CV^2).
    #   IIV CL : 27% CV -> log(1 + 0.27^2) = 0.07037
    #   IIV V1 : 18% CV -> log(1 + 0.18^2) = 0.03189
    #   IIV V2 : 48% CV -> log(1 + 0.48^2) = 0.20724
    # No IIV was estimated on Q (Viberg 2006 Table 4 reports "- -" for Q).
    etalcl ~ 0.07037                                                                                # Viberg 2006 Table 4 (IIV CL 27%)
    etalvc ~ 0.03189                                                                                # Viberg 2006 Table 4 (IIV V1 18%)
    etalvp ~ 0.20724                                                                                # Viberg 2006 Table 4 (IIV V2 48%)

    # Residual error -- proportional only (Viberg 2006 Results paragraph 2:
    # "The residual error was adequately described by only a proportional
    # component."). 15.5% on the linear concentration scale.
    propSd <- 0.155; label("Proportional residual error (fraction)")                                # Viberg 2006 Table 4 (Proportional error 15.5%)
  })

  model({
    # Reciprocal cystatin C transform and centred-linear covariate factors.
    # The CL covariate model is multiplicative across 1/CYSC and WT,
    # per Viberg 2006 Table 4 footnote.
    cysc_inv  <- 1 / CYSC
    cl_factor <- (1 + e_cysc_cl * (cysc_inv - 0.758)) * (1 + e_wt_cl * (WT - 74))
    vc_factor <- (1 + e_wt_vc  * (WT - 74))

    cl <- exp(lcl + etalcl) * cl_factor
    vc <- exp(lvc + etalvc) * vc_factor
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration in central compartment: dose mg / volume L -> mg/L,
    # matching the source bioanalytical reporting scale.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
