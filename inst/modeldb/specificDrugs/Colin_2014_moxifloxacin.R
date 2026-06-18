Colin_2014_moxifloxacin <- function() {
  description <- "Three-compartment population PK model for moxifloxacin in post-bariatric (roux-en-y gastric bypass) volunteers (Colin 2014): linear first-order absorption (no lag, no transit) into a central compartment with two peripheral compartments, allometric scaling on lean body mass (exponent 0.75 on all CL terms and 1 on all volumes, reference LBM 60 kg), and inter-individual variability on ka, central volume, and clearance. Single 400 mg oral and 400 mg 1-h IV infusion doses are fit simultaneously with an implicit bioavailability of 1."
  reference <- paste(
    "Colin P, Eleveld DJ, Struys MMRF, T'Jollyn H, Van Bortel LM,",
    "Ruige J, De Waele J, Van Bocxlaer J, Boussery K.",
    "Moxifloxacin dosing in post-bariatric surgery patients.",
    "Br J Clin Pharmacol. 2014 Jul;78(1):84-90.",
    "doi:10.1111/bcp.12302. PMID 24330006; PMCID PMC4168384.",
    sep = " "
  )
  vignette <- "Colin_2014_moxifloxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Used for allometric scaling on all clearance terms (CL, Q2, Q3) with exponent 0.75, and on all volume terms (V1, V2, V3) with exponent 1, normalised to a reference LBM of 60 kg (Colin 2014 Table 2 footnote). The absorption rate constant ka is the only structural parameter NOT scaled by LBM (Colin 2014 Table 2 footnote: 'All model parameters except ka were centred for a typical subject with a LBM of 60 kg'). LBM in the source dataset was computed from height, total body weight, and sex using the James 1976 equation; care is recommended when applying the model in subjects with BMI > 30 kg/m^2 because the James equation tends to underestimate true LBM in that range (Colin 2014 Discussion).",
      source_name        = "LBM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "25-57 years (median 41)",
    age_median     = "41 years",
    weight_range   = "TBM 57.4-104.0 kg; LBM 41.9-77.6 kg",
    weight_median  = "TBM 78.1 kg; LBM 51.7 kg (LBM reference 60 kg)",
    sex_female_pct = 66.7,
    race_ethnicity = "Not reported (single-centre study at Ghent University Hospital, Belgium)",
    disease_state  = "Post-bariatric (roux-en-y gastric bypass) volunteers at least 6 months after surgery; not actively infected and in generally good condition",
    dose_range     = "Single 400 mg moxifloxacin administered on two occasions per subject in a randomised crossover with a 1-week washout: once as a tablet (oral) and once as a 1-h intravenous infusion",
    regions        = "Belgium (Ghent University Hospital)",
    notes          = "Demographics from Colin 2014 Table 1 (median [range], N=12; 4 male / 8 female). Baseline serum albumin 3.54-4.40 g/dL, height 1.58-1.99 m, CLcr (Cockroft-Gault) 100.4-221.5 mL/min, CLcr (MDRD) 91.9-134.9 mL/min. Source dataset (432 total plasma concentrations from 12 volunteers over a 0-72 h post-dose sampling window per period) originated from De Smet et al. (J Antimicrob Chemother 2012)."
  )

  # Notes on parameter encoding:
  # * Parameters reported in Colin 2014 Table 2 are point estimates for a
  #   typical subject with LBM = 60 kg, with the 95% bootstrap confidence
  #   interval (100 bootstrap samples) shown in brackets.
  # * Allometric exponents (0.75 on CL terms; 1 on volume terms) were
  #   tested but fixed in the final model. Colin 2014 Discussion: removing
  #   the LBM covariate or the allometric exponent of 0.75 on the CL terms
  #   did not improve fit (dAICc = +56.9 and -0.50 respectively); LBM was
  #   superior to TBM (dAICc = +27.3 when using TBM). The exponents are
  #   wrapped in fixed() to preserve this provenance.
  # * Both IV and oral 400 mg moxifloxacin doses were fit simultaneously
  #   (Colin 2014 Methods, PK analysis). The paper does not report a
  #   separate bioavailability parameter in Table 2, consistent with the
  #   authors' prior NCA finding that moxifloxacin oral bioavailability is
  #   unaltered post-bariatric (Introduction, ref 4). The implicit F = 1
  #   used here is a structural assumption documented in the vignette
  #   Assumptions section.
  # * IIV variances reported in Table 2 are omega^2 values on the log
  #   scale; the diagonal variance-covariance matrix was used (Methods,
  #   PK model building) and only ka, V1, and CL carry random effects.
  # * The residual error sigma^2 = 0.03 is the variance of the proportional
  #   component; propSd = sqrt(0.03) = 0.1732 (17.3% CV). The additive
  #   error term used during model exploration was not retained in the
  #   final reported model (Table 2 lists only sigma^2 (Proportional)).
  ini({
    # Structural PK parameters - typical values at LBM = 60 kg.
    lka  <- log(0.95);  label("Absorption rate constant ka (1/h)")                        # Colin 2014 Table 2 ka = 0.95 [95% boot CI 0.72, 1.21]; only structural parameter without LBM scaling
    lcl  <- log(8.60);  label("Clearance from the central compartment CL (L/h)")          # Colin 2014 Table 2 CL x (LBM/60)^0.75 = 8.60 [7.80, 9.70]
    lvc  <- log(47.7);  label("Central volume of distribution V1 (L)")                    # Colin 2014 Table 2 V1 x (LBM/60)^1 = 47.7 [31.6, 78.6]
    lq   <- log(105.3); label("First inter-compartmental clearance Q2 (L/h)")             # Colin 2014 Table 2 Q2 x (LBM/60)^0.75 = 105.3 [55.2, 140.0]
    lvp  <- log(61.5);  label("First peripheral volume of distribution V2 (L)")           # Colin 2014 Table 2 V2 x (LBM/60)^1 = 61.5 [37.6, 75.7]
    lq2  <- log(1.35);  label("Second inter-compartmental clearance Q3 (L/h)")            # Colin 2014 Table 2 Q3 x (LBM/60)^0.75 = 1.35 [1.23, 1.56]
    lvp2 <- log(48.4);  label("Second peripheral volume of distribution V3 (L)")          # Colin 2014 Table 2 V3 x (LBM/60)^1 = 48.4 [34.4, 92.9]

    # Allometric exponents - fixed in the final model (Colin 2014
    # Discussion; ka not scaled).
    allo_cl <- fixed(0.75); label("Allometric exponent on all CL terms (unitless)")  # Colin 2014 Discussion: removing the 0.75 exponent on CL terms gave dAICc = -0.50 (no improvement); exponent fixed in the final model
    allo_v  <- fixed(1.00); label("Allometric exponent on all V terms (unitless)")   # Colin 2014 Table 2 footnote: V terms scale linearly with (LBM/60)

    # IIV - diagonal omega matrix; only ka, V1, and CL carry random
    # effects (Colin 2014 Methods, PK model building).
    etalka ~ 0.24  # Colin 2014 Table 2 omega^2(ka) = 0.24 (variance on log scale)
    etalvc ~ 0.14  # Colin 2014 Table 2 omega^2(V1) = 0.14 (variance on log scale)
    etalcl ~ 0.04  # Colin 2014 Table 2 omega^2(CL) = 0.04 (variance on log scale)

    # Residual error - proportional only.
    propSd <- 0.1732; label("Proportional residual error (fraction)")  # Colin 2014 Table 2 sigma^2 (Proportional) = 0.03; propSd = sqrt(0.03) = 0.1732
  })

  model({
    # Individual PK parameters. ka has IIV but no LBM scaling; vc and cl
    # carry IIV and LBM scaling; the remaining disposition parameters use
    # the typical value scaled by LBM.
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) * (LBM / 60)^allo_cl
    vc  <- exp(lvc + etalvc) * (LBM / 60)^allo_v
    q   <- exp(lq)           * (LBM / 60)^allo_cl
    vp  <- exp(lvp)          * (LBM / 60)^allo_v
    q2  <- exp(lq2)          * (LBM / 60)^allo_cl
    vp2 <- exp(lvp2)         * (LBM / 60)^allo_v

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment disposition with first-order absorption from a
    # depot. Bioavailability F = 1 is implicit (not estimated by Colin
    # 2014). Oral 400 mg doses are administered to depot; the 1-h IV
    # infusion is administered to central with a fixed rate of 400 mg / 1
    # h = 400 mg/h.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Observation - plasma moxifloxacin concentration in mg/L. With dose
    # in mg and vc in L, central / vc returns mg/L directly.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
