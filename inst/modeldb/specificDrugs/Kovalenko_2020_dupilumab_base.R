Kovalenko_2020_dupilumab_base <- function() {
  description <- "Dupilumab primary base population PK model from Kovalenko 2020 (Model 3): 2-compartment with parallel linear + Michaelis-Menten elimination and a 3-transit-compartment SC absorption chain; fit to Phase 3 atopic-dermatitis data with only body weight as a covariate of central volume."
  reference <- "Kovalenko P, Davis JD, Li M, et al. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2020;9(6):756-767. doi:10.1002/cpdd.780"
  vignette <- "Kovalenko_2020_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Model 3 (primary BASE model) from Kovalenko 2020 Table 1 and Supplementary
  # Table S3.  Phase 3 atopic-dermatitis cohort only (R668-AD-1334 / SOLO 1,
  # R668-AD-1416 / SOLO 2, R668-AD-1224 / CHRONOS).  The structural model is
  # the same as Model 1 in Kovalenko_2020_dupilumab.R; the difference is that
  # kcp, kpc, ka, MTT, Vm, Km, and F were FIXED to the values obtained from
  # the rich-data fits (Models 1-2), and only Vc, ke, the correlated
  # between-subject random-effect block on (ln Vc, ln ke), and the residual
  # error were re-estimated on Phase 3 data.  Weight is the only covariate
  # (on Vc) in the base model; albumin / BMI / EASI / race are added in the
  # companion covariate model `Kovalenko_2020_dupilumab_covariate.R`
  # (Model 4).
  #
  # Reference body weight is 75 kg, inherited from Kovalenko 2016
  # (doi:10.1002/psp4.12136, Eq. 1: V2 = theta1 * (WT/75)^theta2).  The 2020
  # paper does not re-state the reference weight; the assumption matches the
  # sibling Kovalenko_2020_dupilumab.R (Model 1) file in this package.
  #
  # IIV: the paper's Methods explicitly define omega as the standard
  # deviation of the between-subject random effect; Supplementary Table S3
  # reports the per-parameter omegas as "sigma of ln(Vc)" = 0.216 and
  # "sigma of ln(ke)" = 0.301 with Corr(ln(ke), ln(Vc)) = -0.373.  The
  # nlmixr2 `~ c(...)` block syntax stores the lower-triangular
  # variance-covariance matrix directly; the inline arithmetic in ini()
  # squares the SDs into variances and multiplies SD * SD * Corr into the
  # covariance.
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume Vc with reference 75 kg.  Reference weight is not restated in Kovalenko 2020; the 75 kg value is inherited from Kovalenko 2016 (doi:10.1002/psp4.12136, Eq. 1) and matches the sibling Kovalenko_2020_dupilumab.R (Model 1) file.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = "Phase 3 cohort only.  The article reports the pooled total across all 16 studies as N = 2115 on study and 2041 on active treatment, with 18,243 of 20,809 samples included.  Phase 3 SOLO 1 (R668-AD-1334) N = 447, SOLO 2 (R668-AD-1416) N = 472, CHRONOS (R668-AD-1224) N = 424 per Supplementary Table S1 (PK analysis set).",
    n_studies      = 3L,
    age_range      = "Adults with moderate-to-severe atopic dermatitis (detailed age breakdown not reported in the main text).",
    age_median     = "Not reported in the main text.",
    weight_range   = "Not reported in the main text.",
    weight_median  = "Not reported in the main text; reference weight 75 kg is the Kovalenko 2016 inheritance.",
    sex_female_pct = "Not reported in the main text.",
    race_ethnicity = "Race was a tested covariate; detailed breakdown not reported in the main text.  Phase 3 AD trials are predominantly White, with Black, Asian, and Other categories represented.",
    disease_state  = "Adults with moderate-to-severe atopic dermatitis (SOLO 1, SOLO 2 monotherapy; CHRONOS concomitant topical corticosteroids).",
    dose_range     = "600 mg SC loading dose on day 1 followed by 300 mg SC qw or q2w for 15 weeks (SOLO 1, SOLO 2) or 51 weeks (CHRONOS).",
    regions        = "Multi-regional Phase 3 programme; see Supplementary Table S1 for per-study geographic coverage.",
    notes          = "Base model (Model 3) is the precursor to the primary covariate model (Model 4).  The structural parameters kcp, kpc, ka, MTT, Vm, Km, and F are FIXED to values obtained from Models 1 and 2 fits on rich Phase 1/2 data (per the stepwise modelling strategy described in the paper's Methods).  Shrinkage in SD of etas for ke and Vc was 22% and 19%, respectively (Results)."
  )

  ini({
    # Estimated structural parameters on Phase 3 data (Supplementary Table S3, Model 3 column)
    lvc  <- log(2.76);   label("central volume at 75 kg (L)")                         # Supp. Table S3 Model 3: Vc = 2.76 L (SE 0.021)
    lkel <- log(0.0448); label("linear elimination rate constant (1/d)")              # Supp. Table S3 Model 3: ke = 0.0448 1/d (SE 0.000490)

    # Fixed structural parameters (carried forward from Models 1 and 2 fits)
    lkcp    <- fixed(log(0.211));  label("central-to-peripheral rate kcp (1/d)")     # Supp. Table S3 Model 3 / Table 1: kcp = 0.211 (fixed; estimated in Model 1)
    lkpc    <- fixed(log(0.310));  label("peripheral-to-central rate kpc (1/d)")     # Supp. Table S3 Model 3 / Table 1: kpc = 0.310 (fixed; derived from Model 1 kcp/Mpc)
    lka     <- fixed(log(0.306));  label("absorption rate ka (1/d)")                  # Supp. Table S3 Model 3 / Table 1: ka  = 0.306 (fixed; estimated in Model 2)
    lmtt    <- fixed(log(0.105));  label("mean transit time MTT (d)")                 # Supp. Table S3 Model 3 / Table 1: MTT = 0.105 (fixed; estimated in Model 1)
    lvmax   <- fixed(log(1.07));   label("maximum target-mediated rate of elimination Vmax (mg/L/d)") # Supp. Table S3 Model 3 / Table 1: Vm = 1.07 (fixed; estimated in Model 2)
    Km      <- fixed(0.01);        label("Michaelis-Menten constant Km (mg/L)")       # Supp. Table S3 Model 3 / Table 1: Km = 0.01 (fixed; carried over from Kovalenko 2016)
    lfdepot <- fixed(log(0.642));  label("subcutaneous bioavailability F (fraction)") # Supp. Table S3 Model 3 / Table 1: F  = 0.642 (fixed; estimated in Model 1)

    # Covariate effect: power of WT/75 on Vc
    e_wt_vc <- 0.919; label("Power exponent of WT/75 on Vc (unitless)")               # Supp. Table S3 Model 3: Vc ~ weight = 0.919 (SE 0.027)

    # Inter-individual variability: correlated (ln Vc, ln ke) block.
    # Supp. Table S3 Model 3 reports omegas as SDs on the log scale:
    #   SD(ln Vc) = 0.216, SD(ln ke) = 0.301
    #   Corr(ln(ke), ln(Vc)) = -0.373
    # Variances: 0.216^2 = 0.046656, 0.301^2 = 0.090601
    # Covariance: -0.373 * 0.216 * 0.301 = -0.024254
    etalvc + etalkel ~ c(0.046656,
                         -0.024254, 0.090601)

    # Residual error (combined proportional + additive on Cc in mg/L)
    propSd <- 0.124; label("Proportional residual error (fraction)")                  # Supp. Table S3 Model 3: sigma_prop (CV%) = 12.4 (SE 0.18)
    addSd  <- 6.17;  label("Additive residual error (mg/L)")                          # Supp. Table S3 Model 3: sigma_add (mg/L) = 6.17 (SE 0.23)
  })
  model({
    # Individual structural parameters (typical reference: 75 kg)
    vc   <- exp(lvc  + etalvc)  * (WT / 75)^e_wt_vc
    kel  <- exp(lkel + etalkel)
    kcp  <- exp(lkcp)
    kpc  <- exp(lkpc)
    ka   <- exp(lka)
    MTT  <- exp(lmtt)
    vmax <- exp(lvmax)

    # Transit-chain rate for 3 transit compartments
    ktr <- (3 + 1) / MTT

    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * (depot    - transit1)
    d/dt(transit2)    <-  ktr * (transit1 - transit2)
    d/dt(transit3)    <-  ktr * transit2 - ka * transit3
    d/dt(central)     <-  ka * transit3 -
                          kel * central -
                          kcp * central + kpc * peripheral1 -
                          central * (vmax / (Km + central / vc))
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1

    # Subcutaneous bioavailability on the depot; IV doses bypass the depot via the event record
    f(depot) <- exp(lfdepot)

    # mg dosed / L central volume gives mg/L concentration; no unit conversion required
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
