Onita_2025_sulbactam <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous sulbactam in 122",
    "pediatric patients (4 weeks to 16.4 years), built by pooled analysis of",
    "690 plasma concentrations digitised from 23 published pediatric sulbactam",
    "PK studies identified by a MEDLINE search (Onita 2025). Clearance carries",
    "two separate allometric power terms, both with a fixed exponent of 0.75:",
    "one on body weight normalised to the cohort median 22.45 kg, and one on",
    "age normalised to the cohort median 8 years; the age term is the novel",
    "feature of this model and drives an order-of-magnitude lower weight-",
    "normalised clearance in infants than in adolescents. Central volume is",
    "linear in body weight (exponent 1); intercompartmental clearance and",
    "peripheral volume carry no covariates. Exponential interindividual",
    "variability is estimated on all four structural parameters and residual",
    "error is proportional. The model also returns the unbound plasma",
    "concentration Cu (71.2% free fraction), which indexes the paper's",
    "60% fT > MIC probability-of-target-attainment analysis against",
    "Acinetobacter baumannii."
  )
  reference <- paste(
    "Onita T, Sano Y, Ikawa K, Ishihara N, Tamaki H, Yano T. Population",
    "pharmacokinetic analysis and pharmacodynamic evaluation of sulbactam in",
    "pediatric patients: dosing suggestions for Acinetobacter baumannii",
    "infections. J Pediatric Infect Dis Soc. 2025;14(5):piaf043.",
    "doi:10.1093/jpids/piaf043.",
    "All structural, covariate, interindividual-variability and residual-error",
    "estimates are from Supplementary Table S1; the free fraction and the",
    "60% fT > MIC target are from the Methods section",
    "'Estimation of the PK/PD Parameter'.",
    sep = " "
  )
  vignette <- "Onita_2025_sulbactam"
  # Doses are in mg and volumes in L, so central/vc is in mg/L. The paper
  # reports MICs and concentrations in ug/mL, which is the same number
  # (1 mg/L = 1 ug/mL), so Cc and Cu compare directly against the CLSI
  # breakpoints without any scaling.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power term on CL with a fixed exponent of 0.75, and a",
        "linear (exponent 1) term on central volume, both normalised to the",
        "cohort median of 22.45 kg (Table S1 footnote a). Cohort range 4-77 kg,",
        "mean 24.7 kg (Table 1). Not a covariate on Q or peripheral volume."
      ),
      source_name        = "BWT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power term on CL with a fixed exponent of 0.75, normalised",
        "to the cohort median of 8 years (Table S1 footnote b). Cohort range",
        "0.083-16.42 years, mean 7.5 years (Table 1). Age acting on CL is the",
        "novel feature of this model relative to prior sulbactam analyses",
        "(Discussion). Not a covariate on any volume term."
      ),
      source_name        = "AGE"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 122L,
    n_studies      = 23L,
    n_observations = 690L,
    age_range      = "0.083-16.42 years",
    age_median     = "8 years",
    age_mean       = "7.5 years (SD 4.0)",
    weight_range   = "4-77 kg",
    weight_median  = "22.45 kg",
    weight_mean    = "24.7 kg (SD 13.4)",
    sex_female_pct = 40.2,
    age_strata     = c(
      "infant (4 weeks to 11 months)" = 10L,
      "child (1-6 years)"             = 44L,
      "pediatric (7-16 years)"        = 68L
    ),
    disease_state  = paste(
      "Not reported. The pooled source publications did not record underlying",
      "disease, treatment indication, or renal function (serum creatinine,",
      "blood urea nitrogen); the paper lists this as its principal limitation."
    ),
    dose_range     = "4.9-41.7 mg/kg IV, given as a bolus or a 0.5-hour or 1-hour infusion",
    regions        = "Not reported; source publications identified by MEDLINE search",
    notes          = paste(
      "Baseline demographics are Table 1 of Onita 2025 (N = 122; sex",
      "male/female/not-applicable 69:49:4, so the female percentage is",
      "49/122 = 40.2%). Concentrations were measured by bioassay in 20 of the",
      "23 pooled studies and by an unclear method in the remaining 3."
    )
  )

  ini({
    # ================================================================
    # Structural parameters. Onita 2025 Table S1, "Fix effects parameter".
    # Values are the typical-individual estimates at the reference
    # covariates (WT = 22.45 kg, AGE = 8 years), where both allometric
    # multipliers equal 1.
    #
    # Cross-check against the Discussion, which states "CL: 0.31 (L/h/kg),
    # V (Vc + Vp): 0.39 (L/kg) as mean body weight 24.7 kg":
    #   7.65 / 24.7          = 0.3097 L/h/kg  -> 0.31  (matches)
    #   (4.75 + 4.88) / 24.7 = 0.3899 L/kg    -> 0.39  (matches)
    # Both identities hold only if the tabulated values are the estimates
    # at the reference covariates, which confirms the parameterisation.
    # ================================================================
    lcl <- log(7.65)
    label("Clearance at WT 22.45 kg and AGE 8 years (L/h)")  # Onita 2025 Table S1: "thetaCL(L/h) 7.65 (5.84) 6.58-8.47"

    lvc <- log(4.75)
    label("Central volume of distribution at WT 22.45 kg (L)")  # Onita 2025 Table S1: "thetaVcentral (L) 4.75 (7.5) 4.13-5.79"

    lq <- log(3.03)
    label("Intercompartmental clearance (L/h)")  # Onita 2025 Table S1: "Q (L/h) = thetaQ  3.03 (11.1) 2.19-3.76"; no covariate on Q

    lvp <- log(4.88)
    label("Peripheral volume of distribution (L)")  # Onita 2025 Table S1: "Vperipheral (L) = thetaVperipheral  4.88 (10.4) 3.88-6.91"; no covariate on Vp

    # ================================================================
    # Covariate-effect exponents. All three are fixed structural
    # assumptions, not estimates: Table S1 prints them inline in the
    # parameter equations with no estimate, RSE or confidence interval of
    # their own, and the Methods state "A standard allometric scaling
    # exponent of 0.75 was included for the clearance and volume of
    # distribution into the model."
    #
    # Table S1 prints:
    #   CL (L/h)       = thetaCL       x (BWT/22.45)^0.75 x (AGE/8)^0.75
    #   Vcentral (L)   = thetaVcentral x (BWT/22.45)
    # so CL carries TWO 0.75 exponents (one per covariate) and central
    # volume carries an exponent of 1, not 0.75. Where the Methods prose
    # ("...for the clearance and volume of distribution") conflicts with
    # the printed equations, the printed equations govern; the Methods'
    # own generic covariate forms agree with Table S1:
    #   CLi(k) = CL(k) x (xi/median(x))^0.75
    #   Vi(k)  = V(k)  x (xi/median(x))
    # ================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on body weight for CL (unitless)")  # Onita 2025 Table S1 equation "(BWT/22.45)^0.75"; Methods: "A standard allometric scaling exponent of 0.75 was included"

    e_age_cl <- fixed(0.75)
    label("Allometric exponent on age for CL (unitless)")  # Onita 2025 Table S1 equation "(AGE/8)^0.75"

    e_wt_vc <- fixed(1)
    label("Exponent on body weight for central volume (unitless)")  # Onita 2025 Table S1 equation "Vcentral (L) = thetaVcentral x (BWT/22.45)" -- printed with no exponent, i.e. 1; matches Methods "Vi(k) = V(k) x (xi/median(x))"

    # ================================================================
    # Reference (normalising) covariate values. Table S1 footnotes give
    # these as the cohort medians in the 122 subjects, not rounded
    # standards.
    # ================================================================
    wtref <- fixed(22.45)
    label("Reference body weight for allometric normalisation (kg)")  # Onita 2025 Table S1 footnote a: "Median as body weight in 122 subjects"

    ageref <- fixed(8)
    label("Reference age for allometric normalisation (years)")  # Onita 2025 Table S1 footnote b: "Median as age in 122 subjects"

    # ================================================================
    # Interindividual variability. Onita 2025 Table S1, "Interindividual
    # variability (exponential error model)". The tabulated numbers are
    # VARIANCES: the Methods define the model as theta_i = theta x
    # exp(eta_i) with eta "normally distributed with a mean of 0 and
    # variance of omega^2", and the Table S1 footnote repeats "eta, random
    # variable which is normally distributed with a mean of zero and
    # variance". They are therefore passed through unsquared. No
    # off-diagonal covariances were reported, so the omega matrix is
    # diagonal.
    # ================================================================
    etalcl ~ 0.229   # Onita 2025 Table S1, row etaCL: estimate 0.229, RSE 25.1%, 95% CI 0.134-0.388
    etalvc ~ 0.223   # Onita 2025 Table S1, row etaVcentral: estimate 0.223, RSE 43.1%, 95% CI 0.0430-0.431
    etalq ~ 0.430    # Onita 2025 Table S1, row etaQ: estimate 0.430, RSE 28.3%, 95% CI 0.184-0.690
    etalvp ~ 0.504   # Onita 2025 Table S1, row etaVperipheral: estimate 0.504, RSE 27.0%, 95% CI 0.275-0.964

    # ================================================================
    # Residual variability. Onita 2025 Table S1, "Residual variability
    # (proportional error model)". As for the etas, the tabulated 0.0699
    # is a VARIANCE (footnote: "epsilon, random error which is normally
    # distributed with a mean of zero and variance"), so the standard
    # deviation nlmixr2 expects is sqrt(0.0699) = 0.2644, i.e. 26.4%
    # proportional error.
    # ================================================================
    propSd <- 0.2644
    label("Proportional residual error (fraction)")  # Onita 2025 Table S1: "epsilon 0.0699 (19.8) 0.0421-0.102"; sqrt(0.0699) = 0.2644

    # ================================================================
    # Pharmacodynamic evaluation inputs. Fixed literature / design
    # constants used by the paper's probability-of-target-attainment
    # analysis, not estimated from the concentration data.
    # ================================================================
    fu <- fixed(0.712)
    label("Unbound fraction of sulbactam in plasma (unitless)")  # Onita 2025 Methods, "Estimation of the PK/PD Parameter": "a value of 71.2% nonprotein-binding rate (the free fraction f) of sulbactam was used"
  })

  model({
    # ------------------------------------------------------------
    # 1. Individual parameters with covariate effects
    #    Onita 2025 Table S1 parameter equations.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / wtref)^e_wt_cl * (AGE / ageref)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / wtref)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # ------------------------------------------------------------
    # 2. Micro-constants
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------
    # 3. ODE system. Onita 2025 Methods, "Population PK Modeling",
    #    reproduced verbatim from the paper (NONMEM ADVAN3):
    #      dXc/dt = -(CL/Vc + Q/Vc) x Xc + Q x Xp/Vp
    #      dXp/dt =  Q x Xc/Vc - Q x Xp/Vp
    #    which is the standard two-compartment system in amounts, with
    #    drug entering the central compartment directly (IV bolus or
    #    infusion; there is no absorption phase and no bioavailability
    #    term in this model).
    # ------------------------------------------------------------
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ------------------------------------------------------------
    # 4. Observation and error.
    #
    #    Cu is the unbound plasma concentration that indexes the paper's
    #    pharmacodynamic target. Onita 2025 Methods, "Estimation of the
    #    PK/PD Parameter": "Fixed-effects parameters were used to
    #    simulate the unbound drug concentrations, where a value of 71.2%
    #    nonprotein-binding rate (the free fraction f) of sulbactam was
    #    used." It is returned as a derived column, not a second
    #    endpoint: only total plasma concentration was measured and
    #    fitted, so no residual error attaches to Cu.
    #
    #    fT > MIC is deliberately NOT accumulated in an ODE state. The
    #    paper computes it by post-processing the simulated profile --
    #    "The time point at which the drug concentrations coincided with
    #    a specific MIC value was determined, and T > MIC was calculated
    #    as the cumulative percentage of 24 hours for different dosing
    #    intervals" -- so the vignette reproduces it the same way, from
    #    Cu on a dense steady-state grid.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
