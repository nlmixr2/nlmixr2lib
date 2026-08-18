Bouazza_2025_prednisolone <- function() {
  description <- paste(
    "One-compartment population PK model for prednisolone (the active metabolite of orally",
    "administered prednisone) in 66 paediatric and adult patients with active systemic lupus",
    "erythematosus (Bouazza 2025). First-order absorption with a lag time that represents the",
    "combined absorption and prednisone-to-prednisolone conversion delay. The disposition is",
    "parameterised on UNBOUND prednisolone (CLu/F and Vu/F, allometrically scaled to body",
    "weight with fixed exponents 0.75 and 1), which linearises the kinetics; the observed TOTAL",
    "plasma concentration is then recovered algebraically through a saturable corticosteroid-",
    "binding-globulin plus linear albumin binding model with constants fixed from Petersen 1983.",
    "Bioavailability decreases as a power function of the administered prednisone dose."
  )
  reference <- paste(
    "Bouazza N, Semeraro M, Lui G, et al. Population pharmacokinetic modelling of prednisolone",
    "in systemic lupus erythematosus patients: Analysis of exposure and disease activity.",
    "Br J Clin Pharmacol. 2025;91(10):2854-2864. doi:10.1002/bcp.70103.",
    "Protein-binding constants (Bmax, K1, Kns) fixed from Petersen HH, Andreassen TK, Breiderhoff T,",
    "et al., as cited by Bouazza 2025 Methods section 2.3 (reference 14)."
  )
  vignette <- "Bouazza_2025_prednisolone"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on both unbound clearance and unbound central volume, with",
        "exponents fixed to 0.75 and 1 respectively (Bouazza 2025 Methods section 2.3,",
        "citing Anderson and Holford). Reference weight 70 kg, as printed in the Table 2",
        "parameter labels 'CLU/F (L.h-1.70 kg-1)' and 'VU/F (L.70 kg-1)'. Cohort median",
        "56.1 kg, IQR 46.8-66.5, range 17.3-113.5 (Bouazza 2025 Table 1)."
      ),
      source_name        = "BW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = "Screened but not retained in the final model (Bouazza 2025 Results section 3.2: 'No other additional covariates resulted in a significant decrease in the OFV')."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. The Discussion notes that corticosteroid-binding globulin is about 20% higher in women, but no sex effect on prednisolone PK was detected in this cohort."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Bouazza 2025 Methods section 2.3, Results section 3.2)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bouazza 2025 Methods section 2.3, Results section 3.2)."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (Schwartz formula under 18 years, MDRD at 18 years and over)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened but not retained. The Discussion notes that only two patients had eGFR below 30 mL/min, limiting the power to detect a renal-function effect."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained as a PK covariate (Bouazza 2025 Methods section 2.3, Results section 3.2)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bouazza 2025 Methods section 2.3, Results section 3.2)."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained (Bouazza 2025 Methods section 2.3, Results section 3.2)."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "prednisone",
      units    = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "prednisolone",
      units    = "nmol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 66,
    n_studies      = 1,
    n_observations = "242 plasma prednisolone concentrations; median 3 per patient (range 1-11); 35 (14%) below the limit of quantification and handled as left-censored data",
    age_range      = "6-63 years",
    age_median     = "23 years (IQR 14.2-31.2)",
    weight_range   = "17.3-113.5 kg",
    weight_median  = "56.1 kg (IQR 46.8-66.5)",
    sex_female_pct = 80.3,
    race_ethnicity = c(
      Caucasian            = 29.2,
      `North African`      = 15.4,
      Turkish              = 1.5,
      Asian                = 20.0,
      `Sub-Saharan African` = 18.5,
      `Central/South American` = 1.5,
      Caribbean            = 10.8,
      Mixed                = 3.1
    ),
    disease_state  = "Active systemic lupus erythematosus (34 juvenile-onset, 32 adult-onset; 51.5% newly diagnosed, 48.5% in relapse; 74.2% renal involvement) meeting ACR or SLICC classification criteria",
    dose_range     = "Oral prednisone at least 0.5 mg/kg/day at initiation; median 0.94 mg/kg/day (IQR 0.65-1.06, range 0.38-2.3)",
    regions        = "France (28 paediatric and adult university-hospital centres)",
    co_medication  = "Methylprednisolone bolus 57.6%, mycophenolate mofetil 53.0%, hydroxychloroquine 87.9%",
    notes          = paste(
      "Prospective, observational, multicentre study (NCT03187743, EUDRACT 2017-002050-36),",
      "enrolment April 2018 to January 2022, 3-month follow-up per patient. Baseline",
      "characteristics from Bouazza 2025 Table 1. Sampling was opportunistic: 5 min to 32 h",
      "after the last prednisone dose, median 3.3 h (IQR 2-6.5 h). Estimation used Monolix",
      "2023R1 with the MCMC-SAEM algorithm."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Bouazza 2025 Table 2.
    #
    # Monolix reports omega as the STANDARD DEVIATION of the random
    # effect, whereas nlmixr2's ini() takes the VARIANCE, so every omega
    # below is squared. Table 2's residual "proportional" value is the
    # proportional SD and carries over directly to propSd.
    # ------------------------------------------------------------------
    ltlag <- log(0.17); label("Absorption lag time, combined absorption and prednisone-to-prednisolone conversion delay (h)")  # Table 2: Tlag = 0.17 h (2% rse)
    lka   <- log(1.19); label("First-order absorption rate constant (1/h)")                                                    # Table 2: ka = 1.19 1/h (1% rse)
    lcl   <- log(54.3); label("Unbound prednisolone apparent clearance CLu/F at 70 kg (L/h)")                                  # Table 2: CLU/F = 54.3 L/h/70 kg (7% rse), 95% CI 48-62
    lvc   <- log(235);  label("Unbound prednisolone apparent central volume Vu/F at 70 kg (L)")                                # Table 2: VU/F = 235 L/70 kg (8% rse), 95% CI 203-274

    # Allometric exponents, fixed a priori (Methods 2.3: "power exponents
    # fixed to 0.75 and 1 for CL and V, respectively").
    e_wt_cl <- fixed(0.75); label("Allometric exponent on unbound clearance (unitless)")        # Table 2 covariate column: (BW/70)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on unbound central volume (unitless)")   # Table 2 covariate column: (BW/70)^1

    # Dose effect on bioavailability. Table 2 fixes the typical F at 1 and
    # modifies it by (Dose/80000)^betaFdose, with Dose the administered
    # prednisone dose in nmol (Table 2 footnote: 80 000 nmol is about
    # 30 mg of prednisone). The exponent is negative, i.e. F falls as the
    # dose rises -- Abstract: "the bioavailability parameter was found to
    # decrease non-linearly with the dose".
    e_dose_fdepot <- -0.28; label("Power exponent of administered prednisone dose on bioavailability (unitless)")  # Table 2: betaFdose = -0.28 (2% rse)

    # ------------------------------------------------------------------
    # Plasma protein binding -- all three constants fixed, not estimated
    # (Table 2 marks each with *; Methods 2.3: "Bmax, K1 and Kns values
    # were set according to Petersen et al.").
    #
    #   Ctot = Cu * (1 + Bmax / (1 + K1 * Cu) + Kns)
    #
    # Bmax and Kns are dimensionless ratios and K1 has units of L/nmol, so
    # K1 * Cu is dimensionless with Cu in nmol/L. The saturable term is
    # corticosteroid-binding globulin (high affinity, low capacity) and the
    # linear Kns term is albumin (low affinity, high capacity).
    # ------------------------------------------------------------------
    lbmax <- fixed(log(6.77)); label("Corticosteroid-binding-globulin binding capacity Bmax (unitless)")  # Table 2: Bmax = 6.77, from Petersen et al.
    kaff  <- fixed(0.0095);    label("Corticosteroid-binding-globulin affinity constant K1 (L/nmol)")     # Table 2: K1 = 0.0095 L/nmol, from Petersen et al.
    kns   <- fixed(0.8);       label("Albumin linear non-saturable binding constant Kns (unitless)")      # Table 2: Kns = 0.8, from Petersen et al.

    # ------------------------------------------------------------------
    # Between-subject variability. Table 2 reports omega on the SD scale;
    # squared here for nlmixr2's variance scale. Only CLu and Bmax carry
    # BSV in the final model -- no BSV is reported on Tlag, ka or Vu.
    # ------------------------------------------------------------------
    etalcl  ~ 0.1225           # Table 2: omega CL,U = 0.35 (17% rse); variance = 0.35^2 = 0.1225
    etalbmax ~ fixed(0.09)     # Table 2: omega Bmx = 0.3, marked * in the table; variance = 0.3^2 = 0.09

    propSd <- 0.42; label("Proportional residual error on total prednisolone concentration (fraction)")  # Table 2: residual variability, proportional = 0.42 (7% rse)
  })

  model({
    # 1. Individual parameters. Allometric scaling to a 70 kg reference.
    tlag <- exp(ltlag)
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc   <- exp(lvc) * (WT / 70)^e_wt_vc
    bmax <- exp(lbmax + etalbmax)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. ODE system. The dose is entered as nmol of prednisone; the 1:1
    #    molar conversion to prednisolone is absorbed into the lag time
    #    and the bioavailability term, so the central compartment holds
    #    nmol of prednisolone.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Lag time and dose-dependent bioavailability. podo(depot) is the
    #    unscaled administered amount of the current dose, in nmol, so the
    #    power term reproduces Table 2's (Dose/80000)^betaFdose without
    #    needing a separate dose covariate column that could drift out of
    #    sync with amt. The typical F is fixed at 1 at the 80 000 nmol
    #    reference dose.
    alag(depot) <- tlag
    f(depot)    <- (podo(depot) / 80000)^e_dose_fdepot

    # 5. Observation. The ODE system is parameterised on UNBOUND
    #    prednisolone, so central/vc is the unbound concentration Cu. The
    #    measured quantity is TOTAL plasma prednisolone, recovered with the
    #    saturable-plus-linear binding model of Methods section 2.3:
    #      Ctot = Cu * (1 + Bmax / (1 + K1 * Cu) + Kns)
    #    The proportional residual error applies to the total concentration,
    #    which is what the LC-MS/MS assay reported.
    Cu <- central / vc
    Cc <- Cu * (1 + bmax / (1 + kaff * Cu) + kns)

    Cc ~ prop(propSd)
  })
}
