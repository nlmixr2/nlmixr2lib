Muntau_2017_sapropterin <- function() {
  description <- "One-compartment population PK model with first-order oral absorption, an absorption lag, linear elimination, and an additive endogenous BH4 baseline for sapropterin dihydrochloride in pediatric patients <4 years with BH4-responsive phenylketonuria or mild hyperphenylalaninemia (Muntau 2017 SPARK trial)."
  reference <- "Muntau AC, Burlina A, Eyskens F, et al. Efficacy, safety and population pharmacokinetics of sapropterin in PKU patients <4 years: results from the SPARK open-label, multicentre, randomized phase IIIb trial. Orphanet Journal of Rare Diseases. 2017;12:47. doi:10.1186/s13023-017-0600-x"
  vignette <- "Muntau_2017_sapropterin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; pediatric range 5-20 kg in the SPARK trial).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent clearance and apparent central volume, normalized to a reference adult weight of 70 kg (Table 4 footnote: 'Reference weight (adult male patient)'). Table 3 reports the exponents (0.839 on CL/F and 0.573 on V/F).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 52L,
    n_studies      = 1L,
    age_range      = "2-47 months",
    age_mean       = "21 months (SD 12)",
    weight_range   = "5-20 kg",
    weight_mean    = "11.3 kg (SD 3.0)",
    sex_female_pct = 46.4,
    disease_state  = "BH4-responsive phenylketonuria (PKU) or mild hyperphenylalaninemia (HPA): per Table 2, 21% classical PKU, 32% mild PKU, 46% mild HPA in the ITT population.",
    dose_range     = "10 mg/kg/day oral sapropterin (could be uptitrated to 20 mg/kg/day at week 4 if Phe tolerance had not increased by >20% vs baseline; only 2 of 27 sapropterin-treated patients escalated).",
    regions        = "Europe (Austria, Belgium, Czech Republic, Germany, Italy, Netherlands, Slovakia, Turkey, United Kingdom; 22 sites in 9 countries).",
    trial          = "SPARK (NCT01376908) -- 26-week open-label, multicentre, randomized phase IIIb study.",
    age_strata     = c("<12 months" = 15L, "12-<24 months" = 18L, "24-<48 months" = 23L),
    notes          = "Baseline demographics from Table 2 (n=56 ITT, n=52 with >=1 PK sample contributing to the popPK analysis). Pharmacokinetic sampling was planned by D-optimization with sparse sampling; plasma samples for endogenous BH4 measurement were collected at baseline and sparsely between weeks 5-12 after oral administration of 10 mg/kg/day. The Discussion notes the one-compartment fit yielded concentration profiles 'virtually identical to those from a two-compartment model evaluated in a previous study' (the Qi 2014 pooled 0-50 years analysis in nlmixr2lib as Qi_2014_sapropterin), supporting the parsimonious one-compartment selection at <4 years."
  )

  ini({
    # Structural parameters -- Muntau 2017 Table 3 final model
    # (reference weight 70 kg adult male per Table 4 footnote; oral dose,
    # apparent parameters).
    lka   <- log(0.234); label("First-order absorption rate (ka, 1/h)")           # Muntau 2017 Table 3: Ka = 0.234 h^-1
    lcl   <- log(2780);  label("Apparent clearance at WT = 70 kg (CL/F, L/h)")    # Muntau 2017 Table 3: CL/F = 2780 L/h
    lvc   <- log(3870);  label("Apparent central volume at WT = 70 kg (V/F, L)")  # Muntau 2017 Table 3: V/F = 3870 L
    ltlag <- log(0.342); label("Absorption lag time (tlag, h)")                   # Muntau 2017 Table 3: LAG = 0.342 h
    lc0   <- log(12.6);  label("Endogenous BH4 baseline plasma concentration (C0, ug/L)")  # Muntau 2017 Table 3: C0 = 12.6 ug/L

    # Body-weight covariate effects: power form normalized to a 70 kg reference.
    e_wt_cl <- 0.839; label("Power exponent of body weight on CL/F (unitless, reference 70 kg)")  # Muntau 2017 Table 3: 'Coefficient describing effect of weight on CL/F' = 0.839
    e_wt_vc <- 0.573; label("Power exponent of body weight on V/F (unitless, reference 70 kg)")   # Muntau 2017 Table 3: 'Coefficient describing effect of weight on V/F' = 0.573

    # Inter-individual variability. Table 3 reports IIV_CL and IIV_V2 as %CV
    # and the Pearson correlation between random effects on CL/F and V/F.
    # Translate %CV to the internal log-scale variance via omega^2 = log(1 + CV^2):
    #   omega^2(CL) = log(1 + 0.2298^2) = 0.0515
    #   omega^2(V ) = log(1 + 0.3256^2) = 0.1008
    #   cov(CL,V)   = 0.134 * sqrt(0.0515 * 0.1008) = 0.0097
    etalcl + etalvc ~ c(0.0515, 0.0097, 0.1008)  # Muntau 2017 Table 3: IIV_CL %CV = 22.98, IIV_V2 %CV = 32.56, Corr(CL,V) = 0.134

    # Residual error. Table 3 reports a single residual error of 65.30 %CV
    # without explicit mention of LTBS / log-additive parameterization;
    # encoded here as a proportional (constant-CV) residual on the linear
    # scale. See the vignette's Assumptions and deviations section.
    propSd <- 0.6530; label("Proportional residual error (fraction CV)")  # Muntau 2017 Table 3: Residual error %CV = 65.30
  })

  model({
    # 1. Derived covariate terms -- power-form WT effects on CL/F and V/F
    #    normalized to a 70 kg reference adult (Muntau 2017 Table 3 + Table 4
    #    footnote).
    wt_cl <- (WT / 70)^e_wt_cl
    wt_vc <- (WT / 70)^e_wt_vc

    # 2. Individual PK parameters.
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * wt_cl
    vc   <- exp(lvc + etalvc) * wt_vc
    tlag <- exp(ltlag)
    c0   <- exp(lc0)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. ODE system (oral dose into depot; apparent parameters).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Absorption lag on the depot compartment.
    alag(depot) <- tlag

    # 6. Observation: predicted concentration = drug-derived + endogenous BH4.
    #    Dose in mg, V/F in L -> central/vc is in mg/L; multiply by 1000 to
    #    express Cc in ug/L, the paper's reporting unit for BH4 plasma
    #    concentration.
    Cc <- (central / vc) * 1000 + c0

    Cc ~ prop(propSd)
  })
}
