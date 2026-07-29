Rovei_1982_theophylline <- function() {
  description <- "One-compartment oral PK model for theophylline tablets (Rovei 1982): first-order absorption with lag time in healthy adult volunteers across single oral doses of 125-500 mg."
  reference <- "Rovei V, Chanoine F, Strolin Benedetti M. Pharmacokinetics of theophylline: a dose-range study. Br J Clin Pharmacol. 1982;14(6):769-778. doi:10.1111/j.1365-2125.1982.tb02035.x"
  vignette <- "Rovei_1982_theophylline"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear scaling on CL and Vc; reference 70 kg. Source paper reports CL and Vd per kg of body weight (Table 3); the per-kg parameterization implies a linear (exponent = 1) weight effect, which is reproduced here by (WT/70)^1 with reference weight 70 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "22-35 years",
    age_mean       = "29 years (SD 4)",
    weight_range   = "48-77 kg",
    weight_mean    = "62 kg (SD 9)",
    sex_female_pct = 50,
    race_ethnicity = c(White = 100),
    disease_state  = "Healthy adult Caucasian non-smokers on a xanthine-free diet; normal hepatic and renal function (Rovei 1982 Table 1, Results page 772).",
    dose_range     = "Single oral 125, 250, 375, 500 mg theophylline tablets (Theodel) in a 4-period cross-over (Rovei 1982 Methods, page 770).",
    regions        = "France / Switzerland (study conducted at Hopital Cantonal of Geneva).",
    notes          = "Plasma sampled 0-48 h post-dose at 13 timepoints; urine collected 0-72 h. Data fitted to a one-compartment open model with first-order absorption and lag time using a Gauss-Newton iterative procedure (G-PHARM, Gomeni & Gomeni 1978)."
  )

  ini({
    # Structural parameters -- pooled across all four doses (Rovei 1982 Table 3, page 774).
    # The paper concludes (page 773, ANOVA on Table 3) that t_lag, t_abs, t_max, t_beta, CL, CL_R, Vd and F are
    # not modified by dose, so across-dose averages are reported here as the typical-value population estimates.
    ltlag <- log(0.09); label("Absorption lag time (hr)")                                              # Rovei 1982 Table 3: across-dose mean of t_lag (0.11, 0.09, 0.09, 0.07 h at 125/250/375/500 mg)
    lka  <- log(1.73); label("First-order absorption rate constant (1/hr)")                           # Rovei 1982 Table 3: ka = ln(2)/t_abs with across-dose mean t_abs = 0.40 h (0.31, 0.44, 0.29, 0.55)
    lcl  <- log(2.94); label("Apparent clearance for a 70 kg adult (L/hr)")                           # Rovei 1982 Table 3: across-dose mean CL = 0.042 L/h/kg * 70 kg = 2.94 L/h
    lvc  <- log(35.7); label("Apparent central volume for a 70 kg adult (L)")                         # Rovei 1982 Table 3: across-dose mean Vd = 0.51 L/kg * 70 kg = 35.7 L

    # IIV approximated from inter-subject SD/mean across the 8 individual fits in Rovei 1982 (Table 3
    # ranges and Results page 773 dose-stratified mean +/- SD); these are not formal popPK omega
    # estimates. Conversion: omega^2 = log(1 + CV^2). Per-parameter rationale recorded in the
    # vignette Source-trace and Assumptions and deviations sections.
    etaltlag ~ log(1 + 0.50^2)
    etalka  ~ log(1 + 0.50^2)
    etalcl  ~ log(1 + 0.25^2)
    etalvc  ~ log(1 + 0.14^2)

    propSd <- 0.10; label("Proportional residual error (fraction)")                                    # Rovei 1982 Methods page 771 / 772: HPLC assay recovery 85 +/- 5% (CV ~6%); 0.10 chosen as a conservative simulation-only default since the paper does not report a residual error model
  })

  model({
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka  + etalka)
    cl   <- exp(lcl  + etalcl) * (WT / 70)
    vc   <- exp(lvc  + etalvc) * (WT / 70)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    lag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
