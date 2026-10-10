Lee_2022_burosumab <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "linear elimination for subcutaneous burosumab (a human anti-FGF23 IgG1",
    "monoclonal antibody) in adult and pediatric (1-12 years) patients with",
    "X-linked hypophosphatemia, coupled to a direct Emax model (Hill",
    "coefficient fixed to 1) for absolute serum phosphorus. Apparent",
    "clearance and volume carry estimated allometric exponents on",
    "time-varying body weight referenced to 70 kg, and body weight also",
    "scales the baseline serum phosphorus E0 (negative exponent) and the",
    "maximal increase Emax (positive exponent)."
  )
  reference <- paste(
    "Lee SK, Gosselin NH, Taylor J, Roberts MS, McKeever K, Shi J.",
    "Population Pharmacokinetics and Pharmacodynamics of Burosumab in Adult",
    "and Pediatric Patients With X-linked Hypophosphatemia.",
    "J Clin Pharmacol. 2022;62(1):87-98. doi:10.1002/jcph.1950"
  )
  vignette <- "Lee_2022_burosumab"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight (time-varying)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying body weight, normalised to 70 kg (Lee 2022 Table 3 and",
        "Table 4 footnotes: 'Effect of time-varying WT was centralized using",
        "70 kg to facilitate the integration of adult data'). Power effects",
        "on CL/F (0.912) and V/F (1.05) in the PK model and on E0 (-0.143)",
        "and Emax (0.168) in the PK-PD model. Time-varying weights were used",
        "to follow the growth of the pediatric patients over treatment",
        "(Results 'Population PK Modeling'). Phoenix covariate wt <- 'WT'",
        "in Supplemental Information 3."
      ),
      source_name = "WT"
    )
  )

  # Covariates screened by Lee 2022 but not retained in the final PK or PK-PD
  # model; no point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (continuous) or age population (adults / children / infants)",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested by stepwise forward inclusion / backward elimination on the",
        "PK parameters and not significant (Results 'Population PK",
        "Modeling'); significant on E0 and Emax in the PK-PD model only",
        "until time-varying WT was added, after which it was dropped",
        "(Results 'Population PK-PD Modeling')."
      )
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status",
      units = "(binary)",
      type = "binary",
      notes = "Tested stepwise on the PK parameters and not significant; screened graphically on the PK-PD parameters and not tested further (Lee 2022 Results)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened graphically on PK and PK-PD parameters (Supplemental Information 5e); no trend, not tested further."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race / ethnicity and country screened graphically (Supplemental Information 5d); no trend, not tested further."
    ),
    FGF23 = list(
      description = "Serum FGF23 (baseline intact, and time-varying total)",
      units = "pg/mL",
      type = "continuous",
      notes = paste(
        "Time-varying total FGF23 (free plus burosumab-bound) was tested",
        "stepwise on CL/F and V/F and was not significant; baseline FGF23",
        "did not appear to affect the PK-PD parameters (Lee 2022 Results).",
        "The Supplemental Information 3 PK-PD code still carries",
        "(fgf/356) and (fgf/103) power terms from an intermediate model;",
        "they are not part of the final model of Table 4."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "burosumab",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "burosumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 277,
    n_studies = 9,
    age_range = "1-12 years (pediatric, 94 subjects) and >17 years (adults, 183 subjects); adult baseline ages 18 to <69 years (Methods 'PK-PD Simulations in Age Groups Without Clinical Data')",
    age_groups = c(infants_1_2y = 6, children_2_12y = 88, adults_gt17y = 183),
    weight_range = "not tabulated; pediatric and adult baseline body weights were used to fit a GAMLSS weight model for the adolescent simulations",
    sex_female_pct = 59.2,
    race_ethnicity = c(White = 84.5, Black = 2.9, Asian = 9.7, Other = 2.9),
    disease_state = "X-linked hypophosphatemia (XLH); baseline serum phosphorus geometric mean 2.6 (infants), 2.4 (children) and 2.0 mg/dL (adults)",
    dose_range = "0.05-1.2 mg/kg subcutaneous, single dose or every 2 or 4 weeks, with intra-subject titration on serum phosphorus",
    regions = "North America and Europe (NCT00830674, NCT01340482, NCT01571596, NCT02163577, NCT02750618, NCT02312687, NCT02526160, NCT02915705, NCT02537431)",
    notes = paste(
      "Demographics from Lee 2022 Table 1 and baseline laboratory values",
      "from Table 2. 2844 measurable serum burosumab concentrations (PK) and",
      "6047 serum phosphorus observations (PD). Intravenous data from the",
      "first single-ascending-dose study were excluded. Sex and race",
      "percentages are pooled across the three age columns of Table 1",
      "(164 of 277 female; 234 White, 8 Black, 27 Asian, 8 Other).",
      "Anti-drug-antibody positive: 24 of 277 (8.7%)."
    )
  )

  ini({
    # ---- Structural PK, Lee 2022 Table 3 (typical values at 70 kg) ----
    lka <- log(0.395)
    label("First-order absorption rate constant (1/day)") # Table 3: Ka = 0.395 1/day (RSE 7.89%); no absorption lag (Results 'Population PK Modeling')
    lcl <- log(0.297)
    label("Apparent clearance CL/F at 70 kg (L/day)") # Table 3: CL/F = 0.297 L/day (RSE 8.41%)
    lvc <- log(9.02)
    label("Apparent volume of distribution V/F at 70 kg (L)") # Table 3: V/F = 9.02 L (RSE 13.9%)

    # ---- Allometric exponents (estimated, not fixed) ----
    e_wt_cl <- 0.912
    label("Power exponent of (WT/70) on CL/F (unitless)") # Table 3: covariate effect of WT on CL/F, (WT/70) exponent = 0.912 (RSE 6.39%)
    e_wt_vc <- 1.05
    label("Power exponent of (WT/70) on V/F (unitless)") # Table 3: covariate effect of WT on V/F, (WT/70) exponent = 1.05 (RSE 11.8%)

    # ---- PK between-subject variability ----
    # Table 3 'BSV (Shrinkage, %)' column on log-normal parameters, converted
    # with omega^2 = log(CV^2 + 1). BSV on Ka is 0% (FIX), so no Ka eta. The
    # text states the BSV was modelled 'without correlation' (Results
    # 'Population PK Modeling') and no covariance is tabulated, so the
    # block(nV, nCl) starting values in the Supplemental Information 3 code
    # are not carried.
    etalcl ~ 0.121227 # Table 3: BSV CL/F = 35.9% (shrinkage 6.15%) -> log(0.359^2 + 1)
    etalvc ~ 0.113694 # Table 3: BSV V/F = 34.7% (shrinkage 24.3%) -> log(0.347^2 + 1)

    # ---- PK residual error ----
    # Supplemental Information 3: observe(CObs = C + CEps * (1 + C * CMixRatio)),
    # i.e. SD = CEps + (CEps * CMixRatio) * C, a linear sum of the additive and
    # proportional SDs (nlmixr2 combined1).
    addSd <- 32.2
    label("Additive residual error on serum burosumab (ng/mL)") # Table 3: additive error = 32.2 ng/mL (RSE 17.7%)
    propSd <- 0.207
    label("Proportional residual error on serum burosumab (fraction)") # Table 3: proportional error = 20.7% (RSE 25.5%)

    # ---- PK-PD sigmoid Emax model, Lee 2022 Table 4 (typical values at 70 kg) ----
    le0 <- log(2.03)
    label("Serum phosphorus with no drug present, E0, at 70 kg (mg/dL)") # Table 4: E0 = 2.03 mg/dL (RSE 1.09%)
    lec50 <- log(4131)
    label("Burosumab concentration at half-maximal effect, EC50 (ng/mL)") # Table 4: EC50 = 4131 ng/mL (RSE 13.3%)
    lemax <- log(1.59)
    label("Maximal increase in serum phosphorus, Emax, at 70 kg (mg/dL)") # Table 4: Emax = 1.59 mg/dL (RSE 5.49%)
    lhill <- fixed(log(1))
    label("Hill coefficient gamma (unitless)") # Table 4: Gamma = 1 (FIX), BSV 0 (FIX); estimate was 0.936 before fixing (Results 'Population PK-PD Modeling')
    e_wt_e0 <- -0.143
    label("Power exponent of (WT/70) on E0 (unitless)") # Table 4: covariate effect of WT on E0, (WT/70) exponent = -0.143 (RSE 8.83%)
    e_wt_emax <- 0.168
    label("Power exponent of (WT/70) on Emax (unitless)") # Table 4: covariate effect of WT on Emax, (WT/70) exponent = 0.168 (RSE 29.1%)

    # ---- PK-PD between-subject variability ----
    # Table 4 'BSV (Shrinkage, %)' column, omega^2 = log(CV^2 + 1); no
    # covariance is tabulated, so the block(nE0, nEmax, nEC50) starting values
    # in the Supplemental Information 3 code are not carried.
    etale0 ~ 0.014297 # Table 4: BSV E0 = 12.0% (shrinkage 11.7%) -> log(0.120^2 + 1)
    etalec50 ~ 0.989541 # Table 4: BSV EC50 = 130% (shrinkage 28.3%) -> log(1.30^2 + 1)
    etalemax ~ 0.249220 # Table 4: BSV Emax = 53.2% (shrinkage 36.6%) -> log(0.532^2 + 1)

    # ---- PK-PD residual error ----
    propSd_serum_phosphorus <- 0.132
    label("Proportional residual error on serum phosphorus (fraction)") # Table 4: proportional error = 13.2% (RSE 1.82%); Supplemental Information 3 observe(EObs = E * (1 + EEps))
  })

  model({
    # ---- 1. Individual PK parameters (Supplemental Information 3 PK code) ----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # ---- 2. One-compartment disposition with first-order SC absorption ----
    # Only subcutaneous data were modelled, so CL/F and V/F are apparent
    # parameters and no bioavailability term is estimated.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 3. Serum burosumab concentration ----
    # central is in mg and vc in L (mg/L = ug/mL); the assay, the additive
    # residual error and EC50 are all in ng/mL, hence the factor of 1000.
    Cc <- 1000 * central / vc

    # ---- 4. Direct sigmoid Emax model for serum phosphorus ----
    # Supplemental Information 3 PK-PD code:
    #   E = E0 + Emax * C^Gam / (EC50^Gam + C^Gam)
    # with E the absolute serum phosphorus (mg/dL) and C the burosumab
    # concentration. The authors fitted the PK-PD model sequentially on
    # individual PK predictions (Methods 'Population PK-PD Modeling'); here
    # the PD is driven by the model's own Cc. No effect compartment or
    # turnover delay: no hysteresis was observed (Discussion).
    e0 <- exp(le0 + etale0) * (WT / 70)^e_wt_e0
    ec50 <- exp(lec50 + etalec50)
    emax <- exp(lemax + etalemax) * (WT / 70)^e_wt_emax
    hill <- exp(lhill)

    serum_phosphorus <- e0 + emax * Cc^hill / (ec50^hill + Cc^hill)

    # ---- 5. Observations ----
    Cc ~ add(addSd) + prop(propSd) + combined1()
    serum_phosphorus ~ prop(propSd_serum_phosphorus)
  })
}
