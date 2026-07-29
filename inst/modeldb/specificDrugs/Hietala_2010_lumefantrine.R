# Population PK model of oral lumefantrine (LUM) in Tanzanian children with
# uncomplicated falciparum malaria treated with the artemether-lumefantrine
# combination Coartem (Hietala 2010, Antimicrob Agents Chemother 54:4780-4788;
# doi:10.1128/AAC.00252-10).

Hietala_2010_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine (LUM) in 50 Tanzanian children",
    "(ages 1-10 years, weights 8-30 kg) with uncomplicated Plasmodium",
    "falciparum malaria treated with the standard six-dose weight-based",
    "Coartem (artemether 20 mg + lumefantrine 120 mg per tablet) regimen at",
    "0, 8, 24, 36, 48, and 60 hours (Hietala 2010). One-compartment",
    "disposition with first-order absorption preceded by an absorption lag",
    "time. The paper tested co-administration with full-fat (3.4%) cow's",
    "milk as a categorical covariate on the PK parameters of LUM; the",
    "effect did not improve the model fit and is not encoded here",
    "(Discussion: 'the resulting number of doses actually administered",
    "with an adequate amount of milk may have been too small to allow the",
    "detection of a difference'). All PK parameters are reported per kg",
    "body weight (linear weight normalisation applied inside model()).",
    sep = " "
  )
  reference <- paste(
    "Friberg Hietala S, Martensson A, Ngasala B, Dahlstrom S, Lindegardh N,",
    "Annerberg A, Premji Z, Farnert A, Gil P, Bjorkman A, Ashton M (2010).",
    "Population pharmacokinetics and pharmacodynamics of artemether and",
    "lumefantrine during combination treatment in children with uncomplicated",
    "falciparum malaria in Tanzania. Antimicrob Agents Chemother",
    "54(11):4780-4788. doi:10.1128/AAC.00252-10.",
    sep = " "
  )
  vignette <- "Hietala_2010_artemether_lumefantrine_malaria"
  units <- list(time = "hour", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (per-kg) weight normalisation of CL / V parameters.",
        "Hietala 2010 reports the lumefantrine PK parameters on a",
        "per-kg basis (Table 2: 'CL/F (ml/h/kg)', 'V/F (liters/kg)'),",
        "so individual CL = lcl * WT, Vc = lvc * WT inside model().",
        "Population mean WT 14 kg (range 8-30 kg, Results)."
      ),
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    MILK = list(
      description        = "Concomitant intake of 200 mL full-fat (3.4%) cow's milk with each dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (water)",
      notes              = paste(
        "Hietala 2010 randomised the 50 patients 1:1 between dosing with",
        "200 mL of full-fat (3.4%) cow's milk and dosing with water (n = 25",
        "per arm). The protocol-mandated 200 mL of milk was completed on",
        "only 43% of intended dose-occasions (Results, 'Demographics and",
        "safety'). The paper tested MILK as a categorical covariate on the",
        "LUM PK parameters; it 'did not explain the variability in LUM",
        "pharmacokinetics and did not result in an improvement of the model'",
        "(Results, 'Pharmacokinetics of artemether, dihydroartemisinin, and",
        "lumefantrine') and was not retained in the final model. Documented",
        "here for provenance but excluded from model()."
      ),
      source_name        = "MILK"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 50L,
    n_studies       = 1L,
    n_observations  = 423L,
    age_range       = "1-10 years (mean 4; Results)",
    weight_range    = "8-30 kg (mean 14; Results)",
    sex_female_pct  = 62,
    disease_state   = paste(
      "Acute uncomplicated Plasmodium falciparum malaria with asexual",
      "parasite density 2,000-200,000 / microL at admission and either",
      "axillary temperature >= 37.5 degC or history of fever within 24 h",
      "(Methods)."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet. Weight-based dosing: 5-14 kg -> 1 tablet/dose,",
      "15-24 kg -> 2 tablets/dose, 25-34 kg -> 3 tablets/dose. Six",
      "doses (oral) at 0, 8, 24, 36, 48, and 60 hours."
    ),
    regions            = "Tanzania (Fukayosi Primary Health Care Centre, Bagamoyo District)",
    trial_registration = "ClinicalTrials.gov NCT00336375",
    notes              = paste(
      "Demographics from Hietala 2010 Results. 423 LUM plasma",
      "concentrations from 50 patients were included in the PK model",
      "(19 LUM samples below limit-of-detection were excluded). Companion",
      "artemether-DHA PK model from the same cohort:",
      "modellib('Hietala_2010_artemether'). Companion parasitemia PD model:",
      "modellib('Hietala_2010_artemether_parasitemia')."
    )
  )

  ini({
    # Structural PK parameters from Hietala 2010 Table 2 ("Estimate (95% CI)" column).
    # The lumefantrine clearance is reported in mL/h/kg (note the milli-litre
    # unit) whereas the volume of distribution is in L/kg; both are
    # normalised by WT inside model(). Absorption lag time tlag and the
    # first-order absorption rate constant ka are individual-scale rates
    # (1/h or h) and are NOT per-kg.
    ltlag <- log(1.92)        ; label("Absorption lag time tlag (h)")                                 # Hietala 2010 Table 2: Lag = 1.92 (95% CI 1.86-1.96)
    lka   <- log(0.82)        ; label("Absorption rate constant ka (1/h)")                            # Hietala 2010 Table 2: ka = 0.82 (95% CI 0.45-1.61)
    lcl   <- log(0.077)       ; label("Apparent oral lumefantrine clearance CL/F (L/h/kg)")           # Hietala 2010 Table 2: CL/F = 77 mL/h/kg = 0.077 L/h/kg (95% CI 52-105 mL/h/kg)
    lvc   <- log(8.9)         ; label("Apparent lumefantrine volume of distribution V/F (L/kg)")      # Hietala 2010 Table 2: V/F = 8.9 (95% CI 6.8-11.7)

    # Inter-individual variability (IIV). Table 2 reports IIV as %CV per the
    # NONMEM log-normal convention; the corresponding internal variance is
    #   omega^2 = log((CV/100)^2 + 1)
    #     ka  CV 156% -> log(1.56^2 + 1) = 1.241237
    #     V/F CV  82% -> log(0.82^2 + 1) = 0.515660
    # IIV on Lag and CL/F was not estimated (Table 2 "NE") and is therefore
    # omitted here.
    etalka ~ 1.241237          # Hietala 2010 Table 2: IIV on ka  = 156% CV (95% CI 126-190)
    etalvc ~ 0.515660          # Hietala 2010 Table 2: IIV on V/F =  82% CV (95% CI  66-102)

    # Residual error. Hietala 2010 reports separate proportional and additive
    # residual errors for LUM (Table 2). The additive component has units of
    # concentration (nM here, matching the SPE / LC / UV quantification used
    # for LUM). The paper provides 95% CIs for both, indicating they were
    # estimated (no "fixed" annotation in Table 2 caption).
    propSd <- 0.46            ; label("Proportional residual SD for lumefantrine plasma concentration (fraction)") # Hietala 2010 Table 2: sigma_prop = 46% (95% CI 41-52)
    addSd  <- 43              ; label("Additive residual SD for lumefantrine plasma concentration (nM)")           # Hietala 2010 Table 2: sigma_add  = 43 nM (95% CI 19-66)
  })

  model({
    # Molecular weight of lumefantrine (g/mol). Used to convert plasma
    # concentration from mg/L (= ug/mL) to nM. Source: lumefantrine free-base
    # MW = 528.94 g/mol (Coartem package insert / DrugBank DB01219).
    mw_lum <- 528.94

    # Individual PK parameters. Linear (per-kg) weight normalisation of CL
    # and V; ka and tlag are individual-scale (not per-kg).
    tlag <- exp(ltlag)
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl)            * WT
    vc   <- exp(lvc + etalvc)   * WT

    kel <- cl / vc

    # One-compartment disposition with first-order absorption and absorption
    # lag time on the depot compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    # Plasma concentration in nM. Dose units mg; vc units L; central / vc has
    # units mg/L. Multiplication by 1e6 / MW converts mg/L to nmol/L = nM.
    Cc <- 1e6 * central / vc / mw_lum

    # Combined proportional + additive residual error on the linear-nM scale.
    Cc ~ add(addSd) + prop(propSd)
  })
}
