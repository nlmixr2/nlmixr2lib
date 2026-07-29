# One-compartment subcutaneous population PK model for anakinra
# (recombinant human IL-1 receptor antagonist) in 87 children and
# adolescents (8 months to 21 years; 4.3 to 83 kg) with systemic-onset
# juvenile idiopathic arthritis (SJIA) and autoinflammatory syndromes
# (Urien 2013, BMC Pharmacology and Toxicology 14:40;
# doi:10.1186/2050-6511-14-40).

Urien_2013_anakinra <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model for subcutaneous",
    "anakinra (recombinant nonglycosylated human IL-1 receptor antagonist)",
    "in 87 children and adolescents (8 months to 21 years, 4.3 to 83 kg)",
    "treated for systemic-onset juvenile idiopathic arthritis (SJIA) and",
    "diverse autoinflammatory syndromes (Urien 2013). First-order",
    "absorption (Ka) into a single central compartment with first-order",
    "elimination; apparent clearance CL/F and apparent volume V/F are",
    "allometrically scaled to body weight with estimated power exponents",
    "(0.47 on CL/F and 0.76 on V/F, reference 70 kg). Inter-individual",
    "variability is reported on CL/F and between-occasion variability on",
    "V/F; no other covariate effect (age, sex, co-administered",
    "anti-inflammatory drugs) was retained.",
    sep = " "
  )
  reference <- paste(
    "Urien S, Bardin C, Bader-Meunier B, Mouy R, Compeyrot-Lacassagne S,",
    "Foissac F, Florkin B, Wouters C, Neven B, Treluyer J-M, Quartier P",
    "(2013). Anakinra pharmacokinetics in children and adolescents with",
    "systemic-onset juvenile idiopathic arthritis and autoinflammatory",
    "syndromes. BMC Pharmacology and Toxicology 14:40.",
    "doi:10.1186/2050-6511-14-40.",
    sep = " "
  )
  vignette <- "Urien_2013_anakinra"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Urien 2013. Used for allometric",
        "scaling of CL/F and V/F with estimated power exponents 0.47",
        "and 0.76 respectively, reference 70 kg (Urien 2013 Table 1",
        "and prose: CL/F = 0.847 * BW^0.47, V/F = 2.581 * BW^0.76,",
        "equivalent to the 70 kg-normalized form in Table 1)."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 87L,
    n_observations  = 148L,
    n_studies       = 2L,
    age_range       = "0.73-21 years (8 months to 21 years; SJIA cohort 2.26-16.8 y, autoinflammatory cohort 0.73-21 y)",
    age_median      = "SJIA cohort 7.6 years; autoinflammatory cohort 8 years (overall median not reported)",
    weight_range    = "4.3-83 kg (SJIA cohort 10-83 kg, autoinflammatory cohort 4.3-60 kg)",
    weight_median   = "21 kg (both cohorts; Results)",
    sex_female_pct  = 36.8,
    race_ethnicity  = NA,
    disease_state   = paste(
      "22 patients with systemic-onset juvenile idiopathic arthritis",
      "(SJIA, from the ANAJIS phase IIB trial) plus 65 patients with",
      "diverse autoinflammatory conditions (20 cryopyrin-associated",
      "periodic syndromes / CAPS, including 14 CINCA/NOMID and 6",
      "Muckle-Wells; 3 mevalonate kinase deficiency; 2 TNF-receptor",
      "associated periodic syndromes; 1 familial mediterranean fever;",
      "remaining genetically undetermined autoinflammatory conditions)."
    ),
    dose_range      = paste(
      "Subcutaneous anakinra once daily. ANAJIS-trial SJIA patients",
      "received 2 mg/kg/day (maximum 100 mg). Autoinflammatory patients",
      "received 2-10 mg/kg/day, with the highest doses in low-weight",
      "CAPS patients who had failed lower doses."
    ),
    regions         = "France (multicentre, Necker-Cochin Inserm CIC and partner sites)",
    notes           = paste(
      "Combined PK dataset from the ANAJIS phase IIB SJIA trial (Quartier",
      "2011 Ann Rheum Dis 70:747-754) and patients with diverse",
      "autoinflammatory conditions subsequently treated at the same",
      "centre. 148 plasma anakinra concentrations were available for",
      "modelling; four below the limit of quantification (40 ng/mL)",
      "were treated as left-censored. Demographics and dose distribution",
      "from the Results, 'Population pharmacokinetic modeling' paragraph.",
      "The paper reports the cohort as 32 girls / 52 boys totalling 87",
      "(32+52 = 84, an apparent typo in the source); sex_female_pct is",
      "computed as 32/87 = 36.8 %."
    )
  )

  ini({
    # Structural PK parameters -- typical 70 kg patient. All values from
    # Urien 2013 Table 1 (estimates with low RSEs).
    lka     <- log(0.38);   label("First-order SC absorption rate Ka (1/h)")                                 # Urien 2013 Table 1: Ka = 0.38 1/h (RSE 19%)
    lcl     <- log(6.24);   label("Apparent clearance CL/F at 70 kg reference (L/h)")                        # Urien 2013 Table 1: CL/F = 6.24 L/h/70 kg (RSE 8%)
    lvc     <- log(65.2);   label("Apparent volume of distribution V/F at 70 kg reference (L)")              # Urien 2013 Table 1: V/F = 65.2 L/70 kg (RSE 12%)

    # Allometric exponents on body weight, estimated (not fixed at the
    # canonical 0.75/1.0). Urien 2013 Discussion explicitly retains the
    # estimated exponents because anakinra is a 153 amino acid peptide
    # whose values diverged from the theoretical allometric defaults.
    e_wt_cl <- 0.47;        label("Allometric exponent on CL/F (unitless)")                                  # Urien 2013 Table 1: beta_CL = 0.47 (RSE 14%)
    e_wt_vc <- 0.76;        label("Allometric exponent on V/F (unitless)")                                   # Urien 2013 Table 1: beta_V = 0.76 (RSE 16%)

    # Inter-individual variability. Monolix v3.2 reports omega (SD of the
    # log-normal random effect) in its "Parameter estimates" column;
    # converted here to the variance scale required by nlmixr2 ini(). The
    # paper places IIV (eta) on CL/F and between-occasion variability
    # (gamma) on V/F only; no eta on Ka or V/F and no gamma on CL/F or Ka
    # are reported (Results, Population pharmacokinetic modeling).
    etalcl ~ 0.0784         # Urien 2013 Table 1: eta_CL/F (omega) = 0.28 (RSE 15 %); variance = 0.28^2 = 0.0784

    # Between-occasion variability on V/F (gamma_V/F = 0.47 in Urien
    # 2013 Table 1, omega/SD scale) is encoded here as an additional
    # eta on V/F so that the simulated population reflects the V/F
    # spread reported in the paper; see vignette Assumptions and
    # deviations for the BOV-as-BSV simulation rationale.
    etalvc ~ 0.2209         # Urien 2013 Table 1: gamma_V/F (omega) = 0.47 (RSE 17 %); variance = 0.47^2 = 0.2209 (encoded as BSV for simulation)

    # Residual error on plasma anakinra concentrations -- additive (no
    # proportional component retained in the final model per Results).
    addSd  <- 0.072;        label("Additive residual error on plasma anakinra (mg/L)")                       # Urien 2013 Table 1: epsilon = 0.072 mg/L (RSE 10 %)
  })

  model({
    # Allometric body-weight scaling (reference 70 kg). Estimated
    # exponents per Table 1 (final model retained the paper-fitted
    # values rather than the canonical 3/4 and 1).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka)

    # Disposition. SC dose enters the depot; first-order absorption into
    # the central compartment; first-order elimination from central.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # Observation: plasma anakinra concentration in mg/L (matches Urien
    # 2013 Figure 1 axis units and the residual-error units in Table 1).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
