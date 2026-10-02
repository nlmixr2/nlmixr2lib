Ali_2022_artemisinin <- function() {
  description <- paste(
    "One-compartment population PK model for oral artemisinin given as a",
    "single dose of the fixed-dose artemisinin-naphthoquine combination to",
    "Tanzanian children (6 years and older) and adults with uncomplicated",
    "Plasmodium falciparum malaria (Ali 2022). Savic transit-compartment",
    "absorption (non-integer NN transit compartments and a separate",
    "first-order absorption rate ka into the central compartment), relative",
    "bioavailability fixed to 1 with between-subject variability, and",
    "allometric scaling by total body weight to a 55 kg reference (exponent",
    "0.75 on CL/F, 1 on V/F). Naphthoquine from the same combination is a",
    "separately fitted model (Ali_2022_naphthoquine)."
  )
  reference <- paste(
    "Ali AM, Gausi K, Jongo SA, Kassim KR, Mkindi C, Simon B, Mtoro AT,",
    "Juma OA, Lweno ON, Gwandu CH, Bakari BM, Mbaga TA, Milando FA, Hamad A,",
    "Shekalaghe SA, Abdulla S, Denti P, Penny MA (2022).",
    "Population Pharmacokinetics of Antimalarial Naphthoquine in Combination",
    "with Artemisinin in Tanzanian Children and Adults: Dose Optimization.",
    "Antimicrobial Agents and Chemotherapy 66(5):e01696-21.",
    "doi:10.1128/aac.01696-21."
  )
  vignette <- "Ali_2022_artemisinin_naphthoquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "artemisinin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "artemisinin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size descriptor on CL/F (exponent 0.75) and V/F (exponent 1), both fixed to the theory-based values (Methods 'Data analysis and pharmacokinetic modeling'), normalised to 55 kg (Table 2 footnote a: 'All clearances and volumes of distribution refer to a patient weighing 55 kg'). Fat-free mass was tested as an alternative descriptor for artemisinin and did not improve the fit (Results 'Pharmacokinetic modeling').",
      source_name = "wt"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = "Screened by stepwise covariate modeling (Methods) and not retained; 'No other available covariate was significant' (Results)."
    ),
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened by stepwise covariate modeling (Methods) and not retained."
    ),
    HGB = list(
      description = "Baseline hemoglobin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened by stepwise covariate modeling (Methods) and not retained."
    ),
    FFM = list(
      description = "Fat-free mass (Janmahasatian)",
      units = "kg",
      type = "continuous",
      notes = "Tested as an alternative to total body weight for artemisinin allometric scaling and did not improve the fit (Results 'Pharmacokinetic modeling')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28L,
    n_studies = 1L,
    age_range = "6.0-56.0 years (enrolled cohort of 29)",
    age_median = "13.1 years (enrolled cohort of 29)",
    weight_range = "20-84 kg (enrolled cohort of 29)",
    weight_median = "32.0 kg (enrolled cohort of 29)",
    sex_female_pct = 48.3,
    race_ethnicity = "Tanzanian (sub-Saharan African)",
    disease_state = "Uncomplicated Plasmodium falciparum malaria",
    dose_range = "Single oral dose of artemisinin-naphthoquine tablets (125 mg artemisinin + 50 mg naphthoquine per tablet), about 20 mg/kg artemisinin (1,000 mg for adults); achieved median 18.5 (IQR 16.7-18.8) mg/kg",
    regions = "Bagamoyo District, Tanzania",
    notes = "Phase IV single-centre randomised study (NCT01930331) conducted in 2014; 29 patients enrolled in the ART-NQ arm in three age groups (6-10 years n = 12, 11-17 years n = 6, >= 18 years n = 11; Table 1). One slow absorber was excluded from the artemisinin analysis, leaving 28 patients and 174 samples drawn at 1, 2, 4, 8, 12 and 18 h after the dose (none below the LLOQ of 1 ng/mL). Demographic summaries above are for the 29 enrolled patients (Table 1); the paper does not report them for the 28-patient artemisinin subset."
  )

  ini({
    # Table 2 (Artemisinin column). All clearances and volumes refer to a
    # patient weighing 55 kg (Table 2 footnote a).
    lcl <- log(66.7)
    label("Apparent clearance CL/F for a 55 kg patient (L/h)") # Table 2 'CL (L/h)' = 66.7 (57.3-78.5)
    lvc <- log(395)
    label("Apparent central volume V/F for a 55 kg patient (L)") # Table 2 'V1 (L)' = 395 (339-446)
    lka <- log(2.11)
    label("First-order absorption rate constant from the absorption compartment (1/h)") # Table 2 'Ka (1/h)' = 2.11 (1.22-3.18)
    lmtt <- log(0.987)
    label("Absorption mean transit time (h)") # Table 2 'MTT (h)' = 0.987 (0.72-1.31)
    lnn <- log(7.53)
    label("Number of absorption transit compartments (unitless, non-integer)") # Table 2 'NN' = 7.53 (5.10-13.7)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability (unitless)") # Table 2 'F' = 1.00 fixed; Methods 'Relative bioavailability (F) was fixed to 1'

    # Allometric exponents fixed to the theory-based values (Methods 'using the
    # suggested exponents of 1 for volumes of distribution and 3/4 for
    # clearance terms'; Table 2 footnote a equation).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)") # Table 2 footnote a: 'CL/F = theta_pop * (wt/55)^0.75 for artemisinin'
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V/F (unitless)") # Table 2 footnote a: 'V/F = theta_pop * (wt/55) for naphthoquine and artemisinin'

    # IIV reported as approximate %CV = sqrt(omega^2) * 100 (Table 2 footnote
    # d), so omega^2 = (CV/100)^2.
    etalcl ~ 0.034596 # Table 2 IIV 'CL' = 18.6 %CV -> 0.186^2
    etalka ~ 0.208849 # Table 2 IIV 'Ka' = 45.7 %CV -> 0.457^2
    etalmtt ~ 0.242064 # Table 2 IIV 'MTT' = 49.2 %CV -> 0.492^2
    etalfdepot ~ 0.168921 # Table 2 IIV 'F' = 41.1 %CV -> 0.411^2

    # Combined additive and proportional residual error (Methods).
    addSd <- fixed(0.20)
    label("Additive residual error (ng/mL)") # Table 2 'Additive error (ng/mL)' = 0.20 fixed; footnote c: 20% of the 1 ng/mL LLOQ
    propSd <- 0.307
    label("Proportional residual error (fraction)") # Table 2 'Proportional error (%)' = 30.7 (26.2-34.7)
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 55)^e_wt_cl
    vc <- exp(lvc) * (WT / 55)^e_wt_vc
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc

    # Savic 2007 transit-compartment absorption (Methods reference 52): the
    # dose passes through nn transit compartments (ktr = (nn + 1) / mtt inside
    # transit()) into the absorption compartment, which empties into central
    # at the separately estimated rate ka (Figure 1).
    d / dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d / dt(central) <- ka * depot - kel * central
    f(depot) <- 0

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
