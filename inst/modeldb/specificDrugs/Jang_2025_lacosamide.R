Jang_2025_lacosamide <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral lacosamide in Korean adults with epilepsy, fit jointly to plasma/serum and saliva concentrations: saliva is the plasma concentration times an estimated saliva:plasma scale factor of 1.44; CYP2C19 poor-metabolizer phenotype, concomitant carbamazepine and phenobarbital, and creatinine clearance on apparent clearance; body surface area on apparent volume"
  reference <- paste(
    "Jang Y, Ahn SJ, Lee S, Kim Y, Lee HS, Lee YK, Shin YW, Lee ST, Jung KH,",
    "Park KI, Jung KY, Cho JY, Lee S, Yu KS, Chu K, Jang IJ, Lee SK.",
    "Saliva-based lacosamide monitoring paves the way toward personalized",
    "epilepsy pharmacotherapy. Sci Rep. 2025;15:19181.",
    "doi:10.1038/s41598-025-04044-x. PMCID: PMC12126552.",
    sep = " "
  )
  vignette <- "Jang_2025_lacosamide"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (estimating equation not stated by the source)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw creatinine clearance in mL/min, NOT BSA-normalized: Supplementary",
        "Table 2 labels the exponent 'Exponent for CLcr (normalized by 107",
        "mL/min)' and the Results print the term as (CLcr/107)^0.579. Jang 2025",
        "does not state how CLcr was computed; the upstream healthy-volunteer",
        "model it starts from (Kim 2012, Pharmacology 89:172) used creatinine",
        "clearance on CL/F, and Cockcroft-Gault is the usual choice in that",
        "group's work, so supply a Cockcroft-Gault value in mL/min. The 107",
        "mL/min reference is the value printed in the covariate equation; the",
        "cohort's CLcr distribution is not tabulated (Table 1 gives only serum",
        "creatinine 0.8 +/- 0.1 mg/dL and a BSA-normalized GFR of 101 +/- 19",
        "mL/min/1.73 m^2).",
        sep = " "
      ),
      source_name = "CLcr"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power-scaled on apparent volume of distribution, referenced to 1.856",
        "m^2 with a FIXED exponent of 1.37 (Supplementary Table 2, theta6 and",
        "theta7; both fixed from the healthy-volunteer model of Kim 2012). The",
        "BSA formula is not stated (DuBois / Mosteller / Haycock unspecified).",
        "Cohort BSA 1.8 +/- 0.2 m^2 (Table 1).",
        sep = " "
      ),
      source_name = "BSA"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP2C19 normal or intermediate metabolizer)",
      notes = paste(
        "1 = CYP2C19 poor metabolizer (*2/*2, *2/*3 or *3/*3; Methods 'Genetic",
        "analysis'); 0 = normal (*1/*1) or intermediate (*1/*2, *1/*3)",
        "metabolizer. The paper estimates one typical CL/F for normal and",
        "intermediate metabolizers pooled (theta1 = 1.78 L/h) and a separate one",
        "for poor metabolizers (theta2 = 1.43 L/h; Supplementary Table 2), so",
        "CYP2C19_IM is deliberately not used. Time-fixed (germline). Only 37",
        "patients were genotyped: 13 normal, 18 intermediate and 6 poor",
        "metabolizers (Results 'Saliva effectively represented ...'); how the",
        "model treated the ungenotyped patients is not stated.",
        sep = " "
      ),
      source_name = "CYP2C19 phenotype"
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes = paste(
        "1 = taking carbamazepine at the time of sampling. Enters CL/F",
        "additively with phenobarbital inside one linear factor,",
        "(1 + 0.771 * CBZ + 0.467 * PB) (Results; Supplementary Table 2 theta3).",
        "The number of carbamazepine users is not reported; 63.1 percent of the",
        "cohort were on any concomitant antiseizure medication (Table 1).",
        sep = " "
      ),
      source_name = "CBZ"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes = paste(
        "1 = taking phenobarbital at the time of sampling. Enters CL/F",
        "additively with carbamazepine inside (1 + 0.771 * CBZ + 0.467 * PB)",
        "(Results; Supplementary Table 2 theta4). The paper does not say whether",
        "primidone was pooled with phenobarbital, and does not report the number",
        "of phenobarbital users.",
        sep = " "
      ),
      source_name = "PB"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lacosamide", units = "mg", specimen = "administration site", verified = TRUE),
    # Saliva is not a compartment: Csaliva is the central (plasma)
    # concentration multiplied by fsaliva (Figure 2A).
    central = list(analyte = "lacosamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 123,
    n_observations = 270,
    n_studies = 1,
    age_range = "41.3 +/- 14.2 years (mean +/- SD)",
    weight_range = "67.2 +/- 14.3 kg (mean +/- SD)",
    bsa_range = "1.8 +/- 0.2 m^2 (mean +/- SD)",
    sex_female_pct = 45.5,
    race_ethnicity = c(Asian = 100),
    disease_state = "epilepsy (91.9 percent focal, 8.1 percent generalized seizures) on stable lacosamide maintenance for at least one month",
    dose_range = "oral lacosamide twice daily; total daily dose median 300 mg (IQR 200-400)",
    renal_function = "serum creatinine 0.8 +/- 0.1 mg/dL; GFR 101 +/- 19 mL/min/1.73 m^2 (Table 1)",
    co_medication = "77 of 123 (63.1 percent) on concomitant antiseizure medications; median 3 antiseizure medications [IQR 2-5]",
    genotype_distribution = c(CYP2C19_NM = 13, CYP2C19_IM = 18, CYP2C19_PM = 6),
    regions = "Republic of Korea (single centre, Seoul National University Hospital)",
    notes = paste(
      "Prospective outpatient study, August 2021 to November 2022 (Methods;",
      "Table 1). The PK dataset holds 109 spot saliva samples, 59 home",
      "self-collected follow-up saliva samples and 102 blood (plasma or serum)",
      "samples from 86 patients, i.e. 270 concentrations (Results 'Population",
      "pharmacokinetic model of lacosamide'). Sparse sampling (median 4.25 h",
      "after the last dose for saliva, 4.67 h for blood) is why ka, V/F and",
      "the BSA exponent on V/F were fixed to the healthy-volunteer model of",
      "Kim 2012. Race recorded as 100 percent Asian because the cohort is a",
      "Korean single-centre population; the paper gives no race breakdown.",
      "genotype_distribution counts the 37 genotyped patients only.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Supplementary Table 2. ka, V/F and the BSA
    # exponent were FIXED to the healthy Korean volunteer model of Kim 2012
    # (Methods 'Pharmacokinetic analysis and modeling': "we fixed the
    # absorption rate (Ka), Vd/F, and the power coefficient of BSA to Vd/F
    # based on values from the published model"); the table prints no RSE
    # for them.
    lka <- fixed(log(6.47)); label("Absorption rate constant ka (1/h)")                                                      # Supp Table 2: theta8 = 6.47 h-1, RSE '-' (fixed)
    lcl <- log(1.78);        label("Apparent clearance CL/F, CYP2C19 normal or intermediate metabolizer, CLcr 107 mL/min (L/h)") # Supp Table 2: theta1 = 1.78 L/h, RSE 3.0% (bootstrap 1.78, 1.68-1.89)
    lvc <- fixed(log(41.6)); label("Apparent volume of distribution V/F at BSA 1.856 m^2 (L)")                                # Supp Table 2: theta6 = 41.6 L, RSE '-' (fixed)

    # Saliva:plasma concentration scale factor, Csaliva = fsaliva * Cc
    # (Figure 2A; Results "The scaling factor of 1.44 well explained the
    # saliva lacosamide concentration").
    lfsaliva <- log(1.44); label("Saliva:plasma concentration scale factor (unitless)")                                       # Supp Table 2: theta9 = 1.44, RSE 4.0% (bootstrap 1.44, 1.33-1.55)

    # Covariate effects. Supplementary Table 2 structure row:
    # CL/F = (theta1 or theta2) x (1 + theta3 x CBZ + theta4 x PB) x (CLcr/107)^theta5
    # V/F  = theta6 x (BSA/1.856)^theta7
    # theta2 is a separately estimated typical CL/F for CYP2C19 poor
    # metabolizers; it is carried here as the ratio theta2/theta1 so the
    # pooled normal + intermediate value stays in lcl. Both published
    # values are written out literally.
    e_cyp2c19_pm_cl <- 1.43 / 1.78; label("Multiplicative factor on CL/F for CYP2C19 poor vs normal/intermediate metabolizers (unitless)") # Supp Table 2: theta2 = 1.43 L/h (RSE 7.4%, bootstrap 1.43, 1.28-1.59) divided by theta1 = 1.78 L/h
    e_conmed_cbz_cl <- 0.771;       label("Fractional increase in CL/F with concomitant carbamazepine (unitless)")                       # Supp Table 2: theta3 = 0.771, RSE 35.4% (bootstrap 0.794, 0.333-1.77)
    e_conmed_pb_cl  <- 0.467;       label("Fractional increase in CL/F with concomitant phenobarbital (unitless)")                       # Supp Table 2: theta4 = 0.467, RSE 45.6% (bootstrap 0.477, 0.112-0.841)
    e_crcl_cl       <- 0.579;       label("Power exponent of (CRCL/107 mL/min) on CL/F (unitless)")                                     # Supp Table 2: theta5 = 0.579, RSE 13.7% (bootstrap 0.576, 0.399-0.771)
    e_bsa_vc        <- fixed(1.37); label("Power exponent of (BSA/1.856 m^2) on V/F (unitless)")                                         # Supp Table 2: theta7 = 1.37, RSE '-' (fixed)

    # Inter-individual variability on CL/F only, reported as 23.8 %CV.
    # Log-normal variance omega^2 = log(1 + 0.238^2) = 0.05510.
    etalcl ~ 0.05510  # Supp Table 2: omega CL/F = 23.8 %CV, RSE 10.3% (bootstrap 23.2, 15.6-29.4)

    # Residual error: separate proportional errors for blood and saliva.
    # The tabulated values are SDs (fractions), not variances: the
    # bootstrap 95% CI relative widths (0.28 blood, 0.22 saliva) match the
    # SD prediction 3.92*sqrt(1/(2N)) (0.27 for N = 102, 0.21 for N = 168)
    # and are half the variance prediction, and the RSEs (7.3%, 6.5%) match
    # 1/sqrt(2N) for an SD.
    propSd         <- 0.270; label("Proportional residual SD for plasma/serum Cc (fraction)") # Supp Table 2: sigma blood,prop = 0.270, RSE 7.3% (bootstrap 0.269, 0.230-0.305)
    propSd_Csaliva <- 0.407; label("Proportional residual SD for saliva Csaliva (fraction)")  # Supp Table 2: sigma saliva,prop = 0.407, RSE 6.5% (bootstrap 0.405, 0.363-0.453)
  })

  model({
    # Individual parameters (Supplementary Table 2 structure rows).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      e_cyp2c19_pm_cl^CYP2C19_PM *
      (1 + e_conmed_cbz_cl * CONMED_CBZ + e_conmed_pb_cl * CONMED_PB) *
      (CRCL / 107)^e_crcl_cl
    vc <- exp(lvc) * (BSA / 1.856)^e_bsa_vc
    fsaliva <- exp(lfsaliva)

    kel <- cl / vc

    # One compartment, first-order absorption and linear elimination
    # (Figure 2A). No bioavailability term: F is absorbed into CL/F and V/F.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma/serum concentration and the saliva concentration scaled from it.
    Cc <- central / vc
    Csaliva <- fsaliva * Cc

    Cc ~ prop(propSd)
    Csaliva ~ prop(propSd_Csaliva)
  })
}
