Ali_2022_naphthoquine <- function() {
  description <- paste(
    "Two-compartment population PK model for oral naphthoquine given as a",
    "single dose of the fixed-dose artemisinin-naphthoquine combination to",
    "Tanzanian children (6 years and older) and adults with uncomplicated",
    "Plasmodium falciparum malaria (Ali 2022). Savic transit-compartment",
    "absorption (non-integer NN transit compartments and a separate",
    "first-order absorption rate ka into the central compartment), relative",
    "bioavailability fixed to 1 with between-subject variability, and",
    "allometric scaling of CL/F by fat-free mass to a 45 kg reference",
    "(exponent 0.75) and of Q/F (exponent 0.75), Vc/F and Vp/F (exponent 1)",
    "by total body weight to a 55 kg reference. Artemisinin from the same",
    "combination is a separately fitted model (Ali_2022_artemisinin)."
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
    depot = list(analyte = "naphthoquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "naphthoquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "naphthoquine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size descriptor on Q/F (exponent 0.75), Vc/F and Vp/F (exponent 1), normalised to 55 kg (Table 2 footnote a). Body weight was slightly better than fat-free mass for the volume parameters (dOFV -10.0 versus -8.9; Results 'Pharmacokinetic modeling').",
      source_name = "wt"
    ),
    FFM = list(
      description = "Fat-free mass at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size descriptor on CL/F (exponent 0.75) normalised to 45 kg, per the Table 2 footnote a equation 'CL/F = theta_pop * (FFM/45)^0.75 for naphthoquine'. The Abstract and Results instead quote the typical clearance for an individual with 44.3 kg fat-free mass; 45 kg is used because it is the printed equation and because it reproduces the paper's own 70 kg rescaling (CL = 52.0 L/h at FFM 56.1 kg; 45 kg gives 52.1 L/h, 44.3 kg gives 52.8 L/h). FFM was 'derived for males and females separately' (Methods, reference 55); the formula is not printed, so downstream users should supply FFM from the Janmahasatian 2005 equations. FFM beat total body weight for clearance (dOFV -14.4 versus -23.2 as printed; Results).",
      source_name = "FFM"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = "Tested and not retained (Results: 'Adding other available covariates (sex, age, fever, hemoglobin, temperature, and hematocrit) did not improve the model fit')."
    ),
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not retained (Results). Sex enters only through the FFM derivation."
    ),
    HGB = list(
      description = "Baseline hemoglobin",
      units = "g/dL",
      type = "continuous",
      notes = "Tested and not retained (Results); a hemoglobin effect on V1 reported in an earlier Papua New Guinea study was not found here (Discussion)."
    ),
    HCT = list(
      description = "Baseline hematocrit",
      units = "%",
      type = "continuous",
      notes = "Tested and not retained (Results)."
    ),
    BODYTEMP = list(
      description = "Body temperature / fever at baseline",
      units = "degC",
      type = "continuous",
      notes = "Tested (as fever and temperature) and not retained (Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 1L,
    age_range = "6.0-56.0 years",
    age_median = "13.1 years",
    weight_range = "20-84 kg",
    weight_median = "32.0 kg",
    sex_female_pct = 48.3,
    race_ethnicity = "Tanzanian (sub-Saharan African)",
    disease_state = "Uncomplicated Plasmodium falciparum malaria",
    dose_range = "Single oral dose of artemisinin-naphthoquine tablets (125 mg artemisinin + 50 mg naphthoquine per tablet), about 8 mg/kg naphthoquine (400 mg for adults); achieved median 7.4 (IQR 6.7-7.5) mg/kg",
    regions = "Bagamoyo District, Tanzania",
    notes = "Phase IV single-centre randomised study (NCT01930331) conducted in 2014; 29 patients in three age groups (6-10 years n = 12, 11-17 years n = 6, >= 18 years n = 11; Table 1). 363 naphthoquine concentrations (median 13 per patient, range 9-13) drawn at 1, 2, 4, 8, 12 and 18 h and on days 4, 7, 14, 21, 28 and 42; 5 (1.4%) below the LLOQ of 0.2 ng/mL, handled with Beal's M6 method. Typical patient for the reported parameters: 55 kg body weight, 45 kg fat-free mass (Table 2 footnote a)."
  )

  ini({
    # Table 2 (Naphthoquine column). Clearances and volumes refer to a patient
    # weighing 55 kg with a fat-free mass of 45 kg (Table 2 footnote a).
    lcl <- log(44.2)
    label("Apparent clearance CL/F for a 45 kg fat-free-mass patient (L/h)") # Table 2 'CL (L/h)' = 44.2 (37.9-50.6)
    lvc <- log(647)
    label("Apparent central volume Vc/F for a 55 kg patient (L)") # Table 2 'V1 (L)' = 647 (394-905)
    lq <- log(601)
    label("Apparent intercompartmental clearance Q/F for a 55 kg patient (L/h)") # Table 2 'Q (L/h)' = 601 (474-707)
    lvp <- log(19100)
    label("Apparent peripheral volume Vp/F for a 55 kg patient (L)") # Table 2 'Vp (L)' = 19100 (16,700-21,700)
    lka <- log(0.108)
    label("First-order absorption rate constant from the absorption compartment (1/h)") # Table 2 'Ka (1/h)' = 0.108 (0.0797-0.136)
    lmtt <- log(1.23)
    label("Absorption mean transit time (h)") # Table 2 'MTT (h)' = 1.23 (0.91-1.723)
    lnn <- log(5.42)
    label("Number of absorption transit compartments (unitless, non-integer)") # Table 2 'NN' = 5.42 (3.56-8.01)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability (unitless)") # Table 2 'F' = 1.00 fixed; Methods 'Relative bioavailability (F) was fixed to 1'

    # Allometric exponents fixed to the theory-based values (Methods 'using the
    # suggested exponents of 1 for volumes of distribution and 3/4 for
    # clearance terms'; Table 2 footnote a equations).
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on CL/F (unitless)") # Table 2 footnote a: 'CL/F = theta_pop * (FFM/45)^0.75 for naphthoquine'
    e_wt_q <- fixed(0.75)
    label("Allometric exponent of body weight on Q/F (unitless)") # Methods: 'exponents of ... 3/4 for clearance terms'; Table 2 footnote a: all clearances except naphthoquine CL scaled by body weight
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # Table 2 footnote a: 'V/F = theta_pop * (wt/55) for naphthoquine and artemisinin'

    # IIV reported as approximate %CV = sqrt(omega^2) * 100 (Table 2 footnote
    # d), so omega^2 = (CV/100)^2.
    etalcl ~ 0.039601 # Table 2 IIV 'CL' = 19.9 %CV -> 0.199^2
    etalka ~ 0.1369 # Table 2 IIV 'Ka' = 37.0 %CV -> 0.370^2
    etalmtt ~ 0.649636 # Table 2 IIV 'MTT' = 80.6 %CV -> 0.806^2
    etalfdepot ~ 0.106929 # Table 2 IIV 'F' = 32.7 %CV -> 0.327^2

    # Combined additive and proportional residual error (Methods).
    addSd <- 0.594
    label("Additive residual error (ng/mL)") # Table 2 'Additive error (ng/mL)' = 0.594 (0.345-0.892)
    propSd <- 0.251
    label("Proportional residual error (fraction)") # Table 2 'Proportional error (%)' = 25.1 (22.2-27.6)
  })
  model({
    cl <- exp(lcl + etalcl) * (FFM / 45)^e_ffm_cl
    vc <- exp(lvc) * (WT / 55)^e_wt_vc
    q <- exp(lq) * (WT / 55)^e_wt_q
    vp <- exp(lvp) * (WT / 55)^e_wt_vc
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Savic 2007 transit-compartment absorption (Methods reference 52): the
    # dose passes through nn transit compartments (ktr = (nn + 1) / mtt inside
    # transit()) into the absorption compartment, which empties into central
    # at the separately estimated rate ka (Figure 1).
    d / dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- 0

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
