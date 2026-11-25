Huang_2026_imeglimin_metformin <- function() {
  description <- paste(
    "Deterministic linearly coupled two-compartment PK model for oral",
    "imeglimin (drug X, bare canonical names) and oral metformin (drug Y,",
    "sibling-drug suffix _metformin) in healthy adult volunteers (Huang",
    "2026, data from Chevalier 2020 and Fouqueray 2020). Same structure as",
    "Huang_2026_metoprolol_captopril: first-order absorption from a",
    "gastrointestinal depot and two-compartment disposition written with",
    "micro rate constants on AMOUNTS, where every rate constant of one drug",
    "is a linear function of the amount of the OTHER drug in the",
    "compartment that the process leaves (Huang 2026 equation 1). The",
    "single-drug rate constants were fitted first to the monotherapy",
    "profiles and the eight interaction slopes then to the combination",
    "profile by particle-swarm optimisation with an L2 penalty. No",
    "inter-individual variability and no residual-error model were",
    "reported; this is a typical-value model."
  )
  reference <- paste(
    "Huang H, Li C, Chen Q, Zhuang C, Yu L, Jin W, Li X. Coupled",
    "pharmacokinetic model unveils drug-drug interactions in plasma",
    "concentration. PLoS One. 2026;21(2):e0339052.",
    "doi:10.1371/journal.pone.0339052 (Supplementary S1 File, S7-S9",
    "Tables). Data: Chevalier C et al. Eur J Drug Metab Pharmacokinet.",
    "2020;45(6):725-33, doi:10.1007/s13318-020-00642-4; Fouqueray P et al.",
    "Clin Pharmacokinet. 2020;59(10):1261-71, doi:10.1007/s40262-020-00886-y."
  )
  vignette <- "Huang_2026_coupled_pk_ddi"

  # Amounts in mg; interaction slopes in 1/(g*h) as printed (S9 Table unit
  # header), so the other drug's amount is divided by 1000 inside model().
  # Concentrations in ng/mL to match Supplementary S7 Table.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "imeglimin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "imeglimin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imeglimin", units = "mg", specimen = "tissue", verified = TRUE),
    depot_metformin = list(analyte = "metformin", units = "mg", specimen = "administration site", verified = TRUE),
    central_metformin = list(analyte = "metformin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_metformin = list(analyte = "metformin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 2L,
    age_range = "18-60 years",
    disease_state = "healthy volunteers",
    dose_range = paste(
      "imeglimin 1500 mg PO and metformin 850 mg PO (Supplementary S7",
      "Table footnote); the combination profile is from repeated dosing",
      "(twice daily on days 1-6, once daily on days 7-12, sampled on day 12)"
    ),
    race_ethnicity = c(White = 100),
    notes = paste(
      "Two datasets as described in Huang 2026 'Data source'. Chevalier",
      "2020: 16 healthy Caucasian men and women aged 19-55 years, BMI",
      "18.5-29.9 kg/m^2; imeglimin alone on day 1 and with cimetidine 400 mg",
      "twice daily on day 8. Fouqueray 2020: 16 healthy Caucasian men aged",
      "18-60 years, BMI 20-29.9 kg/m^2; metformin 850 mg with imeglimin",
      "1500 mg. The model was fitted to the mean profiles in Supplementary",
      "S7 Table (footnote 'Weight is 79 kg'); the combination and",
      "metformin-alone profiles have a non-zero pre-dose concentration",
      "(steady state), the imeglimin-alone profile does not."
    )
  )

  ini({
    # Single-drug rate constants, imeglimin (drug X). Supplementary S8 Table,
    # row [20][21], columns alpha00 / alpha12 / alpha21 / alpha10 (h-1).
    # The column order is read from the table's MathType header cells.
    lka <- log(0.268); label("Imeglimin absorption rate constant alpha00 (1/h)") # S8 Table alpha00 = 0.268
    lk12 <- log(231.984); label("Imeglimin central-to-peripheral rate constant alpha12 (1/h)") # S8 Table alpha12 = 231.984
    lk21 <- log(4.519); label("Imeglimin peripheral-to-central rate constant alpha21 (1/h)") # S8 Table alpha21 = 4.519
    lkel <- log(28.344); label("Imeglimin elimination rate constant alpha10 (1/h)") # S8 Table alpha10 = 28.344

    # Single-drug rate constants, metformin (drug Y). S8 Table, row
    # [20][21], columns beta00 / beta12 / beta21 / beta10 (h-1).
    lka_metformin <- log(0.200); label("Metformin absorption rate constant beta00 (1/h)") # S8 Table beta00 = 0.200
    lk12_metformin <- log(25.194); label("Metformin central-to-peripheral rate constant beta12 (1/h)") # S8 Table beta12 = 25.194
    lk21_metformin <- log(1.299); label("Metformin peripheral-to-central rate constant beta21 (1/h)") # S8 Table beta21 = 1.299
    lkel_metformin <- log(18.969); label("Metformin elimination rate constant beta10 (1/h)") # S8 Table beta10 = 18.969

    # One apparent central volume for the drug pair (S8 Table prints
    # 4721.361 mL under every column); used only for the concentration
    # outputs (Huang 2026 equation 2).
    lvc <- log(4.721361); label("Apparent central volume shared by both drugs (L)") # S8 Table Vc = 4721.361 mL

    # Interaction slopes, effect of metformin on imeglimin (S9 Table, row
    # [20][21]); 1/(g*h).
    ddi_ka <- -0.178; label("Slope of imeglimin ka on metformin depot amount alpha00' (1/(g*h))") # S9 Table alpha00' = -0.178
    ddi_k12 <- 9945; label("Slope of imeglimin k12 on metformin central amount alpha12' (1/(g*h))") # S9 Table alpha12' = 9945
    ddi_k21 <- -2.729; label("Slope of imeglimin k21 on metformin peripheral amount alpha21' (1/(g*h))") # S9 Table alpha21' = -2.729
    ddi_kel <- -45.004; label("Slope of imeglimin kel on metformin central amount alpha10' (1/(g*h))") # S9 Table alpha10' = -45.004

    # Interaction slopes, effect of imeglimin on metformin (S9 Table).
    ddi_ka_metformin <- -0.062; label("Slope of metformin ka on imeglimin depot amount beta00' (1/(g*h))") # S9 Table beta00' = -0.062
    ddi_k12_metformin <- -3.224; label("Slope of metformin k12 on imeglimin central amount beta12' (1/(g*h))") # S9 Table beta12' = -3.224
    ddi_k21_metformin <- -4.141; label("Slope of metformin k21 on imeglimin peripheral amount beta21' (1/(g*h))") # S9 Table beta21' = -4.141
    # Positive as printed; the S1 File text says 'though negative' but then
    # 'Imeglimin facilitates Metformin elimination', which needs a positive slope.
    ddi_kel_metformin <- 1.13; label("Slope of metformin kel on imeglimin central amount beta10' (1/(g*h))") # S9 Table beta10' = 1.13
  })

  model({
    ka <- exp(lka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- exp(lkel)
    ka_metformin <- exp(lka_metformin)
    k12_metformin <- exp(lk12_metformin)
    k21_metformin <- exp(lk21_metformin)
    kel_metformin <- exp(lkel_metformin)
    vc <- exp(lvc)

    # Huang 2026 equation 1, with the other drug's amount converted mg -> g.
    ka_cpl <- ka + ddi_ka * depot_metformin / 1000
    k12_cpl <- k12 + ddi_k12 * central_metformin / 1000
    k21_cpl <- k21 + ddi_k21 * peripheral1_metformin / 1000
    kel_cpl <- kel + ddi_kel * central_metformin / 1000
    ka_cpl_metformin <- ka_metformin + ddi_ka_metformin * depot / 1000
    k12_cpl_metformin <- k12_metformin + ddi_k12_metformin * central / 1000
    k21_cpl_metformin <- k21_metformin + ddi_k21_metformin * peripheral1 / 1000
    kel_cpl_metformin <- kel_metformin + ddi_kel_metformin * central / 1000

    d/dt(depot) <- -ka_cpl * depot
    d/dt(central) <- ka_cpl * depot - k12_cpl * central + k21_cpl * peripheral1 - kel_cpl * central
    d/dt(peripheral1) <- k12_cpl * central - k21_cpl * peripheral1
    d/dt(depot_metformin) <- -ka_cpl_metformin * depot_metformin
    d/dt(central_metformin) <- ka_cpl_metformin * depot_metformin - k12_cpl_metformin * central_metformin +
      k21_cpl_metformin * peripheral1_metformin - kel_cpl_metformin * central_metformin
    d/dt(peripheral1_metformin) <- k12_cpl_metformin * central_metformin - k21_cpl_metformin * peripheral1_metformin

    # Equation 2 (x = Vc * cx) inverted; mg/L * 1000 = ng/mL.
    Cc <- central / vc * 1000
    Cc_metformin <- central_metformin / vc * 1000
  })
}
