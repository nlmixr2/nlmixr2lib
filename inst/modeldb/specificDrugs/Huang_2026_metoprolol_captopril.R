Huang_2026_metoprolol_captopril <- function() {
  description <- paste(
    "Deterministic linearly coupled two-compartment PK model for oral",
    "metoprolol (drug X, bare canonical names) and oral captopril (drug Y,",
    "sibling-drug suffix _captopril) given together to healthy adult",
    "volunteers (Huang 2026, data from Qu 2001). Each drug has first-order",
    "absorption from a gastrointestinal depot and two-compartment",
    "disposition written with micro rate constants on AMOUNTS. Every rate",
    "constant of one drug is a linear function of the amount of the OTHER",
    "drug in the compartment that the process leaves (depot for absorption,",
    "central for distribution out and elimination, peripheral for",
    "redistribution back), e.g. ka_metoprolol = ka + ddi_ka *",
    "depot_captopril (Huang 2026 equation 1); the slopes ddi_* are the",
    "interaction parameters, with negative values inhibiting and positive",
    "values promoting the process. With either drug given alone the slopes",
    "multiply zero and the model reduces to two ordinary independent",
    "two-compartment models. The single-drug rate constants were fitted",
    "first to the monotherapy arms and the eight slopes were then fitted to",
    "the combination arm by particle-swarm optimisation with an L2 penalty",
    "(hierarchical two-step fit). No inter-individual variability and no",
    "residual-error model were reported; this is a typical-value model."
  )
  reference <- paste(
    "Huang H, Li C, Chen Q, Zhuang C, Yu L, Jin W, Li X. Coupled",
    "pharmacokinetic model unveils drug-drug interactions in plasma",
    "concentration. PLoS One. 2026;21(2):e0339052.",
    "doi:10.1371/journal.pone.0339052. Data: Qu F, Zhang X, Song L. Studies",
    "on the pharmacokinetic interaction between metoprolol and captopril in",
    "healthy volunteers. Chin Pharm J. 2001;36(2):111-2."
  )
  vignette <- "Huang_2026_coupled_pk_ddi"
  dosing <- c("depot", "depot_captopril")

  # Amounts are carried in mg. Huang 2026 states amounts in g ('1 g of Cap
  # ... reduces the rate of Met absorption by 1.27 h-1'; Table 1 'Initial
  # amounts (10^3 mg)'), so the interaction slopes are in 1/(g*h) and the
  # other drug's amount is divided by 1000 inside model() to use them as
  # printed. Concentrations are reported in ng/mL to match the source data
  # (Supplementary S6 Table).
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "metoprolol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "metoprolol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "metoprolol", units = "mg", specimen = "tissue", verified = TRUE),
    depot_captopril = list(analyte = "captopril", units = "mg", specimen = "administration site", verified = TRUE),
    central_captopril = list(analyte = "captopril", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_captopril = list(analyte = "captopril", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "22-29 years",
    weight_range = "60-70 kg",
    disease_state = "healthy volunteers",
    dose_range = paste(
      "metoprolol 50 mg PO alone, captopril 50 mg PO alone, and metoprolol",
      "50 mg + captopril 50 mg PO together; single doses, randomised",
      "three-period crossover with 2-week washout"
    ),
    notes = paste(
      "Qu 2001 dataset as described in Huang 2026 'Data source' (n = 12,",
      "height 170-180 cm; three groups of four rotated through the three",
      "treatments). Sampling at 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4,",
      "6, 8 and 12 h. The model was fitted to the mean concentration-time",
      "profiles tabulated in Supplementary S6 Table (footnote 'Weight is 65",
      "kg'). Sex distribution and region not reported."
    )
  )

  ini({
    # Single-drug rate constants, metoprolol (drug X). Huang 2026 Table 7 and
    # Supplementary S8 Table, row [19], columns alpha00 / alpha12 / alpha21 /
    # alpha10. Table 7's unit header reads '10^-3 h-1' but the values are in
    # h-1: S8 Table and the worked example in 'The effect of Cap on Met'
    # (alpha00' * y0(0) / alpha00 = -1.270 * 0.05 / 0.205 = -30.9%) and the
    # Fig 9 fit (peak about 0.47 mg near 1.7 h) all require h-1.
    lka <- log(0.205); label("Metoprolol absorption rate constant alpha00 (1/h)") # Table 7 alpha00 = 0.205
    lk12 <- log(84.192); label("Metoprolol central-to-peripheral rate constant alpha12 (1/h)") # Table 7 alpha12 = 84.192
    lk21 <- log(8.122); label("Metoprolol peripheral-to-central rate constant alpha21 (1/h)") # Table 7 alpha21 = 8.122
    lkel <- log(15.162); label("Metoprolol elimination rate constant alpha10 (1/h)") # Table 7 alpha10 = 15.162

    # Single-drug rate constants, captopril (drug Y). Table 7 / S8 Table,
    # columns beta00 / beta12 / beta21 / beta10.
    lka_captopril <- log(0.716); label("Captopril absorption rate constant beta00 (1/h)") # Table 7 beta00 = 0.716
    lk12_captopril <- log(69.474); label("Captopril central-to-peripheral rate constant beta12 (1/h)") # Table 7 beta12 = 69.474
    lk21_captopril <- log(19.222); label("Captopril peripheral-to-central rate constant beta21 (1/h)") # Table 7 beta21 = 19.222
    lkel_captopril <- log(14.591); label("Captopril elimination rate constant beta10 (1/h)") # Table 7 beta10 = 14.591

    # Apparent central volume, one value for the drug pair (Table 7 prints
    # 4163.137 mL once for both drugs; S8 Table repeats it under every
    # column). Used only to convert the central amounts to concentrations
    # (Huang 2026 equation 2, x = Vc * cx); it does not enter the ODEs.
    lvc <- log(4.163137); label("Apparent central volume shared by both drugs (L)") # Table 7 Vc = 4163.137 mL

    # Interaction slopes, effect of captopril on metoprolol (Huang 2026
    # Table 8 and S9 Table, row [19]); 1/(g*h), multiplied by the captopril
    # amount in g in the compartment named in each label.
    ddi_ka <- -1.270; label("Slope of metoprolol ka on captopril depot amount alpha00' (1/(g*h))") # Table 8 alpha00' = -1.270
    ddi_k12 <- -729.2; label("Slope of metoprolol k12 on captopril central amount alpha12' (1/(g*h))") # Table 8 alpha12' = -729.2
    ddi_k21 <- 585.2; label("Slope of metoprolol k21 on captopril peripheral amount alpha21' (1/(g*h))") # Table 8 alpha21' = 585.2
    ddi_kel <- -4854; label("Slope of metoprolol kel on captopril central amount alpha10' (1/(g*h))") # Table 8 alpha10' = -4854

    # Interaction slopes, effect of metoprolol on captopril (Table 8 / S9).
    ddi_ka_captopril <- -2.354; label("Slope of captopril ka on metoprolol depot amount beta00' (1/(g*h))") # Table 8 beta00' = -2.354
    ddi_k12_captopril <- -9704; label("Slope of captopril k12 on metoprolol central amount beta12' (1/(g*h))") # Table 8 beta12' = -9704
    ddi_k21_captopril <- 4905; label("Slope of captopril k21 on metoprolol peripheral amount beta21' (1/(g*h))") # Table 8 beta21' = 4905
    # Table 8 prints -2481 (95% CI -2611.0, -2351.0) and the text calls it
    # negative ('Met has an inhibitory influence on the elimination of Cap');
    # S9 Table prints +2481 (2351.0, 2611.0). The main-text sign reproduces
    # Fig 11 (captopril peak about 1.68 mg); +2481 gives about 1.53 mg.
    ddi_kel_captopril <- -2481; label("Slope of captopril kel on metoprolol central amount beta10' (1/(g*h))") # Table 8 beta10' = -2481
  })

  model({
    ka <- exp(lka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- exp(lkel)
    ka_captopril <- exp(lka_captopril)
    k12_captopril <- exp(lk12_captopril)
    k21_captopril <- exp(lk21_captopril)
    kel_captopril <- exp(lkel_captopril)
    vc <- exp(lvc)

    # Huang 2026 equation 1: each rate constant is the single-drug value plus
    # a slope times the other drug's amount (mg / 1000 = g) in the
    # compartment the process leaves.
    ka_cpl <- ka + ddi_ka * depot_captopril / 1000
    k12_cpl <- k12 + ddi_k12 * central_captopril / 1000
    k21_cpl <- k21 + ddi_k21 * peripheral1_captopril / 1000
    kel_cpl <- kel + ddi_kel * central_captopril / 1000
    ka_cpl_captopril <- ka_captopril + ddi_ka_captopril * depot / 1000
    k12_cpl_captopril <- k12_captopril + ddi_k12_captopril * central / 1000
    k21_cpl_captopril <- k21_captopril + ddi_k21_captopril * peripheral1 / 1000
    kel_cpl_captopril <- kel_captopril + ddi_kel_captopril * central / 1000

    d/dt(depot) <- -ka_cpl * depot
    d/dt(central) <- ka_cpl * depot - k12_cpl * central + k21_cpl * peripheral1 - kel_cpl * central
    d/dt(peripheral1) <- k12_cpl * central - k21_cpl * peripheral1
    d/dt(depot_captopril) <- -ka_cpl_captopril * depot_captopril
    d/dt(central_captopril) <- ka_cpl_captopril * depot_captopril - k12_cpl_captopril * central_captopril +
      k21_cpl_captopril * peripheral1_captopril - kel_cpl_captopril * central_captopril
    d/dt(peripheral1_captopril) <- k12_cpl_captopril * central_captopril - k21_cpl_captopril * peripheral1_captopril

    # Equation 2 (x = Vc * cx) inverted; mg/L * 1000 = ng/mL.
    Cc <- central / vc * 1000
    Cc_captopril <- central_captopril / vc * 1000
  })
}
