Bousema_2010_malaria_gametocyte_density_tanzania <- function() {
  description <- paste(
    "Two-compartment gametocyte-density model of mature Plasmodium",
    "falciparum gametocytes after efficacious antimalarial treatment in",
    "Tanzanian children with uncomplicated malaria (Bousema 2010, Tanzania",
    "trial 2006, Pfs25 QT-NASBA data). Sequestered gametocytes S are",
    "released into the circulation at rate rho = 1/(mean sequestration",
    "time); circulating gametocytes G decay at rate mu = 1/(mean",
    "circulation time). No new gametocytes are produced from asexual",
    "stages after treatment starts. The model describes the second phase",
    "of carriage only: TIME ZERO IS DAY 3 AFTER THE START OF TREATMENT, and",
    "the initial states are the day-3 circulating density G3 and the day-3",
    "sequestered pool S3 (the drug effect of days 0-2 is folded into",
    "them). The mean circulation time and G3 differ between the ACT arm",
    "(sulfadoxine-pyrimethamine + artesunate; reference) and the ACT-PQ",
    "arm (the same ACT plus a single 0.75 mg/kg primaquine dose on the",
    "third treatment day), selected by TRT_PRIMAQUINE; S3 is shared by",
    "both arms. The patient random effect on overall log density and the",
    "residual SD on log density are not reported and are fixed to 0, so",
    "simulations are typical-value. Sibling models: the Kenya density fit",
    "(Bousema_2010_malaria_gametocyte_density_kenya) and the prevalence",
    "models (Bousema_2010_malaria_gametocyte_prevalence_kenya,",
    "Bousema_2010_malaria_gametocyte_prevalence_tanzania).",
    sep = " "
  )
  reference <- paste(
    "Bousema T, Okell L, Shekalaghe S, Griffin JT, Omar S, Sawa P,",
    "Sutherland C, Sauerwein R, Ghani AC, Drakeley C (2010). Revisiting the",
    "circulation time of Plasmodium falciparum gametocytes: molecular",
    "detection methods to estimate the duration of gametocyte carriage and",
    "the effect of gametocytocidal drugs. Malar J 9:136.",
    "doi:10.1186/1475-2875-9-136.",
    sep = " "
  )
  vignette <- "Bousema_2010_malaria_gametocyte"
  units <- list(
    time = "day (time zero = day 3 after the start of treatment)",
    dosing = "no dosing (treatment arm enters through TRT_PRIMAQUINE)",
    concentration = "gametocytes/uL blood"
  )

  paper_specific_compartments <- c("gametocyte_sequestered", "gametocyte_circulating")

  compartmentData <- list(
    gametocyte_sequestered = list(
      analyte = "Sequestered Plasmodium falciparum gametocytes awaiting release into the circulation",
      units = "gametocytes/uL blood",
      specimen = "tissue",
      verified = TRUE
    ),
    gametocyte_circulating = list(
      analyte = "Circulating mature (stage V) Plasmodium falciparum gametocytes",
      units = "gametocytes/uL blood",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    TRT_PRIMAQUINE = list(
      description = "Primaquine treatment-arm indicator: 1 = sulfadoxine-pyrimethamine + artesunate plus a single primaquine dose (ACT-PQ), 0 = sulfadoxine-pyrimethamine + artesunate alone (ACT, reference)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = ACT without primaquine",
      notes = paste(
        "Randomised arm of the Tanzania trial (Bousema 2010 Methods):",
        "primaquine 0.75 mg/kg once on the third treatment day. Selects the",
        "arm-specific mean circulation time and day-3 circulating density",
        "G3."
      ),
      source_name = "treatment group (ACT / ACT-PQ)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    age_range = "3-15 years (median 5, IQR 3-9; Table 1)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Uncomplicated Plasmodium falciparum malaria, 500-100,000 asexual",
      "parasites/uL at enrolment, Hb >= 8 g/dL; successfully cured of",
      "asexual parasites, no re-infection, and Pfs25 QT-NASBA gametocyte",
      "positive at least once between day 0 and day 28."
    ),
    dose_range = paste(
      "ACT: sulfadoxine-pyrimethamine (25 mg/kg S + 1.25 mg/kg P single",
      "dose) + artesunate 4 mg/kg once daily for 3 days. ACT-PQ: the same",
      "plus primaquine 0.75 mg/kg once on the third treatment day."
    ),
    regions = "Mnyuzi, northeastern Tanzania (EIR ~10 infectious bites per person per year)",
    notes = paste(
      "36 ACT and 41 ACT-PQ children were included in the model fitting",
      "(Table 1). Gametocyte density was measured by Pfs25 QT-NASBA on days",
      "0, 3, 7, 14, 28 and 42; data after day 28 were not used because of",
      "frequent re-infection, and only days 3-28 enter this model.",
      "Detection limit 0.02 gametocytes/uL; values below it were treated as",
      "left-censored in the fit. Fitted with SAS PROC NLMIXED."
    )
  )

  ini({
    lmst <- fixed(log(11)); label("Mean sequestration time per gametocyte, 1/rho (day)") # Table 2 '1/rho' = 1/11 for both trials; Methods: sequestration time fixed to 11 days (best fit within 4-12 days)

    lmct_act <- log(4.61); label("Mean circulation time per gametocyte, 1/mu, ACT arm (day)") # Figure 3 caption 4.61 days (2.92-7.26); Results 4.6 (95% CI 2.9-7.3); Additional file 1 best-fit row
    lmct_actpq <- log(0.53); label("Mean circulation time per gametocyte, 1/mu, ACT-PQ arm (day)") # Figure 3 caption 0.53 days (0.24-1.19); Results 0.5 (95% CI 0.2-1.2); Additional file 1 best-fit row

    lg3_act <- log(1.88); label("Circulating gametocyte density on day 3, G3, ACT arm (gametocytes/uL)") # Table 2 G3 ACT 1.88 (0.37-9.52)
    lg3_actpq <- log(0.07); label("Circulating gametocyte density on day 3, G3, ACT-PQ arm (gametocytes/uL)") # Table 2 G3 ACT-PQ 0.07 (0.01-0.86)

    ls3 <- log(0.0013); label("Sequestered gametocyte pool on day 3 later released, S3, both arms (gametocytes/uL)") # Table 2 S3 Tanzania 0.0013 (0.0001-0.0237); Additional file 1 best-fit row 0.0013 (ln -6.64)

    etalg3 ~ fixed(0) # Methods: patient random effect on overall log density; variance not reported

    expSd <- fixed(0); label("Residual SD of log gametocyte density (natural-log scale)") # Methods: Normal error on log density; SD not reported
  })

  model({
    rho <- 1 / exp(lmst)
    mu <- 1 / exp(lmct_act * (1 - TRT_PRIMAQUINE) + lmct_actpq * TRT_PRIMAQUINE)

    # One patient random effect scales the whole density (both pools), as in
    # the Methods: 'random effects for each patient in overall density'.
    g3 <- exp(lg3_act * (1 - TRT_PRIMAQUINE) + lg3_actpq * TRT_PRIMAQUINE + etalg3)
    s3 <- exp(ls3 + etalg3)

    gametocyte_sequestered(0) <- s3
    gametocyte_circulating(0) <- g3

    d/dt(gametocyte_sequestered) <- -rho * gametocyte_sequestered
    d/dt(gametocyte_circulating) <- rho * gametocyte_sequestered - mu * gametocyte_circulating

    gametocyte_density <- gametocyte_circulating
    gametocyte_density ~ lnorm(expSd)
  })
}
