Bousema_2010_malaria_gametocyte_density_kenya <- function() {
  description <- paste(
    "Two-compartment gametocyte-density model of mature Plasmodium",
    "falciparum gametocytes after efficacious antimalarial treatment in",
    "Kenyan children with uncomplicated malaria (Bousema 2010, Kenya trial",
    "2003-2004, Pfs25 QT-NASBA data). Sequestered gametocytes S are",
    "released into the circulation at rate rho = 1/(mean sequestration",
    "time); circulating gametocytes G decay at rate mu = 1/(mean",
    "circulation time). No new gametocytes are produced from asexual",
    "stages after treatment starts. The model describes the second phase",
    "of carriage only: TIME ZERO IS DAY 3 AFTER THE START OF TREATMENT, and",
    "the initial states are the day-3 circulating density G3 and the day-3",
    "sequestered pool S3 (the drug effect of days 0-2 is folded into",
    "them). The mean circulation time and G3 differ between the non-ACT",
    "arm (sulfadoxine-pyrimethamine + amodiaquine; reference) and the ACT",
    "arm (sulfadoxine-pyrimethamine + artesunate, or",
    "artemether-lumefantrine), selected by TRT_ARTEMISININ. The patient",
    "random effect on overall log density and the residual SD on log",
    "density are not reported and are fixed to 0, so simulations are",
    "typical-value. Sibling models: the Tanzania density fit",
    "(Bousema_2010_malaria_gametocyte_density_tanzania) and the",
    "prevalence models (Bousema_2010_malaria_gametocyte_prevalence_kenya,",
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
    dosing = "no dosing (treatment arm enters through TRT_ARTEMISININ)",
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
    TRT_ARTEMISININ = list(
      description = "Artemisinin-based combination therapy (ACT) arm indicator: 1 = sulfadoxine-pyrimethamine + artesunate or artemether-lumefantrine, 0 = sulfadoxine-pyrimethamine + amodiaquine (non-ACT, reference)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-ACT (sulfadoxine-pyrimethamine + amodiaquine)",
      notes = paste(
        "Bousema 2010 Data analysis pools the SP+AS and AL arms as one ACT",
        "group because both had an indistinguishable effect on gametocyte",
        "carriage. Selects the arm-specific mean circulation time, day-3",
        "circulating density G3 and day-3 sequestered pool S3."
      ),
      source_name = "treatment group (non-ACT / ACT)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 126L,
    n_studies = 1L,
    age_range = "0.5-10 years (median 3, IQR 1-5; Table 1)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Uncomplicated Plasmodium falciparum malaria, 500-100,000 asexual",
      "parasites/uL at enrolment, Hb >= 5 g/dL; successfully cured of",
      "asexual parasites, no re-infection, and Pfs25 QT-NASBA gametocyte",
      "positive at least once between day 0 and day 28."
    ),
    dose_range = paste(
      "Non-ACT: sulfadoxine-pyrimethamine (25 mg/kg S + 1.25 mg/kg P single",
      "dose) + amodiaquine 10 mg/kg once daily for 3 days. ACT:",
      "sulfadoxine-pyrimethamine + artesunate 4 mg/kg once daily for 3",
      "days, or artemether-lumefantrine (10 mg A + 60 mg L per 5 kg twice",
      "daily for 3 days)."
    ),
    regions = "Mbita, western Kenya (EIR ~70 infectious bites per person per year)",
    notes = paste(
      "36 non-ACT and 90 ACT children were included in the model fitting",
      "(Table 1). Gametocyte density was measured by Pfs25 QT-NASBA on days",
      "0, 3, 7, 14 and 28; only days 3-28 enter this model. Detection limit",
      "0.02 gametocytes/uL; values below it were treated as left-censored",
      "in the fit. Fitted with SAS PROC NLMIXED."
    )
  )

  ini({
    lmst <- fixed(log(11)); label("Mean sequestration time per gametocyte, 1/rho (day)") # Table 2 '1/rho' = 1/11 for both trials; Methods: sequestration time fixed to 11 days (best fit within 4-12 days)

    lmct_nonact <- log(6.53); label("Mean circulation time per gametocyte, 1/mu, non-ACT arm (day)") # Figure 2 caption 6.53 days (95% CI 4.84-8.80); Results 6.5; Additional file 1 best-fit row
    lmct_act <- log(5.04); label("Mean circulation time per gametocyte, 1/mu, ACT arm (day)") # Figure 2 caption 5.04 days (95% CI 4.20-6.06); Results 5.0; Additional file 1 best-fit row

    lg3_nonact <- log(1.45); label("Circulating gametocyte density on day 3, G3, non-ACT arm (gametocytes/uL)") # Table 2 G3 non-ACT 1.45 (0.44-4.74)
    lg3_act <- log(0.13); label("Circulating gametocyte density on day 3, G3, ACT arm (gametocytes/uL)") # Table 2 G3 ACT 0.13 (0.03-0.52)

    ls3_nonact <- fixed(log(0.0263)); label("Sequestered gametocyte pool on day 3 later released, S3, non-ACT arm (gametocytes/uL)") # Table 2 S3 non-ACT 0.0263; Additional file 1 best-fit row 0.0263 (ln -3.64)
    ls3_act <- fixed(log(0.0263) - 3); label("Sequestered gametocyte pool on day 3 later released, S3, ACT arm (gametocytes/uL)") # Methods: ACT pool 3 natural-log units below non-ACT; Table 2 S3 ACT 0.0013 (ln -6.64)

    etalg3 ~ fixed(0) # Methods: patient random effect on overall log density; variance not reported

    expSd <- fixed(0); label("Residual SD of log gametocyte density (natural-log scale)") # Methods: Normal error on log density; SD not reported
  })

  model({
    rho <- 1 / exp(lmst)
    mu <- 1 / exp(lmct_nonact * (1 - TRT_ARTEMISININ) + lmct_act * TRT_ARTEMISININ)

    # One patient random effect scales the whole density (both pools), as in
    # the Methods: 'random effects for each patient in overall density'.
    g3 <- exp(lg3_nonact * (1 - TRT_ARTEMISININ) + lg3_act * TRT_ARTEMISININ + etalg3)
    s3 <- exp(ls3_nonact * (1 - TRT_ARTEMISININ) + ls3_act * TRT_ARTEMISININ + etalg3)

    gametocyte_sequestered(0) <- s3
    gametocyte_circulating(0) <- g3

    d/dt(gametocyte_sequestered) <- -rho * gametocyte_sequestered
    d/dt(gametocyte_circulating) <- rho * gametocyte_sequestered - mu * gametocyte_circulating

    gametocyte_density <- gametocyte_circulating
    gametocyte_density ~ lnorm(expSd)
  })
}
