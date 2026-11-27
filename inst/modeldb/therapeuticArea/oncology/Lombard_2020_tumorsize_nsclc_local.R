Lombard_2020_tumorsize_nsclc_local <- function() {
  description <- paste(
    "Stein bi-exponential (exponential growth and decay) tumour growth",
    "inhibition model of the individual-lesion longest diameter in stage",
    "IIIA-IV non-small cell lung cancer treated with cisplatin alone or",
    "gemcitabine plus cisplatin (Lombard 2020), fit to the measurements",
    "made by the treating hospital's radiologists (the 'local' set). The",
    "same lesions read at the centralised review centre were fit",
    "separately and are packaged as Lombard_2020_tumorsize_nsclc_central.",
    "Each lesion is an individual: the baseline size, growth rate and",
    "decay rate carry log-normal inter-lesion variability. There is no",
    "drug input; the two treatment arms were pooled without an arm effect,",
    "so the decay rate is the pooled empirical drug effect."
  )
  reference <- paste(
    "Lombard A, Mistry H, Chapman SC, Gueoguieva I, Aarons L, Ogungbenro K.",
    "Impact of tumour size measurement inter-operator variability on",
    "model-based drug effect evaluation.",
    "Cancer Chemother Pharmacol. 2020;85(4):817-825.",
    "doi:10.1007/s00280-020-04049-5. PMID: 32170415. PMCID: PMC7125250.",
    "Source trial: Sandler AB et al., J Clin Oncol. 2000;18(1):122-130,",
    "doi:10.1200/JCO.2000.18.1.122.",
    sep = " "
  )
  vignette <- "Lombard_2020_tumorsize_nsclc"

  units <- list(
    time = "day",
    dosing = "n/a (no dose input; the pooled cisplatin / gemcitabine-plus-cisplatin effect is absorbed into the empirical decay rate kse)",
    concentration = "cm (the observable TS is the longest diameter of one lesion, not a drug concentration)"
  )

  compartmentData <- list(
    growth = list(
      analyte = "tumour size, growing (treatment-resistant) fraction of the Stein decomposition",
      units = "cm",
      specimen = "tissue",
      verified = TRUE
    ),
    shrink = list(
      analyte = "tumour size, decaying (treatment-sensitive) fraction of the Stein decomposition",
      units = "cm",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 62L,
    n_lesions = 120L,
    n_observations = 357L,
    n_studies = 1L,
    age_range = "not reported in this paper (see Sandler 2000 for the trial demographics)",
    weight_range = "not reported in this paper",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in this paper",
    disease_state = "stage IIIA, IIIB or IV non-small cell lung cancer",
    dose_range = "cisplatin alone or gemcitabine plus cisplatin for up to six cycles (Sandler 2000 phase III protocol; not model inputs)",
    regions = "28 study sites (countries not listed in this paper)",
    notes = paste(
      "Lombard 2020 Methods: a subset of 120 lesions in 62 of the 522 trial patients, selected because the same CT scan was measured both at the treating hospital (local) and at a centralised centre (central); 714 observations = 357 paired measurements, on average 3 (1-5) time points per lesion and a median of 59 (19-116) days between measurements.",
      "Baseline local longest diameter: median 3.0 cm (range 0.25-13.0); median 2 (1-6) lesions per patient.",
      "Measurements followed the WHO criteria; the unidimensional longest diameter was analysed to mimic RECIST. About 10% of measurements were below the 0.5 cm quantification limit and were replaced by 0.25 cm before fitting (no M3-type likelihood).",
      "The model was fitted in NONMEM 7.3 (FOCE with interaction) with each lesion as the individual; trial period August 1995 - February 1997.",
      sep = "\n"
    )
  )

  ini({
    # Lombard 2020 Table 1, 'Local set' column; Eq. 2 gives P_i = theta_P * exp(eta_i).
    lrbase <- log(3.1)
    label("Typical baseline lesion longest diameter TS(0) (cm)") # Table 1: TS(0) = 3.1 cm (5% RSE), Local set
    lkge <- log(0.0012)
    label("Tumour growth rate constant Kg (1/day)") # Table 1: Kg = 0.0012 /d (15% RSE), Local set
    lkse <- log(0.0077)
    label("Tumour decay (drug effect) rate constant Kd (1/day)") # Table 1: Kd = 0.0077 /d (11% RSE), Local set

    # Table 1 'IIV' rows, read as the NONMEM omega^2 variances (see the vignette
    # for why the SD reading is inconsistent with the Figure 3 individual estimates).
    etalrbase ~ 0.31 # Table 1: IIV-TS(0) = 0.31 (14% RSE), Local set
    etalkge ~ 0.20 # Table 1: IIV-Kg = 0.20 (33% RSE), Local set
    etalkse ~ 0.44 # Table 1: IIV-Kd = 0.44 (31% RSE), Local set

    addSd <- 0.43
    label("Additive residual error on the longest diameter (cm)") # Table 1: additive residual error = 0.43 cm (9% RSE), Local set
  })
  model({
    rbase <- exp(lrbase + etalrbase)
    kge <- exp(lkge + etalkge)
    kse <- exp(lkse + etalkse)

    # Lombard 2020 Eq. 1: TS(t) = TS(0) * (exp(Kg * t) + exp(-Kd * t) - 1),
    # written as two exponential sub-states sharing the initial condition TS(0).
    d/dt(growth) <- kge * growth
    d/dt(shrink) <- -kse * shrink
    growth(0) <- rbase
    shrink(0) <- rbase

    TS <- growth + shrink - rbase

    TS ~ add(addSd)
  })
}
