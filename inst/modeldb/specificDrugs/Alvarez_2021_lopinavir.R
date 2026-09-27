Alvarez_2021_lopinavir <- function() {
  description <- paste(
    "One-compartment first-order-absorption population PK model for oral",
    "lopinavir boosted by ritonavir (400/100 mg BID) in 13 hospitalised",
    "adults with Covid-19 (10 in intensive care). Apparent oral clearance",
    "CL/F is inhibited by the per-subject ritonavir trough concentration",
    "through a fixed Imax model, CL/F = CL0/F * (1 - Imax * C / (IC50 + C)),",
    "with ka, Imax and IC50 fixed to the Dickinson 2011 healthy-volunteer",
    "values. IIV on CL/F and V/F; residual error is proportional plus a",
    "fixed additive part (Alvarez 2021)."
  )
  reference <- paste(
    "Alvarez JC, Moine P, Davido B, Etting I, Annane D, Larabi IA, Simon N;",
    "Garches COVID-19 Collaborative Group. Population pharmacokinetics of",
    "lopinavir/ritonavir in Covid-19 patients. Eur J Clin Pharmacol.",
    "2021;77(3):389-397. doi:10.1007/s00228-020-03020-w.",
    "The fixed ka, Imax and IC50 are taken from Dickinson L et al.",
    "Antimicrob Agents Chemother. 2011;55(6):2775-2782 (paper reference 10)."
  )
  vignette <- "Alvarez_2021_lopinavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "lopinavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lopinavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_RTV_CC = list(
      description = "Ritonavir trough plasma concentration (per-subject, time-fixed)",
      units = "mg/L",
      type = "continuous",
      reference_category = "0 (no ritonavir; CL/F equals CL0/F)",
      notes = paste(
        "Measured ritonavir trough concentration (CresRTV) of each patient",
        "on 100 mg ritonavir BID, entering lopinavir CL/F through",
        "CL/F = CL0/F * (1 - Imax * CONMED_RTV_CC / (IC50 + CONMED_RTV_CC))",
        "(Alvarez 2021 Methods Equation 1 and Table 2 footer). The Methods",
        "text writes the unit as mg/mL; the Table 2 legend, the fixed IC50",
        "(0.057 mg/L) and the observed range (< 0.02 to 1.5 mg/L, Results)",
        "are all in mg/L, which is used here. The paper does not report the",
        "cohort median trough."
      ),
      source_name = "CresRTV"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Mean 64 +/- 16 years (Table 1); tested on CL/F and V/F and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Mean 85 +/- 15 kg (Table 1); tested and not retained."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Range 160-190 cm (Table 1); tested and not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Mean 27.9 +/- 5.4 kg/m^2 (Table 1); tested and not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "4 of 13 female (Table 1); tested and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Range 73-478 umol/L at sampling (Table 1); tested and not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Range 31-535 U/L at sampling (Table 1); tested and not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Range 37-297 U/L at sampling (Table 1); tested and not retained."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = paste(
        "Range 34-294 at sampling (Table 1, header printed as U/L; the",
        "Discussion quotes the same values in mg/L); tested and not retained."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 13L,
    n_studies = 1L,
    n_observations = 70L,
    age_range = "42-92 years (mean 64 +/- 16)",
    weight_range = "65-120 kg (mean 85 +/- 15)",
    sex_female_pct = 30.8,
    disease_state = paste(
      "Covid-19 (SARS-CoV-2 confirmed by RT-PCR and/or compatible chest CT);",
      "10 of 13 in intensive care, most intubated and ventilated; 10 of 13 with",
      "creatinine above the laboratory normal range; all with raised CRP."
    ),
    dose_range = paste(
      "Lopinavir/ritonavir (Kaletra) 400/100 mg orally twice daily; tablets",
      "crushed and given by nasogastric tube in ICU patients."
    ),
    regions = "France (Raymond Poincare Hospital, Garches)",
    notes = paste(
      "Retrospective analysis; 1-7 samples per patient (70 lopinavir",
      "concentrations), mostly troughs plus 2-12 h post-dose samples.",
      "Ritonavir troughs ranged from < 0.02 to 1.5 mg/L. NONMEM 7.4 FOCE-I."
    )
  )

  ini({
    lka <- fixed(log(0.572)); label("Absorption rate constant ka (1/h)") # Table 2 'KA (fixed)' = 0.572 1/h (from Dickinson 2011, ref 10)
    lcl <- log(4.88); label("Apparent clearance without ritonavir CL0/F (L/h)") # Table 2 'CL0/F' = 4.88 L/h
    lvc <- log(94.8); label("Apparent volume of distribution V/F (L)") # Table 2 'V/F' = 94.8 L
    ic50 <- fixed(0.057); label("Ritonavir trough concentration giving half of Imax (mg/L)") # Table 2 'IC50 (fixed)' = 0.057 mg/L (from Dickinson 2011)
    imax <- fixed(0.929); label("Maximum fractional inhibition of CL/F by ritonavir (unitless)") # Table 2 'Imax (fixed)' = 0.929 (from Dickinson 2011)

    # Table 2 'Inter-individual variability (omega)' rows, read as NONMEM
    # OMEGA variances of exponential etas (see vignette: a standard-deviation
    # reading cannot reproduce the Figure 4 ECDF percentages).
    etalcl ~ 2.881 # Table 2 IIV 'CL' = 2.881
    etalvc ~ 0.801 # Table 2 IIV 'V' = 0.801

    propSd <- 0.186; label("Proportional residual error (fraction)") # Table 2 RUV 'Proportional' = 0.186
    addSd <- fixed(0.071); label("Additive residual error (mg/L)") # Table 2 RUV 'Additive (fixed)' = 0.071 mg/L
  })
  model({
    # Methods Equation 1 / Table 2 footer:
    # CL/F = CL0/F * [1 - (Imax * CresRTV) / (IC50 + CresRTV)]
    inhib <- imax * CONMED_RTV_CC / (ic50 + CONMED_RTV_CC)
    cl <- exp(lcl + etalcl) * (1 - inhib)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
