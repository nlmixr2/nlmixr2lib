Arshad_2020_fluorouracil <- function() {
  description <- "Semi-physiological PK/PD model of continuous-infusion 5-fluorouracil (5-FU) in adults with gastrointestinal cancer (Arshad 2020). 5-FU follows a two-compartment model with linear elimination; a fixed 85% of 5-FU clearance forms 5-fluoro-5,6-dihydrouracil (5FUH2), which follows a one-compartment model. A single linear body-surface-area effect is shared by the 5-FU and 5FUH2 clearances. 5-FU plasma concentration drives a Friberg-style leukocyte myelosuppression chain (proliferating pool, three transit compartments, circulating cells, feedback exponent fixed at 0.17) through a linear drug effect whose slope differs between 5-FU monotherapy and 5-FU plus cisplatin."
  reference <- "Arshad U, Ploylearmsaeng SA, Karlsson MO, Doroshyenko O, Langer D, Schomig E, Kunze S, Guner SA, Skripnichenko R, Ullah S, Jaehde U, Fuhr U, Jetter A, Taubert M. Prediction of exposure-driven myelotoxicity of continuous infusion 5-fluorouracil by a semi-physiological pharmacokinetic-pharmacodynamic model in gastrointestinal cancer patients. Cancer Chemother Pharmacol. 2020;85(4):711-722. doi:10.1007/s00280-019-04028-5"
  vignette <- "Arshad_2020_fluorouracil"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    WBC = "10^9/L"
  )

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    central_5fuh2 = list(analyte = "5-fluoro-5,6-dihydrouracil", units = "mg", specimen = "plasma", verified = TRUE),
    precursor1 = list(
      analyte = "proliferating leukocyte progenitors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor2 = list(analyte = "maturing leukocytes", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    precursor3 = list(analyte = "maturing leukocytes", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    precursor4 = list(analyte = "maturing leukocytes", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    circ = list(analyte = "total leukocytes", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear centred effect (1 + e_bsa_cl * (BSA - 1.95)) on both the 5-FU and the 5FUH2 clearance; 1.95 m^2 is the cohort median (Table 1). Table 3 footnote a: 'Fractional change in CL per m^2 difference from median BSA value'. The BSA formula (DuBois / Mosteller) is not stated in the paper.",
      source_name = "BSA"
    ),
    CONMED_CISPLATIN = list(
      description = "Concomitant cisplatin during the 5-FU cycle (1 = yes, 0 = 5-FU monotherapy)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the drug-effect slope: 1 -> Slope_comb (5-FU + cisplatin 20 mg/m^2/day, the 14 oesophageal-cancer patients), 0 -> Slope_mono (the 16 colorectal/anal-cancer patients on 5-FU plus radiotherapy). The paper estimated one slope per group rather than a reference plus an offset.",
      source_name = "cisplatin co-medication"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "37-73 years (median 59.5)",
    weight_range = "46-111 kg (median 76)",
    bsa_range = "1.48-2.33 m^2 (median 1.95)",
    sex_female_pct = 100 * 5 / 30,
    race_ethnicity = "Not reported (single German centre)",
    disease_state = "Gastrointestinal cancer: oesophageal (n = 14, 5-FU + cisplatin) and colorectal / rectal / anal (n = 16, 5-FU + radiotherapy). No prior chemotherapy or radiotherapy; Karnofsky performance status >= 70%.",
    dose_range = "5-FU 650 or 1000 mg/m^2/day as a 5-day continuous IV infusion (first cycle only); cisplatin 20 mg/m^2/day for 5 days in the oesophageal-cancer patients.",
    regions = "Germany (University Hospital Cologne, 2002-2005)",
    baseline_wbc = "Median 6.90 x 10^9/L (range 4.68-11.28)",
    notes = "199 5-FU and 251 5FUH2 plasma concentrations; 135 total WBC counts from 29 patients. Estimation with FOCE-I in NONMEM 7.4.2. Demographics from Table 1."
  )

  ini({
    # Values are the bootstrap medians of Arshad 2020 Table 3 (the values the
    # Abstract and Results quote, and the values that reproduce the paper's
    # Figure 3 typical-subject simulation); the NONMEM point estimate of the
    # same row is given in each comment. See the vignette for the evidence.
    lcl <- log(249); label("5-FU total clearance CL5FU at BSA 1.95 m^2 (L/h)") # Table 3 'CL5FU (L/h)' bootstrap 249 (NONMEM 256)
    lvc <- log(5.56); label("5-FU central volume VC,5FU (L)") # Table 3 'VC,5FU (L)' bootstrap 5.56 (NONMEM 5.85)
    lvp <- log(28.5); label("5-FU peripheral volume VP,5FU (L)") # Table 3 'VP,5FU (L)' bootstrap 28.5 (NONMEM 24.0)
    lq <- log(14.8); label("5-FU intercompartmental clearance Q (L/h)") # Table 3 'Q (L/h)' bootstrap 14.8 (NONMEM 17.3)
    e_bsa_cl <- 0.77; label("Linear BSA effect shared by CL5FU and CL5FUH2 (fractional change per m^2 from 1.95 m^2)") # Table 3 'BSA effect (m-2)' bootstrap 0.77 (NONMEM 0.71); Abstract '77%/m2'

    fm <- fixed(0.85); label("Fraction of 5-FU clearance forming 5FUH2 (unitless)") # Table 3 'Fm (%)' 85, Fixed; Methods 'fixed a priori to 0.85'
    lcl_5fuh2 <- log(121); label("5FUH2 clearance CL5FUH2 at BSA 1.95 m^2 (L/h)") # Table 3 'CL5FUH2 (L/h)' bootstrap 121 (NONMEM 124)
    lvc_5fuh2 <- log(96.7); label("5FUH2 volume of distribution VC,5FUH2 (L)") # Table 3 'VC,5FUH2 (L)' bootstrap 96.7 (NONMEM 100)

    lcirc0 <- log(6.86); label("Baseline circulating leukocyte count CIRC0 (10^9/L)") # Table 3 'CIRC0' bootstrap 6.86 (NONMEM 7.16)
    lmtt <- log(281); label("Mean transit time MTT = 4/ktr (h)") # Table 3 'MTT (h)' bootstrap 281 (NONMEM 261)
    lslope_mono <- log(1.17); label("Linear drug-effect slope, 5-FU monotherapy (L/mg)") # Table 3 'Slopemono (L/mg)' bootstrap 1.17 (NONMEM 1.31)
    lslope_comb <- log(2.82); label("Linear drug-effect slope, 5-FU plus cisplatin (L/mg)") # Table 3 'Slopecomb (L/mg)' bootstrap 2.82 (NONMEM 2.10)
    gamma <- fixed(0.17); label("Feedback exponent on (CIRC0/circ) (unitless)") # Table 3 'gamma' 0.17 Fixed; Results 'fixed to a value of 0.17 according to the available literature'

    # IIV: Table 3 reports %CV; converted as omega^2 = log(CV^2 + 1).
    etalcl ~ 0.05155 # Table 3 IIV 'CL5FU' bootstrap 23.0 %CV (NONMEM 24.9)
    etalvc ~ 1.1322 # Table 3 IIV 'VC,5FU' bootstrap 145 %CV (NONMEM 130)
    etalcl_5fuh2 ~ 0.08022 # Table 3 IIV 'CL5FUH2' bootstrap 28.9 %CV (NONMEM 30.5)
    etalvc_5fuh2 ~ 0.30395 # Table 3 IIV 'VC,5FUH2' bootstrap 59.6 %CV (NONMEM 58.9)
    etalcirc0 ~ 0.02654 # Table 3 IIV 'CIRC0' bootstrap 16.4 %CV (NONMEM 16.8)

    # Residual error: Table 3 header 'RUV (sigma2)' -- variances, so SD = sqrt().
    propSd <- sqrt(0.32); label("Proportional residual SD, 5-FU plasma concentration (fraction)") # Table 3 'Proportional error 5FU' sigma2 bootstrap 0.32 (NONMEM 0.36)
    propSd_5fuh2 <- sqrt(0.14); label("Proportional residual SD, 5FUH2 plasma concentration (fraction)") # Table 3 'Proportional error 5FUH2' sigma2 0.14
    propSd_WBC <- sqrt(0.08); label("Proportional residual SD, total WBC count (fraction)") # Table 3 'Proportional error total WBC count' sigma2 0.08
  })
  model({
    # Shared linear BSA effect centred on the cohort median (Table 3 footnote a)
    bsa_eff <- 1 + e_bsa_cl * (BSA - 1.95)

    cl <- exp(lcl + etalcl) * bsa_eff
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    cl_5fuh2 <- exp(lcl_5fuh2 + etalcl_5fuh2) * bsa_eff
    vc_5fuh2 <- exp(lvc_5fuh2 + etalvc_5fuh2)

    circ0 <- exp(lcirc0 + etalcirc0)
    mtt <- exp(lmtt)
    # One slope per treatment group (Table 3 Slopemono / Slopecomb)
    slope <- exp(lslope_mono) * (1 - CONMED_CISPLATIN) + exp(lslope_comb) * CONMED_CISPLATIN
    # Three transit compartments; kprol = ktr = kcirc (Methods; Fig. 1 'MTT = 4/ktr')
    ktr <- 4 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_5fuh2 <- cl_5fuh2 / vc_5fuh2

    Cc <- central / vc
    # 5FUH2 formed on a mass basis (no molecular-weight correction stated)
    Cc_5fuh2 <- central_5fuh2 / vc_5fuh2

    # Linear drug effect on proliferation driven by the 5-FU plasma concentration
    edrug <- slope * Cc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_5fuh2) <- fm * kel * central - kel_5fuh2 * central_5fuh2

    d/dt(precursor1) <- ktr * precursor1 * (1 - edrug) * (circ0 / circ)^gamma - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4
    d/dt(circ) <- ktr * precursor4 - ktr * circ

    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0
    circ(0) <- circ0

    WBC <- circ

    Cc ~ prop(propSd)
    Cc_5fuh2 ~ prop(propSd_5fuh2)
    WBC ~ prop(propSd_WBC)
  })
}
