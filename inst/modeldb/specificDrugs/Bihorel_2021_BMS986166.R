Bihorel_2021_BMS986166 <- function() {
  description <- "Joint population PK model for BMS-986166 (an oral prodrug sphingosine-1-phosphate-1 receptor modulator) and its active phosphorylated metabolite BMS-986166-P in healthy adults (Bihorel 2021). Parent: two-compartment model with first-order absorption, linear elimination and study-specific relative bioavailability. Metabolite: one-compartment model with Michaelis-Menten elimination, fed by complete molar conversion of the parent clearance (FM = 1) plus a virtual pre-systemic BMS-986166-P dose (a fraction F2 of the molar BMS-986166 dose) delivered zero-order over D2 into a metabolite depot and then absorbed first-order. Both bioavailabilities carry study-specific fold shifts for the 12.5% of subjects the analysts identified as having low exposures. The parent PK was fitted first and then fixed while the metabolite parameters were estimated (sequential fit)."
  reference <- "Bihorel S, Singhal S, Shevell D, Sun H, Xie J, Basdeo S, Liu A, Dutta S, Ludwig E, Huang H, Lin K, Fura A, Throup J, Girgis IG. Population Pharmacokinetic Analysis of BMS-986166, a Novel Selective Sphingosine-1-Phosphate-1 Receptor Modulator, and Exposure-Response Assessment of Lymphocyte Counts and Heart Rate in Healthy Participants. Clin Pharmacol Drug Dev. 2021;10(1):8-21. doi:10.1002/cpdd.878. PMCID: PMC7821288."
  vignette <- "Bihorel_2021_BMS986166"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    STUDY_IM018003 = list(
      description = "Study IM018003 (multiple-ascending dose) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = Study IM018001 (single-ascending dose)",
      notes = "Selects the study-specific relative bioavailabilities of BMS-986166 (F1) and of the virtual pre-systemic BMS-986166-P dose (F2), and the study-specific low-exposure fold shifts on both (Bihorel 2021 Table 1 and Figure 1: F1 = F1study x F1study,low; F2 = F2study x F2study,low). F1 in IM018001 is the reference and is fixed to 1. Bihorel 2021 Discussion cautions that the F1 difference may compensate for time-dependent PK changes after repeated dosing rather than being a true bioavailability difference, so the indicator is kept tied to the study rather than to single- vs multiple-dose regimens.",
      source_name = "study (IM018001 / IM018003)"
    ),
    MIX_LOW_EXPOSURE = list(
      description = "Low-exposure subpopulation indicator (1 = subject identified by the analysts as having low BMS-986166 and BMS-986166-P exposures)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = normal-exposure subject",
      notes = "Six of the 48 subjects in the population PK analysis (12.5%) exhibited lower exposures to both analytes (Bihorel 2021 Results, Population PK Model; Discussion). The class was assigned by visual inspection of the data before the fit (Simulations: 'subjects considered to have normal or low (based on visual inspection) exposures') and entered as a fixed subject-level flag, not estimated as a NONMEM $MIXTURE probability. No subject characteristic (age, weight, renal or hepatic function markers, ethnicity) distinguished the two classes. Multiplies F1 by 0.655 (IM018001) or 0.471 (IM018003) and F2 by 1.25 (IM018001) or 1.82 (IM018003). Set to 0 for a typical subject; for population simulation draw MIX_LOW_EXPOSURE ~ Bernoulli(0.125).",
      source_name = "low-exposure subject flag (Table 1 'LOW' subscript)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "BMS-986166", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "BMS-986166", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "BMS-986166", units = "mg", specimen = "tissue", verified = TRUE),
    depot_bms986166p = list(
      analyte = "BMS-986166-P",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_bms986166p = list(analyte = "BMS-986166-P", units = "umol", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 2L,
    age_range = "19-52 years",
    age_median = "33 years",
    weight_range = "58.3-104 kg",
    weight_median = "83.6 kg",
    sex_female_pct = 4.17,
    race_ethnicity = c(White = 54.2, Black = 41.7, AmericanIndianAlaskaNative = 4.17),
    disease_state = "Healthy adult volunteers with normal renal and hepatic function",
    dose_range = "BMS-986166 oral liquid formulation: single doses of 0.75, 2 or 5 mg (Study IM018001) or 0.25, 0.75 or 1.5 mg once daily for 28 days (Study IM018003)",
    regions = "United States (PPD Development, Austin, Texas)",
    notes = "Population PK analysis set: 48 actively treated subjects (24 per study), 1677 BMS-986166 and 1637 BMS-986166-P blood-lysate concentrations, none below the 0.100 ng/mL lower limit of quantification (Bihorel 2021 Results). Demographics from Supplementary Table S2 (overall column: 46 male, 2 female; 26 White, 20 Black or African American, 2 American Indian or Alaska Native; 14 Hispanic). Every subject received a lead-in placebo dose on day -1."
  )

  ini({
    # ------------------------------------------------------------------
    # BMS-986166 (parent). Bihorel 2021 Table 1. Fitted first and then
    # fixed in the combined model (Table 1 'Sequential Model'; Results
    # bullet 'PK of BMS-986166 were fixed to the population mean
    # estimates'). They are estimates of the parent fit, so they are not
    # wrapped in fixed() here.
    # ------------------------------------------------------------------
    lka <- log(0.287); label("BMS-986166 first-order absorption rate constant KA (1/h)") # Table 1 KA 0.287 1/h (RSE 6.65%)
    lcl <- log(2.51); label("BMS-986166 apparent elimination clearance CL (L/h)") # Table 1 CL 2.51 L/h (RSE 5.02%)
    lvc <- log(913); label("BMS-986166 apparent central volume V2 (L)") # Table 1 V2 913 L (RSE 3.63%)
    lq <- log(0.763); label("BMS-986166 apparent distribution clearance Q (L/h)") # Table 1 Q 0.763 L/h (RSE 13.5%)
    lvp <- log(205); label("BMS-986166 apparent peripheral volume V3 (L)") # Table 1 V3 205 L (RSE 11.2%)

    # Study-specific relative bioavailability of BMS-986166 (Figure 1, F1study)
    lfdepot_im018001 <- fixed(log(1)); label("BMS-986166 relative bioavailability F1 in Study IM018001 (reference)") # Table 1 F1 SAD 1.00 FIXED
    lfdepot_im018003 <- log(0.764); label("BMS-986166 relative bioavailability F1 in Study IM018003 (fraction)") # Table 1 F1 MAD 0.764 (RSE 4.87%)
    # Low-exposure fold change in F1 (Figure 1, F1study,low), log scale
    e_mix_low_exposure_fdepot_im018001 <- log(0.655); label("Log fold change in F1 for low-exposure subjects, Study IM018001") # Table 1 F1 SAD,LOW 0.655 (RSE 5.91%)
    e_mix_low_exposure_fdepot_im018003 <- log(0.471); label("Log fold change in F1 for low-exposure subjects, Study IM018003") # Table 1 F1 MAD,LOW 0.471 (RSE 13.6%)

    # ------------------------------------------------------------------
    # BMS-986166-P (active metabolite). Bihorel 2021 Table 1. VMAX and
    # KM are molar (umol/h, uM); footnote a marks them as highly
    # correlated (r^2 >= 0.810).
    # ------------------------------------------------------------------
    fm <- fixed(1); label("Fraction of BMS-986166 clearance converted to BMS-986166-P (FM)") # Table 1 FM 1.00 FIXED
    lvmax_bms986166p <- log(0.649); label("BMS-986166-P apparent maximum elimination rate VMAX (umol/h)") # Table 1 VMAX 0.649 umol/h (RSE 19.4%)
    lkm_bms986166p <- log(0.125); label("BMS-986166-P Michaelis-Menten constant KM (uM)") # Table 1 KM 0.125 uM (RSE 19.5%)
    lvc_bms986166p <- log(38.7); label("BMS-986166-P apparent central volume VM (L)") # Table 1 VM 38.7 L (RSE 7.65%)
    lka_bms986166p <- log(0.382); label("BMS-986166-P first-order absorption rate constant KAM (1/h)") # Table 1 KAM 0.382 1/h (RSE 5.37%)
    ld1_bms986166p <- log(5.56); label("Duration D2 of the zero-order input of the virtual BMS-986166-P dose (h)") # Table 1 D2 5.56 h (RSE 1.05%)

    # Study-specific fraction of the molar dose absorbed as BMS-986166-P (Figure 1, F2study)
    lfdepot_bms986166p_im018001 <- log(0.0771); label("Fraction F2 of the molar dose absorbed as BMS-986166-P, Study IM018001 (fraction)") # Table 1 F2 SAD 0.0771 (RSE 7.79%)
    lfdepot_bms986166p_im018003 <- log(0.0797); label("Fraction F2 of the molar dose absorbed as BMS-986166-P, Study IM018003 (fraction)") # Table 1 F2 MAD 0.0797 (RSE 7.94%)
    # Low-exposure fold change in F2 (Figure 1, F2study,low), log scale
    e_mix_low_exposure_fdepot_bms986166p_im018001 <- log(1.25); label("Log fold change in F2 for low-exposure subjects, Study IM018001") # Table 1 F2 SAD,LOW 1.25 (RSE 16.8%)
    e_mix_low_exposure_fdepot_bms986166p_im018003 <- log(1.82); label("Log fold change in F2 for low-exposure subjects, Study IM018003") # Table 1 F2 MAD,LOW 1.82 (RSE 22.4%)

    # ------------------------------------------------------------------
    # IIV. Table 1 reports magnitudes as %CV of exponential (log-normal)
    # IIV models (Methods); converted with omega^2 = log(1 + CV^2). The
    # paper reports no off-diagonal elements, so the matrix is diagonal.
    # ------------------------------------------------------------------
    etalcl ~ 0.111842 # Table 1 CL IIV 34.4 %CV -> log(1 + 0.344^2)
    etalvc ~ 0.027833 # Table 1 V2 IIV 16.8 %CV -> log(1 + 0.168^2)
    etalka ~ 0.167486 # Table 1 KA IIV 42.7 %CV -> log(1 + 0.427^2)
    etalvmax_bms986166p ~ 0.110003 # Table 1 VMAX IIV 34.1 %CV -> log(1 + 0.341^2)
    etalvc_bms986166p ~ 0.116183 # Table 1 VM IIV 35.1 %CV -> log(1 + 0.351^2)
    etalka_bms986166p ~ 0.065901 # Table 1 KAM IIV 26.1 %CV -> log(1 + 0.261^2)

    # ------------------------------------------------------------------
    # Residual error. Table 1 gives the residual variance and its %CV;
    # %CV = 100 * sqrt(variance) (0.00906 -> 9.52; 0.0174 -> 13.2), i.e.
    # a constant-CV (proportional) error on each analyte.
    # ------------------------------------------------------------------
    propSd <- 0.0952; label("BMS-986166 proportional residual error (fraction)") # Table 1 residual variability 0.00906 (9.52 %CV) -> sqrt(0.00906)
    propSd_bms986166p <- 0.1319; label("BMS-986166-P proportional residual error (fraction)") # Table 1 residual variability 0.0174 (13.2 %CV) -> sqrt(0.0174)
  })

  model({
    # Molecular weights (g/mol) needed for the mg <-> umol conversions.
    # The paper does not print them. Non-paper provenance: PubChem CID
    # 118877516 (BMS-986166, C25H33NO2, 379.5 g/mol) and its phosphate
    # ester BMS-986166-P (C25H34NO5P, 459.5 g/mol); the 80.0 g/mol
    # difference is the HPO3 of the phosphate. The shared m/z 345 product
    # ion monitored for both analytes (Bioanalytical Methods) is
    # consistent with these formulas.
    mw <- 379.5
    mw_bms986166p <- 459.5

    # Study- and exposure-class-specific bioavailabilities (Figure 1:
    # F1 = F1study x F1study,low, F2 = F2study x F2study,low)
    fdepot <- exp(
      lfdepot_im018001 * (1 - STUDY_IM018003) +
        lfdepot_im018003 * STUDY_IM018003 +
        MIX_LOW_EXPOSURE *
          (e_mix_low_exposure_fdepot_im018001 * (1 - STUDY_IM018003) +
            e_mix_low_exposure_fdepot_im018003 * STUDY_IM018003)
    )
    fdepot_bms986166p <- exp(
      lfdepot_bms986166p_im018001 * (1 - STUDY_IM018003) +
        lfdepot_bms986166p_im018003 * STUDY_IM018003 +
        MIX_LOW_EXPOSURE *
          (e_mix_low_exposure_fdepot_bms986166p_im018001 * (1 - STUDY_IM018003) +
            e_mix_low_exposure_fdepot_bms986166p_im018003 * STUDY_IM018003)
    )

    # Parent individual parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # Metabolite individual parameters
    vmax_bms986166p <- exp(lvmax_bms986166p + etalvmax_bms986166p)
    km_bms986166p <- exp(lkm_bms986166p)
    vc_bms986166p <- exp(lvc_bms986166p + etalvc_bms986166p)
    ka_bms986166p <- exp(lka_bms986166p + etalka_bms986166p)
    d1_bms986166p <- exp(ld1_bms986166p)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # BMS-986166-P molar concentration (uM = umol/L) driving its
    # saturable elimination VMAX * C / (KM + C)
    cm_bms986166p <- central_bms986166p / vc_bms986166p

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Virtual pre-systemic BMS-986166-P dose: zero-order input over D2 into
    # this depot, then first-order transfer KAM to the metabolite central
    d/dt(depot_bms986166p) <- -ka_bms986166p * depot_bms986166p
    # Irreversible molar conversion of the parent clearance (mg/h ->
    # umol/h via 1000 / mw) plus the absorbed virtual dose, minus
    # Michaelis-Menten elimination
    d/dt(central_bms986166p) <- ka_bms986166p * depot_bms986166p +
      fm * kel * central * 1000 / mw -
      vmax_bms986166p * cm_bms986166p / (km_bms986166p + cm_bms986166p)

    f(depot) <- fdepot
    # Dose records to depot_bms986166p carry the SAME mg amount as the
    # BMS-986166 dose (rate = -2 so the modelled duration applies);
    # multiplying by 1000 / mw turns it into the molar amount (umol) and
    # F2 takes the absorbed fraction.
    f(depot_bms986166p) <- fdepot_bms986166p * 1000 / mw
    dur(depot_bms986166p) <- d1_bms986166p

    # Blood-lysate concentrations in ng/mL: mg/L x 1000 for the parent;
    # uM x g/mol = ug/L = ng/mL for the metabolite
    Cc <- 1000 * central / vc
    Cc_bms986166p <- cm_bms986166p * mw_bms986166p

    Cc ~ prop(propSd)
    Cc_bms986166p ~ prop(propSd_bms986166p)
  })
}
