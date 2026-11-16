Nakagawa_2026_tisotumabVedotin_tg32mouse_tab <- function() {
  description <- paste(
    "Preclinical (mouse, human FcRn transgenic Tg32).",
    "Two-compartment linear pharmacokinetic model for",
    "total antibody of tisotumab vedotin",
    "(anti-TF ADC; MMAE payload, DAR 4,",
    "cleavable valine-citrulline linker) after a single",
    "1 mg/kg intravenous bolus in male human FcRn transgenic Tg32 mice (Nakagawa 2026).",
    "The ADC does not cross-react with the mouse target antigen, so the",
    "disposition is target-independent and linear. Parameters are the means of",
    "individual two-compartment fits (WinNonlin, n = 3 mice) and are normalised",
    "to body weight: volumes in mL/kg, clearances in mL/day/kg, amounts in ug/kg.",
    "No between-animal variability or residual-error model is reported. Companion",
    "models: Nakagawa_2026_tisotumabVedotin_wtmouse_tab (same analyte, other strain),",
    "Nakagawa_2026_tisotumabVedotin_tg32mouse_adc (total ADC in this strain)",
    "and the same analyte and strain for the seven other ADCs of the panel."
  )
  reference <- paste(
    "Nakagawa K, Haraya K. Human Pharmacokinetic Prediction of Antibody-Drug",
    "Conjugates Using Human FcRn Transgenic Mice. Clin Transl Sci.",
    "2026;19:e70692. doi:10.1111/cts.70692"
  )
  vignette <- "Nakagawa_2026_adc_fcrn_mouse"
  units <- list(time = "day", dosing = "ug/kg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "tisotumab vedotin total antibody", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tisotumab vedotin total antibody", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (human FcRn transgenic Tg32, B6.Cg-Fcgrt Tg(FCGRT)32Dcr/DcrJ)",
    n_subjects = 3,
    n_studies = 1,
    age_range = "5 weeks at purchase",
    sex_female_pct = 0,
    weight_range = "not reported; 20 g was assumed for mice in the allometric scaling",
    disease_state = "healthy, non-tumour-bearing; target-independent disposition (no cross-reactivity with mouse antigens)",
    dose_range = "tisotumab vedotin 1 mg/kg as a single intravenous bolus into the tail vein (10 mL/kg)",
    regions = "Japan",
    notes = paste(
      "Nakagawa 2026 Methods (Pharmacokinetic Study): n = 3 male mice per group,",
      "single tail-vein injection of 1 mg/kg in 10 mL/kg;",
      "jugular-vein blood at 5 min, 4 h, 7 h and 1, 3, 7, 14, 21 and 28 days;",
      "heparinised plasma assayed for total antibody (anti-human IgG sandwich electrochemiluminescence assay; measures conjugated plus unconjugated antibody).",
      "Compartmental and non-compartmental analyses in WinNonlin 8.5.2.4; the",
      "matching NCA summary is Table S1. The same mice supply the companion",
      "total ADC model."
    )
  )

  ini({
    # Means of the three individual two-compartment fits; the paper reports no
    # population (mixed-effects) estimates, so there are no etas.
    lcl <- log(6.2); label("Clearance of total antibody (mL/day/kg)") # Table 3, Tg32 Total antibody CL = 6.2 +/- 0.2 mL/day/kg (mean +/- SD, n = 3)
    lvc <- log(41.8); label("Central volume of distribution of total antibody (mL/kg)") # Table S3, Tg32 Total antibody Vc = 41.8 +/- 1.5 mL/kg
    lq <- log(117.1); label("Intercompartmental clearance of total antibody (mL/day/kg)") # Table S4, Tg32 Total antibody Q = 117.1 +/- 37.6 mL/day/kg
    lvp <- log(43.4); label("Peripheral volume of distribution of total antibody (mL/kg)") # Table S5, Tg32 Total antibody Vp = 43.4 +/- 8.0 mL/kg

    # Not reported: WinNonlin individual fits; no residual-error estimate is given.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)")
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ug/kg over mL/kg gives ug/mL, the assay units of Figures 1 and 2.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
