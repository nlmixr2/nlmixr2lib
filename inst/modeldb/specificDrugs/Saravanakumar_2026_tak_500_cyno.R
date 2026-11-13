Saravanakumar_2026_tak_500_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Population PK/RO/PD model of TAK-500, an anti-CCR2 immune-cell-directed antibody-drug conjugate carrying the STING agonist payload dazostinag (TAK-676), after a single IV bolus of 0.1, 0.5 or 3 mg/kg (Saravanakumar 2026). The conjugate (ADC) and the deconjugated antibody each follow a two-compartment model with linear elimination kel and a shared Michaelis-Menten (CCR2 target-mediated) elimination whose saturation term sums both antibody species; a first-order deconjugation kdec converts ADC to deconjugated antibody. Outputs are the conjugated payload (DAR x ADC, Cc) and total antibody (ADC + deconjugated, Cc_tab), both in nM. Total antibody drives CCR2 receptor occupancy on monocytes (RO, %) through an Emax function above a baseline, and the drug-induced fraction of RO drives monocyte exit from the circulation (mono, cells/uL) to a peripheral pool (mono_periph) that returns them at kout. States hold nmol; dose in nmol of antibody (1 mg = 6.725 nmol, MW 148,700 g/mol)."
  reference <- "Saravanakumar A, Deng R, Dong L, Haridas S, Xia CQ, Appleman VA, Abu-Yousif AO, Piatkov K, Wang H. Population-Based Modeling to Predict Human PK/PD of TAK-500, an Anti-CCR2 Antibody-Drug Conjugate for First-in-Human Study in Cancer Patients. Clin Transl Sci. 2026;19(4):e70522. doi:10.1111/cts.70522"
  vignette <- "Saravanakumar_2026_tak_500"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nM (Cc conjugated TAK-676, Cc_tab total antibody); RO in %; mono in cells/uL"
  )

  # Paper-specific monocyte-trafficking pools (Saravanakumar 2026 Figure 1A
  # 'Mono (central)' / 'Mono (periphery)'), named after the bcell /
  # bcell_periph pair.
  paper_specific_compartments <- c("mono", "mono_periph")

  # Saravanakumar 2026 Methods 'Pharmacokinetic/Pharmacodynamic Model of
  # TAK-500 in Cynomolgus Monkeys': A1/A2 = conjugated antibody (ADC) in the
  # central / peripheral compartments, A3/A4 = deconjugated antibody, A5/A6 =
  # monocytes in the circulation / periphery (Figure 1A). Total antibody and
  # conjugated TAK-676 were measured by LC-MS/MS in plasma.
  compartmentData <- list(
    central = list(analyte = "TAK-500 (conjugated antibody)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "TAK-500 (conjugated antibody)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_dar0 = list(
      analyte = "deconjugated anti-CCR2 antibody",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_dar0 = list(
      analyte = "deconjugated anti-CCR2 antibody",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    mono = list(analyte = "circulating monocytes", units = "cells/uL", specimen = "blood cell", verified = TRUE),
    mono_periph = list(
      analyte = "monocytes that have left the circulation",
      units = "cells/uL (circulating-volume equivalent)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 12,
    n_studies = 1,
    weight_range = "2.6-4.5 kg",
    sex_female_pct = 100,
    disease_state = "Protein-therapeutic-naive healthy female cynomolgus monkeys.",
    dose_range = "Single IV slow bolus of TAK-500 at 0.1, 0.5 or 3 mg/kg (n = 4 per dose).",
    notes = "Saravanakumar 2026 Methods 'Pharmacokinetic/Pharmacodynamic Study of TAK-500 in Cynomolgus Monkeys'. Sampling at predose and 0, 0.83, 1, 6, 24, 72, 168, 240, 336 and 408 h post-dose for total antibody, conjugated TAK-676, CCR2 receptor occupancy and monocyte counts. Fitted in Phoenix NLME (FOCE-ELS). The model has no body-weight covariate: Vc, Vmax and the doses are absolute (per animal), so a mg/kg dose must be multiplied by the animal's weight before conversion to nmol.",
    dar = "Mean drug-to-antibody ratio fixed at 4.6 (Table 2)."
  )

  ini({
    # PK parameters -- Saravanakumar 2026 Table 2 'NLME estimates'.
    lvc <- log(0.113); label("Central volume of distribution (L)") # Table 2 'TV (Vc), L' 0.113 (%CV 4)
    lkel <- log(0.0253); label("Linear elimination rate constant from the central compartment (1/h)") # Table 2 'TV (kel), 1/h' 0.0253 (%CV 17)
    lvmax <- log(0.104); label("Maximum Michaelis-Menten (CCR2-mediated) elimination rate (nmol/h)") # Table 2 'TV (V max), nmole/h' 0.104 (%CV 31)
    lkm <- log(21.58); label("Antibody concentration at half-maximal Michaelis-Menten elimination (nM)") # Table 2 'TV (Km), nM' 21.58 (%CV 26)
    lk12 <- log(0.0177); label("Distribution rate constant from central to peripheral compartment, kpt (1/h)") # Table 2 'TV (kpt), 1/h' 0.0177 (%CV 16)
    lk21 <- log(0.0118); label("Distribution rate constant from peripheral to central compartment, ktp (1/h)") # Table 2 'TV (ktp),1/h' 0.0118 (%CV 11)
    dar <- fixed(4.6); label("Mean drug-to-antibody ratio of TAK-500 (mol TAK-676 per mol antibody)") # Table 2 'TV (DAR)' 4.6 Fixed
    lkdec <- log(0.0324); label("First-order deconjugation rate constant of the ADC to deconjugated antibody (1/h)") # Table 2 'TV (kdec), 1/h' 0.0324 (%CV 6); Results: deconjugation half-life 21 h

    # PD parameters -- Saravanakumar 2026 Table 2 'PD parameters'.
    lrbase_mono <- log(167); label("Baseline circulating monocyte count (cells/uL)") # Table 2 'TV (Monocytes), cell/uL' 167 (%CV 19)
    lk12_mono <- log(1.96); label("Maximal rate constant of monocyte exit from the circulation at full drug-induced occupancy, kin (1/h)") # Table 2 'TV (kin), 1/h' 1.96 (%CV 26)
    lk21_mono <- log(0.055); label("Rate constant of monocyte return from the periphery to the circulation, kout (1/h)") # Table 2 'TV (kout), 1/h' 0.055 (%CV 32)
    lrbase_ro <- log(7); label("Baseline (assay background) CCR2 receptor occupancy, RO0 (%)") # Table 2 'TV (RO_base), %' 7 (%CV 8)
    lemax_ro <- log(59); label("Maximal drug-induced CCR2 receptor occupancy above baseline, ROmax (%)") # Table 2 'TV (ROmax), %' 59 (%CV 6)
    lec50 <- log(5.31); label("Total-antibody concentration giving half of ROmax, EC50 (nM)") # Table 2 'TV (EC50), nM' 5.31 (%CV 20)
    lgamma <- log(1.24); label("Power exponent r on the drug-induced receptor-occupancy fraction driving monocyte exit (unitless)") # Table 2 'TV (r)' 1.24 (%CV 35)

    # IIV -- Table 2 'Intersubject variability, eta' column, read as the
    # Phoenix Omega variance x 100 (see the vignette for the evidence).
    etalkel ~ 0.06 # Table 2 'TV (kel)' eta 6 (shrinkage 18%)
    etalrbase_mono ~ 0.07 # Table 2 'TV (Monocytes)' eta 7 (shrinkage 7%)

    # Residual error -- proportional, Table 2 (percent -> fraction).
    propSd <- 0.31; label("Proportional residual error, conjugated payload (fraction)") # Table 2 'Proportional error model (conjugated payload) %' 31
    propSd_tab <- 0.33; label("Proportional residual error, total antibody (fraction)") # Table 2 'Proportional error (total Ab) %' 33
    propSd_RO <- 0.57; label("Proportional residual error, receptor occupancy (fraction)") # Table 2 'Proportional error model (receptor occupancy), %' 57
    propSd_mono <- 0.63; label("Proportional residual error, monocyte count (fraction)") # Table 2 'Proportional error (Monocytes), %' 63
  })

  model({
    vc <- exp(lvc)
    kel <- exp(lkel + etalkel)
    vmax <- exp(lvmax)
    km <- exp(lkm)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kdec <- exp(lkdec)

    rbase_mono <- exp(lrbase_mono + etalrbase_mono)
    k12_mono <- exp(lk12_mono)
    k21_mono <- exp(lk21_mono)
    rbase_ro <- exp(lrbase_ro)
    emax_ro <- exp(lemax_ro)
    ec50 <- exp(lec50)
    gamma <- exp(lgamma)

    # Equations 1-4: Michaelis-Menten elimination of each antibody species,
    # saturated by the SUM of conjugated and deconjugated antibody (amounts).
    mm_denom <- km * vc + central + central_dar0
    d/dt(central) <- -kel * central - vmax * central / mm_denom - k12 * central + k21 * peripheral1 - kdec * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_dar0) <- -kel * central_dar0 - vmax * central_dar0 / mm_denom - k12 * central_dar0 + k21 * peripheral1_dar0 + kdec * central
    d/dt(peripheral1_dar0) <- k12 * central_dar0 - k21 * peripheral1_dar0

    # Equations 5-6: total antibody (s2) and conjugated payload (s1), nM.
    Cc_tab <- (central + central_dar0) / vc
    Cc <- dar * central / vc

    # Equations 9-10: receptor occupancy (%) driven by total antibody, and
    # the drug-induced occupancy fraction s4.
    RO <- Cc_tab / (ec50 + Cc_tab) * emax_ro + rbase_ro
    ro_frac <- (RO - rbase_ro) * 0.01

    # Equations 11-12: monocyte trafficking out of and back into the
    # circulation. The periphery holds no drug-displaced monocytes at baseline.
    d/dt(mono) <- k21_mono * mono_periph - k12_mono * ro_frac^gamma * mono
    d/dt(mono_periph) <- k12_mono * ro_frac^gamma * mono - k21_mono * mono_periph
    mono(0) <- rbase_mono
    mono_periph(0) <- 0

    # Equations 7, 8, 13, 14: proportional residual error on every output.
    Cc ~ prop(propSd)
    Cc_tab ~ prop(propSd_tab)
    RO ~ prop(propSd_RO)
    mono ~ prop(propSd_mono)
  })
}
