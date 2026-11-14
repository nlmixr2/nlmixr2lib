Tan_2026_fluorouracil <- function() {
  description <- "Two-compartment population PK model of intravenous 5-fluorouracil (5-FU) with Michaelis-Menten elimination from the central compartment in adults with gastrointestinal cancer receiving bolus plus continuous infusion or continuous-infusion-only 5-FU (FOLFOX, FOLFIRI, FOLFIRINOX, FLOT) (Tan 2026). Body surface area enters Vmax as a power function normalised to 1.94 m^2 with an estimated exponent; inter-individual variability on Vmax only; combined proportional + (fixed, negligible) additive residual error. A cumulative plasma-AUC state mirrors the authors' AUC compartment."
  reference <- "Tan Z, Sancho-Araiz A, Voller S, Peeters SLJ, Manten T, Deenen MJ, Maring JG, Bet PM, Mathot RAA, Knibbe CAJ, Moes DJAR. Advancing Precision Dosing of 5-FU: Population PK Model Development, Limited Sampling Strategies, and Fit-for-Use Application. Clin Pharmacokinet. 2026;65(7):1027-1042. doi:10.1007/s40262-026-01645-1"
  vignette <- "Tan_2026_fluorouracil"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (BSA / 1.94)^e_bsa_vmax on Vmax; 1.94 m^2 is the cohort median (Results 'The median BSA was 1.94 m2'; Table 2 footnote 'TVVmax typical value of Vmax for a BSA of 1.94 m2'). The BSA formula is not stated in the paper.",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    auc_central = list(analyte = "5-fluorouracil", units = "mg*h/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 140L,
    n_studies = 4L,
    age_range = "IQR 56-70 years",
    age_median = "63 years",
    weight_range = "mean 78.64 kg (SD 16.84)",
    bsa_range = "mean 1.94 m^2 (SD 0.25); median 1.94 m^2",
    sex_female_pct = 45,
    race_ethnicity = "Not reported (Dutch centres)",
    disease_state = "Gastrointestinal cancer: colorectal 42.1%, pancreatic 34.3%, gastric / gastro-oesophageal junction 17.9%, other GI 5.7%. DPYD wild type 71.4% (the CombiGaPP study enrolled DPYD variant carriers on genotype-guided reduced doses).",
    dose_range = "5-FU as part of FOLFOX, FOLFIRI, FOLFIRINOX or FLOT: bolus median 385 mg/m^2 (IQR 293-400) plus a 46-48 h continuous infusion of median 2340 mg/m^2 (IQR 1790-2422), or a 24 h continuous infusion (FLOT, no bolus); 1-11 cycles per patient.",
    regions = "Netherlands (Amsterdam UMC, Catharina Hospital Eindhoven, Isala Hospital Zwolle)",
    notes = "Pooled CombiGaPP (intensive + sparse arms), PERFU, FUUT and FLOT prospective cohorts; 435 5-FU plasma concentrations (none below the quantification limit), median 2 samples per patient (range 1-17). Samples stabilised ex vivo with the DPD inhibitor gimeracil. Demographics from Table 1; NONMEM 7.4.4 FOCE-I."
  )

  ini({
    lvc <- log(19.7); label("Central volume of distribution V1 (L)") # Table 2 'V1 (L)' 19.7; ESM control stream THETA(3)
    lvp <- log(15); label("Peripheral volume of distribution V2 (L)") # Table 2 'V2 (L)' 15; ESM THETA(4)
    lq <- log(52.4); label("Intercompartmental clearance Q (L/h)") # Table 2 'Q (L/h)' 52.4; ESM THETA(5)
    lvmax <- log(1900); label("Maximum elimination rate Vmax at BSA 1.94 m^2 (mg/h)") # Table 2 'TVVmax (mg/h)' 1900; ESM THETA(6)
    lkm <- log(8.04); label("Michaelis-Menten constant Km (mg/L)") # Table 2 'KM (mg/L)' 8.04; ESM THETA(7)
    e_bsa_vmax <- 1.14; label("Power exponent of (BSA / 1.94) on Vmax (unitless)") # Table 2 'EXPO' 1.14; ESM THETA(8)

    etalvmax ~ 0.10 # Table 2 'ETA_Vmax (omega2)' 0.10; ESM $OMEGA 0.10

    # ESM $ERROR: W = SQRT(PROPER**2 * IPRED**2 + ADDER**2) with $SIGMA 1 FIX,
    # so THETA(1) and THETA(2) are standard deviations (combined2 form).
    propSd <- 0.34; label("Proportional residual SD (fraction)") # Table 2 'Proportional error (%)' 34%; ESM THETA(1)
    addSd <- fixed(0.0001); label("Additive residual SD (mg/L)") # ESM THETA(2) '(0.0001) FIX ; additive error'
  })
  model({
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    vmax <- exp(lvmax + etalvmax) * (BSA / 1.94)^e_bsa_vmax
    km <- exp(lkm)

    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    d/dt(central) <- -vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Cumulative plasma AUC (ESM control stream compartment 3, DADT(3) = C1)
    d/dt(auc_central) <- Cc

    Cc ~ add(addSd) + prop(propSd)
  })
}
