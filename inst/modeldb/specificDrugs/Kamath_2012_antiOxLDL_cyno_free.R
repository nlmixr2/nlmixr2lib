Kamath_2012_antiOxLDL_cyno_free <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Two-compartment intravenous",
    "pharmacokinetic model for free anti-oxLDL, a human recombinant IgG1",
    "antibody to malondialdehyde-oxidized ApoB-100 (oxidized LDL), in healthy",
    "young cynomolgus monkeys given a single 10 mg/kg IV bolus (Kamath 2012).",
    "Serum anti-oxLDL was measured by Assay A, a ligand-binding ELISA",
    "(MDA-ApoB100 capture) that detects free antibody (unbound and partially",
    "bound). Kamath 2012 fitted a two-compartment model with IV bolus input",
    "and first-order elimination to each animal separately in WinNonlin and",
    "reported the mean of the four individual estimates (Table 1); the",
    "central volume, clearance and Vss are those means and the",
    "intercompartmental clearance is back-calculated so the terminal",
    "half-life equals the mean T1/2beta. All parameters are per kg body",
    "weight (dose in mg/kg, volumes in L/kg, clearances in L/day/kg). No",
    "between-animal variability or residual-error model is available. The",
    "companion total-antibody (Assay B) model is",
    "modellib('Kamath_2012_antiOxLDL_cyno_total')."
  )
  reference <- paste(
    "Kamath AV, Williams SP, Bullens S, Cowan KJ, Stenberg Y, Cherry SR,",
    "Rendig S, Kukis DL, Griesemer C, Damico-Beyer LA, Bunting S (2012).",
    "Pharmacokinetics and Biodistribution of a Human Monoclonal Antibody to",
    "Oxidized LDL in Cynomolgus Monkey Using PET Imaging. PLoS One",
    "7(9):e45116. doi:10.1371/journal.pone.0045116"
  )
  vignette <- "Kamath_2012_antiOxLDL"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "anti-oxLDL (free)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "anti-oxLDL (free)", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 4,
    n_studies = 1,
    age_range = "3 years 2 months to 7 years 10 months",
    sex_female_pct = 50,
    disease_state = "normal healthy young cynomolgus monkeys",
    dose_range = "anti-oxLDL approximately 10 mg/kg single IV bolus (per-animal AUCinf x CL implies 9.8-10.4 mg/kg; mixture of 64Cu-DOTA-labelled and unlabelled antibody, about 1.5 mCi/animal)",
    regions = "United States (UC Davis Center for Molecular and Genomic Imaging)",
    notes = paste(
      "Kamath 2012 Methods 'Study design in cynomolgus monkeys': two male and",
      "two female monkeys, serum sampled predose and at 5 min, 1 h, 24 h, 4,",
      "8, 15 and 29 days post-dose; lower limit of quantification 400 ng/mL.",
      "Body weights are not reported. One animal (2M) had a low-level",
      "anti-therapeutic-antibody signal already present predose, judged not",
      "drug-induced and of minimal impact on PK. Concentrations of the",
      "oxidized-LDL antigen were not measured. No race/ethnicity data apply."
    )
  )

  ini({
    # Kamath 2012 Table 1 (Assay A) reports per-animal and mean CL, V1, Vss,
    # T1/2alpha and T1/2beta in mL-based per-kg units; values are converted
    # to L/kg here (divide by 1000) so mg/kg / (L/kg) gives ug/mL.
    lcl <- fixed(log(14.2 / 1000)); label("Clearance (L/day/kg)") # Table 1 Assay A, Mean CL 14.2 mL/day/kg
    lvc <- fixed(log(57.3 / 1000)); label("Central volume of distribution (L/kg)") # Table 1 Assay A, Mean V1 57.3 mL/kg
    # Vss - V1 with Table 1 Assay A means: 152 - 57.3 = 94.7 mL/kg.
    lvp <- fixed(log(94.7 / 1000)); label("Peripheral volume of distribution (L/kg)") # derived: Table 1 Assay A Mean Vss 152 minus Mean V1 57.3 mL/kg
    # Q is not printed. With CL, V1 and V2 fixed it is the only free
    # parameter; 23.2829 mL/day/kg makes the beta half-life equal the Table 1
    # mean T1/2beta 9.40 day. The implied alpha half-life is 0.839 day against
    # the reported mean 0.831 day. The same back-solve on each animal's own
    # CL, V1, Vss and T1/2beta returns that animal's T1/2alpha to rounding
    # (vignette 'Source trace').
    lq <- fixed(log(23.2829 / 1000)); label("Intercompartmental clearance (L/day/kg)") # derived: back-solved from Table 1 Assay A Mean T1/2beta 9.40 day

    # Each animal was fitted separately in WinNonlin; the weighting scheme
    # and residual error are not reported.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)") # not reported in Kamath 2012
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
