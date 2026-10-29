Huygens_2024_antiSarsCov2IgG <- function() {
  description <- "Two-compartment population PK model for passively transferred anti-SARS-CoV-2 spike IgG (binding-antibody titre) after intravenous convalescent plasma (ConvP) or hyperimmune anti-SARS-CoV-2 globulin (COVIg) in seronegative immunocompromised adults (Huygens 2024)"
  reference <- paste(
    "Huygens S, Preijers T, Swaneveld FH, Kleine Budde I, GeurtsvanKessel CH,",
    "Koch BCP, Rijnders BJA. Dosing of Convalescent Plasma and Hyperimmune",
    "Anti-SARS-CoV-2 Immunoglobulins: A Phase I/II Dose-Finding Study.",
    "Clin Pharmacokinet. 2024;63(4):497. doi:10.1007/s40262-024-01351-w.",
    sep = " "
  )
  vignette <- "Huygens_2024_antiSarsCov2IgG"
  # The CLIA (LIAISON SARS-CoV-2 TrimericS IgG) titre is expressed in WHO
  # binding antibody units (BAU). The dose is the product titre (BAU/mL)
  # times the infused volume (mL); V1 and V2 are in mL, so central / vc is
  # directly the serum titre in BAU/mL.
  units <- list(time = "h", dosing = "BAU", concentration = "BAU/mL")

  compartmentData <- list(
    central = list(analyte = "anti-SARS-CoV-2 spike IgG", units = "BAU", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "anti-SARS-CoV-2 spike IgG", units = "BAU", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    LBM = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the central volume, (LBM / 70)^0.538, Huygens 2024 Table 2 footnote equation for V1. The reference value is 70 kg as printed in that equation, not the cohort median that the Methods describe for continuous covariates generically. The paper does not name the lean-body-weight formula and does not tabulate the cohort lean body weight; Table 1 gives median total body weight 82 kg, median height 178 cm and 52.4% male. The base-model V1 (3340 mL) equals 3790 * (55/70)^0.538, consistent with a cohort-typical lean body weight near 55 kg.",
      source_name = "LBW"
    ),
    FORM_SARSCOV2IGG_HYPERIMMUNE = list(
      description = "Antibody product indicator: 1 = hyperimmune anti-SARS-CoV-2 globulin (COVIg), 0 = convalescent plasma (ConvP)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (convalescent plasma)",
      notes = "Multiplicative theta^COV effect on the peripheral volume, V2 = 1640 * 1.99^COVIg, Huygens 2024 Table 2 footnote equation and Table 2 row 'V2 - if using COVIg'. The Table 2 footnote prose and the Results / Discussion text describe the factor as applying to ConvP instead; the printed equation is followed here, and it is the reading that reproduces the base-model V2 of 2700 mL as the cohort geometric mean (30 COVIg and 12 ConvP subjects) and the longer COVIg terminal half-life in supplement S7. See the vignette for the full reconciliation.",
      source_name = "COVIg"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 42L,
    n_studies = 1L,
    age_range = "Median 61 years (IQR 51-64); subgroup ranges 31-74 years",
    weight_range = "Median 82 kg; subgroup ranges 48-145 kg",
    height_range = "Median 178 cm (IQR 167-183)",
    sex_female_pct = 47.6,
    race_ethnicity = "Caucasian 88.1%, African 4.8%, Asian 4.8%",
    disease_state = "Immunocompromised adults without anti-SARS-CoV-2 spike antibodies at baseline despite vaccination (rituximab 61.9%, solid-organ transplant 16.7%, ocrelizumab 9.5%, haematopoietic stem-cell transplant 4.8%, CVID 2.4%, other 4.8%), uninfected at enrolment",
    dose_range = "Single IV infusion: ConvP 300 or 600 mL (neutralising titre 500 or 910 IU/mL; median CLIA titre 3230 or 3070 BAU/mL) or COVIg 150 or 300 mL (910 IU/mL; median CLIA titre 3985 BAU/mL)",
    regions = "Netherlands (Erasmus MC, Rotterdam)",
    n_observations = 326L,
    notes = "Phase I/II open-label non-randomised dose-finding trial NL9379 (April 2021 to April 2022), stopped early when Omicron became dominant. 44 enrolled, 42 analysed (12 ConvP: arms P2 n = 3, P3 n = 3, P4 n = 6; 30 COVIg: arms C3 n = 11, C4 n = 19). Serum CLIA titres at baseline, 24 h, 48 h and 1, 2, 4, 6, 8, 12, 18 and 24 weeks or until seronegative (< 33.8 BAU/mL). Huygens 2024 Table 1, Fig. 1 and supplement S6."
  )

  ini({
    # Final-model estimates, Huygens 2024 Table 2 ('Final model' column) and
    # the parameter equations printed in the Table 2 footnote.
    lcl <- log(9.1); label("Clearance (mL/h)") # Table 2 final CL = 9.1 mL/h
    lvc <- log(3790); label("Central volume of distribution at a lean body weight of 70 kg (mL)") # Table 2 final V1 = 3790 mL
    lq <- log(27.0); label("Inter-compartmental clearance (mL/h)") # Table 2 final Q = 27.0 mL/h
    lvp <- log(1640); label("Peripheral volume of distribution for convalescent plasma (mL)") # Table 2 final V2 = 1640 mL

    e_lbm_vc <- 0.538; label("Power exponent of lean body weight on the central volume (unitless)") # Table 2 'V1 - lean body weight' = 0.538
    e_form_sarscov2igg_hyperimmune_vp <- 1.99; label("Multiplicative factor on the peripheral volume for hyperimmune globulin versus convalescent plasma (ratio)") # Table 2 'V2 - if using COVIg' = 1.99

    # IIV reported as %CV; omega^2 = log(CV^2 + 1). No IIV on Q.
    etalcl ~ 0.17920 # Table 2 final IIV CL 44.3% CV: log(0.443^2 + 1)
    etalvc ~ 0.071882 # Table 2 final IIV V1 27.3% CV: log(0.273^2 + 1)
    etalvp ~ 0.081823 # Table 2 final IIV V2 29.2% CV: log(0.292^2 + 1)

    addSd <- 3.75; label("Additive residual error (BAU/mL)") # Table 2 final additive residual error = 3.75 BAU/mL
    propSd <- 0.0968; label("Proportional residual error (fraction)") # Table 2 final proportional residual error = 9.68% CV
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (LBM / 70)^e_lbm_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * e_form_sarscov2igg_hyperimmune_vp^FORM_SARSCOV2IGG_HYPERIMMUNE

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in BAU, volume in mL -> BAU/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
