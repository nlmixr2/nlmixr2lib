Leber_2023_everolimus_sheep <- function() {
  description <- "Preclinical (sheep). One-compartment population PK model for orally administered everolimus in healthy ewes, with first-order absorption, an absorption lag, and a saturable CytoSorb hemoadsorption clearance arm whose clearance decreases linearly with the drug amount already adsorbed on the cartridge (Leber 2023)."
  reference <- "Leber B, Liebchen U, Rohrhofer L, Weber J, Klaus T, Scheier J, Sucher R, Stiegler P. Pharmacokinetics of immunosuppressive agents during hemoperfusion in a sheep model. Front Med (Lausanne). 2023;10:1258661. doi:10.3389/fmed.2023.1258661"
  vignette <- "Leber_2023_immunosuppressants_hemoperfusion"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    HEMOADSORB_ACTIVE = list(
      description        = "Hemoadsorption-active indicator (1 while a CytoSorb cartridge is in the extracorporeal circuit, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hemoadsorption cartridge in circuit)",
      notes              = "Time-varying within subject. Gates the saturable CytoSorb adsorption clearance arm: cl_hemoadsorption is added to the systemic CL/F only while HEMOADSORB_ACTIVE = 1. In the source study the cartridge ran for a single 6-hour session (Methods, Catheter implantation and extracorporeal circulation); control animals received the identical extracorporeal circuit with no cartridge, i.e. HEMOADSORB_ACTIVE = 0 throughout. Setting it to 0 for the whole record reproduces the sham-circuit / no-device arm. The source paper has no data column for this indicator; it encodes the intervention-versus-control group assignment described in Methods, Drug administration.",
      source_name        = NULL
    ),
    BFR = list(
      description        = "Blood flow rate through the extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Held at approximately 120 mL/min throughout the study (Methods, Catheter implantation and extracorporeal circulation), so it is a study-wide constant rather than a fitted covariate. It enters ONLY the post-adsorber observation equation, never the clearance model: Leber 2023 Methods defines the cross-adsorber extraction denominator as FL = blood flow x (1 - hematocrit), the plasma / serum flow through the adsorber. Carried in the canonical mL/min and converted to L/h inside model() (x 0.06). Meaningful only when HEMOADSORB_ACTIVE = 1.",
      source_name        = "blood flow"
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not measured per animal. Leber 2023 Methods states haematocrit \"was calculated as 33% according to the standard value in sheep\" (reference 25), so 33 is the study-wide constant. Enters ONLY the post-adsorber observation equation, as the plasma fraction 1 - HCT/100 that converts blood flow into the plasma / serum flow FL through the adsorber. Same role as in ButraguenoLaiseca_2025_teicoplanin.",
      source_name        = NULL
    )
  )

  compartmentData <- list(
    depot          = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "everolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    adsorbed       = list(analyte = "everolimus", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "sheep (ewe)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "5 years (all animals)",
    weight_range   = "approximately 85 kg",
    sex_female_pct = 100,
    disease_state  = "Healthy (no induced pathology). Animals were acclimatised for 14 days, then dosed to steady state against clinical target blood levels before a single extracorporeal session.",
    dose_range     = "Everolimus 3-8.25 mg orally twice daily, individually titrated to a whole-blood target of 3-8 ng/mL (Table 1). Given in combination with mycophenolate mofetil 1 g twice daily and prednisolone 10 mg (group 3).",
    regions        = "Single-centre preclinical study, Medical University of Graz, Austria (Austrian Committee for Animal Trials approval 2020-0.437.202).",
    notes          = "Fifteen 5-year-old ewes of approximately 85 kg were studied across five drug-combination groups (Methods, Animals). For each immunosuppressant combination five animals were allocated to the CytoSorb intervention group and three to the sham-circuit control group (Methods, Drug administration), giving n = 8 for the everolimus group. Four additional pilot animals were used to characterise the oral absorption of ciclosporin, tacrolimus, mycophenolate mofetil and everolimus, and their data were included in the analysis. Everolimus was given orally because no intravenous preparation exists (Methods, Drug administration). Blood was sampled from the extracorporeal circuit immediately before and at 30, 90, 250 and 330 min after the start of the adsorber / sham procedure; everolimus was measured in whole blood by LC-MS/MS (Methods, Blood samples and laboratory analysis)."
  )

  ini({
    # Structural parameters. The paper fits the model in two steps: step 1
    # estimates the disposition model ignoring hemoperfusion (Supplementary
    # Table 1); step 2 adds the CytoSorb adsorption pathway with "all
    # parameters of step 1 except clearance ... fixed" (Methods,
    # Pharmacokinetic parameters calculations). CL/F below is therefore the
    # re-estimated step-2 value and every other step-1 parameter is fixed().
    lcl     <- log(130);          label("Apparent oral clearance (CL/F, L/h)")                   # Leber 2023 Supplementary Table 2 (final adsorption model)
    lvc     <- fixed(log(3100));  label("Apparent central volume of distribution (V1/F, L)")     # Leber 2023 Supplementary Table 1 (fixed in the step-2 adsorption model)
    lka     <- fixed(log(0.232)); label("First-order oral absorption rate constant (KA, 1/h)")   # Leber 2023 Supplementary Table 1 (fixed in the step-2 adsorption model)
    ltlag   <- fixed(log(4.14));  label("Absorption lag time (ALAG, h)")                         # Leber 2023 Supplementary Table 1 (fixed in the step-2 adsorption model)

    # CytoSorb saturable adsorption sub-model (Leber 2023 Eq. 2):
    #   CL_CytoSorb(t) = CLmax * (1 - A_CytoSorb(t) / Amax)
    # Both parameters are estimated in step 2.
    lclmax_hemoadsorption <- log(3.23);    label("Maximum CytoSorb hemoadsorption clearance at zero saturation (CLmax, L/h)")  # Leber 2023 Supplementary Table 2
    lamax_hemoadsorption  <- log(0.0163);  label("Maximum drug amount adsorbable by the CytoSorb cartridge (Amax, mg)")        # Leber 2023 Supplementary Table 2

    # Interindividual variability, carried from the step-1 model and held
    # fixed in step 2. Supplementary Table 1 reports CV%; the exponential
    # random-effect model implies omega^2 = log(CV^2 + 1).
    #   CL/F : 16.4% CV -> log(0.164^2 + 1) = 0.0265407
    etalcl ~ fixed(0.0265407)   # Leber 2023 Supplementary Table 1 (IIV CL/F, 16.4% CV)

    # Residual error on the pre-adsorber (inlet) samples, which are the
    # systemic whole-blood concentrations this model predicts.
    propSd <- 0.158; label("Proportional residual error, pre-adsorber samples (fraction)")  # Leber 2023 Supplementary Table 2 (Prop. error pre, 15.8% CV)
    propSd_Cpostfilter <- 0.09; label("Proportional residual error, post-adsorber samples (fraction)")  # Leber 2023 Supplementary Table 2 (Prop. error post, 9.0% CV)
  })
  model({
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    alag_d <- exp(ltlag)

    clmax_hemoadsorption <- exp(lclmax_hemoadsorption)
    amax_hemoadsorption  <- exp(lamax_hemoadsorption)

    # Saturable CytoSorb adsorption clearance (Leber 2023 Eq. 2). The arm is
    # self-limiting: as the cartridge load adsorbed approaches Amax the
    # clearance falls to zero, so the adsorbed amount asymptotes to Amax and
    # never exceeds it. Gated off entirely when no cartridge is in circuit.
    cl_hemoadsorption <- HEMOADSORB_ACTIVE * clmax_hemoadsorption *
      (1 - adsorbed / amax_hemoadsorption)

    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - (cl + cl_hemoadsorption) * central / vc
    d/dt(adsorbed)       <-  cl_hemoadsorption * central / vc

    alag(depot) <- alag_d

    # Plasma / serum flow through the adsorber (Leber 2023 Methods):
    # FL = blood flow x (1 - hematocrit). BFR is carried in the canonical
    # mL/min and converted to L/h (x 0.06); at the study's 120 mL/min and
    # 33% haematocrit this is 4.824 L/h. The floor keeps the ratio finite on
    # a record with the circuit off.
    plasma_flow_filter <- max(BFR * 0.06 * (1 - HCT / 100), 0.001)
    filter_extraction  <- cl_hemoadsorption / plasma_flow_filter

    # Two observables, matching the paired goodness-of-fit panels of
    # Supplementary Figure 9 ("Observations preCyto [mg/L]" and
    # "Observations postCyto [mg/L]"), each against its own individual
    # predictions. Cc is the inlet (systemic) concentration; Cpostfilter is
    # the outlet concentration after the cartridge, obtained by inverting the
    # paper's own cross-adsorber clearance definition CL = (Ci - Co)/Ci x FL
    # (Methods, Determination of clearance and elimination by the adsorber).
    # Dose in mg and vc in L give central / vc directly in mg/L, the unit the
    # Supplementary Figure 9 axes carry.
    Cc          <- central / vc
    Cpostfilter <- Cc * max(1 - filter_extraction, 0.001)

    Cc          ~ prop(propSd)
    Cpostfilter ~ prop(propSd_Cpostfilter)
  })
}
