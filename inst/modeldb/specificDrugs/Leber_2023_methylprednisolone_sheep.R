Leber_2023_methylprednisolone_sheep <- function() {
  description <- "Preclinical (sheep). Two-compartment population PK model for intravenous methylprednisolone in healthy ewes, with linear elimination, combined proportional and additive residual error, and a saturable CytoSorb hemoadsorption clearance arm whose clearance decreases linearly with the drug amount already adsorbed on the cartridge (Leber 2023)."
  reference <- "Leber B, Liebchen U, Rohrhofer L, Weber J, Klaus T, Scheier J, Sucher R, Stiegler P. Pharmacokinetics of immunosuppressive agents during hemoperfusion in a sheep model. Front Med (Lausanne). 2023;10:1258661. doi:10.3389/fmed.2023.1258661"
  vignette <- "Leber_2023_immunosuppressants_hemoperfusion"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    HEMOADSORB_ACTIVE = list(
      description        = "Hemoadsorption-active indicator (1 while a CytoSorb cartridge is in the extracorporeal circuit, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hemoadsorption cartridge in circuit)",
      notes              = "Time-varying within subject. Gates the saturable CytoSorb adsorption clearance arm: cl_hemoadsorption is added to the systemic CL only while HEMOADSORB_ACTIVE = 1. In the source study the cartridge ran for a single 6-hour session (Methods, Catheter implantation and extracorporeal circulation); the two methylprednisolone control animals received the identical extracorporeal circuit with no cartridge, i.e. HEMOADSORB_ACTIVE = 0 throughout. Setting it to 0 for the whole record reproduces the sham-circuit / no-device arm. The source paper has no data column for this indicator; it encodes the intervention-versus-control group assignment described in Methods, Drug administration.",
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
    central        = list(analyte = "methylprednisolone", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1    = list(analyte = "methylprednisolone", units = "mg", specimen = "serum", verified = TRUE),
    adsorbed    = list(analyte = "methylprednisolone", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "sheep (ewe)",
    n_subjects     = 6L,
    n_studies      = 1L,
    age_range      = "5 years (all animals)",
    weight_range   = "approximately 85 kg",
    sex_female_pct = 100,
    disease_state  = "Healthy (no induced pathology). The methylprednisolone arm simulated intraoperative corticosteroid administration during cardiopulmonary bypass with concurrent hemoadsorption (a heart-transplant scenario).",
    dose_range     = "Methylprednisolone intravenously (Table 1 and Methods, Drug administration): two intervention sheep received 1 g followed by a second 1 g dose 1.5 h later; two intervention sheep received a single increased 1.5 g starting dose; two control sheep received 2 x 1 g without an adsorber in the circuit.",
    regions        = "Single-centre preclinical study, Medical University of Graz, Austria (Austrian Committee for Animal Trials approval 2020-0.437.202).",
    notes          = "Methylprednisolone was the only drug given intravenously and was studied alone (group 5) rather than in an immunosuppressant combination. Six ewes contributed: four in the CytoSorb intervention group across two dosing scenarios and two sham-circuit controls (Table 1; Methods, Drug administration; Figure 3). Samples were taken at 10, 30, 60, 90, 120 and 180 min after the first methylprednisolone dose, from the extracorporeal circuit before (inlet) and after (outlet) the adsorber. Methylprednisolone was quantified by UPLC-MS/MS at University Hospital Reims, France (Methods, Blood samples and laboratory analysis). NOTE: Table 1 describes the control animals as receiving 2 x 1 g while the Methods prose states they received 1 g; the Table is used here. See vignette Assumptions and deviations."
  )

  ini({
    # Structural parameters. The paper fits the model in two steps: step 1
    # estimates the disposition model ignoring hemoperfusion (Supplementary
    # Table 1); step 2 adds the CytoSorb adsorption pathway with "all
    # parameters of step 1 except clearance ... fixed" (Methods,
    # Pharmacokinetic parameters calculations). CL below is therefore the
    # re-estimated step-2 value and every other step-1 parameter is fixed().
    #
    # Methylprednisolone was given intravenously, so these are true CL and V
    # rather than the apparent CL/F and V/F of the orally dosed drugs.
    lcl   <- log(75.1);         label("Systemic clearance (CL, L/h)")                              # Leber 2023 Supplementary Table 2 (final adsorption model)
    lvc   <- fixed(log(21.0));  label("Central volume of distribution (V1, L)")                    # Leber 2023 Supplementary Table 1 (fixed in the step-2 adsorption model)
    lvp   <- fixed(log(34.9));  label("Peripheral volume of distribution (V2, L)")                 # Leber 2023 Supplementary Table 1 (fixed in the step-2 adsorption model)
    lq    <- fixed(log(31.8));  label("Intercompartmental clearance (Q, L/h)")                     # Leber 2023 Supplementary Table 1 (fixed in the step-2 adsorption model)

    # CytoSorb saturable adsorption sub-model (Leber 2023 Eq. 2):
    #   CL_CytoSorb(t) = CLmax * (1 - A_CytoSorb(t) / Amax)
    # Both parameters are estimated in step 2. The paper flags that this CLmax
    # exceeds the 7.2 L/h extracorporeal blood flow, although the confidence
    # interval (6.16-10.25 L/h) covers it (Results). See vignette Assumptions
    # and deviations.
    lclmax_hemoadsorption <- log(8.21);  label("Maximum CytoSorb hemoadsorption clearance at zero saturation (CLmax, L/h)")  # Leber 2023 Supplementary Table 2
    lamax_hemoadsorption  <- log(53.4);  label("Maximum drug amount adsorbable by the CytoSorb cartridge (Amax, mg)")        # Leber 2023 Supplementary Table 2

    # Interindividual variability, carried from the step-1 model and held
    # fixed in step 2. Supplementary Table 1 reports CV%; the exponential
    # random-effect model implies omega^2 = log(CV^2 + 1).
    #   CL : 28.1% CV -> log(0.281^2 + 1) = 0.0759985
    etalcl ~ fixed(0.0759985)   # Leber 2023 Supplementary Table 1 (IIV CL, 28.1% CV)

    # Combined proportional plus additive residual error on the pre-adsorber
    # (inlet) samples, which are the systemic serum concentrations this model
    # predicts. Supplementary Table 2 labels the additive term "[mg]"; it is
    # read here as mg/L, the concentration unit the model predicts, because an
    # additive residual error must carry the units of the dependent variable.
    propSd <- 0.323;  label("Proportional residual error, pre-adsorber samples (fraction)")  # Leber 2023 Supplementary Table 2 (Prop. error pre, 32.3% CV)
    addSd  <- 0.0042; label("Additive residual error, pre-adsorber samples (mg/L)")          # Leber 2023 Supplementary Table 2 (Add. error pre)
    propSd_Cpostfilter <- 0.262;  label("Proportional residual error, post-adsorber samples (fraction)")  # Leber 2023 Supplementary Table 2 (Prop. error post, 26.2% CV)
    addSd_Cpostfilter  <- 0.129;  label("Additive residual error, post-adsorber samples (mg/L)")          # Leber 2023 Supplementary Table 2 (Add. error post)
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    clmax_hemoadsorption <- exp(lclmax_hemoadsorption)
    amax_hemoadsorption  <- exp(lamax_hemoadsorption)

    # Saturable CytoSorb adsorption clearance (Leber 2023 Eq. 2). The arm is
    # self-limiting: as the cartridge load adsorbed approaches Amax the
    # clearance falls to zero, so the adsorbed amount asymptotes to Amax and
    # never exceeds it. Gated off entirely when no cartridge is in circuit.
    cl_hemoadsorption <- HEMOADSORB_ACTIVE * clmax_hemoadsorption *
      (1 - adsorbed / amax_hemoadsorption)

    kel <- (cl + cl_hemoadsorption) / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)        <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)    <-  k12 * central - k21 * peripheral1
    d/dt(adsorbed)    <-  cl_hemoadsorption * central / vc

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
    #
    # NOTE: methylprednisolone is the one drug whose CLmax (8.21 L/h) exceeds
    # the plasma flow FL (4.824 L/h), so the extraction ratio starts above 1
    # and the 0.001 floor BINDS at the onset of hemoperfusion, pinning the
    # outlet at 0.1% of the inlet until the cartridge saturates enough to
    # bring the ratio below 1. The paper raises the same point itself ("the
    # estimated CLmax of MP was higher than the blood flow, but the
    # confidence interval included the value of 7.2 L/h") and Supplementary
    # Figure 4 shows the MP outlet dropping by roughly an order of magnitude
    # at the first on-cartridge sample. See the vignette Errata.
    #
    # Dose in mg and vc in L give central / vc directly in mg/L, the unit the
    # Supplementary Figure 9 axes carry.
    Cc          <- central / vc
    Cpostfilter <- Cc * max(1 - filter_extraction, 0.001)

    Cc          ~ prop(propSd) + add(addSd)
    Cpostfilter ~ prop(propSd_Cpostfilter) + add(addSd_Cpostfilter)
  })
}
