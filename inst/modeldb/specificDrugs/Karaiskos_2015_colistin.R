Karaiskos_2015_colistin <- function() {
  description <- "Population PK model for colistimethate sodium (CMS, prodrug) and colistin (active polymyxin formed by in vivo hydrolysis) in critically ill adults after a 9 MU CMS loading dose. CMS distributes through four compartments representing two states of the prodrug (CMS1 = more fully sulfomethylated, CMS2 = partially sulfomethylated derivatives); each state has central and peripheral compartments sharing volumes Vc and Vp but distinct inter-compartmental clearances Q1 and Q2. The same nonrenal clearance CL_NR drives the first-order hydrolysis CMS1 -> CMS2 (in both central and peripheral, with the same rate constant) and CMS2 -> colistin (central only); CMS1 and CMS2 central compartments are additionally cleared by renal clearance proportional to creatinine clearance. Colistin disposition follows a one-compartment model with apparent clearance and volume (CL/fm, V/fm) scaled to the unknown fraction of administered CMS converted to colistin. Measured colistimethate concentration is the sum of CMS1 and CMS2 central concentrations."
  reference <- paste(
    "Karaiskos I, Friberg LE, Pontikis K, Ioannidis K, Tsagkari V, Galani L,",
    "Kostakou E, Baziaka F, Paskalis C, Koutsoukou A, Giamarellou H.",
    "Colistin population pharmacokinetics after application of a loading dose of",
    "9 MU colistin methanesulfonate in critically ill patients.",
    "Antimicrob Agents Chemother. 2015;59(12):7240-7248.",
    "doi:10.1128/AAC.00554-15.",
    sep = " "
  )
  vignette <- "Karaiskos_2015_colistin"
  paper_specific_compartments <- c("central_cms2", "peripheral1_cms2")

  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance computed from ideal body weight (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate on CMS clearance. Source column CL_CR. Computed by the Cockcroft-Gault equation using ideal body weight (IBW) per Karaiskos 2015 Materials and Methods. Values are capped at 150 mL/min in the analysis dataset; the cap is reproduced in model() via `if (CRCL > 150) CRCL_cap <- 150`. The paper's covariate formula CL_R (L/h) = Sl_CRCL * CR_CL (L/h) requires conversion from mL/min to L/h (multiply by 60/1000). The renal CL effect applies to both CMS1 and CMS2 central compartments with the same slope. Stored under canonical CRCL per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization).",
      source_name        = "CL_CR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 3L,
    age_range      = "18-86 years",
    age_median     = "56.2 years (current study mean)",
    weight_range   = "50-120 kg (current study actual BW)",
    weight_median  = NA_character_,
    sex_female_pct = 42,
    race_ethnicity = "Not reported (Greek ICU population)",
    disease_state  = "Critically ill adults with suspected or microbiologically documented extensively drug-resistant (XDR) Gram-negative infections (Acinetobacter baumannii, Pseudomonas aeruginosa, Klebsiella pneumoniae, Citrobacter spp.); ICU patients NOT on renal replacement therapy. Indications include ventilator-associated pneumonia and tracheobronchitis, bacteremia, complicated intra-abdominal infection, complicated urinary tract infection, catheter-related bloodstream infection, and necrotizing fasciitis. Mean APACHE II score 18.4, mean serum albumin 2.8 g/dL.",
    dose_range     = "9 MU CMS loading dose (approximately 270 mg colistin base activity; approximately 413 umol CMS) over 30 min or 1 h IV infusion, followed by 4.5 MU q12h maintenance commenced 24 h after the loading dose; in patients with creatinine clearance < 60 mL/min the maintenance dose was reduced per the modified Garonzik formula (daily maintenance colistin dose IU = CLCR/10 + 2).",
    regions        = "Greece (Hygeia General Hospital and Sotiria Chest Diseases and General Hospital, Athens)",
    renal_function = "Median creatinine clearance 92.1 mL/min on day 1 (capped at 150 mL/min in the PK analysis); range across the current cohort 29-220 mL/min on day 1.",
    notes          = "The PK analysis pools 19 new ICU patients (current study; 9 MU loading) with 28 patients from two earlier studies of the same group (3 MU q8h and 6 MU loading + 3 MU q8h regimens) for a total of 47 patients and 1144 observed concentrations. Baseline demographics for the current cohort per Karaiskos 2015 Table 1. The packaged model assumes the new-study dosing convention (F1 = 1, all dose into CMS1 central); users reproducing the earlier-study cohorts can split each dose between central and central_cms2 using the source-reported F1 = 0.892 and F = 0.610 fractions (see the validation vignette for the recipe)."
  )

  ini({
    # Structural parameters - Karaiskos 2015 Table 2 final-model column for the
    # simultaneous fit of the new-study and earlier-study data.

    # CMS nonrenal clearance: drives BOTH hydrolysis steps (CMS1 -> CMS2 in
    # central and peripheral, and CMS2 -> colistin in central) at first-order
    # rate constant CL_NR / Vc (the "same rate constant" reported for both
    # central and peripheral conversion). Reported as CL_NR,CMS in the paper.
    lcl_nonren <- log(5.84);    label("CMS nonrenal / hydrolysis clearance CL_NR,CMS (L/h)") # Karaiskos 2015 Table 2: CL_NR,CMS = 5.84 L/h (11% RSE)

    # Renal CL slope on creatinine clearance for CMS1 and CMS2 (same slope
    # applied to both central compartments). The paper formula
    # CL_R (L/h) = Sl_CRCL * CR_CL (L/h) makes Sl_CRCL unitless when CrCL is
    # in L/h. The value also captures CrCL-dependent IIV via a shared eta
    # with lcl_nonren in model() (paper reports identical 16% IIV / 40% IOV on
    # both CL_NR and Sl_CRCL, consistent with a single shared eta on a common
    # CMS-clearance scale factor).
    e_crcl_cl_renal <- 0.541;   label("CMS renal CL slope on CrCL (L/h per L/h; unitless)") # Karaiskos 2015 Table 2: Sl_CRCL = 0.541 (16% RSE); note Table 2 footnote also reports 0.510 when CrCL is not capped at 150 mL/min

    lvc        <- log(1.42);    label("CMS central volume of distribution Vc (L; shared by CMS1 and CMS2)") # Karaiskos 2015 Table 2: V1 = 1.42 L (13% RSE)
    lvp        <- log(12.5);    label("CMS peripheral volume of distribution Vp (L; shared by CMS1 and CMS2)") # Karaiskos 2015 Table 2: V2 = 12.5 L (10% RSE)
    lq         <- log(550);     label("CMS1 inter-compartmental clearance Q1 (L/h)") # Karaiskos 2015 Table 2: Q1 = 550 L/h (31% RSE)
    lq_cms2    <- log(7.75);    label("CMS2 inter-compartmental clearance Q2 (L/h)") # Karaiskos 2015 Table 2: Q2 = 7.75 L/h (11% RSE)

    # Colistin elimination CL/fm (apparent clearance scaled to unknown
    # fraction of CMS converted to colistin). Renal clearance of colistin was
    # not significantly different from zero, so colistin elimination is a
    # single CL term.
    lcl_col    <- log(4.99);    label("Colistin apparent elimination clearance CL/fm (L/h)") # Karaiskos 2015 Table 2: CL_NR,CMS column, colistin = 4.99 L/h (25% RSE)
    lvc_col    <- log(80.4);    label("Colistin apparent central volume of distribution V/fm (L)") # Karaiskos 2015 Table 2: V1 column, colistin = 80.4 L (11% RSE)

    # Inter-individual variability. The source paper applies a SHARED eta to
    # both CL_NR,CMS and Sl_CRCL (Table 2 reports identical 16% IIV / 40% IOV
    # for both parameters, with matching RSE), so model() applies etalcl_nonren
    # to both lcl_nonren AND e_crcl_cl_renal. The colistin clearance IIV in the
    # source is constructed as 4.52 * eta_CMS (perfect correlation), and the
    # 71% CV reported in the table is derived as 4.52 * 16% via the
    # scaling-factor footnote b. This extraction simplifies that to an
    # independent etalcl_col with the same marginal variance (71% CV) and
    # documents the loss of correlation in the vignette Assumptions section.
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl_nonren ~ 0.02528  # log(0.16^2 + 1); 16% CV on CMS CL (shared with Sl_CRCL via model())
    etalcl_col    ~ 0.40546  # log(0.71^2 + 1); 71% CV on colistin CL/fm

    # Residual error. Source applies separate additive + proportional errors
    # to CMS and colistin observations (Table 2). The source ALSO links the
    # two outputs by sharing one component of the residual error (since both
    # are determined from the same sample); this extraction does not
    # reproduce the shared component and treats the residuals as independent,
    # which is documented in the vignette Assumptions section.
    propSd     <- 0.157;  label("CMS proportional residual error (fraction)") # Karaiskos 2015 Table 2 CMS column: proportional error 0.157 (6.7% RSE; 15.7%)
    addSd      <- 0.159;  label("CMS additive residual error (umol/L)")       # Karaiskos 2015 Table 2 CMS column: additive 0.159 umol/L (10% RSE)
    propSd_col <- 0.0884; label("Colistin proportional residual error (fraction)") # Karaiskos 2015 Table 2 colistin column: proportional error 0.0884 (13% RSE; 8.84%)
    addSd_col  <- 0.0629; label("Colistin additive residual error (umol/L)") # Karaiskos 2015 Table 2 colistin column: additive 0.0629 umol/L (14% RSE)
  })

  model({
    # Creatinine-clearance capping (Karaiskos 2015 Results: "creatinine
    # clearance was capped at 150 mL/min" in the PK analysis) and conversion
    # from data-side mL/min to model-side L/h.
    CRCL_cap <- CRCL
    if (CRCL_cap > 150) CRCL_cap <- 150
    crcl_lph <- CRCL_cap * 60 / 1000

    # Individual structural parameters. The shared eta on CL_NR and Sl_CRCL
    # reproduces the Karaiskos 2015 Table 2 finding that both parameters
    # carry the same 16% IIV (and 40% IOV in the source; IOV is omitted in
    # this extraction).
    cl_nonren     <- exp(lcl_nonren + etalcl_nonren)
    sl_crcl_i     <- e_crcl_cl_renal * exp(etalcl_nonren)
    cl_renal_cms  <- sl_crcl_i * crcl_lph
    vc            <- exp(lvc)
    vp            <- exp(lvp)
    q1            <- exp(lq)
    q2            <- exp(lq_cms2)
    cl_col        <- exp(lcl_col + etalcl_col)
    vc_col        <- exp(lvc_col)

    # Micro-constants. The hydrolysis CMS1 -> CMS2 first-order rate constant
    # k_hyd = CL_NR / Vc is applied identically in the central AND peripheral
    # CMS1 compartments per Karaiskos 2015 Results: "conversion from CMS1 to
    # CMS2 also was allowed to occur in the peripheral compartment, with the
    # same rate constant as that in the central compartment." The same k_hyd
    # is applied to the CMS2 central -> colistin step; CMS2 peripheral does
    # not contribute to colistin formation (Figure 2 caption / labels show
    # the CL_NR,CMS/Vc arrow only from CMS2,c to colistin).
    k_hyd        <- cl_nonren / vc
    k_ren_cms    <- cl_renal_cms / vc
    k12_cms1     <- q1 / vc
    k21_cms1     <- q1 / vp
    k12_cms2     <- q2 / vc
    k21_cms2     <- q2 / vp
    k_col        <- cl_col / vc_col

    # ODE system. Compartments and their roles:
    #   central          = CMS1 central (dose entry; F1 = 1 in the current
    #                      study). Cleared renally and by hydrolysis to CMS2c;
    #                      exchanges with CMS1 peripheral via Q1.
    #   peripheral1      = CMS1 peripheral. Hydrolyses to CMS2 peripheral at
    #                      the same rate constant k_hyd.
    #   central_cms2     = CMS2 central. Receives hydrolysis influx from CMS1
    #                      central; cleared renally and by hydrolysis to
    #                      colistin; exchanges with CMS2 peripheral via Q2.
    #   peripheral1_cms2 = CMS2 peripheral. Receives hydrolysis influx from
    #                      CMS1 peripheral; exchanges with CMS2 central via Q2.
    #   central_col      = colistin central. Receives hydrolysis influx from
    #                      CMS2 central; eliminated at rate k_col.
    d/dt(central) <- -k_hyd * central - k_ren_cms * central -
                      k12_cms1 * central + k21_cms1 * peripheral1
    d/dt(peripheral1) <- k12_cms1 * central - k21_cms1 * peripheral1 -
                          k_hyd * peripheral1
    d/dt(central_cms2) <- k_hyd * central - k_hyd * central_cms2 -
                           k_ren_cms * central_cms2 -
                           k12_cms2 * central_cms2 +
                           k21_cms2 * peripheral1_cms2
    d/dt(peripheral1_cms2) <- k_hyd * peripheral1 +
                               k12_cms2 * central_cms2 -
                               k21_cms2 * peripheral1_cms2
    d/dt(central_col) <- k_hyd * central_cms2 - k_col * central_col

    # Observations. CMS measurement is the SUM of CMS1 and CMS2 central
    # concentrations (Karaiskos 2015 Results: "The measured CMS (C_CMS) was
    # the sum of the predicted concentrations in the two central
    # compartments (C_CMS1 plus C_CMS2)").
    Cc     <- (central + central_cms2) / vc
    Cc_col <- central_col / vc_col

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_col ~ add(addSd_col) + prop(propSd_col)
  })
}
