`Leuppi-Taegtmeyer_2019_colistin` <- function() {
  description <- "Six-structural-compartment population PK model for colistimethate sodium (CMS, prodrug) and colistin (Col, active metabolite) during continuous renal replacement therapy (CRRT) in critically ill adults. CMS converts to Col by first-order metabolism (CL1M); Col is cleared by metabolism (CL2M). Both CMS and Col exchange with a CRRT filter compartment (driven by blood flow QBL and patient hematocrit HCT) and a downstream cartridge / effluent compartment (driven by effluent flow QEFF and species-specific sieving coefficients SC_CMS and SC_COL). Filter and cartridge priming volumes are fixed device constants (V_filter = 0.2 L, V_cart = 0.3 L). Bioavailability factor F = 1155.5/1749.8 on the dose compartment converts mg of CMS sodium administered into mg colistin-base equivalents (the molar-mass ratio of colistin to CMS sodium); central concentrations are reported as mg/L colistin-base equivalents."
  reference <- paste(
    "Leuppi-Taegtmeyer A, Decosterd LA, Osthoff M, Mueller NJ, Buclin T, Corti N.",
    "Multicenter Population Pharmacokinetic Study of Colistimethate Sodium and",
    "Colistin Dosed as in Normal Renal Function in Patients on Continuous Renal",
    "Replacement Therapy. Antimicrob Agents Chemother. 2019;63(2):e01957-18.",
    "doi:10.1128/AAC.01957-18.",
    "DDMORE Foundation Model Repository: DDMODEL00000295.",
    sep = " "
  )
  vignette <- "Leuppi-Taegtmeyer_2019_colistin"
  paper_specific_compartments <- c("filter", "filter_col", "cart", "cart_col", "effl", "effl_col")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  ddmore_id <- "DDMODEL00000295"
  replicate_of <- NULL

  covariateData <- list(
    HCT = list(
      description        = "Hematocrit, expressed as a fraction (0-1) inside the model. Determines the plasma-fraction scaling on the CRRT filter exchange rates: K13 / K31 / K24 / K42 each scale with (1 - HCT). Source column HT in the DDMORE bundle is encoded as a fraction; canonical HCT in inst/references/covariate-columns.md is documented in % units (0-100). When using HCT-in-% data, divide by 100 in model() OR supply the column already as a fraction (0-1).",
      units              = "fraction (0-1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The .mod treats HT == 0 as missing and substitutes a default of 0.25 (the population-typical hematocrit assumed for critically ill CRRT patients). The same fallback is reproduced in model() with `if (HCT <= 0) HCT_use <- 0.25`. Source column name in the DDMORE bundle: HT.",
      source_name        = "HT"
    ),
    QBL = list(
      description        = "CRRT blood-flow rate read from the dialysis machine setting at each event. Drives the CMS / Col mass-transfer rates K13, K31, K24, K42 between systemic circulation and the CRRT filter. Source units in the DDMORE bundle are mL/min; the model converts to L/h via the constant 60/1000.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying CRRT operational parameter. Not in the canonical covariate register (inst/references/covariate-columns.md); CRRT-specific. Source column name: QBL.",
      source_name        = "QBL"
    ),
    QEFF = list(
      description        = "CRRT effluent flow rate read from the dialysis machine setting at each event. Drives the rates K35, K46 (filter -> cartridge) and K50, K60 (cartridge -> effluent reservoir). Source units in the DDMORE bundle are mL/h; the model converts to L/h via the constant 1/1000.",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying CRRT operational parameter. Not in the canonical covariate register; CRRT-specific. Source column name: QEFF.",
      source_name        = "QEFF"
    )
  )

  population <- list(
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "adults (critically ill, ICU); detailed age statistics not in the DDMORE bundle and the linked publication is not on disk",
    weight_range   = "adult body-weight range; not in the DDMORE bundle",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Critically ill adult ICU patients receiving continuous renal replacement therapy (CRRT) for acute kidney injury, treated with colistimethate sodium (CMS) for systemic Gram-negative infection.",
    dose_range     = "Loading dose 9 MIU CMS (~720 mg of CMS sodium / ~300 mg colistin equivalents) IV, then 3 MIU q8h maintenance (~240 mg CMS sodium / ~100 mg colistin equivalents per dose). Dosing follows the standard normal-renal-function regimen (the study question was whether CRRT removal mandates dose adjustment).",
    regions        = "Multicenter Switzerland (Basel, Lausanne, Zurich) per the DDMORE-bundle Model_Accommodations file; specific site mix not enumerated.",
    crrt_modality  = "Continuous renal replacement therapy. The model abstracts the CRRT circuit as a series of two compartments (filter and cartridge) with priming volumes V_filter = 0.2 L and V_cart = 0.3 L, and operational flow rates QBL (blood) and QEFF (effluent) supplied per event. Sieving coefficients SC_CMS and SC_COL govern membrane permeability for the two analytes.",
    notes          = "n_subjects = 10 per the DDMORE bundle's Model_Accommodations.txt and the .rdf model-has-description-long field ('CMS and colistin pharmacokinetics were assessed prospectively in 10 critically ill patients requiring CRRT. Extensive pharmacokinetic sampling was performed on treatment day 1, 3 and 5 after administration of a loading dose of CMS, followed by a maintenance dosage every eight hours.'). The linked publication (Leuppi-Taegtmeyer 2019, AAC) is not present on the operator's disk, so detailed baseline demographics (age, weight, sex, study sites, comorbidities) could not be read from Table 1; populate these fields if the paper PDF becomes available."
  )

  ini({
    # Final estimates from Output_real_CMS_colistin_PK_CRRT.lst (DDMORE
    # bundle DDMODEL00000295), FINAL PARAMETER ESTIMATE block at lines
    # 650-661 (THETA), 665-680 (OMEGA diagonals). MINIMIZATION SUCCESSFUL
    # at OBJV = 300.337 (lines 576, 645). The .mod $THETA initial values
    # at .mod lines 116-123 happen to match these estimates because the
    # DDMORE re-run was a posthoc EVALUATION (MAXEVAL=1 ... POSTHOC) at
    # the previously fitted point estimates.

    # CMS structural PK
    lcl       <- log(2.31);  label("CMS metabolic clearance to colistin CL1M (L/h)")        # .lst TH 1 = 2.31E+00; .mod L116
    lvc       <- log(12.1);  label("CMS systemic distribution volume V1 (L)")               # .lst TH 2 = 1.21E+01; .mod L117

    # CMS sieving coefficient (dimensionless, bounded [0, 1.2] per .mod L118)
    sc_cms    <- c(0, 1.05, 1.2); label("CMS sieving coefficient SC_CMS through CRRT filter (dimensionless, 0..1.2)") # .lst TH 3 = 1.05E+00; .mod L118

    # Colistin structural PK
    lcl_col   <- log(1.93);  label("Colistin metabolic-elimination clearance CL2M (L/h)")   # .lst TH 4 = 1.93E+00; .mod L119
    lvc_col   <- log(70.1);  label("Colistin systemic distribution volume V2 (L)")          # .lst TH 5 = 7.01E+01; .mod L120

    # Colistin sieving coefficient (dimensionless, bounded [0, 1.2] per .mod L121)
    sc_col    <- c(0, 0.454, 1.2); label("Colistin sieving coefficient SC_COL through CRRT filter (dimensionless, 0..1.2)") # .lst TH 6 = 4.54E-01; .mod L121

    # Inter-individual variability (diagonal $OMEGA, .mod L125-129).
    # NONMEM stores omega^2 on the internal log scale because the .mod
    # parameterizes CL1M = TVCL1M * EXP(ETA(1)) etc.
    etalcl     ~ 0.270; label("IIV variance on log CMS CL1M (~55% CV)")           # .lst OMEGA(1,1) = 2.70E-01; .mod L126
    etalvc     ~ 0.133; label("IIV variance on log CMS V1 (~37% CV)")             # .lst OMEGA(2,2) = 1.33E-01; .mod L127
    etalcl_col ~ 0.304; label("IIV variance on log colistin CL2M (~59% CV)")      # .lst OMEGA(3,3) = 3.04E-01; .mod L128
    etalvc_col ~ 0.247; label("IIV variance on log colistin V2 (~52% CV)")        # .lst OMEGA(4,4) = 2.47E-01; .mod L129

    # Residual error. NONMEM .mod $ERROR (lines 98-103) implements a single
    # combined error model shared across both observation compartments
    # (CMS at CMT=1 and colistin at CMT=2):
    #     W = SQRT(THETA(7)^2 * IPRED^2 + THETA(8)^2)
    #     Y = IPRED + W * EPS(1)
    # with SIGMA(1,1) FIX 1, so THETA(7) and THETA(8) carry the per-record
    # proportional and additive residual SDs respectively. Here those two
    # values are duplicated to give per-output residual-error parameters
    # (Cc for CMS, Cc_col for colistin) per the nlmixr2lib multi-output
    # convention; the duplication is a naming-only choice and reproduces
    # the source's shared-error behaviour exactly.
    propSd     <- 0.222; label("CMS proportional residual SD (fraction)")              # .lst TH 7 = 2.22E-01; .mod L122
    addSd      <- 0.459; label("CMS additive residual SD (mg/L colistin equivalent)")  # .lst TH 8 = 4.59E-01; .mod L123
    propSd_col <- 0.222; label("Colistin proportional residual SD (fraction; shared with CMS in source)")    # source uses shared THETA(7); see ini comment above
    addSd_col  <- 0.459; label("Colistin additive residual SD (mg/L colistin equivalent; shared with CMS in source)") # source uses shared THETA(8); see ini comment above
  })

  model({
    # =====================================================================
    # Verbatim translation of Executable_CMS_colistin_PK_CRRT.mod
    # ($PK lines 28-83; $DES lines 86-96; $ERROR lines 98-101). The model
    # has 10 NONMEM compartments; the last two (CMSAUC, COLAUC) are
    # cumulative-AUC bookkeeping integrators that are not needed at the
    # nlmixr2lib level (PKNCA in the validation vignette computes AUC).
    # The 8 structural compartments are kept.
    # =====================================================================

    # ----- Hematocrit fallback: HT == 0 in the source dataset is treated
    # as a missing-value sentinel and replaced by the typical CRRT-cohort
    # value 0.25 (.mod L35-41 / L53-59). The data-side encoding in the
    # DDMORE simulated CSV represents missing as "." (mapped to 0 by
    # NONMEM's `IGNORE=#` reader); we reproduce the 0-treated-as-missing
    # convention here.
    HCT_use <- HCT
    if (HCT_use <= 0) HCT_use <- 0.25

    # ----- Individual structural parameters (CL ~ log-normal IIV)
    cl     <- exp(lcl + etalcl)            # CMS -> Col conversion CL (L/h)
    vc     <- exp(lvc + etalvc)            # CMS systemic volume (L)
    cl_col <- exp(lcl_col + etalcl_col)    # Colistin elimination CL (L/h)
    vc_col <- exp(lvc_col + etalvc_col)    # Colistin systemic volume (L)

    # ----- CRRT device constants (priming volumes; .mod L34, L43, L52, L61).
    v_filter <- 0.2   # CMS filter compartment volume (L); also used for Col filter (.mod V3 = V4 = 0.2)
    v_cart   <- 0.3   # CMS cartridge compartment volume (L); also used for Col cartridge (.mod V5 = V6 = 0.3)

    # ----- Operational flow conversions (.mod L36-40, L54-58, L44-45, L62-63)
    # QBL is mL/min in the source dataset; convert to L/h via *60/1000.
    # QEFF is mL/h in the source dataset; convert to L/h via /1000.
    qbl_lph  <- QBL  * 60 / 1000   # blood flow (L/h)
    qeff_lph <- QEFF        / 1000 # effluent flow (L/h)

    # ----- CMS rate constants (CRRT exchange and elimination)
    # CMS to filter and back (.mod L36-37, L39-40):
    k13 <- qbl_lph * (1 - HCT_use) / vc        # CMS central -> CMS filter (1/h)
    k31 <- qbl_lph * (1 - HCT_use) / v_filter  # CMS filter -> CMS central (1/h)
    # CMS to colistin metabolic conversion (.mod L33):
    k12 <- cl / vc                              # CMS central -> Col central (1/h)
    # CMS filter -> CMS cartridge (sieving-modulated; .mod L44):
    k35 <- sc_cms * qeff_lph / v_filter        # CMS filter -> CMS cartridge (1/h)
    # CMS cartridge -> effluent waste (.mod L45):
    k50 <- qeff_lph / v_cart                    # CMS cartridge -> effluent (1/h)

    # ----- Col rate constants (CRRT exchange and elimination)
    # Col to filter and back (.mod L54-55, L57-58):
    k24 <- qbl_lph * (1 - HCT_use) / vc_col    # Col central -> Col filter (1/h)
    k42 <- qbl_lph * (1 - HCT_use) / v_filter  # Col filter -> Col central (1/h)
    # Col metabolic elimination (.mod L51):
    k20 <- cl_col / vc_col                      # Col central -> elimination (1/h)
    # Col filter -> Col cartridge (.mod L62):
    k46 <- sc_col * qeff_lph / v_filter        # Col filter -> Col cartridge (1/h)
    # Col cartridge -> effluent waste (.mod L63):
    k60 <- qeff_lph / v_cart                    # Col cartridge -> effluent (1/h)

    # =====================================================================
    # ODE system (.mod $DES L86-96; AUC integrators DADT(9), DADT(10) dropped)
    # =====================================================================
    d/dt(central)     <- -k12 * central - k13 * central + k31 * filter                              # CMS central (.mod DADT(1))
    d/dt(central_col) <- -k24 * central_col - k20 * central_col + k12 * central + k42 * filter_col   # Col central (.mod DADT(2))
    d/dt(filter)      <-  k13 * central - k31 * filter - k35 * filter                                 # CMS filter (.mod DADT(3))
    d/dt(filter_col)  <-  k24 * central_col - k42 * filter_col - k46 * filter_col                     # Col filter (.mod DADT(4))
    d/dt(cart)        <-  k35 * filter - k50 * cart                                                   # CMS cartridge (.mod DADT(5))
    d/dt(cart_col)    <-  k46 * filter_col - k60 * cart_col                                           # Col cartridge (.mod DADT(6))
    d/dt(effl)        <-  k50 * cart                                                                  # cumulative CMS effluent (.mod DADT(7))
    d/dt(effl_col)    <-  k60 * cart_col                                                              # cumulative Col effluent (.mod DADT(8))

    # =====================================================================
    # Bioavailability factor on the CMS dose compartment.
    # F1 = 1155.5/1749.8 (.mod L76) is the molar-mass ratio of colistin to
    # colistimethate sodium; it converts the administered AMT (mg of CMS
    # sodium) into mg colistin-base equivalents in the modelled CMS
    # central compartment, so all internal concentrations are reported as
    # mg/L colistin-base equivalents (the standard reporting convention
    # for colistin / CMS PK).
    # =====================================================================
    f(central) <- 1155.5 / 1749.8

    # =====================================================================
    # Observations and residual error (.mod $ERROR L98-101).
    # C1 = A(1)/V1 (CMS) and C2 = A(2)/V2 (Col) are the predicted
    # plasma concentrations at the two NONMEM observation compartments
    # CMT = 1 and CMT = 2. Both share the same residual error structure
    # in the source; we declare per-output names per nlmixr2lib
    # multi-output convention while initialising to the same values.
    # =====================================================================
    Cc     <- central     / vc
    Cc_col <- central_col / vc_col
    Cc     ~ add(addSd)     + prop(propSd)
    Cc_col ~ add(addSd_col) + prop(propSd_col)
  })
}
