Savic_2017_cladribine <- function() {
  description <- "Population PK model for cladribine (CdA) in patients with relapsing-remitting multiple sclerosis (Savic 2017): three-compartment disposition with first-order oral absorption, separate fasted vs fed (or unknown food-state) absorption parameters, renal clearance proportional to Cockcroft-Gault creatinine clearance, and a multiplicative non-renal-clearance effect of concomitant subcutaneous interferon beta-1a coadministration."
  reference <- "Savic RM, Novakovic AM, Ekblom M, Munafo A, Karlsson MO. Population Pharmacokinetics of Cladribine in Patients with Multiple Sclerosis. Clin Pharmacokinet. 2017 Oct;56(10):1245-1253. doi:10.1007/s40262-017-0516-6. PMID: 28255848; PMCID: PMC5591346."
  vignette <- "Savic_2017_cladribine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation in raw mL/min (not BSA-normalized). Stored under canonical CRCL per inst/references/covariate-columns.md (CRCL accepts raw Cockcroft-Gault mL/min when the source paper does not BSA-normalize). Reference value 105.2 mL/min corresponds to the paper's typical patient CLCR = 6.31 L/h. Linear effect on renal clearance: cl_renal_typ * CRCL / 105.2.",
      source_name        = "CLCR"
    ),
    FED = list(
      description        = "Fed (or unknown food-state) vs fasted at oral dose administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Source paper Table 4 reports a single set of absorption parameters for the 'unknown/fed state' (Ka = 1.03 1/h, F = 0.4) distinct from the fasted state (Ka = 1.08 1/h, F = 0.456). The 'unknown' category covers dose records in studies 26486 and 25643 that did not record food state; these records share parameters with the fed state. For the packaged model, FED = 1 covers both confirmed-fed and unknown-food-state records.",
      source_name        = "FED"
    ),
    CONMED_IFNB1A = list(
      description        = "Concomitant subcutaneous interferon beta-1a (Rebif) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant IFN beta-1a)",
      notes              = "Time-varying per subject; Savic 2017 study 26486 alternated between cladribine monotherapy and cladribine + IFN beta-1a periods. Multiplicative effect on non-renal clearance: cl_nonrenal *= (1 + e_ifn_clnr * CONMED_IFNB1A), with e_ifn_clnr = 0.21 (21% increase in non-renal CL during IFN beta-1a coadministration).",
      source_name        = "IFNB1A"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 173L,
    n_studies      = 4L,
    age_range      = "19-65 years",
    age_median     = "40 years",
    weight_range   = "48.5-116.1 kg",
    weight_median  = "69.2 kg",
    sex_female_pct = 65.9,
    race_ethnicity = "Not reported in detail; multinational phase I and phase III studies in patients with relapsing-remitting multiple sclerosis",
    disease_state  = "Relapsing-remitting multiple sclerosis",
    dose_range     = "Phase I: 3 mg IV infusion over 1 h, 10 mg single oral tablet, or 1.75 mg/kg orally over 8 weeks +/- subcutaneous IFN beta-1a (8.8/22/44 ug three times weekly). Phase III (CLARITY, study 25643): cumulative oral 3.5 or 5.25 mg/kg over 2 years (one or two 10 mg tablets daily for 4-5 days in weeks 1 and 5 of each year, +/- additional dosing weeks 9 and 13 of year 1).",
    regions        = "Multinational (CLARITY phase III; NCT00213135)",
    renal_function = "CLCR median 107.9 mL/min, range 49.6-244.4 mL/min (Cockcroft-Gault, raw)",
    co_medication  = "Phase I drug-drug interaction study (26486) coadministered subcutaneous IFN beta-1a (Rebif); 16 of 173 patients received IFN beta-1a coadministration",
    notes          = "Pooled analysis of four clinical studies: 25803 (phase I IV+oral crossover, n=16), 26127 (phase I food-effect crossover, n=16), 26486 (phase I cladribine+IFN beta-1a DDI, n=16), and 25643 (CLARITY phase III, n=125). Demographics from Table 3 of the source paper. The original dataset contained 4790 records; 45% were excluded (mostly BLQ, 96.5% of exclusions), leaving 2619 plasma and urine concentration records for population analysis."
  )

  # Notes on implementation choices (see vignette 'Assumptions and deviations'
  # for full justification):
  # * Parent-only model: the source paper jointly fits cladribine and its
  #   metabolite 2-chloroadenine (CAde). The packaged model implements the
  #   cladribine parent kinetics only; CAde compartment is omitted because
  #   Table 4 reports CAde apparent CL/F_m and V/F_m corrected for the
  #   (unidentifiable) fraction metabolized, so CAde concentrations cannot be
  #   simulated on an absolute scale without auxiliary information.
  # * Single shared BSV on Q3, Q4, V3, V4: Table 4 reports one variance
  #   (0.0365) for these four distribution parameters jointly, consistent with
  #   a single ETA applied multiplicatively to all four in the source NONMEM
  #   control stream. Implemented here as one shared 'etadist' added to each
  #   of the four log-transformed distribution parameters.
  # * Logit-scale BSV on bioavailability F: the source paper notes 'Variance
  #   on a logit scale' for BSV(F) = 0.223. Implemented as 'etalogitf' added
  #   inside the logit transform so F stays bounded in (0, 1).
  # * Fed-state absorption delay: the source paper uses a transit-compartment
  #   chain (N = 2.24, MTT = 0.910 h) only in the fed state. The packaged
  #   model approximates this with a lag time alag(depot) = MTT * FED. The
  #   approximation captures the absorption-delay magnitude but loses the
  #   gamma-distributed input-rate shape; the source paper notes the transit
  #   model improved fit over a pure lag time, so for downstream simulation
  #   needing the precise absorption curve, users should consider a
  #   transit() override.
  # * Phase III-specific 0.319 h lag time omitted as a study-specific
  #   artefact that is not generalizable to other simulation populations.
  # * Single proportional residual error magnitude (35.3% from study 25643,
  #   phase III, largest cohort) chosen over the per-study/per-route values
  #   reported in Table 4. The between-subject variability in residual error
  #   (0.159) reported in Table 4 is also omitted from the packaged model.
  ini({
    # Absorption parameters (Savic 2017 Table 4)
    lka           <- log(1.08); label("Absorption rate constant, fasted (1/h)")               # Table 4: Ka = 1.08 1/h (RSE 21.14%)
    e_fed_ka      <- log(1.03 / 1.08); label("FED effect on log-Ka (unitless)")              # Table 4: Ka (unknown/fed state) = 1.03 1/h vs fasted 1.08 1/h
    etalka        ~ 0.102                                                                   # Table 4: BSV(Ka) = 0.102 variance

    # Bioavailability on logit scale (Savic 2017 Table 4 footnote d:
    # "Variance on a logit scale")
    logitf        <- log(0.456 / (1 - 0.456)); label("Logit of oral bioavailability, fasted (unitless)") # Table 4: F = 0.456 (RSE 7.03%)
    e_fed_logitf  <- log(0.4 / (1 - 0.4)) - log(0.456 / (1 - 0.456)); label("FED effect on logit-F (unitless)") # Table 4: F (unknown/fed) = 0.4 vs fasted 0.456
    etalogitf     ~ 0.223                                                                   # Table 4: BSV(F) = 0.223 variance on logit scale (footnote d)

    # Fed-state absorption-delay lag time (transit-compartment approximation;
    # see header notes)
    lmtt_fed      <- log(0.910); label("Lag time for fed-state absorption (h)")              # Table 4: Mean transit time (fed state) = 0.910 h (RSE 11.03%)

    # Renal clearance: CLR proportional to CLCR. CLR (L/h) = CLR_typ * CRCL / 105.2
    # where CLR_typ = 22.2 L/h is the typical CLR at the paper's reference
    # patient CLCR = 6.31 L/h = 105.2 mL/min. Source paper expresses this as
    # CLR (L/h) = coefficient * CLCR (L/h) with a dimensionless coefficient =
    # 3.52; either parameterisation gives the same fit.
    lcl_renal_typ <- log(22.2); label("Renal clearance at typical CRCL (L/h)")               # Table 4: CLR coefficient = 3.52 dimensionless; CLR at CLCR=6.31 L/h = 22.2 L/h (RSE 9.26%)

    # Non-renal clearance (Savic 2017 Table 4)
    lcl_nonrenal  <- log(23.4); label("Non-renal clearance, monotherapy (L/h)")              # Table 4: CLNR = 23.4 L/h (RSE 9.58%)
    e_ifn_clnr    <- 0.21; label("Fractional increase in CLNR with concomitant IFN beta-1a") # Table 4: IFN beta-1a fold-increase in CLNR = 1.21
    etalcl_nonrenal ~ 0.00574                                                                # Table 4: BSV(CLNR) = 0.00574 variance

    # Central volume (Savic 2017 Table 4)
    lvc           <- log(44.0); label("Central volume of distribution V1 (L)")               # Table 4: Central volume = 44.0 L (RSE 22.77%)
    etalvc        ~ 0.209                                                                    # Table 4: BSV(V) = 0.209 variance

    # Distribution (intercompartmental) parameters (Savic 2017 Table 4); a
    # single shared eta governs BSV for Q3, Q4, V3, V4 (Table 4 reports one
    # BSV row "Q3,Q4,V3,V4" = 0.0365). Paper terminology Q3/V3 = central <-> peripheral1,
    # Q4/V4 = central <-> peripheral2; mapped to canonical lq/lvp and lq2/lvp2.
    lq            <- log(14.3); label("Intercompartmental clearance to peripheral1 (L/h)")    # Table 4: Q3 = 14.3 L/h (RSE 7.73%)
    lq2           <- log(53.7); label("Intercompartmental clearance to peripheral2 (L/h)")    # Table 4: Q4 = 53.7 L/h (RSE 19.06%)
    lvp           <- log(347);  label("Peripheral1 volume of distribution (L)")               # Table 4: V3 = 347 L (RSE 6.07%)
    lvp2          <- log(89.5); label("Peripheral2 volume of distribution (L)")               # Table 4: V4 = 89.5 L (RSE 7.97%)
    # Shared eta on the distribution parameters Q, Q2, Vp, Vp2; named etalq
    # for convention matching (Table 4 reports one BSV row covering all four).
    etalq         ~ 0.0365                                                                    # Table 4: BSV(Q3,Q4,V3,V4) = 0.0365 variance (single shared eta on q/q2/vp/vp2)

    # Residual error (proportional; Savic 2017 Table 4)
    propSd        <- 0.353; label("Proportional residual error, fraction")                    # Table 4: RUV plasma (oral; study 25643) = 35.3% (phase III, largest cohort)
  })
  model({
    # Individual absorption parameters
    ka          <- exp(lka + e_fed_ka * FED + etalka)
    logitf_ind  <- logitf + e_fed_logitf * FED + etalogitf
    f_oral      <- 1 / (1 + exp(-logitf_ind))
    tlag_depot  <- exp(lmtt_fed) * FED

    # Individual clearances
    cl_renal     <- exp(lcl_renal_typ) * (CRCL / 105.2)
    cl_nonrenal  <- exp(lcl_nonrenal + etalcl_nonrenal) * (1 + e_ifn_clnr * CONMED_IFNB1A)
    cl           <- cl_renal + cl_nonrenal

    # Individual disposition parameters (shared distribution eta etalq)
    vc  <- exp(lvc  + etalvc)
    q   <- exp(lq   + etalq)
    q2  <- exp(lq2  + etalq)
    vp  <- exp(lvp  + etalq)
    vp2 <- exp(lvp2 + etalq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3-compartment disposition with first-order oral absorption
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central - k13 * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot)    <- f_oral
    alag(depot) <- tlag_depot

    # Cladribine plasma concentration. central (mg) / vc (L) = mg/L =
    # ug/mL; multiply by 1000 for ng/mL to match the source paper's
    # bioanalytical units.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
