Fournier_2018_amoxicillin <- function() {
  description <- paste(
    "Two-compartment IV population PK model for amoxicillin in adult ICU",
    "burn patients hospitalized at a Swiss tertiary-care centre, with",
    "Cockcroft-Gault creatinine clearance as a linear covariate on CL",
    "(centered at 110 mL/min) and body weight as a linear (allometric",
    "exponent 1) covariate on the central volume V1 (centered at 70 kg)",
    "(Fournier 2018)."
  )
  reference <- paste(
    "Fournier A, Goutelle S, Que Y-A, Eggimann P, Pantet O, Sadeghipour F,",
    "Voirol P, Csajka C. Population pharmacokinetic study of",
    "amoxicillin-treated burn patients hospitalized at a Swiss",
    "tertiary-care center. Antimicrob Agents Chemother. 2018;62(9):e00505-18.",
    "doi:10.1128/AAC.00505-18."
  )
  vignette <- "Fournier_2018_amoxicillin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Cockcroft-Gault creatinine clearance (raw mL/min, NOT",
        "BSA-normalized)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CL_CR. Computed by the Cockcroft-Gault equation in",
        "raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under",
        "the canonical CRCL column per inst/references/covariate-columns.md",
        "(CRCL accepts raw mL/min when the source paper does not apply BSA",
        "normalization; document the assay form per model). Reference value",
        "110 mL/min (population average CRCL, Fournier 2018 Table 3",
        "footnote). The effect is applied to CL as a divisively-centered",
        "(median-normalized) term: TV_CL = exp(lcl) * (1 + e_crcl_cl *",
        "(CRCL/110 - 1)). The literal Table 3 footnote",
        "TV_CL = CL * [1 + theta * (CRCL - 110)] with theta = 0.57 yields",
        "physically absurd CL across the simulated CRCL range (negative",
        "below ~108 mL/min, 711 L/h at CRCL = 200); the divisively-centered",
        "interpretation matches the paper's discussion-level CL ratios",
        "across the 15-200 mL/min range and follows the Delattre 2010",
        "amikacin precedent. Cohort median CRCL was 128 mL/min (IQR",
        "65-150) with augmented CRCL (>150 mL/min) in 25% of patients",
        "(Fournier 2018 Table 1 and Results paragraph 1)."
      ),
      source_name        = "CL_CR"
    ),
    WT = list(
      description        = paste(
        "Actual body weight at the time of amoxicillin sampling."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BW. Time-varying within the amoxicillin therapy",
        "window in principle, but Fournier 2018 reports limited",
        "intraindividual change (median -2.3%; minimum -11.4%, maximum",
        "+10.6%) over the course of therapy (Results paragraph 1).",
        "Reference 70 kg (Fournier 2018 Table 3 footnote: TVV1 = V1 *",
        "(BW/70)). Applied as a linear (allometric exponent 1) scalar on",
        "V1 only. Body weight on admission ranged from 60 to 132 kg",
        "(Results paragraph 1; Table 1 median 72.4 kg, IQR 67.0-83.6)."
      ),
      source_name        = "BW"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (binary; 1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Cohort 16 male / 5 female (Fournier 2018 Table 1: 76.2% male).",
        "Screened in the stepwise covariate analysis (Methods 'Population",
        "PK model building'); not retained in the final model."
      )
    ),
    BSA_BURN = list(
      description = "Total burnt body surface area (percentage of TBSA).",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Cohort median TBSA 23% (IQR 12.5-44; Fournier 2018 Table 1).",
        "Screened, not retained in the final model."
      )
    ),
    ALBUMIN = list(
      description = "Serum albumin concentration.",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened in the stepwise covariate analysis (Methods 'Population",
        "PK model building'); not retained in the final model. Concentration",
        "value not tabulated in the paper."
      )
    ),
    CREAT = list(
      description = "Serum creatinine concentration.",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened in the stepwise covariate analysis (Methods 'Population",
        "PK model building'); not retained in the final model (CRCL was",
        "retained instead)."
      )
    ),
    WT_ADM = list(
      description = "Body weight on admission to the burn ICU.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Median admission BW 72.4 kg (IQR 67.0-83.6; range 60-132;",
        "Fournier 2018 Table 1 and Results paragraph 1). Screened",
        "alongside time-varying actual BW (used as WT); not retained as",
        "an independent covariate."
      )
    ),
    WT_GAIN = list(
      description = paste(
        "Body weight gain since admission, defined as max(0, BW - BW_ADM)."
      ),
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Hinge-style covariate tested in the forward selection (Methods",
        "'Population PK model building'); not retained in the final model."
      )
    ),
    WT_LOSS = list(
      description = paste(
        "Body weight loss since admission, defined as max(0, BW_ADM - BW)."
      ),
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Hinge-style covariate tested in the forward selection (Methods",
        "'Population PK model building'); not retained in the final model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "16-93 years (mean 50.1, SD 24.3)",
    weight_range   = "60-132 kg admission (median 72.4, IQR 67.0-83.6)",
    sex_female_pct = 23.8,
    disease_state  = paste(
      "Adult patients with severe burns hospitalized at the Burn Centre",
      "of a Swiss tertiary-care intensive care unit (Centre Hospitalier",
      "Universitaire Vaudois, Lausanne). All patients received a course",
      "of intravenous amoxicillin (alone or in combination with",
      "clavulanic acid) for treatment of an infection during the first",
      "weeks of hospitalization. Median total burnt body surface area",
      "23% (IQR 12.5-44); 76.2% had inhalation lesions; mean SAPS II",
      "35.9 (SD 18.6); median Ryan score 1 (IQR 1-2); median ICU length",
      "of stay 23 days (IQR 13.0-39.5); ICU mortality 9.5%."
    ),
    renal_function = paste(
      "Median Cockcroft-Gault creatinine clearance 128 mL/min (IQR",
      "65-150); 25% of patients had CRCL > 150 mL/min (augmented renal",
      "clearance)."
    ),
    dose_range     = paste(
      "Intravenous amoxicillin 1 to 2 g every 6 to 8 hours per",
      "manufacturer recommendations in patients with normal renal",
      "function; infusion duration 30 min for the first dose and 2 h",
      "(amoxicillin) or 1 h (amoxicillin-clavulanic acid) starting from",
      "the second dose per local guidelines. Dosage adjusted for renal",
      "insufficiency (eGFR < 30 mL/min: 500 mg-2 g q8-12h; eGFR < 15",
      "mL/min: 750 mg-2 g q24h)."
    ),
    regions        = "Switzerland (single centre, CHUV Lausanne Burn Centre)",
    enrollment     = "Prospective consecutive enrollment, October 2013-October 2016",
    trial_id       = "ClinicalTrials.gov NCT01965340",
    notes          = paste(
      "Baseline demographics per Fournier 2018 Table 1. 185 amoxicillin",
      "plasma concentrations from 21 burn patients; a rich kinetic",
      "profile (samples at 0, 1, 2, 3, 4, 5 h after end of infusion) was",
      "obtained on one occasion for 18 of 21 patients (1 patient had two",
      "rich profiles); the remaining samples were trough (days 2, 4, 6,",
      "8) and random (days 6, 8) levels. Plasma free fraction assumed",
      "82% for amoxicillin (Methods 'Dosing simulations', citing",
      "reference 57)."
    )
  )

  ini({
    # Structural disposition parameters at the reference covariate state
    # (CRCL = 110 mL/min, BW = 70 kg). Source: Fournier 2018 Table 3
    # "Covariate model mean estimate (RSE, %)" column.
    lcl <- log(13.6); label("Amoxicillin clearance at CRCL = 110 mL/min (L/h)")           # Table 3 row CL
    lvc <- log(9.73); label("Central volume V1 at BW = 70 kg (L)")                         # Table 3 row V1
    lvp <- log(17.6); label("Peripheral volume V2 (L)")                                    # Table 3 row V2
    lq  <- log(20.1); label("Intercompartmental clearance Q (L/h)")                        # Table 3 row Q

    # Fractional covariate effect of Cockcroft-Gault CRCL on CL, centered
    # at the population average CRCL of 110 mL/min:
    #   TV_CL = exp(lcl) * (1 + e_crcl_cl * (CRCL/110 - 1))
    # Source: Fournier 2018 Table 3 row "theta_CLCR_CL" (covariate model
    # column) and Table 3 footnote. NOTE: The literal Table 3 footnote
    # equation TVCL = CL * [1 + theta * (CRCL - 110)] yields negative or
    # physically absurd CL across most of the CRCL range with theta = 0.57
    # (e.g. -374 L/h at CRCL = 60 and 711 L/h at CRCL = 200). The
    # divisively-centered (median-normalized) form below interprets the
    # CRCL term as the relative change in CRCL from the population mean and
    # yields plausible CL values across the entire 15-200 mL/min span the
    # paper simulated (19.94 L/h at CRCL = 200; 6.91 L/h at CRCL = 15).
    # See the vignette "Assumptions and deviations" section for the full
    # reconciliation. The policy precedent is Delattre 2010 amikacin
    # (inst/modeldb/specificDrugs/Delattre_2010_amikacin.R), where the same
    # absurd-as-printed centered-linear equation was operator-resolved to
    # the divisively-centered form.
    e_crcl_cl <- 0.57; label("Fractional CL slope per (CRCL/110 - 1)")                     # Table 3 row theta_CLCR_CL

    # Inter-individual variability on CL only. Fournier 2018 Table 3 row
    # "Random effect omega CL (% CV)" reports omega_CL = 37.3% CV under
    # an exponential (lognormal) eta model (Methods 'Population PK model
    # building'). omega^2 = log(CV^2 + 1) = log(0.373^2 + 1) = 0.13028.
    # IIV was not retained on V1, V2, or Q (Results paragraph 3: "the
    # others had fixed, estimated values").
    etalcl ~ 0.13028   # Table 3 row omega CL (37.3% CV converted to log-normal variance)

    # Residual error: combined additive + proportional model (Results
    # paragraph 3, "best described by a combined additive and proportional
    # residual error model").
    propSd <- 0.37;  label("Proportional residual error (fraction)")                       # Table 3 row "Proportional error (%)"
    addSd  <- 0.08;  label("Additive residual error (mg/L)")                               # Table 3 row "Additive error (mg/L)"
  })

  model({
    # Individual CL: typical value with the divisively-centered
    # (median-normalized) CRCL effect, then multiplied by lognormal IIV.
    # See ini() comment on e_crcl_cl for the interpretation rationale.
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL / 110 - 1))

    # Individual V1: typical value with linear (allometric exponent 1)
    # body-weight scaling, no IIV.
    vc <- exp(lvc) * (WT / 70)

    # V2 and Q: typical values only, no IIV.
    vp <- exp(lvp)
    q  <- exp(lq)

    # Micro-constants for the two-compartment model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: IV infusion enters the central compartment directly
    # (no depot; infusion rate is set on the dose record).
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # Observation. Dose in mg / vc in L -> mg/L, matching the paper's
    # observed amoxicillin concentrations (Cmax up to ~200 mg/L,
    # Fig. 2 caption).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
