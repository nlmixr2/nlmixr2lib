Leshinsky_2017_caspofungin_cat <- function() {
  description <- "Preclinical (cat). Two-compartment population PK model with first-order linear elimination from the central compartment for intravenous caspofungin acetate in healthy adult cats (Leshinsky 2017)"
  reference <- "Leshinsky J, McLachlan A, Foster DJR, Norris R, Barrs VR. Pharmacokinetics of caspofungin acetate to guide optimal dosing in cats. PLoS ONE. 2017;12(6):e0178783. doi:10.1371/journal.pone.0178783"
  vignette <- "Leshinsky_2017_caspofungin_cat"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with reference weight 4 kg (cohort mean per Leshinsky 2017 Methods 'Animals'; Table 1 footnote 'CL * (WT/4)^0.75'). Fixed exponents 0.75 on CL and Q, 1 on V1 and V2.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "cat (domestic shorthair)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "mean 3.94 +/- 1.84 years (healthy adult)",
    weight_range   = "mean 4.4 +/- 0.56 kg",
    sex_female_pct = 50.0,
    disease_state  = "Healthy desexed domestic shorthair cats from a research colony.",
    dose_range     = "1 mg/kg IV infusion over 1 h; single dose (n=8) and once-daily multiple dose for 7 days (n=6, 3M/3F).",
    regions        = "Australia (Sydney; Eurofins SCEC research colony).",
    scope_note     = "Preclinical-only popPK fit. Filed under inst/modeldb/pharmacokinetics/ rather than specificDrugs/ because nlmixr2lib's specificDrugs tier is reserved for human drugs.",
    notes          = "Baseline demographics per Leshinsky 2017 Methods 'Animals' (page on Materials and methods). 8 adult desexed domestic shorthair cats (4F/4M), mean age 3.94 +/- 1.84 years and mean body weight 4.4 +/- 0.56 kg. The multi-dose phase used a subset of 6 cats (3M/3F)."
  )

  ini({
    # Structural parameters. Leshinsky 2017 Table 1 (final FOCE estimates).
    # Reference weight WT = 4 kg (cohort mean per Methods 'Animals' and
    # Table 1 footnote). Concentration in the central compartment is
    # Cc = central / vc, with dose in mg and volumes in L -> Cc in mg/L
    # (= ug/mL).
    lcl <- log(0.0175); label("Clearance for the reference 4 kg cat (L/h)")          # Leshinsky 2017 Table 1 (CL = 17.5 mL/h, %RSE 7)
    lvc <- log(0.214);  label("Central volume of distribution for the reference 4 kg cat (L)")   # Leshinsky 2017 Table 1 (V1 = 214 mL, %RSE 8.6)
    lvp <- log(0.143);  label("Peripheral volume of distribution for the reference 4 kg cat (L)") # Leshinsky 2017 Table 1 (V2 = 143 mL, %RSE 8.3)
    lq  <- log(0.150);  label("Intercompartmental clearance Q for the reference 4 kg cat (L/h)") # Leshinsky 2017 Table 1 (Q = 150 mL/h, %RSE 20)

    # Allometric exponents fixed by the authors (Leshinsky 2017 Methods 'General
    # modeling strategy': 'exponent fixed of 0.75 for clearance parameters and 1
    # for volumes', and Table 1 row labels 'CL * (WT/4)^0.75', 'V1 * (WT/4)',
    # 'V2 * (WT/4)', 'Q * (WT/4)^0.75'). Reference weight 4 kg.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of WT on CL and Q (unitless; fixed)") # Leshinsky 2017 Methods / Table 1 footnote
    e_wt_vc_vp <- fixed(1.0);  label("Allometric exponent of WT on V1 and V2 (unitless; fixed)") # Leshinsky 2017 Methods / Table 1 footnote

    # IIV. Leshinsky 2017 Table 1 reports an exponential PPV model (eq 2,
    # Methods 'General modeling strategy'). Only CL had an independent eta;
    # V1's PPV was encoded as a scale factor applied to the CL random effect
    # ('scale factor for the PPV on V1 as a function of the PPV on CL', Results
    # 'Pharmacokinetic modeling'). Equivalent NONMEM-style coding:
    #   ETA_CL = ETA(1);  ETA_V1 = scale_etalvc * ETA(1)
    # Convert 18% CV on CL via omega^2 = log(CV^2 + 1) = log(1 + 0.18^2)
    #         = 0.03194.
    etalcl ~ 0.03194    # Leshinsky 2017 Table 1 (BSV on CL = 18% CV, 0% shrinkage, %RSE 21)
    scale_etalvc <- 1.05; label("Scale factor: V1 BSV = scale_etalvc * CL random effect (unitless)") # Leshinsky 2017 Table 1 (V/F BSV scale factor = 1.05, %RSE 8.4)

    # Residual error: combined proportional + additive (Leshinsky 2017 Methods
    # 'General modeling strategy' eq 3 and Results / Table 1).
    propSd <- 0.078; label("Proportional residual error (fraction)") # Leshinsky 2017 Table 1 (proportional %CV = 7.8, %RSE 18)
    addSd  <- 0.61;  label("Additive residual error SD (ug/mL)")     # Leshinsky 2017 Table 1 (additive SD = 0.61 ug/mL, %RSE 20)
  })

  model({
    # Individual PK parameters with fixed allometric WT scaling on
    # CL / Q (exponent 0.75) and V1 / V2 (exponent 1.0), centred at the
    # reference 4 kg cat. V1's PPV is the CL random effect scaled by
    # scale_etalvc (preserves the high CL-V1 correlation reported in
    # Leshinsky 2017 Results 'Pharmacokinetic modeling' without estimating an
    # independent eta on V1).
    cl <- exp(lcl + etalcl)                * (WT / 4)^e_wt_cl_q
    vc <- exp(lvc + scale_etalvc * etalcl) * (WT / 4)^e_wt_vc_vp
    q  <- exp(lq)                          * (WT / 4)^e_wt_cl_q
    vp <- exp(lvp)                         * (WT / 4)^e_wt_vc_vp

    # Two-compartment IV model with first-order linear elimination from the
    # central compartment (Leshinsky 2017 Results 'Pharmacokinetic modeling':
    # 'The two-compartment linear model provided the best description'). Dose
    # is delivered directly into central via IV infusion; no depot or
    # absorption process.
    Cc <- linCmt()
    Cc ~ prop(propSd) + add(addSd)
  })
}
