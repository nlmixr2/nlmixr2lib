Plock_2017_TAK648_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). PK/PD disease-progression model for HbA1c in",
    "female db/db mice treated once daily for 28 days with the",
    "phosphodiesterase-4 (PDE4) inhibitor TAK-648 (1, 3 or 10 mg/kg/day),",
    "roflumilast (0.3, 1 or 3 mg/kg/day) or its active metabolite",
    "roflumilast N-oxide (3 mg/kg/day) (Plock 2017). HbA1c rises linearly",
    "from its pre-dose baseline H0 with slope SLOPE; drug slows the rise",
    "by the factor (1 - EFF), where EFF = EMAX * tPDE4i / (PDE50 + tPDE4i).",
    "The exposure metric tPDE4i (total PDE4 inhibition) is the free",
    "steady-state average plasma concentration of each active compound",
    "divided by its in-vitro PDE4 IC50, summed over roflumilast and",
    "roflumilast N-oxide; it is computed here from the daily dose and the",
    "paper's typical apparent clearances (AUC = Dose / CL, F fixed to 1).",
    "One PDE50 is estimated for TAK-648 and one for the joint",
    "roflumilast + N-oxide exposure, sharing a single random effect. The",
    "mouse compartmental PK models are not reproducible from the paper",
    "(volumes and absorption rates are not printed); only their apparent",
    "clearances, which fully determine the steady-state AUC the PD model",
    "consumes, are used. Time zero is the pre-treatment HbA1c measurement",
    "(study day -1)."
  )
  reference <- paste(
    "Plock N, Vollert S, Mayer M, Hanauer G, Lahu G.",
    "Pharmacokinetic/Pharmacodynamic Modeling of the PDE4 Inhibitor TAK-648",
    "in Type 2 Diabetes: Early Translational Approaches for Human Dose",
    "Prediction. Clin Transl Sci. 2017;10(3):185-193.",
    "doi:10.1111/cts.12436."
  )
  vignette <- "Plock_2017_TAK648_mouse"

  units <- list(
    time = "day",
    dosing = "(mg/kg/day, entered as the DOSE_* covariates; the model has no dose events)",
    concentration = "% (HbA1c, percent of total haemoglobin; the output is not a drug concentration)"
  )

  # PDE50 for TAK-648 and for roflumilast(-N-oxide) carry the same ETA(3)
  # in the published $PRED, so one eta multiplies both.
  paper_specific_etas <- c("etalec50")

  covariateData <- list(
    DOSE_TAK648_MGKG = list(
      description = "Once-daily oral TAK-648 dose of the animal's treatment group",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "0 for animals not receiving TAK-648. Converted to the steady-state AUC over one 24 h dosing interval as 1000 * DOSE / CL (ug*h/L, CL in L/h/kg, F = 1; Plock 2017 Eq. 2).",
      source_name = "dose group (the published $PRED reads a precomputed tPDE4i column TPDE)"
    ),
    DOSE_ROFLUMILAST_MGKG = list(
      description = "Once-daily oral roflumilast dose of the animal's treatment group",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "0 for animals not receiving roflumilast. Produces both roflumilast and roflumilast N-oxide exposure through the apparent clearances estimated after roflumilast dosing.",
      source_name = "dose group (the published $PRED reads a precomputed tPDE4i column TPDE)"
    ),
    DOSE_ROFLUMILASTNOXIDE_MGKG = list(
      description = "Once-daily oral roflumilast N-oxide dose of the animal's treatment group",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "0 for animals not receiving roflumilast N-oxide. Produces both roflumilast N-oxide and (back-converted) roflumilast exposure through the apparent clearances estimated after N-oxide dosing.",
      source_name = "dose group (the published $PRED reads a precomputed tPDE4i column TPDE)"
    )
  )

  population <- list(
    species = "mouse (db/db, BKS.Cg-m+/+Lepr db/BomTac)",
    n_subjects = 80L,
    n_studies = 1L,
    sex_female_pct = 100,
    age_range = "7 weeks at treatment start",
    disease_state = "Leptin-receptor-deficient db/db mouse model of type 2 diabetes",
    dose_range = "Once-daily oral gavage for 28 days: vehicle (4% methocel); TAK-648 1, 3 or 10 mg/kg; roflumilast 0.3, 1 or 3 mg/kg; roflumilast N-oxide 3 mg/kg",
    regions = "Germany",
    notes = "Plock 2017 Supplementary Methods and Supplementary Table 1: 8 PD groups of 10 animals each (80 PD animals) with HbA1c on study days -1, 14 and 28 and two sparse PK samples (day 20 and day 28); 3 satellite PK groups of 6 animals (TAK-648 10 mg/kg, roflumilast 3 mg/kg, roflumilast N-oxide 3 mg/kg) dosed for 34 days with 8 samples each. The PK models were fitted to PK and PD animals jointly; the PK/PD model to the 80 PD animals."
  )

  ini({
    # PD model, Plock 2017 Table 2 (final PK/PD model) and Supplementary
    # PK/PD model code ($PRED)
    lrbase <- log(4.35); label("HbA1c on study day -1, H0 (%)") # Table 2 H0 = 4.35 % (RSE 1.25%); THETA(4)
    lslope <- log(0.138); label("Linear rate of HbA1c increase, SLOPE (%/day)") # Table 2 slope = 0.138 %/DAY (RSE 4.80%); THETA(5)
    lemax <- log(1.00); label("Maximum fractional inhibition of the HbA1c rise, EMAX (unitless)") # Table 2 Emax = 1.00 (RSE 6.27%); THETA(1)
    lec50_tak648 <- log(1.10); label("tPDE4i of TAK-648 giving half-maximal effect, PDE50 TAK-648 (unitless)") # Table 2 PDE50,TAK-648 = 1.10 (RSE 31.0%); THETA(3)
    lec50_roflumilast <- log(27.1); label("Joint roflumilast + roflumilast N-oxide tPDE4i giving half-maximal effect, PDE50 roflumilast(-N-oxide) (unitless)") # Table 2 PDE50,roflumilast(-N-oxide) = 27.1 (RSE 25.2%); THETA(2)

    # Exposure link: apparent clearances of the mouse PK models (Results,
    # 'TAK-648 pharmacokinetic' and 'Pharmacokinetic of roflumilast and
    # roflumilast-N-oxide'), fixed because they come from the separate PK
    # step. Per-kg units follow Eq. 2; 1000 / 2.97 = 337 ug*h/L reproduces
    # Table 1 for TAK-648 1 mg/kg.
    lcl_tak648 <- fixed(log(2.97)); label("Apparent clearance of TAK-648, steady state (L/h/kg)") # Results: 'typical clearance was estimated to be 2.97 L/h'
    lcl_roflumilast <- fixed(log(2.34)); label("Apparent roflumilast clearance after roflumilast dosing (L/h/kg)") # Results: roflumilast CL 2.34 L/h for roflumilast dosing
    lcl_noxide <- fixed(log(3.52)); label("Apparent roflumilast N-oxide clearance after roflumilast dosing (L/h/kg)") # Results: roflumilast-N-oxide CL 3.52 L/h (roflumilast dosing)
    lcl_roflumilast_noxidedose <- fixed(log(41.8)); label("Apparent roflumilast clearance after roflumilast N-oxide dosing (L/h/kg)") # Results: roflumilast CL 41.8 L/h for roflumilast-N-oxide dosing
    lcl_noxide_noxidedose <- fixed(log(1.09)); label("Apparent roflumilast N-oxide clearance after roflumilast N-oxide dosing (L/h/kg)") # Results: roflumilast-N-oxide CL 1.09 L/h (roflumilast-N-oxide dosing)

    # tPDE4i constants, Eq. 3 and the footnote paragraph below Table 1
    fu_tak648 <- fixed(0.200); label("Unbound fraction of TAK-648 in mouse plasma (unitless)") # Results below Table 1: fu 20.0% for TAK-648
    fu_roflumilast <- fixed(0.037); label("Unbound fraction of roflumilast in mouse plasma (unitless)") # Results below Table 1: fu 3.7% for roflumilast
    fu_noxide <- fixed(0.127); label("Unbound fraction of roflumilast N-oxide in mouse plasma (unitless)") # Results below Table 1: fu 12.7% for roflumilast-N-oxide
    ic50_tak648 <- fixed(5.77); label("In-vitro PDE4 IC50 of TAK-648 (ug/L)") # Results below Table 1: IC50 5.77 ug/L for TAK-648
    ic50_roflumilast <- fixed(0.282); label("In-vitro PDE4 IC50 of roflumilast (ug/L)") # Results below Table 1: IC50 0.282 ug/L for roflumilast
    ic50_noxide <- fixed(0.629); label("In-vitro PDE4 IC50 of roflumilast N-oxide (ug/L)") # Results below Table 1: IC50 0.629 ug/L for roflumilast-N-oxide
    tau <- fixed(24); label("Dosing interval (h)") # Eq. 3; once-daily dosing (Methods)

    # IIV: exponential, Table 2 'Between-subject variability' CV%,
    # omega^2 = log(1 + CV^2)
    etalrbase ~ 0.016001 # Table 2 BSV_H0 12.7 CV% (ETA(1))
    etalslope ~ 0.044764 # Table 2 BSV_slope 21.4 CV% (ETA(2))
    etalec50 ~ 0.49277 # Table 2 BSV_PDE50 79.8 CV% (ETA(3), shared by both PDE50 values)

    # Residual error: Y = A1C + A1C*EPS(1)
    propSd <- 0.0387; label("Proportional residual error (fraction)") # Table 2 residual variability proportional 3.87 CV%
  })

  model({
    # Steady-state AUC over one dosing interval, Eq. 2 (ug*h/L; dose in
    # ug/kg = 1000 * mg/kg, CL in L/h/kg, F = 1)
    auc_tak648 <- 1000 * DOSE_TAK648_MGKG / exp(lcl_tak648)
    auc_roflumilast <- 1000 * DOSE_ROFLUMILAST_MGKG / exp(lcl_roflumilast) +
      1000 * DOSE_ROFLUMILASTNOXIDE_MGKG / exp(lcl_roflumilast_noxidedose)
    auc_noxide <- 1000 * DOSE_ROFLUMILAST_MGKG / exp(lcl_noxide) +
      1000 * DOSE_ROFLUMILASTNOXIDE_MGKG / exp(lcl_noxide_noxidedose)

    # Total PDE4 inhibition, Eq. 3; roflumilast and its N-oxide are added
    tpde4i_tak648 <- auc_tak648 * fu_tak648 / (ic50_tak648 * tau)
    tpde4i_roflumilast <- auc_roflumilast * fu_roflumilast / (ic50_roflumilast * tau) +
      auc_noxide * fu_noxide / (ic50_noxide * tau)

    rbase <- exp(lrbase + etalrbase)
    slope <- exp(lslope + etalslope)
    emax <- exp(lemax)
    ec50_tak648 <- exp(lec50_tak648 + etalec50)
    ec50_roflumilast <- exp(lec50_roflumilast + etalec50)

    # Eq. 4. Each animal received one compound, for which
    # EMAX*CE/(PDE50 + CE) equals EMAX*x/(1 + x) with x = CE/PDE50. Writing
    # it through the potency-normalised sum x lets both compound classes
    # share one expression; for a mixture it is an additive extension the
    # paper did not test.
    xrel <- tpde4i_tak648 / ec50_tak648 + tpde4i_roflumilast / ec50_roflumilast
    eff <- emax * xrel / (1 + xrel)

    Hba1c <- rbase + slope * (1 - eff) * t
    Hba1c ~ prop(propSd)
  })
}
