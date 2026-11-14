Hosawi_2026_hydroxyurea <- function() {
  description <- paste(
    "Two-compartment oral population PK model of hydroxyurea in lactating",
    "women (Hosawi 2026, HELPS study NCT02990598; 16 women, single 1000 mg",
    "oral dose) with first-order absorption after a lag time, linear",
    "elimination from the central compartment and a reversible breast-milk",
    "compartment that exchanges with plasma through a single inter-",
    "compartmental clearance (no elimination from milk). Apparent clearance",
    "and central volume are allometrically scaled to body weight (fixed",
    "exponents 0.75 and 1, reference 70 kg). Proportional residual error on",
    "plasma and additive residual error on breast-milk concentrations.",
    sep = " "
  )
  reference <- paste(
    "Hosawi A, Irie K, Dong M, Ware J, Ware R, Mizuno T. Model-Informed",
    "Evaluation of Hydroxyurea Exposure During Lactation. Clin Pharmacol",
    "Ther. 2026. doi:10.1002/cpt.70394. PMCID: PMC13388475.",
    sep = " "
  )
  vignette <- "Hosawi_2026_hydroxyurea"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "hydroxyurea", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "hydroxyurea", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hydroxyurea", units = "mg", specimen = "tissue", verified = TRUE),
    milk = list(analyte = "hydroxyurea", units = "mg", specimen = "milk", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Maternal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of apparent clearance (exponent 0.75) and",
        "central volume (exponent 1), reference 70 kg (Methods 'Covariate",
        "analysis'; Table 1 units 'L/h/70 kg' and 'L/70 kg'). Inclusion",
        "lowered the OFV by 11.9 (Results). Q, V2, the plasma-milk",
        "clearance and the milk volume are reported unscaled in Table 1.",
        sep = " "
      ),
      source_name = "body weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    n_observations = "101 plasma + 93 breast-milk hydroxyurea concentrations",
    sex_female_pct = 100,
    disease_state = paste(
      "Lactating women enrolled in the Hydroxyurea Exposure in Lactation",
      "Study (HELPS, NCT02990598). Their underlying clinical status is not",
      "described in this paper.",
      sep = " "
    ),
    dose_range = "Single oral dose of 1000 mg hydroxyurea (approximately 15-20 mg/kg)",
    weight_range = paste(
      "Not tabulated; a 1000 mg dose at approximately 15-20 mg/kg implies",
      "roughly 50-67 kg",
      sep = " "
    ),
    lactation_stage = "Breastfeeding women; postpartum stage not reported",
    notes = paste(
      "Plasma and milk sampled pre-dose and serially to 8-25 h post-dose",
      "(Methods 'Patients and sampling'); two women contributed milk",
      "samples only. Assayed by validated HPLC. Fit in NONMEM 7.5.0 with",
      "FOCE-I; plasma parameters were estimated first, then the milk",
      "transfer parameters.",
      sep = " "
    )
  )

  ini({
    lka <- log(2.5); label("First-order absorption rate constant, KA (1/h)") # Table 1 'KA (/h)' = 2.5 (RSE 36.3%)
    ltlag <- log(0.2); label("Absorption lag time (h)") # Table 1 'LAG Time (h)' = 0.2 (RSE 14.2%)
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (unitless)") # Table 1 'F' = 1 Fixed
    lcl <- log(8.57); label("Apparent clearance CL/F at 70 kg (L/h)") # Table 1 'CL (L/h/70 kg)' = 8.57 (RSE 8.65%)
    lvc <- log(39.0); label("Apparent central volume V/F at 70 kg (L)") # Table 1 'V(L/70 kg)' = 39.0 (RSE 11.5%)
    lq <- log(4.4); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 1 'Q (L/h)' = 4.4 (RSE 18.3%)
    lvp <- log(19.1); label("Apparent peripheral volume V2/F (L)") # Table 1 'V2 (L)' = 19.1 (RSE 19.8%)
    lq_milk <- log(2.3); label("Apparent plasma-milk inter-compartmental clearance (L/h)") # Table 1 'CL plasma-milk (L/h)' = 2.3 (RSE 17.8%)
    lvmilk <- log(0.79); label("Apparent breast-milk compartment volume, VM (L)") # Table 1 'VM (L)' = 0.79 (RSE 6.80%)

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Methods 'Covariate analysis': theoretical exponent 0.75, standardized to 70 kg
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V (unitless)") # Methods 'Covariate analysis': theoretical exponent 1, standardized to 70 kg

    # IIV: Table 1 reports CV%; Methods Eq. 2 CV% = 100 * sqrt(exp(omega^2) - 1),
    # so omega^2 = log(1 + CV^2). No IIV on Q, V2 or F.
    etalcl ~ 0.0723914 # Table 1 'IIVCL' = 27.4 CV%
    etalvc ~ 0.164606 # Table 1 'IIVV' = 42.3 CV%
    etalka ~ 0.832735 # Table 1 'IIVKA' = 114 CV%
    etaltlag ~ 0.174767 # Table 1 'IIVLAG time' = 43.7 CV%
    etalq_milk ~ 0.23688 # Table 1 'IIVCL plasma-milk' = 51.7 CV%
    etalvmilk ~ 0.0524255 # Table 1 'IIVVM' = 23.2 CV%

    propSd <- 0.15; label("Proportional residual SD, plasma (fraction)") # Table 1 'Prop.Err.(CV%)' = 15.0 (plasma)
    addSd_Cmilk <- 1.12; label("Additive residual SD, breast milk (mg/L)") # Table 1 'Add.Err. (ug/mL)' = 1.12 (milk)
  })

  model({
    # Individual parameters (Methods Eq. 1, exponential IIV; allometric
    # scaling on CL and V only, per the Table 1 units and footnote).
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    q_milk <- exp(lq_milk + etalq_milk)
    vmilk <- exp(lvmilk + etalvmilk)

    Cc <- central / vc
    Cp <- peripheral1 / vp
    Cmilk <- milk / vmilk

    # Figure 1: depot -> central (KA, lag, F); central <-> peripheral (Q);
    # central <-> breast milk (CL plasma-milk); elimination (CL) from central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cc - q * (Cc - Cp) - q_milk * (Cc - Cmilk)
    d/dt(peripheral1) <- q * (Cc - Cp)
    d/dt(milk) <- q_milk * (Cc - Cmilk)

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose mg, volumes L -> concentrations in mg/L (numerically ug/mL, the
    # units of the paper).
    Cc ~ prop(propSd)
    Cmilk ~ add(addSd_Cmilk)
  })
}
