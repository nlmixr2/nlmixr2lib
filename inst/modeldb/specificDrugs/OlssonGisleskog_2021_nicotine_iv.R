OlssonGisleskog_2021_nicotine_iv <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous nicotine in",
    "healthy adult smokers (Olsson Gisleskog 2021; 80 subjects, 4 studies,",
    "0.028 mg/kg infused over 10 min). Clearance and inter-compartmental",
    "flows are allometrically scaled by (WT/70)^0.75 and the three volumes",
    "by (WT/70)^1. Residual nicotine from pre-study smoking is described by",
    "a virtual 1-mg bolus into the central compartment at the start of the",
    "36-h washout, whose bioavailability (the pre-washout nicotine dose,",
    "4.90 mg) carries its own IIV; the record is flagged with",
    "VIRTUAL_DOSE = 1. This IV model supplies the disposition parameters",
    "(and their IIV) that the paper's oral, buccal and transdermal models",
    "hold fixed.",
    sep = " "
  )
  reference <- paste(
    "Olsson Gisleskog PO, Perez Ruixo JJ, Westin A, Hansson AC, Soons PA.",
    "Nicotine Population Pharmacokinetics in Healthy Smokers After",
    "Intravenous, Oral, Buccal and Transdermal Administration.",
    "Clin Pharmacokinet. 2021;60(4):541-561.",
    "doi:10.1007/s40262-020-00960-5",
    sep = " "
  )
  vignette <- "OlssonGisleskog_2021_nicotine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Olsson Gisleskog 2021 Sect. 2.2: nicotine concentrations were assayed in
  # plasma by GC with nitrogen-sensitive or chemical-ionisation MS detection.
  compartmentData <- list(
    central = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling to a 70-kg reference: exponent 0.75 on CL, Q2",
        "and Q3 and 1 on V1, V2 and V3 (Sect. 2.3.1; control stream",
        "WTCL/WTQ = (WT/70)**0.75, WTV = (WT/70)**1)."
      ),
      source_name = "WT"
    ),
    VIRTUAL_DOSE = list(
      description = "Virtual pre-washout smoking-dose record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (a real nicotine dose)",
      notes = paste(
        "Set to 1 on the dose record that represents residual nicotine",
        "from pre-study smoking: a virtual bolus of amt = 1 into central",
        "at the start of the washout (36 h before the IV infusion in the",
        "IV studies), whose bioavailability is the estimated pre-washout",
        "nicotine dose in mg. 0 on every other record, including the IV",
        "infusion. Needed only in this IV model because the real dose also",
        "enters central; the extravascular sibling models dose central",
        "with the virtual bolus only, so they need no flag. Omit the",
        "virtual dose (and set the flag to 0) to simulate a nicotine-naive",
        "subject."
      ),
      source_name = "PREDOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 4L,
    age_range = "20-76 years",
    age_median = "36 years",
    weight_range = "49.0-99.0 kg",
    weight_median = "72.9 kg",
    sex_female_pct = 46.2,
    race_ethnicity = c(White = 97.5, Asian = 2.5),
    disease_state = "healthy adult smokers",
    dose_range = "0.028 mg/kg nicotine intravenous infusion over 10 min",
    regions = "Sweden (Universities of Lund, Linkoping and Gothenburg; studies 1993-2012)",
    notes = paste(
      "Olsson Gisleskog 2021 Table 1 (IV row) and Table 2 (IV column). The",
      "IV dataset excluded subjects with hepatic or renal impairment",
      "(control stream IGNORE=(HEPAT.GT.0) and IGNORE=(RENAL.GT.0)). 173",
      "observations (12%) below the limit of quantification were handled",
      "with the M3 method (LAPLACIAN estimation). The paper pooled 930",
      "subjects across 29 studies and seven formulations; this file is the",
      "IV subset."
    )
  )

  ini({
    # Disposition (Table 3). The full-precision values are the ones the
    # oral, buccal and transdermal control streams (ESM model code) carry
    # as $THETA ... FIX, which match Table 3 to its printed 3 significant
    # figures.
    lcl <- log(67.4136)
    label("Clearance CL for a 70-kg subject (L/h)") # Table 3 'CL (L/h)' 67.4; ESM $THETA 67.4136
    lvc <- log(117.373)
    label("Central volume V1 for a 70-kg subject (L)") # Table 3 'V1 (L)' 117; ESM $THETA 117.373
    lq <- log(38.615)
    label("Inter-compartmental flow Q2 to peripheral1 for a 70-kg subject (L/h)") # Table 3 'Q2 (L/h)' 38.6; ESM $THETA 38.615
    lvp <- log(130.372)
    label("Peripheral volume V2 for a 70-kg subject (L)") # Table 3 'V2 (L)' 130; ESM $THETA 130.372
    lq2 <- log(216.29)
    label("Inter-compartmental flow Q3 to peripheral2 for a 70-kg subject (L/h)") # Table 3 'Q3 (L/h)' 216; ESM $THETA 216.29
    lvp2 <- log(53.4189)
    label("Peripheral volume V3 for a 70-kg subject (L)") # Table 3 'V3 (L)' 53.4; ESM $THETA 53.4189

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on CL, Q2 and Q3 (unitless)") # Sect. 2.3.1 'allometric exponent for CL and intercompartmental flows was fixed to 0.75'
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on V1, V2 and V3 (unitless)") # Sect. 2.3.1 'and to 1.0 for Vn'

    lfcentral <- log(4.90)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 3 'Pre-washout nicotine dose (mg)' 4.90; ESM output TH7 4.90E+00

    # IIV (Table 3 CV% = 100*sqrt(exp(omega^2)-1); variances from the ESM
    # output and the downstream control streams' $OMEGA ... FIX blocks).
    etalcl ~ 0.0705245 # Table 3 'IIV CL' 27.0% CV; ESM $OMEGA IIV_CL 0.0705245
    # Block over V1, V3, V2 in the control-stream order ETA(2), ETA(3),
    # ETA(4): V1 is central (vc), V3 the 53.4-L peripheral (vp2), V2 the
    # 130-L peripheral (vp). The V1-V2 covariance is 0 in the block.
    etalvc + etalvp2 + etalvp ~ c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    ) # Table 3 IIV V1 68.1%, IIV V3 75.4%, IIV V2 230%, covariances V1/V3 -0.231 and V2/V3 0.553; ESM $OMEGA BLOCK(3)
    etalfcentral ~ 0.519 # Table 3 'IIV pre washout nicotine dose' 82.5% CV; ESM output OMEGA(5,5) 5.19E-01

    # Residual error: W = SQRT(addRUV**2 + propRUV**2*IPRED**2) (ESM $ERROR).
    propSd <- 0.0926
    label("Proportional residual error (fraction)") # Table 3 'Proportional residual error' 0.0926
    addSd <- 0.212
    label("Additive residual error (ng/mL)") # Table 3 'Additive residual error (ng/mL)' 0.212
  })
  model({
    # 1. Allometric scaling to 70 kg (ESM $PK WTCL, WTV, WTQ).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 2. Individual parameters. Q2 and Q3 carry no IIV (ESM $PK Q2 = TVQ2,
    #    Q3 = TVQ3).
    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. Three-compartment disposition (ESM ADVAN11 TRANS4).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Residual pre-study nicotine: F1 = THETA(7)*EXP(ETA(5)) on the
    #    PREDOSE record, F1 = 1 on the IV infusion (ESM $PK).
    fcentral <- exp(lfcentral + etalfcentral)
    f(central) <- (1 - VIRTUAL_DOSE) + VIRTUAL_DOSE * fcentral

    # 5. Observation: amounts in mg, volumes in L; S1 = V1/1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
