Patel_2025_radamts13 <- function() {
  description <- "Two-compartment population PK model of plasma ADAMTS13 activity in patients with congenital thrombotic thrombocytopenic purpura (cTTP) receiving intravenous recombinant ADAMTS13 (rADAMTS13; TAK-755) or plasma-based therapy (PBT). Zero-order IV infusion into the central compartment with first-order linear elimination. Body weight is carried a priori as a power function centred on the 68.7 kg median with fixed allometric exponents (0.75 on clearances, 1.0 on volumes). The only retained covariate effect is treatment type, entered as a relative ADAMTS13 activity (relative bioavailability) multiplier on the delivered dose: 1 for rADAMTS13 (reference), 39.0% lower for plasma / solvent-detergent-treated plasma, and 93.3% lower for plasma-derived FVIII:VWF concentrates. Pooled analysis of three rADAMTS13 trials (65 patients, 2,462 quantifiable activity samples). This model supplies the exposure metric that drives the four companion exposure-response models from the same paper."
  reference <- paste(
    "Patel M, Xu H, Barriere O, Diderichsen P, Patwari P, Zhu AZX,",
    "Marier JF, Peyret T, Wang LT, Mellgard B, Wang W, Bhattacharya I.",
    "Use of PopPK and E-R Analyses toward Explaining Causal Link Between",
    "ADAMTS13 in Recombinant vs. Plasma-Based Therapies and Clinical",
    "Effects in cTTP. Clin Pharmacol Ther. 2025;118(4):813-822.",
    "doi:10.1002/cpt.3720",
    sep = " "
  )
  vignette <- "Patel_2025_radamts13_exposure_response"

  # ADAMTS13 activity is a functional (enzyme-activity) readout dosed and
  # reported in international units. With dose amounts in IU and volumes in L
  # the model's natural concentration unit is IU/L, which is also the unit the
  # source paper uses for the additive residual error (Table 1, "Additive
  # (IU/L): 79.9"). The paper reports every activity summary (Cmax, Cave,
  # EC50) in IU/mL, where 1 IU/mL = 1000 IU/L = 100% of normal ADAMTS13
  # activity. The four companion exposure-response models take their CAV
  # covariate in IU/mL, so a solve of this model must be divided by 1000
  # before it is used to drive them; the validation vignette does this
  # explicitly.
  units <- list(time = "h", dosing = "IU", concentration = "IU/L")

  compartmentData <- list(
    central     = list(analyte = "ADAMTS13 activity", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ADAMTS13 activity", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in this extraction. Carried a priori (not by covariate selection) as a power function centred on the 68.7 kg population median body weight, with allometric exponents fixed at 0.75 on clearance and intercompartmental clearance and 1.0 on both volumes (paper Methods, 'body weight was incorporated a priori in the model using a power function centered to the median body weight (68.7 kg). Fixed allometric exponents were used'; Table 1 rows 'x (WT/68.7)^0.75' and 'x (WT/68.7)^1.0'). Body weight also sets the administered amount, because both rADAMTS13 (40 IU/kg) and PBT (10 IU/kg) are dosed per kilogram. Observed range 18.3-130.0 kg (Table 2).",
      source_name        = "WT"
    ),
    TRT_PBT = list(
      description        = "Plasma / solvent-detergent-treated-plasma treatment-arm indicator (1 = the infusion on this dose record is fresh frozen plasma, solvent/detergent-treated plasma or an equivalent plasma infusion; 0 = it is not).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = rADAMTS13 (the reference treatment, whose relative ADAMTS13 activity is fixed to 1)",
      notes              = "Required input. Per-dose-record covariate, not a subject-level one: the pivotal phase III trial is a crossover in which each patient receives rADAMTS13 in one period and PBT in another, so the indicator changes within a subject. Mutually exclusive with TRT_PDFVIII_VWF; both are 0 on an rADAMTS13 record. Enters as a relative-ADAMTS13-activity (relative bioavailability) multiplier on the delivered dose, not as an effect on clearance or volume.",
      source_name        = "PBT"
    ),
    TRT_PDFVIII_VWF = list(
      description        = "Plasma-derived factor VIII / von Willebrand factor concentrate treatment-arm indicator (1 = the infusion on this dose record is a pdFVIII:VWF concentrate; 0 = it is not).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = rADAMTS13 (the reference treatment, whose relative ADAMTS13 activity is fixed to 1)",
      notes              = "Required input. Per-dose-record covariate; mutually exclusive with TRT_PBT. Separated from TRT_PBT because pdFVIII:VWF concentrates carry far less ADAMTS13 than plasma itself (93.3% vs 39.0% lower than rADAMTS13), which the paper attributes to the variable ADAMTS13 content across PBT preparations (paper Methods, 'the type of PBT (fresh frozen plasma, solvent/detergent-treated plasma, or pdFVIII:VWF concentrates) as different PBTs can have highly variable rADAMTS13 content').",
      source_name        = "FVIII:VWF concentrates"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 65L,
    n_studies      = 3L,
    age_range      = "0-70 years at enrolment (both phase III trials enrolled ages 0-70). Age-group split of the 65-patient PK analysis set (Table 2): <6 years 4 (6.2%), 6 to <12 years 4 (6.2%), 12 to <18 years 6 (9.2%), >=18 years 51 (78.5%).",
    weight_range   = "18.3-130.0 kg (Table 2, overall)",
    weight_median  = "68.7 kg (Table 2, overall median; also the allometric reference weight)",
    sex_female_pct = 60.0,
    race_ethnicity = c(White = 61.5, Black_African_American = 3.1, Asian = 16.9, Multiple = 1.5, Missing = 16.9),
    disease_state  = "Congenital thrombotic thrombocytopenic purpura (cTTP), an ultra-rare hereditary ADAMTS13 deficiency diagnosed by ADAMTS13 activity <10% of normal in the absence of an acquired inhibitor.",
    dose_range     = "rADAMTS13 40 IU/kg IV once weekly (Q1W) or once every 2 weeks (Q2W) for prophylaxis (the phase I dose-escalation study also studied 5, 20 and 40 IU/kg); PBT approximately 10 IU/kg IV.",
    regions        = "Multinational (NCT02216084 phase I; NCT03393975 phase III crossover; NCT04683003 phase IIIb continuation). Per-region breakdown not reported; country of origin is reported only as Chinese 3 (4.6%) and Japanese 7 (10.8%).",
    notes          = "Baseline characteristics are in Table 2 of the paper. The PopPK analysis set is 65 unique patients contributing 2,462 samples with measurable ADAMTS13 activity across the three trials (paper Results, 'PopPK analysis'). ADAMTS13 activity was measured with the FRETS-VWF73 assay; samples below the limit of quantitation were set to missing rather than imputed, so the model carries no BLQ handling and no endogenous ADAMTS13 baseline term. Number of PK samples was very limited in patients aged <12 years (paper Discussion, limitations)."
  )

  ini({
    # ---- Structural parameters, for the reference 68.7 kg patient ----
    # Table 1 reports clearances and volumes on the linear scale; they are
    # log-transformed here per library convention.
    lcl <- log(0.0398); label("Clearance (L/h) for a 68.7 kg patient")                             # Table 1, "Clearance, L/h" = 0.0398 (RSE 10.7%; bootstrap median 0.0402, 95% CI 0.0288-0.0506)
    lvc <- log(2.69);   label("Central volume of distribution (L) for a 68.7 kg patient")          # Table 1, "Central volume of distribution, L" = 2.69 (RSE 4.99%; bootstrap median 2.71, 95% CI 2.35-3.00)
    lq  <- log(0.0456); label("Intercompartmental clearance (L/h) for a 68.7 kg patient")          # Table 1, "Peripheral clearance, L/h" = 0.0456 (RSE 12.1%; bootstrap median 0.0487)
    lvp <- log(3.71);   label("Peripheral volume of distribution (L) for a 68.7 kg patient")       # Table 1, "Peripheral volume of distribution, L" = 3.71 (RSE 51.4%; bootstrap median 3.43, 95% CI 1.04-19.9)

    # ---- Allometric exponents (held fixed by the authors) ----
    # Paper Methods: "Fixed allometric exponents were used to describe the
    # effect of body weight on clearance and central volume of distribution".
    # Table 1 prints the exponents as literal superscripts (0.75 and 1.0) with
    # no RSE, bootstrap median or CI, which is the table's signal for a fixed
    # parameter.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on clearance and intercompartmental clearance (unitless)")  # Table 1, "x (WT/68.7)^0.75" rows on Clearance and Peripheral clearance
    e_wt_vc  <- fixed(1.0);  label("Allometric exponent on central and peripheral volume (unitless)")               # Table 1, "x (WT/68.7)^1.0" rows on both volumes

    # ---- Treatment-type effect on relative ADAMTS13 activity ----
    # Table 1 prints this block as a multiplier on "Relative ADAMTS13
    # activity", fixed to 1 for rADAMTS13 and multiplied by (1 - theta) for
    # the two PBT categories. The estimated thetas are therefore fractional
    # reductions; the bootstrap medians are reported as the signed quantities
    # (-0.374 and -0.915), confirming the sign convention used here.
    e_trt_pbt  <- -0.390; label("Fractional change in relative ADAMTS13 activity for plasma / S-D-treated plasma vs rADAMTS13")   # Table 1, "x (1-0.390) if PBT" (RSE 6.14%; bootstrap median -0.374, 95% CI -0.427 to -0.270)
    e_trt_pdvw <- -0.933; label("Fractional change in relative ADAMTS13 activity for pdFVIII:VWF concentrates vs rADAMTS13")      # Table 1, "x (1-0.933) if FVIII:VWF concentrates" (RSE 2.06%; bootstrap median -0.915, 95% CI -0.986 to -0.828)

    # ---- Inter-individual variability ----
    # Table 1 reports IIV as a percent coefficient of variation, with the
    # back-transform given in its own footnote a: "IIV is presented as a %
    # coefficient of variation derived as (100 x [exp(omega^2)-1]^0.5)".
    # Inverting that footnote gives the variances used here:
    #   omega^2 = log(1 + CV^2)
    #   CL: log(1 + 0.363^2) = 0.1237819
    #   Vc: log(1 + 0.254^2) = 0.0625202
    # The paper reports no IIV on Q or Vp (Table 1 leaves those cells blank),
    # and reports no eta correlations, so the OMEGA matrix is diagonal.
    etalcl ~ 0.1237819  # Table 1, "IIV, % (RSE%)" on Clearance = 36.3% (RSE 45.8%), shrinkage 14.3%; converted via footnote a
    etalvc ~ 0.0625202  # Table 1, "IIV, % (RSE%)" on Central volume = 25.4% (RSE 20.9%), shrinkage 4.9%; converted via footnote a

    # ---- Residual unexplained variability ----
    # Combined additive + proportional. The additive term is reported in IU/L
    # (Table 1 labels it "Additive (IU/L): 79.9") and is used verbatim here
    # because this model's concentration unit is IU/L; expressed in the
    # paper's reporting unit it is 0.0799 IU/mL, about 7% of the ~1.1 IU/mL
    # rADAMTS13 peak.
    addSd  <- 79.9;  label("Additive residual error (IU/L)")             # Table 1, "Error model, Additive (IU/L): 79.9" (RSE 16.4%; bootstrap median 75.9, 95% CI 41.0-108.0)
    propSd <- 0.204; label("Proportional residual error (fraction)")     # Table 1, "Error model, Proportional (Fraction): 0.204" (RSE 14.7%; bootstrap median 0.203, 95% CI 0.124-0.275)
  })

  model({
    # 1. Relative ADAMTS13 activity of the administered product.
    #    Table 1 parameterises this as a multiplicative chain on a reference
    #    value fixed to 1 for rADAMTS13. TRT_PBT and TRT_PDFVIII_VWF are
    #    mutually exclusive, so the product form below is identical to the
    #    additive-indicator form; the product is written to mirror the
    #    table's printed "x (1 - theta)" structure.
    frel <- (1 + e_trt_pbt * TRT_PBT) * (1 + e_trt_pdvw * TRT_PDFVIII_VWF)

    # 2. Individual parameters, allometrically scaled to the 68.7 kg reference.
    cl <- exp(lcl + etalcl) * (WT / 68.7)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 68.7)^e_wt_vc
    q  <- exp(lq)           * (WT / 68.7)^e_wt_cl
    vp <- exp(lvp)          * (WT / 68.7)^e_wt_vc

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Zero-order input is supplied by the event table
    #    (rate= or dur= on the dose record), not by the model, so the
    #    infusion enters `central` directly.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Relative ADAMTS13 activity acts as a relative bioavailability on the
    #    delivered dose amount.
    f(central) <- frel

    # 6. Observation and error. Cc is plasma ADAMTS13 activity in IU/L;
    #    divide by 1000 to obtain the IU/mL used throughout the paper and by
    #    the companion exposure-response models.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
