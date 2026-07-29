Janssen_2017_cabazitaxel <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous cabazitaxel in",
    "10 men with metastatic castration-resistant prostate cancer, with",
    "individual midazolam clearance as a CYP3A metabolic-phenotype",
    "covariate on cabazitaxel clearance (Janssen 2017 Table 2B 'Metabolic",
    "phenotype model'). The covariate enters via a linear deviation from",
    "the population midazolam clearance reference of 26 L/h. The",
    "individual midazolam clearance is sourced from the companion Janssen",
    "2017 midazolam model (see modellib('Janssen_2017_midazolam'))."
  )
  reference <- paste(
    "Janssen A, Verkleij CPM, van der Vlist A, Mathijssen RHJ,",
    "Bloemendal HJ, ter Heine R (2017).",
    "Towards better dose individualisation: metabolic phenotyping to",
    "predict cabazitaxel pharmacokinetics in men with prostate cancer.",
    "Br J Cancer 116(10):1312-1317.",
    "doi:10.1038/bjc.2017.91.",
    "The CYP3A4 covariate column carries the individual midazolam",
    "clearance produced by the companion model; see",
    "modellib('Janssen_2017_midazolam').",
    sep = " "
  )
  vignette <- "Janssen_2017_cabazitaxel_phenotyping"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CYP3A4 = list(
      description        = "Individual midazolam clearance used as a CYP3A metabolic-phenotype probe",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant per-subject covariate. The column carries each",
        "patient's empirical-Bayes individual midazolam clearance",
        "(L/h) from the companion Janssen 2017 one-compartment IV",
        "midazolam model (modellib('Janssen_2017_midazolam')). The",
        "covariate enters cabazitaxel CL via a linear deviation",
        "(Janssen 2017 Equation (1)):",
        "CL_cabazitaxel = CLbase + y_clearance * (CLmdz_i - 26),",
        "with CLbase = 119 L/h, y_clearance = 1.71, and the population",
        "midazolam clearance reference 26 L/h taken from Table 2A. The",
        "canonical CYP3A4 column accepts probe-substrate-derived CYP3A",
        "activity scores in paper-specific units (the Ter Heine 2014",
        "tamoxifen extraction uses dextromethorphan-probe CYP3A4/5 CL",
        "in ng/L; this extraction uses midazolam CL in L/h). The probe",
        "in Janssen 2017 is midazolam alone (a near-pure CYP3A",
        "substrate), so the column captures CYP3A activity (CYP3A4 +",
        "CYP3A5 are not distinguished). Janssen 2017 Methods 'Modelling'",
        "investigated linear, power and exponential covariate functions",
        "and reported that the linear form best described the",
        "midazolam-CL / cabazitaxel-CL relationship (Results)."
      ),
      source_name        = "CLmdz_i"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "65-77 years (median 67)",
    weight_range   = "Body surface area 1.76-2.34 m^2 (median 1.95-1.96 m^2); body weight not reported in the paper",
    sex_female_pct = 0,
    disease_state  = "Metastatic castration-resistant prostate cancer; all previously treated with docetaxel",
    dose_range     = "Cabazitaxel 25 mg/m^2 by intravenous infusion as part of routine palliative care (median absolute dose 46.25 mg, range 38-50 mg per Table 1)",
    regions        = "Meander Medical Center, Amersfoort, The Netherlands (single-centre)",
    notes          = paste(
      "Eight plasma samples per cabazitaxel infusion at 0, 30, 60, 120,",
      "240, 360, 600 min post-infusion plus one sample at 24 h",
      "(Methods). Bioanalysis by validated LC-MS/MS with LOQ = 1 ng/mL;",
      "inter- and intra-run precision < 9%, accuracy within 12% of",
      "nominal (Bioanalytics). Hepatic function (bilirubin) not",
      "decreased in any patient. Estimation with NONMEM 7.3.0, first-",
      "order conditional estimation with interaction. Body surface area",
      "for dose calculation per the DuBois and DuBois formula (Methods).",
      "Sex split: 100% male (advanced prostate cancer cohort)."
    )
  )

  ini({
    # Final-model parameter values from Janssen 2017 Table 2B "Metabolic
    # phenotype model" (the column that includes midazolam clearance as a
    # linear covariate on cabazitaxel CL). The "Base model" column in the
    # same table is the development-stage model without the covariate and
    # is not extracted here (default policy: extract the final/headline
    # model per references/replicate-author-structure.md).
    lcl <- log(119.0)
    label("Cabazitaxel clearance baseline CLbase at population midazolam CL of 26 L/h (L/h)")  # Table 2B metab-phenotype: CLbase = 119 L/h, RSE 29.6%
    lvc <- log(142.0)
    label("Cabazitaxel central volume of distribution V1 (L)")                                  # Table 2B metab-phenotype: V1 = 142 L, RSE 35.1%
    lvp <- log(2090.0)
    label("Cabazitaxel peripheral volume of distribution V2 (L)")                               # Table 2B metab-phenotype: V2 = 2090 L, RSE 30.2%
    lq  <- log(220.0)
    label("Cabazitaxel inter-compartmental clearance Q (L/h)")                                  # Table 2B metab-phenotype: Q = 220 L/h, RSE 47.9%

    # Linear covariate effect of individual midazolam clearance on
    # cabazitaxel clearance. Janssen 2017 Equation (1):
    #   CL = CLbase + y_clearance * (CLmdz_i - CLmdz_pop)
    # where CLmdz_pop = 26 L/h (the population midazolam clearance from
    # Table 2A). The slope is linear on the L/h scale, not log-scale or
    # power-form; the linear form was selected over power and exponential
    # alternatives (Methods 'Modelling'; Results).
    e_cyp3a4_cl <- 1.71
    label("Linear slope of cabazitaxel CL on (midazolam CL - 26 L/h) (L/h per L/h)")            # Table 2B metab-phenotype: y_clearance = 1.71, RSE 69%

    # IIV is reported as %CV; convert to internal log-scale variance via
    # omega^2 = log(1 + CV^2).
    #   CL : 10.2% -> log(1 + 0.102^2) = 0.01035
    # The Final model retains IIV on CL only. Vc, Vp, and Q do not carry
    # IIV in Table 2B (no entries in the IIV section other than CL_CBZ).
    # Note the high RSE (294%) and shrinkage (42%) reported by the
    # authors -- the residual IIV after the covariate is small and
    # imprecisely estimated; see vignette Errata.
    etalcl ~ 0.01035   # Table 2B metab-phenotype: IIV CL_CBZ 10.2% CV (RSE 294.2%, shrinkage 41.7%)

    # Proportional residual error.
    propSd <- 0.333
    label("Proportional residual error for cabazitaxel plasma concentration (fraction)")        # Table 2B metab-phenotype: residual error 33.3% CV (RSE 19.7%, shrinkage 2.3%)
  })

  model({
    # Janssen 2017 Equation (1): typical cabazitaxel CL deviates linearly
    # from CLbase as the individual midazolam clearance (CYP3A4 column,
    # L/h) deviates from the population reference 26 L/h. The IIV (etalcl)
    # then applies multiplicatively to the covariate-adjusted typical CL.
    cl_typ <- exp(lcl) + e_cyp3a4_cl * (CYP3A4 - 26)
    cl     <- cl_typ * exp(etalcl)
    vc     <- exp(lvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    # Two-compartment IV-infusion disposition. Cabazitaxel was given by
    # intravenous infusion at 25 mg/m^2; supply the infusion rate or
    # duration via the event table (rate / dur columns) and the dose into
    # central.
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
