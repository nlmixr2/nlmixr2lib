Bauer_2023_vonicogAlfa <- function() {
  description <- paste(
    "Two-compartment population PK model for von Willebrand factor:ristocetin",
    "cofactor activity (VWF:RCo) after intravenous recombinant von Willebrand",
    "factor (rVWF, vonicog alfa) coupled to an indirect-response PK/PD model",
    "for endogenous factor VIII activity (FVIII:C), from Bauer 2023.",
    "VWF:RCo disposition is linear two-compartment with first-order",
    "elimination from the central compartment plus an additive endogenous",
    "background activity E_VWF, fixed at half the assay LLOQ (0.5 IU/dL) for",
    "von Willebrand disease (VWD) type 3. CL and Q are allometrically scaled",
    "by body weight with a fixed exponent of 0.75, Vc and Vp with a fixed",
    "exponent of 1, both referenced to 75 kg; Vc additionally decreases with",
    "hematocrit as (HCT/40)^-0.334. FVIII:C is a turnover pool with zero-order",
    "production kin = FVIII0 * kout and first-order removal kout (15.9 1/h)",
    "that VWF:RCo inhibits through an Imax function,",
    "1 - Imax * VWF:RCo / (IC50 + VWF:RCo), so that rising VWF:RCo protects",
    "FVIII from clearance; FVIII0 is the theoretical baseline FVIII:C at",
    "VWF:RCo = 0 IU/dL and scales with hematocrit as (HCT/39)^-0.571.",
    "Fitted to 1664 VWF:RCo samples from 79 patients across four studies",
    "(VWD types 1, 2, 3 and severe hemophilia A) and 686 FVIII:C samples from",
    "the 41 patients of the two phase 3 VWD studies."
  )
  reference <- paste(
    "Bauer A, Friberg-Hietala S, Smania G, Wolfsegger M.",
    "Pharmacokinetic-pharmacodynamic comparison of recombinant and",
    "plasma-derived von Willebrand factor in patients with von Willebrand",
    "disease type 3.",
    "J Blood Med. 2023;14:399-411.",
    "doi:10.2147/JBM.S395845.",
    sep = " "
  )
  vignette <- "Bauer_2023_vonWillebrandFactor"

  # `fviii` is the FVIII:C turnover pool of the indirect-response PD model
  # (Supplementary Figure 1). It is a measured plasma analyte with its own
  # production / removal balance rather than an abstract effect site, so it is
  # declared as a paper-specific compartment rather than reusing `effect`.
  paper_specific_compartments <- c("fviii")

  units <- list(
    time          = "h",
    dosing        = "IU (VWF:RCo activity)",
    concentration = "IU/dL (VWF:RCo activity; FVIII:C activity)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "vonicog alfa", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "vonicog alfa", units = NA_character_, specimen = "plasma", verified = FALSE),
    fviii       = list(analyte = "endogenous factor VIII activity (FVIII:C)", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight. Allometric structural covariate on all four",
        "disposition parameters with fixed exponents and a 75 kg reference:",
        "CL and Q scale as (WT/75)^0.75, Vc and Vp as (WT/75)^1",
        "(Bauer 2023 Supplementary Results, PK Model for pdVWF/FVIII, the",
        "equations that define the shared rVWF parameterisation). The",
        "exponents were estimated at 0.54 (CL) and 0.33 (volumes) but",
        "retained at the canonical fixed values because the estimates were",
        "not in line with previous experience (Holford 1996). Median (range)",
        "in the VWF:RCo PK analysis population was 74 (43.8-145) kg",
        "(Bauer 2023 Table 2)."
      ),
      source_name        = "WT"
    ),
    HCT = list(
      description        = "Hematocrit -- packed red blood cell volume fraction.",
      units              = "% (volume fraction x 100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline hematocrit. Two effects, both power functions. (1) On the",
        "central volume of distribution, (HCT/40)^-0.334: a higher red-cell",
        "volume fraction leaves a smaller plasma volume, and VWF:RCo activity",
        "is confined to plasma, so Vc falls as HCT rises (Bauer 2023 Table 4,",
        "'Hematocrit effect on Vc'; Supplementary Results Vc equation).",
        "(2) On the theoretical baseline FVIII:C, (HCT/39)^-0.571",
        "(Bauer 2023 Table 4, 'Hematocrit effect on baseline FVIII';",
        "Supplementary Results baseline FVIII equation). The paper reports",
        "hematocrit in L/L with reference values 0.4 L/L (Vc) and 0.39 L/L",
        "(baseline FVIII); this model file uses the canonical percent unit,",
        "so the references become 40 % and 39 %. The power ratio is",
        "numerically identical under either unit as long as the supplied HCT",
        "column and the reference share a unit. Median (range) in the",
        "VWF:RCo PK analysis population was 0.417 (0.310-0.480) L/L, i.e.",
        "41.7 (31.0-48.0) % (Bauer 2023 Table 2)."
      ),
      source_name        = "Hematocrit"
    )
  )

  # Covariates screened in the stepwise covariate-modelling procedure
  # (Bauer 2023 Supplementary Table 1) but not retained in the final model.
  # Documented here for provenance; not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL, Vc, Vp and on FVIII0, IC50 (Supplementary Table 1);",
        "not retained. Median (range) 35 (18-70) years in the VWF:RCo PK",
        "population and 36 (18-70) years in the FVIII PK/PD population",
        "(Table 2)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL, Vc, Vp and on FVIII0, IC50 (Supplementary Table 1);",
        "not retained. 33/79 (42 %) female in the VWF:RCo PK population and",
        "21/41 (51 %) female in the FVIII PK/PD population (Table 2)."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was screened on CL only (Supplementary Table 1) and not",
        "retained. 73/79 (92 %) White in the VWF:RCo PK population (Table 2)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was screened on CL only (Supplementary Table 1) and not",
        "retained. 6/79 (8 %) Asian in the VWF:RCo PK population (Table 2)."
      )
    ),
    BLOOD_GROUP_O = list(
      description = "ABO blood group O indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL, Vc, Vp and on FVIII0, IC50 (Supplementary Table 1);",
        "not retained. Blood-group composition of the analysis population is",
        "not reported in Bauer 2023."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 79L,
    n_studies      = 4L,
    age_range      = "18-70 years",
    age_median     = "35 years",
    weight_range   = "43.8-145 kg",
    weight_median  = "74 kg",
    sex_female_pct = 42,
    race_ethnicity = c(White = 92, Asian = 8),
    disease_state  = paste(
      "Adults with severe von Willebrand disease (type 1 n=5, type 2 n=7,",
      "type 3 n=57) or severe hemophilia A with FVIII activity < 1 % (n=10).",
      "The hemophilia A cohort was retained in the VWF:RCo PK model because,",
      "once individual endogenous VWF:RCo levels were accounted for, no PK",
      "differences between VWD types and hemophilia A remained."
    ),
    dose_range     = "2-80 IU/kg VWF:RCo intravenously (rVWF, vonicog alfa)",
    regions        = "North America, Europe, Australia, Japan, India, Taiwan, Turkey",
    n_observations = 1664L,
    n_observations_pd = 686L,
    biomarkers     = c(
      "VWF:RCo -- von Willebrand factor:ristocetin cofactor activity (IU/dL)",
      "FVIII:C -- factor VIII activity by one-stage clotting assay (IU/dL)"
    ),
    notes          = paste(
      "VWF:RCo PK model: 1664 samples from 79 patients pooled across",
      "NCT00816660 (phase 1, dose escalation, 29 patients / 479 samples),",
      "NCT01410227 (phase 3 on-demand, 31 / 694), NCT02283268 (phase 3",
      "surgery, 11 / 241) and EudraCT 2011-004314-42 (phase 1 hemophilia A",
      "proof of concept, 10 / 250); Bauer 2023 Table 1. FVIII PK/PD model:",
      "686 FVIII:C samples from the 41 patients of the two phase 3 studies",
      "only, because patients in the phase 1 studies received concomitant",
      "exogenous FVIII concentrates. Patients with inhibitory or binding",
      "antibodies to VWF or FVIII (n=7), observations associated with a",
      "bleeding event (n=4) and observations after coadministration of",
      "plasma-derived replacement therapy (n=308) were excluded. VWF:RCo LLOQ",
      "was 1 IU/dL in NCT00816660 and 8 IU/dL in the other three studies;",
      "FVIII:C LLOQ was 1 IU/dL. Samples below LLOQ were retained and modelled",
      "with the likelihood-based M3 method. Estimation used FOCEI in NONMEM",
      "7.3.0 / 7.4.0; the PK/PD model was fitted sequentially, conditioned on",
      "the individual VWF:RCo PK parameters. The parameter values encoded here",
      "are the VWD type 3 parameterisation: E_VWF is the type 3 endogenous",
      "level and its interindividual variability was set to 0 for that group.",
      "E_VWF values for VWD types 1 and 2 and for severe hemophilia A were",
      "estimated as structural covariates but are not reported in the paper",
      "or its supplement."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # VWF:RCo disposition -- typical values for a 75 kg patient with a
    # hematocrit of 40 % (Bauer 2023 Table 4, rVWF Model column).
    # -----------------------------------------------------------------
    lcl <- log(2.10)  ; label("Clearance of VWF:RCo activity (dL/h)")                     # Bauer 2023 Table 4: CL = 2.10 dL/h, RSE 8.10 %
    lvc <- log(43.5)  ; label("Central volume of distribution for VWF:RCo (dL)")           # Bauer 2023 Table 4: Vc = 43.5 dL, RSE 4.98 %
    lq  <- log(2.29)  ; label("Intercompartmental clearance of VWF:RCo activity (dL/h)")   # Bauer 2023 Table 4: Q = 2.29 dL/h, RSE 12.2 %
    lvp <- log(15.8)  ; label("Peripheral volume of distribution for VWF:RCo (dL)")        # Bauer 2023 Table 4: Vp = 15.8 dL, RSE 4.47 %

    # Endogenous VWF:RCo background activity, fixed at half the 1 IU/dL
    # assay LLOQ for VWD type 3. Kept on the linear scale so that a user
    # can set bl_vwf = 0 to model a complete absence of endogenous VWF.
    bl_vwf <- fixed(0.500); label("Endogenous VWF:RCo activity E_VWF, VWD type 3 (IU/dL)") # Bauer 2023 Table 4: E_VWF, VWD type 3 = 0.500 IU/dL (FIX); Supplementary Results: "the typical endogenous level for this population was set to half the LLOQ (0.5 IU/dL)"

    # -----------------------------------------------------------------
    # Allometric body-weight exponents, fixed at the canonical values
    # with a 75 kg reference. Shared exponents: one value on CL and Q,
    # one value on Vc and Vp.
    # -----------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.750); label("Allometric exponent of (WT/75) on CL and Q (unitless)")  # Bauer 2023 Table 4: WT effect on CL and Q = 0.750 (FIX)
    e_wt_vc_vp <- fixed(1.00) ; label("Allometric exponent of (WT/75) on Vc and Vp (unitless)") # Bauer 2023 Table 4: WT effect on Vc and Vp = 1.00 (FIX)

    # Hematocrit on Vc, estimated in the rVWF model.
    e_hct_vc <- -0.334; label("Power exponent of (HCT/40) on Vc (unitless)")                    # Bauer 2023 Table 4: Hematocrit effect on Vc = -0.334, RSE 37.5 %

    # -----------------------------------------------------------------
    # FVIII:C indirect-response PD model (Bauer 2023 Table 4, rVWF
    # Model column, FVIII PK/PD model block). Precision for this block
    # was assessed by bootstrap and is not reported per parameter
    # (Table 4 footnote a).
    # -----------------------------------------------------------------
    lrbase <- log(0.500); label("Theoretical baseline FVIII:C at VWF:RCo = 0 IU/dL, FVIII0 (IU/dL)") # Bauer 2023 Table 4: Baseline FVIII = 0.500 IU/dL; footnote d: "the theoretical baseline FVIII:C assuming that no VWF is affecting kout"
    lkout  <- log(15.9) ; label("First-order removal rate constant of FVIII:C, kout (1/h)")          # Bauer 2023 Table 4: kout = 15.9 1/h
    limax  <- log(0.998); label("Maximum fractional inhibition of kout by VWF:RCo, Imax (unitless)") # Bauer 2023 Table 4: Imax = 0.998
    lec50  <- log(0.0658); label("VWF:RCo activity giving 50 % inhibition of kout, IC50 (IU/dL)")    # Bauer 2023 Table 4: IC50 = 0.0658 IU/dL

    e_hct_rbase <- -0.571; label("Power exponent of (HCT/39) on baseline FVIII:C (unitless)")        # Bauer 2023 Table 4: Hematocrit effect on baseline FVIII = -0.571

    # -----------------------------------------------------------------
    # IIV -- exponential (log-normal multiplicative), P_i = TVP * exp(eta)
    # per Supplementary Methods, Stochastic Model. Table 4 reports the
    # magnitudes as CV; convert to the log-scale variance with the
    # log-normal identity omega^2 = log(1 + CV^2).
    # -----------------------------------------------------------------
    etalcl   ~ 0.149110  # log(1 + 0.401^2); Bauer 2023 Table 4: IIV CL, CV = 0.401, RSE 9.88 %
    etalvc   ~ 0.081823  # log(1 + 0.292^2); Bauer 2023 Table 4: IIV Vc, CV = 0.292, RSE 10.7 %
    etalrbase ~ 0.083988 # log(1 + 0.296^2); Bauer 2023 Table 4: IIV baseline FVIII, CV = 0.296
    etalec50 ~ 0.296947  # log(1 + 0.588^2); Bauer 2023 Table 4: IIV IC50, CV = 0.588

    # -----------------------------------------------------------------
    # Residual unexplained variability -- combined additive and
    # proportional on each output (Supplementary Methods, Stochastic
    # Model: "the residual unexplained variability was described by a
    # combined additive and proportional error model").
    # -----------------------------------------------------------------
    propSd <- 0.138; label("Proportional residual error on VWF:RCo (fraction)")  # Bauer 2023 Table 4: Proportional RUV, CV = 0.138, RSE 1.26 %
    addSd  <- 2.80 ; label("Additive residual error on VWF:RCo (IU/dL)")         # Bauer 2023 Table 4: Additive RUV = 2.80 IU/dL, RSE 2.08 %

    propSd_Cfviii <- 0.190; label("Proportional residual error on FVIII:C (fraction)") # Bauer 2023 Table 4: FVIII PK/PD model, Proportional RUV, CV = 0.190
    addSd_Cfviii  <- 1.55 ; label("Additive residual error on FVIII:C (IU/dL)")        # Bauer 2023 Table 4: FVIII PK/PD model, Additive RUV = 1.55 IU/dL
  })

  model({
    # -----------------------------------------------------------------
    # 1. Reference values for the covariate models.
    #    Body weight 75 kg and hematocrit 40 % (0.4 L/L) for the
    #    disposition parameters; hematocrit 39 % (0.39 L/L) for the
    #    baseline FVIII:C (Bauer 2023 Supplementary Results).
    # -----------------------------------------------------------------
    wtRef     <- 75  # kg
    hctRef    <- 40  # %, reference for the Vc covariate model
    hctRefBl  <- 39  # %, reference for the baseline FVIII:C covariate model

    wtCl <- (WT / wtRef)^e_wt_cl_q
    wtV  <- (WT / wtRef)^e_wt_vc_vp

    # -----------------------------------------------------------------
    # 2. Individual VWF:RCo disposition parameters. IIV is supported on
    #    CL and Vc only (Supplementary Results, PK Model for
    #    pdVWF/FVIII: "Interindividual variability terms were supported
    #    on CL and Vc").
    # -----------------------------------------------------------------
    cl <- exp(lcl + etalcl) * wtCl
    q  <- exp(lq)           * wtCl
    vc <- exp(lvc + etalvc) * wtV * (HCT / hctRef)^e_hct_vc
    vp <- exp(lvp)          * wtV

    # -----------------------------------------------------------------
    # 3. Micro-constants.
    # -----------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # -----------------------------------------------------------------
    # 4. VWF:RCo disposition -- linear two-compartment with first-order
    #    elimination from the central compartment. rVWF is given
    #    intravenously, so dose enters `central` directly.
    # -----------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed VWF:RCo activity is the exogenous contribution plus the
    # endogenous background E_VWF (Supplementary Methods: "a component
    # describing endogenous VWF:RCo activity levels").
    Cc <- central / vc + bl_vwf

    # -----------------------------------------------------------------
    # 5. FVIII:C indirect-response model. VWF:RCo inhibits the
    #    first-order removal of FVIII:C through an Imax function
    #    (Supplementary Methods:
    #    VWF effect = 1 - Imax * VWF:RCo / (IC50 + VWF:RCo)),
    #    and kin follows from the steady-state condition
    #    kin = FVIII0 * kout.
    # -----------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase) * (HCT / hctRefBl)^e_hct_rbase
    kout  <- exp(lkout)
    imax  <- exp(limax)
    ec50  <- exp(lec50 + etalec50)

    vwfEffect <- 1 - imax * Cc / (ec50 + Cc)
    kin       <- rbase * kout

    d/dt(fviii) <- kin - kout * vwfEffect * fviii

    # Pre-dose initial condition. FVIII0 is defined at VWF:RCo = 0 IU/dL,
    # but a VWD type 3 patient carries the endogenous E_VWF activity, so
    # the drug-free steady state is kin / (kout * VWF effect at E_VWF)
    # rather than FVIII0 itself. Solving d/dt(fviii) = 0 at Cc = bl_vwf
    # gives the expression below; for the typical patient this is
    # 4.17 IU/dL, consistent with the low but non-zero FVIII:C of type 3
    # VWD. The paper does not state the initial condition explicitly
    # (see the vignette's Assumptions and deviations section).
    fviii(0) <- rbase / (1 - imax * bl_vwf / (ec50 + bl_vwf))

    Cfviii <- fviii

    # -----------------------------------------------------------------
    # 6. Observations and residual error.
    # -----------------------------------------------------------------
    Cc     ~ add(addSd) + prop(propSd)
    Cfviii ~ add(addSd_Cfviii) + prop(propSd_Cfviii)
  })
}
