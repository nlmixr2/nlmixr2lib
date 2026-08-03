Bauer_2023_vwfFviii_humateP <- function() {
  description <- paste(
    "Two-compartment population PK model for von Willebrand factor:ristocetin",
    "cofactor activity (VWF:RCo) after intravenous plasma-derived von",
    "Willebrand factor / factor VIII concentrate (pdVWF/FVIII; Humate-P,",
    "VWF:RCo/FVIII:C 2.4:1) coupled to an indirect-response PK/PD model for",
    "factor VIII activity (FVIII:C), from Bauer 2023. The structure is the one",
    "developed for recombinant VWF (see modellib('Bauer_2023_vonicogAlfa'))",
    "refitted to the plasma-derived product: linear two-compartment VWF:RCo",
    "disposition with first-order elimination from the central compartment",
    "plus an additive endogenous background E_VWF fixed at half the assay LLOQ",
    "(0.5 IU/dL) for von Willebrand disease (VWD) type 3, allometric body",
    "weight on CL and Q (exponent 0.75) and on Vc and Vp (exponent 1) with a",
    "75 kg reference, and (HCT/40)^-0.334 on Vc. FVIII:C is a turnover pool",
    "whose first-order removal kout is inhibited by VWF:RCo through",
    "1 - Imax * VWF:RCo / (IC50 + VWF:RCo). Because the product delivers FVIII",
    "as well as VWF, a volume of distribution for FVIII (V FVIII = 32.9 dL) is",
    "added so that the administered FVIII:C dose enters the FVIII pool; the",
    "elimination of plasma-derived FVIII is assumed identical to that of",
    "endogenous FVIII and is therefore carried by kout. The system-specific",
    "parameters FVIII0, kout and the hematocrit effect on FVIII0 are fixed to",
    "the recombinant-VWF model estimates. VWF:RCo clearance is roughly twice",
    "that of recombinant VWF (4.14 vs 2.10 dL/h), giving a 1.76-fold shorter",
    "mean residence time. Fitted to 281 VWF:RCo and FVIII:C samples from 20",
    "patients with VWD type 3."
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
  # (Supplementary Figure 1). It is a measured plasma analyte that also
  # receives the plasma-derived FVIII dose, rather than an abstract effect
  # site, so it is declared as a paper-specific compartment rather than
  # reusing `effect`.
  paper_specific_compartments <- c("fviii")

  units <- list(
    time          = "h",
    dosing        = "IU (VWF:RCo activity into central; FVIII:C activity into fviii)",
    concentration = "IU/dL (VWF:RCo activity; FVIII:C activity)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight. Allometric structural covariate with fixed",
        "exponents and a 75 kg reference: CL and Q scale as (WT/75)^0.75,",
        "Vc and Vp as (WT/75)^1 (Bauer 2023 Supplementary Results, PK Model",
        "for pdVWF/FVIII equations). A formal covariate analysis was not",
        "possible for the pdVWF/FVIII model given the small number of",
        "patients, so the mechanistically plausible covariates identified in",
        "the rVWF model (body weight and hematocrit) were applied unchanged",
        "(Bauer 2023 Discussion). Median (range) 66.3 (43.8-132) kg",
        "(Bauer 2023 Table 3)."
      ),
      source_name        = "WT"
    ),
    HCT = list(
      description        = "Hematocrit -- packed red blood cell volume fraction.",
      units              = "% (volume fraction x 100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline hematocrit. Two effects, both power functions and both",
        "fixed to the rVWF model estimates. (1) On the central volume of",
        "distribution, (HCT/40)^-0.334: a higher red-cell volume fraction",
        "leaves a smaller plasma volume, and VWF:RCo activity is confined to",
        "plasma, so Vc falls as HCT rises (Bauer 2023 Table 4, 'Hematocrit",
        "effect on Vc', FIX in the pdVWF/FVIII column). (2) On the",
        "theoretical baseline FVIII:C, (HCT/39)^-0.571 (Bauer 2023 Table 4,",
        "'Hematocrit effect on baseline FVIII', FIX). The paper reports",
        "hematocrit in L/L with reference values 0.4 L/L (Vc) and 0.39 L/L",
        "(baseline FVIII); this model file uses the canonical percent unit,",
        "so the references become 40 % and 39 %. The power ratio is",
        "numerically identical under either unit as long as the supplied HCT",
        "column and the reference share a unit. Median (range) 0.406",
        "(0.334-0.483) L/L, i.e. 40.6 (33.4-48.3) % (Bauer 2023 Table 3)."
      ),
      source_name        = "Hematocrit"
    )
  )

  # Covariates screened in the rVWF stepwise covariate-modelling procedure
  # (Bauer 2023 Supplementary Table 1) but not retained. A formal covariate
  # analysis was not feasible for the pdVWF/FVIII model at all (n=20;
  # Bauer 2023 Discussion). Documented here for provenance; not referenced
  # in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the rVWF model on CL, Vc, Vp and on FVIII0, IC50",
        "(Supplementary Table 1); not retained, and not testable in the",
        "pdVWF/FVIII model. Median (range) 30 (18-60) years (Table 3)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the rVWF model on CL, Vc, Vp and on FVIII0, IC50",
        "(Supplementary Table 1); not retained. 11/20 (55 %) female",
        "(Table 3)."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was screened on CL in the rVWF model (Supplementary Table 1)",
        "and not retained. The pdVWF/FVIII cohort was 20/20 (100 %) White",
        "(Table 3), so the covariate is not estimable here."
      )
    ),
    BLOOD_GROUP_O = list(
      description = "ABO blood group O indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the rVWF model on CL, Vc, Vp and on FVIII0, IC50",
        "(Supplementary Table 1); not retained. Blood-group composition of",
        "the analysis population is not reported in Bauer 2023."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "18-60 years",
    age_median     = "30 years",
    weight_range   = "43.8-132 kg",
    weight_median  = "66.3 kg",
    sex_female_pct = 55,
    race_ethnicity = c(White = 100),
    disease_state  = "Adults with von Willebrand disease type 3 (20/20, 100 %)",
    dose_range     = paste(
      "Single dose of 50 IU/kg VWF:RCo pdVWF/FVIII (Humate-P), which at a",
      "VWF:RCo/FVIII:C ratio of 2.4:1 delivers about 20.8 IU/kg FVIII:C"
    ),
    regions        = "Europe, North America",
    n_observations = 281L,
    biomarkers     = c(
      "VWF:RCo -- von Willebrand factor:ristocetin cofactor activity (IU/dL)",
      "FVIII:C -- factor VIII activity by one-stage clotting assay (IU/dL)"
    ),
    notes          = paste(
      "Single crossover cohort (n=22 randomised) of the phase 1 dose-escalation",
      "study NCT00816660, in which patients with VWD type 3 received either",
      "rVWF plus recombinant FVIII (octocog alfa, 50 IU/kg VWF:RCo /",
      "38.5 IU/kg FVIII:C) or pdVWF/FVIII (Humate-P, 50 IU/kg VWF:RCo,",
      "VWF:RCo/FVIII:C 2.4:1). Only the pdVWF/FVIII arm contributes to this",
      "model; measurements after rVWF plus rFVIII were excluded, as were",
      "screening and end-of-study (30 days) samples because of uncertain",
      "dosing history. Two of the 22 patients were excluded (one with no",
      "PK/PD measurements after pdVWF/FVIII, one who dropped out before",
      "dosing), leaving 20 patients and 281 samples collected up to 96 h",
      "post-infusion (Bauer 2023 Results and Table 3). Of the 281 VWF:RCo",
      "samples, 232 were above and 49 below the 1 IU/dL LLOQ; all 281 FVIII:C",
      "samples were above LLOQ. Samples below LLOQ were retained and modelled",
      "with the likelihood-based M3 method. Estimation used FOCEI in NONMEM",
      "7.3.0 / 7.4.0, starting from the rVWF model; the PK/PD model was",
      "fitted sequentially, conditioned on the individual VWF:RCo PK",
      "parameters. Interindividual variability was estimated with",
      "considerable uncertainty in this small cohort (RSE 38.3 % on IIV CL,",
      "versus 9.9 % in the rVWF model). The paper cautions that the findings",
      "are specific to this one concentrate and may not extend to other",
      "pdVWF products, in particular to those containing very little FVIII."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # VWF:RCo disposition -- typical values for a 75 kg patient with a
    # hematocrit of 40 % (Bauer 2023 Table 4, pdVWF/FVIII Model column).
    # -----------------------------------------------------------------
    lcl <- log(4.14)  ; label("Clearance of VWF:RCo activity (dL/h)")                     # Bauer 2023 Table 4: CL = 4.14 dL/h, RSE 12.7 %
    lvc <- log(47.0)  ; label("Central volume of distribution for VWF:RCo (dL)")           # Bauer 2023 Table 4: Vc = 47.0 dL, RSE 16.7 %
    lq  <- log(4.47)  ; label("Intercompartmental clearance of VWF:RCo activity (dL/h)")   # Bauer 2023 Table 4: Q = 4.47 dL/h, RSE 12.3 %
    lvp <- log(19.3)  ; label("Peripheral volume of distribution for VWF:RCo (dL)")        # Bauer 2023 Table 4: Vp = 19.3 dL, RSE 3.12 %

    # Endogenous VWF:RCo background activity, fixed at half the 1 IU/dL
    # assay LLOQ for VWD type 3. Kept on the linear scale so that a user
    # can set bl_vwf = 0 to model a complete absence of endogenous VWF.
    bl_vwf <- fixed(0.500); label("Endogenous VWF:RCo activity E_VWF, VWD type 3 (IU/dL)") # Bauer 2023 Table 4: E_VWF, VWD type 3 = 0.500 IU/dL (FIX); Supplementary Results: "the pre-dose VWF:RCo was below the LLOQ in 18 patients and 1 IU/dL in two patients, and so fixing E_VWF to half the LLOQ at 0.5 IU/dL was considered appropriate"

    # -----------------------------------------------------------------
    # Allometric body-weight exponents, fixed at the canonical values
    # with a 75 kg reference. Shared exponents: one value on CL and Q,
    # one value on Vc and Vp.
    # -----------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.750); label("Allometric exponent of (WT/75) on CL and Q (unitless)")  # Bauer 2023 Table 4: WT effect on CL and Q = 0.750 (FIX)
    e_wt_vc_vp <- fixed(1.00) ; label("Allometric exponent of (WT/75) on Vc and Vp (unitless)") # Bauer 2023 Table 4: WT effect on Vc and Vp = 1.00 (FIX)

    # Hematocrit on Vc, fixed to the rVWF model estimate because the
    # pdVWF/FVIII cohort was too small for a covariate analysis.
    e_hct_vc <- fixed(-0.334); label("Power exponent of (HCT/40) on Vc (unitless)")              # Bauer 2023 Table 4: Hematocrit effect on Vc = -0.334 (FIX in the pdVWF/FVIII model)

    # -----------------------------------------------------------------
    # FVIII:C indirect-response PD model (Bauer 2023 Table 4, pdVWF/FVIII
    # Model column, FVIII PK/PD model block). The system-specific
    # parameters FVIII0, kout and the hematocrit effect on FVIII0 "can be
    # expected to be similar across all VWD type 3 patients and
    # independent of drug product" and were fixed to the rVWF estimates
    # (Supplementary Results, PK/PD Model for pdVWF/FVIII); Table 4
    # footnote c marks them as fixed in this model.
    # -----------------------------------------------------------------
    lrbase <- fixed(log(0.500)); label("Theoretical baseline FVIII:C at VWF:RCo = 0 IU/dL, FVIII0 (IU/dL)") # Bauer 2023 Table 4: Baseline FVIII = 0.500 IU/dL (FIX)
    lkout  <- fixed(log(15.9)) ; label("First-order removal rate constant of FVIII:C, kout (1/h)")          # Bauer 2023 Table 4: kout = 15.9 1/h (FIX)
    limax  <- log(0.994)       ; label("Maximum fractional inhibition of kout by VWF:RCo, Imax (unitless)") # Bauer 2023 Table 4: Imax = 0.994, RSE 0.0378 %
    lec50  <- log(0.0577)      ; label("VWF:RCo activity giving 50 % inhibition of kout, IC50 (IU/dL)")     # Bauer 2023 Table 4: IC50 = 0.0577 IU/dL, RSE 20.4 %

    e_hct_rbase <- fixed(-0.571); label("Power exponent of (HCT/39) on baseline FVIII:C (unitless)")        # Bauer 2023 Table 4: Hematocrit effect on baseline FVIII = -0.571 (FIX)

    # Volume of distribution for FVIII, added so that the plasma-derived
    # FVIII:C dose delivered by the concentrate translates into an FVIII
    # activity increment. Paper-mechanistic volume term for which the
    # standard `vc` shape does not apply, hence the canonical `vd` name.
    lvd <- log(32.9); label("Volume of distribution for FVIII, V FVIII (dL)")                              # Bauer 2023 Table 4: FVIII V = 32.9 dL, RSE 1.65 %; Results: "in line with a previous report (32.8 dL)"

    # -----------------------------------------------------------------
    # IIV -- exponential (log-normal multiplicative), P_i = TVP * exp(eta)
    # per Supplementary Methods, Stochastic Model. Table 4 reports the
    # magnitudes as CV; convert to the log-scale variance with the
    # log-normal identity omega^2 = log(1 + CV^2).
    # -----------------------------------------------------------------
    etalcl    ~ 0.053754 # log(1 + 0.235^2); Bauer 2023 Table 4: IIV CL, CV = 0.235, RSE 38.3 %
    etalvc    ~ 0.130264 # log(1 + 0.373^2); Bauer 2023 Table 4: IIV Vc, CV = 0.373, RSE 51.7 %
    etalrbase ~ 0.045188 # log(1 + 0.215^2); Bauer 2023 Table 4: IIV baseline FVIII, CV = 0.215, RSE 44.4 %
    etalec50  ~ 0.257556 # log(1 + 0.542^2); Bauer 2023 Table 4: IIV IC50, CV = 0.542, RSE 29.4 %

    # -----------------------------------------------------------------
    # Residual unexplained variability -- combined additive and
    # proportional on each output (Supplementary Methods, Stochastic
    # Model: "the residual unexplained variability was described by a
    # combined additive and proportional error model").
    # -----------------------------------------------------------------
    propSd <- 0.0479; label("Proportional residual error on VWF:RCo (fraction)") # Bauer 2023 Table 4: Proportional RUV, CV = 0.0479, RSE 14.1 %
    addSd  <- 1.48  ; label("Additive residual error on VWF:RCo (IU/dL)")        # Bauer 2023 Table 4: Additive RUV = 1.48 IU/dL, RSE 12.4 %

    propSd_Cfviii <- 0.172; label("Proportional residual error on FVIII:C (fraction)") # Bauer 2023 Table 4: FVIII PK/PD model, Proportional RUV, CV = 0.172, RSE 2.50 %
    addSd_Cfviii  <- 2.91 ; label("Additive residual error on FVIII:C (IU/dL)")        # Bauer 2023 Table 4: FVIII PK/PD model, Additive RUV = 2.91 IU/dL, RSE 15.5 %
  })

  model({
    # -----------------------------------------------------------------
    # 1. Reference values for the covariate models.
    #    Body weight 75 kg and hematocrit 40 % (0.4 L/L) for the
    #    disposition parameters; hematocrit 39 % (0.39 L/L) for the
    #    baseline FVIII:C (Bauer 2023 Supplementary Results).
    # -----------------------------------------------------------------
    wtRef    <- 75  # kg
    hctRef   <- 40  # %, reference for the Vc covariate model
    hctRefBl <- 39  # %, reference for the baseline FVIII:C covariate model

    wtCl <- (WT / wtRef)^e_wt_cl_q
    wtV  <- (WT / wtRef)^e_wt_vc_vp

    # -----------------------------------------------------------------
    # 2. Individual VWF:RCo disposition parameters. IIV supported on CL
    #    and Vc only (Supplementary Results).
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
    #    elimination from the central compartment. pdVWF/FVIII is given
    #    intravenously, so the VWF:RCo dose enters `central` directly.
    # -----------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed VWF:RCo activity is the exogenous contribution plus the
    # endogenous background E_VWF.
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
    vd    <- exp(lvd)

    vwfEffect <- 1 - imax * Cc / (ec50 + Cc)
    kin       <- rbase * kout

    d/dt(fviii) <- kin - kout * vwfEffect * fviii

    # The FVIII pool is carried in activity units (IU/dL), so the
    # plasma-derived FVIII:C dose (IU) delivered by the concentrate is
    # divided by V FVIII on the way in. This is the role the paper
    # assigns to V FVIII: "To account for the distribution of the pdFVIII
    # dose, a volume of distribution of FVIII was introduced into the
    # model, V FVIII" (Supplementary Results). Elimination of
    # plasma-derived FVIII is assumed identical to that of endogenous
    # FVIII and is therefore carried by kout, with no separate term.
    f(fviii) <- 1 / vd

    # Pre-dose initial condition. FVIII0 is defined at VWF:RCo = 0 IU/dL,
    # but a VWD type 3 patient carries the endogenous E_VWF activity, so
    # the drug-free steady state is kin / (kout * VWF effect at E_VWF)
    # rather than FVIII0 itself. Solving d/dt(fviii) = 0 at Cc = bl_vwf
    # gives the expression below. The paper does not state the initial
    # condition explicitly (see the vignette's Assumptions and deviations
    # section).
    fviii(0) <- rbase / (1 - imax * bl_vwf / (ec50 + bl_vwf))

    Cfviii <- fviii

    # -----------------------------------------------------------------
    # 6. Observations and residual error.
    # -----------------------------------------------------------------
    Cc     ~ add(addSd) + prop(propSd)
    Cfviii ~ add(addSd_Cfviii) + prop(propSd_Cfviii)
  })
}
