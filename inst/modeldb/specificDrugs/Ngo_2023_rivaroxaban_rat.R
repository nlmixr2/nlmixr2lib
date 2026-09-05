Ngo_2023_rivaroxaban_rat <- function() {
  description <- "Preclinical (rat). Two-compartment population PK model for rivaroxaban (RIV) in rats with a carbamazepine (CBZ) drug-drug-interaction covariate, as reproduced in Ngo 2023 Table 1 (column 'Estimated in Rats') from the original fit in Ngo 2020. Absorption is a parallel mixed zero- and first-order process: a fraction F1 of the dose enters a depot and is absorbed first-order with rate Ka (lymphatic uptake of the lipophilic drug), while the complementary fraction 1 - F1 enters the central compartment directly as a zero-order input of duration D2 after a lag Alag2 (solubility-limited transfer across the enterocyte capillary network and hepatic portal vein). Six days of carbamazepine pretreatment, a strong CYP and P-gp / BCRP inducer, acts through two fractional-change covariate terms: it increases apparent clearance CL/F by 211 % and prolongs the zero-order absorption time D2 by 33.9 %, together reducing the observed rivaroxaban AUC by 57.9 % (3088.5 to 1299.3 ng*h/mL) and Cmax by 43.3 % (540.0 to 306.0 ng/mL). Disposition parameters are reported per kilogram of body weight, so the model multiplies them by the WT covariate. The sibling model Ngo_2023_rivaroxaban_human carries the allometric extrapolation of this fit to a 60 kg human."
  reference   <- "Ngo LT, Yun H-y, Chae J-w. Application of the Population Pharmacokinetics Model-Based Approach to the Prediction of Drug-Drug Interaction between Rivaroxaban and Carbamazepine in Humans. Pharmaceuticals. 2023;16(5):684. doi:10.3390/ph16050684 (Table 1, Figure 2, Equations 1-2 and Sections 2.1, 4.1-4.2). The estimates were originally published in Ngo LT, Yang S, Tran QT, Kim SK, Yun H, Chae J. Effects of Carbamazepine and Phenytoin on Pharmacokinetics and Pharmacodynamics of Rivaroxaban. Pharmaceutics. 2020;12(11):1040. doi:10.3390/pharmaceutics12111040, which Ngo 2023 reproduces in full; this extraction is transcribed from Ngo 2023."
  vignette    <- "Ngo_2023_rivaroxaban"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Ngo 2023 Table 1 reports every rat disposition parameter in per-kilogram units (CL/F and Q/F in L/h/kg, Vc/F and Vp/F in L/kg), so multiplying by WT is the unit conversion those units denote rather than an estimated allometric relationship -- the exponent is exactly 1 for both clearances and volumes, and no reference weight or centring is involved. Ngo 2023 Section 4.3 assumes an average rat body weight of 0.25 kg when converting these per-kilogram values for the interspecies extrapolation, which is the natural default for simulating this model. Rats were individually weighed and dosed on a mg/kg basis (Ngo 2023 Section 4.1).",
      source_name        = "Body weight (implicit in the L/h/kg and L/kg parameter units)"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator: 1 = rat pretreated with carbamazepine 45 mg/kg twice daily for 6 consecutive days before the rivaroxaban dose (the 'test group'), 0 = rivaroxaban administered alone (the 'control group')",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (rivaroxaban alone; the 'control group' of Ngo 2023)",
      notes              = "Time-fixed: pretreatment ran for 6 days before the single rivaroxaban dose on Day 7, so induction is already established and stable across the 24 h PK sampling window. Enters the model twice, both as fractional-change (1 + coefficient * CONMED_CBZ) terms rather than log-multiplicative terms: on apparent clearance CL/F (Ngo 2023 Equation 1, coefficient CL_CBZ = 2.11) and on the zero-order absorption duration D2 (Ngo 2023 Equation 2, coefficient D2_CBZ = 0.339). Ngo 2023 Section 4.2 states that carbamazepine was also tested on the first-order absorption rate constant Ka but was not retained there. On Day 7 rivaroxaban was given 30 minutes after the morning carbamazepine dose.",
      source_name        = "Pretreatment with CBZ (index variable; 'control group' vs 'test group')"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "rat",
    n_subjects     = 12L,
    n_studies      = 1L,
    weight_median  = "0.25 kg assumed as the average rat body weight (Ngo 2023 Section 4.3); individual weights were measured but are not reported in Ngo 2023",
    disease_state  = "Healthy rats, randomised into two groups of six: a control group receiving rivaroxaban alone and a test group pretreated with carbamazepine",
    dose_range     = "Single oral rivaroxaban 3 mg/kg on Day 7; carbamazepine 45 mg/kg orally twice daily on Days 1-7 in the test group (placebo in the control group). The 3 mg/kg rivaroxaban dose was raised above the 2.0 mg/kg body-surface-area equivalent of the 20 mg human dose to keep concentrations within assay sensitivity once carbamazepine reduced exposure (Ngo 2023 Discussion).",
    regions        = "Republic of Korea (Chungnam National University, Daejeon)",
    notes          = "Plasma sampled before the rivaroxaban dose and at 0.25, 0.5, 1, 2, 4, 8, 10 and 24 h after it; water and food available ad libitum in both groups (Ngo 2023 Section 4.1). Animal study approved by the Animal Ethics Committee of Chungnam National University (No. 2019012A-CNU-193, 27 December 2019). Fitted in NONMEM 7.3.0 driven by PsN 4.4.0 through Pirana (Ngo 2023 Section 4.2). The rat plasma unbound fraction of rivaroxaban is 1.3 % versus about 6.5 % in humans, a difference Ngo 2023 deliberately did NOT correct for when extrapolating (Ngo 2023 Discussion, 'limitations'). Ngo 2023 describes only the carbamazepine arm of the original study; Ngo 2020 also studied phenytoin, which is out of scope here."
  )

  ini({
    # ====================================================================
    # Structural PK parameters -- Ngo 2023 Table 1, column 'Estimated in
    # Rats' (originally Ngo 2020). Clearances and volumes are reported per
    # kilogram of body weight and are multiplied by WT in model().
    # ====================================================================
    lcl <- log(0.609); label("Apparent clearance CL/F from the central compartment per kg body weight, without carbamazepine (L/h/kg)")  # Ngo 2023 Table 1: CL/F control group = 0.609 L/h/kg, RSE 36.9 %
    lvc <- log(0.701); label("Apparent central volume of distribution Vc/F per kg body weight (L/kg)")                                   # Ngo 2023 Table 1: Vc/F = 0.701 L/kg, RSE 44.2 %
    lq  <- log(0.665); label("Apparent inter-compartmental clearance Q/F per kg body weight (L/h/kg)")                                   # Ngo 2023 Table 1: Q/F = 0.665 L/h/kg, RSE 21.4 %
    lvp <- log(5.60);  label("Apparent peripheral volume of distribution Vp/F per kg body weight (L/kg)")                                # Ngo 2023 Table 1: Vp/F = 5.60 L/kg, RSE 38.1 %

    lka     <- log(2.31);  label("First-order absorption rate constant Ka of the depot arm (1/h)")                                       # Ngo 2023 Table 1: Ka = 2.31 1/h in rats, RSE 33.6 %; the human sibling model instead uses 0.97 1/h from Mueck 2007
    ld2     <- log(6.62);  label("Zero-order absorption duration D2 into the central compartment without carbamazepine (h)")             # Ngo 2023 Table 1 and Equation 2: D2 = 6.62 h in the control group (no RSE printed)
    ltlag   <- log(0.501); label("Lag time Alag2 before the zero-order absorption arm begins (h)")                                       # Ngo 2023 Table 1: Alag2 = 0.501 h (no RSE printed)
    lfdepot <- log(0.260); label("Fraction F1 of the dose absorbed by the first-order (depot) arm (unitless)")                           # Ngo 2023 Table 1: F1 = 0.260, RSE 9.10 %; the complementary 1 - F1 = 0.740 takes the zero-order route

    # ====================================================================
    # Carbamazepine drug-drug-interaction coefficients. Ngo 2023 codes
    # carbamazepine pretreatment as a 0/1 index variable entering as a
    # FRACTIONAL CHANGE, not as a log-multiplicative shift:
    #     Equation 1 / 3:  CL/F = TVCL * (1 + CL_CBZ)
    #     Equation 2:      D2   = TVD2 * (1 + D2_CBZ)
    # ====================================================================
    e_conmed_cbz_cl <- 2.11;  label("Fractional increase in apparent clearance CL/F with carbamazepine pretreatment (unitless)")             # Ngo 2023 Equation 1 and Section 2.1: CL_CBZ = 2.11 (a 211 % increase); check 0.609 * (1 + 2.11) = 1.894 L/h/kg, the Table 1 test-group value (RSE 69.7 %)
    e_conmed_cbz_d2 <- 0.339; label("Fractional increase in the zero-order absorption duration D2 with carbamazepine pretreatment (unitless)")  # Ngo 2023 Equation 2 and Section 2.1: D2_CBZ = 0.339 (a 33.9 % increase, RSE 17.2 %); 6.62 * (1 + 0.339) = 8.86 h vs the 8.84 h printed in Table 1 -- see vignette Errata

    # ====================================================================
    # Inter-individual variability. Ngo 2023 Table 1 reports IIV as a
    # percentage, read as a coefficient of variation and converted to the
    # log-scale variance by the lognormal identity omega^2 = log(CV^2 + 1).
    # ====================================================================
    etalcl ~ log(0.490^2 + 1)  # Ngo 2023 Table 1: IIV for CL/F = 49.0 % CV (RSE 35.0 %) -> omega^2 = 0.2151
    etalvc ~ log(0.470^2 + 1)  # Ngo 2023 Table 1: IIV for Vc/F = 47.0 % CV (RSE 33.0 %) -> omega^2 = 0.1996

    # ====================================================================
    # Residual unexplained variability -- combined additive and
    # proportional, on the ng/mL scale of the observed plasma
    # concentrations.
    # ====================================================================
    addSd  <- 13.6;  label("Additive residual error on plasma rivaroxaban concentration (ng/mL)")        # Ngo 2023 Table 1: additive error 13.6 ng/mL (RSE 35.8 %)
    propSd <- 0.232; label("Proportional residual error on plasma rivaroxaban concentration (fraction)")  # Ngo 2023 Table 1: proportional error 23.2 % (RSE 17.2 %)
  })

  model({
    # ====================================================================
    # Individual PK parameters. Clearances and volumes are per-kilogram
    # quantities, so each is multiplied by body weight to give the
    # absolute L/h and L that the ODEs below require. The carbamazepine
    # terms are applied on the NATURAL scale as (1 + coefficient *
    # indicator), matching Ngo 2023 Equations 1-3 exactly; writing them as
    # exp(l... + coef * CONMED_CBZ) would be a log-multiplicative form the
    # paper did not use and would not reproduce the printed test-group
    # values.
    # ====================================================================
    cl     <- exp(lcl + etalcl) * WT * (1 + e_conmed_cbz_cl * CONMED_CBZ)
    vc     <- exp(lvc + etalvc) * WT
    q      <- exp(lq) * WT
    vp     <- exp(lvp) * WT
    ka     <- exp(lka)
    d2     <- exp(ld2) * (1 + e_conmed_cbz_d2 * CONMED_CBZ)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot)

    # ====================================================================
    # Two-compartment disposition with parallel mixed-order absorption
    # (Ngo 2023 Figure 2). The oral dose leaves the GI tract by two routes
    # that both terminate in the central compartment:
    #
    #   depot   : fraction F1 = 0.260, absorbed first-order at rate Ka,
    #             representing lymphatic uptake of rivaroxaban carried in
    #             the lipid vehicle. No lag on this arm.
    #   central : fraction 1 - F1 = 0.740, delivered as a zero-order input
    #             of duration D2 beginning after the lag Alag2,
    #             representing solubility-limited transfer through the
    #             capillary network around the enterocytes and the hepatic
    #             portal vein.
    #
    # A user therefore encodes each oral administration as TWO dose
    # records at the same time: one with cmt = "depot" and one with
    # cmt = "central". The f() multipliers split the dose; dur(central)
    # imposes the zero-order duration and alag(central) the lag. Because
    # every disposition parameter is apparent (divided by the true
    # bioavailability F), F1 and 1 - F1 sum to 1 and describe only how the
    # absorbed dose is partitioned between the two routes, not how much of
    # the dose is absorbed.
    # ====================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    f(depot)      <- fdepot
    f(central)    <- 1 - fdepot
    dur(central)  <- d2
    alag(central) <- tlag

    # Plasma rivaroxaban concentration. Doses are in mg and vc is in L, so
    # central / vc is mg/L = ug/mL; the factor 1000 converts to the ng/mL
    # units in which Ngo 2023 reports every concentration, Cmax and AUC.
    Cc <- (central / vc) * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
