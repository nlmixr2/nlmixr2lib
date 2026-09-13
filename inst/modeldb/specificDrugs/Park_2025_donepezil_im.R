Park_2025_donepezil_im <- function() {
  description <- "Two-compartment population PK model for GB-5001A, a long-acting intramuscular donepezil formulation, with three-phase absorption: two lagged parallel first-order depots plus a simultaneous zero-order input into central, in healthy Korean adult men"
  reference <- paste(
    "Park YC, Seol E, Lee J, Hong JH, Jung J-G, Sunwoo J.",
    "Pharmacokinetic Evaluation of GB-5001, a Long-Acting Injectable Formulation",
    "of Donepezil, in Healthy Korean Participants: Population Pharmacokinetics",
    "with Phase 1 Study. Pharmaceutics. 2025;17(12):1517.",
    "doi:10.3390/pharmaceutics17121517.",
    "Structural model Figure 2; parameter estimates Table 6.",
    "The earlier GB-5001 formulation modelled by this group's reference 14 is",
    "distributed as modellib('Khwarg_2024_donepezil_im').",
    sep = " "
  )
  vignette <- "Park_2025_donepezil"
  # Each IM administration is entered as three parallel dose records (Fig. 2:
  # NONMEM depot compartments 3 and 4, plus the zero-order arm into central).
  # Declared explicitly because the registry's default detection only
  # recognises depot and central.
  dosing <- c("depot", "depot2", "central")
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Sex, age, body weight and CYP2D6 metabolizer phenotype were carried in the
  # modelling dataset (Park 2025 Section 2.10) but no covariate effect is
  # retained in the final model: Table 6 reports structural, IIV and residual
  # parameters only. They are documented here rather than in covariateData so
  # the paper's covariate screen keeps its provenance without declaring
  # covariates that model() never references.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Carried in the modelling dataset (Section 2.10); no weight effect retained in the final model (Table 6)."
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Carried in the modelling dataset (Section 2.10); no age effect retained in the final model (Table 6)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Carried in the modelling dataset (Section 2.10) but not estimable:",
        "all 50 enrolled participants were male (Table 1), so SEXF = 0 throughout.",
        sep = " "
      )
    ),
    CYP2D6_EM = list(
      description = "CYP2D6 extensive-metabolizer phenotype indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Phenotyped in Part B only (Cohorts E and F; Section 3.9): 5/8 and 7/8",
        "extensive metabolizers respectively. No CYP2D6 effect retained in the",
        "final model (Table 6); the paper states the sample was too small for a",
        "meaningful comparison (Discussion).",
        sep = " "
      )
    ),
    CYP2D6_IM = list(
      description = "CYP2D6 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Phenotyped in Part B only (Cohorts E and F; Section 3.9): 3/8 and 1/8",
        "intermediate metabolizers respectively. CYP2D6 poor metabolizers were",
        "excluded from Part B by protocol (Section 2.4). No CYP2D6 effect",
        "retained in the final model (Table 6).",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32,
    n_studies = 1,
    age_range = "19-55 years (eligibility); IM cohort means 29.22 (SD 6.82), 30.75 (SD 10.21) and 31.25 (SD 6.02) years at 70, 140 and 280 mg",
    weight_range = "at least 55 kg (eligibility); IM cohort means 74.28 (SD 8.90), 72.69 (SD 13.39) and 70.98 (SD 7.21) kg at 70, 140 and 280 mg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy",
    dose_range = "70, 140 and 280 mg single intramuscular (right ventrogluteal) dose of GB-5001A",
    regions = "Republic of Korea (Clinical Trials Center, Chungnam National University Hospital, Daejeon); NCT06127368; sponsor G2GBIO, Cheongju",
    notes = paste(
      "Healthy Asian men aged 19-55 years with body weight at least 55 kg and BMI",
      "18.5-30.0 kg/m2, from the open-label active-controlled dose-escalation",
      "phase 1 study NCT06127368. CYP2D6 poor metabolizers were excluded from",
      "Part B (Section 2.4).",
      "This model file covers the GB-5001A intramuscular arms only (Cohorts A, E",
      "and F at 70, 140 and 280 mg). Section 3.8 states that 32 participants were",
      "included in the modelling analysis, pooling those IM cohorts with the oral",
      "Aricept Cohort D; the per-cohort PK analysis counts in Table 2 sum to 31",
      "(8 + 8 + 7 IM and 8 oral) and the randomised counts to 33, so the stated",
      "total is not fully reconcilable from the published tables.",
      "The companion oral Aricept 10 mg two-compartment model (Figure S1) is NOT",
      "distributed: Table 6 reports GB-5001A parameters only and the sole oral",
      "estimate published anywhere in the paper is ka = 0.7 1/h (Discussion).",
      "The 6-O-desmethyl donepezil metabolite and the AChE-inhibition",
      "pharmacodynamics were analysed non-compartmentally (Tables S2 and 3) and",
      "have no published structural model.",
      sep = " "
    )
  )

  ini({
    # Absorption -- Park 2025 Table 6. The IM dose is split three ways (Fig. 2):
    # fraction F3 into a lagged first-order depot, fraction F4 into a second
    # lagged first-order depot, and the remaining 1 - F3 - F4 = 0.148 as a
    # zero-order input delivered directly into the central compartment over D1
    # starting at the time of dose. F1 is named in the Figure 2 caption as the
    # fraction absorbed by the zero-order pathway but is not tabulated; it is
    # the complement of the two estimated depot fractions, which the paper's own
    # steady-state simulation confirms (see the vignette source-trace table:
    # AUCtau = dose / CL requires the three fractions to sum to one).
    lka <- log(0.00173); label("Absorption rate constant from the first IM depot (KA3, 1/h)") # Table 6: KA3 = 0.00173 1/h (SE 0.000141; 95% CI 0.00145-0.00201)
    lka2 <- log(0.0113); label("Absorption rate constant from the second IM depot (KA4, 1/h)") # Table 6: KA4 = 0.0113 1/h (SE 0.00481; 95% CI 0.00187-0.0207)
    lfdepot <- log(0.657); label("Fraction of the IM dose entering the first depot (F3, fraction)") # Table 6: F3 = 0.657 (SE 0.0172; 95% CI 0.623-0.691)
    lfdepot2 <- log(0.195); label("Fraction of the IM dose entering the second depot (F4, fraction)") # Table 6: F4 = 0.195 (SE 0.0142; 95% CI 0.167-0.223)
    ltlag <- log(299); label("Absorption lag time of the first IM depot (ALAG3, h)") # Table 6: ALAG3 = 299 h (SE 6.97; 95% CI 285.456-312.544)
    ltlag2 <- log(1130); label("Absorption lag time of the second IM depot (ALAG4, h)") # Table 6: ALAG4 = 1130 h (SE 27.7; 95% CI 1075.708-1184.292)
    ld1 <- log(497); label("Duration of the zero-order absorption input into central (D1, h)") # Table 6: D1 = 497 h (SE 5.9; 95% CI 485.436-508.564)

    # Disposition -- Park 2025 Table 6. Q was fixed at 100 L/h in the final model
    # (Table 6 'Q  100 FIX'); Section 2.10 explains that estimating Q was
    # unstable because the flip-flop kinetics of the formulation plus sparse
    # early sampling leave the distribution phase uncharacterised, and that 100
    # L/h was chosen from preliminary fits of this dataset with reference to the
    # Q = 185 L/h of the earlier formulation (reference 14).
    lcl <- log(9.57); label("Clearance (CL, L/h)") # Table 6: CL = 9.57 L/h (SE 0.57; 95% CI 8.453-10.687)
    lvc <- log(58.2); label("Central volume of distribution (V1, L)") # Table 6: V1 = 58.2 L (SE 47.9; 95% CI -35.684-152.084, i.e. imprecise; see Section 3.8)
    lq <- fixed(log(100)); label("Inter-compartmental clearance (Q, L/h)") # Table 6: Q = 100 L/h FIXED (Section 2.10)
    lvp <- log(2050); label("Peripheral volume of distribution (V2, L)") # Table 6: V2 = 2050 L (SE 303; 95% CI 1456.12-2643.88)

    # IIV -- Park 2025 Table 6 reports the log-scale variance directly; the
    # Table 6 footnote gives the derived CV as CV(%) = sqrt(exp(omega) - 1) * 100.
    # Values below are those variances: sqrt(exp(0.0945) - 1) = 31.5% for CL and
    # sqrt(exp(6.97) - 1) = 3264% for V1. No covariance between the two etas is
    # reported, so they are encoded as independent. The 31.5% CL CV is confirmed
    # exactly by the paper's own steady-state simulation, which reports 31.5% and
    # 31.6% CV on AUCtau and Cav,ss at 140 and 280 mg (Section 3.8) -- both of
    # which scale as 1/CL and so carry the CL CV unchanged.
    etalcl ~ 0.0945 # Table 6: IIV CL = 0.0945 (SE 0.0334), a log-scale variance
    etalvc ~ 6.97 # Table 6: IIV V1 = 6.97 (SE 1.69), a log-scale variance; very imprecise, matching the V1 point estimate whose 95% CI spans zero

    # Residual error -- Table 6 reports the additive term in linear ng/mL and
    # labels the proportional term 'a proportional residual error represented as
    # CV', so both are standard deviations rather than NONMEM $SIGMA variances.
    # Equation (1) as typeset (Y = F + sqrt(F^2 * theta1^2 * theta2^2 * eps)) is
    # not a usable error model -- it takes the square root of a signed normal
    # deviate and multiplies rather than adds the two components. Section 2.10's
    # prose ('a single residual error term (eps) ... with weights assigned based
    # on a proportional component (theta1) and an additive component (theta2),
    # thereby implementing a combined error structure') describes the standard
    # NONMEM combined form W = sqrt((theta1 * F)^2 + theta2^2), Y = F + W * eps,
    # which is what is encoded here.
    propSd <- 0.271; label("Proportional residual error (fraction)") # Table 6: eps2 (proportional) = 0.271 as CV (SE 0.0174; 95% CI 0.237-0.305)
    addSd <- 0.0713; label("Additive residual error (ug/L)") # Table 6: eps1 (additive) = 0.0713 ng/mL = 0.0713 ug/L (SE 0.0207; 95% CI 0.0307-0.112)
  })

  model({
    # Individual parameters
    ka <- exp(lka)
    ka2 <- exp(lka2)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)
    fdepot2 <- exp(lfdepot2)
    tlag <- exp(ltlag)
    tlag2 <- exp(ltlag2)
    d1 <- exp(ld1)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Three-phase IM absorption into a two-compartment disposition model
    # (Park 2025 Fig. 2). Each IM administration is encoded by the USER as
    # THREE parallel dose records in the event table, all with the same amt:
    #   (a) cmt = "depot"   -- first-order arm, fraction F3, lag ALAG3
    #   (b) cmt = "depot2"  -- first-order arm, fraction F4, lag ALAG4
    #   (c) cmt = "central", rate = -2 -- zero-order arm, fraction 1 - F3 - F4,
    #       modelled duration D1 supplied by dur(central) below.
    # The f() multipliers perform the dose split, so each record carries the
    # full nominal dose. See the vignette for a worked event table.
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag
    f(depot2) <- fdepot2
    alag(depot2) <- tlag2
    f(central) <- 1 - fdepot - fdepot2
    dur(central) <- d1

    # Plasma donepezil concentration. Dose is in mg and vc in L, so central/vc is
    # mg/L; the factor 1000 converts to the ug/L (= ng/mL) scale used in
    # Park 2025 Tables 2 and 6 and for the additive residual error.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
