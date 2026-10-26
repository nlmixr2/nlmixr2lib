Purohit_2023_ritlecitinib <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption for oral ritlecitinib (PF-06651600, a covalent JAK3/TEC-family kinase inhibitor) in healthy participants and patients with rheumatoid arthritis. Apparent clearance and apparent central volume carry allometric body-weight scaling with exponents fixed at 0.75 and 1 referenced to 75 kg; apparent clearance is lower in rheumatoid arthritis patients and at total daily doses above 100 mg, and a high-fat meal slows absorption. Interindividual variability on CL/F and Vc/F is correlated, and the log-scale residual error differs between healthy participants and rheumatoid arthritis patients. This is the interim (fit-for-purpose) model used to simulate in silico healthy-participant reference cohorts for the ritlecitinib hepatic and renal impairment studies."
  reference <- paste(
    "Purohit V, Huh Y, Wojciechowski J, Plotka A, Salts S, Antinew J,",
    "Dimitrova A, Nicholas T.",
    "Leveraging Prior Healthy Participant Pharmacokinetic Data to Evaluate the",
    "Impact of Renal and Hepatic Impairment on Ritlecitinib Pharmacokinetics.",
    "AAPS J. 2023;25(3):32. doi:10.1208/s12248-023-00792-8. PMCID PMC10047454.",
    "Parameter estimates from Supplementary Table B; model equations from the",
    "final-model NONMEM control stream (run12.mod) in Supplementary Text S1.",
    sep = " "
  )
  vignette <- "Purohit_2023_ritlecitinib"
  paper_specific_residual_sds <- c("expSdHv", "expSdRa")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling referenced to 75 kg with exponents fixed at 0.75 on",
        "CL/F and 1 on Vc/F (control stream $PK: TVCL = ... * (BWT/75)**0.75 and",
        "TVV2 = THETA(3) * (BWT/75)). Supplementary Text S1 'Structural Model",
        "Development' says Q/F and Vp/F were also scaled, but the final-model",
        "control stream and the 'Final Model' paragraph scale only CL/F and Vc/F;",
        "the control stream is followed here."
      ),
      source_name = "BWT"
    ),
    DIS_RA = list(
      description = "Rheumatoid arthritis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant)",
      notes = paste(
        "Derived from the NONMEM PROT column: FLG2 = 0 for PROT = 1006, the",
        "phase II rheumatoid arthritis study NCT02969044, and 1 for the two phase",
        "I healthy-participant studies. DIS_RA = 1 - FLG2. Scales CL/F by 0.625",
        "(Supplementary Table B 'RA on CL/F') and selects the RA residual-error",
        "magnitude (control stream $ERROR uses EPS(2) when FLG2 = 0)."
      ),
      source_name = "PROT"
    ),
    FED_HIGHFAT = list(
      description = "High-fat meal at dosing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, or food status not recorded)",
      notes = paste(
        "Derived from the NONMEM FOOD column: FLG3 = 0 when FOOD is non-zero,",
        "FED_HIGHFAT = 1 - FLG3. Supplementary Text S1 'Covariate Model",
        "Development': 'The effect of food on the first-order absorption rate",
        "constant (ka) was assessed in participants who had high-fat meals",
        "recorded. Participants where food status was not captured were grouped",
        "with the fasted participants.' The fed records come from the high-fat",
        "meal arm of the relative-bioavailability study NCT02684760. Multiplies",
        "ka by (1 - 0.613)."
      ),
      source_name = "FOOD"
    ),
    DOSE_HIGH = list(
      description = "Total-daily-dose-above-100-mg indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (total daily dose of 100 mg or less)",
      notes = paste(
        "Control stream $PK: TRT = DOSE*DFREQ (total daily dose, mg);",
        "FLG1 = 0 IF (TRT.GT.100), and DOSE_HIGH = 1 - FLG1. A 100 mg BID",
        "regimen therefore has DOSE_HIGH = 1 and a 100 mg QD regimen has",
        "DOSE_HIGH = 0. Multiplies CL/F by 0.666 (Supplementary Table B 'Dose on",
        "CL/F'); Supplementary Text S1 attributes it to apparent non-linear CL/F",
        "seen in noncompartmental analyses. In the analysis data the >100 mg/day",
        "regimens were the 200-800 mg single doses, the 200 and 400 mg QD and",
        "100 and 200 mg BID multiple-dose arms, and the 200 mg QD RA arm."
      ),
      source_name = "DOSE, DFREQ"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Evaluated graphically against the individual random effects",
        "(Supplementary Text S1 'Covariate Model Development'); 'There was no",
        "apparent difference in clearance with respect to age.' Not in the final",
        "model."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ritlecitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ritlecitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ritlecitinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 113L,
    n_studies = 3L,
    n_observations = 1519L,
    disease_state = "Healthy participants (two phase I studies) and patients with moderate to severe active rheumatoid arthritis with an inadequate response to methotrexate (one phase II study).",
    dose_range = "Single oral doses 5-800 mg; 50-400 mg QD or 100-200 mg BID for 14 days; 50 mg tablet or oral solution fasted and tablet after a high-fat meal; 200 mg QD for 8 weeks in rheumatoid arthritis.",
    studies = c(
      "NCT02309827 (phase I SAD/MAD, healthy participants, n = 80 enrolled)",
      "NCT02684760 (phase I relative bioavailability and food effect, healthy participants, n = 14)",
      "NCT02969044 (phase II, rheumatoid arthritis, n = 70 randomized including placebo)"
    ),
    notes = paste(
      "Supplementary Text S1 and Table A. 1519 ritlecitinib concentrations from",
      "113 participants. Baseline demographics of the modelling population are",
      "not reported; the enrolment counts above include placebo recipients, who",
      "contribute no concentrations."
    )
  )

  ini({
    # Structural parameters (typical values at 75 kg, healthy participant,
    # fasted, total daily dose <= 100 mg). All from Supplementary Table B,
    # which reports linear-scale estimates; the control stream estimates the
    # THETAs on the linear scale (ADVAN4 TRANS4).
    lka <- log(3.48); label("First-order absorption rate constant (1/h)") # Table B 'ka' = 3.48 (95% CI 3.09-3.87); THETA(1)
    lcl <- log(78.1); label("Apparent clearance CL/F at 75 kg (L/h)") # Table B 'CL/F' = 78.1 (70.3-86.0); THETA(2)
    lvc <- log(119); label("Apparent central volume Vc/F at 75 kg (L)") # Table B 'V/F' = 119 (99.7-138); THETA(3)
    lq <- log(0.424); label("Apparent intercompartmental clearance Q/F (L/h)") # Table B 'Q' = 0.424 (0.339-0.509); THETA(4)
    lvp <- log(3.54); label("Apparent peripheral volume Vp/F (L)") # Table B 'Vp/F' = 3.54 (3.12-3.97); THETA(5)

    # Allometric exponents, not estimated (Supplementary Text S1 and the
    # main-text Population Pharmacokinetic Model paragraph).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Text S1: exponents 'fixed to 0.75 and 1'; control stream (BWT/75)**0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F (unitless)") # Text S1; control stream TVV2 = THETA(3)*(BWT/75)

    # Covariate effects. The control stream applies the dose and RA effects
    # as ratio multipliers on CL/F (THETA(7) and THETA(10) replace a factor of
    # 1); they are carried here on the log scale so exp(e * indicator) gives
    # the printed ratio for indicator = 1 and 1 for indicator = 0.
    e_dose_high_cl <- log(0.666); label("Log ratio of CL/F at total daily dose above 100 mg (unitless)") # Table B 'Dose on CL/F' = 0.666 (0.641-0.691); THETA(7)
    e_dis_ra_cl <- log(0.625); label("Log ratio of CL/F in rheumatoid arthritis patients (unitless)") # Table B 'RA on CL/F' = 0.625 (0.564-0.739); THETA(10)
    e_fed_highfat_ka <- -0.613; label("Fractional change in ka with a high-fat meal (unitless)") # Table B 'Food on ka' = -0.613 (-0.649 to -0.577); THETA(9), TVKA = THETA(1)*(1 + FED*THETA(9))

    # IIV: $OMEGA BLOCK(2). Table B's footnote calls these rows variances,
    # but they are read as NONMEM's correlation-matrix print: SDs on the
    # diagonal (0.156, 0.280) and a correlation off it (0.173). Only this
    # reading reproduces the spread of the paper's own Table IV simulation
    # (taken as variances, the trial-level AUC0-24 90% PI is ~2.2-fold too
    # wide), and it matches the healthy-participant IIV SD on CL/F of 0.201
    # in the Wojciechowski_2023_ritlecitinib_base sibling model.
    # Variances: 0.156^2 = 0.024336, 0.280^2 = 0.0784;
    # covariance: 0.173 * 0.156 * 0.280 = 0.00755664.
    etalcl + etalvc ~ c(
      0.024336,
      0.00755664, 0.0784
    ) # Table B 'IIV-CL/F' = 0.156, 'Cov CL/F-V/F' = 0.173, 'IIV-V/F' = 0.280, read as SD / correlation / SD

    # Residual error. Concentrations were log-transformed both sides,
    # Y = LOG(F) + EPS, with EPS(1) for healthy participants and EPS(2) for
    # RA patients. These rows are read as the labelled $SIGMA variances
    # (initial values 0.15 and 0.7 in $SIGMA), so the nlmixr2 lnorm() SD is
    # their square root; as SDs the simulated Table IV Cmax medians fall
    # about 15 percent short.
    expSdHv <- sqrt(0.166); label("Log-scale residual SD in healthy participants") # Table B 'RUV-HP' = 0.166 (variance); SD = 0.407
    expSdRa <- sqrt(0.860); label("Log-scale residual SD in rheumatoid arthritis patients") # Table B 'RUV-RA' = 0.860 (variance); SD = 0.927
  })

  model({
    # Individual parameters (control stream $PK).
    ka <- exp(lka) * (1 + e_fed_highfat_ka * FED_HIGHFAT)
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl *
      exp(e_dose_high_cl * DOSE_HIGH) * exp(e_dis_ra_cl * DIS_RA)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    # The control stream carries an absorption lag ALAG1 = THETA(6) +
    # FED*THETA(8) with both THETAs set to 0 and not estimated, so there is
    # no lag and none is coded here.

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; S2 = V2/1000 in the control
    # stream reports ng/mL.
    Cc <- 1000 * central / vc

    expSdPop <- expSdHv * (1 - DIS_RA) + expSdRa * DIS_RA
    Cc ~ lnorm(expSdPop)
  })
}
