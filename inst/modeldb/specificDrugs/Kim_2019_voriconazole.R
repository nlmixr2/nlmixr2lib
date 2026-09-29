Kim_2019_voriconazole <- function() {
  description <- "Three-compartment population pharmacokinetic model for intravenous and oral voriconazole in Korean healthy volunteers and patients (Kim 2019), with first-order absorption, an absorption lag time, a logit-scale absolute bioavailability, and a mechanistic auto-inhibition compartment that makes clearance time-dependent: CL(t) = CL0 * [RCLF + (1 - RCLF) * IC50 / (IC50 + Cinh)], so clearance decays from CL0 toward the non-inhibitable floor CL0 * RCLF as drug accumulates in the inhibition compartment. CYP2C19 metabolizer phenotype acts on both CL0 and RCLF, body weight on CL0, the first peripheral volume and the first inter-compartmental clearance, and a CTCAE grade >= 3 liver function abnormality on CL0."
  reference <- paste(
    "Kim Y, Rhee SJ, Park WB, Yu KS, Jang IJ, Lee S.",
    "A personalized CYP2C19 phenotype-guided dosing regimen of",
    "voriconazole using a population pharmacokinetic analysis.",
    "J Clin Med. 2019 Feb 10;8(2):227. doi:10.3390/jcm8020227.",
    "PMCID PMC6406770.",
    sep = " "
  )
  vignette <- "Kim_2019_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "voriconazole", units = "mg/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight; power-function covariate on CL0, on the first peripheral volume of distribution (V3) and on the first inter-compartmental clearance (Q2)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline body weight. Kim 2019 Methods 'Population pharmacokinetic analysis' states that continuous covariates were tested 'using power functions normalized to their median values or generally accepted typical value (e.g., 70 kg for body weight)', so the reference weight is 70 kg. Exponents (Kim 2019 Table 3): 0.595 on CL, 2.2 on V3, 2.56 on Q2. The V3 and Q2 exponents are far above the allometric 0.75 / 1 range; they are applied exactly as published.",
      source_name = "BW"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator; exponential effect on CL0 and on the non-inhibitable clearance fraction",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer, the EM reference group)",
      notes = "1 = CYP2C19 intermediate metabolizer. Kim 2019 Methods 'Study Population' classifies phenotypes per the CPIC guideline: EM = *1/*1; IM = *1/*2, *1/*3, *2/*17; PM = *2/*2, *2/*3, *3/*3. One subject genotyped as a rapid metabolizer (*1/*17) was analysed as an EM. Paired with CYP2C19_PM; both 0 is the EM reference.",
      source_name = "CYP2C19 phenotype (IM)"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator; exponential effect on CL0 and on the non-inhibitable clearance fraction",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer, the EM reference group)",
      notes = "1 = CYP2C19 poor metabolizer (*2/*2, *2/*3, *3/*3 per Kim 2019 Methods 'Study Population'). Paired with CYP2C19_IM; both 0 is the EM reference.",
      source_name = "CYP2C19 phenotype (PM)"
    ),
    HEPIMP_SEV = list(
      description = "Severe liver function abnormality indicator; exponential effect on CL0",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CTCAE liver function abnormality grade < 3)",
      notes = "Classification scheme is the Common Terminology Criteria for Adverse Events version 4.0 liver function abnormality grade (Kim 2019 Methods 'Population pharmacokinetic analysis' and Table 1), NOT the NCI ODWG or Child-Pugh schemes that other HEPIMP_SEV models use. 1 = CTCAE grade >= 3; 0 = grade < 3. Six of 193 subjects (grade 3 n = 5, grade 4 n = 1) were grade >= 3 (Kim 2019 Table 1), so the effect is supported by a small subgroup; Kim 2019 Discussion notes that 'further evaluation is needed due to limited data available from only 6 patients with grade >= 3 hepatic abnormality'.",
      source_name = "Liver function abnormality grade (CTCAE v4.0)"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator; selects which of the two published additive residual-error magnitudes applies to an observation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient)",
      notes = "1 = healthy volunteer (Kim 2019 studies 1-4, intensively sampled), 0 = patient (study 5, sparsely sampled TDM data). Kim 2019 Methods 'Population pharmacokinetic analysis' states that residual error models were 'tested for the healthy volunteers' and the patients' data independently', and Table 2 reports separate additive errors of 0.208 mg/L (healthy subjects) and 0.799 mg/L (patients). Time-fixed per subject.",
      source_name = "study population (healthy subjects vs patients)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 193L,
    n_studies = 5L,
    age_range = "18-80 years",
    age_median = "34 years (healthy subjects 26, patients 59)",
    weight_range = "40.8-88.5 kg",
    weight_median = "66.0 kg (healthy subjects 70.3, patients 59.4)",
    sex_female_pct = 15.0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Pooled analysis of 93 Korean healthy volunteers (studies 1-4, intensively sampled) and 100 Korean patients receiving voriconazole for suspected or proven invasive fungal infection under therapeutic drug monitoring (study 5, sparsely sampled). CYP2C19 phenotype: extensive metabolizer 75 (39%), intermediate metabolizer 70 (36%), poor metabolizer 48 (25%). CTCAE v4.0 liver function abnormality grade 0 in 165 (85.5%), grade 1 in 9, grade 2 in 13, grade 3 in 5, grade 4 in 1. Co-medications: proton pump inhibitors 22 (11.4%), steroids 9 (4.7%); neither was a significant covariate.",
    dose_range = "Study 1: 200 mg IV single dose, then single and multiple oral 200 mg q12h. Study 2: 400 mg oral single dose. Study 3: 200 mg IV single dose. Study 4: 200 mg IV single dose then 200 mg oral single dose. Study 5: loading 6 mg/kg IV or 400 mg oral q12h on day 1, then TDM-based maintenance 4 mg/kg IV or 200 mg oral q12h.",
    regions = "Republic of Korea (all five studies conducted at Seoul National University Hospital)",
    lab_range = "Aspartate aminotransferase 21 (7-377) U/L; alanine aminotransferase 21 (4-363) U/L (median (range), Kim 2019 Table 1)",
    notes = "1,828 plasma voriconazole concentrations (1,579 from healthy volunteers, 249 from patients). Demographics from Kim 2019 Table 1 and Table S1; values are median (range). NONMEM 7.3.0 FOCE-I on log-transformed concentrations; the model was built on the healthy-volunteer data and then extended with the patient data."
  )

  # Implementation notes (the vignette 'Assumptions and deviations' section
  # carries the full justification and the arbitration evidence for each):
  #
  # * Auto-inhibition structure. Kim 2019 Equations 1 and 2 give
  #     CL  = CL0 * INH
  #     INH = RCLF + (1 - RCLF) * (1 - Cinh / (IC50 + Cinh))
  #   which is algebraically
  #     INH = RCLF + (1 - RCLF) * IC50 / (IC50 + Cinh),
  #   i.e. INH = 1 when Cinh = 0 and INH -> RCLF as Cinh grows. RCLF is
  #   therefore the fraction of CL0 that cannot be inhibited, and the
  #   equivalent Imax form has Imax = 1 - RCLF. The packaged model uses the
  #   IC50 / (IC50 + Cinh) spelling because it is numerically identical and
  #   avoids the cancelling 1 - x/(a+x) round trip.
  #
  # * Inhibition-compartment ODE. Kim 2019 does not print the differential
  #   equation for the inhibition compartment. Figure 1 draws a single
  #   dashed arrow from the central compartment into the inhibition
  #   compartment labelled KIC and no return or elimination arrow, and
  #   Table 2 names KIC the "rate constant INTO inhibition compartment" and
  #   reports no second rate constant. The model therefore uses the
  #   degenerate pure-integrator effect-compartment form
  #     d/dt(effect) <- ke0 * Cc
  #   which is the form the canonical `ke0` register entry
  #   (inst/references/parameter-names.md, "Degenerate pure-integrator
  #   form") directs extractions to encode under the `ke0` name rather than
  #   minting a new canonical. The alternative Sheiner equilibration form
  #   ke0 * (Cc - effect) was simulated against Supplementary Figure S1 and
  #   differs by at most ~1% in CL over the published 72 h window, because
  #   ke0 * t = 0.14 there; the choice is not load-bearing. The vignette
  #   shows both.
  #
  # * IC50 is fixed at 0.01 mg/L (Table 2, "0.01 FIX"). Under the
  #   pure-integrator form Cinh = ke0 * AUC, so ke0 and IC50 enter only
  #   through their ratio and one of the two must be fixed for the other to
  #   be identifiable; fixing IC50 is what makes the reported ke0 (RSE
  #   14.9%) estimable.
  #
  # * Bioavailability. Kim 2019 Methods states F was "estimated using a
  #   logit model", so F is carried as logitfdepot with logit(0.876) =
  #   1.955 and the reported 84.4% CV IIV is applied on the logit scale.
  #   An exponential IIV of that size on the natural scale would put a
  #   large share of subjects above F = 1.
  #
  # * CYP2C19 effect on RCLF -- Table 3 versus the Discussion. Table 3
  #   gives the RCLF effects as IM = -0.51 and PM = -0.44, so
  #   RCLF_IM = 0.162 * exp(-0.51) = 0.0973 and
  #   RCLF_PM = 0.162 * exp(-0.44) = 0.1043, matching the "decreased by
  #   approximately 36-40% (0.097-0.104)" range printed in Results 3.2 and
  #   in the Discussion. The Discussion additionally states that the
  #   time-dependent CL "changes from 45.3 to 7.3 L/h in CYP2C19 EMs, from
  #   37.6 to 3.9 L/h in IMs, and from 21.5 to 2.1 L/h in PMs". The EM pair
  #   is exactly CL0 * RCLF = 45.3 * 0.162 = 7.34, but the IM and PM values
  #   are only reproduced by pairing each CL0 with the OTHER phenotype's
  #   RCLF (37.61 * 0.1043 = 3.92 and 21.48 * 0.0973 = 2.09); Table 3 as
  #   printed gives 3.66 and 2.24. Supplementary Figure S1 arbitrates in
  #   favour of Table 3: digitising its three median CL-versus-time curves
  #   and simulating this model at the patients' median weight reproduces
  #   the IM and PM curves to 2-3% RMS across 3-60 h under Table 3 as
  #   printed, whereas the swapped assignment splits the bias between the
  #   two phenotypes (IM +15%, PM +6%). The mechanism is unambiguous: PM
  #   subjects have the higher exposure, hence the larger Cinh, hence the
  #   inhibition factor closer to its floor, which only the Table 3
  #   assignment can reconcile with the measured curves. Table 3 is used;
  #   the Discussion sentence transposes IM and PM.
  #
  # * Liver function effect wording. Table 3 gives -0.75 on CL, i.e.
  #   exp(-0.75) = 0.472, so CL falls to 47.2% of its value (a 52.8%
  #   reduction). Kim 2019 Results 3.2 and Discussion describe this as a
  #   "47% reduction"; the 47% figure is the remaining fraction, not the
  #   reduction. The coefficient is used as printed.
  #
  # * IIV scale. Kim 2019 Table 2 reports inter-individual variability as
  #   %CV with the footnote "Standard error given on the variance scale",
  #   and reports the off-diagonal OMEGA elements as covariances (footnote
  #   "Standard error of the covariance estimate"). The packaged model uses
  #   the exact log-normal conversion omega^2 = log(CV^2 + 1) for the
  #   diagonal and the published covariances verbatim for the off-diagonal,
  #   which gives a positive-definite 4x4 block (minimum eigenvalue 0.005).
  #
  # * Residual error. Table 2 reports two additive errors, 0.208 mg/L for
  #   healthy subjects and 0.799 mg/L for patients, fitted independently.
  #   Both are carried and selected per observation by DIS_HEALTHY.
  paper_specific_residual_sds <- c("addSdHealthy", "addSdPatient")

  ini({
    # ----- Absorption (Kim 2019 Table 2) -----
    lka <- log(1.23)
    label("Absorption rate constant Ka (1/h)") # Kim 2019 Table 2: Ka = 1.23 1/h (RSE 15.4%)
    ltlag <- log(0.237)
    label("Absorption lag time ALAG1 (h)") # Kim 2019 Table 2: ALAG1 = 0.237 h (RSE 1.8%)
    logitfdepot <- 1.955085
    label("Logit of absolute oral bioavailability F1 (unitless)") # Kim 2019 Table 2: F1 = 0.876 (RSE 2.3%), estimated with a logit model; logit(0.876) = 1.955

    # ----- Disposition (Kim 2019 Table 2) -----
    lvc <- log(35.7)
    label("Central volume of distribution V2 at 70 kg (L)") # Kim 2019 Table 2: V2 = 35.7 L (RSE 15.7%)
    lcl <- log(45.3)
    label("Uninhibited clearance CL0 at 70 kg, CYP2C19 EM (L/h)") # Kim 2019 Table 2: CL = 45.3 L/h (RSE 5.8%); the instantaneous clearance is cl * inh
    lvp <- log(58.9)
    label("First peripheral volume of distribution V3 at 70 kg (L)") # Kim 2019 Table 2: V3 = 58.9 L (RSE 6.2%)
    lq <- log(10.9)
    label("Inter-compartmental clearance Q2 to peripheral1 at 70 kg (L/h)") # Kim 2019 Table 2: Q2 = 10.9 L/h (RSE 8.0%)
    lvp2 <- log(25.4)
    label("Second peripheral volume of distribution V4 (L)") # Kim 2019 Table 2: V4 = 25.4 L (RSE 16.7%)
    lq2 <- log(54.6)
    label("Inter-compartmental clearance Q3 to peripheral2 (L/h)") # Kim 2019 Table 2: Q3 = 54.6 L/h (RSE 45.4%)

    # ----- Auto-inhibition of clearance (Kim 2019 Table 2, Equations 1-2) -----
    lfcl_noinh <- log(0.162)
    label("Fraction of CL0 that cannot be inhibited, RCLF, CYP2C19 EM (unitless)") # Kim 2019 Table 2: RCLF = 0.162 (RSE 9.7%)
    lic50 <- fixed(log(0.01))
    label("Inhibition-compartment concentration giving 50% of maximum CL inhibition, IC50 (mg/L)") # Kim 2019 Table 2: IC50 = 0.01 FIX
    lke0 <- log(0.002)
    label("Rate constant into the inhibition compartment, KIC (1/h)") # Kim 2019 Table 2: KIC = 0.002 1/h (RSE 14.9%)

    # ----- Covariate effects (Kim 2019 Table 3) -----
    e_wt_cl <- 0.595
    label("Body-weight power exponent on CL0 (unitless)") # Kim 2019 Table 3: body weight exponent for CL = 0.595 (RSE 31.8%)
    e_wt_vp <- 2.2
    label("Body-weight power exponent on V3 (unitless)") # Kim 2019 Table 3: body weight exponent for V3 = 2.2 (RSE 20.0%)
    e_wt_q <- 2.56
    label("Body-weight power exponent on Q2 (unitless)") # Kim 2019 Table 3: body weight exponent for Q2 = 2.56 (RSE 18.1%)
    e_cyp2c19_im_cl <- -0.186
    label("Exponential CYP2C19 intermediate-metabolizer effect on CL0 (unitless)") # Kim 2019 Table 3: IM effect for CL = -0.186 (RSE 29.5%); 45.3 * exp(-0.186) = 37.6 L/h
    e_cyp2c19_pm_cl <- -0.746
    label("Exponential CYP2C19 poor-metabolizer effect on CL0 (unitless)") # Kim 2019 Table 3: PM effect for CL = -0.746 (RSE 10.9%); 45.3 * exp(-0.746) = 21.5 L/h
    e_hepimp_sev_cl <- -0.75
    label("Exponential CTCAE grade >= 3 liver function abnormality effect on CL0 (unitless)") # Kim 2019 Table 3: liver function grade >= 3 effect for CL = -0.75 (RSE 49.3%); exp(-0.75) = 0.472
    e_cyp2c19_im_fcl_noinh <- -0.51
    label("Exponential CYP2C19 intermediate-metabolizer effect on RCLF (unitless)") # Kim 2019 Table 3, 'Effect on RCLF': IM = -0.51 (RSE 27.5%); 0.162 * exp(-0.51) = 0.0973
    e_cyp2c19_pm_fcl_noinh <- -0.44
    label("Exponential CYP2C19 poor-metabolizer effect on RCLF (unitless)") # Kim 2019 Table 3, 'Effect on RCLF': PM = -0.44 (RSE 42.3%); 0.162 * exp(-0.44) = 0.1043

    # ----- Inter-individual variability (Kim 2019 Table 2) -----
    # Diagonal: exact log-normal conversion omega^2 = log(CV^2 + 1) from the
    # published %CV. Off-diagonal: the published covariance estimates, used
    # verbatim (Table 2 footnote b, 'Standard error of the covariance
    # estimate'). Block order V2, CL, V3, Q2 matches the Table 2 listing.
    etalvc + etalcl + etalvp + etalq ~ c(
      0.149802,
      0.0116, 0.044778,
      -0.0117, -0.0119, 0.041560,
      -0.0734, 0.008, 0.0345, 0.079683
    ) # Table 2: CV 40.2 / 21.4 / 20.6 / 28.8 %; covariances V2-CL 0.0116, V2-V3 -0.0117, V2-Q2 -0.0734, CL-V3 -0.0119, CL-Q2 0.008, V3-Q2 0.0345
    etalka ~ 0.571479 # Table 2: 'IIV for Ka (% CV)' = 87.8; log(1 + 0.878^2) = 0.571479
    etalogitfdepot ~ 0.537859 # Table 2: 'IIV for F1 (% CV)' = 84.4, applied on the logit scale; log(1 + 0.844^2) = 0.537859
    etalfcl_noinh ~ 0.259233 # Table 2: 'IIV for RCLF (% CV)' = 54.4; log(1 + 0.544^2) = 0.259233

    # ----- Residual error (Kim 2019 Table 2) -----
    addSdHealthy <- 0.208
    label("Additive residual SD, healthy volunteers (mg/L)") # Kim 2019 Table 2: additive error for healthy subjects = 0.208 mg/L (RSE 8.4%)
    addSdPatient <- 0.799
    label("Additive residual SD, patients (mg/L)") # Kim 2019 Table 2: additive error for patients = 0.799 mg/L (RSE 6.7%)
  })
  model({
    # ----- Reference weight -----
    # Kim 2019 Methods: continuous covariates entered as power functions
    # normalized to a "generally accepted typical value (e.g., 70 kg for
    # body weight)".
    ref_wt <- 70

    # ----- Individual parameters -----
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp) * (WT / ref_wt)^e_wt_vp
    q <- exp(lq + etalq) * (WT / ref_wt)^e_wt_q
    vp2 <- exp(lvp2)
    q2 <- exp(lq2)

    # Uninhibited clearance CL0: body weight power function, exponential
    # CYP2C19 phenotype effects (EM = reference) and an exponential effect
    # of a CTCAE grade >= 3 liver function abnormality.
    cl0 <-
      exp(lcl + etalcl) *
      (WT / ref_wt)^e_wt_cl *
      exp(
        e_cyp2c19_im_cl * CYP2C19_IM +
          e_cyp2c19_pm_cl * CYP2C19_PM +
          e_hepimp_sev_cl * HEPIMP_SEV
      )

    # Non-inhibitable clearance fraction RCLF: exponential CYP2C19
    # phenotype effects on the same reference (EM).
    fcl_noinh <-
      exp(lfcl_noinh + etalfcl_noinh) *
      exp(
        e_cyp2c19_im_fcl_noinh * CYP2C19_IM +
          e_cyp2c19_pm_fcl_noinh * CYP2C19_PM
      )
    ic50 <- exp(lic50)
    ke0 <- exp(lke0)

    # ----- Plasma concentration and the inhibition factor -----
    Cc <- central / vc
    # Kim 2019 Equation 2, written in its algebraically identical
    # IC50 / (IC50 + Cinh) form: inh = 1 when effect = 0 and inh ->
    # fcl_noinh as effect grows.
    inh <- fcl_noinh + (1 - fcl_noinh) * ic50 / (ic50 + effect)
    cl <- cl0 * inh

    # ----- ODE system -----
    d/dt(depot) <- -ka * depot
    d/dt(central) <-
      ka * depot -
      cl * Cc -
      q * (Cc - peripheral1 / vp) -
      q2 * (Cc - peripheral2 / vp2)
    d/dt(peripheral1) <- q * (Cc - peripheral1 / vp)
    d/dt(peripheral2) <- q2 * (Cc - peripheral2 / vp2)
    # Inhibition compartment, degenerate pure-integrator effect-compartment
    # form (Kim 2019 Figure 1: a single KIC arrow into the compartment, no
    # return or elimination arrow; Table 2 lists no second rate constant).
    # The state holds the driving concentration Cinh in mg/L.
    d/dt(effect) <- ke0 * Cc

    alag(depot) <- tlag
    f(depot) <- fdepot

    # ----- Observation and population-specific residual error -----
    # Kim 2019 Table 2 reports additive residual errors fitted separately
    # for the healthy-volunteer and the patient data.
    addSd <- addSdHealthy * DIS_HEALTHY + addSdPatient * (1 - DIS_HEALTHY)
    Cc ~ add(addSd)
  })
}
