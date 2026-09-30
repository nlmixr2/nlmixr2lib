Taneja_2019_glpg1690 <- function() {
  description <- "Two-compartment population PK with first-order absorption and dose-dependent apparent clearance, coupled to an effect-compartment Imax model for plasma lysophosphatidic acid (LPA) C18:2 reduction, for the autotaxin inhibitor GLPG1690 (ziritaxestat) in healthy volunteers and patients with idiopathic pulmonary fibrosis"
  reference <- "Taneja A, Desrivot J, Diderichsen PM, Blanque R, Allamasey L, Fagard L, Fieuw A, Van der Aar E, Namour F. Population Pharmacokinetic and Pharmacodynamic Analysis of GLPG1690, an Autotaxin Inhibitor, in Healthy Volunteers and Patients with Idiopathic Pulmonary Fibrosis. Clin Pharmacokinet. 2019;58(9):1175-1188. doi:10.1007/s40262-019-00755-3"
  vignette <- "Taneja_2019_glpg1690"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_GLPG1690_MGD = list(
      description = "Total daily GLPG1690 dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the linear decrease of apparent clearance with dose",
        "(Fig. 3: CL_i(DOSE) = TVCL_i * (1 - CLSLP * DOSE)).",
        "Total DAILY dose, not the per-administration amount: Table 4 reports",
        "identical steady-state AUC and identical 95% CI for the QD and BID",
        "rows of every total daily dose, which is only possible if both arms",
        "share one CL. Studied range 20-1500 mg/day. The linear form is",
        "extrapolation-unsafe: it reaches CL = 0 at 2410 mg/day, well above",
        "the 1500 mg/day maximum studied."
      ),
      source_name = "DOSE"
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral suspension), the comparator formulation in the first-in-human study",
      notes = paste(
        "1 = capsule, 0 = oral suspension. Absorption rate is 16.2% higher for",
        "the capsule (Sect. 3.3). Relative bioavailability is NOT a function of",
        "formulation in this model: F1 is anchored at 1 for both arms and only",
        "rifampicin coadministration moves it."
      ),
      source_name = "Formulation"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin (rifampin) coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant rifampicin)",
      notes = paste(
        "1 = GLPG1690 given after 10 days of rifampicin 600 mg once daily",
        "(study IND130687). Steady-state induction phase only, so the",
        "unqualified canonical is used rather than the",
        "CONMED_RIFAMPICIN_SD / _MD pair. Acts on absorption rate (+117%) and",
        "relative bioavailability (9.94% of the rifampicin-free value), not on",
        "clearance; the authors attribute this to intestinal P-gp induction",
        "(Sect. 4)."
      ),
      source_name = "DDI"
    ),
    DIS_IPF = list(
      description = "Idiopathic pulmonary fibrosis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = paste(
        "1 = patient with idiopathic pulmonary fibrosis (proof-of-concept study",
        "NCT02738801), 0 = healthy volunteer. Raises apparent central volume by",
        "559% and inflates both residual-error variances. The authors attribute",
        "the volume effect to the sparse PoC sampling schedule, which took no",
        "samples near Tmax and so under-estimated Cmax (Sect. 4) - it is a",
        "study-design artefact rather than a physiological difference."
      ),
      source_name = "Health status (HV vs IPF)"
    ),
    LPAC182_BL = list(
      description = "Baseline plasma lysophosphatidic acid C18:2",
      units = "peak-area ratio (LPA C18:2 / LPA C17:0 internal standard)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Unitless peak-area ratio, not a molar concentration: the LC-MS/MS",
        "method is non-quantitative because of endogenous LPA, so results are",
        "expressed relative to the LPA C17:0 internal standard (Sect. 2.2).",
        "Sets the PD model's baseline AND enters as a covariate on both IC50",
        "(power form) and Imax (linear on the probit scale), each centred on the",
        "typical value 0.36 used for the Sect. 2.3.9 typical-patient",
        "simulations. Observed range 0.104-1.305 across studies (Table 2)."
      ),
      source_name = "LPA C18:2 BL"
    ),
    OCC = list(
      description = "First-in-human single-ascending-dose occasion number",
      units = "(count)",
      type = "categorical",
      reference_category = "1 (the subject's first SAD occasion); 0 on records from every other study",
      notes = paste(
        "Volunteers in the SAD part received ascending doses on up to three",
        "occasions separated by a one-week washout and were treated as",
        "independent subjects in the analysis (Sect. 2.3.7). Imax is lower on",
        "the later occasions. Set OCC = 0 on all MAD, drug-drug-interaction and",
        "proof-of-concept records so both indicators zero out, following the",
        "Oosten 2016 fentanyl precedent."
      ),
      source_name = "Occasion"
    ),
    STUDY_FIH_MAD = list(
      description = "First-in-human multiple-ascending-dose study-part indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (every other study part)",
      notes = paste(
        "1 = record from the multiple-ascending-dose part of the first-in-human",
        "study NCT02179502 (150 mg BID, 600 mg QD or 1000 mg QD for 14 days).",
        "Used only in the product STUDY_FIH_MAD * DAY14, which gates the",
        "higher Imax the authors estimated at day 14 and later of that study",
        "part. The gate is required because the 12-week proof-of-concept study",
        "also samples past day 14 but does not carry this effect."
      ),
      source_name = "Study part (FiH MAD)"
    ),
    DAY14 = list(
      description = "Day-14-or-later landmark indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (before day 14 of dosing)",
      notes = paste(
        "1 = the observation falls on or after day 14 of dosing. Enters only as",
        "the product STUDY_FIH_MAD * DAY14 (see STUDY_FIH_MAD). Derive as",
        "as.integer(time_since_first_dose_days >= 14)."
      ),
      source_name = "Day 14 and later"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Sect. 2.3.3, demographics on PK and PK/PD parameters) but not retained in the final model; no estimate reported. Observed range 64.20-110 kg (Table 2)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened as a demographic covariate (Sect. 2.3.3) but not retained. Observed range 21-79 years (Table 2)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened as a demographic covariate (Sect. 2.3.3) but not retained. Observed range 20.0-39.1 kg/m^2 (Table 2)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened as a demographic covariate (Sect. 2.3.3) but not retained.",
        "Only the proof-of-concept study enrolled women (43.75% of its analysed",
        "subjects); the authors report higher median baseline LPA C18:2 in women",
        "(0.47 vs 0.36) and note that sex is confounded with the unequal",
        "representation, so no firm conclusion is drawn (Sect. 4, ESM Fig. S4 /",
        "Table S1). The sex signal therefore reaches the model only through",
        "LPAC182_BL."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "GLPG1690", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "GLPG1690", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "GLPG1690", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "GLPG1690", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 81,
    n_studies = 3,
    age_range = "21-79 years",
    weight_range = "64.20-110 kg",
    sex_female_pct = 8.6,
    disease_state = "healthy volunteers (2 studies) and patients with idiopathic pulmonary fibrosis (1 study)",
    dose_range = "20-1500 mg single oral dose (suspension or capsule); 150 mg BID and 600 or 1000 mg QD for 14 days; 600 mg QD for 12 weeks",
    regions = "Belgium (first-in-human), USA (drug-drug interaction), 14 sites (proof-of-concept)",
    notes = paste(
      "Pooled from NCT02179502 (first-in-human SAD 16 + MAD 24 male volunteers),",
      "IND130687 (rifampicin drug-drug interaction, 18 male volunteers) and",
      "NCT02738801 (proof-of-concept, 23 patients with IPF, of whom 17 received",
      "GLPG1690). Baseline demographics in Table 2. Subject counts differ from",
      "NONMEM ID counts because each SAD occasion was entered as an independent",
      "subject (48 SAD IDs from 16 volunteers). 1348 GLPG1690 concentrations and",
      "894 LPA C18:2 records were analysed; 11% (166/1514) of concentrations were",
      "below the 1.00 ng/mL limit of quantification and were excluded. Only the",
      "proof-of-concept study enrolled women (7 of the 16 analysed subjects)."
    )
  )

  ini({
    # -- Structural PK (Table 3, 'Estimate (%RSE)' column; footnote a = log-transformed) --
    lcl <- 3.15; label("Apparent clearance (CL/F, L/h)") # Table 3 CL: exp(3.15) = 23.3 L/h
    lvc <- 2.65; label("Apparent central volume of distribution (VP2/F, L)") # Table 3 VP2: exp(2.65) = 14.1 L
    lka <- -1.51; label("Absorption rate constant (KA, 1/h)") # Table 3 KA: exp(-1.51) = 0.222 /h
    lq <- 0.130; label("Apparent inter-compartmental clearance (Q23/F, L/h)") # Table 3 Q23: exp(0.130) = 1.14 L/h
    lvp <- 2.53; label("Apparent peripheral volume of distribution (VP3/F, L)") # Table 3 VP3: exp(2.53) = 12.5 L

    # -- Structural PD (Table 3) --
    probitimax <- 1.33; label("Maximal fractional reduction of plasma LPA C18:2 (Imax, probit scale)") # Table 3 Imax, footnote b probit-transformed: pnorm(1.33) = 0.908
    lic50 <- 4.74; label("Effective concentration producing half-maximal LPA C18:2 reduction (IC50, ng/mL)") # Table 3 IC50: exp(4.74) = 114 ng/mL
    lke0 <- -2.32; label("Effect-compartment equilibration rate constant (KEO, 1/h)") # Table 3 KEO: exp(-2.32) = 0.0982 /h
    lcprel <- -1.20; label("Weighting of effect-site relative to plasma concentration in the effective concentration (CPREL, unitless)") # Table 3 CPREL: exp(-1.20) = 0.301 (paper prints 0.303)

    # -- Covariate effects (Table 3) --
    e_dose_glpg1690_mgd_cl <- 0.000415; label("Fractional decrease in apparent clearance per mg of total daily dose (1/mg)") # Table 3 CL(DOSE): -0.0415% change/mg, i.e. -4.15% per 100 mg
    e_dis_ipf_vc <- 5.59; label("Fractional increase in apparent central volume in patients with IPF (unitless)") # Table 3 'VP2 in IPF patients': 559% change
    e_form_capsule_ka <- 0.162; label("Fractional increase in absorption rate for the capsule vs the oral suspension (unitless)") # Table 3 'KA (capsule)': 16.2% change
    e_conmed_rifampicin_ka <- 1.17; label("Fractional increase in absorption rate with concomitant rifampicin (unitless)") # Table 3 'KA with DDI': 117% change
    e_conmed_rifampicin_fdepot <- -2.31; label("Log relative bioavailability with concomitant rifampicin (unitless)") # Table 3 'Relative bioavailability F1 with DDI', footnote a: exp(-2.31) = 0.0993
    e_lpac182_bl_ic50 <- -0.473; label("Power exponent of baseline LPA C18:2 on IC50 (unitless)") # Table 3 'IC50 (LPA C18:2 BL)', footnote c: 1.1^-0.473 - 1 = -4.41% per +10%
    e_lpac182_bl_imax <- 0.5805; label("Linear effect of baseline LPA C18:2 on Imax, probit scale (per unit peak-area ratio)") # Table 3 'Imax (LPA C18:2 BL)', footnote d: pnorm(1.33 + 0.5805*0.1) = 91.7% at baseline +0.1
    e_occ2_imax <- -0.36; label("Shift in Imax on the second single-ascending-dose occasion, probit scale (unitless)") # Table 3 'Imax at occasion 2 of the FIH study, SAD part': pnorm(1.33 - 0.36) = 83.3%
    e_occ3_imax <- -0.4504; label("Shift in Imax on the third single-ascending-dose occasion, probit scale (unitless)") # Table 3 'Imax at occasion 3 of the FIH study, SAD part': pnorm(1.33 - 0.4504) = 81.0%
    e_day14_imax <- 0.2043; label("Shift in Imax from day 14 of the multiple-ascending-dose study part, probit scale (unitless)") # Table 3 'Imax at day 14 and later of the FIH study, MAD part': pnorm(1.33 + 0.2043) = 93.7%
    e_dis_ipf_ruv <- 0.209; label("Log variance ratio of the GLPG1690 residual error in patients with IPF (unitless)") # Table 3 'GLPG1690 RUV in patients with IPF': exp(0.209) = 1.233, i.e. 23.3% change
    e_dis_ipf_ruv_lpaC182 <- 0.500; label("Log variance ratio of the LPA C18:2 residual error in patients with IPF (unitless)") # Table 3 'Plasma LPA C18:2 RUV in patients with IPF': exp(0.500) = 1.649, i.e. 64.9% change

    # -- Between-subject variability (Table 3, 'BSV variance (%RSE)' column) --
    # Correlation 0.357 / sqrt(0.294 * 0.877) = 0.703, matching the printed 70.3%.
    etalcl + etalvc ~ c(0.294, 0.357, 0.877) # Table 3 CL, 'CL x VP2 covariance' and VP2 BSV variances
    etalka ~ 0.0194 # Table 3 KA BSV variance (sqrt = 13.9%cv)
    etalic50 ~ 0.129 # Table 3 IC50 BSV variance (sqrt = 35.9%cv)
    etaprobitimax ~ 0.0324 # Table 3 Imax BSV variance, on the probit scale

    # -- Residual unexplained variability --
    # Table 3 reports variances; the paper's %cv column is their square root.
    expSd <- 0.447214; label("Log-scale residual SD for GLPG1690 plasma concentration (unitless)") # Table 3 'GLPG1690 RUV': sqrt(0.200) = 0.447 (44.8%cv)
    expSd_lpaC182 <- 0.276043; label("Log-scale residual SD for plasma LPA C18:2 (unitless)") # Table 3 'Plasma LPA C18:2 RUV': sqrt(0.0762) = 0.276 (27.6%cv)
  })

  model({
    # 1. Individual PK parameters.
    #    Apparent clearance falls LINEARLY with total daily dose (Fig. 3:
    #    CL_i(DOSE) = TVCL_i * (1 - CLSLP * DOSE)). Verified against Table 4:
    #    at 600 mg/day, 23.3 * (1 - 0.000415 * 600) = 17.5 L/h and
    #    600 / 17.5 = 34.3 ug*h/mL, the printed steady-state AUC.
    cl <- exp(lcl + etalcl) * (1 - e_dose_glpg1690_mgd_cl * DOSE_GLPG1690_MGD)
    vc <- exp(lvc + etalvc) * (1 + e_dis_ipf_vc * DIS_IPF)
    ka <- exp(lka + etalka) *
      (1 + e_form_capsule_ka * FORM_CAPSULE) *
      (1 + e_conmed_rifampicin_ka * CONMED_RIFAMPICIN)
    q <- exp(lq)
    vp <- exp(lvp)
    ke0 <- exp(lke0)
    cprel <- exp(lcprel)

    # 2. Individual PD parameters. Both baseline-LPA effects are centred on the
    #    typical baseline 0.36 used for the Sect. 2.3.9 typical-patient
    #    simulations, which reproduces the paper's own statements: at
    #    LPAC182_BL = 0.36, Imax = pnorm(1.33) = 90.8% ("91%", Sect. 3.3), and a
    #    +10% baseline gives 1.1^-0.473 = -4.41% on IC50 (Table 3 footnote c).
    ic50 <- exp(lic50 + etalic50) * (LPAC182_BL / 0.36)^e_lpac182_bl_ic50
    imax <- pnorm(
      probitimax + etaprobitimax +
        e_lpac182_bl_imax * (LPAC182_BL - 0.36) +
        e_occ2_imax * (OCC == 2) +
        e_occ3_imax * (OCC == 3) +
        e_day14_imax * STUDY_FIH_MAD * DAY14
    )

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (Sect. 3.2 display equations, reproduced verbatim).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # GLPG1690 plasma concentration. Dose is in mg and volumes in L, so
    # central / vc is mg/L; the factor 1000 converts to the ng/mL the paper
    # reports throughout (LLOQ 1.00 ng/mL, IC50 114 ng/mL, Table 4 Cmax).
    Cc <- 1000 * central / vc

    # Hypothetical effect compartment; the state holds a CONCENTRATION (Cpe),
    # matching the paper's d Cpe / dt = KEO * (Cplasma - Cpe).
    d/dt(effect) <- ke0 * (Cc - effect)

    # 5. Bioavailability. F1 is anchored at 1 without rifampicin; rifampicin
    #    coadministration drops it to exp(-2.31) = 9.94% (Sect. 3.3).
    f(depot) <- exp(e_conmed_rifampicin_fdepot * CONMED_RIFAMPICIN)

    # 6. Observations.
    #    Effective concentration: CEFF = Cplasma + CPREL * Cpe (Sect. 3.2). This
    #    is the display equation, which is a WEIGHTED SUM rather than a convex
    #    combination; the Fig. 3 legend's looser wording ("weighting factor
    #    between plasma and effect compartment concentration") is not followed
    #    where it conflicts with the printed equation.
    ceff <- Cc + cprel * effect
    lpaC182 <- LPAC182_BL * (1 - imax * ceff / (ceff + ic50))

    # Both observations were fitted log-transformed with additive residual
    # error (Sect. 3.2), i.e. log-normal on the linear scale. A study effect
    # multiplies each residual VARIANCE by exp(theta) in patients with IPF, so
    # the SD carries exp(theta / 2).
    sdCc <- expSd * exp(0.5 * e_dis_ipf_ruv * DIS_IPF)
    sdLpaC182 <- expSd_lpaC182 * exp(0.5 * e_dis_ipf_ruv_lpaC182 * DIS_IPF)

    Cc ~ lnorm(sdCc)
    lpaC182 ~ lnorm(sdLpaC182)
  })
}
