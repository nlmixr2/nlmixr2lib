Hanzel_2022_vedolizumab <- function() {
  description <- paste(
    "Two-compartment population PK model for vedolizumab (humanised",
    "anti-alpha4-beta7 integrin IgG1 monoclonal antibody) with parallel linear",
    "and Michaelis-Menten elimination in adults with active Crohn's disease",
    "(LOVE-CD trial), with interindividual and inter-occasion variability",
    "(induction / maintenance) on linear clearance and time-varying effects of",
    "serum albumin, antibodies to vedolizumab and anti-TNF-naive status on",
    "linear clearance, plus the sequential first-order discrete-time Markov",
    "model for endoscopic remission (SES-CD < 4) and dropout at weeks 26 and",
    "52 driven by the individual predicted week-22 trough concentration",
    "(Hanzel 2022). Q, Vp, Km and Vmax were held at the Rosario 2015 values."
  )
  reference <- paste(
    "Hanzel J, Dreesen E, Vermeire S, Lowenberg M, Hoentjen F, Bossuyt P,",
    "Clasquin E, Baert FJ, D'Haens GR, Mathot R. Pharmacokinetic-Pharmacodynamic",
    "Model of Vedolizumab for Targeting Endoscopic Remission in Patients With",
    "Crohn Disease: Posthoc Analysis of the LOVE-CD Study. Inflamm Bowel Dis.",
    "2022;28:689-699. doi:10.1093/ibd/izab143 (PMC9071095). A corrigendum",
    "(doi:10.1093/ibd/izab270; PMC9071094) corrects reference 36 only and does",
    "not change any parameter value."
  )
  vignette <- "Hanzel_2022_vedolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vedolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vedolizumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Linear effect on linear clearance centred on the cohort median 41 g/L: CL multiplier = 1 + (ALB - 41) * e_alb_cl (Hanzel 2022 Equation 1 and Appendix $PK 'TVCL = THETA(3) * (1 + (ALBUMIN - 41)*THETA(9)) ...'). Lowering albumin from 41 to 28 g/L raises CL by 26%.",
      source_name = "ALBUMIN"
    ),
    ADA_POS = list(
      description = "Antibodies to vedolizumab (ATV) detected at the sampling time",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ATV-negative)",
      notes = "Time-varying ('ATV i,j equals 1 if ATV were positive at time j'); measured with a drug-sensitive assay before each infusion. Multiplicative effect on linear clearance e_ada_cl^ADA_POS (+89% when positive). 4 of 108 patients (10 of 737 samples) were ATV-positive.",
      source_name = "ADA"
    ),
    PRIOR_TNF = list(
      description = "Previous exposure to anti-TNF agents",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (previously exposed to anti-TNF agents; the paper's typical patient)",
      notes = "Time-constant. The source column is TNFNAIVE (1 = never exposed to anti-TNF agents), so TNFNAIVE = 1 - PRIOR_TNF and the paper's multiplier enters as e_tnfnaive_cl^(1 - PRIOR_TNF): anti-TNF-naive patients have 0.755-fold the linear clearance of previously exposed patients. The typical value lcl therefore refers to a PRIOR_TNF = 1 patient, as in the paper. 96 of 108 patients (89%) were previously exposed.",
      source_name = "TNFNAIVE"
    ),
    OCC = list(
      description = "Dosing occasion for the inter-occasion variability on linear clearance: 1 = induction, 2 = maintenance",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- decomposed into occasion indicators that select the per-occasion eta on CL",
      notes = "Two occasions per Hanzel 2022 Methods: induction (up to the week-14 dose, i.e. time < 98 days) and maintenance (from the week-14 dose onward, time >= 98 days). Any value other than 1 or 2 switches the IOV term off.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 108L,
    n_studies = 1L,
    age_median = "36 years (IQR 28-46)",
    weight_median = "71 kg (IQR 61-82)",
    sex_female_pct = 69,
    disease_state = "Active Crohn's disease (CDAI > 220 with mucosal ulceration at baseline ileocolonoscopy); baseline CDAI median 261 (IQR 238-312), SES-CD median 12 (IQR 7-17), disease duration median 9 years.",
    dose_range = "Vedolizumab 300 mg IV at weeks 0, 2 and 6 and every 8 weeks thereafter through week 52; 68 patients (63%) received an additional infusion at week 10 for an insufficient CDAI response.",
    regions = "Netherlands, Belgium",
    prior_tnf_pct = 89,
    albumin_median = "41 g/L (IQR 38-43)",
    ada_positive_pct = 3.7,
    notes = "LOVE-CD (NCT02646683), a phase 4 prospective open-label multicentre trial. 108 of 110 enrolled patients had at least one quantifiable vedolizumab sample (737 trough samples, all drawn before an infusion). Concomitant corticosteroids 41%, immunomodulators 20% at baseline. Endoscopic remission (SES-CD < 4, central reading) in 36 (33%) at week 26 and 40 (37%) at week 52. Demographics from Hanzel 2022 Table 1."
  )

  ini({
    # PK structural parameters -- Hanzel 2022 Table 2, final model. Typical
    # patient: albumin 41 g/L, ATV-negative, previously exposed to anti-TNF.
    lcl <- log(0.215); label("Linear clearance CL_L for the typical patient (L/day)") # Table 2 final: CL_L = 0.215 L/d (RSE 3%)
    lvc <- log(4.92); label("Central volume of distribution V1 (L)") # Table 2 final: V1 = 4.92 L (RSE 4%)
    lvp <- fixed(log(1.65)); label("Peripheral volume of distribution V2 (L)") # Table 2: V2 = 1.65 FIX (from Rosario 2015)
    lq <- fixed(log(0.12)); label("Intercompartmental clearance Q (L/day)") # Table 2: Q = 0.12 FIX (from Rosario 2015)
    lkm <- fixed(log(0.964)); label("Michaelis-Menten constant Km (mg/L)") # Table 2: Km = 0.964 FIX (from Rosario 2015)
    lvmax <- fixed(log(0.265)); label("Maximum Michaelis-Menten elimination rate Vm (mg/day)") # Table 2: Vm = 0.265 FIX (from Rosario 2015)

    # Covariate effects on linear clearance -- Table 2 final and Equation 1.
    e_alb_cl <- -0.020; label("Linear albumin effect on CL_L per g/L from 41 g/L (1/(g/L))") # Table 2: Albumin on CL_L = -0.020 (RSE 29%)
    e_ada_cl <- 1.89; label("Multiplier on CL_L when antibodies to vedolizumab are present (unitless)") # Table 2: Antidrug antibodies on CL_L = 1.89 (RSE 26%)
    e_tnfnaive_cl <- 0.755; label("Multiplier on CL_L for anti-TNF-naive patients (unitless)") # Table 2: No previous biologic exposure on CL = 0.755 (RSE 9%)

    # IIV and IOV on CL_L. Table 2 footnote: 'The coefficient of variation was
    # calculated as the square root of variance', so omega^2 = CV^2.
    etalcl ~ 0.068644 # Table 2 final: IIV CL_L 26.2 CV% -> 0.262^2
    # The Appendix control stream also carries an eta on V1 ('IIV-V1'), but
    # neither Table 2 nor the Results text reports it (IIV and IOV 'were
    # estimated for linear clearance'); carried as a zero-variance slot.
    etalvc ~ fixed(0) # Appendix $PK V1 = TVV1 * EXP(ETA(2)); variance not reported in Table 2
    # IOV: two occasions (induction / maintenance), one shared variance
    # ($OMEGA BLOCK(1) ... $OMEGA BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.023104 # Table 2 final: IOV CL_L 15.2 CV% -> 0.152^2
    etaiov_cl_2 ~ fixed(0.023104) # SAME as occasion 1 (Appendix $OMEGA BLOCK(1) SAME)

    # Residual error -- Table 2 final; Appendix $ERROR
    # W = SQRT(THETA(1)**2 + IPRED*IPRED*THETA(2)**2) with $SIGMA 1 held at 1.
    addSd <- 0.469; label("Additive residual error (mg/L)") # Table 2 final: additive error = 0.469 mg/L (RSE 63%)
    propSd <- 0.189; label("Proportional residual error (fraction)") # Table 2 final: proportional error = 0.189 (RSE 8%)

    # Markov endoscopic-remission model -- Table 2 final (Markov block) and
    # Appendix $PRED. States: 0 = no endoscopic remission, 1 = endoscopic
    # remission (SES-CD < 4), 2 = dropout (absorbing). Driver: individual
    # predicted vedolizumab concentration at week 22 (IPRED22).
    emax_01 <- fixed(0.7); label("Maximum probability of the 0 -> 1 (remission) transition (fraction)") # Table 2: Emax01 = 70 FIX; Appendix THETA(1) 0.7 FIX
    lec50_01 <- log(20.0); label("Week-22 concentration giving half of emax_01 for the 0 -> 1 transition (mg/L)") # Table 2 final: EC50,01 = 20.0 mg/L (RSE 23%)
    emax_02 <- fixed(1); label("Maximum probability of the 0 -> 2 (dropout) transition (fraction)") # Table 2: Emax02 = 1 FIX; Appendix THETA(3) 1 FIX
    let50_02 <- log(515); label("Assessment day giving half of emax_02 for the 0 -> 2 transition (day)") # Table 2 final: ET50,02 = 515 days (RSE 21%)
    emax_10 <- fixed(1); label("Maximum probability of the 1 -> 0 (loss of remission) transition (fraction)") # Table 2 final: Emax10 = 1 FIX
    lec50_10 <- log(1.78); label("Week-22 concentration halving the 1 -> 0 transition probability (mg/L)") # Table 2 final: EC50,10 = 1.78 mg/L (RSE 55%)
    emax_12 <- fixed(1); label("Maximum probability of the 1 -> 2 (dropout from remission) transition (fraction)") # Table 2 final: Emax12 = 1 FIX
    lec50_12 <- log(0.47); label("Week-22 concentration halving the 1 -> 2 transition probability (mg/L)") # Table 2 final: EC50,12 = 0.47 mg/L (RSE 94%)
  })

  model({
    # Occasion indicators for IOV on CL_L (induction = 1, maintenance = 2).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # Equation 1 / Appendix $PK: linear albumin effect centred on 41 g/L,
    # power-form multipliers for anti-TNF-naive status and ATV positivity.
    cl <- exp(lcl + etalcl + iov_cl) *
      (1 + (ALB - 41) * e_alb_cl) *
      e_tnfnaive_cl^(1 - PRIOR_TNF) *
      e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    km <- exp(lkm)
    vmax <- exp(lvmax)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Appendix $DES: two-compartment model, IV infusion into central, parallel
    # linear and Michaelis-Menten elimination from central.
    Cc <- central / vc
    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      vmax * Cc / (km + Cc)
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Markov model (Appendix $PRED; final-model parameterisation of Table 2).
    # The exposure driver is IPRED22, the individual predicted trough
    # concentration at week 22 (day 154, before the week-22 infusion). The
    # probabilities below use the CURRENT Cc, so they reproduce the paper only
    # when read at t = 154 days before any week-22 dose.
    ec50_01 <- exp(lec50_01)
    et50_02 <- exp(let50_02)
    ec50_10 <- exp(lec50_10)
    ec50_12 <- exp(lec50_12)
    # 0 -> 1: Emax in the week-22 concentration.
    p01 <- emax_01 * Cc / (ec50_01 + Cc)
    # 0 -> 2: Emax in the assessment day, DAYS = VISIT * 7 + 0.01.
    p02_wk26 <- emax_02 * (26 * 7 + 0.01) / (et50_02 + 26 * 7 + 0.01)
    p02_wk52 <- emax_02 * (52 * 7 + 0.01) / (et50_02 + 52 * 7 + 0.01)
    # 1 -> 0 and 1 -> 2: inhibitory Emax in the week-22 concentration with
    # Emax = 1 (see the vignette for why the inhibitory form is used).
    p10 <- emax_10 * (1 - Cc / (ec50_10 + Cc))
    p12 <- emax_12 * (1 - Cc / (ec50_12 + Cc))
    # Appendix $PRED transition fractions: remission is resolved first and
    # dropout competes only among those not transitioning, e.g.
    # P(0 -> 2) = P02 * (1 - P01) and P(0 -> 0) = 1 - P01 - P02 * (1 - P01).
    # Marginal state probabilities, all patients in state 0 at baseline:
    prob_endorem_wk26 <- p01
    prob_dropout_wk26 <- p02_wk26 * (1 - p01)
    prob_noendorem_wk26 <- 1 - prob_endorem_wk26 - prob_dropout_wk26
    prob_endorem_wk52 <- prob_endorem_wk26 * (1 - p10) * (1 - p12) +
      prob_noendorem_wk26 * p01
    prob_dropout_wk52 <- prob_dropout_wk26 +
      prob_endorem_wk26 * p12 * (1 - p10) +
      prob_noendorem_wk26 * p02_wk52 * (1 - p01)
    prob_noendorem_wk52 <- 1 - prob_endorem_wk52 - prob_dropout_wk52

    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
