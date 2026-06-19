Ngamprasertwong_2016_propofol_sheep <- function() {
  description <- "Preclinical (sheep). Maternal-fetal population PK model of propofol in mid-gestational pregnant Dorset ewes (Ngamprasertwong 2016; N = 8 ewe-fetus pairs at 110-125 days gestation; term ~147-150 days). Two-compartment maternal disposition (central + peripheral1) linked to a single fetal compartment via a reversible inter-compartmental clearance QM-F; fetal clearance was tested but estimated near zero (<0.001 L/min, RSE >100%) and set to zero in the final model. Maternal clearance scales with heart rate via the normalised power model CL = theta1 * (HR/158)^theta2; no other covariate (gestational age, body weight, blood pressure, uterine blood flow) reached statistical significance. Inter-individual variability was estimated on CL and QM-F; IIV on Vc, Q, Vp, and VFetus was fixed to zero in the source and is omitted here. Residual error is purely proportional, with separate variances for maternal-ewe and fetal observations."
  reference <- paste(
    "Ngamprasertwong P, Dong M, Niu J, Venkatasubramanian R, Vinks AA,",
    "Sadhasivam S.",
    "Propofol pharmacokinetics and estimation of fetal propofol exposure",
    "during mid-gestational fetal surgery: a maternal-fetal sheep model.",
    "PLoS ONE 2016;11(1):e0146563.",
    "doi:10.1371/journal.pone.0146563.",
    sep = " "
  )
  vignette <- "Ngamprasertwong_2016_propofol_sheep"
  paper_specific_compartments <- c("fetus")

  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    HR = list(
      description        = "Maternal heart rate during the propofol observation window",
      units              = "beats/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Encoded as the per-ewe median heart rate over the propofol-infusion observation window per Ngamprasertwong 2016 Methods (no separate intra-individual HR trajectory was modelled). Power covariate effect normalised to a reference of 158 beats/min as written in the Table 2 equation CL = theta1 * (HR/158)^theta2; the paper's Results text quotes 4.17 L/min as the typical CL 'in a typical ewe with a median heart rate of 135 beats/min', so the population reference HR (158) in the equation differs from the typical-subject median HR (135) reported in the narrative -- the table equation is authoritative for the model.",
      source_name        = "HR"
    )
  )

  population <- list(
    species        = "sheep (pregnant Dorset ewe)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "gestational age 110-125 days (term ~147-150 days, mid-gestation)",
    age_median     = "gestational age 116.5 days",
    weight_range   = "60-82 kg",
    weight_median  = "72.5 kg",
    sex_female_pct = 100,
    disease_state  = "Singleton pregnant Dorset ewes under general anesthesia with propofol + remifentanil + desflurane, instrumented for chronic maternal-fetal physiologic monitoring (femoral arterial / venous catheters, umbilical and bilateral uterine flow probes). Used as a mid-gestational fetal-surgery anesthesia model.",
    dose_range     = "Induction: propofol 3 mg/kg IV bolus + succinylcholine 1.5 mg/kg via the maternal femoral vein. Anesthesia phase 1 (0-60 min): propofol 450 ug/kg/min IV continuous infusion + remifentanil 0.5 ug/kg/min. Anesthesia phase 2 (60-150 min): propofol reduced to 75 ug/kg/min + remifentanil 0.25 ug/kg/min + 1.5 MAC desflurane (10.2 percent end-tidal). Propofol infusion stopped at 150 min; final blood sample at 180 min.",
    regions        = "United States (Cincinnati Children's Hospital Medical Center).",
    n_observations = "160 propofol plasma measurements (80 paired ewe-fetus simultaneous draws) at 5, 15, 25, 60, 75, 100, 110, 150, and 180 min after the start of the propofol infusion. Sampling times were derived from a D-optimal design (Ngamprasertwong 2016 Methods / Pharmacokinetics study).",
    notes          = "Propofol was assayed by LC-MS/MS (iC42 Integrated Solutions in Clinical Research and Development, University of Colorado, Denver). LLOQ = 0.05 ng/mL; intra- and inter-assay variability < 10 percent. Anesthesia / instrumentation protocol: protocol 0D03027, Cincinnati Children's Hospital Committee for Animal Care. The data underlying the modelling are deposited under doi:10.5281/zenodo.35672 (open access). The 'typical ewe' in the Results narrative has a median heart rate of 135 beats/min and body weight 71.6 kg (mean across the cohort); the Table 2 equation uses 158 beats/min as the population-level reference for the (HR/158) normalisation."
  )

  ini({
    # Final-model fixed-effect parameter estimates (Ngamprasertwong 2016
    # Table 2). NONMEM v7.2.0 FOCE-I; structural model is a 2-compartment
    # maternal disposition (central + peripheral1) plus a single fetal
    # compartment exchanging with the maternal central via QM-F; fetal
    # clearance was tested but driven to ~zero with RSE > 100 percent
    # and was dropped from the final model.

    # Maternal clearance (CL) is parameterised via the normalised power
    # model CL = theta1 * (HR/158)^theta2 (Table 2). theta1 is the
    # reference-subject CL at HR = 158 beats/min.
    lcl     <- log(4.17);    label("Maternal clearance at HR=158 bpm (CL, L/min)")           # Table 2: theta1 = 4.17 (RSE 8.4 percent)
    lvc     <- log(37.7);    label("Maternal central volume of distribution (Vc, L)")         # Table 2: Vc = 37.7 (RSE 8.0 percent)
    lq      <- log(1.22);    label("Maternal inter-compartmental clearance (Q, L/min)")       # Table 2: Q = 1.22 (RSE 21.6 percent)
    lvp     <- log(60.8);    label("Maternal peripheral volume of distribution (Vp, L)")      # Table 2: Vp = 60.8 (RSE 21.9 percent)
    lqmf    <- log(0.0138);  label("Maternal-fetal inter-compartmental clearance (QM-F, L/min)") # Table 2: QM-F = 0.0138 (RSE 26.5 percent); units row in the source PDF lists '(min)' which is a typo for 'L/min' per the Fig 3 caption
    lvfetus <- log(0.144);   label("Fetal volume of distribution (VFetus, L)")                # Table 2: VFetus = 0.144 (RSE 6.5 percent)

    # Covariate effect (Table 2 equation footer): CL = theta1 * (HR/158)^theta2.
    # Reference HR = 158 beats/min is the population-level normaliser written
    # into the Table 2 equation; the per-subject median HR used as the typical
    # value in the Results narrative is 135 beats/min (see covariateData[[HR]]$notes).
    e_hr_cl <- 0.764;        label("Power exponent of (HR/158) on maternal CL (theta2, unitless)") # Table 2: theta2 = 0.764 (RSE 28.3 percent)

    # IIV (log-normal). NONMEM reports the %CV directly under
    # 'Inter-individual variability (% CV)' in Table 2; the internal
    # omega^2 is computed via omega^2 = log(CV^2 + 1).
    # IIVs on Vc (omega2^2), Q (omega3^2), Vp (omega4^2), and VFetus
    # (omega6^2) were '0 FIX' (held at zero during estimation), so the
    # corresponding etas are omitted entirely (Przybylowski 2015
    # propofol convention).
    etalcl  ~ 0.046429    # Table 2: omega1^2 (CL) = 21.8 %CV; omega^2 = log(0.218^2 + 1) = 0.046429
    etalqmf ~ 0.366170    # Table 2: omega5^2 (QM-F) = 66.5 %CV; omega^2 = log(0.665^2 + 1) = 0.366170

    # Residual error. Ngamprasertwong 2016 Methods Eq (2) describes a
    # combined proportional + additive error model on linear-scale
    # plasma concentration, but Table 2 reports only the proportional
    # components for ewe and fetus (additive components were dropped
    # in the final model). NONMEM reports the proportional component as
    # %CV, which is numerically equal to the proportional-residual SD.
    propSd        <- 0.260;  label("Proportional residual SD for maternal ewe Cc (fraction)")    # Table 2: sigma^2 prop ewe = 26.0 %CV (RSE 26.9 percent)
    propSd_Cfetus <- 0.218;  label("Proportional residual SD for fetal Cfetus (fraction)")       # Table 2: sigma^2 prop fetus = 21.8 %CV (RSE 32.1 percent)
  })
  model({
    # Individual parameters. Maternal clearance follows the Table 2
    # normalised power covariate model; the remaining structural
    # parameters carry no covariate effect (their %CV IIVs were fixed
    # to zero in the source).
    cl      <- exp(lcl + etalcl) * (HR / 158)^e_hr_cl
    vc      <- exp(lvc)
    q       <- exp(lq)
    vp      <- exp(lvp)
    qmf     <- exp(lqmf + etalqmf)
    vfetus  <- exp(lvfetus)

    # Micro-rate constants for the explicit ODE form. The maternal-fetal
    # transfer rate is reversible: drug crosses from maternal central to
    # fetus at rate (qmf / vc) * central and returns at rate (qmf / vfetus)
    # * fetus. No fetal clearance term is included because the source
    # paper found it negligible (<0.001 L/min, RSE >100 percent) and
    # dropped it from the final model.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k1f <- qmf / vc
    kf1 <- qmf / vfetus

    # IV bolus + continuous IV infusion are administered into `central`
    # via the events table. The fetus compartment has no dose input.
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k1f * central + kf1 * fetus
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(fetus)       <-  k1f * central - kf1 * fetus

    # Observations. Maternal plasma concentration (Cc) is the canonical
    # central-compartment output; fetal plasma concentration (Cfetus)
    # is a paper-named non-PK output keyed to the `fetus` physiologic
    # compartment.
    Cc     <- central / vc
    Cfetus <- fetus   / vfetus

    Cc     ~ prop(propSd)
    Cfetus ~ prop(propSd_Cfetus)
  })
}
