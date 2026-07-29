Zhao_2010_mycophenolic_acid <- function() {
  description <- "Two-compartment population PK model for mycophenolic acid (MPA, the active moiety delivered as oral mycophenolate mofetil MMF) in children with idiopathic nephrotic syndrome (Zhao 2010). First-order absorption (ka = 5.16 1/h) with absorption lag time (tlag = 0.215 h) into a central compartment. Apparent oral clearance CL/F (typical value 9.7 L/h at the cohort medians WT = 23.5 kg, ALB = 38.6 g/L) is modeled with two covariates: a power effect of body weight on CL/F with exponent 0.753 referenced to 23.5 kg (close to allometric but estimated, not fixed), and an unusual linear-in-ratio effect of serum albumin in the form CL/F = q1 * (WT/23.5)^q2 * [1 - q3 * (ALB/38.6)] with q1 = 22.5 L/h, q2 = 0.753, q3 = 0.570 (higher serum albumin reduces apparent CL/F, consistent with stronger MPA-albumin binding in nephrotic patients with restored albumin). Apparent central V1/F = 22.3 L; apparent peripheral V2/F was fixed at 250 L (estimation between 100 and 600 L was non-identifiable; the fixed value lies in the range reported for adult transplant cohorts). Apparent inter-compartment clearance Q/F = 18.8 L/h. Exponential inter-individual variability is estimated on lag time, V1/F, Q/F, and CL/F (no IIV on ka or V2/F). A proportional residual error (44.6%) on MPA plasma concentration completes the model. Dosing in this packaged form is in mg of MMF; the MMF-to-MPA hydrolysis is implicit in the apparent bioavailability F."
  reference <- paste(
    "Zhao W, Elie V, Baudouin V, Bensman A, Andre JL, Brochard K,",
    "Broux F, Cailliez M, Loirat C, Jacqz-Aigrain E.",
    "Population pharmacokinetics and Bayesian estimator of mycophenolic",
    "acid in children with idiopathic nephrotic syndrome.",
    "Br J Clin Pharmacol. 2010;69(4):358-366.",
    "doi:10.1111/j.1365-2125.2010.03615.x.",
    sep = " "
  )
  vignette <- "Zhao_2010_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the pharmacokinetic occasion (time-varying across the M1 and M6 PK days).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort medians 23.6 kg at M1 and 23.3 kg at M6 (Table 1); the model uses the pooled-cohort median 23.5 kg as the reference weight per Results 'The median body weight ... in the population were 23.5 kg'. Enters CL/F as a power term with exponent q2 = 0.753.",
      source_name        = "weight"
    ),
    ALB = list(
      description        = "Serum albumin concentration measured on the pharmacokinetic day.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort means 36.0 g/L (M1) and 39.6 g/L (M6), Table 1; the model uses the pooled-cohort median 38.6 g/L as the reference per Results 'The median ... albumin in the population were ... 38.6 g l-1'. Enters CL/F in the unusual linear-in-ratio form [1 - q3 * (ALB/38.6)] with q3 = 0.570 -- this is NOT the typical power form (ALB/ref)^exponent. At ALB = 38.6 the factor evaluates to (1 - 0.570) = 0.430 so that the typical-CL anchor q1 = 22.5 L/h yields CL = 9.675 L/h at the cohort medians, matching the abstract value 9.7 L/h. Higher albumin reduces CL/F (biologically sensible: more bound MPA).",
      source_name        = "albumin"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 23L,
    n_profiles     = 41L,
    n_observations = 285L,
    n_studies      = 1L,
    age_range      = "2.9-14.9 years (mean 7.4 +/- 3.9 years)",
    age_median     = "5.4 years (M6 occasion); 5.2 years (M1 occasion)",
    weight_range   = "14.0-83.2 kg (mean 29.9 +/- 18.0 kg)",
    weight_median  = "23.5 kg (pooled cohort median; M1 23.6 kg, M6 23.3 kg)",
    sex_female_pct = 21.7,
    race_ethnicity = "Not reported in source paper.",
    disease_state  = "Pediatric patients with steroid-dependent idiopathic nephrotic syndrome (INS) who had relapsed despite cyclophosphamide therapy. All patients were on background prednisone at the PK day (mean 47 mg per 2 days at M1, 13 mg per 2 days at M6).",
    dose_range     = "MMF 1200 mg/m^2/day given orally BID; per-dose amounts 300-1000 mg MMF (mean 563 mg M1, 540 mg M6).",
    regions        = "Six pediatric nephrology centres in France (Paris Robert-Debre, Paris Trousseau, Nancy, Toulouse, Rouen, Marseille).",
    sampling_window= "Two PK occasions per patient (M1 ~ day 30, M6 ~ day 180); samples at pre-dose, 0.5, 1, 2, 4, 8, 12 h post-dose. 41 profiles from 23 patients (18 had both M1 and M6; 3 had M1 only; 2 had M6 only).",
    albumin_range  = "26.5-45.5 g/L (M1); 25.6-44.0 g/L (M6)",
    crcl_range     = "All patients had Schwartz-formula creatinine clearance > 25 mL/min (range 87.6-249.5 mL/min); renal function did not influence MPA CL/F in this cohort.",
    notes          = "Covariates screened but not retained in the final model: age, height, body surface area, sex, creatinine clearance, urine protein, haemoglobin, cholesterol, alkaline phosphatase, ASAT, ALAT, prednisone dose, time after start of therapy. BSA, age, haemoglobin, and alkaline phosphatase showed initial correlation with CL/F in graphical screening (Maitre method) but did not survive forward + backward covariate selection (forward P < 0.05, backward P < 0.01); only WT and ALB were retained. Body weight on CL/F was kept in preference to body surface area (delta-OFV -23.6 forward, +14.5 backward vs -19.6 forward for BSA which dropped out in the backward step). Patient baseline characteristics from Table 1; final-model parameters from Table 3; covariate selection from Table 2."
  )

  ini({
    # Final-model parameter estimates from Table 3 ('Population
    # pharmacokinetic parameters of MMF') of Zhao 2010. Standard errors
    # (% SE) reported in Table 3 are not encoded as fixed() flags because
    # they are precision summaries on estimated point estimates, not
    # fixed-during-estimation flags. The only fixed parameter in the
    # final model is V2/F per the Results text 'V2/F was fixed to 250 l'.

    # Absorption parameters. Lag time is estimated and carries IIV;
    # ka is estimated but with no IIV per the Results text 'Interindividual
    # variability ... was then estimated for V1/F, Q/F, CL/F and lag time'
    # (Ka excluded explicitly).
    ltlag   <- log(0.215);        label("Absorption lag time (h)")                                                    # Table 3 lag time = 0.215 h
    lka     <- log(5.16);         label("First-order absorption rate constant ka (1/h)")                              # Table 3 Ka = 5.16 1/h

    # Disposition: 2-compartment apparent oral PK. V2/F was fixed to 250 L
    # per Results 'After testing different values for V2/F, between 100 and
    # 600. V2/F was fixed to 250 l, as this gave the best estimation for
    # both pharmacokinetic parameters and residual variability and was in
    # reasonable agreement the reported values in transplant patients
    # ranging from 137 l to 496 l'.
    lvc     <- log(22.3);         label("Apparent central volume V1/F (L)")                                           # Table 3 V1/F = 22.3 L
    lvp     <- fixed(log(250));   label("Apparent peripheral volume V2/F (L), fixed at 250 L")                        # Table 3 V2/F = 250 L (fixed)
    lq      <- log(18.8);         label("Apparent inter-compartment clearance Q/F (L/h)")                             # Table 3 Q/F = 18.8 L/h

    # Apparent oral clearance. The paper writes
    #   CL/F = q1 * (WT/23.5)^q2 * [1 - q3 * (ALB/38.6)]
    # with q1 = 22.5 L/h, q2 = 0.753, q3 = 0.570 (Table 3). lcl is the
    # log of q1 (the typical-CL anchor before the ALB factor multiplies in
    # 0.430 at reference ALB to give 22.5 * 0.430 = 9.675 L/h, matching the
    # abstract value 9.7 L/h).
    lcl     <- log(22.5);         label("Apparent oral clearance anchor q1 (L/h); CL/F = q1 * (WT/23.5)^q2 * [1-q3*(ALB/38.6)]")  # Table 3 q1 = 22.5 L/h

    # Covariate effects. e_wt_cl is the power exponent (q2) of the WT/23.5
    # ratio on CL/F. e_alb_cl is the linear coefficient (q3) of the
    # ALB/38.6 ratio inside the [1 - q3 * (ALB/38.6)] factor on CL/F --
    # note this is NOT a power form; see vignette Assumptions and
    # deviations for the rationale (the paper writes it this way in
    # Equation just above Table 2).
    e_wt_cl   <- 0.753;           label("Power exponent of WT/23.5 on CL/F (unitless)")                               # Table 3 q2 = 0.753
    e_alb_cl  <- 0.570;           label("Linear coefficient of ALB/38.6 in [1 - q3 * (ALB/38.6)] on CL/F (unitless)") # Table 3 q3 = 0.570

    # Inter-individual variability. Source paper reports IIV as % CV
    # (Table 3 'Interindividual variability'); translated here to the
    # log-normal log-scale variance via omega^2 = log(1 + CV^2). Per
    # Results 'Interindividual variability ... was then estimated for V1/F,
    # Q/F, CL/F and lag time' -- so ka and V2/F have NO IIV. Intraindividual
    # variability (between-occasion) was deliberately not modeled per
    # Results 'Intraindividual variability was not included in the model
    # because of imprecision in the estimates of variance'.
    etaltlag ~ 0.257                                                                                                  # Table 3 IIV lag time = 54.0% CV -> log(1 + 0.540^2)
    etalvc   ~ 0.472                                                                                                  # Table 3 IIV V1/F     = 79.9% CV -> log(1 + 0.799^2)
    etalq    ~ 0.285                                                                                                  # Table 3 IIV Q/F      = 57.6% CV -> log(1 + 0.576^2)
    etalcl   ~ 0.0474                                                                                                 # Table 3 IIV CL/F     = 22.0% CV -> log(1 + 0.220^2)

    # Residual variability. Source paper text 'A proportional model was
    # used for the residual random effects'; Table 3 reports 44.6% as the
    # 'Residual proportional' magnitude. Encoded as a proportional SD
    # (fractional) on the linear-space concentration.
    propSd   <- 0.446;            label("Proportional residual SD (fraction of Cc)")                                  # Table 3 Residual proportional = 44.6%
  })

  model({
    # Individual parameters. Log-normal IIV on ltlag, lvc, lq, lcl only;
    # ka and V2/F (= exp(lvp)) carry no etas because the source paper did
    # not estimate IIV on those parameters.
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka)
    vc   <- exp(lvc   + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq    + etalq)

    # CL/F: power-of-WT * [1 - linear-of-ALB] structure per Table 3.
    # The minus sign in [1 - q3 * (ALB/38.6)] is part of the model
    # equation in the paper (not a sign-of-coefficient subtlety): the
    # paper estimated q3 as a positive 0.570 inside this functional form,
    # so the term decreases CL/F with rising ALB and yields 0.430 at the
    # reference ALB = 38.6 g/L.
    cl <- exp(lcl + etalcl) * (WT / 23.5)^e_wt_cl * (1 - e_alb_cl * (ALB / 38.6))

    # Micro-constants and ODEs for the 2-compartment depot model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Absorption lag on the depot compartment.
    alag(depot) <- tlag

    # Observation: MPA plasma concentration with proportional residual.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
