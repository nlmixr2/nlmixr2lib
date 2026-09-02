Downes_2022_tobramycin <- function() {
  description <- "Two-compartment IV population PK model for tobramycin in hospitalized children less than 5 years of age with cystic fibrosis treated for a pulmonary exacerbation (Downes 2022). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 70 kg) while central and peripheral volumes scale linearly with body weight (exponent 1, reference 70 kg). Clearance additionally carries a power effect of age (exponent 0.136, reference 2.7 years), a power effect of bedside-Schwartz estimated GFR (exponent 0.246, reference 128 mL/min/1.73 m^2), and a 29.2% reduction with concomitant vancomycin. Between-subject variability is a full block on clearance and central volume; residual variability is proportional."
  reference <- "Downes KJ, Grim A, Shanley L, Rubenstein RC, Zuppa AF, Gastonguay MR. A Pharmacokinetic Analysis of Tobramycin in Patients Less than Five Years of Age with Cystic Fibrosis: Assessment of Target Attainment with Extended-Interval Dosing through Simulation. Antimicrob Agents Chemother. 2022;66(5):e02377-21. doi:10.1128/aac.02377-21"
  vignette <- "Downes_2022_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Tobramycin was given as a 30-minute IV infusion
  # (Materials and Methods, "Study design"), so the dose enters `central`
  # directly and there is no depot state. `central` is verified: Downes 2022
  # "Data collection" states that tobramycin was measured in SERUM in the CHOP
  # Chemistry Laboratory by competitive immunoassay (VITROS TOBRA reagent).
  # `peripheral1` is left unverified because the paper never states what matrix
  # the peripheral distribution compartment represents; "plasma" follows the
  # repository default for a mathematical peripheral compartment and is not a
  # paper-sourced claim.
  compartmentData <- list(
    central     = list(analyte = "tobramycin", units = "mg", specimen = "serum",  verified = TRUE),
    peripheral1 = list(analyte = "tobramycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Downes 2022 Results, 'Model development': clearance and Q were allometrically scaled for weight to 0.75 and normalized to 70 kg, while central (V1) and peripheral (V2) volumes were scaled linearly for weight, also normalized to 70 kg. Table 1 gives a first-course median of 11.2 kg (IQR 8.2-14.7) and an all-courses median of 12.4 kg (IQR 9.7-14.8); the 70 kg normalizer is the conventional adult reference weight and is far outside the observed range, so the reported thetas are extrapolated adult-equivalent values and only the weight-normalized quantities quoted in the text (0.252 L/hr/kg^0.75 for CL and 0.308 L/kg for V1) are within the fitted range.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at the start of the tobramycin course",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Included a priori on clearance based on prior published pediatric cystic-fibrosis tobramycin models (Downes 2022 Materials and Methods, 'Covariate selection'). Evaluated both as a Hill function and as an exponential covariate normalized to the population median; the retained form in Table 2 is the power (exponential-on-log) form (AGE/2.7)^0.136. The 2.7-year normalizer is the value printed in the Table 2 parameterization footnote and is close to but not identical with the Table 1 all-courses baseline median of 2.6 years (IQR 1.2-4); the difference is consistent with 2.7 being the median across analysis records rather than across first courses (first-course median 2.2 years, IQR 0.8-3.8). Eligibility restricted the cohort to less than 5 years of age, and Downes 2022 Discussion warns explicitly that the parameterization may differ beyond the observed age range.",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the bedside Schwartz equation, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Included a priori on clearance (Downes 2022 Materials and Methods, 'Covariate selection'), calculated by the bedside Schwartz equation, and allowed to vary over the first 48 h of therapy when creatinine was measured more than once, so the covariate is time-varying in the source analysis. The 128 mL/min/1.73 m^2 normalizer is the value printed in the Table 2 parameterization footnote and sits between the Table 1 first-course median of 126 (IQR 110-149) and all-courses median of 129 (IQR 110-148), again consistent with a median taken across analysis records. Renal impairment (eGFR < 60 mL/min/1.73 m^2) was an exclusion criterion and only 3 of 58 first courses (5.2%) had an eGFR below 90, so the cohort is essentially normal-to-supranormal in renal function and the model carries no information about renal impairment. Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates; the assay form here is the creatinine-based bedside Schwartz estimate.",
      source_name        = "GFR"
    ),
    CONMED_VANCOMYCIN = list(
      description        = "Concomitant intravenous vancomycin coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant vancomycin)",
      notes              = "Dichotomized (Y/N) as a TIME-VARYING covariate, allowed to vary over the course of tobramycin therapy (Downes 2022 Materials and Methods, 'Covariate selection'). Entered the model by forward selection from an exploratory screen of five nephrotoxic comedications and was the only one retained. Multiplicative factor 0.708^CONMED_VANCOMYCIN on clearance, i.e. a 29.2% (95% CI 3.7-55.7%) reduction in tobramycin CL when vancomycin is coadministered. Downes 2022 flags this as a novel and fragile finding: only 5 of 58 patients over 8 of 111 courses received concomitant vancomycin, the %RSE on the theta is 19.1, and the bootstrap 95% CI (0.462-1.24) crosses 1. The authors' own Monte Carlo simulations deliberately set all simulated patients to CONMED_VANCOMYCIN = 0 so as not to over-state the effect (Results, 'Target attainment'), and the Discussion states further studies are needed. Users reproducing the paper's target-attainment analysis should likewise hold this at 0.",
      source_name        = "VAN"
    )
  )

  # Screened in the Downes 2022 exploratory nephrotoxin analysis (Materials and
  # Methods, 'Covariate selection') but NOT retained in the final model, so
  # these are documentation only and are not referenced in model(). Each was
  # dichotomized (Y/N) as a time-varying covariate and offered to a forward
  # selection requiring a >= 6.63 drop in OFV; only concomitant vancomycin met
  # that threshold. Names follow the auto-approved CONMED_<agent> family; only
  # CONMED_VANCOMYCIN is registered in inst/references/covariate-columns.md as
  # an actually-used covariate, alongside the four screened-only siblings.
  covariatesDataExcluded <- list(
    CONMED_TICARCLAV = list(
      description        = "Concomitant ticarcillin-clavulanate coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant ticarcillin-clavulanate)",
      notes              = "Screened on tobramycin CL in the exploratory nephrotoxin forward selection and not retained. Listed in Table 1 footnote c as one of the concurrent nephrotoxic medications counted for the baseline characteristics.",
      source_name        = "ticarcillin/clavulanate"
    ),
    CONMED_TMPSMX = list(
      description        = "Concomitant trimethoprim-sulfamethoxazole coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant trimethoprim-sulfamethoxazole)",
      notes              = "Screened on tobramycin CL in the exploratory nephrotoxin forward selection and not retained. Listed in Table 1 footnote c as one of the concurrent nephrotoxic medications counted for the baseline characteristics.",
      source_name        = "trimethoprim-sulfamethoxazole"
    ),
    CONMED_PIPTAZO = list(
      description        = "Concomitant piperacillin-tazobactam coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant piperacillin-tazobactam)",
      notes              = "Screened on tobramycin CL in the exploratory nephrotoxin forward selection and not retained. Listed in Table 1 footnote c as one of the concurrent nephrotoxic medications counted for the baseline characteristics.",
      source_name        = "piperacillin-tazobactam"
    ),
    CONMED_NSAID = list(
      description        = "Concomitant non-steroidal anti-inflammatory drug coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant NSAID)",
      notes              = "Screened on tobramycin CL in the exploratory nephrotoxin forward selection and not retained. Listed in Table 1 footnote c as one of the concurrent nephrotoxic medications counted for the baseline characteristics. Downes 2022 additionally tested the nephrotoxin exposures as a class effect, grouping the agents as none versus >= 1 agent and as 0-1 versus >= 2 agents; neither composite was retained either.",
      source_name        = "NSAID"
    ),
    CONMED_ACYCLOVIR = list(
      description        = "Concomitant acyclovir coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant acyclovir)",
      notes              = "Counted toward the concurrent-nephrotoxic-medication tally reported in Table 1 (footnote c) but NOT screened individually on clearance: the Materials and Methods forward-selection list names only vancomycin, ticarcillin/clavulanate, trimethoprim-sulfamethoxazole, piperacillin-tazobactam and NSAIDs. Recorded here so the composition of the Table 1 nephrotoxin count is auditable.",
      source_name        = "acyclovir"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 58L,
    n_courses        = 111L,
    n_studies        = 1L,
    n_concentrations = 224L,
    age_range        = "less than 5 years by eligibility; Table 1 first-course median 2.2 years (IQR 0.8-3.8), all-courses median 2.6 years (IQR 1.2-4)",
    age_median       = "2.2 years at first course (IQR 0.8-3.8)",
    weight_range     = "first-course IQR 8.2-14.7 kg (full range not reported)",
    weight_median    = "11.2 kg at first course (IQR 8.2-14.7)",
    height_median    = "85.8 cm at first course (IQR 70.5-95.5)",
    sex_female_pct   = 43,
    disease_state    = "Cystic fibrosis with a pulmonary exacerbation requiring hospitalization and intravenous tobramycin. Excluded: estimated GFR < 60 mL/min/1.73 m^2 by bedside Schwartz, extracorporeal membrane oxygenation, postmenstrual age < 44 weeks, and concurrent nebulized tobramycin. Concurrent nephrotoxic medications at tobramycin initiation in 34.5% of first courses (1 agent 31.0%, 2 or more 3.4%).",
    dose_range       = "Clinician-chosen intravenous regimens; the CHOP formulary starting dose for this age group was 3.3 mg/kg every 8 h as a 30-minute infusion, and the observed first-course median was 3.2 mg/kg/dose (IQR 3.1-3.3). The model's simulation application in Downes 2022 projects 10-15 mg/kg once-daily extended-interval regimens, which the cohort did NOT receive.",
    regions          = "United States (Children's Hospital of Philadelphia, Philadelphia PA)",
    renal_function   = "Bedside-Schwartz eGFR first-course median 126 mL/min/1.73 m^2 (IQR 110-149); only 3 of 58 (5.2%) below 90 mL/min/1.73 m^2. Serum creatinine median 0.3 mg/dL (IQR 0.2-0.3).",
    notes            = "Retrospective single-center analysis of standard-of-care therapeutic-drug-monitoring data collected 2011-03-01 to 2018-09-01 (Downes 2022 Table 1; Materials and Methods). 61 patients / 115 courses were screened and 4 courses in 3 patients excluded (2 in premature infants, 2 with no concentrations), leaving 58 patients over 111 courses; 35 patients contributed one course, 13 two, and 10 three or more. Only concentrations drawn within the first 48 h of a course were used. Of 228 collected concentrations 4 were excluded as mistimed, leaving 224, of which 53 (23.7%) were below the 0.6 mg/L limit of quantification and were handled by the Beal M3 likelihood method. Because sampling was routine TDM peaks and troughs, Q and V2 were not identifiable from these data alone: Downes 2022 digitized the individual concentration-time profiles of 6 richly sampled adult cystic-fibrosis patients from Figure 2 of an earlier publication with WebPlotDigitizer, fit them to obtain prior distributions, and estimated the pediatric model by MAP-Bayesian penalized likelihood with informative priors on Q and V2. Fit in NONMEM 7.4 with the PDx-Pop 5.2.1 interface; covariate selection by backward elimination from a full model at a critical OFV change of 6.63, then an exploratory forward selection for nephrotoxic comedications at the same threshold. Inter-occasion variability on CL was evaluated and rejected (it raised the AIC, left the CL and V1 inter-individual variances and the point estimates unchanged, reduced residual variance only from 0.334 to 0.328, and inflated every %RSE). Bootstrap n = 1000 overall, plus n = 500 stratified by vancomycin receipt, by age band and by eGFR band, with all point estimates inside the bootstrap 95% CIs."
  )

  ini({
    # Structural parameters (Downes 2022 Table 2). The reference subject weighs
    # 70 kg (the conventional allometric reference, far above this cohort), is
    # 2.7 years old, has a bedside-Schwartz eGFR of 128 mL/min/1.73 m^2 and is
    # not receiving concomitant vancomycin. The text quotes the corresponding
    # weight-normalized values, which reproduce exactly from these thetas:
    #   CL 6.10 / 70^0.75 = 0.2521 L/hr/kg^0.75  (text: 0.252, 95% CI 0.233-0.271)
    #   V1 21.6 / 70      = 0.3086 L/kg          (text: 0.308, 95% CI 0.264-0.353)
    #   Q  4.73 / 70^0.75 = 0.1954 L/hr/kg^0.75  (text: 0.195, 95% CI 0.171-0.219)
    #   V2 6.69 / 70      = 0.0956 L/kg          (text: 0.096, 95% CI 0.081-0.110)
    lcl <- log(6.10); label("Clearance at WT=70 kg, AGE=2.7 y, CRCL=128 mL/min/1.73 m^2, no vancomycin (L/h)")  # Downes 2022 Table 2: theta CL 6.10 L/hr/70 kg (%RSE 3.89)
    lvc <- log(21.6); label("Central volume at WT=70 kg (L)")                                                   # Downes 2022 Table 2: theta V1 21.6 L/70 kg (%RSE 7.41)
    lq  <- log(4.73); label("Intercompartmental clearance at WT=70 kg (L/h)")                                    # Downes 2022 Table 2: theta Q 4.73 L/hr/70 kg (%RSE 6.24)
    lvp <- log(6.69); label("Peripheral volume at WT=70 kg (L)")                                                 # Downes 2022 Table 2: theta V2 6.69 L/70 kg (%RSE 7.55)

    # Allometric exponents. Downes 2022 Results ('Model development') states
    # that CL and Q "were allometrically scaled for weight to 0.75, normalized
    # to 70-kg, while central (V1) and peripheral (V2) volume parameters were
    # scaled linearly for weight, normalized to 70-kg", citing the standard
    # allometry reference. The Table 2 parameterization footnote prints
    # [WT/70]^0.75 on CL and Q and a bare (WT/70) on V1 and V2, with no
    # exponent row, no %RSE and no bootstrap CI, so both are structural
    # assumptions rather than estimates.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT/70) for CL and Q (unitless)")   # Downes 2022 Table 2 footnote: [WT/70]^0.75 on TVCL and TVQ
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT/70) for Vc and Vp (unitless)")  # Downes 2022 Table 2 footnote: (WT/70) on TVV1 and TVV2

    # Covariate effects on clearance, all estimated (each has a %RSE and a
    # bootstrap 95% CI in Table 2), so none is wrapped in fixed().
    e_age_cl  <- 0.136; label("Power exponent on (AGE/2.7) for CL (unitless)")   # Downes 2022 Table 2: theta AGE 0.136 (%RSE 12.1)
    e_crcl_cl <- 0.246; label("Power exponent on (CRCL/128) for CL (unitless)")  # Downes 2022 Table 2: theta GFR 0.246 (%RSE 32.4)

    # Exponentiated-switch form, matching the Table 2 footnote's (theta VAN ^
    # VAN): the factor multiplies CL when the indicator is 1 and is inert when
    # it is 0. 1 - 0.708 = 0.292 reproduces the 29.2% CL reduction quoted in
    # the Results text.
    e_conmed_vancomycin_cl <- 0.708; label("Multiplicative factor on CL when CONMED_VANCOMYCIN = 1")  # Downes 2022 Table 2: theta VAN 0.708 (%RSE 19.1)

    # Between-subject variability. Downes 2022 Materials and Methods
    # ('Base model') specifies exponential BSV, and Table 2 footnote c states
    # that a FULL BLOCK covariance matrix was used for the CL and V1 random
    # effects. Table 2 reports the omega elements directly on the variance
    # scale, with the %CV in parentheses computed as sqrt(omega^2) rather than
    # the log-normal sqrt(exp(omega^2) - 1): sqrt(0.0296) = 0.172 and
    # sqrt(0.0633) = 0.252 match the printed 17.2% and 25.2% exactly, so the
    # tabulated numbers are used as-is with no CV-to-variance conversion. The
    # printed off-diagonal correlation 0.751 is reproduced by
    # 0.0325 / sqrt(0.0296 * 0.0633) = 0.7508, confirming that 0.0325 is the
    # covariance and not a variance. The block is positive definite
    # (determinant 0.0296 * 0.0633 - 0.0325^2 = 8.17e-4 > 0).
    # Downes 2022 Table 2: omega 1,1-CL = 0.0296 (%RSE 17.8);
    # omega 1,2-CL:V1 = 0.0325 with r = 0.751 (%RSE 32.3);
    # omega 2,2-V1 = 0.0633 (%RSE 28.9).
    etalcl + etalvc ~ c(0.0296,
                        0.0325, 0.0633)

    # No etalq / etalvp: Downes 2022 Table 2 footnote c states "Interindividual
    # random effects fixed to 0 for Q and V2", and Results ('Model
    # development') explains why -- "Given poor estimation of between-subject
    # random effects for Q and V2 during base model development, they were
    # fixed to 0 for these parameters." A variance of exactly 0 is equivalent
    # to omitting the eta, and omitting it keeps the OMEGA block non-singular.

    # Residual variability. Materials and Methods ('Base model'): "Residual
    # variability (RV) was estimated using a proportional error model." Table 2
    # reports it on the variance scale as 0.334 with 57.8% CV in parentheses,
    # and sqrt(0.334) = 0.5779 recovers that 57.8%, so the tabulated value is a
    # variance and nlmixr2's propSd is its square root.
    propSd <- 0.578; label("Proportional residual error (fraction)")  # Downes 2022 Table 2: residual variability 0.334 (57.8% CV, %RSE 6.89); 0.578 = sqrt(0.334)
  })
  model({
    # Individual PK parameters, transcribed from the Downes 2022 Table 2
    # parameterization footnote:
    #   TVCL = thetaCL * ([WT/70]^0.75) * (AGE/2.7)^thetaAGE
    #                  * (GFR/128)^thetaGFR * (thetaVAN^VAN)
    #   TVV1 = thetaV1 * (WT/70)
    #   TVQ  = thetaQ  * ([WT/70]^0.75)
    #   TVV2 = thetaV2 * (WT/70)
    # The footnote as typeset renders the AGE and GFR terms as
    # "(AGE/2.7 ^ u AGE)" and "(GFR/128 ^ u GFR)", i.e. with the closing
    # parenthesis displaced past the exponent, exactly as the WT term is
    # correctly typeset "([WT/70] ^ 0.75)". The intended grouping is
    # (AGE/2.7)^thetaAGE and (GFR/128)^thetaGFR, which is also the form
    # Materials and Methods ('Covariate selection') describes -- "an
    # exponential covariate normalized to the population median" -- and is the
    # only reading under which the four weight-normalized typical values
    # quoted in the Results text reproduce from the Table 2 thetas.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * (AGE / 2.7)^e_age_cl *
      (CRCL / 128)^e_crcl_cl * e_conmed_vancomycin_cl^CONMED_VANCOMYCIN
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L, so central/vc is mg/L, which is the mg/L
    # unit Downes 2022 reports serum tobramycin in and is numerically
    # identical to ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
