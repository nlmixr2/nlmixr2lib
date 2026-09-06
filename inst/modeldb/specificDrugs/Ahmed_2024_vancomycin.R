Ahmed_2024_vancomycin <- function() {
  description <- "One-compartment intravenous population PK model for vancomycin in Sudanese adult inpatients, developed from routine therapeutic-drug-monitoring peak and trough concentrations at a single hospital in Khartoum. Clearance is a median-centered power function of creatinine clearance; volume of distribution carries no covariates. IMPORTANT: the CRCL column is NOT a conventional creatinine clearance in mL/min -- the source computes it as Bjornson's creatinine production rate divided by serum creatinine with the body-weight and 14.4 unit factors omitted, which deflates it about five-fold relative to a Cockcroft-Gault value; see covariateData$CRCL."
  reference   <- "Ahmed KA, Ibrahim A, Gonzalez D, Nur AO. Population Pharmacokinetics and Model-Based Dose Optimization of Vancomycin in Sudanese Adult Patients with Renal Impairment. Drug Des Devel Ther. 2024;18:81-95. doi:10.2147/DDDT.S432439"
  vignette    <- "Ahmed_2024_vancomycin"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the source: the single disposition
  # compartment holds vancomycin and is sampled as plasma concentration by the
  # EMIT immunoassay, calibration range 2-50 mg/L (Methods, "Vancomycin
  # Assay"). Both peak (1 h post-dose) and trough (30 min pre-dose) samples
  # inform the fit (Methods, "Dosing and Sample Collection").
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance as computed by Ahmed 2024 -- Bjornson creatinine production rate divided by serum creatinine, WITHOUT the body-weight multiplier or the 14.4 unit-conversion divisor of the published Bjornson method",
      units              = "mL/min as labelled by the source; dimensionally (mg/kg/24 h)/(mg/dL) as actually computed",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "READ THIS BEFORE SUPPLYING A CRCL VALUE TO THIS MODEL. Ahmed 2024 Methods, 'Patients and Data Collection',",
        "Equations 1-3: CLcr = Rcr / SCR, with Rcr(males) = 27 - 0.173 * age and Rcr(females) = 25 - 0.175 * age, both in",
        "mg/kg/24 h, and SCR in mg/dL. The published Bjornson method (Bjornsson TD, Clin Pharmacokinet 1979) completes the",
        "calculation as CrCl = Rcr * body weight / (SCr * 14.4) to reach mL/min; Ahmed 2024 states explicitly that 'the data",
        "lacked patient weight and height' and therefore stops at Rcr / SCr, while still labelling the column mL/min. The",
        "omitted factor is weight / 14.4, i.e. about 4.9 for a 70 kg adult, so this column runs roughly five-fold BELOW a",
        "Cockcroft-Gault or complete-Bjornson creatinine clearance for the same patient. The arithmetic is confirmed by the",
        "paper's own Table 1: the median male aged 65 with SCr 1.2 mg/dL gives (27 - 0.173 * 65) / 1.2 = 13.1, against the",
        "reported cohort median of 12.7. Because e_crcl_cl scales CL as a power of this column, feeding it a conventional",
        "mL/min creatinine clearance would inflate typical CL by about 4.9^0.49 = 2.2-fold. Users must reproduce the source's",
        "own calculation, not substitute a standard renal-function estimate. This also explains why a general medical-ward",
        "cohort appears almost uniformly renally impaired (median 12.7, IQR 5.52-25.78, range 1.4-107.5): the column is",
        "systematically deflated, not the cohort uniformly anuric. Enters CL as the median-centered power term",
        "(CRCL / 12.7)^0.49 -- see the model() block for why 12.7 is used and what the source does and does not print.",
        sep = " "
      ),
      source_name        = "CLcr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at the time of the therapeutic-drug-monitoring episode",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and significant in univariate forward inclusion (Table 2 model 4, dOFV -16.9, p < 0.05) but not",
        "retained: once CLcr was in the model, adding age gave dOFV -0.06 (Table 2 model 8, p > 0.05). Age nonetheless enters",
        "the model indirectly, because it is an input to the Bjornson Rcr equations that generate CRCL.",
        sep = " "
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened on CL as log-transformed SCr and significant (Table 2 model 2, dOFV -43.73, p < 0.05), but CLcr gave the",
        "larger drop in both OFV (-49.69) and BSV on CL (omega 0.46 versus 0.49), so the CLcr model was selected as final",
        "(Results, 'Model Development'). The two are not independent -- SCr is the denominator of the CRCL calculation.",
        "Cohort median 1.2 mg/dL (IQR 0.7-2.5, range 0.2-9.5; Table 1).",
        sep = " "
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Screened on both V (Table 2 model 6, dOFV -0.86) and CL (model 7, dOFV -0.30) and rejected at p > 0.05 in each case;",
        "also rejected when added to the CLcr model (model 10, dOFV -1.07). Cohort median 2.5 g/dL (IQR 2.1-2.9,",
        "range 1.4-4.1; Table 1). Reported here in the source's g/dL, not the canonical register unit g/L.",
        sep = " "
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and significant in univariate forward inclusion (Table 2 model 5, dOFV -28.84, p < 0.05) but not",
        "retained once CLcr was in the model (Table 2 model 9, dOFV -0.34, p > 0.05). Cohort median 52 mg/dL",
        "(IQR 26.1-89, range 5-215; Table 1, where it is headed 'Blood urea').",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = NULL,
      type        = "categorical",
      notes       = paste(
        "Screened on CL. The Discussion reports that adding sex to CL reduced the OFV by 6.03 points (p > 0.05) at the forward",
        "step and that it was removed at backward elimination, so it does not appear in Table 2 at all and no coefficient is",
        "published. Cohort 66 male / 33 female (Table 1). As with age, sex enters indirectly: it selects between the two",
        "Bjornson Rcr equations that generate CRCL.",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 99L,
    n_observations = 194L,
    n_studies      = 1L,
    age_range      = "18-90 years; median 65 (IQR 50-75) (Table 1)",
    age_median     = "65 years",
    weight_range   = "NOT RECORDED. Ahmed 2024 states that patient weight and height were absent from the medical records, which is why the Cockcroft-Gault equation could not be used and why the CRCL column omits its weight factor.",
    sex_female_pct = 33,
    race_ethnicity = "Not reported beyond nationality; single-country cohort of Sudanese adults. The paper motivates the study by noting that no vancomycin population PK information existed for Sudanese patients.",
    disease_state  = "Adult inpatients receiving intravenous vancomycin at a single hospital. Patients under 18 years, pregnant patients, and patients on renal replacement therapy were excluded, so the model carries no information about dialysis. As reported by the source the cohort is dominated by renal impairment (median CLcr 12.7), but see covariateData$CRCL: the renal-function column is computed without its weight and unit factors and is therefore about five-fold deflated relative to a conventional creatinine clearance, so the apparent severity of impairment is partly an artefact of the covariate definition.",
    renal_function = "CLcr median 12.7 (IQR 5.52-25.78, range 1.4-107.5) in the source's own units (Table 1); serum creatinine median 1.2 mg/dL (IQR 0.7-2.5, range 0.2-9.5)",
    dose_range     = "500-1000 mg intravenously every 12 h, infused over 60 min (Methods, 'Dosing and Sample Collection'); administered one to two times per day (Results, 'Patient Characteristics')",
    regions        = "Sudan (Aliaa Specialist Hospital, Khartoum)",
    notes          = "Retrospective single-centre observational cohort of patients treated between August 2016 and January 2019. 194 concentrations from 99 patients: 129 troughs (median 16.22 mg/L, IQR 11.1-26.53) drawn 30 min before the next dose and 65 peaks (median 29.55 mg/L, IQR 22.38-37.36) drawn 1 h after the dose, all from the fourth dose onwards. No samples were below the limit of quantitation. Assay: enzyme-multiplied immunoassay (EMIT), calibration range 2-50 mg/L. Fitted in MonolixSuite 2020R1 by SAEM. Baseline demographics are Table 1 of Ahmed 2024."
  )

  ini({
    # All estimates are the "VALUE" (Monolix) column of Ahmed 2024 Table 3.
    # The model was fitted in MonolixSuite 2020R1 with log-normally distributed
    # individual parameters and an exponential BSV model (Methods, "Population
    # Pharmacokinetic Analysis"), so the population values below are entered on
    # the log scale.

    lvc <- log(65);   label("Volume of distribution (L)")  # Ahmed 2024 Table 3, V_pop = 65 L (SE 6.12, RSE 9.41%; bootstrap median 66.78, 95% CI 62.49-71.18)
    lcl <- log(2.02); label("Clearance at the reference CRCL of 12.7 (L/h)")  # Ahmed 2024 Table 3, CL_pop = 2.02 L/h (SE 0.13, RSE 6.40%; bootstrap median 2.004, 95% CI 1.92-2.08)

    # Covariate effect. Monolix names it beta_CL_logtCLcr: the slope of
    # log(CL) on the log-transformed, centered covariate logtCLcr, which is
    # identical to a power exponent on the centered ratio. See the model()
    # block for the two conflicting forms the paper prints and why the
    # centered one is used.
    e_crcl_cl <- 0.49; label("Power exponent on (CRCL / 12.7) for CL (unitless)")  # Ahmed 2024 Table 3, beta_CL_logtCLcr = 0.49 (SE 0.064, RSE 13%; bootstrap median 0.488, 95% CI 0.443-0.541)

    # Inter-individual variability. Table 3 heads this block "Standard
    # Deviation of the Random Effects" and Monolix reports omega as an SD, so
    # the tabulated 0.39 / 0.46 are SDs of the log-normal etas and the
    # variances entered here are their squares. The 0.46 is corroborated in
    # the Results text, which tracks omega for CL falling "from 0.66 to 0.46"
    # when the CLcr effect was added.
    etalvc ~ 0.1521  # Ahmed 2024 Table 3, omega_V  = 0.39 (SD) -> variance 0.39^2 = 0.1521 (RSE 19.3%; bootstrap median 0.42, 95% CI 0.35-0.49)
    etalcl ~ 0.2116  # Ahmed 2024 Table 3, omega_Cl = 0.46 (SD) -> variance 0.46^2 = 0.2116 (RSE 12.5%; bootstrap median 0.45, 95% CI 0.41-0.49)

    # Residual error. Both proportional and combined models were tested
    # (Methods); Table 3 reports only the proportional term b under the
    # heading "Standard Deviation of the Proportional Error", so the final
    # model is proportional-only. Monolix's proportional error model is
    # y = f + b * f * e with e ~ N(0, 1), so b maps directly onto propSd.
    propSd <- 0.28; label("Proportional residual error (fraction)")  # Ahmed 2024 Table 3, b = 0.28 (SE 0.023, RSE 8.26%; bootstrap median 0.28, 95% CI 0.26-0.29)
  })

  model({
    # Covariate model on clearance.
    #
    #   CL_i = 2.02 * (CRCL / 12.7)^0.49 * exp(eta_CL)
    #   V_i  = 65 * exp(eta_V)
    #
    # THE CENTERING IS DELIBERATE AND THE PAPER PRINTS TWO CONFLICTING FORMS.
    # Ahmed 2024 Methods, "Population Pharmacokinetic Analysis" (p. 83), gives
    # the form that was actually fitted:
    #
    #   log(CL) = log(CL_pop) + beta_CL_logtCLcr * logtCLcr + eta_CL
    #   where   logtCLcr = log(CLcr / mean(CLcr))
    #   i.e.    CL = CL_pop * (CLcr / mean(CLcr))^beta_CL_logtCLcr * exp(eta_CL)
    #
    # Equation 5 on p. 85 and the Table 3 footnote both drop the denominator
    # and print CL = CL_pop * CLcr^0.49 * exp(eta). That uncentered form is
    # arithmetically impossible: it returns 2.02 * 12.7^0.49 = 6.9 L/h at the
    # cohort median and 2.02 * 54.5^0.49 = 14.2 L/h at the top simulated CLcr
    # group, against the 2.22 and 4.28 L/h that the paper's own Table 4 reports
    # as the median simulated clearance for those groups. The centered form is
    # therefore used, per the Methods equation.
    #
    # THE CENTERING CONSTANT ITSELF IS NOT PRINTED. The Methods define the
    # reference as mean(CLcr), a value the paper never reports; the only
    # central-tendency statistic published for the covariate is the median,
    # 12.7 (Table 1), which is the value used here. Back-solving the paper's
    # Table 4 median simulated clearances against CL_pop = 2.02 and the 0.49
    # exponent implies an effective reference of about 11.8 -- consistent to
    # within 2% across all five independent CLcr groups, and close to what a
    # Monolix log-transform centered on the mean of log(CLcr) would give (the
    # geometric mean, which for a right-skewed covariate sits just below the
    # median). Using 12.7 rather than that back-solved 11.8 leaves typical CL
    # about 4% low against Table 4 across the whole covariate range; that is
    # accepted here in preference to fitting a constant to the validation
    # target. The vignette quantifies the residual offset.
    cl <- exp(lcl + etalcl) * (CRCL / 12.7)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    # One-compartment disposition with first-order elimination. Ahmed 2024
    # Results, "Model Development": "A one-compartment model with first-order
    # elimination best characterized vancomycin's PK"; one- and two-compartment
    # models were both evaluated (Methods).
    kel <- cl / vc

    # Doses enter the central compartment directly as intravenous infusions.
    # The source infused each dose over 60 min (Methods, "Dosing and Sample
    # Collection"); the infusion duration is a property of the event table
    # (rate / dur) rather than of the model.
    d/dt(central) <- -kel * central

    # Vancomycin plasma concentration in mg/L (doses in mg, volume in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
