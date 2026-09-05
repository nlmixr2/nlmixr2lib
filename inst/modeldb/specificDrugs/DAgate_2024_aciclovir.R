DAgate_2024_aciclovir <- function() {
  description <- "One-compartment IV population PK model for aciclovir in term and pre-term neonates with suspected systemic (herpes simplex virus) infection. Total clearance is the sum of an allometrically-scaled, post-menstrual-age-maturing residual clearance (tubular secretion plus metabolism) and the individual creatinine clearance (glomerular filtration) entering at unit slope; central volume scales linearly with body weight and is 2.67-fold higher in the presence of systemic infection."
  reference   <- "D'Agate S, Ruiz Gabarre D, Della Pasqua O. Population pharmacokinetics and dose rationale for aciclovir in term and pre-term neonates with herpes. Pharmacol Res Perspect. 2024;12(3):e1193. doi:10.1002/prp2.1193"
  vignette    <- "DAgate_2024_aciclovir"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against D'Agate 2024 Section 4.1 ("one-
  # compartment open model with zero-order infusion and first-order
  # elimination") and Figure 1 (aciclovir PLASMA concentration versus time).
  compartmentData <- list(
    central = list(analyte = "aciclovir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 1.37 kg, the population median stated beneath D'Agate 2024 Methods Eq. 2 and Eq. 3 ('1.37 kg is the median weight for the population'); Table 1 reports a median weight of 1420 g (range 373-5720 g) over the 32 enrolled infants. Enters CL as (WT/1.37)^0.75 and V as (WT/1.37)^1, both exponents pre-defined and fixed (Methods Section 2.3 item 1 for CL, and 'The effect of body weight on V was evaluated considering a pre-defined fixed allometric exponent with value 1' for V). The source records a single weight per infant alongside birth weight, so treat it as baseline unless a time-varying record is available.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Post-menstrual age (gestational age plus post-natal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months: D'Agate 2024 states post-menstrual age in weeks throughout, and the maturation half-time PMA50 = 50 weeks is only meaningful on that scale (see the PAGE entry in inst/references/covariate-columns.md, which provides for this). Table 1 median 31 weeks (range 25-41). Drives the maturation fraction PMA^HILL / (PMA50^HILL + PMA^HILL) on the residual clearance; Methods Section 2.2 states the sigmoidal function is used 'based on the assumption that PMA can be considered a proxy for maturation-related changes in renal clearance'.",
      source_name        = "PMA"
    ),
    CRCL = list(
      description        = "Individual creatinine clearance (Schwartz formula; raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. D'Agate 2024 Methods Section 2.3 item 3 computes CLCR with the Schwartz formula CLCR = k * HT / SCr, with k = 0.33 when gestational age is below 36 weeks and k = 0.45 otherwise, HT the subject height in cm and SCr the serum creatinine in mg/L; 'SCr was measured multiple times for each subject; so, CLCR was introduced as a time-varying covariate.' The Schwartz value is BSA-normalized (mL/min/1.73 m^2), and the paper then states 'The values obtained were converted from mL/min/1.73 m^2 to L/h using the individual body surface area of the subject as calculated with Gehan and George formula' -- i.e. the model consumes the INDIVIDUAL, non-BSA-normalized creatinine clearance. This column therefore carries the individual (de-normalized) value in raw mL/min, matching the raw-mL/min variant of the canonical CRCL column documented in inst/references/covariate-columns.md (Delattre 2010, Georges 2009 precedent); model() applies the 0.06 mL/min -> L/h conversion. To construct it from a BSA-normalized Schwartz estimate use CRCL = CLCR_schwartz * BSA / 1.73. It enters total clearance ADDITIVELY at unit slope (Table 2: 'CL (L/h) = theta1*(WT/1.37)^0.75*(PMA/(theta3 + PMA)) + CLCR'), representing the glomerular-filtration arm; the Discussion reports it to be about 25% of total aciclovir clearance.",
      source_name        = "CLCR"
    ),
    DIS_INFECT_ACTIVE = list(
      description        = "Systemic (herpes simplex virus) infection indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no evidence of systemic infection)",
      notes              = "D'Agate 2024 Methods Section 2.3 used STUDY as the discrete proxy for disease status because individual baseline viral load was unavailable: 'The use of study as a proxy for the disease status provided a suitable alternative to the lack of individual details on baseline viral load' (Results Section 4.1). Study 2 (multicentre, infants 23-34 weeks gestational age with suspected systemic HSV infection) had the larger proportion of positive virological findings and the lower observed concentrations, hence the larger volume of distribution; Study 1 is the reference. The paper's Table 2 abbreviation list names the column DIS (disease status). Set to 1 for the infected group, which multiplies V by 2.67. Note that the paper's own virtual-cohort simulations (Table 3) correspond to DIS = 0 -- see the model vignette.",
      source_name        = "DIS"
    )
  )

  # Screened during covariate model building but NOT retained in the final
  # model (D'Agate 2024 Methods Section 2.2 lists the candidate set; the
  # Discussion states 'Post-natal age, race and co-medication use were also
  # tested as potentially influential covariates but were not found to be
  # statistically significant'). No point estimates are reported for these,
  # so they are documentation only and are not referenced in model().
  covariatesDataExcluded <- list(
    PNA = list(
      description = "Post-natal age",
      units       = "days",
      type        = "continuous",
      notes       = "Table 1 median 3 days (range 1-30). Screened; not statistically significant (Discussion)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1 reports 17 female / 15 male among the 32 enrolled infants. Screened; not retained."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1 reports White/Black/Asian = 20/11/1. Race screened as a candidate covariate; not statistically significant (Discussion)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1 reports White/Black/Asian = 20/11/1. Race screened as a candidate covariate; not statistically significant (Discussion)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Table 1 median 0.9 mg/dL (range 0.3-1.8). Screened as a covariate in its own right; it enters the final model only indirectly, as an input to the Schwartz CRCL column."
    ),
    WT_BIRTH = list(
      description = "Birth weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Table 1 median 1295 g (range 420-4840 g). Screened alongside current weight; current weight was retained."
    ),
    CONMED_VASOPRESSIN = list(
      description = "Concomitant vasopressin use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1 reports 1 subject. Co-medication use screened; not statistically significant (Discussion)."
    ),
    CONMED_DOPAMINE = list(
      description = "Concomitant dopamine use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1 reports 4 subjects. Co-medication use screened; not statistically significant (Discussion)."
    ),
    CONMED_EPINEPHRINE = list(
      description = "Concomitant epinephrine use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1 reports 7 subjects. Co-medication use screened; not statistically significant (Discussion)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28L,
    n_studies      = 2L,
    age_range      = "Gestational age 23-40 weeks (median 30); post-menstrual age 25-41 weeks (median 31); post-natal age 1-30 days (median 3)",
    ga_range       = "23-40 weeks (median 30)",
    pma_range      = "25-41 weeks (median 31)",
    weight_range   = "0.373-5.720 kg (median 1.420 kg)",
    weight_median  = "1.42 kg",
    sex_female_pct = 53.1,
    race_ethnicity = c(White = 62.5, Black = 34.4, Asian = 3.1),
    disease_state  = "Term and pre-term neonates and young infants with suspected systemic infection; Study 2 specifically enrolled infants with suspected systemic herpes simplex virus infection.",
    renal_function = "Serum creatinine median 0.9 mg/dL (range 0.3-1.8). Post-hoc Schwartz creatinine clearance median 15.8 mL/min/1.73 m^2 (90% CI 9.4-40.0), against a post-hoc total clearance of 71.0 mL/min/1.73 m^2 (90% CI 26.2-257.6), so glomerular filtration accounts for roughly one quarter of total aciclovir clearance (Results Section 4.1).",
    co_medication  = "Vasopressin 1 subject, dopamine 4 subjects, epinephrine 7 subjects (Table 1, of the 32 enrolled).",
    dose_range     = "Intravenous aciclovir; the currently recommended neonatal regimen is 20 mg/kg every 8 h (60 mg/kg/day).",
    regions        = "United States (data from Sampson et al.; Study 1 single-centre, Study 2 multicentre)",
    notes          = "Demographics are D'Agate 2024 Table 1, reported for the 32 enrolled infants. Ninety-two plasma samples were collected; 9 were excluded before model development as contaminated or drawn during infusion, leaving 83 samples from 28 infants in the final data set (Results Section 4). The population metadata n_subjects therefore records the 28 infants contributing to the fit, while Table 1 percentages are over the 32 enrolled."
  )

  ini({
    # Structural parameters. Reference subject: WT = 1.37 kg (population
    # median weight) with the maturation and creatinine-clearance terms
    # evaluated at the individual's own PMA and CLCR.
    lcl <- log(0.748); label("Residual (non-glomerular-filtration) clearance at 1.37 kg and full maturation (L/h)")  # D'Agate 2024 Table 2 theta1 = 0.748 L/h (RSE 18%; bootstrap median 0.730, 95% CI 0.429-1.04); Methods Eq. 3 names it theta_CL,res
    lvc <- log(1.93);  label("Central volume of distribution at 1.37 kg without systemic infection (L)")              # D'Agate 2024 Table 2 theta2 = 1.93 L (RSE 39%; bootstrap median 1.89, 95% CI 0.823-3.50)

    # Allometric exponents. Both are pre-defined and held constant rather
    # than estimated: Methods Section 2.3 item 1 ("the use of pre-defined
    # fixed allometric exponents (0.75 for CL)") and, for V, "The effect of
    # body weight on V was evaluated considering a pre-defined fixed
    # allometric exponent with value 1." Neither carries an RSE in Table 2.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on the residual clearance (unitless)")  # D'Agate 2024 Methods Section 2.3 item 1; the exponent 0.75 also appears explicitly in the Table 2 CL equation
    e_wt_vc <- fixed(1);    label("Allometric exponent on the central volume (unitless)")      # D'Agate 2024 Methods Section 2.3, paragraph on V; the Table 2 V equation uses (WT/1.37) with no exponent, i.e. exponent 1

    # Maturation of the residual clearance, PMA^HILL / (PMA50^HILL + PMA^HILL).
    # Table 2 prints the final-model equation with this term collapsed to
    # (PMA / (theta3 + PMA)), i.e. HILL = 1; the paper reports no estimate,
    # standard error or bootstrap interval for HILL anywhere.
    pma_tm50 <- fixed(50); label("Post-menstrual age at 50% of mature residual clearance (weeks)")  # D'Agate 2024 Table 2 theta3 = 50, footnote b "Fixed to literature value"
    pma_hill <- fixed(1);  label("Hill coefficient of the post-menstrual-age maturation function (unitless)")  # D'Agate 2024 Table 2 final-model CL equation is printed as (PMA/(theta3 + PMA)), i.e. the Methods Eq. 3 sigmoid with HILL = 1; no HILL estimate is reported

    # Creatinine clearance enters total clearance additively at UNIT slope:
    # Table 2 writes '... + CLCR' with no coefficient, the paper having
    # converted CLCR to L/h beforehand (Methods Section 2.3 item 3). This
    # coefficient is therefore purely the mL/min -> L/h unit conversion
    # (60 min/h / 1000 mL/L = 0.06) applied to the CRCL column, not an
    # estimated covariate effect.
    e_crcl_cl <- fixed(0.06); label("Glomerular-filtration contribution to CL per unit CRCL (L/h per mL/min)")  # D'Agate 2024 Table 2 CL equation '+ CLCR' at unit slope; 0.06 = 60/1000 converts the mL/min CRCL column to L/h

    # Systemic infection on V. Table 2 prints 'V (L) = theta2*(WT/1.37)*DIS*theta4',
    # which taken literally would set V to zero for a non-infected subject.
    # The Abstract states the intended form: "Population estimate for volume
    # of distribution was 1.93 L with systemic infection increasing this
    # value by almost 3-fold (2.67 times higher)" -- i.e. theta4^DIS.
    e_infect_vc <- 2.67; label("Fold change in central volume with systemic infection (unitless)")  # D'Agate 2024 Table 2 theta4 = 2.67 (RSE 38%; bootstrap median 2.74, 95% CI 1.46-6.43); orientation from the Abstract

    # IIV. Table 2 reports Omega values on the log scale as VARIANCES; the
    # footnote c confirms it -- "Value in parentheses represents the
    # inter-individual variability of the PK parameters calculated as the
    # square root of (e^Omega - 1) x 100%" -- and sqrt(exp(0.451) - 1) =
    # 75.5% and sqrt(exp(0.646) - 1) = 95.3% reproduce the printed CV%
    # exactly. The reported correlation 0.442 / sqrt(0.451 * 0.646) = 0.819
    # likewise reproduces the printed rho of 82%.
    etalcl + etalvc ~ c(0.451,
                        0.442, 0.646)  # D'Agate 2024 Table 2: Omega1 = 0.451 (75.5% CV, RSE 13%), eta CL-V covariance 0.442 (rho 82%, RSE 42.3%), Omega2 = 0.646 (95.3% CV, RSE 23%)

    # Residual error: proportional only. Table 2 reports sigma1 = 0.143 with
    # (37.8%) alongside; sqrt(0.143) = 0.3782, so the printed value is the
    # variance and the linear-scale proportional SD is its square root.
    propSd <- sqrt(0.143); label("Proportional residual SD on Cc (fraction)")  # D'Agate 2024 Table 2 sigma1 = 0.143 (37.8% CV, RSE 17%) -> propSd = sqrt(0.143) = 0.3782
  })

  model({
    # 1. Derived covariate terms.
    #    Maturation of the residual clearance with post-menstrual age
    #    (D'Agate 2024 Methods Eq. 3). PAGE is carried in WEEKS for this
    #    model -- see covariateData$PAGE$notes.
    maturation_cl <- PAGE^pma_hill / (pma_tm50^pma_hill + PAGE^pma_hill)

    #    Systemic infection multiplies V by e_infect_vc (2.67-fold). Written
    #    as an exponentiated indicator so that DIS_INFECT_ACTIVE = 0 leaves
    #    V unchanged; see the ini() note on the Table 2 typographical form.
    infect_vc <- e_infect_vc^DIS_INFECT_ACTIVE

    # 2. Individual PK parameters.
    #    Total clearance is the sum of the size- and maturation-scaled
    #    residual clearance (tubular secretion plus metabolism, per the
    #    Discussion) and the individual creatinine clearance converted from
    #    mL/min to L/h. The log-normal eta multiplies the SUM, matching the
    #    NONMEM CL = TVCL * exp(eta) convention and the Delattre 2010 /
    #    Georges 2009 additive-CLCR precedents in this package.
    cl <- (exp(lcl) * (WT / 1.37)^e_wt_cl * maturation_cl +
             e_crcl_cl * CRCL) * exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 1.37)^e_wt_vc * infect_vc

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. ODE system. One compartment; doses are given as zero-order
    #    intravenous infusions through the event table (Results Section 4.1:
    #    "one-compartment open model with zero-order infusion and
    #    first-order elimination").
    d/dt(central) <- -kel * central

    # 5. Observation. Dose in mg and vc in L give plasma aciclovir in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
