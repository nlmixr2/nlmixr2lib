Faessel_2019_pevonedistat <- function() {
  description <- "Two-compartment population PK model for intravenous pevonedistat (TAK-924, MLN4924) in adults with solid tumours or haematological malignancies, with body-surface-area scaling of clearance and volumes, an albumin effect on distributional clearance, and reduced clearance on concomitant carboplatin plus paclitaxel"
  reference <- paste(
    "Faessel HM, Mould DR, Zhou X, Faller DV, Sedarati F, Venkatakrishnan K.",
    "Population pharmacokinetics of pevonedistat alone or in combination with",
    "standard of care in patients with solid tumours or haematological",
    "malignancies. Br J Clin Pharmacol. 2019;85(11):2568-2579.",
    "doi:10.1111/bcp.14078",
    sep = " "
  )
  vignette <- "Faessel_2019_pevonedistat"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling on (BSA / 1.73). A single exponent (1.33) is shared by CL and Q,",
        "and a second single exponent (1.39) is shared by Vc and Vp; Faessel 2019 Table 3",
        "footnotes a and b state that each pair was constrained to be equal after the",
        "unconstrained model produced high standard errors on the BSA effects.",
        "Observed range 1.38-3 m^2 (Results, 'Justification of BSA-based dosing');",
        "per-study medians 1.84-2.01 m^2 (Table 2). The paper does not state which BSA",
        "formula the clinical databases used, so it is recorded as unspecified.",
        sep = " "
      ),
      source_name = "BSA"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling on (ALB / 40) applied to Q only; albumin was not retained on CL,",
        "so it changes the shape of the concentration-time profile without changing AUC",
        "(Faessel 2019 section 4.3). Observed range 20-50 g/L; over that range Q moves",
        "from 8.01 to 30.27 L/h. Per-study medians 35-39.6 g/L (Table 2); albumin was",
        "not reported for studies C15009 and C15010.",
        sep = " "
      ),
      source_name = "ALB"
    ),
    CONMED_CARBOPLATIN = list(
      description = "Concomitant carboplatin coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving concomitant carboplatin)",
      notes = paste(
        "Faessel 2019 estimated a single covariate effect for the carboplatin + paclitaxel",
        "doublet (study C15010 arm 2), not for either agent alone, so the model applies the",
        "effect to the product CONMED_CARBOPLATIN * CONMED_PACLITAXEL. In the analysis",
        "dataset the two agents were only ever given together, so the product reproduces",
        "the paper's indicator exactly while keeping each agent on its own canonical column.",
        sep = " "
      ),
      source_name = "ConCarboTax"
    ),
    CONMED_PACLITAXEL = list(
      description = "Concomitant paclitaxel coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving concomitant paclitaxel)",
      notes = paste(
        "See CONMED_CARBOPLATIN: the effect is estimated for the doublet and is applied to",
        "the product of the two indicators. Concomitant azacitidine, docetaxel and",
        "gemcitabine were tested and did not affect pevonedistat CL (Faessel 2019",
        "sections 3 and 4.2), so those comedications carry no covariate in this model.",
        sep = " "
      ),
      source_name = "ConCarboTax"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Identified early in covariate model development as an important predictor of CL",
        "but supplanted by BSA, which is highly correlated with weight and is the basis of",
        "clinical dosing (Faessel 2019 section 4.4). Not in the final model.",
        sep = " "
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault, capped at 150 mL/min for covariate testing)",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Met the forward-selection criterion (P < .05) in univariate testing but was removed",
        "during backward elimination (P < .01); CrCL explained no more of the IIV in CL than",
        "BSA already did (Faessel 2019 section 4.3). Not in the final model.",
        "NOTE the units: Faessel 2019 Table 2 footnote a reports raw Cockcroft-Gault",
        "clearance in mL/min, whereas the canonical CRCL column is BSA-normalized",
        "(mL/min/1.73 m^2). Because this covariate was screened and rejected it never enters",
        "model(), but any future promotion into covariateData must normalize first.",
        sep = " "
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested on CL; not influential over the range represented (Faessel 2019 section 4.3). Not in the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested on CL; not influential over the range represented (Faessel 2019 section 4.3). Not in the final model."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Tested on CL; not influential. The dataset excluded patients with bilirubin above",
        "1.5 x ULN per protocol, so moderate and severe hepatic impairment are unrepresented",
        "(Faessel 2019 section 4.3). Not in the final model.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units = "y",
      type = "continuous",
      notes = "Tested (range 23-90 years); no effect on any PK parameter (Faessel 2019 sections 3 and 4.4). Not in the final model."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Tested; no effect on pevonedistat CL (Faessel 2019 Figure 3 panel A). Not in the final model."
    ),
    # The three CONMED_* names below are well-formed members of the auto-approved
    # CONMED_<INN> family but are deliberately NOT added to
    # inst/references/covariate-columns.md: Faessel 2019 screened and rejected all
    # three, they are never referenced in model(), and registering a canonical for a
    # negative screen would be register pollution. Register them if a later paper
    # retains one.
    CONMED_AZACITIDINE = list(
      description = "Concomitant azacitidine coadministration",
      units = "(binary)",
      type = "binary",
      notes = "Tested; did not enter the model as a covariate on CL (Faessel 2019 sections 3 and 4.2). Not in the final model."
    ),
    CONMED_DOCETAXEL = list(
      description = "Concomitant docetaxel coadministration",
      units = "(binary)",
      type = "binary",
      notes = "Tested; did not enter the model as a covariate on CL (Faessel 2019 sections 3 and 4.2). Not in the final model."
    ),
    CONMED_GEMCITABINE = list(
      description = "Concomitant gemcitabine coadministration",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Identified as a possible covariate in univariate testing but eliminated in the",
        "stepwise backward removal (Faessel 2019 section 3). Not in the final model.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    central = list(analyte = "pevonedistat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pevonedistat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 346,
    n_studies = 6,
    age_range = "23-90 years",
    age_median = "62.1 years (mean)",
    weight_range = "43.5-180 kg",
    weight_median = "80.3 kg (mean)",
    sex_female_pct = 41,
    race_ethnicity = c(Caucasian = 87),
    disease_state = "advanced solid tumours or haematological malignancies (AML, MDS, lymphoma, multiple myeloma, melanoma)",
    dose_range = "15-278 mg/m^2 as a 1-hour intravenous infusion, across 6 dosing schedules in 21- or 28-day cycles",
    bsa_range = "1.38-3 m^2",
    albumin_range = "20-50 g/L",
    renal_function = "46% normal, 36% mild impairment (CrCL 60-89 mL/min), 18% moderate (CrCL 30-59 mL/min), 1 patient severe at baseline",
    hepatic_function = "total bilirubin <= ULN and ALT/AST <= 2.5 x ULN per protocol exclusion criteria; no moderate or severe hepatic impairment represented",
    co_medication = "single agent, or with azacitidine (C15009) or docetaxel / carboplatin + paclitaxel / gemcitabine (C15010)",
    notes = paste(
      "Faessel 2019 Tables 1 and 2. 346 adult patients were enrolled across studies",
      "C15001, C15002, C15003, C15005, C15009 and C15010; 11 had no PK data, leaving 335",
      "evaluable subjects contributing 3768 concentration observations. Per-study median",
      "BSA ranged 1.84-2.01 m^2 and median serum albumin 35-39.6 g/L (albumin was not",
      "reported for C15009 or C15010).",
      sep = " "
    )
  )

  ini({
    # Structural parameters - typical values for the reference patient
    # (BSA 1.73 m^2, albumin 40 g/L, not receiving carboplatin + paclitaxel;
    # Faessel 2019 Results paragraph following Table 3).
    lcl <- log(31.5)
    label("Clearance (L/h)") # Faessel 2019 Table 3, 'CL (L/h)' = 31.5 (SE 6.4%; bootstrap 95% CI 27.3-35)
    lvc <- log(117)
    label("Central volume of distribution (L)") # Faessel 2019 Table 3, 'Vc (L)' = 117 (SE 8.4%; bootstrap 95% CI 96.7-136)
    lq <- log(21.9)
    label("Intercompartmental clearance (L/h)") # Faessel 2019 Table 3, 'Q (L/h)' = 21.9 (SE 6.6%; bootstrap 95% CI 19.7-24.6)
    lvp <- log(122)
    label("Peripheral volume of distribution (L)") # Faessel 2019 Table 3, 'Vp (L)' = 122 (SE 3.9%; bootstrap 95% CI 112-131)

    # Covariate effects. The two BSA exponents are each shared by a pair of
    # parameters, exactly as Faessel 2019 Table 3 footnotes a and b describe.
    e_bsa_cl_q <- 1.33
    label("Shared power exponent on (BSA / 1.73) for CL and Q (unitless)") # Faessel 2019 Table 3, 'BSA on CL' = 'BSA on Q' = 1.33 (footnote a; SE 12.3%)
    e_bsa_vc_vp <- 1.39
    label("Shared power exponent on (BSA / 1.73) for Vc and Vp (unitless)") # Faessel 2019 Table 3, 'BSA on Vc' = 'BSA on Vp' = 1.39 (footnote b; SE 12.5%)
    e_alb_q <- 1.45
    label("Power exponent on (ALB / 40) for Q (unitless)") # Faessel 2019 Table 3, 'Albumin on Q' = 1.45 (SE 34.7%)
    e_conmed_carboplatin_paclitaxel_cl <- -0.441
    label("Fractional change in CL on concomitant carboplatin + paclitaxel (unitless)") # Faessel 2019 Table 3, 'Carboplatin + paclitaxel on CL' = -0.441 (SE 8.8%); enters as (1 + theta * indicator)

    # IIV. Faessel 2019 Table 3 reports interindividual variability as %CV of a
    # log-normally distributed parameter; omega^2 = log(1 + CV^2).
    # CL is deliberately OUTSIDE the covariance block: the authors removed it
    # from the full 4-parameter block because the condition number rose from
    # 12.51 to 99.93, and the reduced block converged with a condition number of
    # 16.38 (Faessel 2019 Results, second paragraph).
    etalcl ~ 0.125712 # Faessel 2019 Table 3, 'IIV_CL (%CV)' = 36.6; log(1 + 0.366^2) = 0.125712
    # Correlated block on Vc, Q and Vp. Off-diagonals are
    # corr * omega_i * omega_j with omega_Vc = 0.354559, omega_Q = 0.558485,
    # omega_Vp = 0.357288.
    etalvc + etalq + etalvp ~ c(
      0.125712,
      0.136829, 0.311905,
      0.119206, 0.102763, 0.127655
    ) # Faessel 2019 Table 3: 'IIV_Vc(%CV)' = 36.6, 'IIV_Q (%CV)' = 60.5, 'IIV_Vp (%CV)' = 36.9; 'Corr (Vc, Q)' = 0.691, 'Corr (Vc,Vp)' = 0.941, 'Corr (Q,Vp)' = 0.515

    # Residual error. Faessel 2019 section 2.3: 'Residual variability was
    # modelled using the log transform both sides approach with an additive
    # error model', i.e. log-normal (exponential) error on the linear scale.
    expSd <- 0.323360
    label("Log-normal residual error (SD on the log scale)") # Faessel 2019 Table 3, 'Residual error (%CV)' = 33.2; sqrt(log(1 + 0.332^2)) = 0.323360
  })

  model({
    # 1. Derived covariate terms (Faessel 2019 Table 3, displayed equations)
    bsaRatio <- BSA / 1.73 # reference BSA 1.73 m^2
    albRatio <- ALB / 40 # reference albumin 40 g/L
    carboTax <- CONMED_CARBOPLATIN * CONMED_PACLITAXEL # doublet indicator

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * bsaRatio^e_bsa_cl_q *
      (1 + e_conmed_carboplatin_paclitaxel_cl * carboTax)
    vc <- exp(lvc + etalvc) * bsaRatio^e_bsa_vc_vp
    q <- exp(lq + etalq) * bsaRatio^e_bsa_cl_q * albRatio^e_alb_q
    vp <- exp(lvp + etalvp) * bsaRatio^e_bsa_vc_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (intravenous infusion into the central compartment)
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
