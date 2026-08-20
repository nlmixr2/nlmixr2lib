Chen_2024_rituximab <- function() {
  description <- "Two-compartment population PK of rituximab coupled to a CD19+ B-lymphocyte turnover PD model with saturable stimulation of B-cell loss, in pediatric patients with frequent-relapsing or steroid-dependent nephrotic syndrome"
  reference <- paste(
    "Chen Y, Shen Q, Xiong Y, Dong M, Xu H, Li Z (2024).",
    "Using real-world data to inform dosing strategies of rituximab for pediatric",
    "patients with frequent-relapsing or steroid-dependent nephrotic syndrome:",
    "a prospective pharmacokinetic-pharmacodynamic study.",
    "Front Pharmacol 14:1319744. doi:10.3389/fphar.2023.1319744.",
    "The PK sub-model (structural parameters and the BSA covariate model) was",
    "developed in, and is fixed from, Chen Y, Shen Q, Dong M, Xiong Y, Xu H, Li Z",
    "(2021). Population pharmacokinetics of rituximab in pediatric patients with",
    "frequent-relapsing or steroid-dependent nephrotic syndrome.",
    "Front Pharmacol 12:725665. doi:10.3389/fphar.2021.725665.",
    sep = " "
  )
  vignette <- "Chen_2024_rituximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Model structure (Chen 2024 Eqs 1-4, reproduced verbatim from the PDF):
  #   (1) C1     = A(1) / V1
  #   (2) dA1/dt = k21 * A(2) - k12 * A(1) - CL * C1
  #   (3) dA2/dt = -k21 * A(2) + k12 * A(1)
  #   (4) dA3/dt = Kin - Kout * (1 + EMAX * C1 / (EC50 + C1)) * A(3)
  # A3 is the CD19+ lymphocyte count. Chen 2024 Methods: "In absence of
  # rituximab, baseline CD19+ count is given by kin/kout, from which we derived
  # kin by estimating baseline and kout as parameters. The baseline parameter
  # was used for initializing CD19+ lymphocytes compartment." Hence
  # kin = rbase * kout and bcell(0) = rbase.

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "rituximab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "rituximab", units = "mg", specimen = "plasma", verified = FALSE),
    bcell       = list(analyte = "CD19+ B-lymphocytes", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference (normalising) BSA is 0.9 m^2, the median BSA of the",
        "14-patient cohort (Chen 2021 Table 1). CL scales as (BSA/0.9)^1.26 and",
        "V1 scales linearly as (BSA/0.9); Q and V2 carry no covariate",
        "(Chen 2021 Table 2). The BSA computation formula is unspecified in",
        "both publications. BSA also sets the dose, which is prescribed in",
        "mg/m^2 (375 mg/m^2, capped at 500 mg per infusion).",
        sep = " "
      ),
      source_name        = "BSA"
    )
  )

  # Covariates screened against the PD parameters in Chen 2024 but not retained:
  # "None of the covariates investigated exhibited a statistically significant
  # impact on the depletion of CD19+ lymphocytes" (Chen 2024 Results). The
  # screen also covered disease-state descriptors with no canonical column name
  # (age at onset, nephrotic-syndrome duration, proteinuria status and
  # proteinuria-to-creatinine ratio); those are recorded in population$notes
  # rather than here.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as a PD covariate in Chen 2024; not retained. The authors",
        "attribute the null result to the limited sample size and note that",
        "age-related decreases in CD19+ B cells are reported elsewhere",
        "(Morbach 2010).",
        sep = " "
      )
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = paste(
        "Screened as a PD covariate in Chen 2024; not retained. The cohort was",
        "13 boys / 1 girl (Chen 2021 Table 1).",
        sep = " "
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    TCHOL = list(
      description = "Serum total cholesterol",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance (Schwartz equation)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    CYSC = list(
      description = "Cystatin C",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened as a PD covariate in Chen 2024; not retained."
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 14,
    n_studies     = 1,
    age_range     = "3.0-15.6 years",
    age_median    = "6.8 years",
    weight_range  = "15.0-96.5 kg",
    weight_median = "23.2 kg",
    bsa_range     = "0.6-2.1 m^2",
    bsa_median    = "0.9 m^2",
    sex_female_pct = 7.1,
    race_ethnicity = c(Asian = 100),
    disease_state = "frequent-relapsing or steroid-dependent nephrotic syndrome (5 FRNS, 7 SDNS, 2 FRNS/SDNS; 8 minimal change disease, 3 focal segmental glomerulosclerosis, 3 not biopsied)",
    dose_range    = "375 mg/m^2 IV (maximum 500 mg per infusion) once weekly for up to 2 weeks; 11 of 14 patients received 2 infusions, 3 received 1 infusion",
    regions       = "China (single centre, Children's Hospital of Fudan University, Shanghai)",
    n_observations = "72 rituximab concentrations (Chen 2021); 102 CD19+ lymphocyte counts (Chen 2024)",
    notes = paste(
      "Baseline demographics are in Chen 2021 Table 1; the CD19+ lymphocyte",
      "summary is in Chen 2024 Table 1 (baseline CD19+ median 548.0",
      "[258.4-701.6] x10^6/L; median time to nadir 37.5 days; median time to",
      "B-cell recovery 154.5 days). Race is not reported in either paper;",
      "recorded as Asian on the basis of the single-centre Shanghai cohort.",
      "Median infusion duration was 4.75 h (range 3-6.5 h). Anti-rituximab",
      "antibodies were negative in all 14 patients. Disease-state covariates",
      "screened against the PD parameters but with no canonical register entry:",
      "age at onset (median 3.9 years), nephrotic-syndrome duration (median",
      "31.5 months), proteinuria status (3 positive / 11 negative) and",
      "proteinuria-to-creatinine ratio (median 0.1). None were retained.",
      sep = " "
    )
  )

  ini({
    # --- PK sub-model -------------------------------------------------------
    # All PK structural parameters were fixed in the sequential PK-PD fit;
    # Chen 2024 Table 2 marks each "(fixed)" with footnote a, "Parameters
    # derived from our previous publication (Chen et al., 2021)".
    lcl <- fixed(log(8.69)); label("Clearance at the reference BSA (CL, mL/h)")     # Chen 2024 Table 2 (= Chen 2021 Table 2 theta1)
    lvc <- fixed(log(1.86)); label("Central volume at the reference BSA (V1, L)")   # Chen 2024 Table 2 (= Chen 2021 Table 2 theta2)
    lq <- fixed(log(7.5)); label("Inter-compartmental clearance (Q, mL/h)")         # Chen 2024 Table 2 (= Chen 2021 Table 2 theta3)
    lvp <- fixed(log(1.9)); label("Peripheral volume (V2, L)")                      # Chen 2024 Table 2 (= Chen 2021 Table 2 theta4)

    e_bsa_cl <- fixed(1.26); label("Power exponent of BSA on CL (unitless)")        # Chen 2024 Table 2 theta_BSA~CL (= Chen 2021 Table 2 theta5)
    e_bsa_vc <- fixed(1); label("Power exponent of BSA on V1 (unitless)")           # Chen 2021 Table 2: V1 = theta2 * (BSA/0.9), i.e. structurally linear in BSA

    # --- PD sub-model -------------------------------------------------------
    # Estimated in the sequential PK-PD fit (Chen 2024 Table 2, PD block).
    lemax <- log(99.6); label("Maximum fractional increase in the CD19+ elimination rate (EMAX, unitless)")  # Chen 2024 Table 2, RSE 84.7%
    lec50 <- log(5.87); label("Rituximab concentration at half-maximal effect (EC50, ug/mL)")                # Chen 2024 Table 2, RSE 50.3%
    lrbase <- log(395); label("Baseline CD19+ lymphocyte count (BSLN, 10^6/L)")                              # Chen 2024 Table 2, RSE 124.8%
    lkout <- log(0.051); label("CD19+ lymphocyte elimination rate constant (KOUT, 1/day)")                   # Chen 2024 Table 2, RSE 136.3%

    # --- Inter-individual variability --------------------------------------
    # Chen 2024 Table 2 reports IIV as CV%: CL 37.0%, V1 25.4%, BSLN 65.8%,
    # KOUT 65.4%. Exponential (log-normal) IIV, so omega^2 = log(CV^2 + 1).
    # Q, V2, EMAX and EC50 are reported "NE" (not estimated) and carry no IIV.
    etalcl ~ 0.128305; label("IIV on CL")            # Chen 2024 Table 2: CV 37.0%; log(0.370^2 + 1)
    etalvc ~ 0.062520; label("IIV on V1")            # Chen 2024 Table 2: CV 25.4%; log(0.254^2 + 1)
    etalrbase ~ 0.359745; label("IIV on BSLN")       # Chen 2024 Table 2: CV 65.8%; log(0.658^2 + 1)
    etalkout ~ 0.356076; label("IIV on KOUT")        # Chen 2024 Table 2: CV 65.4%; log(0.654^2 + 1)

    # --- Residual error -----------------------------------------------------
    # "Residual variability in both serum concentration and CD19+ lymphocytes
    # was characterized using a proportional error model" (Chen 2024 Results).
    propSd <- 0.18; label("Proportional residual error on rituximab concentration (fraction)")     # Chen 2024 Table 2 sigma_PK, RSE 30.2%
    propSd_Bcell <- 0.84; label("Proportional residual error on CD19+ lymphocyte count (fraction)") # Chen 2024 Table 2 sigma_PD, RSE 49.4%
  })

  model({
    # 1. Reference / unit-conversion constants
    bsaRef <- 0.9    # m^2; median BSA of the Chen 2021 cohort, the normalising
                     # value in CL = theta1 * (BSA/0.9)^theta5 (Chen 2021 Table 2)
    mlhPerLd <- 24 / 1000  # (mL/h) -> (L/day): x24 h/day, /1000 mL/L

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (BSA / bsaRef)^e_bsa_cl * mlhPerLd  # L/day
    vc <- exp(lvc + etalvc) * (BSA / bsaRef)^e_bsa_vc             # L
    q <- exp(lq) * mlhPerLd                                       # L/day
    vp <- exp(lvp)                                                # L

    emax <- exp(lemax)
    ec50 <- exp(lec50)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    kin <- rbase * kout   # Chen 2024 Methods: baseline = kin/kout, so kin = BSLN * KOUT

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    # PK (Chen 2024 Eqs 2-3); CL * C1 is written here as kel * central, which
    # is identical because C1 = central / vc.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Rituximab concentration driving the PD effect (Chen 2024 Eq 1).
    Cc <- central / vc   # mg/L == ug/mL

    # PD: CD19+ B-lymphocyte turnover with saturable stimulation of the
    # elimination rate (Chen 2024 Eq 4).
    d/dt(bcell) <- kin - kout * (1 + emax * Cc / (ec50 + Cc)) * bcell
    bcell(0) <- rbase

    # 5. Observations and error
    Bcell <- bcell   # CD19+ lymphocyte count, 10^6/L
    Cc ~ prop(propSd)
    Bcell ~ prop(propSd_Bcell)
  })
}
