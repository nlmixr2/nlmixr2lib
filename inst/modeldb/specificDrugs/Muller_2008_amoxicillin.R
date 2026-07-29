Muller_2008_amoxicillin <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous amoxicillin",
    "in pregnant women before, during and immediately after labour, with",
    "labour-state binary indicators reducing the peripheral volume V2",
    "during active labour (-13.7%) and the immediate postpartum period",
    "(-29.5%) relative to before labour (Muller 2008)."
  )
  reference <- paste(
    "Muller AE, Dorr PJ, Mouton JW, De Jongh J, Oostvogel PM,",
    "Steegers EAP, Voskuyl RA, Danhof M. The influence of labour on the",
    "pharmacokinetics of intravenously administered amoxicillin in",
    "pregnant women. Br J Clin Pharmacol. 2008;66(6):866-874.",
    "doi:10.1111/j.1365-2125.2008.03292.x."
  )
  vignette <- "Muller_2008_amoxicillin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    LABOR_ACTIVE = list(
      description        = paste(
        "Binary indicator that the sample was collected while the subject",
        "is in active labour (uterine contractions with progressive",
        "cervical dilatation per the delivering physician's vaginal exam,",
        "Muller 2008 Methods 'Patients'). Time-varying within a single",
        "subject's PK window as she transitions from before-labour",
        "through active labour to delivery."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in active labour; before-labour or postpartum samples)",
      notes              = paste(
        "Paired with LABOR_POSTPARTUM in a 3-level encoding: both 0",
        "= before labour (reference), LABOR_ACTIVE = 1 / LABOR_POSTPARTUM",
        "= 0 = during active labour, LABOR_ACTIVE = 0 /",
        "LABOR_POSTPARTUM = 1 = immediate postpartum (within ~28 h of",
        "delivery)."
      ),
      source_name        = "LABOUR"
    ),
    LABOR_POSTPARTUM = list(
      description        = paste(
        "Binary indicator that the sample was collected in the immediate",
        "postpartum period. In Muller 2008, postpartum samples were",
        "obtained up to ~28 h after delivery (the additional 1 g",
        "postpartum dose was administered 1.5-3.8 h after child birth;",
        "Results paragraph 1)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in immediate postpartum; before-labour or during-labour samples)",
      notes              = paste(
        "Paired with LABOR_ACTIVE; reference state (both 0) is before",
        "labour. Distinct from TPP (continuous weeks-postpartum) which",
        "models multi-week recovery of pregnancy-induced physiology;",
        "LABOR_POSTPARTUM captures only the acute peri-delivery state."
      ),
      source_name        = "POSTPARTUM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at study entry.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Reported in Muller 2008 Table 1 (79.0 +/- 13.6 kg; range 53-107,",
        "n=33). Screened in the forward-addition covariate analysis on",
        "the parameters with random effects (CL, V1, V2) but not retained",
        "as an independent covariate in the final model (Methods",
        "'Pharmacokinetic analysis'; Results paragraph 3)."
      )
    ),
    BMI = list(
      description = "Body mass index at study entry.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Reported in Muller 2008 Table 1 (28.8 +/- 5.0 kg/m^2; range 18-38, n=33). Screened, not retained."
    ),
    GA = list(
      description = "Gestational age at the time of the PK study.",
      units       = "weeks",
      type        = "continuous",
      notes       = "Reported in Muller 2008 Table 1 (35.9 +/- 2.3 weeks; range 30.0-40.6, n=34). Screened, not retained."
    ),
    OEDEMA = list(
      description = "Semi-quantitative oedema score scored 0 (none) / 1 (around the ankle) / 2 (up to the knee) / 3 (above the knee). Time-fixed at study entry.",
      units       = "(ordinal 0-3)",
      type        = "count",
      notes       = paste(
        "Reported in Muller 2008 Table 1 with cohort distribution",
        "21 / 10 / 2 across scores 0 / 1 / 2 (no patient scored 3).",
        "Muller 2008 reports a significant effect of OEDEMA on the central",
        "volume V1 (Results paragraph 3: 'The volume of distribution",
        "increased with an increasing amount of oedema. (Figure 1b)'),",
        "but the relationship is presented graphically only (Figure 1B,",
        "V1 vs amount of oedema). No numeric coefficient, slope, or",
        "categorical-stratum estimate is reported in Table 2 or the",
        "narrative. The effect is therefore not implemented in this",
        "model; the typical V1 = 8.7 L from Table 2 is carried verbatim.",
        "See the vignette's 'Assumptions and deviations' section for the",
        "downstream consequences."
      )
    ),
    CREAT = list(
      description = "Serum creatinine at study entry.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Reported in Muller 2008 Table 1 (45.1 +/- 8.2 umol/L; range 30-74, n=34). Screened on CL with both raw value and as input to Cockcroft-Gault / MDRD GFR estimators; not retained."
    ),
    CRCL_CG = list(
      description = "Estimated creatinine clearance from the Cockcroft-Gault formula.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL, not retained (Muller 2008 Results paragraph 3: 'Using the serum creatinine concentration and the estimated creatinine clearance calculated with the CG-formula, no influence on the CL of amoxicillin was seen.')."
    ),
    GFR_MDRD = list(
      description = "Estimated glomerular filtration rate from the modified MDRD equation.",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and reported as having a small but statistically",
        "significant inverse association with CL (Muller 2008 Results",
        "paragraph 3, Figure 1C). The paper explicitly judges this",
        "finding implausible ('This is unlikely, and therefore illustrates",
        "the importance of the use of validated formulae in special",
        "patient groups') and notes that 'Pharmacokinetic parameter",
        "estimates did not change after the implementation of GFR",
        "calculated with the MDRD formula'. The effect is therefore not",
        "retained in the published final model and is not implemented",
        "here."
      )
    ),
    TWIN = list(
      description = "Twin (vs singleton) pregnancy indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Cohort 31 singleton + 3 twin pregnancies (Muller 2008 Results paragraph 1). Screened, not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "20-38 years (maternal age; mean 29.0 +/- 5.5)",
    weight_range   = "53-107 kg (mean 79.0 +/- 13.6, n=33)",
    ga_range       = "30.0-40.6 weeks gestational age at the time of the PK study",
    sex_female_pct = 100,
    disease_state  = paste(
      "Pregnant women treated with amoxicillin for prevention of",
      "neonatal group B streptococcus (GBS) disease per local guidelines",
      "(indications: preterm premature rupture of the membranes,",
      "rupture of the membranes >18 h, prematurity, fever >37.8 C,",
      "bacteriuria, prior child with invasive GBS disease).",
      "Cohort comprises 17 patients sampled only before the onset of",
      "labour, 9 sampled only during labour, and 8 sampled both before",
      "and during labour; 8 of the 34 patients additionally received the",
      "elective postpartum dose. 31 singleton and 3 twin pregnancies."
    ),
    dose_range     = paste(
      "Initial 2 g amoxicillin IV infusion over 30 min (50 mg/mL),",
      "followed by 1 g amoxicillin IV infusion over 15 min every 4 h",
      "until delivery; one optional 1 g postpartum dose administered 4 h",
      "after the last antepartum dose."
    ),
    regions        = "Netherlands (single centre, Medical Centre Haaglanden, the Hague)",
    notes          = paste(
      "Subset of the 34-patient cohort overlaps with the prior Muller",
      "study of amoxicillin in PPROM (n=416 samples shared, all collected",
      "before the onset of labour; Methods 'Patients'). Demographic table",
      "is Muller 2008 Table 1; 898 total plasma samples (550",
      "before-labour / 187 during-labour / 161 postpartum)."
    )
  )

  ini({
    # Structural disposition parameters reported as the typical value
    # for the reference state (LABOR_ACTIVE = 0 and LABOR_POSTPARTUM = 0,
    # i.e., before the onset of labour). Source: Muller 2008 Table 2
    # 'Structural model parameters', Mean (model) column.
    lcl  <- log(21.1);  label("Clearance (L/h)")                                # Table 2 row CL
    lvc  <- log(8.7);   label("Central volume V1 (L)")                           # Table 2 row V1
    lvp  <- log(11.8);  label("First peripheral volume V2 (L) at reference state") # Table 2 row V2
    lvp2 <- log(20.5);  label("Second peripheral volume V3 (L)")                  # Table 2 row V3
    lq   <- log(21.9);  label("Intercompartmental clearance Q1 between V1 and V2 (L/h)")  # Table 2 row Q1
    lq2  <- log(1.5);   label("Intercompartmental clearance Q2 between V1 and V3 (L/h)")  # Table 2 row Q2

    # Covariate effects on V2 (peripheral1) for the three-level labour
    # state, encoded as an additive-fractional model with before-labour
    # as the reference category:
    #   V2_indiv = TVV2 * (1 + e_labor_active_vp     * LABOR_ACTIVE
    #                        + e_labor_postpartum_vp * LABOR_POSTPARTUM)
    # Source: Muller 2008 Results paragraph 4 and Discussion paragraph 1
    # ('Compared with women before the onset of labour, V2 was decreased
    # with 13.7% during labour and 29.5% in the immediate postpartum
    # period').
    e_labor_active_vp     <- -0.137;  label("LABOR_ACTIVE fractional effect on V2 (unitless)")      # Results paragraph 4
    e_labor_postpartum_vp <- -0.295;  label("LABOR_POSTPARTUM fractional effect on V2 (unitless)")  # Results paragraph 4

    # Inter-individual variability. The variances below are the
    # 'Mean (model) x 10^-2' entries of Muller 2008 Table 2's
    # 'Variance model parameters' block, i.e., omega^2 on the log scale.
    # Equivalent lognormal CV%:
    #   sqrt(exp(0.042) - 1) = 20.7% (Table 2 reports CV 19.8% bootstrap)
    #   sqrt(exp(0.051) - 1) = 22.9% (Table 2 reports CV 23.1% bootstrap)
    #   sqrt(exp(0.097) - 1) = 31.9% (Table 2 reports CV 31.6% bootstrap)
    # Muller 2008 Results paragraph 3 reports that 'Correlations between
    # the random parameters for interindividual variability were found
    # and implemented in the model', but the numeric off-diagonal
    # entries are not given in Table 2 or elsewhere on disk; a diagonal
    # OMEGA is used here and the gap is flagged in the vignette's
    # 'Assumptions and deviations' section.
    etalcl ~ 0.042   # Table 2 row 'IIV in CL'  (variance, log scale)
    etalvc ~ 0.051   # Table 2 row 'IIV in V1'  (variance, log scale)
    etalvp ~ 0.097   # Table 2 row 'IIV in V2'  (variance, log scale)

    # Residual error. Table 2 reports a proportional residual variance of
    # 0.046 (Mean (model) x 10^-2 = 4.6); the corresponding proportional
    # CV is sqrt(0.046) = 0.2145 (~21.5%, matching the 21.5% bootstrap
    # CV reported in Table 2).
    propSd <- sqrt(0.046);  label("Proportional residual SD on Cc (fraction)")   # Table 2 row 'Residual variability'
  })

  model({
    # Apply the labour-state effect on the typical V2 before adding eta.
    # The reference state (LABOR_ACTIVE = 0 AND LABOR_POSTPARTUM = 0) is
    # before labour; the during-labour and postpartum strata are paired
    # binary indicators (both = 1 is invalid by construction).
    labor_vp_factor <- 1 + e_labor_active_vp     * LABOR_ACTIVE +
                            e_labor_postpartum_vp * LABOR_POSTPARTUM

    # Individual disposition parameters.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp) * labor_vp_factor
    vp2 <- exp(lvp2)
    q   <- exp(lq)
    q2  <- exp(lq2)

    # Micro-constants for the three-compartment model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODE system: IV dose enters the central compartment directly via the
    # event table (no depot; infusion rate is set on the dose record).
    d/dt(central)     <- -(kel + k12 + k13) * central +
                           k21 * peripheral1 +
                           k31 * peripheral2
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-   k13 * central - k31 * peripheral2

    # Observation. Dose in mg / vc in L -> mg/L, matching the paper's
    # reported peak amoxicillin concentrations (e.g., 97.4 mg/L after a
    # 2 g infusion; Muller 2008 Results paragraph 1).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
