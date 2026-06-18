AbdelJalil_2013_tacrolimus <- function() {
  description <- "One-compartment population PK model for oral tacrolimus in paediatric liver transplant recipients (Abdel Jalil 2013), with first-order absorption (ka fixed at the literature value 4.5 1/h), an apparent volume of distribution fixed at the literature value 30 L/kg, allometric (WT/13.2 kg)^0.75 scaling on apparent clearance with the theory-based exponent fixed, multiplicative exponential effects of time post-transplantation (days) and CYP3A5*1 carrier status on CL/F, exponential (log-normal) inter-individual variability on CL/F, and proportional residual error."
  reference   <- "Abdel Jalil MH, Hawwa AF, McKiernan PJ, Shields MD, McElnay JC. Population pharmacokinetic and pharmacogenetic analysis of tacrolimus in paediatric liver transplant patients. Br J Clin Pharmacol. 2014;77(1):130-140. doi:10.1111/bcp.12174"
  vignette    <- "AbdelJalil_2013_tacrolimus"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within the year-long observation window per the paper's Methods. Used twice in model(): (i) allometric (WT/13.2)^0.75 scaling on CL/F with the theory-based exponent fixed and the reference set to the study median weight (13.2 kg per Abdel Jalil 2013 Methods 'Selecting the base model'); (ii) linear scaling of V/F because the literature value used to fix V/F is expressed per kg (30 L/kg per Abdel Jalil 2013 Methods 'Selecting the base model'). Study median weight 13.2 kg (range 6.06-69.95 kg).",
      source_name        = "WT"
    ),
    POD = list(
      description        = "Days post-transplantation (paper's TPT, time post-transplantation in days)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject; rises from 0 on the day of liver transplantation. Enters the apparent-clearance equation as a non-centred exponential term: cl_typ *= exp(-0.00158 * POD). The negative coefficient encodes the paper's finding that tacrolimus apparent oral clearance decreases exponentially over the first post-transplant year. Observations were collected over POD 15-364 days (Abdel Jalil 2013 Table 2 demographic summary); model extrapolation outside this window is not supported by the data.",
      source_name        = "TPT"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype derived from rs776746 (CYP3A5 6986A>G; the paper's CYP3A5*3 SNP). In the Abdel Jalil 2013 paediatric cohort the genotype distribution was *1/*1 = 2 (4.7%), *1/*3 = 9 (20.9%), *3/*3 = 32 (74.4%); CYP3A5_EXPR = 1 for the 11 *1 carriers, 0 for the 32 nonexpressers. Multiplicative exponential effect on CL/F: cl_typ *= exp(0.4282 * CYP3A5_EXPR), so expressers have CL/F that is exp(0.4282) = 1.53-fold higher than nonexpressers at the same weight and post-transplant day.",
      source_name        = "CYP3A5"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 43L,
    n_observations   = 628L,
    age_range        = "0.65-17.56 years",
    age_median       = "5 years (mean)",
    weight_range     = "6.06-69.95 kg",
    weight_median    = "21.6 kg (mean); 13.2 kg (study median, used as the allometric reference)",
    sex_female_pct   = 44.2,
    race_ethnicity   = c(White = 72.1, Asian = 20.9, Black = 7.0),
    disease_state    = "Paediatric liver transplant recipients on tacrolimus-based immunosuppression in combination with low-dose steroids; tacrolimus administered orally as capsules or oral suspension.",
    dose_range       = "0.4-30 mg/day total tacrolimus orally (0.03-1.23 mg/kg/day); typical dosing twice daily aiming for therapeutic trough concentrations.",
    regions          = "United Kingdom (single-centre paediatric liver-transplant cohort at Birmingham Children's Hospital).",
    cyp3a5_distribution = "*1/*1 n = 2 (4.7%); *1/*3 n = 9 (20.9%); *3/*3 n = 32 (74.4%). Hardy-Weinberg equilibrium observed.",
    pod_range        = "15-364 days post-transplantation (Abdel Jalil 2013 Table 2).",
    notes            = "Retrospective analysis of tacrolimus pre-dose trough concentrations from therapeutic drug monitoring records; concentrations measured by enzyme immunoassay (Abbott IMx; LLQ 2 ng/mL, linear to 30 ng/mL). Because only pre-dose troughs were available, the volume of distribution and absorption rate constant were not identifiable from the data and were fixed to literature values (V/F = 30 L/kg, ka = 4.5 1/h). The reference weight 13.2 kg is the study median, not 70 kg, so the typical-value clearance estimate theta_CL = 12.92 L/h corresponds to a 13.2-kg paediatric subject (a CYP3A5 nonexpresser at the transplant date)."
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a power covariate on CL/F during forward inclusion (Abdel Jalil 2013 Table 4). Was significant univariately (-2LLD = -52.225) but did not improve the model after weight (allometric) and TPT were included (-2LLD = -0.092 in the second round) and was therefore dropped. Not retained in the final model."
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Tested additively on CL/F (Abdel Jalil 2013 Table 4 'Gender x theta_2'); -2LLD = -0.493, not significant, not retained in the final model."
    ),
    RACE_ASIAN = list(
      description        = "Asian or African-Caribbean race indicator (paper-defined pooling)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian)",
      notes              = "Tested additively on CL/F (Abdel Jalil 2013 Table 4 'Race x theta_2'; coded as 1 for Asian and African-Caribbean patients, 0 for Caucasian patients). Univariate -2LLD = -5.017 (significant), but did not survive forward inclusion after TPT and CYP3A5 (3rd-round -2LLD = -1.935). Not retained in the final model."
    ),
    TX_FULL = list(
      description        = "Whole-liver-graft indicator (vs partial / split graft)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (partial / split liver graft)",
      notes              = "Tested as a multiplicative-on-CL/F effect (Abdel Jalil 2013 Table 4 'Type x theta_2'); -2LLD = -1.71, not significant, not retained in the final model."
    ),
    ABCB1_C1236T = list(
      description        = "ABCB1 C1236T variant carrier indicator (>= 1 variant T allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C)",
      notes              = "Tested additively on CL/F (Abdel Jalil 2013 Table 4 'C1236T x theta_2'). Univariate -2LLD = -21.486 (significant), 2nd-round -2LLD = -0.05 after TPT. Not retained in the final model."
    ),
    ABCB1_G2677T = list(
      description        = "ABCB1 G2677T variant carrier indicator (>= 1 variant T allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (G/G)",
      notes              = "Tested additively on CL/F (Abdel Jalil 2013 Table 4 'G2677T x theta_2'). Univariate -2LLD = -10.156, 2nd-round -2LLD = -1.335. Not retained in the final model."
    ),
    ABCB1_C3435T = list(
      description        = "ABCB1 C3435T variant carrier indicator (>= 1 variant T allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C)",
      notes              = "Tested additively on CL/F (Abdel Jalil 2013 Table 4 'C3435T x theta_2'). Univariate -2LLD = -21.192, 2nd-round -2LLD = -0.796. Not retained in the final model."
    ),
    ABCB1_HAP_TTT = list(
      description        = "ABCB1 T-T-T haplotype (C1236T-G2677T-C3435T) carrier indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other haplotype)",
      notes              = "Tested additively on CL/F (Abdel Jalil 2013 Table 4 'HAP x theta_2'); -2LLD = -1.012, not significant, not retained in the final model."
    ),
    CSD = list(
      description        = "Concomitant corticosteroid dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as an exponential covariate on CL/F (Abdel Jalil 2013 Table 4 'CSD x theta_2'). Univariate -2LLD = -37.99 (significant), 2nd-round -2LLD = -0.006. Not retained in the final model. A moderate inverse correlation between TPT and CSD (r = -0.405) explains why CSD was significant univariately."
    )
  )

  ini({
    # Structural PK -- Abdel Jalil 2013 Table 5 final-model estimates plus
    # the two literature-fixed parameters described in Methods 'Selecting the
    # base model'. Time in hours; CL/F in L/h; V/F in L/kg (multiplied by WT
    # inside model() to get the individual apparent volume in L).
    lka <- fixed(log(4.5)) ; label("Absorption rate constant ka (1/h)")                                                                   # Abdel Jalil 2013 Methods: ka fixed to 4.5 /h (literature value [reference 9])
    lvc <- fixed(log(30))  ; label("Apparent volume of distribution V/F per kg body weight (L/kg)")                                       # Abdel Jalil 2013 Methods: V/F fixed to 30 L/kg (literature value [references 7, 8])
    lcl <- log(12.92)      ; label("Apparent oral clearance CL/F at WT = 13.2 kg, CYP3A5 nonexpresser, day of transplant POD = 0 (L/h)") # Abdel Jalil 2013 Table 5 theta_1 = 12.92 L/h (4.8% RSE)

    # Covariate effects on CL/F -- Abdel Jalil 2013 Results final equation:
    # CL/F_i = theta_1 * (WT/13.2)^0.75 * exp(theta_2 * TPT) * exp(theta_3 * CYP3A5_EXPR).
    e_wt_cl          <- fixed(0.75) ; label("Allometric exponent of (WT/13.2 kg) on CL/F (unitless; fixed)")               # Abdel Jalil 2013 Methods: allometric exponent fixed to 0.75 (theory-based)
    e_pod_cl         <- -0.00158    ; label("Time-post-transplant exponential coefficient on CL/F (per day)")              # Abdel Jalil 2013 Table 5 theta_2 = -0.00158 /day (13.8% RSE)
    e_cyp3a5_expr_cl <- 0.4282      ; label("CYP3A5 expresser exponential coefficient on CL/F (unitless)")                 # Abdel Jalil 2013 Table 5 theta_3 = 0.4282 (25.5% RSE)

    # Inter-individual variability -- Abdel Jalil 2013 Table 5 reports
    # omega^2 (variance on log-scale of clearance) = 0.16 (19.8% RSE).
    # Equivalent approximate CV is sqrt(omega^2) ~ 40% as stated in the
    # paper's Results.
    etalcl ~ 0.16                                                                                                          # Abdel Jalil 2013 Table 5 omega^2 = 0.16

    # Residual unexplained variability -- Abdel Jalil 2013 Table 5 reports
    # sigma^2 = 0.125 (11.9% RSE) on a proportional error model (the
    # additive arm of the combined error structure was driven to ~0 during
    # estimation and dropped; Abdel Jalil 2013 Results 'Population
    # pharmacokinetic analysis'). Proportional CV = sqrt(0.125) ~ 35.4%,
    # consistent with the paper's quoted 35.4% residual variability.
    propSd <- sqrt(0.125) ; label("Proportional residual error (fraction)")                                                # Abdel Jalil 2013 Table 5 sigma^2 = 0.125
  })

  model({
    # Individual PK parameters with the Abdel Jalil 2013 covariate equation.
    # Reference subject: WT = 13.2 kg (study median, NOT 70 kg), CYP3A5
    # nonexpresser (CYP3A5_EXPR = 0), day-of-transplant (POD = 0). V/F is
    # fixed per kg body weight so the apparent volume scales linearly with
    # WT; ka is a literature-fixed constant.
    ka <- exp(lka)
    vc <- exp(lvc) * WT
    cl <- exp(lcl + etalcl) *
          (WT / 13.2) ^ e_wt_cl *
          exp(e_pod_cl * POD) *
          exp(e_cyp3a5_expr_cl * CYP3A5_EXPR)
    kel <- cl / vc

    # One-compartment oral PK with first-order absorption and elimination.
    # Dose lands in `depot`; the paper did not estimate bioavailability and
    # absorbs F into the apparent parameters CL/F and V/F throughout.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Tacrolimus assay reports whole-blood concentrations in ng/mL. With
    # dose in mg and vc in L, central/vc has units mg/L = ug/mL; multiply
    # by 1000 to convert to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
