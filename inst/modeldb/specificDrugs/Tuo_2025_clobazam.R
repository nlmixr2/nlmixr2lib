Tuo_2025_clobazam <- function() {
  description <- "Joint parent-plus-metabolite population PK model for oral clobazam and its active metabolite N-desmethylclobazam (norclobazam) in Chinese children with refractory epilepsy (Tuo 2025). Tandem one-compartment disposition: a first-order absorption depot feeds a one-compartment parent, whose entire elimination is routed into a one-compartment metabolite that is then cleared. Absorption was not identifiable from the opportunistic trough-dominated sampling, so Ka is fixed at 1.99 1/h from Jullien 2015; the dosage conversion fraction Fm from clobazam to N-desmethylclobazam was likewise not estimable, so the metabolite clearance and volume are apparent with respect to Fm. Fixed allometric body-weight exponents (0.75 on both clearances, 1 on both volumes, 70 kg reference) scale all four disposition parameters, and CYP2C19 metabolizer phenotype shifts the metabolite clearance only (intermediate and poor metabolizers relative to a normal-metabolizer reference), which is why CYP2C19 poor metabolizers accumulate N-desmethylclobazam without a matching rise in parent exposure. Between-subject variability was estimated on the two clearances only."
  reference   <- "Tuo Y, Yu X, Li S, Wang J, Liu M, Song X, Ma J, Wang Y, Liu Z, Sun D. Population Pharmacokinetics and Model-Informed Precision Dosing of Clobazam Based on the Developmental and Genetic Characteristics of Children with Epilepsy. Pharmaceutics. 2025;17(7):813. doi:10.3390/pharmaceutics17070813"
  vignette    <- "Tuo_2025_clobazam"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte  = "clobazam",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "clobazam",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_ndmclb = list(
      # Holds an Fm-scaled amount: the source could not estimate the
      # clobazam-to-N-desmethylclobazam dosage conversion fraction Fm, so
      # both the metabolite clearance and its volume are reported apparent
      # with respect to Fm. The predicted metabolite CONCENTRATION is
      # nonetheless the true one, because the same Fm divides the state and
      # the volume it is divided by. See the model() block for the algebra.
      analyte  = "N-desmethylclobazam",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling of all four disposition parameters against a 70 kg reference with exponents FIXED at 0.75 (both clearances) and 1 (both volumes); Tuo 2025 Equations (9)-(12) print the exponents inline and Table 2 reports no exponent parameter, so they were not estimated. The cohort weight range is 6.60-73.00 kg (median 20.00 kg), so the 70 kg reference sits at the extreme upper edge of the observed data and the reported typical values are extrapolated adult-standardised values rather than values observed at 70 kg.",
      source_name        = "Weight"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer; *1/*1 -- both CYP2C19_IM = 0 and CYP2C19_PM = 0)",
      notes              = "1 = subject has CYP2C19 IM phenotype (*1/*2, *1/*3, *2/*17 or *3/*17 in Tuo 2025); 0 = otherwise. Cohort distribution: NM 39.81% (41/103), IM 43.69% (45/103), PM 14.56% (15/103), RM 1.94% (2/103). Affects the metabolite clearance only. The two CYP2C19 rapid metabolizers (*1/*17) were excluded from the covariate analysis for lack of sample size (Tuo 2025 Results 3.1) and no ultrarapid metabolizers (*17/*17) were observed, so neither phenotype has an estimated effect and both fall into the CYP2C19_IM = 0, CYP2C19_PM = 0 reference cell by default; users simulating RM or UM subjects should treat them as an explicit extrapolation.",
      source_name        = "CYP2C19 genotype (NMs / IMs / PMs)"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal, intermediate, or rapid metabolizer)",
      notes              = "1 = subject has CYP2C19 PM phenotype (*2/*2, *2/*3 or *3/*3 in Tuo 2025); 0 = otherwise. Paired with `CYP2C19_IM` to encode the three-level NM (reference) / IM / PM phenotype with two binary indicators. Affects the metabolite clearance only: the paper reports mean post-hoc CL_N-CLB/Fm of 0.46, 0.34 and 0.13 L/h in NMs, IMs and PMs (a 71.7% reduction in PMs vs NMs) against no significant difference in parent CL/F across the three groups.",
      source_name        = "CYP2C19 genotype (NMs / IMs / PMs)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis and not retained; Tuo 2025 Results 3.2 reports that age had no significant impact on the PK parameters of clobazam or N-desmethylclobazam, which the Discussion attributes to the cohort being purely pediatric (0.85-16.75 years) rather than the combined pediatric-plus-adult range of Tolbert 2016."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (Tuo 2025 Results 3.2)."
    ),
    BSA = list(
      description = "Body surface area.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened and not retained (Tuo 2025 Results 3.2); body weight was the size descriptor carried into the final model."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (modified Schwartz formula).",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Screened as part of the renal-function panel and not retained (Tuo 2025 Results 3.2). The Discussion notes the number of patients with impaired renal function was too small to draw firm conclusions."
    ),
    ALB = list(
      description = "Serum albumin.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as part of the hepatic-function panel and not retained (Tuo 2025 Results 3.2)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as part of the hepatic-function panel and not retained (Tuo 2025 Results 3.2)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as part of the hepatic-function panel and not retained (Tuo 2025 Results 3.2)."
    ),
    CONMED_VALPROIC_ACID = list(
      description = "Concomitant valproic acid indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (Tuo 2025 Results 3.2 and Supplementary Figure S1); 80.58% of the cohort received valproic acid. The Discussion notes this agrees with prior modelling analyses that found no clinically relevant antiepileptic drug-drug interaction with clobazam."
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (Tuo 2025 Results 3.2 and Supplementary Figure S1); 28.16% of the cohort."
    ),
    DIET_KETOGENIC = list(
      description = "Adherence to a ketogenic diet.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (Tuo 2025 Results 3.2 and Discussion); only 6.80% of the cohort adhered to a ketogenic diet, which the authors judged too few to resolve an effect despite a published case report of a 42% fall in clobazam and N-desmethylclobazam concentrations on diet initiation."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 3435C>T (rs1045642) genotype.",
      units       = "(genotype)",
      type        = "categorical",
      notes       = "Genotyped and screened; no significant effect on clobazam or N-desmethylclobazam PK (Tuo 2025 Results 3.2 and Supplementary Figure S2). Two further ABCB1 SNPs (rs1128503 1236C>T, rs2032582 2677G>T/A), two CYP3A4 SNPs (rs2740574 *1B, rs2242480 *1G) and six GABA-receptor SNPs (rs2279020, rs279858, rs11503014, rs2229944, rs211014, rs211037) were screened with the same negative result; they are represented here by this single entry rather than one entry each because none carries an estimated coefficient."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 103L,
    n_studies      = 1L,
    n_observations = "156 plasma samples yielding 302 analyte concentrations (154 clobazam + 148 N-desmethylclobazam). Sampling depth per patient: 68 sampled once, 21 twice, 12 three times, 1 four times and 1 six times. Assay quantitative ranges 3-1200 ug/L (clobazam) and 40-16000 ug/L (N-desmethylclobazam) by HPLC-MS/MS.",
    age_range      = "0.85-16.75 years",
    age_median     = "5.46 years (mean 5.94, SD 3.15)",
    weight_range   = "6.60-73.00 kg",
    weight_median  = "20.00 kg (mean 22.94, SD 10.48)",
    sex_female_pct = 43.7,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Pediatric refractory epilepsy: Lennox-Gastaut syndrome, Dravet syndrome, infantile spasms and other refractory epilepsies. 75.73% of patients were taking three or more antiepileptic drugs; 80.58% received concomitant valproic acid, 28.16% lamotrigine, 24.27% perampanel, 18.45% levetiracetam and 17.48% topiramate, and 6.80% adhered to a ketogenic diet.",
    dose_range     = "Oral clobazam tablets, dosed twice daily when the total dose exceeded 5 mg. Starting dose 5 mg for patients weighing 30 kg or less and 10 mg above 30 kg, then individually titrated on efficacy and tolerability.",
    regions        = "China (single centre: Wuhan Children's Hospital, Tongji Medical College, Huazhong University of Science and Technology; enrolment December 2022 to March 2024).",
    genotype       = "CYP2C19 phenotype: normal metabolizers 41 (39.81%), intermediate 45 (43.69%), poor 15 (14.56%), rapid 2 (1.94%). No ultrarapid metabolizers were observed. All genotype frequencies were consistent with Hardy-Weinberg equilibrium.",
    notes          = "Demographics from Tuo 2025 Table 1. Prospective single-centre opportunistic-sampling study using scavenged residual blood drawn for routine biochemistry during safety follow-up, so the sampling is trough-dominated and carries essentially no information on the absorption or distribution phases -- the reason Ka was fixed and a one-compartment rather than two-compartment parent disposition was selected. Reported therapeutic trough ranges applied by the authors are 30-300 ug/L for clobazam and 300-3000 ug/L for N-desmethylclobazam, with laboratory alert levels of 500 ug/L and 5000 ug/L respectively; these are adult-derived targets carried over to the pediatric setting because no pediatric-specific ranges exist."
  )

  ini({
    # ---------------------------------------------------------------
    # All final-model estimates are Tuo 2025 Table 2 ("Final Model,
    # Estimate" column), cross-checked against the printed final-model
    # Equations (8)-(12) on the same page. Bootstrap medians and 95% CIs
    # from the same table (1000 resamples) are quoted alongside each
    # value; every estimate lies inside its bootstrap CI with bias below
    # 5%, and all relative standard errors are below 40%.
    #
    # Typical values are reported standardised to a 70 kg adult, which
    # is at the extreme upper edge of the 6.60-73.00 kg cohort. The
    # cohort median 20 kg child has cl = 5.66 * (20/70)^0.75 = 2.21 L/h.
    # ---------------------------------------------------------------

    # ----- Structural parameters: clobazam (parent) -----
    # Ka was NOT estimated. The opportunistic sampling was almost
    # entirely in the elimination phase, so the authors fixed Ka to the
    # value from Jullien et al. (reference [19] of Tuo 2025), the prior
    # pediatric clobazam popPK model.
    lka <- fixed(log(1.99)); label("Absorption rate constant Ka (1/h)")                                  # Tuo 2025 Table 2 "Ka (h-1) = 1.99 (fixed)" and Equation (8); value carried from Jullien 2015 (ref [19])
    lcl <- log(5.66);        label("Apparent clearance CL_CLB/F of clobazam at 70 kg (L/h)")             # Tuo 2025 Table 2: CL_CLB/F = 5.66, RSE 10.53%, bootstrap median 5.63 (95% CI 4.09-7.26); Equation (10)
    lvc <- log(92.07);       label("Apparent central volume V_CLB/F of clobazam at 70 kg (L)")           # Tuo 2025 Table 2: V_CLB/F = 92.07, RSE 31.15%, bootstrap median 93.08 (95% CI 40.42-176.27); Equation (9)

    # ----- Structural parameters: N-desmethylclobazam (metabolite) -----
    # Both are apparent with respect to Fm, the clobazam-to-N-CLB
    # dosage conversion fraction, which the source states was not
    # estimable ("Since the amount of CLB converted to N-CLB was not
    # clear, Fm was not estimated in this study", Tuo 2025 Methods 2.5).
    lcl_ndmclb <- log(1.01); label("Apparent clearance CL_N-CLB/Fm of N-desmethylclobazam at 70 kg in CYP2C19 normal metabolizers (L/h)")  # Tuo 2025 Table 2: CL_N-CLB/Fm = 1.01, RSE 11.83%, bootstrap median 1.01 (95% CI 0.71-1.30); Equation (12)
    lvc_ndmclb <- log(1.84); label("Apparent central volume V_N-CLB/Fm of N-desmethylclobazam at 70 kg (L)")                                # Tuo 2025 Table 2: V_N-CLB/Fm = 1.84, RSE 29.66%, bootstrap median 1.87 (95% CI 0.72-2.76); Equation (11)

    # ----- Allometric body-weight exponents -----
    # FIXED, not estimated: Equations (9)-(12) print the exponents
    # inline as literal constants (1.0 on the volumes, 0.75 on the
    # clearances) and Table 2 carries no exponent row, no RSE and no
    # bootstrap CI for them. Reference weight 70 kg, likewise printed in
    # the equations rather than estimated.
    e_wt_cl        <- fixed(0.75); label("Allometric exponent on clobazam CL/F (unitless)")                    # Tuo 2025 Equation (10): (Weight/70)^0.75
    e_wt_vc        <- fixed(1);    label("Allometric exponent on clobazam V/F (unitless)")                     # Tuo 2025 Equation (9): (Weight/70)^1.0
    e_wt_cl_ndmclb <- fixed(0.75); label("Allometric exponent on N-desmethylclobazam CL/Fm (unitless)")        # Tuo 2025 Equation (12): (Weight/70)^0.75
    e_wt_vc_ndmclb <- fixed(1);    label("Allometric exponent on N-desmethylclobazam V/Fm (unitless)")         # Tuo 2025 Equation (11): (Weight/70)^1.0

    # ----- CYP2C19 phenotype effects on metabolite clearance -----
    # Equation (12) applies the genotype term as exp(theta_CYP2C19,genotype)
    # multiplying CL_N-CLB/Fm, i.e. a LOG-ADDITIVE shift. Because
    # `lcl_ndmclb` is already on the log scale these coefficients are
    # added there directly, which reproduces the published equation
    # exactly. The normal-metabolizer level is theta = 0.00 (fixed) --
    # the reference cell where both indicators are 0 -- so it needs no
    # parameter of its own. CYP2C19 was NOT retained on the parent
    # clearance (Tuo 2025 Results 3.2: "the CYP2C19 genotype only
    # significantly affected CL/Fm").
    e_cyp2c19_im_cl_ndmclb <- -0.25; label("CYP2C19 intermediate-metabolizer log-additive shift on N-desmethylclobazam CL/Fm (unitless)")  # Tuo 2025 Table 2: theta_CYP2C19,IM = -0.25, RSE 39.03%, bootstrap median -0.24 (95% CI -0.53 to -0.02); Equation (12)
    e_cyp2c19_pm_cl_ndmclb <- -1.30; label("CYP2C19 poor-metabolizer log-additive shift on N-desmethylclobazam CL/Fm (unitless)")          # Tuo 2025 Table 2: theta_CYP2C19,PM = -1.30, RSE 15.01%, bootstrap median -1.28 (95% CI -1.67 to -0.89); Equation (12)

    # ----- Between-subject variability -----
    # Exponential IIV, Tuo 2025 Equation (4): P_i = theta * exp(eta_i),
    # eta_i ~ N(0, omega^2). Table 2 reports the two omega^2 rows as
    # percentages, i.e. omega^2 = 15.94% = 0.1594 and 45.35% = 0.4535
    # (equivalently 41.6% and 75.8% CV). These are VARIANCES expressed
    # as percentages, not CVs -- confirmed against the paper's own
    # Monte Carlo output rather than assumed. Reading them as variances
    # reproduces the Table 3 probability-of-target-attainment column for
    # the 20 kg normal-metabolizer / 0.25 mg/kg BID row (predicted 3.9%
    # and 0.2% of subjects above 300 and 500 ug/L against the published
    # 5.1% and 0.1%); reading them as CVs predicts 0.0005% and 7e-11%,
    # i.e. four to ten orders of magnitude too small. The metabolite row
    # is confirmed the same way against the N-desmethylclobazam columns.
    # No IIV was estimated on Ka or on either volume -- Tuo 2025
    # Discussion paragraph 2 states the simplified model was reached by
    # "reducing model compartments, freezing parameter values, and
    # ignoring interindividual variability for certain parameters".
    # Table 2 reports no correlation block, so the two etas are
    # independent.
    etalcl        ~ 0.1594;  label("IIV on clobazam CL/F (variance on the log scale)")                   # Tuo 2025 Table 2: omega^2 CL_CLB/F = 15.94%, RSE 25.85%, bootstrap median 15.97 (95% CI 7.31-24.63); eta-shrinkage 9.73%
    etalcl_ndmclb ~ 0.4535;  label("IIV on N-desmethylclobazam CL/Fm (variance on the log scale)")        # Tuo 2025 Table 2: omega^2 CL_N-CLB/Fm = 45.35%, RSE 20.20%, bootstrap median 44.96 (95% CI 25.28-64.64); eta-shrinkage 17.42%

    # ----- Residual unexplained variability -----
    # Proportional on the linear scale, Tuo 2025 Equation (5):
    # Y = IPRED * (1 + eps), eps ~ N(0, sigma^2). Table 2 reports the
    # rows as "sigma", i.e. the SD itself and not its square, so these
    # are 33% and 53% proportional error. Separate values per analyte.
    propSd        <- 0.33; label("Proportional residual SD for clobazam (fraction)")                     # Tuo 2025 Table 2: sigma_CLB = 0.33, RSE 8.52%, bootstrap median 0.33 (95% CI 0.27-0.37)
    propSd_ndmclb <- 0.53; label("Proportional residual SD for N-desmethylclobazam (fraction)")          # Tuo 2025 Table 2: sigma_N-CLB = 0.53, RSE 8.60%, bootstrap median 0.53 (95% CI 0.44-0.61); eps-shrinkage 33.17%
  })

  model({
    # ----- Covariate model -----
    # Continuous covariates enter as a power of the ratio to the cohort
    # reference (Tuo 2025 Equation (6)); categorical covariates enter
    # log-additively (Equation (7)). Both are instantiated here by
    # Equations (9)-(12).
    wtr <- WT / 70

    cl        <- exp(lcl + etalcl) * wtr^e_wt_cl
    vc        <- exp(lvc) * wtr^e_wt_vc
    cl_ndmclb <- exp(lcl_ndmclb + etalcl_ndmclb +
                       e_cyp2c19_im_cl_ndmclb * CYP2C19_IM +
                       e_cyp2c19_pm_cl_ndmclb * CYP2C19_PM) * wtr^e_wt_cl_ndmclb
    vc_ndmclb <- exp(lvc_ndmclb) * wtr^e_wt_vc_ndmclb

    ka <- exp(lka)

    # ----- ODE system -----
    # Tuo 2025 Equations (1)-(3):
    #   dAa/dt      = -Ka * Aa * F
    #   dA_CLB/dt   =  Ka * Aa * F - CL_CLB * A_CLB / V_CLB
    #   dA_N-CLB/dt =  CL_CLB * A_CLB * Fm / V_CLB - CL_N-CLB * A_N-CLB / V_N-CLB
    #
    # Neither F (clobazam oral bioavailability) nor Fm (the clobazam-to-
    # N-CLB dosage conversion fraction) was estimated, so the system is
    # written in the apparent parameterisation the paper reports. Divide
    # the parent equations through by F and the metabolite equation by
    # (F * Fm): every state becomes a scaled amount, and the scaling
    # cancels in each observation because the SAME factor divides the
    # clearance, the volume and the state.
    #
    #   depot          holds Aa                (dosed amount, F-scaled)
    #   central        holds A_CLB / F
    #   central_ndmclb holds A_N-CLB / (F * Fm)
    #
    # so `cl` = CL_CLB/F, `vc` = V_CLB/F, `cl_ndmclb` = CL_N-CLB/Fm and
    # `vc_ndmclb` = V_N-CLB/Fm are exactly the four apparent quantities
    # of Table 2, and `Cc` / `Cc_ndmclb` below are the true plasma
    # concentrations. This is the same convention the register records
    # for the sibling `ndmima` (N-desmethyl-imatinib) models.
    #
    # The parent's ENTIRE elimination is routed into the metabolite: the
    # paper writes no separate non-metabolic clobazam elimination arm,
    # and the Fm factor of Equation (3) -- which is what would otherwise
    # split it -- has been absorbed into the metabolite state scaling.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - cl * (central / vc)
    d/dt(central_ndmclb) <-  cl * (central / vc) -
                             cl_ndmclb * (central_ndmclb / vc_ndmclb)

    # ----- Observations -----
    Cc        <- central / vc
    Cc_ndmclb <- central_ndmclb / vc_ndmclb

    Cc        ~ prop(propSd)
    Cc_ndmclb ~ prop(propSd_ndmclb)
  })
}
