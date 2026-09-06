NassarSheikhRashid_2024_adalimumab <- function() {
  description <- paste(
    "Two-compartment population PK model with sequential zero- then",
    "first-order subcutaneous absorption and linear elimination for",
    "adalimumab in 50 children with juvenile idiopathic arthritis (JIA)",
    "treated at a single Dutch centre, fitted to 78 therapeutic-drug-",
    "monitoring samples drawn during routine clinical care. The structural",
    "model and every disposition parameter are carried over unchanged from",
    "the best-fitting literature model, Kang 2020 model 3 (the rheumatoid",
    "arthritis maintenance-phase model); only apparent clearance, its",
    "covariate effects, the IIV on clearance and the residual error were",
    "re-estimated on the JIA data. Standard allometric scaling on body",
    "weight is applied to CL/F, Q/F (exponent 0.75) and V1/F, V2/F",
    "(exponent 1) at a 70 kg reference. Apparent clearance additionally",
    "carries anti-drug antibodies (fixed, 108% higher CL), concomitant",
    "methotrexate (28% lower CL), C-reactive protein (power form) and",
    "comorbid uveitis during adalimumab treatment (44% higher CL) - the",
    "last of these a novel association reported by this paper. Together",
    "these covariates cut inter-patient variability in CL/F from 58.6% to",
    "28.0% CV. Residual error is additive only. Upstream model:",
    "modellib('Kang_2020_adalimumab_phase3_extension').",
    sep = " "
  )
  reference <- paste(
    "Nassar-Sheikh Rashid A, Hooijberg F, Bergkamp SC, Gruppen MP,",
    "Kuijpers TW, Nurmohamed M, Rispens T, Wolbink G, van den Berg JM,",
    "Schonenberg-Meinema D, Mathot RAA. Population pharmacokinetics of",
    "adalimumab in juvenile idiopathic arthritis patients: a retrospective",
    "cohort study using clinical care data.",
    "Paediatr Drugs. 2024;26(4):441-450. doi:10.1007/s40272-024-00629-7.",
    "Parameter estimates from Table 2; covariate functional forms from",
    "Equations 1 and 2 and Section 3.3. Disposition parameters and their",
    "IIV are fixed to Kang J, Eudy-Byrne RJ, Mondick J, Knebel W,",
    "Jayadeva G, Liesenfeld K-H. Br J Clin Pharmacol. 2020;86(11):",
    "2274-2285. doi:10.1111/bcp.14330, model 3; see",
    "modellib('Kang_2020_adalimumab_phase3_extension').",
    sep = " "
  )
  vignette <- "NassarSheikhRashid_2024_adalimumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "adalimumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "adalimumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "adalimumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalised to a 70 kg reference. Standard allometric scaling was applied to CL/F, V1/F, Q/F and",
        "V2/F with the exponents held at 0.75 (clearances) and 1 (volumes) rather than estimated",
        "(Nassar-Sheikh Rashid 2024 Section 3.3 'Standard allometric scaling was used ... with a reference",
        "weight of 70 kg', and the Table 2 row labels (WT/70)^0.75 and (WT/70)^1). Median body weight in the",
        "JIA cohort was 49 kg (IQR 29.4-59.8 kg), so the 70 kg reference sits well above the study population",
        "and the typical CL/F of 0.374 L/day is an extrapolated adult-sized value rather than an observed",
        "paediatric one. Body weight also determines the approved dose band (20 mg every other week for 10-30",
        "kg, 40 mg every other week above 30 kg), which is why Figure 2A shows the 30-40 kg stratum reaching",
        "roughly twice the concentrations of the above-60 kg stratum.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    ADA_POS = list(
      description        = "Anti-drug-antibody positivity during adalimumab treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative; the reference covariate value for the typical CL/F per Equation 1)",
      notes              = paste(
        "Screened as 'ever detection of ADA (yes/no)' (Section 2.3), so the column is subject-level rather",
        "than sample-level. ADA were only assayed when the adalimumab concentration fell below 5 mg/L, using",
        "a drug-sensitive antigen-binding test; ADA were detected in 9 of 50 patients (18%, Table 1). The",
        "covariate effect could not be identified stably on these data and was FIXED at 2.08 (Section 3.3).",
        "That value is not a JIA estimate: it was back-derived from Kang 2020 model 3, whose reference group",
        "is ADA-POSITIVE at a titre of 16 and in which ADA-negative patients clear 35% more slowly (a factor",
        "of 0.654, i.e. about 1/1.53). Because Kang additionally carried ADA titre as a separate covariate,",
        "this paper doubled that effect ('the covariate effect of ADA was effectively doubled'). It is the",
        "fractional INCREASE that is doubled, not the ratio: inverting 0.65 gives a 53% higher clearance in",
        "ADA-positive patients, and 2 x 53% = 106%, i.e. a factor of about 2.08 - which is exactly the '108%",
        "increase' the paper quotes. Doubling the ratio itself would instead give 2 x 1.53 = 3.06, which is",
        "not the reported value. The orientation is therefore inverted relative to",
        "modellib('Kang_2020_adalimumab_phase3_extension'), which stores an ADA-negative indicator against an",
        "ADA-positive reference; here the canonical ADA_POS orientation is used directly with an ADA-negative",
        "reference, matching Equation 1. The authors themselves flag the imputation as a limitation (Section",
        "4), because assay-to-assay differences in ADA detection may misstate the true effect size.",
        sep = " "
      ),
      source_name        = "ADA"
    ),
    CONMED_MTX = list(
      description        = "Concomitant methotrexate use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant methotrexate)",
      notes              = paste(
        "Screened as 'use of methotrexate (yes/no)' (Section 2.3) and retained in the final model. 39 of 50",
        "patients (78%) were on methotrexate (Table 1). Methotrexate lowers apparent clearance by 28% (Table",
        "2 theta_MTX = 0.720), taking the typical 70 kg CL/F from 0.374 to 0.269 L/day - the paper quotes",
        "0.27 L/day in Section 4. The authors attribute the effect to methotrexate suppressing ADA formation.",
        "The JIA data set could not resolve a methotrexate dose-response, so the column is a plain yes/no",
        "indicator with no dose term.",
        sep = " "
      ),
      source_name        = "MTX"
    ),
    CRP = list(
      description        = "C-reactive protein, a marker of inflammatory disease activity",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Standard (not high-sensitivity) assay. Time-varying: the paper applies last-observation-carried-",
        "forward to missing time-varying covariates and imputes the population median otherwise (Section",
        "2.5). Normalised to a 0.6 mg/L reference in the power form of Equation 2; the reference is the",
        "population median rounded from the 0.55 mg/L reported in Table 1 (the value 0.6 is what appears in",
        "the Table 2 row label and in the Section 3.3 worked example). CRP values in this cohort are very low",
        "- median 0.55 mg/L, IQR 0.3-1.2 mg/L - which is why the authors ran a sensitivity analysis dropping",
        "CRP entirely (Supplementary Table 3); it changed CL/F and the other covariate effects negligibly and",
        "worsened the OFV by only 10.5 points. Beware extrapolating this covariate far outside the observed",
        "range: the paper's own worked example only goes as far as 10 mg/L, where CL/F is 59% higher than at",
        "the reference. This 0.6 mg/L reference is much lower than the references used in the adult",
        "inflammatory-disease adalimumab and IBD models registered elsewhere in the library (3-15 mg/L).",
        sep = " "
      ),
      source_name        = "CRP"
    ),
    DIS_UVEITIS = list(
      description        = "Active uveitis during adalimumab treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no active uveitis during adalimumab treatment)",
      notes              = paste(
        "Screened as 'uveitis during treatment with adalimumab (yes/no)' (Section 2.3). Eight of 50 patients",
        "(16%) had active uveitis, three of whom also had active arthritis (Table 1). Patients with",
        "concomitant uveitis have 44% higher apparent clearance (Table 2 theta_UV = 1.44). This is the",
        "paper's novel finding and matches the clinical observation that JIA patients with uveitis often need",
        "escalated (up to weekly) adalimumab dosing. The authors caution that ascertainment was not random -",
        "concentrations were mostly measured when treatment appeared clinically ineffective, which may have",
        "preferentially sampled uveitis patients with low levels (Section 4).",
        sep = " "
      ),
      source_name        = "UVEITIS"
    )
  )

  # Covariates the paper screened in its stepwise forward-backward regression
  # but did not retain in the final model (Section 2.3 lists the full screen).
  # They are recorded for provenance only and are deliberately absent from
  # model(); checkModelConventions() treats this list as documentation.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as 'sex' (Section 2.3) and not retained. 36 of 50 patients (72%) were female (Table 1).",
        "Sex was a clearance covariate in the Ternant 2015 adult rheumatoid arthritis model that this paper",
        "evaluated, but it did not reach the OFV drop of 3.84 in the JIA data.",
        sep = " "
      ),
      source_name = "SEX"
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened as an alternative size descriptor to body weight (Section 2.3) and not retained; body",
        "weight with standard allometric exponents was kept instead. No BSA summary is reported in Table 1.",
        sep = " "
      ),
      source_name = "BSA"
    ),
    CONMED_IMMUNOMOD = list(
      description = "Any concomitant immunosuppressive medication",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as 'use of concomitant immunosuppressive medication (yes/no)' (Section 2.3) and not",
        "retained; the methotrexate-specific indicator CONMED_MTX was retained instead. 45 of 50 patients",
        "(90%) were on some concomitant medication: methotrexate 78%, steroids 12%, azathioprine 6%,",
        "leflunomide 2% (Table 1). Note this paper's composite is broader than the registered",
        "CONMED_IMMUNOMOD definition (purine analogue or methotrexate) because it also counts leflunomide and",
        "steroids.",
        sep = " "
      ),
      source_name = "IMMUNOSUPPRESSIVE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50,
    n_studies      = 1,
    n_observations = 78,
    age_range      = "under 18 years by inclusion criterion; no minimum or maximum reported",
    age_median     = "mean 11.8 years (SD 3.9)",
    weight_range   = "not reported; IQR 29.4-59.8 kg",
    weight_median  = "49 kg",
    sex_female_pct = 72,
    race_ethnicity = "not reported",
    disease_state  = paste(
      "juvenile idiopathic arthritis; 68% with active joint inflammation and 16% with active uveitis at the",
      "time of drug measurement",
      sep = " "
    ),
    dose_range     = paste(
      "subcutaneous adalimumab at the approved European JIA dose: 20 mg every other week for 10-30 kg and 40",
      "mg every other week above 30 kg. Most patients received 40 mg every other week.",
      sep = " "
    ),
    regions        = "Netherlands (single centre: Emma Children's Hospital, Amsterdam UMC)",
    notes          = paste(
      "Retrospective single-centre chart review of routine therapeutic drug monitoring, screened January",
      "2013 to August 2023 (Section 2.1). Baseline demographics are Table 1. Sampling was extremely sparse:",
      "78 concentrations from 50 patients, with 33 patients contributing only one sample, and samples were",
      "drawn at random times within the dosing interval rather than to a PK schedule. Serum adalimumab was",
      "measured by a Sanquin ELISA with an LLOQ of 0.01 mg/L; 7 samples (9.0%) were below the detection",
      "limit and were set to 0.005 mg/L. Median observed adalimumab level was 12.0 mg/L (IQR 6.12-15.8).",
      "This sparseness is why only CL/F, its covariates, IIV on CL/F and the residual error were estimated -",
      "everything else is fixed from Kang 2020 model 3 (Section 3.3). eta-shrinkage on CL was 27%, and a",
      "1000-sample bootstrap completed successfully in 885 runs (88.5%).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. All disposition parameters are apparent
    # (CL/F, V1/F, Q/F, V2/F) because only subcutaneous data were
    # available. The reference subject is a 70 kg, ADA-negative,
    # methotrexate-naive, uveitis-free patient with CRP 0.6 mg/L
    # (Nassar-Sheikh Rashid 2024 Equations 1-2 and Table 2 row labels).
    #
    # Every parameter marked with * in Table 2 was fixed to Kang 2020
    # model 3 and is wrapped in fixed() here. Kang reports in hours, this
    # paper in days; the Table 2 values below are the paper's own day-scale
    # numbers and reproduce Kang's exactly (Q/F 0.0709 /h x 24 = 1.70
    # L/day; ka 0.0109 /h x 24 = 0.262 /day; D1 2.91 h / 24 = 0.121 day).
    # ------------------------------------------------------------------
    lcl <- log(0.374)        ; label("Apparent clearance CL/F (L/day)")                            # Table 2: CL/F, L/day (WT/70)^0.75 = 0.374 (RSE 10.0%; bootstrap median 0.372, 95% CI 0.303-0.547)
    lvc <- fixed(log(2.67))  ; label("Apparent central volume of distribution V1/F (L)")           # Table 2: V1/F*, L (WT/70)^1 = 2.67, fixed to Kang 2020 model 3
    lq  <- fixed(log(1.7))   ; label("Apparent inter-compartmental clearance Q/F (L/day)")         # Table 2: Q/F*, L/day (WT/70)^0.75 = 1.7, fixed to Kang 2020 model 3
    lvp <- fixed(log(3.94))  ; label("Apparent peripheral volume of distribution V2/F (L)")        # Table 2: V2/F*, L (WT/70)^1 = 3.94, fixed to Kang 2020 model 3
    lka <- fixed(log(0.262)) ; label("First-order absorption rate constant ka (1/day)")            # Table 2: ka*, day^-1 = 0.262, fixed to Kang 2020 model 3
    ld1 <- fixed(log(0.121)) ; label("Duration of the zero-order input into the depot D1 (day)")   # Table 2: D1*, day = 0.121, fixed to Kang 2020 model 3

    # No intravenous arm was studied and Table 2 reports every disposition
    # parameter in its apparent form, so bioavailability is unidentifiable
    # and is carried at 1 as a structural anchor rather than a paper value.
    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (unitless)")                  # Table 2 reports CL/F, V1/F, Q/F and V2/F throughout; F not estimated

    # Allometric exponents on body weight, held at the canonical 0.75 / 1.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on (WT/70) for CL/F (unitless)")           # Table 2 row label (WT/70)^0.75; Section 3.3 "Standard allometric scaling was used ... with a reference weight of 70 kg"
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on (WT/70) for Q/F (unitless)")            # Table 2 row label (WT/70)^0.75; Section 3.3
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on (WT/70) for V1/F (unitless)")           # Table 2 row label (WT/70)^1; Section 3.3
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on (WT/70) for V2/F (unitless)")           # Table 2 row label (WT/70)^1; Section 3.3

    # ------------------------------------------------------------------
    # Covariate model on CL/F. Equation 1 is the categorical power form
    # TVP = theta_p * theta_cov^COV, Equation 2 the continuous power form
    # TVP = theta_p * (COV/refCOV)^theta_cov. Both are additive on the log
    # scale, so each categorical multiplier below is stored as its natural
    # log and the continuous exponent as-is.
    # ------------------------------------------------------------------
    e_ada_pos_cl     <- fixed(log(2.08)) ; label("Effect of anti-drug-antibody positivity on log CL/F (log-scale)")  # Table 2: (theta_ADA)^ADA* = 2.08, fixed (Section 3.3: fixed because of limited power to identify it; derived from Kang 2020 model 3 by doubling the fractional increase implied by the 0.654 ADA-negative factor: 1/0.65 = 1.53, and 2 x 53% = 106%, i.e. the "108% increase" quoted in Section 3.3)
    e_conmed_mtx_cl  <- log(0.720)       ; label("Effect of concomitant methotrexate on log CL/F (log-scale)")       # Table 2: (theta_MTX)^MTX = 0.720 (RSE 10.8%; bootstrap median 0.723, 95% CI 0.472-0.903) = 28% lower CL
    e_crp_cl         <- 0.165            ; label("Power exponent on (CRP/0.6) for CL/F (unitless)")                  # Table 2: (CRP/0.6)^theta_CRP = 0.165 (RSE 31.6%; bootstrap median 0.175, 95% CI 0.0421-0.323); Section 3.3 worked example: (10/0.6)^0.165 = 1.59, a 59% increase at CRP 10 mg/L
    e_dis_uveitis_cl <- log(1.44)        ; label("Effect of active uveitis during treatment on log CL/F (log-scale)") # Table 2: (theta_UV)^UVEITIS = 1.44 (RSE 16.4%; bootstrap median 1.48, 95% CI 0.987-2.94) = 44% higher CL

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 2 reports %CV, so the internal
    # variances below are omega^2 = log(CV^2 + 1). Each fixed value
    # reproduces its Kang 2020 model 3 variance exactly, which is the
    # cross-check that the CV convention is the right one:
    #   V1/F 89.2% -> 0.585   Q/F 165% -> 1.31
    #   V2/F 40.3% -> 0.150   ka   74.2% -> 0.439
    #
    # Only the CL/F variance was estimated here: Section 3.3 states that
    # "interindividual variability (IIV) was estimated for CL while keeping
    # the full omega variance-covariance matrix for the other parameters as
    # described previously [Kang 2020]". The four fixed etas are therefore a
    # single BLOCK(4) carried over intact, and etalcl sits outside it.
    #
    # NON-PAPER PROVENANCE: Table 2 prints only the diagonal %CV values.
    # The six off-diagonal covariances below are the corresponding 4x4
    # sub-block of Kang 2020 Table 4's BLOCK(5), read from the registered
    # upstream extraction modellib('Kang_2020_adalimumab_phase3_extension').
    # They are what "the full omega variance-covariance matrix ... as
    # described previously" refers to; dropping them would contradict the
    # paper's own description of the model. The resulting 4x4 is positive
    # definite (smallest eigenvalue 0.064). Ordering is lower-triangular
    # row-major in the block order V1, Q, V2, ka.
    # ------------------------------------------------------------------
    etalcl ~ 0.0754785                                        # Table 2: IIV CL/F = 28.0% CV -> log(1 + 0.280^2) = 0.0754785 (bootstrap median 25.8%, 95% CI 7.57-42.6%); down from 58.6% in the covariate-free model
    etalvc + etalq + etalvp + etalka ~ fixed(c(
      0.585,
     -0.205,  1.31,
     -0.138,  0.287, 0.150,
      0.0438, 0.559, 0.119, 0.439
    ))                                                        # Table 2 diagonals IIV V1/F* 89.2%, Q/F* 165%, V2/F* 40.3%, ka* 74.2%; off-diagonals from Kang 2020 Table 4 COV rows via modellib('Kang_2020_adalimumab_phase3_extension')

    # ------------------------------------------------------------------
    # Residual error: additive only (Section 3.3 "residual variability was
    # modelled using an additive error model"). Table 2 prints the VARIANCE
    # 13.0 mg^2/L^2 with its square root alongside as "SD 3.61 mg/L"; the
    # SD is what nlmixr2 wants. This residual is large next to a typical
    # trough of ~7 mg/L, which the authors call out as the main obstacle to
    # using the model directly for therapeutic drug monitoring (Section 4).
    # ------------------------------------------------------------------
    addSd <- 3.61 ; label("Additive residual error (mg/L)")   # Table 2: Additive, mg/L = 13.0 (SD 3.61 mg/L); bootstrap median variance 12.0, 95% CI 3.54-23.4
  })

  model({
    # 1. Individual parameters. Equation 1 (categorical power form) and
    #    Equation 2 (continuous power form) are both additive on the log
    #    scale. Only CL/F carries covariates beyond allometry; the sparse
    #    JIA data could not support covariates elsewhere (Section 3.3).
    cl <- exp(lcl + etalcl +
                e_ada_pos_cl * ADA_POS +
                e_conmed_mtx_cl * CONMED_MTX +
                e_dis_uveitis_cl * DIS_UVEITIS +
                e_crp_cl * log(CRP / 0.6)) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Sequential zero- then first-order subcutaneous absorption, carried
    #    over from Kang 2020 model 3: the dose enters the depot as a
    #    zero-order input over D1 days and then transfers to central
    #    first-order at ka. Dose records must set rate = -2 so that rxode2
    #    takes the input duration from dur(depot).
    f(depot)   <- exp(lfdepot)
    dur(depot) <- d1

    # 5. Observation. central is in mg and vc in L, so central / vc is mg/L,
    #    the unit in which Table 2 reports the additive residual SD.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
