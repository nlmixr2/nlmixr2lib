Francis_2019_rifapentine <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order elimination and Savic transit-compartment absorption (NN = 10.2, MTT = 1.47 h) followed by a separate first-order absorption step (ka = 0.814 1/h) for oral rifapentine in 326 southern African adults with drug-susceptible pulmonary tuberculosis pooled from two trials (RIFAQUIN, 900 mg twice weekly or 1,200 mg once weekly in the continuation phase; Daily RPE, 450 or 600 mg daily in the intensive phase). Clearance is allometrically scaled on fat-free mass (reference 46 kg, exponent 0.75 fixed) and apparent volume on total body weight (reference 56 kg, exponent 1 fixed), giving typical values CL = 1.33 L/h and V = 25 L. Four covariate effects are carried: AADAC rs1803155 AA homozygotes have 10.4% lower clearance (the paper's headline pharmacogenetic finding), the RIFAQUIN 1,200-mg-once-weekly arm has 13.2% lower clearance, HIV-positive patients have 21.9% lower bioavailability, and Daily RPE participants have 23.3% lower bioavailability than RIFAQUIN participants (attributed to non-standardised meals). Bioavailability is fixed to 1 in the RIFAQUIN HIV-negative reference. Inter-individual variability is on CL (23.0% CV) and V (12.8% CV); the paper's inter-occasion variability on ka (48.9% CV), MTT (37.4% CV) and NN (20.3% CV) is carried as single-draw random effects because the source does not enumerate occasions - see the vignette Assumptions and deviations. Residual variability is combined additive (0.247 mg/L) plus proportional (9.56%)."
  reference <- paste(
    "Francis J, Zvada SP, Denti P, Hatherill M, Charalambous S, Mungofa S,",
    "Dawson R, Dorman S, Gupte N, Wiesner L, Jindani A, Harrison TS,",
    "Olagunju A, Egan D, Owen A, McIlleron HM. (2019).",
    "A population pharmacokinetic analysis shows that arylacetamide",
    "deacetylase (AADAC) gene polymorphism and HIV infection affect the",
    "exposure of rifapentine.",
    "Antimicrob Agents Chemother 63(4):e01964-18.",
    "doi:10.1128/AAC.01964-18.",
    sep = " "
  )
  vignette <- "Francis_2019_rifapentine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "Rifapentine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "Rifapentine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass, the size descriptor used for allometric scaling of rifapentine clearance.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 46 kg per Francis 2019 Table 2 footnote b ('The typical",
        "values of clearance and volume of distribution reported for a patient with",
        "a body weight of 56 kg and FFM of 46 kg') and the Results sentence 'In a",
        "typical patient (FFM, 46 kg; weight, 56 kg), the values of clearance and",
        "volume of distribution were 1.33 liters/h and 25 liters'. The Abstract",
        "instead reports 45 kg, which is also the Table 1 overall cohort median;",
        "the 46 kg reading is used here because it appears in the parameter table's",
        "own footnote and in the Results body. The discrepancy moves typical",
        "clearance by (46/45)^0.75 = 1.7 percent - see the vignette Assumptions and",
        "deviations. The exponent is fixed at 0.75 per Methods 'Allometric scaling",
        "was applied to clearance (CL) and the volume of distribution (V) to adjust",
        "for the effect of body size, as described by Anderson and Holford'. FFM was",
        "selected over total body weight for clearance on a 23-point OFV improvement.",
        "The paper does not state which FFM equation produced the data column;",
        "Janmahasatian et al. (Clin Pharmacokinet 2005;44:1051-1065) is the usual",
        "choice in this literature, and the cohort's Table 1 median weight of 56 kg",
        "against a median FFM of 45 kg is the only in-paper calibration point."
      ),
      source_name = "FFM"
    ),
    WT = list(
      description = "Total body weight, the size descriptor used for allometric scaling of rifapentine apparent volume of distribution.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 56 kg, the Table 1 overall cohort median, per Francis 2019",
        "Table 2 footnote b. Exponent fixed at 1 per the Anderson and Holford",
        "allometric convention cited in Methods. Total body weight beat fat-free",
        "mass as the size descriptor for volume (delta-OFV 20, P < 0.001), the",
        "opposite of the clearance result. Observed range 38 to 94 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    HIV_POS = list(
      description = "HIV-positive comorbidity indicator (1 = HIV-infected, 0 = HIV-uninfected).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (HIV-negative)",
      notes = paste(
        "59 of 326 patients (18.1 percent) were HIV-positive (Francis 2019 Table 1).",
        "Multiplicative shift on bioavailability, F = 1 * (1 + e_hiv_pos_fdepot *",
        "HIV_POS); Table 2 row 'Effect of HIV+ on F (%)' = -21.9 (95% CI -33.2 to",
        "-6.64), i.e. HIV-positive patients have 21.9 percent lower bioavailability.",
        "The paper notes the data were not sufficient to identify interactions with",
        "the individual concomitant antiretroviral drugs."
      ),
      source_name = "HIV"
    ),
    SNP_AADAC_RS1803155_HOM = list(
      description = "AADAC rs1803155 (G>A) homozygous-variant indicator (1 = AA genotype, 0 = GA or GG genotype).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (GA or GG genotype)",
      notes = paste(
        "Arylacetamide deacetylase is the enzyme that deacetylates rifapentine to",
        "its primary metabolite 25-desacetyl rifapentine. 106 of 162 genotyped",
        "patients (65.4 percent) were AA; GA 53 (32.7 percent) and GG 3 (1.85",
        "percent), variant A allele frequency 0.82 (Francis 2019 Table 3).",
        "Multiplicative shift on clearance, CL = ... * (1 +",
        "e_snp_aadac_rs1803155_hom_cl * SNP_AADAC_RS1803155_HOM); Table 2 row",
        "'AADAC rs1803155 (AA) effect on CL (%)' = -10.4 (95% CI -17.3 to -3.53).",
        "GA and GG were first estimated as separate groups but their effects were",
        "similar and neither was statistically significant, so they were pooled into",
        "the reference category (Results paragraph 3); only 3 patients were GG.",
        "Genotype was missing for 164 of 326 patients and was imputed in the source",
        "analysis by mixture modelling (Methods, MIX method); a downstream user must",
        "supply an observed or imputed 0/1 value per subject."
      ),
      source_name = "AADAC rs1803155"
    ),
    DOSE_HIGH = list(
      description = "Highest-dose-cohort indicator (1 = RIFAQUIN 1,200 mg once-weekly arm, 0 = any other arm).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (RIFAQUIN 900 mg twice weekly, Daily RPE 450 mg daily, or Daily RPE 600 mg daily)",
      notes = paste(
        "Per-paper threshold: DOSE_HIGH = 1 for the 125 patients in the RIFAQUIN",
        "1,200-mg-once-weekly continuation-phase arm (Francis 2019 Table 1), the top",
        "of the 450 / 600 / 900 / 1,200 mg dose range in the pooled analysis.",
        "Multiplicative shift on clearance, CL = ... * (1 + e_dose_high_cl *",
        "DOSE_HIGH); Table 2 row 'Effect of group on 1,200-mg dose in RIFAQUIN study",
        "on CL (%)' = -13.2 (95% CI -22.8 to -4.36). The Discussion attributes the",
        "reduction to less pronounced autoinduction under once-weekly rather than",
        "more frequent dosing: 'The reduced dosing frequency in this group may have",
        "led to reduced autoinduction and, thus, increased exposure.' The indicator",
        "therefore encodes a dose-and-frequency arm, not a dose level alone; it is",
        "confounded with dosing interval by design and must not be read as a",
        "concentration-dependent saturation of clearance."
      ),
      source_name = "1,200-mg dose group in RIFAQUIN"
    ),
    STUDY_DAILY_RPE = list(
      description = "Daily RPE study cohort indicator (1 = Daily RPE trial, 0 = RIFAQUIN trial).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (RIFAQUIN trial)",
      notes = paste(
        "85 of 326 patients (44 at 450 mg and 41 at 600 mg daily) came from the Daily",
        "RPE two-stage activity-safety study of daily rifapentine",
        "(ClinicalTrials.gov NCT00814671); the remaining 241 came from the phase III",
        "RIFAQUIN trial (ISRCTN44153044). Multiplicative shift on bioavailability,",
        "F = 1 * (1 + e_study_daily_rpe_fdepot * STUDY_DAILY_RPE); Table 2 row",
        "'Effect of Daily RPE study on F (%)' = -23.3 (95% CI -35.6 to -9.25).",
        "The Discussion attributes the difference to food: RIFAQUIN gave a",
        "standardised light meal (two hard-boiled eggs with bread) 15 min before",
        "every dose, whereas Daily RPE advised patients to take the dose with food",
        "but neither standardised nor recorded the meal. The indicator is therefore",
        "a proxy for concomitant food intake, not an intrinsic trial effect."
      ),
      source_name = "study (Daily RPE vs RIFAQUIN)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 326L,
    n_studies = 2L,
    age_range = "18-80 years",
    age_median = "32 years",
    weight_range = "38-94 kg",
    weight_median = "56 kg",
    ffm_range = "27-62 kg",
    ffm_median = "45 kg",
    sex_female_pct = 33.1,
    race_ethnicity = "not reported by race; all participants were enrolled in southern Africa and the 162 genotyped patients were all from South African sites.",
    disease_state = "Drug-susceptible pulmonary tuberculosis; 59 of 326 (18.1 percent) HIV-coinfected.",
    dose_range = paste(
      "RIFAQUIN continuation phase: 1,200 mg rifapentine once weekly with 400 mg",
      "moxifloxacin for 4 months (n = 125), or 900 mg rifapentine twice weekly with",
      "400 mg moxifloxacin twice weekly for 2 months (n = 116); doses taken with",
      "240 ml water 15 min after a standardised light meal of two hard-boiled eggs",
      "with bread. Daily RPE intensive phase: 450 mg (n = 44) or 600 mg (n = 41)",
      "rifapentine daily replacing 600 mg rifampin; patients advised to dose with",
      "food, but no standardised meal was provided and food intake was not recorded."
    ),
    regions = "Western Cape and Gauteng, South Africa; Harare, Zimbabwe.",
    sampling_design = paste(
      "1,144 rifapentine plasma concentrations from 326 patients entered the",
      "analysis (1,151 measured, 7 below the LLOQ of 0.156 mg/L omitted).",
      "RIFAQUIN sampling took place during the 4th month of treatment: rich",
      "(predose and 1, 2, 3, 5, 7, 10, 12, 26 and 50 h post dose) or sparse (about",
      "2, 5 and 24 or 48 h post dose). Daily RPE sampling took place about 1 month",
      "after starting therapy: intensive (predose and 0.75, 1.5, 3.5, 5, 12 and 24 h",
      "post dose) or sparse (0.5-2 h and 5-8 h post dose). Pharmacogenetic data were",
      "available for 162 of 326 patients (49.7 percent)."
    ),
    assay = paste(
      "Validated LC-MS/MS assay (Division of Clinical Pharmacology, University of",
      "Cape Town) with protein-precipitation extraction and rifaximin as internal",
      "standard; quadratic calibration weighted by 1/concentration over 0.156 to",
      "40.0 mg/L; interbatch accuracy 103.9, 102.8 and 97.5 percent at the low,",
      "medium and high quality-control levels; LLOQ 0.156 mg/L."
    ),
    notes = paste(
      "Baseline demographics are Francis 2019 Table 1. Genotype and allele",
      "frequencies for the five screened SNPs (SLCO1B1 rs2306283 and rs4149032,",
      "NR1I2 rs2472677 and rs1523130, AADAC rs1803155) are Table 3; only AADAC",
      "rs1803155 was retained in the final model. Estimation used NONMEM 7.4.2 with",
      "FOCE INTER; parameter precision came from sampling importance resampling",
      "(SIR, n = 1,000), which is the source of the 95 percent confidence intervals",
      "quoted in the parameter comments below."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters - Francis 2019 Table 2, 'Final parameter estimates
    # for rifapentine population pharmacokinetic model'. Typical values are for
    # the reference patient of Table 2 footnote b: body weight 56 kg, FFM 46 kg.
    # Time is in hours, doses in mg, concentrations in mg/L throughout, matching
    # the source.
    # =========================================================================
    lcl <- log(1.33)
    label("Rifapentine apparent oral clearance CL/F at FFM = 46 kg (L/h)")
    # Francis 2019 Table 2 row 'CL (liters/h)' = 1.33 (95% CI 1.14, 1.54).

    lvc <- log(25)
    label("Rifapentine apparent central volume V/F at WT = 56 kg (L)")
    # Francis 2019 Table 2 row 'V (liters)' = 25 (95% CI 21.9, 28.4).

    lka <- log(0.814)
    label("Rifapentine first-order absorption rate constant ka (1/h)")
    # Francis 2019 Table 2 row 'ka (h-1)' = 0.814 (95% CI 0.568, 1.26). This is
    # the rate out of the depot that terminates the transit chain; it is
    # estimated separately from the transit rate ktr = (NN + 1) / MTT.

    lmtt <- log(1.47)
    label("Rifapentine absorption mean transit time MTT (h)")
    # Francis 2019 Table 2 row 'MTT (h)' = 1.47 (95% CI 1.20, 1.78).

    lnn <- log(10.2)
    label("Rifapentine transit-chain shape parameter NN (unitless)")
    # Francis 2019 Table 2 row 'NN' = 10.2 (95% CI 6.70, 14.0), the number of
    # hypothetical transit compartments. Non-integer; the rxode2 transit()
    # closed form evaluates the gamma density at non-integer n.

    lfdepot <- fixed(log(1))
    label("Rifapentine reference oral bioavailability F (unitless anchor of 1)")
    # Francis 2019 Table 2 row 'F' = '1 (fixed)'. The reference is an
    # HIV-negative RIFAQUIN participant; the HIV and study effects below are
    # applied multiplicatively to this anchor.

    # =========================================================================
    # Allometric size scaling - Francis 2019 Methods, 'Allometric scaling was
    # applied to clearance (CL) and the volume of distribution (V) to adjust for
    # the effect of body size, as described by Anderson and Holford'. The
    # Anderson-Holford convention fixes the exponents at 0.75 for flow terms and
    # 1 for volume terms; neither exponent is reported with an estimate or
    # confidence interval in Table 2, confirming they were not estimated.
    # =========================================================================
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent on CL for fat-free mass, (FFM/46)^0.75 (unitless)")
    # Francis 2019 Methods, Pharmacokinetic analysis, final paragraph. FFM was
    # the best size descriptor for clearance (delta-OFV 93 vs no scaling, and 23
    # points better than total body weight).

    e_wt_vc <- fixed(1)
    label("Allometric exponent on V for total body weight, (WT/56)^1 (unitless)")
    # Francis 2019 Methods, Pharmacokinetic analysis, final paragraph. Total body
    # weight was the best size descriptor for volume (delta-OFV 20, P < 0.001).

    # =========================================================================
    # Covariate effects - Francis 2019 Table 2, lower block. All four are
    # reported as percent changes and are applied in the NONMEM fractional form
    # theta_typical * (1 + coefficient * indicator), so a reported -21.9 percent
    # becomes a coefficient of -0.219.
    # =========================================================================
    e_snp_aadac_rs1803155_hom_cl <- -0.104
    label("Fractional change in CL for AADAC rs1803155 AA vs GA/GG (unitless)")
    # Francis 2019 Table 2 row 'AADAC rs1803155 (AA) effect on CL (%)' = -10.4
    # (95% CI -17.3, -3.53); Results 'patients homozygous for the AADAC rs1803155
    # AA polymorphism were found to have a 10.4% lower clearance of rifapentine
    # than subjects that were rs1803155 GG or GA (delta-OFV 6.2; P = 0.013)'.

    e_dose_high_cl <- -0.132
    label("Fractional change in CL for the RIFAQUIN 1,200 mg once-weekly arm (unitless)")
    # Francis 2019 Table 2 row 'Effect of group on 1,200-mg dose in RIFAQUIN study
    # on CL (%)' = -13.2 (95% CI -22.8, -4.36); Results 'clearance reduced by
    # 13.2% compared to the clearance for the patients in the other dose groups
    # (delta-OFV 17; P < 0.001)'.

    e_hiv_pos_fdepot <- -0.219
    label("Fractional change in bioavailability for HIV-positive patients (unitless)")
    # Francis 2019 Table 2 row 'Effect of HIV + on F (%)' = -21.9 (95% CI -33.2,
    # -6.64); Results 'Patients infected with HIV were found to have a 21.9%
    # lower bioavailability (delta-OFV 42; P < 0.001)'.

    e_study_daily_rpe_fdepot <- -0.233
    label("Fractional change in bioavailability for Daily RPE vs RIFAQUIN (unitless)")
    # Francis 2019 Table 2 row 'Effect of Daily RPE study on F (%)' = -23.3
    # (95% CI -35.6, -9.25); Results 'the bioavailability of rifapentine in the
    # Daily RPE study was 23.3% lower than that in the RIFAQUIN study
    # (delta-OFV 59; P < 0.001)'.

    # =========================================================================
    # Random effects - Francis 2019 Table 2 'Variability' columns, reported as
    # percent CV with 95 percent confidence intervals. Methods: 'A lognormal
    # distribution was assumed for IIV and IOV', so the internal-scale variance
    # is omega^2 = log(1 + CV^2).
    #
    # Table 2 labels the CL and V terms IIV and the ka, MTT and NN terms IOV.
    # The source never states how many pharmacokinetic occasions each patient
    # contributed, and nlmixr2lib has no standard OCC encoding for an
    # unenumerated occasion design, so the three IOV terms are carried here as
    # single-draw subject-level random effects. That reproduces the published
    # variance magnitudes exactly for a one-occasion simulation, which is what
    # both trials' sampling schedules describe (RIFAQUIN sampled during month 4,
    # Daily RPE at about month 1). A user simulating more than one occasion per
    # subject must redraw etalka, etalmtt and etalnn per occasion while holding
    # etalcl and etalvc fixed. See the vignette Assumptions and deviations.
    # =========================================================================
    # Table 2 'CL (liters/h)' variability 23.0% CV (IIV) (95% CI 17.7, 28.6);
    # omega^2 = log(1 + 0.230^2) = 0.0515483.
    etalcl ~ 0.0515483

    # Table 2 'V (liters)' variability 12.8% CV (IIV) (95% CI 8.8, 17.4);
    # omega^2 = log(1 + 0.128^2) = 0.0162512.
    etalvc ~ 0.0162512

    # Table 2 'ka (h-1)' variability 48.9% CV (IOV) (95% CI 36.4, 59.8);
    # omega^2 = log(1 + 0.489^2) = 0.2144023.
    etalka ~ 0.2144023

    # Table 2 'MTT (h)' variability 37.4% CV (IOV) (95% CI 28.3, 48.6);
    # omega^2 = log(1 + 0.374^2) = 0.1309195.
    etalmtt ~ 0.1309195

    # Table 2 'NN' variability 20.3% CV (IOV) (95% CI 14.9, 26.4);
    # omega^2 = log(1 + 0.203^2) = 0.0403825.
    etalnn ~ 0.0403825

    # =========================================================================
    # Residual unexplained variability - Francis 2019 Table 2, combined additive
    # plus proportional model (Methods: 'a combined additive and proportional
    # model for the residual unexplained variability (RUV) was evaluated').
    # =========================================================================
    propSd <- 0.0956
    label("Proportional residual standard deviation (fraction)")
    # Francis 2019 Table 2 row 'Proportional residual error (%)' = 9.56
    # (95% CI 7.09, 13.2).

    addSd <- 0.247
    label("Additive residual standard deviation (mg/L)")
    # Francis 2019 Table 2 row 'Additive residual error (mg/liter)' = 0.247
    # (95% CI 0.143, 0.401).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual disposition parameters.
    #    Allometric size scaling on both terms, then the two multiplicative
    #    clearance covariate effects in the NONMEM fractional form
    #    theta * (1 + coefficient * indicator).
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (FFM / 46)^e_ffm_cl *
      (1 + e_snp_aadac_rs1803155_hom_cl * SNP_AADAC_RS1803155_HOM) *
      (1 + e_dose_high_cl * DOSE_HIGH)
    vc <- exp(lvc + etalvc) * (WT / 56)^e_wt_vc

    # -----------------------------------------------------------------------
    # 2. Individual absorption parameters. The Savic transit chain is
    #    parameterised by MTT and NN with transit rate ktr = (NN + 1) / MTT,
    #    which is what rxode2's transit() computes internally; the final step
    #    out of the depot runs at the separately estimated ka. Francis 2019
    #    cites Savic et al. as reference 31 for the transit model, and the
    #    sibling rifapentine model from the same group and the same
    #    (NN ~ 10, MTT ~ 1.5 h, separate ka) parameterisation is
    #    Zvada_2010_rifapentine.R.
    # -----------------------------------------------------------------------
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    nn <- exp(lnn + etalnn)

    # -----------------------------------------------------------------------
    # 3. Bioavailability. The reference (F = 1) is an HIV-negative RIFAQUIN
    #    participant; the two effects are fractional shifts off that anchor.
    # -----------------------------------------------------------------------
    fdepot <- exp(lfdepot) *
      (1 + e_hiv_pos_fdepot * HIV_POS) *
      (1 + e_study_daily_rpe_fdepot * STUDY_DAILY_RPE)

    # -----------------------------------------------------------------------
    # 4. ODE system. transit() supplies the gamma-density input rate into the
    #    depot and carries the bioavailability, so f(depot) <- 0 suppresses the
    #    ordinary dose bolus that would otherwise enter the depot a second time.
    # -----------------------------------------------------------------------
    kel <- cl / vc

    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- 0

    # -----------------------------------------------------------------------
    # 5. Observation. Doses are in mg and volumes in L, so central / vc is
    #    already in mg/L, the unit of the source's reported concentrations.
    # -----------------------------------------------------------------------
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
