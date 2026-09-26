# One-compartment parent-plus-metabolite population PK of oral sunitinib and
# its active metabolite N-desethyl sunitinib (SU12662) in Chinese patients with
# renal-cell carcinoma, with body weight and ABCB1 rs2032582 genotype on the
# apparent clearance of SU12662 (Zhang 2018, Oncotarget 9(18):14109-14123;
# doi:10.18632/oncotarget.23881).

Zhang_2018_sunitinib <- function() {
  description <- paste(
    "One-compartment parent-plus-metabolite population PK model for oral",
    "sunitinib and its active metabolite N-desethyl sunitinib (SU12662) in 53",
    "Chinese patients with renal-cell carcinoma dosed on the 4-weeks-on /",
    "2-weeks-off schedule (Zhang 2018). Sunitinib is absorbed first order into",
    "a single central compartment and eliminated with apparent oral clearance",
    "Clp/F; that elimination flux is the sole input to a single SU12662",
    "compartment, which is eliminated with apparent clearance Clm/F. Because",
    "only oral data were available, the fraction of sunitinib converted to",
    "SU12662 is not identifiable and is folded into the apparent metabolite",
    "clearance and volume, so the metabolite parameters are apparent values",
    "conditioned on complete conversion. Sunitinib and SU12662 were fitted",
    "simultaneously in Phoenix NLME. Body weight (power 0.538, reference 68.3",
    "kg) and the ABCB1 rs2032582 genotype entered as a six-level linear",
    "proportional effect (AT reference) on Clm/F only; the paper reports that",
    "no covariate (weight, age, sex or genotype) was retained on the parent",
    "Clp/F or V/F. Inter-individual variability was estimated (eta shrinkages",
    "are tabulated) but no omega values were published, so every eta is",
    "carried at fixed(0); residual error is combined proportional plus",
    "additive, reported separately for sunitinib and SU12662.",
    sep = " "
  )
  reference <- paste(
    "Zhang Y, Mai H, Guo G, Bi G, Hao G, Li Y, Wang X, Cheng L, Wang J,",
    "Dong R, Liu Z, Chen L, Qu H (2018). Association analysis of SNPs present",
    "in plasma with adverse events and population pharmacokinetics in Chinese",
    "sunitinib treated patients with renal cell carcinoma.",
    "Oncotarget 9(18):14109-14123. doi:10.18632/oncotarget.23881.",
    sep = " "
  )
  vignette <- "Zhang_2018_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at baseline. Enters only the apparent clearance of the",
        "metabolite SU12662, as a power function normalised to the cohort",
        "average body weight of 68.3 kg stated in Methods 'Model",
        "development' (equation 4); the exponent 0.538 was estimated, not",
        "fixed at an allometric value (Table 6 dClmdBW). Zhang 2018 does not",
        "tabulate the weight distribution; Table 1 reports body surface area",
        "(median 1.86 m2, quartiles 1.74-1.94) instead."
      ),
      source_name = "weight (kg)"
    ),
    SNP_ABCB1_RS2032582_TG = list(
      description = "ABCB1 rs2032582 TG genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 in combination with the other four rs2032582 indicators = 0, i.e. the AT genotype group (Z3 = 0)",
      notes = paste(
        "1 = subject's reported rs2032582 genotype string is 'TG'; 0",
        "otherwise. One of the five indicators that encode Zhang 2018's",
        "six-level Z3 covariate (Z3 = 0 AT, 1 TG, 2 GG, 3 TT, 4 AG, 5 GT;",
        "final equation, page 14114). AT is the reference. Note that the",
        "source distinguishes 'TG' (Z3 = 1) from 'GT' (Z3 = 5) and assigns",
        "them materially different coefficients (0.314 vs 0.0456) even though",
        "they denote the same unordered heterozygous genotype; the allele",
        "order in the source's genotype strings is therefore load-bearing and",
        "must be reproduced as reported."
      ),
      source_name = "Z3 = 1 (TG)"
    ),
    SNP_ABCB1_RS2032582_GG = list(
      description = "ABCB1 rs2032582 GG genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 in combination with the other four rs2032582 indicators = 0, i.e. the AT genotype group (Z3 = 0)",
      notes = paste(
        "1 = subject's reported rs2032582 genotype string is 'GG'; 0",
        "otherwise. Z3 = 2 in Zhang 2018's final equation. Table 2 reports 10",
        "of 53 subjects (18.87%) as GG."
      ),
      source_name = "Z3 = 2 (GG)"
    ),
    SNP_ABCB1_RS2032582_TT = list(
      description = "ABCB1 rs2032582 TT genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 in combination with the other four rs2032582 indicators = 0, i.e. the AT genotype group (Z3 = 0)",
      notes = paste(
        "1 = subject's reported rs2032582 genotype string is 'TT'; 0",
        "otherwise. Z3 = 3 in Zhang 2018's final equation. Table 2 pools TT",
        "with AA and TA into a single 17-subject (32.08%) row, so the TT-only",
        "count is not recoverable from the publication."
      ),
      source_name = "Z3 = 3 (TT)"
    ),
    SNP_ABCB1_RS2032582_AG = list(
      description = "ABCB1 rs2032582 AG genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 in combination with the other four rs2032582 indicators = 0, i.e. the AT genotype group (Z3 = 0)",
      notes = paste(
        "1 = subject's reported rs2032582 genotype string is 'AG'; 0",
        "otherwise. Z3 = 4 in Zhang 2018's final equation. rs2032582 is",
        "tri-allelic (G / T / A) in this cohort, which is why genotypes",
        "containing an A allele appear at all; Table 2 reports an A allele",
        "frequency of 0.12."
      ),
      source_name = "Z3 = 4 (AG)"
    ),
    SNP_ABCB1_RS2032582_GT = list(
      description = "ABCB1 rs2032582 GT genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 in combination with the other four rs2032582 indicators = 0, i.e. the AT genotype group (Z3 = 0)",
      notes = paste(
        "1 = subject's reported rs2032582 genotype string is 'GT'; 0",
        "otherwise. Z3 = 5 in Zhang 2018's final equation. Distinct from the",
        "'TG' stratum (Z3 = 1) in the source's coding; see the",
        "SNP_ABCB1_RS2032582_TG notes for why the allele order is reproduced",
        "rather than collapsed."
      ),
      source_name = "Z3 = 5 (GT)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_su12662 = list(
      analyte = "N-desethyl sunitinib (SU12662)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 53L,
    n_observations = 127L,
    n_studies = 1L,
    age_range = "19-71 years (Results, 'Baseline characteristic of Chinese patients with RCC')",
    age_median = "54 years (Results text); Table 1 reports 51.65 years with quartiles 45-61",
    weight_range = "not reported; Methods states a cohort average body weight of 68.3 kg",
    weight_median = "68.3 kg (cohort average used as the covariate reference, Methods equation 4)",
    bsa_median = "1.86 m2 (quartiles 1.74-1.94; Table 1)",
    sex_female_pct = 24.53,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "Histologically confirmed renal-cell carcinoma; 75.47% had a prior",
      "nephrectomy, ECOG performance status 0 in 37 and 1 in 16 subjects,",
      "45.28% with one metastatic site, 15.09% with two and 20.75% with three",
      "or more (Table 1)."
    ),
    dose_range = paste(
      "Oral sunitinib 50 mg once daily for 4 weeks followed by 2 weeks off",
      "(schedule 4/2), with protocol-permitted reductions for toxicity. Daily",
      "dose over the first four cycles: 50 mg in 64.15%, 37.5 mg in 15.09%",
      "and 25 mg in 20.75% of subjects (Table 1); 35.85% had a dose reduction",
      "after cycle 1, 2 or 3."
    ),
    regions = "China (Academy of Military Medical Sciences Affiliated Hospital, Beijing)",
    notes = paste(
      "127 plasma samples from 53 subjects enrolled March 2014 to January",
      "2016. PK samples were drawn on day 15 +/- 1 of treatment, after the",
      "steady state the authors considered established at 14 days; sampling",
      "is therefore trough-dominated and sparse, which the Discussion",
      "acknowledges as the reason the model differs from richer published",
      "sunitinib models. Sunitinib and SU12662 plasma concentrations were",
      "measured by a validated LC-MS/MS method (intra- and inter-day",
      "precision < 1.64% and accuracy within +/-5.69% for sunitinib; < 12.4%",
      "and +/-12.74% for SU12662). Estimation used Phoenix NLME 1.4 with",
      "FOCE; the final model was evaluated by 1000-run nonparametric",
      "bootstrap (Table 7) and visual predictive check."
    )
  )

  ini({
    # Structural parameters -- Zhang 2018 Table 6, 'Final model' block.
    # Volumes are reported in mL and clearances in mL/h; those units are
    # preserved here and converted to ng/mL at the observation lines.

    lka <- log(0.0117)
    label("Sunitinib first-order absorption rate constant ka (1/h)") # Table 6 final model: tvKa = 0.0117 1/h (RSE 11.4%; bootstrap 95% CI 0.00892-0.0140)

    lvc <- log(99438)
    label("Sunitinib apparent central volume of distribution V/F (mL)") # Table 6 final model: tvV = 99438 mL (RSE 43.4%; bootstrap 95% CI 4459-191880)

    lcl <- log(24576)
    label("Sunitinib apparent oral clearance Clp/F (mL/h)") # Table 6 final model: tvClp = 24576 mL/h (RSE 6.21%; bootstrap 95% CI 21769-27036)

    lvc_su12662 <- log(916641)
    label("SU12662 apparent volume of distribution V2/F (mL)") # Table 6 final model: tvV2 = 916641 mL (RSE 50.7%; bootstrap 95% CI 415037-1569500)

    lcl_su12662 <- log(53614)
    label("SU12662 apparent clearance Clm/F at 68.3 kg and the AT reference genotype (mL/h)") # Table 6 final model and the 'final equation' intercept: tvClm = 53614 mL/h (RSE 9.9%; bootstrap 95% CI 39962-64261)

    # Covariate effects on Clm/F only. Zhang 2018 Results: 'No covariates
    # (BW, age, sex and SNPs) had a relationship with Vd/F and CL/F
    # parameters', so the parent disposition carries no covariates.
    e_wt_cl_su12662 <- 0.538
    label("Power of body weight / 68.3 kg on SU12662 apparent clearance (unitless)") # Table 6 dClmdBW = 0.538 (RSE 63.7%; bootstrap 95% CI 0.00691-1.41); reference weight 68.3 kg from Methods equation 4

    # ABCB1 rs2032582 genotype effects, entered as the source's linear
    # proportional form (1 - d * indicator) with AT (Z3 = 0) as reference.
    e_snp_abcb1_rs2032582_tg_cl_su12662 <- 0.314
    label("Proportional reduction in SU12662 apparent clearance for the ABCB1 rs2032582 TG genotype (fraction)") # Table 6 dClmdZ31 = 0.314 (RSE 60.5%); final equation term (1 - 0.314 * (Z3 = 1))

    e_snp_abcb1_rs2032582_gg_cl_su12662 <- 0.269
    label("Proportional reduction in SU12662 apparent clearance for the ABCB1 rs2032582 GG genotype (fraction)") # Table 6 dClmdZ32 = 0.269 (RSE 68.4%); final equation term (1 - 0.269 * (Z3 = 2))

    e_snp_abcb1_rs2032582_tt_cl_su12662 <- 0.308
    label("Proportional reduction in SU12662 apparent clearance for the ABCB1 rs2032582 TT genotype (fraction)") # Table 6 dClmdZ33 = 0.308 (RSE 62.1%); final equation term (1 - 0.308 * (Z3 = 3))

    e_snp_abcb1_rs2032582_ag_cl_su12662 <- 0.0368
    label("Proportional reduction in SU12662 apparent clearance for the ABCB1 rs2032582 AG genotype (fraction)") # Table 6 dClmdZ34 = 0.0368 (RSE 19.8%); final equation term (1 - 0.0368 * (Z3 = 4))

    e_snp_abcb1_rs2032582_gt_cl_su12662 <- 0.0456
    label("Proportional reduction in SU12662 apparent clearance for the ABCB1 rs2032582 GT genotype (fraction)") # Table 6 dClmdZ35 = 0.0456 (RSE 7.58%); final equation term (1 - 0.0456 * (Z3 = 5))

    # Inter-individual variability. Zhang 2018 Table 6 reports an eta
    # shrinkage for each of the five structural parameters (tvKa 0.386,
    # tvV 0.773, tvV2 0.975, tvClp 0.317, tvClm 0.163), which proves the
    # final model carried an eta on each, but no omega / variance / CV is
    # published anywhere in the paper -- the 'CV%' column is the relative
    # standard error of the estimate, as the bootstrap confidence intervals
    # in Table 7 confirm. The variances are therefore carried at fixed(0)
    # rather than invented; see the vignette 'Assumptions and deviations'.
    etalka ~ fixed(0) # Table 6 final model reports only eta shrinkage 0.386 for tvKa; no omega published
    etalvc ~ fixed(0) # Table 6 final model reports only eta shrinkage 0.773 for tvV; no omega published
    etalcl ~ fixed(0) # Table 6 final model reports only eta shrinkage 0.317 for tvClp; no omega published
    etalvc_su12662 ~ fixed(0) # Table 6 final model reports only eta shrinkage 0.975 for tvV2; no omega published
    etalcl_su12662 ~ fixed(0) # Table 6 final model reports only eta shrinkage 0.163 for tvClm; no omega published

    # Residual error -- combined proportional plus additive, separately for
    # sunitinib and SU12662 (Methods equation 3). Phoenix names the
    # proportional standard deviation of observation k '<C>MultStdev' and the
    # additive standard deviation of its epsilon 'stdev<k-1>', so tvCMultStdev
    # pairs with stdev0 (sunitinib, C) and tvC2MultStdev with stdev1 (SU12662,
    # C2). See the vignette for why the two components are read as
    # independent standard deviations rather than as a Phoenix mixed-ratio.
    propSd <- 0.31
    label("Proportional residual error for sunitinib (fraction)") # Table 6 final model: tvCMultStdev = 0.31 (RSE 13%; bootstrap 95% CI 0.169-0.338)

    addSd <- 0.0751
    label("Additive residual error for sunitinib (ng/mL)") # Table 6 final model: stdev0 = 0.0751 (RSE 12.9%; bootstrap 95% CI 0.0000526-6.33 -- very poorly identified)

    propSd_su12662 <- 0.242
    label("Proportional residual error for SU12662 (fraction)") # Table 6 final model: tvC2MultStdev = 0.242 (RSE 14.5%; bootstrap 95% CI 0.194-0.315)

    addSd_su12662 <- 1.23
    label("Additive residual error for SU12662 (ng/mL)") # Table 6 final model: stdev1 = 1.23 (RSE 37.5%; bootstrap 95% CI 0.00464-1.61)
  })

  model({
    # 1. Individual parameters. The parent carries no covariates (Results:
    #    weight, age, sex and genotype were all rejected on V/F and CL/F).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. SU12662 apparent clearance -- the paper's 'final equation'. Body
    #    weight enters as a power term normalised to 68.3 kg; the ABCB1
    #    rs2032582 genotype enters as the product of five linear proportional
    #    terms, of which at most one indicator is 1 for any subject, so the
    #    product collapses to a single multiplier. All five indicators are 0
    #    for the AT reference genotype (Z3 = 0).
    vc_su12662 <- exp(lvc_su12662 + etalvc_su12662)
    cl_su12662 <- exp(lcl_su12662 + etalcl_su12662) *
      (WT / 68.3)^e_wt_cl_su12662 *
      (1 - e_snp_abcb1_rs2032582_tg_cl_su12662 * SNP_ABCB1_RS2032582_TG) *
      (1 - e_snp_abcb1_rs2032582_gg_cl_su12662 * SNP_ABCB1_RS2032582_GG) *
      (1 - e_snp_abcb1_rs2032582_tt_cl_su12662 * SNP_ABCB1_RS2032582_TT) *
      (1 - e_snp_abcb1_rs2032582_ag_cl_su12662 * SNP_ABCB1_RS2032582_AG) *
      (1 - e_snp_abcb1_rs2032582_gt_cl_su12662 * SNP_ABCB1_RS2032582_GT)

    # 3. Micro-constants.
    kel <- cl / vc
    kel_su12662 <- cl_su12662 / vc_su12662

    # 4. ODE system. Sunitinib is absorbed first order from the depot and
    #    eliminated at kel; that whole elimination flux is the input to the
    #    SU12662 compartment. Because only oral data were fitted, neither the
    #    fraction of sunitinib converted to SU12662 nor its bioavailability is
    #    identifiable, so both are folded into the apparent Clm/F and V2/F and
    #    the conversion is written as complete. The resulting steady-state
    #    identity Cm,avg = Dose / (Clm/F * tau) is exactly what the paper's
    #    observed SU12662 trough levels support (see the vignette).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_su12662) <- kel * central - kel_su12662 * central_su12662

    # 5. Observations. Amounts are in mg and volumes in mL, so amount/volume
    #    is mg/mL; multiplying by 1e6 gives ng/mL, the units of every
    #    concentration reported in Zhang 2018 (Table 4, Figure 1).
    Cc <- 1e6 * central / vc
    Cc_su12662 <- 1e6 * central_su12662 / vc_su12662

    Cc ~ prop(propSd) + add(addSd)
    Cc_su12662 ~ prop(propSd_su12662) + add(addSd_su12662)
  })
}
