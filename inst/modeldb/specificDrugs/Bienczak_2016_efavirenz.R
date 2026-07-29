Bienczak_2016_efavirenz <- function() {
  description <- paste(
    "Two-compartment population PK model for oral efavirenz in African",
    "children (Bienczak 2016), with Savic 2007 transit-compartment",
    "absorption (NN = 25 fixed transit compartments and a separate",
    "first-order absorption step ka from the depot to the central",
    "compartment), oral bioavailability fixed to 1 (no intravenous data),",
    "Anderson-Holford allometric scaling of all clearance and volume",
    "parameters to a 15.4 kg reference child (exponents 0.75 on CL and Q,",
    "1.0 on Vc and Vp), and a composite CYP2B6 516G>T (rs3745274) | 983T>C",
    "(rs28399499) SNP-vector effect on apparent oral clearance that",
    "distinguishes six metabolic subgroups (516GG|983TT extensive",
    "metabolizer reference, 516GG|983TC and 516GT|983TT intermediate,",
    "516TT|983TT and 516GT|983TC slow, 516GG|983CC ultra-slow). Encoded as",
    "log-ratio multiplicative shifts on the 516GG|983TT EM reference so",
    "the single etalcl IIV applies uniformly on the log-CL scale across",
    "all six SNP-vector subgroups."
  )
  reference <- paste(
    "Bienczak A, Cook A, Wiesner L, Olagunju A, Mulenga V, Kityo C,",
    "Kekitiinwa A, Owen A, Walker AS, Gibb DM, McIlleron H, Burger D,",
    "Denti P (2016).",
    "The impact of genetic polymorphisms on the pharmacokinetics of",
    "efavirenz in African children.",
    "British Journal of Clinical Pharmacology 82(1):185-198.",
    "doi:10.1111/bcp.12934.",
    sep = " "
  )
  vignette <- "Bienczak_2016_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline / time-varying body weight. Drives the Anderson-Holford allometric scaling of apparent oral clearance CL/F and inter-compartmental clearance Q/F (exponent 0.75) and central and peripheral volume Vc/F and Vp/F (exponent 1.0), with the structural typical-value parameters reported in Bienczak 2016 Table 2 corresponding to the cohort median 15.4 kg (Table 2 footnote: 'All clearance and volume parameters scaled allometrically to median weight of 15.4 kg').",
      source_name        = "WT"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Bienczak 2016 Methods 'Genotyping' paragraph 1 -- the rs3745274 variant is jointly tested with rs28399499 (983T>C) to define the six observed composite SNP-vector subgroups on which the final CL/F effect was estimated. Cohort allele-genotype frequencies (n = 162 genotyped of 169 enrolled; Bienczak 2016 Table 1): GG 40%, GT 41%, TT 19%; minor-allele frequency 0.39.",
      source_name        = "CYP2B6 516G>T (rs3745274)"
    ),
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description        = "Count of CYP2B6 c.983T>C (rs28399499, p.I328T) C-alleles per subject (0/1/2). 0 = TT homozygous wild-type, 1 = TC heterozygous, 2 = CC homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Bienczak 2016 Methods 'Genotyping' paragraph 1 -- combined with rs3745274 to define the six observed CYP2B6 SNP-vector subgroups. Cohort allele-genotype frequencies (n = 162 genotyped of 169 enrolled; Bienczak 2016 Table 1): TT 86%, TC 14%, CC 1%; minor-allele frequency 0.07. The 983C variant (CYP2B6*18) is virtually absent in European-ancestry populations and reaches appreciable frequency only in sub-Saharan African cohorts; Bienczak 2016 is one of the first paediatric studies to quantify the 983CC homozygote (ultra-slow metabolizer) CL/F effect.",
      source_name        = "CYP2B6 983T>C (rs28399499)"
    )
  )

  covariatesDataExcluded <- list(
    STUDY_ARROW = list(
      description = "Trial-of-origin indicator (ARROW vs CHAPAS-3 reference)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Bienczak 2016 retained a significant trial-of-origin effect on the absorption rate constant ka (1.6-fold larger in ARROW, dOFV = 37.9, P < 0.001) and mean transit time MTT (1.4-fold longer in ARROW, dOFV = 21.4, P < 0.001), attributed by the paper to formulation differences between trials (CHAPAS-3 used only the double-scored 600 mg paediatric tablet; ARROW used a mix of 50, 100, and 200 mg capsules and half / whole 600 mg tablets). This packaged model encodes the CHAPAS-3 reference absorption parameters (MTT = 0.82 h, ka = 0.79 /h, Bienczak 2016 Table 2 column 'CHAPAS-3') and does NOT carry the trial-of-origin covariate, so the model is faithful to the CHAPAS-3 subset but slightly under-predicts the absorption rate (and slightly over-predicts the absorption transit time) for the ARROW subset. The ARROW typical values (MTT = 1.17 h, ka = 1.27 /h) are reported in Bienczak 2016 Table 2 column 'ARROW' for downstream users who need them. See vignette Assumptions and deviations for the rationale (avoiding a new STUDY_ARROW canonical covariate registration)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 169L,
    n_studies      = 2L,
    age_range      = "2.1-13.8 years (paediatric)",
    age_median     = "4.7 years",
    weight_range   = "7.8-30.0 kg",
    weight_median  = "15.5 kg (cohort) / 15.4 kg (allometric-scaling reference, Bienczak 2016 Table 2 footnote)",
    sex_female_pct = 52.7,
    race_ethnicity = "African (all 169 patients black African; recruited from Uganda and Zambia per Bienczak 2016 Table 1 footnote)",
    disease_state  = "HIV-1 infection on once-daily efavirenz-based combination antiretroviral therapy (paediatric). Companion nucleoside reverse-transcriptase inhibitors (NRTI backbone) and tuberculosis co-treatment status were tested as covariates on efavirenz PK but were not retained in the final model (Bienczak 2016 Results 'Population pharmacokinetics' paragraph 4: 'No other covariate ... was found to significantly improve the model fit').",
    dose_range     = "Once-daily oral efavirenz dosed by WHO weight-band guidelines. CHAPAS-3 used a new paediatric double-scored 600 mg tablet (provided by Cipla Pharmaceuticals) split to give 200, 300, 400 or 600 mg doses by weight band (Bienczak 2016 Table 4 dosing table). ARROW used 50, 100 and 200 mg capsules and half or whole 600 mg tablets per modified WHO 2006 paediatric recommendations. Drug was taken either in the morning or at night; on intensive PK days observed in the clinic, dosing was advised to be in the morning for 4-6 weeks prior to sampling.",
    regions        = "Uganda and Zambia (sub-Saharan Africa).",
    notes          = "Pooled cohort from two paediatric HIV trials: CHAPAS-3 (Children with HIV in Africa - Pharmacokinetics and Adherence/Acceptability of Simple antiretroviral regimens; 128 children with intensive + sparse PK) and ARROW (Anti-Retroviral Research for Watoto; 41 children with intensive PK). Final analysis dataset had 2086 efavirenz concentration measurements (611 intensive ARROW, 474 intensive CHAPAS-3, 1002 sparse CHAPAS-3) after exclusion of 22 samples (Bienczak 2016 Table 1). Genotypes were available for 162 of 169 children (95.9%); the remaining 7 were assigned via mixture-model imputation using cohort frequencies (Bienczak 2016 Results 'Demographic results and samples' paragraph 1). The mixture-model imputation is NOT encoded in this nlmixr2lib model -- this file assumes known genotype and the user supplies the two SNP allele-count covariate columns directly. Six observed CYP2B6 516G>T|983T>C SNP-vector subgroups (Bienczak 2016 Table 2 / Table 3): 516GG|983TT (n = 56, 33%), 516GT|983TT (n = 59, 35%), 516GG|983TC (n = 10, 6%), 516TT|983TT (n = 31, 18%), 516GT|983TC (n = 12, 7%), 516GG|983CC (n = 1, 1%). The remaining three combinations (516GT|983CC, 516TT|983TC, 516TT|983CC) were not observed."
  )

  ini({
    # =========================================================================
    # Structural disposition (Bienczak 2016 Table 2 'Final parameter
    # estimates'). All clearance and volume parameters reported at the
    # cohort median 15.4 kg (Table 2 footnote: 'All clearance and volume
    # parameters scaled allometrically to median weight of 15.4 kg.'). The
    # reference for CL/F here is the 516GG|983TT extensive-metabolizer
    # phenotype (CL/F = 6.94 L/h, Table 2 'CL 516GG|983TT' row), and the
    # other five observed SNP-vector subgroups enter as log-ratio
    # multiplicative shifts in model() below.
    # =========================================================================
    lcl_GG_TT <- log(6.94)
    label("CL/F at the 516GG|983TT extensive-metabolizer reference (L/h), allometrically scaled to a 15.4 kg child")
    # Bienczak 2016 Table 2 row 'CL 516GG|983TT = 6.94 (6.47-7.61)' (bootstrap median, 5th-95th)

    lvc <- log(64.1)
    label("Apparent central volume of distribution Vc/F (L), allometrically scaled to a 15.4 kg child")
    # Bienczak 2016 Table 2 row 'Vc = 64.1 (49.1-73.3)'

    lq <- log(17.1)
    label("Apparent inter-compartmental clearance Q/F (L/h), allometrically scaled to a 15.4 kg child")
    # Bienczak 2016 Table 2 row 'Q = 17.1 (14.1-20.9)'

    lvp <- log(92.2)
    label("Apparent peripheral volume of distribution Vp/F (L), allometrically scaled to a 15.4 kg child")
    # Bienczak 2016 Table 2 row 'Vp = 92.2 (80.1-112.7)'

    # =========================================================================
    # Transit-compartment absorption (Savic 2007, Bienczak 2016 reference
    # [39]). The model uses NN = 25 transit compartments fixed at the
    # bootstrap median estimate, with mean transit time MTT and absorption
    # rate constant ka reported separately per trial. This packaged file
    # encodes the CHAPAS-3 reference values (MTT = 0.82 h, ka = 0.79 /h)
    # and does NOT carry the trial-of-origin (CHAPAS-3 vs ARROW) covariate;
    # see covariatesDataExcluded$STUDY_ARROW notes and the vignette
    # Assumptions and deviations for the rationale and the ARROW typical
    # values that downstream users can apply if they need them.
    # =========================================================================
    lmtt <- log(0.82)
    label("Mean transit time MTT through the NN-compartment Savic chain (h), CHAPAS-3 reference")
    # Bienczak 2016 Table 2 row 'MTT CHAPAS-3 = 0.82 (0.69-0.96)'

    lka <- log(0.79)
    label("First-order absorption rate ka from the depot to the central compartment (1/h), CHAPAS-3 reference")
    # Bienczak 2016 Table 2 row 'Ka CHAPAS-3 = 0.79 (0.37-0.95)'

    nn_fix <- fixed(25)
    label("Number of Savic-style transit compartments NN (integer, unitless)")
    # Bienczak 2016 Table 2 row 'NN = 25.0 (17.7-35.1)' -- bootstrap median used as the fixed value; rxode2 transit() supports continuous NN but the integer-rounded median preserves the published structural model with negligible numerical impact at this NN

    lfdepot <- fixed(log(1))
    label("Oral bioavailability F/F (unitless; FIXED to 1 because no intravenous data were available)")
    # Bienczak 2016 Table 2 row 'BIO = 1 (FIXED)'; paper Results 'Population pharmacokinetics' paragraph 1: 'PK parameters were estimated relative to oral bioavailability whose typical value was fixed to one due to lack of intravenous data.'

    # =========================================================================
    # Allometric exponents. Fixed by convention (Anderson and Holford,
    # reference 43 of the paper); the paper Methods 'Covariates' paragraph
    # 1: 'Allometric scaling was added to the model at an early development
    # stage as previously suggested [43].'
    # =========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F and Q/F with body weight (unitless)")
    # Bienczak 2016 Methods 'Covariates' paragraph 1 (Anderson and Holford 2008, ref [43])

    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on Vc/F and Vp/F with body weight (unitless)")
    # Bienczak 2016 Methods 'Covariates' paragraph 1 (Anderson and Holford 2008, ref [43])

    # =========================================================================
    # CYP2B6 SNP-vector effects on CL/F (Bienczak 2016 Table 2 / Table 3).
    # Encoded as log-ratio multiplicative effects with the 516GG|983TT
    # extensive-metabolizer phenotype (CL/F = 6.94 L/h) as the structural
    # reference. The five non-reference observed combinations enter via
    # mutually-exclusive indicator products built in model() from the two
    # SNP allele-count covariate columns. The paper's reported per-group
    # CL/F typical values are reproduced exactly: e.g. for 516GT|983TT the
    # encoded log-effect e_GT_TT_cl = log(4.90/6.94) = -0.3478 gives the
    # 29.4% lower CL/F reported in Bienczak 2016 Table 2 / Discussion
    # paragraph 2 ('presence of one variant allele in 516G>T causes
    # clearance to drop by 34%').
    # =========================================================================
    e_GG_TC_cl <- log(3.93 / 6.94)
    label("Log-ratio of 516GG|983TC CL/F vs the 516GG|983TT EM reference (unitless)")
    # Bienczak 2016 Table 2 row 'CL 516GG|983TC = 3.93 (2.61-5.65)'; log(3.93/6.94) = -0.5689

    e_GG_CC_cl <- log(0.74 / 6.94)
    label("Log-ratio of 516GG|983CC CL/F vs the 516GG|983TT EM reference (unitless)")
    # Bienczak 2016 Table 2 row 'CL 516GG|983CC = 0.74 (0.72-0.75)'; log(0.74/6.94) = -2.2384

    e_GT_TT_cl <- log(4.90 / 6.94)
    label("Log-ratio of 516GT|983TT CL/F vs the 516GG|983TT EM reference (unitless)")
    # Bienczak 2016 Table 2 row 'CL 516GT|983TT = 4.90 (4.40-5.46)'; log(4.90/6.94) = -0.3478

    e_GT_TC_cl <- log(1.36 / 6.94)
    label("Log-ratio of 516GT|983TC CL/F vs the 516GG|983TT EM reference (unitless)")
    # Bienczak 2016 Table 2 row 'CL 516GT|983TC = 1.36 (0.97-1.76)'; log(1.36/6.94) = -1.6296

    e_TT_TT_cl <- log(1.92 / 6.94)
    label("Log-ratio of 516TT|983TT CL/F vs the 516GG|983TT EM reference (unitless)")
    # Bienczak 2016 Table 2 row 'CL 516TT|983TT = 1.92 (1.52-2.33)'; log(1.92/6.94) = -1.2854

    # =========================================================================
    # Inter-individual variability (Bienczak 2016 Table 2 'Random Effects
    # (ETA)' column). Source reports approximate %CV on the SD scale
    # (footnote: 'Expressed as approximate %CV on SD scale sqrt(ETA) * 100');
    # nlmixr2lib uses variances on the log scale, converted via
    # omega^2 = log(1 + CV^2). The source distinguishes BSV (between-
    # subject) from BOV (between-occasion). nlmixr2lib has no idiomatic
    # encoding for BOV separate from BSV; per the convention used in
    # Bienczak_2016_nevirapine.R and Svensson_2018_bedaquiline.R, BOV is
    # dropped where a BSV term is reported on the same parameter and
    # folded in as a BSV-equivalent where only BOV is reported.
    # Specifically:
    #   - CL/F:  BSV 36.9% kept; BOV 26.6% dropped (separate BSV reported)
    #   - F:     BSV 42.2% kept; BOV 50.5% dropped (separate BSV reported)
    #   - MTT:   only BOV 78.0% reported -> folded as BSV-equivalent
    #   - ka:    only BOV 57.7% reported -> folded as BSV-equivalent
    # =========================================================================
    etalcl     ~ 0.12779
    # Bienczak 2016 Table 2 BSV CL/F = 36.9%; omega^2 = log(1 + 0.369^2) = 0.12779
    etalfdepot ~ 0.16415
    # Bienczak 2016 Table 2 BSV BIO = 42.2%; omega^2 = log(1 + 0.422^2) = 0.16415
    etalmtt    ~ 0.47521
    # Bienczak 2016 Table 2 BOV MTT = 78.0% folded as BSV-equivalent; omega^2 = log(1 + 0.780^2) = 0.47521
    etalka     ~ 0.28738
    # Bienczak 2016 Table 2 BOV Ka = 57.7% folded as BSV-equivalent; omega^2 = log(1 + 0.577^2) = 0.28738

    # =========================================================================
    # Residual error (Bienczak 2016 Table 2 'Error model (SIGMA)' rows).
    # The paper additionally reports a 2-fold increase in the residual
    # error for sparse PK samples (Table 2 row 'Increased error for
    # sparse data = 2x (1.7x - 2.5x)'); this per-occasion / per-record
    # residual-error scaling is dropped here, leaving the typical
    # (intensive PK) residual-error magnitude. See vignette Assumptions
    # and deviations.
    # =========================================================================
    addSd <- 0.101
    label("Additive residual error (mg/L)")
    # Bienczak 2016 Table 2 row 'Additive error (mg/L) = 0.101 (0.067-0.131)'

    propSd <- 0.0672
    label("Proportional residual error (fraction)")
    # Bienczak 2016 Table 2 row 'Proportional error (%) = 6.72 (5.20-7.90)'
  })

  model({
    # --- 1. CYP2B6 SNP-vector indicators -----------------------------------
    # Six observed combinations of (516G>T T-count, 983T>C C-count) define
    # the metabolic subgroups in Bienczak 2016 Table 2. Built as mutually-
    # exclusive products of allele-count equality checks. The 516GG|983TT
    # reference is implicit: when all five non-reference indicators are 0,
    # the log-CL/F equals lcl_GG_TT.
    is_GG_TC <- (SNP_CYP2B6_RS3745274_T_COUNT == 0) * (SNP_CYP2B6_RS28399499_C_COUNT == 1)
    is_GG_CC <- (SNP_CYP2B6_RS3745274_T_COUNT == 0) * (SNP_CYP2B6_RS28399499_C_COUNT == 2)
    is_GT_TT <- (SNP_CYP2B6_RS3745274_T_COUNT == 1) * (SNP_CYP2B6_RS28399499_C_COUNT == 0)
    is_GT_TC <- (SNP_CYP2B6_RS3745274_T_COUNT == 1) * (SNP_CYP2B6_RS28399499_C_COUNT == 1)
    is_TT_TT <- (SNP_CYP2B6_RS3745274_T_COUNT == 2) * (SNP_CYP2B6_RS28399499_C_COUNT == 0)

    # Log-additive composite genotype effect on log-CL/F (516GG|983TT is
    # the reference; all five terms are zero for the reference subject).
    log_cl_geno <- e_GG_TC_cl * is_GG_TC +
                   e_GG_CC_cl * is_GG_CC +
                   e_GT_TT_cl * is_GT_TT +
                   e_GT_TC_cl * is_GT_TC +
                   e_TT_TT_cl * is_TT_TT

    # --- 2. Individual PK parameters ---------------------------------------
    # Apparent oral clearance: 516GG|983TT typical CL/F times the
    # log-additive genotype factor, allometrically scaled to body weight
    # (reference 15.4 kg), with BSV on the log scale.
    cl <- exp(lcl_GG_TT + log_cl_geno + etalcl) * (WT / 15.4)^e_wt_cl

    # Central volume of distribution: allometric on body weight; no IIV
    # reported in Bienczak 2016 Table 2.
    vc <- exp(lvc) * (WT / 15.4)^e_wt_vc

    # Inter-compartmental clearance and peripheral volume, allometric on
    # body weight; no IIV reported in Bienczak 2016 Table 2.
    q  <- exp(lq)  * (WT / 15.4)^e_wt_cl
    vp <- exp(lvp) * (WT / 15.4)^e_wt_vc

    # --- 3. Absorption parameters ------------------------------------------
    # MTT and ka use the CHAPAS-3 reference typical values; BOV folded as
    # BSV-equivalent (see ini() comment block on inter-individual
    # variability).
    mtt    <- exp(lmtt + etalmtt)
    ka     <- exp(lka  + etalka)
    nn     <- nn_fix
    fdepot <- exp(lfdepot + etalfdepot)

    # --- 4. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp

    # --- 5. ODE system -----------------------------------------------------
    # Dose is given to the depot compartment in the dataset (cmt = "depot",
    # amt = dose). The Savic 2007 analytical Erlang transit chain
    # rxode2 builtin transit(nn, mtt, fdepot) returns the per-time input
    # rate to depot from the most recent dose, applying the bioavailability
    # factor fdepot. f(depot) <- 0 below suppresses the standard bolus so
    # transit() is the only input pathway, and depot then empties at
    # first-order rate ka into central. The two-compartment disposition
    # is parameterised with rate constants kel = CL/Vc, k23 = Q/Vc, and
    # k32 = Q/Vp.
    d/dt(depot)      <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central)    <-  ka * depot - kel * central - k23 * central + k32 * peripheral1
    d/dt(peripheral1) <-                              k23 * central - k32 * peripheral1

    # Suppress the dose bolus on depot so the analytical transit() chain
    # is the only input pathway from the most recent dose.
    f(depot) <- 0

    # --- 6. Observation and residual error ---------------------------------
    # Dose in mg / volume in L -> concentration in mg/L (the paper's
    # reported unit). Combined additive (0.101 mg/L) and proportional
    # (6.72%) residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
