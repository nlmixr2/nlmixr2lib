Chala_2023_efavirenz <- function() {
  description <- "One-compartment population pharmacokinetic-pharmacogenetic model with first-order absorption for oral efavirenz in antiretroviral-naive HIV-1-infected Ethiopian children aged 3-16 years (Chala 2023). Apparent oral clearance CL/F carries five multiplicative covariate factors: allometric body weight (exponent fixed to 0.75, reference 22 kg), CYP2B6*6 (c.516G>T, rs3745274) heterozygous and homozygous genotype factors, an ABCB1 c.4036A>G (rs3842) A-allele-carrier factor, a genotype-gated autoinduction step (CL/F rises from week 12 in CYP2B6*1/*1 and from week 8 in CYP2B6*1/*6, with no change in *6/*6), and a two-class latent mixture in which 7.5% of children form a subpopulation with 3.5-fold lower CL/F. Apparent volume V/F scales linearly with weight (exponent fixed to 1). Absorption rate ka was not identifiable from the sparse design and is fixed to 0.776 1/h; interindividual variability on V/F and ka was likewise not identifiable and was fixed to zero, leaving CL/F as the only random effect."
  reference <- paste(
    "Chala A, Kitabi EN, Ahmed JH, Tadesse BT, Chaka TE, Makonnen E,",
    "Aklillu E.",
    "Genetic and non-genetic factors influencing efavirenz population",
    "pharmacokinetics among human immunodeficiency virus-1-infected",
    "children in Ethiopia.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(6):783-794.",
    "doi:10.1002/psp4.12951."
  )
  vignette <- "Chala_2023_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot   = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometrically scaled on both CL/F and V/F with the reference weight fixed at the cohort median, 22 kg (Chala 2023 Results 'Population pharmacokinetic model' step 3: 'CLi = CLpop * (WT/22)^0.75, Vi = Vpop * (WT/22)'; supplement Appendix S1 control stream `TVCL = THETA(1) * ... * (WT/22)**0.75` and `TVV = THETA(2)*(WT/22)`). Both exponents were FIXED to the theoretical allometric values 0.75 and 1 rather than estimated, which the Discussion (paragraph 3) states was done deliberately because the pediatric weight range was narrow. Cohort distribution: median 22.05 kg, IQR 16.8-28.25 kg (Table 1).",
      source_name        = "WT"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type = CYP2B6*1/*1, 1 = GT heterozygous = CYP2B6*1/*6, 2 = TT homozygous variant = CYP2B6*6/*6. Chala 2023 genotypes and reports the variant as 'CYP2B6*6' throughout and identifies it with the 516G>T SNP explicitly in the Discussion ('the CYP2B6*6 (c.516G>T), which is the most common single-nucleotide polymorphism in Sub-Saharan Africans'), so *6 status maps one-to-one onto 516G>T genotype in this cohort.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (germline genotype). The paper encodes two mutually-exclusive binary indicators as multiplicative fold-factors on CL/F, which the canonical count column reconstructs deterministically: *1/*6 = (count == 1), *6/*6 = (count == 2); *1/*1 (count == 0) is the reference. Supplement Appendix S1 control stream: `IF(CP2BS6.LE.1) CP2B6S1S1 = 1; IF(CP2BS6.EQ.2) CP2B6S1S6 = 1; IF(CP2BS6.EQ.3) CP2B6S6S6 = 1; CP2B6CL = THETA(6)**CP2B6S1S6 * THETA(7)**CP2B6S6S6`. The same column ALSO gates the autoinduction term (see T_FIRSTDOSE): the week-12 step applies only to *1/*1 and the week-8 step only to *1/*6, with no time effect in *6/*6. Cohort distribution (Table 1, n = 100): *1/*1 45 (45%), *1/*6 45 (45%), *6/*6 8 (8%), missing 2 (2%). The control stream maps missing genotype (coded 0 or 1) onto the *1/*1 reference via `.LE.1`.",
      source_name        = "CP2BS6 (1 = *1/*1, 2 = *1/*6, 3 = *6/*6)"
    ),
    SNP_ABCB1_RS3842_A_CARRIER = list(
      description        = "Binary carrier indicator for the ABCB1 (MDR1 / P-glycoprotein) c.4036A>G (rs3842, 3' untranslated region) polymorphism, pooled as Chala 2023 pools it: 1 = subject carries at least one A allele (heterozygous G/A or homozygous A/A); 0 = homozygous G/G. Note that this is the OPPOSITE pooling to the register's older `SNP_ABCB1_RS3842` column, which is a G-allele-carrier indicator with A/A as reference (Mukonzo 2009); the two cannot be derived from one another because the heterozygote falls in the '1' group under both.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous G/G)",
      notes              = "Time-fixed per subject (germline genotype). A-allele carriers have 1.45-fold HIGHER apparent oral clearance than G/G homozygotes (Chala 2023 Table 2 row 'Fold of typical CL for subjects with ABCB1.rs3842 G/A or A/A genotype'). Supplement Appendix S1 control stream: `IF(ABRS3842.EQ.1) ABCB1RS3842 = 1; ABCB1RS3842CL = THETA(8)**ABCB1RS3842`. Cohort distribution (Table 1, n = 100): G/G 17 (17%), G/A or A/A 80 (80%), missing 3 (3%). Direction check against the other rs3842 model in the registry: Mukonzo 2009 (Ugandan adults) found G-allele carriers had 25.7% HIGHER relative bioavailability, i.e. higher exposure; Chala 2023 finds G/G homozygotes have the LOWER clearance, i.e. also higher exposure. The two papers therefore agree that the G allele tracks with higher efavirenz exposure -- they differ only in which genotypes they pool against which reference, which is why a second canonical column is needed rather than a sign flip on the first.",
      source_name        = "ABRS3842 (1 = G/A or A/A; reference G/G)"
    ),
    T_FIRSTDOSE = list(
      description        = "Time elapsed since the first efavirenz dose of the antiretroviral-therapy course (treatment duration), NOT time after the most recent dose",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the genotype-gated autoinduction step on CL/F. Divided by 168 in model() to convert the canonical hours to the weeks-on-treatment clock the paper uses. Two mutually-exclusive step effects: CL/F is 12.35% higher from week 12 onward in CYP2B6*1/*1 subjects, and 22.98% higher from week 8 onward in CYP2B6*1/*6 subjects; CYP2B6*6/*6 subjects show no change with treatment duration (Chala 2023 Table 2 and Results 'Population pharmacokinetic model' paragraph 5). The supplement Appendix S1 control stream implements this on the DISCRETE sampling grid actually observed -- `CLWKCP2B61 = ((1+THETA(9))**WK12 * (1+THETA(9))**WK24)**CP2B6S1S1` with `WK12 = (WEEK.EQ.12)` and `WK24 = (WEEK.GE.24)`, and correspondingly `WK8 = (WEEK.EQ.8)` for the *1/*6 arm -- because PK samples were only drawn at weeks 0/1, 4, 8, 12, 24 and 48. On that grid the discrete coding and the >= threshold coding used here are identical; the >= form is what Table 2 ('from >= 12-weeks', 'from >= 8 weeks'), Equation 1 and the abstract describe, and it is the only form that is well defined at the intermediate times a simulation may request. See the vignette Errata.",
      source_name        = "WEEK (weeks since start of efavirenz-based ART; sampled values 0, 4, 8, 12, 24, 48)"
    ),
    MIX_SLOW_ELIM_EFV = list(
      description        = "Per-subject latent mixture-model class indicator from the Chala 2023 NONMEM $MIX block. 1 = subject belongs to the minority slow-eliminator subpopulation whose typical CL/F is 0.2829 times (about 3.5-fold lower than) the main population; 0 = subject belongs to the majority main population (the reference). Not a measured clinical covariate -- the assignment is a latent-class index, and no measured covariate in the screen identified it.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (majority main population, 92.5% of the source cohort)",
      notes              = "Estimated population probability of MIX_SLOW_ELIM_EFV = 1 is 0.07511 (RSE 71%; Chala 2023 Table 2 'Proportion of an unknown subpopulation among the studied population' 0.075, bootstrap median 0.09, 95% CI 0.014-0.272; supplement Appendix S1 `$MIX NSPOP=2; P(1) = THETA(11)` with THETA(11) final estimate 0.07511 in Table S3 run32). For typical-value simulation set MIX_SLOW_ELIM_EFV = 0 (the recommended default) to reproduce the main population; for population simulation draw MIX_SLOW_ELIM_EFV ~ Bernoulli(0.075) per subject, which is what the paper itself did for its Figure 3 Monte-Carlo simulations (Methods 'Population pharmacokinetic analysis' final paragraph: 'The virtual subjects were also randomly assigned to an unknown subpopulation at the probability of 8%'). A bimodal histogram of the CL/F EBEs motivated the mixture (Results paragraph 6, Figure S1); the Discussion (paragraph 3) speculates the subpopulation 'might represent unstudied covariates like other CYP2B6 and CYP2A6 polymorphisms', neither of which was genotyped.",
      source_name        = "MIXNUM (NONMEM $MIX class index; subpopulation 1 = the slow-CL class)"
    )
  )

  # Covariates screened during covariate-model building but NOT retained in the
  # final model (Chala 2023 Methods 'Population pharmacokinetic analysis'
  # paragraph 3 lists the full screen; Results 'Population pharmacokinetic
  # model' paragraph 4 reports the outcomes quoted below). Documented so the
  # provenance of the screen is preserved; none is referenced in model().
  covariatesDataExcluded <- list(
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 c.3435C>T (rs1045642) mutant-allele carrier indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a third genotype covariate on CL/F after CYP2B6*6 and rs3842; dropped because the improvement in fit was not significant. Results paragraph 4: 'Stepwise additions of CYP2B6*6, ABCB1c.4036A>G, and ABCB1c.3435C>T genotypes as covariates for CL resulted in -28.5, -12.1, and -1.2 units decrease in OFV, respectively. Therefore CYP2B6*6 and ABCB1c.4036A>G were retained in the model.' Cohort distribution (Table 1): C/C 65 (65%), T/C or T/T 32 (32%), missing 3 (3%)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male) (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Results paragraph 4: 'ETA versus covariate plots indicated sex, eGFR, and PTB could be potential covariates, but their inclusion into the model did not improve the fit.' No point estimate is published, so no effect can be encoded. Cohort distribution (Table 1): male 58 (58%), female 42 (42%)."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (screened, not retained)",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Showed a visual trend in the ETA-versus-covariate plots but did not improve the fit when added (Results paragraph 4). No point estimate is published. Cohort distribution (Table 1): median 95.453, IQR 70.7-115.7 mL/min/1.73m^2. Efavirenz is cleared by hepatic CYP2B6 metabolism, so a renal-function effect had no strong mechanistic prior."
    ),
    TB_POS = list(
      description = "Active pulmonary tuberculosis disease indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Showed a visual trend in the ETA-versus-covariate plots but did not improve the fit when added (Results paragraph 4). No point estimate is published. Cohort distribution (Table 1): no 83 (83%), yes 15 (15%), missing 2 (2%)."
    ),
    AGE = list(
      description = "Age (screened, not retained)",
      units       = "year",
      type        = "continuous",
      notes       = "Listed in the Methods covariate screen ('demographics (age and sex)') but not reported as retained or as showing a trend. Body weight carried the size effect. Cohort distribution (Table 1): median 9 years, IQR 6-13 years."
    ),
    CYP3A5_STAR1_HET = list(
      description = "CYP3A5 one functional (*1) allele indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the metabolizing-enzyme genotype panel screened on CL/F (Methods paragraph 3: 'genotypes of EFV metabolizing enzymes (CYP2B6*6, CYP3A5 [*3, *6, *7], and UGT2B7c.372G>A)'). Not retained; no point estimate is published. Cohort distribution (Table 1, CYP3A5 number of functional *1 alleles): zero 54 (54%), one 38 (38%), two 6 (6%), missing 2 (2%)."
    ),
    CYP3A5_STAR1_HOM = list(
      description = "CYP3A5 two functional (*1) alleles indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Companion indicator to CYP3A5_STAR1_HET; see that entry. Cohort distribution (Table 1): two functional alleles in 6 (6%) of subjects."
    ),
    UGT2B7_372AG = list(
      description = "UGT2B7 c.372G>A carrier indicator (A/G or A/A versus G/G reference) (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the metabolizing-enzyme genotype panel screened on CL/F (Methods paragraph 3). Not retained; no point estimate is published. Cohort distribution (Table 1): G/G 28 (28%), A/G or A/A 69 (69%), missing 3 (3%). Chala 2023 pools heterozygotes with variant homozygotes, so this documentation entry uses the register's carrier-style `UGT2B7_372AG` name for the pooled non-reference group rather than the strict heterozygote-only semantics that name carries elsewhere; the entry is documentation only and is not referenced in model()."
    ),
    SLCO1B1_HAP15_HET = list(
      description = "SLCO1B1 functional-allele indicators, screened via SLCO1B1*1B, SLCO1B1*5 and SLCO1B1 rs4149032 (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the drug-transporter genotype panel screened on CL/F (Methods paragraph 3: 'drug transporters (ABCB1c.3435C>T, ABCB1 c.4036A>G [rs3842], and SLCO1B1*1B, SLCO1B1*5, SLCO1B1.rs4149032)'). Not retained; no point estimate is published. Cohort distribution (Table 1, SLCO1B1 number of functional *1 alleles): zero 73 (73%), one or two 24 (24%), missing 3 (3%)."
    ),
    ALB = list(
      description = "Plasma albumin (screened, not retained)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "One of the baseline laboratory covariates screened on CL/F (Methods paragraph 3). Not retained; no point estimate is published. Cohort distribution (Table 1): median 3.8 mg/dL, IQR 3.4-4.2. The full baseline laboratory screen also covered AST, ALT, ALP, total cholesterol, urea, total bilirubin, eGFR, hemoglobin, hematocrit, LDL, HDL, triglycerides, viral load and CD4 count; and the baseline clinical screen covered isoniazid preventive therapy, hepatitis C infection, hepatitis B surface antigen, cotrimoxazole prophylaxis, pulmonary tuberculosis, WHO HIV clinical stage and ART regimen. None was retained and none has a published point estimate, so they are summarised here and in population$notes rather than given one documentation entry each."
    ),
    CD4_ABS = list(
      description = "Baseline absolute CD4 count (screened, not retained)",
      units       = "cells/dL",
      type        = "continuous",
      notes       = "Screened on CL/F as a baseline laboratory covariate (Methods paragraph 3). Not retained; no point estimate is published. Cohort distribution (Table 1): median 330 cells/dL, IQR 200.5-671."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100,
    n_studies      = 1,
    age_range      = "3-16 years (enrolment criterion); most subjects below 12 years",
    age_median     = "9 years (IQR 6-13)",
    weight_range   = "interquartile range 16.8-28.25 kg (full range not reported; the dosing weight bands in Table S1 span 7.5 kg to above 40 kg)",
    weight_median  = "22.05 kg (IQR 16.8-28.25)",
    sex_female_pct = 42,
    race_ethnicity = c(Ethiopian = 100),
    disease_state  = "Combination-antiretroviral-therapy-naive children with HIV-1 infection, starting efavirenz-based combination ART. Relatively normal liver and renal function and immunocompetent CD4 counts at baseline; 15% had active pulmonary tuberculosis, 69% were on cotrimoxazole prophylaxis, and WHO clinical stage was 1/2/3/4 in 40/22/30/5% of subjects.",
    dose_range     = "Weight-band efavirenz dosing per the SUSTIVA label, once daily. Doses actually received by weight band (Table S1): 7.5-15 kg 200 mg (66.7%), 250 mg (13.3%), 300 mg (20%); 15-20 kg 250 mg (55.6%), 300 mg (44.4%); 20-25 kg 250 mg (5.6%), 300 mg (61.1%), 400 mg (33.3%); 25-32.5 kg 350 mg (13.6%), 400 mg (22.7%), 600 mg (63.6%); 32.5-40 kg 400 mg (7.1%), 600 mg (92.9%); above 40 kg 600 mg (100%). A few individuals received higher doses than the label recommends.",
    regions        = "Ethiopia (seven hospital ART centres in the Oromia and Southern Nations, Nationalities and Peoples regional states)",
    cyp2b6_freq    = "CYP2B6*6 (c.516G>T, rs3745274): *1/*1 45 (45%), *1/*6 45 (45%), *6/*6 8 (8%), missing 2 (2%). For its Monte-Carlo simulations the paper used 45%/45%/10% (Methods final paragraph).",
    abcb1_freq     = "ABCB1 c.4036A>G (rs3842): G/G 17 (17%), G/A or A/A 80 (80%), missing 3 (3%). For its Monte-Carlo simulations the paper used 20% G/G and 80% G/A or A/A. ABCB1 c.3435C>T: C/C 65 (65%), T/C or T/T 32 (32%), missing 3 (3%).",
    notes          = "554 efavirenz plasma concentrations from 100 children, collected over a span of one year (Results 'Population pharmacokinetic model' paragraph 1, Figure 1). Two-arm sampling design: 13 children gave rich samples at 0, 2.5, 16 and 24 h after the FIRST dose (9 of them repeated at week 8), and 87 children gave a single mid-dose sample between 8 and 16 h post dose at weeks 4, 8, 12, 24 and 48. The study was planned for 30 rich and 90 sparse subjects but only 13 and 87 consented, and the resulting data sparseness is why interindividual variability could not be estimated for V/F or ka. Efavirenz was assayed by HPLC-MS/MS with LLOQ 15.78 ng/mL over a calibration range of 15.78-15,783.75 ng/mL. NONMEM 7.5.0 with Pirana 2.9.9, PsN 5.2.6 and R 4.1.2; FOCE-I estimation, ADVAN2 TRANS2, model qualified by goodness-of-fit plots, prediction-corrected VPC and bootstrap. Covariate screen (Methods paragraph 3): genotypes of CYP2B6*6, CYP3A5 (*3, *6, *7), UGT2B7 c.372G>A, ABCB1 c.3435C>T, ABCB1 c.4036A>G (rs3842), SLCO1B1*1B, SLCO1B1*5 and SLCO1B1 rs4149032; demographics age and sex; baseline laboratory AST, ALT, ALP, total cholesterol, urea, total bilirubin, eGFR, albumin, hemoglobin, hematocrit, LDL, HDL, triglycerides, viral load and CD4 count; and baseline clinical isoniazid preventive therapy, hepatitis C infection, hepatitis B surface antigen, cotrimoxazole prophylaxis, pulmonary tuberculosis, WHO HIV clinical stage and ART regimen. Interoccasion variability on CL/F was explored (occasions = weeks 0, 4, 8, 12 and >= 24) and found to be negligible. Food effect was not assessed; the study could not monitor food co-administration compliance over the one-year follow-up (Discussion 'limitations' paragraph)."
  )

  ini({
    # ---- Structural PK parameters (Chala 2023 Table 2 final population PK model) ----
    # Reference subject: 22 kg (cohort median weight), CYP2B6*1/*1, ABCB1
    # rs3842 G/G, week 1 of treatment, main (non-mixture) population.
    lcl <- log(4.30)   ; label("Apparent oral clearance CL/F (L/h) at the 22 kg / CYP2B6*1/*1 / rs3842 G/G / week-1 / main-population reference")  # Chala 2023 Table 2 'CL (L/h)' 4.30 (RSE 13%; bootstrap median 4.30, 95% CI 3.20-5.30); Equation 1 prints '4.3(CLpop)'. Table S3 run32 (the same final run) prints 4.29 -- see vignette Errata.
    lvc <- log(123.80) ; label("Apparent volume of distribution V/F (L) at the 22 kg reference")  # Chala 2023 Table 2 'V c (L)' 123.80 (RSE 10%; bootstrap median 125.24, 95% CI 103.30-163.80)
    lka <- fixed(log(0.776)) ; label("First-order absorption rate constant ka (1/h)")  # Chala 2023 Table 2 'K a (/h) 0.78' with footnote 'a Fixed to this value'; the unrounded 0.776 is given in Results paragraph 7 ('K a were 4.3 L/h, 124 L, and 0.776/h'), the abstract, and the supplement Appendix S1 control stream `$THETA (0.776) FIX ; KA`. Results paragraph 3 step 5: a sensitivity analysis fixing ka 10-fold lower or 5-fold higher worsened the fit, so ka was fixed at its step-4 estimate.

    # ---- Allometric exponents on body weight, reference 22 kg ----
    # Both FIXED to theoretical values rather than estimated (Discussion
    # paragraph 3: 'the allometric scaling of CL and V was implemented with
    # allometric exponents fixed to theoretical values of 0.75 and 1 for CL
    # and V, respectively'). They appear as literal constants in the
    # supplement Appendix S1 control stream, not as $THETA entries.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F (unitless)")  # Chala 2023 Results paragraph 3 step 3 'CLi = CLpop * (WT/22)^0.75'; Appendix S1 `(WT/22)**0.75`
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on V/F (unitless)")   # Chala 2023 Results paragraph 3 step 3 'Vi = Vpop * (WT/22)'; Appendix S1 `TVV = THETA(2)*(WT/22)`

    # ---- Covariate effects on CL/F ----
    # Chala 2023 Equation 1 writes CL/F as a product of factors, each raised to
    # a 0/1 indicator power so it collapses to 1 when the indicator is 0:
    #   CLi = 4.3 * (WT/22)^0.75 * 0.72^f1 * 0.28^f2 * 1.45^f3
    #             * 1.12^f4 * 1.23^f5 * 0.28^f6 * exp(eta)
    # The genotype and mixture factors are FOLDS of the typical value; the two
    # autoinduction thetas are PROPORTIONAL INCREASES entered as (1 + theta).
    # Point estimates below are taken from supplement Table S3 run32 -- the
    # final run, carrying all twelve thetas -- which prints the same values as
    # main-text Table 2 to more significant figures.
    e_cyp2b6_6het_cl <- 0.7245 ; label("Fold of typical CL/F for CYP2B6*1/*6 heterozygotes (unitless)")  # Chala 2023 Table 2 'Fraction of typical CL for subjects with CYP2B6*1/*6' 0.72 (RSE 11%; bootstrap median 0.73, 95% CI 0.59-0.89); Table S3 run32 CP2B6S1S6 = 0.7245. Abstract: clearance reduced by 28% in *1/*6.
    e_cyp2b6_6hom_cl <- 0.2823 ; label("Fold of typical CL/F for CYP2B6*6/*6 homozygotes (unitless)")     # Chala 2023 Table 2 'Fraction of typical CL for subjects with CYP2B6*6/*6' 0.28 (RSE 17%; bootstrap median 0.28, 95% CI 0.2-0.46); Table S3 run32 CP2B6S6S6 = 0.2823. Abstract: clearance reduced by 72% in *6/*6.
    e_snp_abcb1_rs3842_a_cl <- 1.452 ; label("Fold of typical CL/F for ABCB1 rs3842 G/A or A/A carriers (unitless)")  # Chala 2023 Table 2 'Fold of typical CL for subjects with ABCB1.rs3842 G/A or A/A genotype' 1.45 (RSE 12%; bootstrap median 1.47, 95% CI 1.15-1.87); Table S3 run32 ABCB1RS3842 = 1.452

    # Genotype-gated autoinduction. Efavirenz induces its own metabolism, and
    # the induction appears at a different treatment week in each CYP2B6*6
    # genotype (Results paragraph 5): significant from week 12 in *1/*1, from
    # week 8 in *1/*6, and absent in *6/*6.
    e_autoind_cyp2b6_6wt_cl  <- 0.1235 ; label("Proportional increase in CL/F from week 12 onward in CYP2B6*1/*1 (fraction)")  # Chala 2023 Table 2 'Proportional increase in CL from >= 12-weeks for subjects with CYP2B6*1/*1' 0.12 (RSE 71%; bootstrap median 0.12, 95% CI 0.001-0.27); Table S3 run32 CP2B6S1S1WEEK12 = 0.1235. Appendix S1: `CLWKCP2B61 = ((1+THETA(9))**WK12 * (1+THETA(9))**WK24)**CP2B6S1S1`.
    e_autoind_cyp2b6_6het_cl <- 0.2298 ; label("Proportional increase in CL/F from week 8 onward in CYP2B6*1/*6 (fraction)")   # Chala 2023 Table 2 'Proportional increase in CL from >= 8weeks for subjects with CYP2B6*1/*6' 0.23 (RSE 39%; bootstrap median 0.23, 95% CI 0.07-0.40); Table S3 run32 CP2B6S1S6WEEK8 = 0.2298. Appendix S1: `CLWKCP2B62 = ((1+THETA(10))**WK8 * (1+THETA(10))**WK12 * (1+THETA(10))**WK24)**CP2B6S1S6`.

    # Latent two-class mixture on CL/F. The estimated class probability itself
    # (0.07511) is not a model() term; it is recorded in
    # covariateData$MIX_SLOW_ELIM_EFV$notes so a simulation can draw the class.
    e_mix_slow_elim_efv_cl <- 0.2829 ; label("Fold of typical CL/F for the latent slow-eliminator subpopulation (unitless)")  # Chala 2023 Table 2 'Fraction of typical of CL for the subpopulation compared to the studied population' 0.28 (RSE 36%; bootstrap median 0.29, 95% CI 0.13-0.535); Table S3 run32 MIXPOP = 0.2829. Appendix S1: `CLMIX = THETA(12)**MIXPOP1`.

    # ---- IIV ----
    # CL/F is the ONLY random effect in the final model. Results paragraph 7:
    # 'The estimate of IIV of CL was 35.4%. For V and Ka, IIVs could not be
    # estimated with reasonable precision and were therefore fixed to 0.'
    # The Appendix S1 $OMEGA block confirms: `0.227 ; IIVCL / 0 FIX ; IIVV /
    # 0.000001 FIX ; IIVKA`. The two zero-variance etas are omitted here
    # rather than written as `~ fixed(0)`, which would make OMEGA singular
    # and break rxode2's Cholesky sampler.
    #
    # Methods paragraph 2 defines the CV-to-variance relationship used by the
    # paper: CV% = 100 * sqrt(exp(omega^2) - 1), so omega^2 = log(CV^2 + 1).
    etalcl ~ 0.118127  # Chala 2023 Table 2 'Interindividual variability for CL (%CV)' 35.4 (RSE 13%); Table S3 run32 IIVCL = 0.3541 i.e. 35.41% CV. log(0.3541^2 + 1) = 0.118127. ETA shrinkage for CL was 20.7%. Scale check: the SAME Table 2 row reports its bootstrap as 0.105 (95% CI 0.05-0.169), which is on the VARIANCE scale, and the back-transformed 0.118127 sits inside that interval -- whereas reading 0.3541 as a variance (CV 65%) would fall far outside it. That independently confirms 0.3541 is a CV and not an omega^2.

    # ---- Residual error (combined proportional + additive) ----
    # Appendix S1 $ERROR: `PROP = THETA(4); ADD = THETA(5);
    # W = SQRT(ADD**2 + PROP**2*IPRED**2); Y = IPRED + W*ERR(1)` with
    # `$SIGMA 1 FIX`. Main-text Table 2 reports only the proportional term;
    # the additive term is FIXED at 0.00028 mg/L (0.28 ng/mL) in the control
    # stream and acts as a numerical stabiliser -- it is ~56-fold below the
    # assay LLOQ of 15.78 ng/mL and is negligible at therapeutic
    # concentrations of 1-4 ug/mL. It is carried here for fidelity to the
    # published control stream.
    propSd <- 0.4969   ; label("Proportional residual error (fraction)")  # Chala 2023 Table 2 'Proportional residual error (%CV)' 50% (RSE 4%; bootstrap median 0.49, 95% CI 0.45-0.54); Table S3 run32 PROP = 0.4969
    addSd  <- fixed(0.00028) ; label("Additive residual error (ug/mL)")   # Chala 2023 supplement Appendix S1 `$THETA (0.00028) FIX ; ADD`; not reported in main-text Table 2
  })

  model({
    # 1. Decompose the canonical 0/1/2 CYP2B6 516G>T allele count into the
    #    three CYP2B6*6 genotype indicators the paper uses. *1/*1 (count 0)
    #    is the reference for the fold-factors but is itself the gate for the
    #    week-12 autoinduction step, so all three are needed.
    #    Appendix S1: IF(CP2BS6.LE.1) *1/*1; IF(CP2BS6.EQ.2) *1/*6;
    #    IF(CP2BS6.EQ.3) *6/*6.
    wt_6  <- (SNP_CYP2B6_RS3745274_T_COUNT == 0)
    het_6 <- (SNP_CYP2B6_RS3745274_T_COUNT == 1)
    hom_6 <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)

    # 2. Treatment duration on the paper's weeks-on-treatment clock. The
    #    canonical T_FIRSTDOSE column is in hours, so divide by 24 * 7 = 168.
    wk <- T_FIRSTDOSE / 168

    # 3. Genotype-gated autoinduction indicators: the week-12 step applies
    #    only to CYP2B6*1/*1 and the week-8 step only to CYP2B6*1/*6, so each
    #    is the product of a genotype indicator and a treatment-week
    #    threshold. CYP2B6*6/*6 has no time effect.
    ind_autoind_wt  <- wt_6  * (wk >= 12)
    ind_autoind_het <- het_6 * (wk >= 8)

    # 4. Typical CL/F, written as the product of factors of Chala 2023
    #    Equation 1, each raised to its 0/1 indicator power so it collapses
    #    to 1 when the indicator is 0.
    tvcl <- exp(lcl) * (WT / 22)^e_wt_cl *
      e_cyp2b6_6het_cl^het_6 *
      e_cyp2b6_6hom_cl^hom_6 *
      e_snp_abcb1_rs3842_a_cl^SNP_ABCB1_RS3842_A_CARRIER *
      (1 + e_autoind_cyp2b6_6wt_cl)^ind_autoind_wt *
      (1 + e_autoind_cyp2b6_6het_cl)^ind_autoind_het *
      e_mix_slow_elim_efv_cl^MIX_SLOW_ELIM_EFV

    # 5. Individual PK parameters. CL/F is the only parameter with IIV.
    ka <- exp(lka)
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc) * (WT / 22)^e_wt_vc

    # 6. Micro-constants
    kel <- cl / vc

    # 7. ODE system: one compartment with first-order absorption (ADVAN2 TRANS2)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 8. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
