Milosheska_2016_lamotrigine <- function() {
  description <- "One-compartment first-order-absorption parent + one-compartment metabolite population pharmacokinetic model for oral lamotrigine (LTG) and its N-2-glucuronide (LTG-glu) in 100 adult Slovenian epilepsy patients on stable mono- or adjunctive therapy (Milosheska 2016). Complete conversion of parent to metabolite is assumed (LTG-N-5-glucuronide is a minor route, < 10% urinary excretion of unchanged drug per the paper's Discussion citing Ref [5]). Parent apparent oral clearance (CL/F) carries a power effect of total body weight and additive-in-fraction effects of smoking, concomitant enzyme-inducing antiepileptic drugs (carbamazepine, phenobarbital, or phenytoin, pooled as CONMED_EIAED), concomitant UGT2B7 inhibitors (valproic acid or sertraline, pooled as CONMED_UGT_INH), Cockcroft-Gault estimated creatinine clearance (deviation from 110 mL/min), and two UGT2B7 SNP genotype categoricals (-161C>T rs7668258 and 372A>G rs28365063). Parent apparent volume (V/F) carries a linear deviation-from-reference weight effect. Metabolite apparent clearance (CL_LTG-glu / F_metab) carries a power weight and linear Cockcroft-Gault CLcr effect; metabolite apparent volume (V_LTG-glu / F_metab) is estimated as a typical value only (no IIV supported by the sparse metabolite data)."
  reference <- paste(
    "Milosheska D, Lorber B, Vovk T, Kastelic M, Dolzan V, Grabnar I.",
    "Pharmacokinetics of lamotrigine and its metabolite N-2-glucuronide:",
    "Influence of polymorphism of UDP-glucuronosyltransferases and drug",
    "transporters. Br J Clin Pharmacol. 2016 Sep;82(3):399-411.",
    "doi:10.1111/bcp.12984.",
    sep = " "
  )
  vignette <- "Milosheska_2016_lamotrigine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight. Enters parent CL/F as a power model with reference 70 kg (Milosheska 2016 Table 4 row 'Body weight' on CL, exponent 0.938); parent V/F as a linear deviation-from-70-kg model (Table 4 row 'Body weight' on V, coefficient 0.0181 per kg); and metabolite CL_LTG-glu / F_metab as a power model with the same 70 kg reference (Discussion / 'Population pharmacokinetic analysis of the lamotrigine-N-2-glucuronide' section, exponent 1.01).",
      source_name        = "Wt"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-smoker)",
      notes              = "Multiplicative effect on parent CL/F: cl *= (1 + 0.340 * SMOKE) per Milosheska 2016 Table 4 row 'Cigarette smoking'. Missing smoking status was imputed to 'non-smoker' in the primary analysis; multiple imputation (paper Discussion 'the mean estimate ... was 0.418') was in good agreement with the imputation-to-non-smoker estimate.",
      source_name        = "Tob"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant enzyme-inducing antiepileptic drug indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no EIAED coadministration)",
      notes              = "Pooled indicator covering carbamazepine, phenobarbital, or phenytoin per Milosheska 2016 Table 4 footnote 'Carbamazepine, phenobarbital or phenytoin (Ind)'. Multiplicative effect on parent CL/F: cl *= (1 + 0.546 * CONMED_EIAED). Oral contraceptives were also observed to induce CL but were excluded from the final model because only n = 2 patients were on them (paper Results paragraph 3).",
      source_name        = "Ind"
    ),
    CONMED_UGT_INH = list(
      description        = "Concomitant UGT-inhibitor coadministration indicator (pooled)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no UGT-inhibitor coadministration)",
      notes              = "Pooled indicator covering valproic acid (VPA, n = 13) or sertraline (n = 2) per Milosheska 2016 Table 4 footnote '*Valproic acid or sertraline (Inh)'. Both drugs are hypothesised to inhibit LTG glucuronidation via competitive inhibition of UGT2B7 (paper Discussion paragraph 4); the paper's small sertraline sample size prevented separate estimation. Multiplicative effect on parent CL/F: cl *= (1 - 0.579 * CONMED_UGT_INH).",
      source_name        = "Inh"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault estimate of creatinine clearance (NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Raw Cockcroft-Gault CLcr in mL/min (not BSA-normalized) per Milosheska 2016 Methods 'Patients and blood sampling' paragraph 2. Linear deviation-from-110 mL/min effect on parent CL/F: cl *= (1 + 0.00328 * (CRCL - 110)), Table 4 row 'CL cr'; and on metabolite CL_LTG-glu / F_metab: cl_gluc *= (1 + 0.00759 * (CRCL - 110)), Discussion 'a decrease by 0.759% per 1 ml deviation from a standard CLcr of 110 ml min-1'. Reference value 110 mL/min is the paper's chosen center (~median in-cohort CLcr was 107.6 mL/min per Table 1). Cohort range 40.84 - 246 mL/min (Table 1).",
      source_name        = "CLcr"
    ),
    UGT2B7_M161CT = list(
      description        = "UGT2B7 -161C>T (rs7668258) heterozygous C/T indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C wild-type). Effective reference-vs-non-reference contrast is realised only when the paired non-reference indicator UGT2B7_M161TT is also 0; both indicators equal 0 encodes the C/C homozygous wild-type reference stratum for the three-level -161C>T genotype categorical.",
      notes              = "Multiplicative effect on parent CL/F: cl *= (1 - 0.0358 * UGT2B7_M161CT) per Milosheska 2016 Table 4 row 'UGT2B7 -161C>T genotype CT vs CC'. Population frequency in the paper's cohort: 48.9% of 94 genotyped patients (Table 2; the paired C/C reference stratum was 22.3%, T/T homozygous variant 28.7%).",
      source_name        = "UGT2B7 -161 C/T"
    ),
    UGT2B7_M161TT = list(
      description        = "UGT2B7 -161C>T (rs7668258) homozygous T/T indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C wild-type). Effective reference-vs-non-reference contrast is realised only when the paired non-reference indicator UGT2B7_M161CT is also 0.",
      notes              = "Multiplicative effect on parent CL/F: cl *= (1 - 0.204 * UGT2B7_M161TT) per Milosheska 2016 Table 4 row 'UGT2B7 -161C>T genotype TT vs CC'. Population frequency: 28.7% of 94 genotyped patients (Table 2).",
      source_name        = "UGT2B7 -161 T/T"
    ),
    UGT2B7_372AG = list(
      description        = "UGT2B7 372A>G (rs28365063) heterozygous A/G indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (A/A wild-type). Effective reference-vs-non-reference contrast is realised only when the paired non-reference indicator UGT2B7_372GG is also 0; both indicators equal 0 encodes the A/A homozygous wild-type reference stratum for the three-level 372A>G genotype categorical.",
      notes              = "Multiplicative effect on parent CL/F: cl *= (1 + 0.194 * UGT2B7_372AG) per Milosheska 2016 Table 4 row 'UGT2B7 372 A > G genotype AG vs AA'. Population frequency: 25.3% of 99 genotyped patients (Table 2; the paired A/A reference stratum was 72.7%, G/G homozygous variant 2.0%).",
      source_name        = "UGT2B7 372 A/G"
    ),
    UGT2B7_372GG = list(
      description        = "UGT2B7 372A>G (rs28365063) homozygous G/G indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (A/A wild-type). Effective reference-vs-non-reference contrast is realised only when the paired non-reference indicator UGT2B7_372AG is also 0.",
      notes              = "Multiplicative effect on parent CL/F: cl *= (1 + 1.17 * UGT2B7_372GG) per Milosheska 2016 Table 4 row 'UGT2B7 372 A > G genotype GG vs AA'. Only n = 2 subjects were GG homozygous in the cohort (2.0% of 99 genotyped); the effect magnitude is by far the largest single-variant effect but the 95% CI is wide (0.448, 2.47).",
      source_name        = "UGT2B7 372 G/G"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained in the final model. Milosheska 2016 Results paragraph 3 lists significant continuous covariates as weight and CLcr on CL and weight on V; age was tested but not carried into the final CL or V models."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained in the final model. Milosheska 2016 Methods 'Model development' lists gender among the tested categorical covariates; Results paragraph 3 does not include gender among the retained final-model covariates."
    ),
    IBW = list(
      description = "Devine ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened univariately but not retained; total body weight (WT) was preferred in the final model. Milosheska 2016 Methods 'Model development' lists IBW among continuous covariates considered for inclusion."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "microkat/L",
      type        = "continuous",
      notes       = "Screened but not retained. Chronic hepatic disease was an exclusion criterion; AST tested for residual hepatic effect."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "microkat/L",
      type        = "continuous",
      notes       = "Screened but not retained. Chronic hepatic disease was an exclusion criterion; ALT tested for residual hepatic effect."
    ),
    UGT1A4_70CC = list(
      description = "UGT1A4 70C>A (rs6755571) homozygous C/C indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested and found not associated with LTG PK (Milosheska 2016 Results paragraph 3 'What This Study Adds': 'No significant association of lamotrigine pharmacokinetics with UGT1A4, SLC22A1 and ABCB1 polymorphisms was observed'). Population frequency in the cohort (Table 2): CC = 93.8%, CA = 6.2%, AA = 0%."
    ),
    SLC22A1_1222 = list(
      description = "SLC22A1 1222G>A (rs628031) genotype (screened, not retained)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Tested but not retained; per Discussion paragraph 8 the observed CL differences across GG / GA / AA strata were small and non-significant (P = 0.600). Population frequencies (Table 2): GG = 39.8%, GA = 45.9%, AA = 14.3%."
    ),
    ABCB1_2677 = list(
      description = "ABCB1 2677G>T/A (rs2032582) genotype (screened, not retained)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Entered forward-inclusion step but dropped during backward elimination (P > 0.01) per Results paragraph 3. Population frequencies (Table 2): GG = 42.6%, GT+GA = 44.7%, TT = 12.8%."
    ),
    ABCB1_3435 = list(
      description = "ABCB1 3435C>T (rs1045642) genotype (screened, not retained)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not retained. Population frequencies (Table 2): CC = 28.0%, CT = 56.0%, TT = 16.0%."
    ),
    ABCB1_1236 = list(
      description = "ABCB1 1236C>T (rs1128503) genotype (screened, not retained)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not retained. Population frequencies (Table 2): CC = 40.0%, CT = 49.0%, TT = 11.0%."
    ),
    ABCB1_HAP_TTT = list(
      description = "ABCB1 T-T-T haplotype indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Milosheska 2016 constructed ABCB1 1236C>T / 2677G>T/A / 3435C>T haplotypes and tested them on CL; no significant association was retained (Results paragraph 4). T-T-T haplotype frequency 53%, C-G-C 23%, C-G-T 14%."
    ),
    CONMED_BIRTHCONTROL = list(
      description = "Concomitant oral-contraceptive use indicator (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Observed as a trend of higher CL in patients on oral contraceptives (Results paragraph 3) but excluded from formal testing because only n = 2 patients used them; not in the final model."
    ),
    DOSE_LTG = list(
      description = "LTG daily dose",
      units       = "mg/day",
      type        = "continuous",
      notes       = "Entered forward-inclusion step but dropped during backward elimination (P > 0.01) per Results paragraph 3."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 1L,
    age_range      = "20.7 - 80.4 years (median 39.8)",
    age_median     = "39.8 years",
    weight_range   = "50 - 124 kg (median 70)",
    weight_median  = "70 kg",
    height_range   = "1.52 - 1.90 m (median 1.67)",
    bmi_range      = "18.7 - 41.3 kg/m^2 (median 24.5)",
    sex_female_pct = 70.0,
    race_ethnicity = c(Slovenian = 100),
    disease_state  = "Adult epilepsy on stable lamotrigine treatment (monotherapy in 54 patients or adjunctive therapy in 46 patients with carbamazepine, oxcarbazepine, valproic acid, phenytoin, phenobarbital, levetiracetam, topiramate, lacosamide, zonisamide, pregabalin, clonazepam, and/or clobazam) for at least 2 months; chronic renal / hepatic disease and pregnancy were exclusion criteria.",
    dose_range     = "50 - 600 mg/day oral lamotrigine (median 200 mg/day); dosing regimens QD, BID, or TID at steady state.",
    regions        = "Slovenia (single-centre; Department of Neurology, University Medical Centre Ljubljana).",
    co_medication  = "Inducers (Ind, retained pooled): carbamazepine n = 4, oxcarbazepine n = 6 (not counted as inducer, no CL association observed), phenytoin n = 2, phenobarbital n = 2, oral contraceptives n = 2 (excluded due to small N). Inhibitors (Inh, retained pooled): valproic acid n = 13, sertraline n = 2. Enzyme-inducer group defined by the paper as CBZ + PB + PHT (see Table 4 footnote).",
    n_observations = "195 lamotrigine and 195 lamotrigine-N-2-glucuronide plasma concentrations (up to two samples per patient, approximating trough and peak at steady state).",
    notes          = "Baseline demographics from Milosheska 2016 Table 1. Genotyping counts from Table 2 (94 - 99 patients genotyped depending on the SNP; the paper lists percentages, and this list captures the corresponding n via `round(pct * total)` in the accompanying vignette). Population priors from a meta-analysis of 11 prior LTG popPK studies (paper Methods 'Meta-analysis of previous population pharmacokinetic studies with LTG') were used via NONMEM PRIOR functionality to stabilise Ka and V estimation given the sparse (up to two samples per patient) design."
  )

  ini({
    # ----- Structural parameters -----
    # Parent lamotrigine (LTG) is fit to 195 plasma concentrations with
    # prior stabilisation from a meta-analysis of 11 prior LTG popPK
    # studies. The metabolite N-2-glucuronide (LTG-glu) is added in a
    # second-step joint model that assumes complete metabolic conversion.
    # Parent parameters below are taken from Milosheska 2016 Table 4
    # (the paper's parent-only final model); the paper's Discussion
    # states that the LTG parameters obtained in the joint parent +
    # metabolite fit "were very similar to those obtained with the
    # parent drug model" but does not tabulate them separately, so we
    # use Table 4 values for the LTG side (see vignette Assumptions).
    lka       <- log(1.96)   ; label("Absorption rate constant Ka (1/h)")                          # Milosheska 2016 Table 4 row 'Ka': 1.96 (95% CI 1.72 - 2.24)
    lcl       <- log(2.40)   ; label("Apparent parent oral clearance CL/F at 70 kg, non-smoker, no comed, CLcr = 110 (L/h)")  # Milosheska 2016 Table 4 row 'CL': 2.40 (95% CI 2.30 - 2.50)
    lvc       <- log(76.2)   ; label("Apparent parent volume of distribution V/F at 70 kg (L)")   # Milosheska 2016 Table 4 row 'V': 76.2 (95% CI 66.9 - 84.2)

    # Metabolite (LTG-N-2-glucuronide) apparent CL and V, from the
    # 'Population pharmacokinetic analysis of the lamotrigine-N-2-
    # glucuronide' section (the paper's joint-model narrative). Values
    # reported in prose only (no table). Because mass balance across
    # the parent -> metabolite arrow is written without an explicit
    # molecular-weight correction, CL_LTG-glu and V_LTG-glu are
    # "apparent" values that fold together the MW ratio and any
    # fraction-metabolised deviation from 1 (the paper assumes 100%
    # conversion; see Discussion 'renal elimination of unchanged LTG
    # is insignificant ... conversion to LTG-N-5-glucuronide is a minor
    # route').
    lcl_gluc  <- log(3.16)   ; label("Apparent metabolite clearance CL_LTG-glu / F_metab at 70 kg, CLcr = 110 (L/h)")  # Milosheska 2016 Results 'Population PK analysis of the lamotrigine-N-2-glucuronide' paragraph 1: 'Typical value of CL_LTG-glu was estimated at 3.16 l h-1'
    lvc_gluc  <- log(110)    ; label("Apparent metabolite volume of distribution V_LTG-glu / F_metab (L)")  # Milosheska 2016 Results, same paragraph: 'Distribution volume of the metabolite (V LTG-glu) was estimated at 110 l'

    # ----- Covariate effects on parent CL -----
    # Encoded uniformly as (1 + e * COV) so all coefficient signs match
    # Table 4. The two indicator groups written with a minus sign in
    # the Table 4 formula (Inh, UGT2B7 -161 CT / TT) carry negative
    # e_ coefficients.
    e_conmed_ugt_inh_cl <- -0.579   ; label("Effect of pooled UGT inhibitor (VPA or sertraline) on parent CL/F (multiplicative, unitless)")  # Milosheska 2016 Table 4 row 'Co-treatment with inhibitors*': -0.579 (95% CI -0.674, -0.483)
    e_wt_cl             <-  0.938   ; label("Allometric power exponent on parent CL/F (unitless)")                                            # Milosheska 2016 Table 4 row 'Body weight' on CL: 0.938 (95% CI 0.558, 1.32) -- the 95% CI includes the theoretical value 0.75
    e_smoke_cl          <-  0.340   ; label("Effect of current smoking on parent CL/F (multiplicative, unitless)")                            # Milosheska 2016 Table 4 row 'Cigarette smoking': 0.340 (95% CI 0.147, 0.590)
    e_conmed_eiaed_cl   <-  0.546   ; label("Effect of concomitant enzyme-inducing AED (CBZ / PB / PHT) on parent CL/F (multiplicative, unitless)")  # Milosheska 2016 Table 4 row 'Co-treatment with inducers': 0.546 (95% CI 0.114, 1.21)
    e_ugt2b7_m161ct_cl  <- -0.0358  ; label("Effect of UGT2B7 -161 CT vs CC on parent CL/F (multiplicative, unitless)")                       # Milosheska 2016 Table 4 row 'UGT2B7 - 161C>T genotype CT vs CC': -0.0358 (95% CI -0.161, 0.098)
    e_ugt2b7_m161tt_cl  <- -0.204   ; label("Effect of UGT2B7 -161 TT vs CC on parent CL/F (multiplicative, unitless)")                       # Milosheska 2016 Table 4 row 'UGT2B7 - 161C>T genotype TT vs CC': -0.204 (95% CI -0.336, -0.0364)
    e_ugt2b7_372ag_cl   <-  0.194   ; label("Effect of UGT2B7 372 AG vs AA on parent CL/F (multiplicative, unitless)")                        # Milosheska 2016 Table 4 row 'UGT2B7 372 A > G genotype AG vs AA': 0.194 (95% CI 0.0331, 0.392)
    e_ugt2b7_372gg_cl   <-  1.17    ; label("Effect of UGT2B7 372 GG vs AA on parent CL/F (multiplicative, unitless)")                        # Milosheska 2016 Table 4 row 'UGT2B7 372 A > G genotype GG vs AA': 1.17 (95% CI 0.448, 2.47)
    e_crcl_cl           <-  0.00328 ; label("Linear effect of CRCL - 110 on parent CL/F (multiplicative, per mL/min)")                        # Milosheska 2016 Table 4 row 'CL cr': 0.00328 (95% CI 0.000550, 0.00640)

    # ----- Covariate effects on parent V -----
    e_wt_vc <- 0.0181 ; label("Linear effect of WT - 70 on parent V/F (multiplicative, per kg)")   # Milosheska 2016 Table 4 row 'Body weight' on V: 0.0181 (95% CI 0.0102, 0.0239)

    # ----- Covariate effects on metabolite CL_gluc -----
    # The metabolite arm carries a power weight and linear CLcr effect
    # per the Results narrative in the 'Population pharmacokinetic
    # analysis of the lamotrigine-N-2-glucuronide' section.
    e_wt_cl_gluc   <- 1.01    ; label("Allometric power exponent on metabolite CL_LTG-glu / F_metab (unitless)")                # Milosheska 2016 Results paragraph 1: 'patient weight (power model, exponent 1.01)'
    e_crcl_cl_gluc <- 0.00759 ; label("Linear effect of CRCL - 110 on metabolite CL_LTG-glu / F_metab (multiplicative, per mL/min)")  # Milosheska 2016 Results paragraph 1: 'linear model, a decrease of 0.759% per 1 ml min-1 decrease in CLcr'

    # ----- Interindividual variability -----
    # Milosheska 2016 Table 4 reports IIV as CV%; the internal
    # log-normal variance is omega^2 = log(1 + CV^2).
    etalka       ~ 0.4090   # Table 4 row 'IIV Ka (%)': 71.1; log(1 + 0.711^2) = 0.4090
    etalcl       ~ 0.1040   # Table 4 row 'IIV CL (%)': 33.1; log(1 + 0.331^2) = 0.1040
    etalvc       ~ 0.0867   # Table 4 row 'IIV V (%)':  30.1; log(1 + 0.301^2) = 0.0867
    etalcl_gluc  ~ 0.1603   # Metabolite text: 'Unexplained IIV on CL_LTG-glu was 41.7%'; log(1 + 0.417^2) = 0.1603
    # No etalvc_gluc: the paper explicitly reports 'Due to sparse
    # concentration measurements and no prior data on pharmacokinetics
    # of the metabolite we were not able to estimate the IIV of
    # V LTG-glu' (Results paragraph 1 of the metabolite section).

    # ----- Residual error -----
    propSd      <- 0.180   ; label("Proportional residual error on parent LTG concentration (fraction)")            # Milosheska 2016 Table 4 row 'Residual variability (%)': 18.0 (95% CI 10.4, 23.5)
    propSd_gluc <- 0.138   ; label("Proportional residual error on metabolite LTG-glu concentration (fraction)")   # Milosheska 2016 Results, metabolite section paragraph 1: 'The residual (intra-individual) variability of LTG-glu concentration was 13.8%'
  })

  model({
    # ----- Reference values -----
    ref_wt   <- 70                 # Milosheska 2016 Table 4 formula uses (Wt/70) as the weight reference
    ref_crcl <- 110                # Milosheska 2016 Table 4 formula uses (CLcr - 110) as the CRCL reference

    # ----- Composite parent CL covariate multiplier -----
    # Reproduces the Table 4 footnote equation
    # CL = 2.40 * (1 - 0.579*Inh) * (Wt/70)^0.938 * (1 + 0.340*Tob) *
    #      (1 + 0.546*Ind) * (1 - 0.0358*UGT2B7_-161CT) *
    #      (1 - 0.204*UGT2B7_-161TT) * (1 + 0.194*UGT2B7_372AG) *
    #      (1 + 1.17*UGT2B7_372GG) * (1 + 0.00328*(CLcr - 110)).
    cov_cl <- (1 + e_conmed_ugt_inh_cl * CONMED_UGT_INH) *
              (WT / ref_wt)^e_wt_cl *
              (1 + e_smoke_cl        * SMOKE) *
              (1 + e_conmed_eiaed_cl * CONMED_EIAED) *
              (1 + e_ugt2b7_m161ct_cl * UGT2B7_M161CT) *
              (1 + e_ugt2b7_m161tt_cl * UGT2B7_M161TT) *
              (1 + e_ugt2b7_372ag_cl  * UGT2B7_372AG) *
              (1 + e_ugt2b7_372gg_cl  * UGT2B7_372GG) *
              (1 + e_crcl_cl          * (CRCL - ref_crcl))

    # ----- Composite parent V covariate multiplier -----
    # V(l) = 76.2 * (1 + 0.0181 * (Wt - 70)).
    cov_vc <- 1 + e_wt_vc * (WT - ref_wt)

    # ----- Composite metabolite CL covariate multiplier -----
    # CL_LTG-glu / F_metab = 3.16 * (Wt/70)^1.01 * (1 + 0.00759 * (CLcr - 110)).
    cov_cl_gluc <- (WT / ref_wt)^e_wt_cl_gluc *
                   (1 + e_crcl_cl_gluc * (CRCL - ref_crcl))

    # ----- Individual parameters -----
    ka       <- exp(lka       + etalka)
    cl       <- exp(lcl       + etalcl)      * cov_cl
    vc       <- exp(lvc       + etalvc)      * cov_vc
    cl_gluc  <- exp(lcl_gluc  + etalcl_gluc) * cov_cl_gluc
    vc_gluc  <- exp(lvc_gluc)                # no IIV per paper text

    # ----- Driving concentrations -----
    Cc      <- central      / vc
    Cc_gluc <- central_gluc / vc_gluc

    # ----- ODE system -----
    # Parent absorption + linear elimination. Metabolite formation rate
    # equals the parent's total elimination flux (mass per time) under
    # the paper's complete-conversion assumption; metabolite is then
    # eliminated linearly. Following the Abduljalil 2009 clarithromycin
    # + 14-OH-CLA precedent, this ODE is written in mass units without
    # an explicit MW correction; CL_LTG-glu and V_LTG-glu are already
    # "apparent" values that absorb the MW ratio and any deviation of
    # the true metabolic conversion fraction from 1.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot     - cl      * Cc
    d/dt(central_gluc) <-  cl  * Cc       - cl_gluc * Cc_gluc

    # ----- Observations and residual error -----
    Cc      ~ prop(propSd)
    Cc_gluc ~ prop(propSd_gluc)
  })
}
