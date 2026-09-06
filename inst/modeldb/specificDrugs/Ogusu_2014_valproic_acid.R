Ogusu_2014_valproic_acid <- function() {
  description <- paste0(
    "Sequential population PK / PK-PD model for valproic acid (VPA) and ",
    "the probability of gamma-glutamyltransferase (gamma-GT) elevation in ",
    "Japanese patients with epilepsy on long-term therapy (Ogusu 2014; ",
    "n = 237 with 827 steady-state therapeutic-drug-monitoring ",
    "concentrations for the PK layer, n = 169 with 472 gamma-GT ",
    "measurements for the PK-PD layer, ages 2.2-52.2 years, 100-2600 ",
    "mg/day). PK is one-compartment first-order oral absorption with a ",
    "3.00 h absorption lag fixed by the authors; the daily VPA dose ",
    "enters both Vd/F and CL/F as a power term, and female sex, ",
    "carbamazepine, phenobarbital, phenytoin and clobazam coadministration ",
    "each scale CL/F multiplicatively (Ogusu 2014 equations 5-8, Table 2). ",
    "PD is a logistic regression for the probability that serum gamma-GT ",
    "exceeds its age- and sex-stratified upper limit of normal, with the ",
    "logit driven linearly by the individual steady-state daily AUC and ",
    "shifted by complication with intellectual disability and the SOD2 ",
    "Val16Ala (rs4880) Val/Val genotype (equation 9, Table 3). The PD ",
    "layer consumes the PK layer's AUC, so the two are one coupled model. ",
    "Fitted in NONMEM 7.2.0 (FOCE with ADVAN2 TRANS2 for PK; Laplacian ",
    "for the PK-PD layer). NOTE: the source article is under a PLOS ONE ",
    "Expression of Concern on funding-disclosure grounds only -- no ",
    "concern is raised about the data, methods or results, and it is not ",
    "a retraction; see the reference field and the vignette Errata."
  )
  reference <- paste(
    "Ogusu N, Saruwatari J, Nakashima H, Noai M, Nishimura M, Deguchi M,",
    "Oniki K, Yasui-Furukori N, Kaneko S, Ishitsu T, Nakagawa K.",
    "Impact of the superoxide dismutase 2 Val16Ala polymorphism on the",
    "relationship between valproic acid exposure and elevation of",
    "gamma-glutamyltransferase in patients with epilepsy: a population",
    "pharmacokinetic-pharmacodynamic analysis.",
    "PLoS One. 2014;9(11):e111066. doi:10.1371/journal.pone.0111066.",
    "EXPRESSION OF CONCERN: The PLOS ONE Editors. Expression of Concern:",
    "Impact of the Superoxide Dismutase 2 Val16Ala Polymorphism on the",
    "Relationship between Valproic Acid Exposure and Elevation of",
    "gamma-Glutamyltransferase in Patients with Epilepsy: A Population",
    "Pharmacokinetic-Pharmacodynamic Analysis.",
    "PLoS One. 2023;18(1):e0279162. doi:10.1371/journal.pone.0279162.",
    "The stated grounds are funding disclosure only (the study was",
    "part-funded by the Smoking Research Foundation, which had received",
    "tobacco-industry support, contrary to the journal's 2010 policy on",
    "funding from tobacco companies). The notice raises no concern about",
    "the data, methods or results and is not a retraction.",
    sep = " "
  )
  vignette <- "Ogusu_2014_valproic_acid"
  units <- list(
    time          = "h",
    dosing        = "mg (single administration of the divided daily dose)",
    concentration = "Cc in ug/mL; prob_ggt_elevation is a probability (0-1)"
  )

  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    DOSE_VPA_MGD = list(
      description        = "Patient's own total daily dose of valproic acid, summed across the daily administrations and NOT normalised by body weight. Drives Vd/F and CL/F as a power term in the PK layer (Ogusu 2014 equations 6 and 7) and the logit slope in the PD layer (equation 9).",
      units              = "mg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters every term as the ratio DOSE_VPA_MGD/1000, i.e. the printed equations normalise the daily dose to 1000 mg/day; that normalisation is written explicitly inside model() rather than folded into the coefficients. PK cohort mean 934.3 +/- 540.2 mg/day, range 100-2600; PK-PD cohort mean 903.8 +/- 502.7 mg/day, range 100-2600 (Ogusu 2014 Table 1). Must be strictly positive. Table 4 tabulates model predictions over 400-1200 mg/day.",
      source_name        = "Dose (daily VPA dose, mg/day)"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. Multiplies CL/F as a power-of-binary factor (Ogusu 2014 equation 7).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; equation 7 states 'female = 1, male = 0')",
      notes              = "PK cohort 100/237 women (42.2%); PK-PD cohort 67/169 women (39.6%) (Ogusu 2014 Table 1). The 0.917 coefficient makes female CL/F 8.3% lower than male, which Ogusu 2014 Discussion attributes to lower UDP-glucuronosyltransferase activity in women.",
      source_name        = "Gender"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator; 1 = carbamazepine co-administered, 0 = not. Multiplies CL/F as a power-of-binary factor (Ogusu 2014 equation 7).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = "Carbamazepine induces CYP2C9, CYP2C19 and the UGTs (Ogusu 2014 Discussion), raising VPA CL/F by 19%. 190/827 PK records (23.0%) and 94/472 PK-PD records (19.9%) were on carbamazepine (Table 1; the Table 1 co-administration percentages are per record, not per patient). Retained in the PK layer but deliberately EXCLUDED from the final PK-PD layer for collinearity with intellectual disability -- see the DIS_INTELLDIS_MODSEV notes.",
      source_name        = "CBZ"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital indicator; 1 = phenobarbital co-administered, 0 = not. Multiplies CL/F as a power-of-binary factor (Ogusu 2014 equation 7).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes              = "Enzyme-inducing AED; raises VPA CL/F by 12%. 73/827 PK records (8.8%) and 28/472 PK-PD records (5.9%) (Ogusu 2014 Table 1). Excluded from the final PK-PD layer for the same collinearity reason as carbamazepine.",
      source_name        = "PB"
    ),
    CONMED_PHT = list(
      description        = "Concomitant phenytoin indicator; 1 = phenytoin co-administered, 0 = not. Multiplies CL/F as a power-of-binary factor (Ogusu 2014 equation 7).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes              = "The strongest of the three induction effects; raises VPA CL/F by 43%. 88/827 PK records (10.6%) and 45/472 PK-PD records (9.5%) (Ogusu 2014 Table 1). Excluded from the final PK-PD layer for the same collinearity reason as carbamazepine.",
      source_name        = "PHT"
    ),
    CONMED_CLB = list(
      description        = "Concomitant clobazam indicator; 1 = clobazam co-administered, 0 = not. Multiplies CL/F as a power-of-binary factor (Ogusu 2014 equation 7).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant clobazam)",
      notes              = "Clobazam is NOT an enzyme inducer, and its coefficient has the opposite sign to the CBZ / PB / PHT terms estimated in the same equation: it LOWERS VPA CL/F by 9.4%. Ogusu 2014 Discussion cites a prior paediatric report of the same direction and states the mechanism 'remains unknown'. 128/827 PK records (15.5%) and 90/472 PK-PD records (19.0%) (Table 1).",
      source_name        = "CLB"
    ),
    DIS_INTELLDIS_MODSEV = list(
      description        = "Complication with moderate or severe intellectual disability; 1 = present, 0 = absent. Shifts the logit intercept of the gamma-GT elevation model additively (Ogusu 2014 equation 9).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no intellectual disability; equation 9 states 'intellectual disability = 1 if an intellectual disability was present, and was otherwise 0')",
      notes              = "The largest single effect in the PD model (+3.62 on the logit). 131/237 (55.3%) of the PK cohort and 97/169 (57.4%) of the PK-PD cohort (Ogusu 2014 Table 1). COLLINEARITY: Ogusu 2014 reports that intellectual disability was significantly associated with an increased number of co-administered CBZ, PB and PHT (P < 0.0001) and therefore removed those three indicators from the final PK-PD model to reduce multicollinearity, so this coefficient is partly a proxy for enzyme-inducing polytherapy. This is why the PD layer carries no CONMED_* terms although the PK layer does.",
      source_name        = "Intellectual disability"
    ),
    SNP_SOD2_RS4880_TT = list(
      description        = "SOD2 Val16Ala (rs4880) homozygous wild-type indicator; 1 = Val/Val, 0 = Val/Ala or Ala/Ala. Shifts the logit intercept of the gamma-GT elevation model additively (Ogusu 2014 equation 9).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Val/Ala heterozygotes pooled with Ala/Ala homozygous-variant carriers; equation 9 states 'SOD2 Val/Val genotype = 1, SOD2 Val/Ala or Ala/Ala genotype = 0')",
      notes              = "NOTE THE INVERTED POLARITY relative to the usual variant-carrier convention: the indicator flags the WILD-TYPE homozygote, which is the higher-risk stratum here. rs4880 is c.47T>C, so the T allele codes Val (wild-type) and C codes Ala (variant). Genotype frequencies in the 169-patient PK-PD cohort: Val/Val 77.6%, Val/Ala 20.7%, Ala/Ala 1.7%; 16Ala allele frequency 12.7%, consistent with Hardy-Weinberg (Ogusu 2014 Results). The two strata were pooled because only 2 patients were Ala/Ala -- a cell-count decision, not a genetic-model claim.",
      source_name        = "SOD2 Val/Val genotype"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at the observation.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a covariate on the PD SLOPE and found significant during forward inclusion, but removed from the final PK-PD model because age, body weight and daily VPA dose 'significantly correlated with each other (P < 0.05)' and the authors reduced the set to break the multicollinearity (Ogusu 2014 Results, PK-PD Model). No point estimate is printed anywhere on disk. Cohort: 17.2 +/- 8.3 years, range 2.2-52.2 (PK) and 18.0 +/- 7.8, range 3.0-52.2 (PK-PD); 93.3% were 30 years or younger."
    ),
    WT = list(
      description = "Body weight at the observation.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on the PD SLOPE (and, per Methods, on Vd/F in the PK layer) but not retained in either final model, for the same age/weight/dose collinearity reason as AGE. No point estimate is printed. Cohort: 48.8 +/- 20.9 kg, range 9.6-120.5 (PK) and 51.0 +/- 20.1, range 13.0-120.5 (PK-PD) (Ogusu 2014 Table 1). Its absence is a genuine feature of this model, not an omission: a paediatric-to-adult VPA model with no weight term is unusual and is noted in the vignette Errata."
    ),
    DUR_VPA_THERAPY = list(
      description = "Duration of valproic acid therapy up to the observation.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on the PD SLOPE, significant during forward inclusion but NOT significant during backward elimination, and therefore removed from the full covariate model (Ogusu 2014 Results, PK-PD Model; Table S2). No point estimate is printed. Mean follow-up 3.2 +/- 4.0 years (PK) and 6.6 +/- 5.1 years (PK-PD)."
    ),
    SNP_CYP2C9_STAR3 = list(
      description = "CYP2C9*3 allele carrier indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F in the PK forward-inclusion analysis and found not statistically significant (Ogusu 2014 Results, Table S1). CYP2C9*1/*1 93.7%, *1/*3 6.3%, *3/*3 0%. No point estimate exists."
    ),
    SNP_CYP2C19_PHENOTYPE = list(
      description = "CYP2C19 metabolizer phenotype derived from *2 and *3 (homozygous EM / heterozygous EM / PM).",
      units       = "(category)",
      type        = "categorical",
      notes       = "Tested on CL/F in the PK forward-inclusion analysis and found not statistically significant (Ogusu 2014 Results, Table S1). Homozygous EM 35.9%, heterozygous EM 47.2%, PM 16.9%. No point estimate exists."
    ),
    SNP_GSTM1_NULL = list(
      description = "GSTM1 null-genotype indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on the PD logit and NOT identified as a statistically significant covariate (Ogusu 2014 Table S2, Discussion). 56.8% null among the 169 PK-PD patients. No point estimate exists."
    ),
    SNP_GSTT1_NULL = list(
      description = "GSTT1 null-genotype indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on the PD logit and NOT identified as a statistically significant covariate (Ogusu 2014 Table S2, Discussion). 47.9% null among the 169 PK-PD patients. No point estimate exists."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 237L,
    n_studies      = 1L,
    n_observations = "827 steady-state VPA concentrations from 237 patients (PK layer) and 472 serum gamma-GT measurements from 169 patients (PK-PD layer)",
    age_range      = "PK 17.2 +/- 8.3 years, range 2.2-52.2; PK-PD 18.0 +/- 7.8 years, range 3.0-52.2; 93.3% were 30 years or younger (Ogusu 2014 Table 1, Figure S2)",
    weight_range   = "PK 48.8 +/- 20.9 kg, range 9.6-120.5; PK-PD 51.0 +/- 20.1 kg, range 13.0-120.5 (Ogusu 2014 Table 1)",
    sex_female_pct = 42.2,
    race_ethnicity = c(Japanese = 100),
    disease_state  = "epilepsy on long-term valproic acid therapy; seizure locus generalized 46.8% / partial 50.2% / unidentified 3.0%; seizure aetiology idiopathic 27.0% / symptomatic 31.7% / cryptogenic 41.3%; 55.3% had a comorbid moderate or severe intellectual disability",
    dose_range     = "valproic acid 100-2600 mg/day orally (mean 934.3 +/- 540.2 mg/day), given in 1-3 divided doses; 45.7% of PK records were monotherapy",
    regions        = "Japan (single-centre retrospective therapeutic-drug-monitoring cohort)",
    notes          = paste0(
      "Retrospective therapeutic-drug-monitoring data, sparsely sampled: ",
      "the interval between the last dose and the sample was distributed ",
      "over the full 24 h. Standard errors could not be determined for the ",
      "final PK model (Ogusu 2014 Discussion, limitation 1), which is why ",
      "Table 2 reports no RSE column and precision comes from a 1000-run ",
      "stratified nonparametric bootstrap (965 successful minimisations ",
      "for PK, 997 for PK-PD). gamma-GT elevation is defined against the ",
      "age- and sex-stratified upper limit of normal; gamma-GT values from ",
      "the first 6 months of VPA therapy were excluded to avoid transient ",
      "post-initiation elevation. Baseline gamma-GT 51.3 +/- 68.6 IU/L ",
      "(PK) and 48.3 +/- 65.1 IU/L (PK-PD)."
    )
  )

  ini({
    # ==================================================================
    # PK LAYER -- Ogusu 2014 equations 5-8 and Table 2 ("NONMEM Final
    # Estimates" column). No standard errors were estimable for the final
    # PK model (Discussion, limitation 1); the precision quoted in the
    # trailing comments is the bootstrap 95% CI from Table 2.
    #
    # UNITS: Dose is the daily VPA dose in mg/day and enters every PK term
    # as (Dose/1000); CL/F is L/h, Vd/F is L, Ka is 1/h, ALAG is h.
    # ==================================================================

    lka   <- log(0.109) ; label("Absorption rate constant Ka (1/h)")                       # Ogusu 2014 eq 5 / Table 2: 0.109 (bootstrap median 0.104, 95% CI 0.0317-0.748)
    lvc   <- log(21.4)  ; label("Apparent volume of distribution Vd/F at a 1000 mg/day daily dose (L)")  # Ogusu 2014 eq 6 / Table 2: 21.4 (bootstrap median 21.2, 95% CI 8.01-68.1)
    lcl   <- log(0.559) ; label("Apparent oral clearance CL/F at a 1000 mg/day daily dose in an untreated male (L/h)")  # Ogusu 2014 eq 7 / Table 2: 0.559 (bootstrap median 0.558, 95% CI 0.520-0.589)

    # ALAG is reported as "3.00 (Fixed)" in the Table 2 NONMEM column, with
    # no bootstrap median and no CI -- the fixed() wrapper records that.
    ltlag <- fixed(log(3.00)) ; label("Absorption lag time ALAG (h)")                      # Ogusu 2014 eq 8 / Table 2: "3.00 (Fixed)"

    # ----- Covariate effects on the PK parameters (eq 6, eq 7) --------
    # Daily dose enters as a POWER of (Dose/1000) on both Vd/F and CL/F.
    e_dose_vc  <- 1.52  ; label("Power exponent of the normalised daily VPA dose (Dose/1000) on Vd/F (unitless)")  # Ogusu 2014 eq 6 / Table 2 "Dose on Vd/F": 1.52 (bootstrap median 1.46, 95% CI 0.457-3.74)
    e_dose_cl  <- 0.596 ; label("Power exponent of the normalised daily VPA dose (Dose/1000) on CL/F (unitless)")  # Ogusu 2014 eq 7 / Table 2 "Dose on CL/F": 0.596 (bootstrap median 0.592, 95% CI 0.525-0.652)

    # The five categorical CL/F effects are POWER-OF-BINARY factors: the
    # printed equation writes each coefficient raised to the 0/1 indicator
    # (e.g. 0.917^female), so the coefficient is the multiplicative ratio
    # in the indicated group relative to its reference.
    e_sexf_cl <- 0.917 ; label("Ratio of CL/F in women relative to men (unitless)")                             # Ogusu 2014 eq 7 / Table 2 "Gender on CL/F": 0.917 (bootstrap median 0.915, 95% CI 0.847-0.998)
    e_cbz_cl  <- 1.19  ; label("Ratio of CL/F with concomitant carbamazepine relative to without (unitless)")   # Ogusu 2014 eq 7 / Table 2 "CBZ on CL/F": 1.19 (bootstrap median 1.20, 95% CI 1.13-1.29)
    e_pb_cl   <- 1.12  ; label("Ratio of CL/F with concomitant phenobarbital relative to without (unitless)")   # Ogusu 2014 eq 7 / Table 2 "PB on CL/F": 1.12 (bootstrap median 1.12, 95% CI 1.03-1.24)
    e_pht_cl  <- 1.43  ; label("Ratio of CL/F with concomitant phenytoin relative to without (unitless)")       # Ogusu 2014 eq 7 / Table 2 "PHT on CL/F": 1.43 (bootstrap median 1.43, 95% CI 1.31-1.56)
    e_clb_cl  <- 0.906 ; label("Ratio of CL/F with concomitant clobazam relative to without (unitless)")        # Ogusu 2014 eq 7 / Table 2 "CLB on CL/F": 0.906 (bootstrap median 0.907, 95% CI 0.852-0.970)

    # ==================================================================
    # PD LAYER -- Ogusu 2014 equation 9 and Table 3 ("NONMEM Final
    # Estimates (RSE, %)" column).
    #
    # FORM. Equation 9 typesets its two categorical indicators as
    # SUPERSCRIPTS ("3.62^{Intellectual disability}"), which is a PLOS
    # display-equation rendering artifact, not an exponent. Equation 4 --
    # the paper's own statement of how a PD covariate enters -- is
    # ADDITIVE ("PD parameter = theta_p + theta_cov * covariate"), and
    # only the additive reading reproduces Table 4 (see the vignette's
    # 36-anchor gate; the exponent reading cannot even reproduce the
    # reference group, missing it by 11.3 percentage points on average).
    #
    # SIGN. Table 3 prints the intercept as "26.63" and its bootstrap
    # median as "26.68" with a 95% CI of "211.6-24.82". Those are
    # -6.63, -6.68 and (-11.6, -4.82): the PDF renders a minus sign as
    # the digit 2. A positive intercept is impossible here -- it would
    # put the reference-group probability above 99.8% against Table 4's
    # 5.2-14.8%.
    # ==================================================================

    logit_ref <- -6.63 ; label("Logit of the gamma-GT elevation probability at zero VPA exposure in the reference group (no intellectual disability, SOD2 Val/Ala or Ala/Ala) (unitless logit)")  # Ogusu 2014 eq 9 / Table 3 "Base": -6.63 (RSE 17.5%; bootstrap median -6.68, 95% CI -11.6 to -4.82)

    e_intelldis_logit <- 3.62 ; label("Additive shift in the gamma-GT elevation logit for complication with moderate or severe intellectual disability, versus no intellectual disability (unitless logit)")  # Ogusu 2014 eq 9 / Table 3 "Intellectual disability on BASE": 3.62 (RSE 28.4%; bootstrap median 3.62, 95% CI 1.87-8.94)
    e_sod2_vv_logit   <- 1.96 ; label("Additive shift in the gamma-GT elevation logit for the SOD2 rs4880 Val/Val genotype, versus Val/Ala or Ala/Ala (unitless logit)")  # Ogusu 2014 eq 9 / Table 3 "SOD2 genotype on BASE": 1.96 (RSE 44.1%; bootstrap median 2.02, 95% CI 0.33-4.42)

    # Table 3 lists ONLY "Dose on SLOPE" -- there is no separate base
    # SLOPE parameter anywhere in the paper, so the exposure slope IS
    # (Dose/1000)^1.55 with no multiplier. That is dimensionally
    # consistent only if the AUC entering equation 9 is expressed in
    # mg*h/mL (= g*h/L), i.e. the daily AUC divided by 1000; the paper
    # never states the unit. See the vignette's 36-anchor gate, where the
    # ug*h/mL reading saturates every cell at 100%.
    e_dose_slope <- 1.55 ; label("Power exponent of the normalised daily VPA dose (Dose/1000) giving the slope of the gamma-GT elevation logit on the steady-state daily AUC (per mg*h/mL)")  # Ogusu 2014 eq 9 / Table 3 "Dose on SLOPE": 1.55 (RSE 19.4%; bootstrap median 1.56, 95% CI 0.83-2.29)

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY
    #
    # SCALE. Table 2 and Table 3 each print a "NONMEM Final Estimates"
    # column and a "Bootstrap Median" column that look irreconcilable on
    # every variance row. They are the same quantities on different
    # scales: the NONMEM column is the VARIANCE and the bootstrap column
    # is the SD. sqrt() of the NONMEM value reproduces the bootstrap
    # median to three significant figures on five of the six rows:
    #   ALAG  sqrt(4.48e-9) = 6.69e-5 vs 6.69e-5
    #   Ka    sqrt(7.77e-7) = 8.81e-4 vs 8.82e-4
    #   Vd/F  sqrt(1.83e-7) = 4.28e-4 vs 4.28e-4
    #   CL/F  sqrt(0.0587)  = 0.242   vs 0.241
    #   sigma sqrt(0.0617)  = 0.248   vs 0.247
    #   logit sqrt(12.3)    = 3.507   vs 3.48
    # nlmixr2 takes variances here, so the NONMEM column is what is
    # encoded. The logit variance in particular is load-bearing: reading
    # 3.48 as the variance instead misses Table 4 by 6.5 percentage
    # points on average (see the vignette gate).
    #
    # The Ka, Vd/F and ALAG variances are effectively zero (1e-7 to 1e-9);
    # they are kept at the published values rather than dropped so the
    # file is a faithful transcription. In practice only CL/F and the
    # logit carry meaningful between-subject spread.
    # ==================================================================

    # Single quotes only in these trailing comments: an ini() line with no
    # label() whose comment contains a double quote is promoted into a label
    # and breaks the parse (tests/testthat/test-checkModelConventions.R).
    etaltlag ~ 4.48e-9  # Ogusu 2014 Table 2 'v2 on ALAG' (NONMEM column, a variance); bootstrap SD 6.69e-5
    etalka   ~ 7.77e-7  # Ogusu 2014 Table 2 'v2 on Ka'   (NONMEM column, a variance); bootstrap SD 8.82e-4
    etalvc   ~ 1.83e-7  # Ogusu 2014 Table 2 'v2 on Vd/F' (NONMEM column, a variance); bootstrap SD 4.28e-4
    etalcl   ~ 0.0587   # Ogusu 2014 Table 2 'v2 on CL/F' (NONMEM column, a variance); bootstrap SD 0.241, 95% CI 0.202-0.281
    etalogit_ref ~ 12.3 # Ogusu 2014 Table 3 'v2 on logit (Pr)' (NONMEM column, a variance; SD 3.51); RSE 43.4%, bootstrap SD 3.48, 95% CI 2.29-10.7

    # ==================================================================
    # RESIDUAL ERROR
    #
    # The paper contradicts itself on the form. Methods: variability "was
    # best described using a proportional error model". Results
    # (Population PK Model): "The best residual error model was an
    # additive model." Table 2's own row label settles it -- "s2
    # (proportional error)" -- and a physical check agrees: 0.248 read as
    # an ADDITIVE SD would be 0.248 ug/mL against a cohort mean VPA
    # concentration of 68.15 ug/mL, i.e. 0.36%, which is impossible for a
    # routine immunoassay; 24.8% proportional is ordinary for sparse TDM.
    # Encoded as proportional; the contradiction is in the vignette Errata.
    # ==================================================================

    propSd <- 0.248 ; label("Proportional residual error SD for VPA concentration (fraction)")  # Ogusu 2014 Table 2 "s2 (proportional error)" = 0.0617 (a variance); SD = sqrt(0.0617) = 0.248, matching the bootstrap median 0.247 (95% CI 0.224-0.268)

    # The PD endpoint is a Bernoulli outcome (gamma-GT elevated: yes/no),
    # so the source model has no residual-error parameter for it -- the
    # likelihood IS the logistic probability. This placeholder exists only
    # so prob_ggt_elevation can be emitted as an observation in rxode2;
    # it is not a published value. Same device as
    # Fukae_2024_valemetostat_anc_decrease.R.
    addSd_prob_ggt_elevation <- fixed(0.001) ; label("Placeholder additive residual SD on the gamma-GT elevation probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ---- Normalised daily dose ---------------------------------------
    # Every dose term in equations 6, 7 and 9 is written as (Dose/1000)
    # with Dose the DAILY dose in mg/day, not the amount in an individual
    # administration. Computed once here.
    dose_ratio <- DOSE_VPA_MGD / 1000

    # ---- Individual PK parameters (eq 5-8) ---------------------------
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    vc   <- exp(lvc + etalvc) * dose_ratio^e_dose_vc
    cl   <- exp(lcl + etalcl) * dose_ratio^e_dose_cl *
      e_sexf_cl^SEXF *
      e_cbz_cl^CONMED_CBZ *
      e_pb_cl^CONMED_PB *
      e_pht_cl^CONMED_PHT *
      e_clb_cl^CONMED_CLB

    # ---- One-compartment oral disposition (NONMEM ADVAN2 TRANS2) -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl / vc * central

    alag(depot) <- tlag

    Cc <- central / vc

    # ---- Steady-state daily exposure driving the PD layer ------------
    # Equation 9's AUC is "the individual AUC value of VPA that was
    # simulated based on the population PK analysis". At steady state the
    # AUC over one day equals the daily dose divided by CL/F, exactly --
    # no integration is required and none is implied by the paper.
    #
    # With DOSE_VPA_MGD in mg/day and cl in L/h this is mg*h/L = ug*h/mL.
    auc_ss <- DOSE_VPA_MGD / cl

    # Equation 9 needs the AUC in mg*h/mL (see the e_dose_slope comment in
    # ini()); the /1000 conversion is written out here rather than folded
    # into the coefficient so both quantities stay inspectable.
    auc_ss_mghml <- auc_ss / 1000

    # ---- Logistic model for gamma-GT elevation (eq 9) ----------------
    slope_ggt <- dose_ratio^e_dose_slope

    logit_ggt_elevation <- logit_ref +
      e_intelldis_logit * DIS_INTELLDIS_MODSEV +
      e_sod2_vv_logit * SNP_SOD2_RS4880_TT +
      slope_ggt * auc_ss_mghml +
      etalogit_ref

    prob_ggt_elevation <- expit(logit_ggt_elevation)

    # ---- Observations -------------------------------------------------
    Cc ~ prop(propSd)
    prob_ggt_elevation ~ add(addSd_prob_ggt_elevation)
  })
}
