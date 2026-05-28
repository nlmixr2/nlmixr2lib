Majekodunmi_2017_HIV_HCV_CD4_recovery <- function() {
  description <- paste0(
    "Longitudinal disease-progression / immune-reconstitution model for ",
    "age-standardised CD4 T-cell counts (z-scores) in HIV-infected ",
    "children receiving antiretroviral therapy (ART), with HIV/HCV ",
    "coinfection slowing the recovery rate (Majekodunmi 2017, fit to ",
    "401 children -- 355 HIV monoinfected and 46 HIV/HCV coinfected -- ",
    "from the Ukraine Paediatric HIV Cohort Study and the EPPICC ",
    "HIV/HCV coinfection study across 8 European countries). The age-",
    "standardised CD4 z score is modelled as an asymptotic recovery ",
    "curve z(t) = asy + (int - asy) * exp(-krec * t), where t is duration ",
    "on ART, int is the pre-ART z score, asy is the long-term z score, ",
    "and krec (the paper's symbol c) is the per-subject recovery rate ",
    "(1/year; ln(2)/krec is the time to half the total recovery from int ",
    "to asy; renamed from the paper's c to avoid shadowing R's built-in ",
    "c() combine function). Age at start of ART (centred at 4.3 years ",
    "-- the all-cohort median per Table 2 footnote) shifts both asy and ",
    "int (younger children start higher and reach higher long-term ",
    "levels). EPPICC enrollment country shifts int with Ukraine as the ",
    "implicit reference. HCV coinfection is a multiplicative fractional ",
    "reduction on krec: krec_coinf = krec_mono * (1 + e_hcv_pos_krec * ",
    "HCV_POS) with e_hcv_pos_krec = -0.77, so coinfected children ",
    "recover at krec = 1.55 * 0.23 = 0.357 /year (half-time ~2 years) ",
    "versus 1.55 /year for monoinfected (half-time ~0.45 year). Disease-",
    "progression model with no drug dosing -- the population were on ",
    "combination ART (most commonly lamivudine + zidovudine + lopinavir/",
    "ritonavir) and the model's t = 0 is the time of ART initiation."
  )
  reference <- paste(
    "Majekodunmi AO, Thorne C, Malyuta R, Volokha A, Callard RE,",
    "Klein NJ, Lewis J; on behalf of The European Paediatric HIV/HCV",
    "Co-infection Study group in the European Pregnancy and Paediatric",
    "HIV Cohort Collaboration and the Ukraine Paediatric HIV Cohort",
    "Study in EuroCoord. Modelling CD4 T Cell Recovery in Hepatitis C",
    "and HIV Co-infected Children Receiving Antiretroviral Therapy.",
    "Pediatr Infect Dis J. 2017 May;36(5):e123-e129.",
    "doi:10.1097/INF.0000000000001478.",
    sep = " "
  )
  vignette <- "Majekodunmi_2017_HIV_HCV_CD4_recovery"
  units <- list(
    time = "year",
    dosing = "n/a (disease-progression / immune-reconstitution model with no drug dosing; ART is implicit in the population)",
    concentration = "(age-standardised CD4 T-cell count z-score, unitless)"
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age at start of ART (years). Time-fixed per subject in this model (the same baseline value is used for every observation row); the covariate shifts both the pre-ART z score (int) and the long-term z score (asy) via linear-additive forms centred on the all-cohort median 4.3 years.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centring value 4.3 years from the Table 2 footnote ('The reference case is a Ukrainian child starting ART 4.3 years of age (median age in the dataset)'). Used in the model as (AGE - 4.3) so that AGE = 4.3 reproduces the typical-value asy = -1.07 and int = -2.42. Pre-ART z score decreases by 0.29 units per year older at ART start (Int:age = -0.29) and long-term z score decreases by 0.11 units per year older (Asy:age = -0.11). Source range: ART initiation between approximately 0 and 18 years (Table 1: median 4.40 y for monoinfected IQR 1.73-7.07; median 3.12 y for coinfected IQR 1.31-5.65).",
      source_name = "age at start of ART"
    ),
    HCV_POS = list(
      description = "HCV coinfection indicator. 1 = HIV/HCV coinfected (HCV antibody positive and/or two positive HCV RNA detections on separate visits >= 3 months apart per the Methods 'Definitions' section); 0 = HIV monoinfected. Time-fixed per subject. Used as a multiplicative fractional reduction on the recovery-rate constant krec: krec = (krec_pop + etakrec) * (1 + e_hcv_pos_krec * HCV_POS), with e_hcv_pos_krec = -0.77, so coinfected subjects recover at 23% of the monoinfected rate.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Paper 'Definitions' section: HCV antibody positive and/or >= 2 positive HCV RNA detections on separate visits >= 3 months apart. Subjects with known spontaneous HCV clearance were excluded. The covariate enters krec (Table 2 row 'C:Coinf'); pre-ART z score (int) and long-term z score (asy) were tested for HCV effect in stepwise selection but not retained (paper Results: 'no statistically significant effect of HCV coinfection on either pre-ART or long-term CD4 z scores').",
      source_name = "HCV status"
    ),
    REGION_POLAND = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled at the Polish EPPICC cohort site (Medical University Warsaw / Regional Hospital of Infectious Disease cohort), 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference (all REGION_* = 0 implies Ukraine). Additive shift on int of +0.44 (Table 2 row 'Int:Poland'). Cohort size in source: 1 of 46 coinfected (no monoinfected children from Poland).",
      source_name = "EPPICC cohort"
    ),
    REGION_RUSSIA = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled at the Russian EPPICC cohort site (Republican Hospital of Infectious Diseases, St Petersburg), 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference. Additive shift on int of +0.69 (Table 2 row 'Int:Russia'). Cohort size in source: 17 of 46 coinfected.",
      source_name = "EPPICC cohort"
    ),
    REGION_SWITZERLAND = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled at the Swiss Mother and Child HIV Cohort Study (MoCHiV), 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference. Additive shift on int of +0.02 (Table 2 row 'Int:Switzerland'; essentially zero, retained because all 8 country indicators were carried as a fixed-effect block per the EPPICC covariate analysis). Cohort size in source: 3 of 46 coinfected.",
      source_name = "EPPICC cohort"
    ),
    REGION_UK = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled in the UK Collaborative HIV Paediatric Study (CHIPS), 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference. Additive shift on int of -17.5 (Table 2 row 'Int: United Kingdom') -- the magnitude is implausibly large for a CD4 z score effect and is anchored on a UK cohort of only 2 of 46 coinfected children, likely a small-sample artifact; the value is reproduced verbatim per the published table. See the validation vignette's Assumptions and deviations section.",
      source_name = "EPPICC cohort"
    ),
    REGION_SPAIN = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled at the Spanish Paediatric HIV Network (CoRISpe; Madrid and Barcelona), 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference. Additive shift on int of +2.89 (Table 2 row 'Int:Spain'). Cohort size in source: 2 of 46 coinfected.",
      source_name = "EPPICC cohort"
    ),
    REGION_GERMANY = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled in the German Competence Network on HIV-infected Children, 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference. Additive shift on int of +0.34 (Table 2 row 'Int:Germany'). Cohort size in source: 1 of 46 coinfected.",
      source_name = "EPPICC cohort"
    ),
    REGION_ITALY = list(
      description = "EPPICC enrollment-country indicator: 1 = subject enrolled in the Italian Register for HIV-infection in Children, 0 = otherwise. Time-fixed per subject.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Encoded with REGION_UKRAINE as the implicit reference. Additive shift on int of -3.63 (Table 2 row 'Int:Italy'); like the UK estimate this is anchored on a very small subgroup (2 of 46 coinfected). Cohort size in source: 2 of 46 coinfected.",
      source_name = "EPPICC cohort"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 401L,
    n_studies = 2L,
    age_range = "0-18 years at ART initiation; whole-cohort median 4.3 y (Table 2 footnote); HIV monoinfected median 4.40 y (IQR 1.73-7.07); HIV/HCV coinfected median 3.12 y (IQR 1.31-5.65)",
    age_median = "4.3 years at ART initiation (all-cohort median, Table 2 footnote)",
    weight_range = "(not a covariate in this model; not tabulated in the source)",
    weight_median = "(not tabulated in the source)",
    sex_female_pct = 52.1,
    race_ethnicity = "(not tabulated in the source -- multi-country European pediatric cohort, predominantly White by demography)",
    disease_state = "HIV-infected children on combination antiretroviral therapy (ART); 46 of 401 (11.5%) also HCV coinfected (HCV antibody positive and/or >= 2 HCV RNA detections >= 3 months apart). 38% of monoinfected and 13% of coinfected had AIDS-defining diagnosis at ART start (Table 1). Most common ART regimen across both cohorts was lamivudine + zidovudine + kaletra (lopinavir/ritonavir); 33% of monoinfected and 49% of coinfected on a 3-drug ART regimen.",
    dose_range = "(not applicable; disease-progression model with no drug dosing -- ART is implicit in the population)",
    regions = "Ukraine (Ukraine Paediatric HIV Cohort, 355 monoinfected + 18 coinfected; all monoinfected children from Ukraine), and 7 other EPPICC sites in Europe for the HIV/HCV coinfection study: Russia (17), Switzerland (3), Spain (2), UK (2), Italy (2), Germany (1), Poland (1)",
    n_observations = "Median 5 CD4 measurements per child (IQR 3-7); monoinfected IQR 4-7; coinfected median 4 (IQR 3-6). Median follow-up on ART: 4.2 years overall (IQR 2.7-5.3), 4.1 y monoinfected (IQR 2.7-5.2), 5.1 y coinfected (IQR 3.1-5.6).",
    notes = "Two pooled cohorts: (a) Ukraine Paediatric HIV Cohort Study (HIV monoinfected children, established January 2011) and (b) EPPICC HIV/HCV coinfection sub-study (8 European countries including Ukraine). Children with known spontaneous HCV clearance were excluded. Children with fewer than 2 CD4 measurements were excluded. Sex breakdown: monoinfected 171/355 male (48.2%), 184/355 female (51.8%); coinfected 20/46 male (43.5%), 25/46 female (54.3%), 1/46 unknown. Sex was tested but not retained as a model covariate. HIV vertical transmission was the dominant route (98% monoinfected, 98% coinfected). Pre-ART HIV viral load and AIDS status were tested as covariates and not retained in the final model. Demographics from Table 1. Software: NONMEM for model fitting; R and Wolfram Mathematica for post-processing (Methods 'Software and Algorithms')."
  )

  ini({
    # ====================================================================
    # Asymptotic-recovery CD4 z-score model (Majekodunmi 2017 Eq. 1, Figure
    # 1 schematic). Final multivariate model parameter estimates from Table
    # 2 ('Multivariate model' column). The reference case is a Ukrainian
    # child starting ART at 4.3 years of age (Table 2 footnote).
    # ====================================================================

    # ---------- Asymptote (long-term CD4 z score) ----------
    asy <- -1.07
    label("Typical long-term CD4 z score after ART, for a Ukrainian child starting ART at 4.3 y (z-score units)")
    # Majekodunmi 2017 Table 2 Asy 'Estimate' = -1.07

    e_age_asy <- -0.11
    label("Linear-additive effect of age at ART start on asy, per year above 4.3 y (z-score units / year)")
    # Majekodunmi 2017 Table 2 Asy:age 'Estimate' = -0.11

    # ---------- Intercept (pre-ART CD4 z score) ----------
    # Renamed from the paper's symbol 'int' to avoid the C reserved word
    # int (rxode2 transpiles to C); the eta pairs as 'etaintercept'.
    intercept <- -2.42
    label("Typical pre-ART CD4 z score, for a Ukrainian child starting ART at 4.3 y (z-score units)")
    # Majekodunmi 2017 Table 2 Int 'Estimate' = -2.42

    e_age_intercept <- -0.29
    label("Linear-additive effect of age at ART start on int, per year above 4.3 y (z-score units / year)")
    # Majekodunmi 2017 Table 2 Int:age 'Estimate' = -0.29

    e_region_poland_intercept <- 0.44
    label("Additive shift on int for Polish EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int:Poland 'Estimate' = 0.44

    e_region_russia_intercept <- 0.69
    label("Additive shift on int for Russian EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int:Russia 'Estimate' = 0.69

    e_region_switzerland_intercept <- 0.02
    label("Additive shift on int for Swiss EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int:Switzerland 'Estimate' = 0.02

    e_region_uk_intercept <- -17.5
    label("Additive shift on int for UK EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int: United Kingdom 'Estimate' = -17.5 (n = 2; reproduced verbatim, magnitude likely a small-sample artifact -- see vignette Assumptions)

    e_region_spain_intercept <- 2.89
    label("Additive shift on int for Spanish EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int:Spain 'Estimate' = 2.89

    e_region_germany_intercept <- 0.34
    label("Additive shift on int for German EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int:Germany 'Estimate' = 0.34

    e_region_italy_intercept <- -3.63
    label("Additive shift on int for Italian EPPICC cohort (z-score units; Ukraine reference)")
    # Majekodunmi 2017 Table 2 Int:Italy 'Estimate' = -3.63

    # ---------- Recovery-rate constant krec (paper's c) ----------
    # NOT log-transformed: the paper reports c on the linear scale (Table 2
    # Estimate 1.55) with an additive 'Variance of REs' of 0.39, i.e.
    # additive normal IIV on the linear scale (not log-normal CV). Keeping
    # the linear-scale parameterisation reproduces the source faithfully;
    # see Assumptions and deviations in the vignette for the downstream
    # implication (individual krec can occasionally simulate negative
    # under additive normal IIV with SD 0.62; the typical-value trajectory
    # is unaffected). Renamed from the paper's symbol c to avoid shadowing
    # R's built-in c() combine function.
    krec <- 1.55
    label("Typical CD4 z-score recovery-rate constant for an HIV-monoinfected child (1/year)")
    # Majekodunmi 2017 Table 2 c 'Estimate' = 1.55 (ln(2)/1.55 = 0.45 year half-recovery time, matching the paper text '5 months (0.45 years)')

    e_hcv_pos_krec <- -0.77
    label("Fractional shift on krec for HCV coinfection (multiplicative; krec_coinf = krec_mono * (1 + e_hcv_pos_krec))")
    # Majekodunmi 2017 Table 2 C:Coinf 'Estimate' = -0.77; verifies against paper text 'reduced recovery rate of 0.357 per year' (1.55 * (1 - 0.77) = 0.357) and '2 years' coinfected half-recovery time (ln(2)/0.357 = 1.94 yr)

    # ---------- Inter-individual variability ----------
    # Variance of REs from Table 2 (column 'Variance of REs'). Encoded as
    # additive normal etas on the natural (z-score) scale for asy and int;
    # additive normal on the linear (1/year) scale for krec.
    etaasy  ~ 1.83  # Majekodunmi 2017 Table 2 Variance of REs (Asy) = 1.83
    etaintercept ~ 4.78  # Majekodunmi 2017 Table 2 Variance of REs (Int) = 4.78
    etakrec ~ 0.39  # Majekodunmi 2017 Table 2 Variance of REs (c)   = 0.39

    # ---------- Residual error ----------
    # The paper reports 'Residual error 1.43 +/- 0.16'. The 1.43 is the SD
    # (additive on the z-score scale): the dependent variable in the
    # NONMEM fit was 'CD4 z scores' (Methods 'Mixed-effects Modeling') and
    # the residual epsilon_ij in Eq. 1 sits on the z-score scale together
    # with z_ij, int_i, and asy_i; the magnitude (1.43 z-score units) is
    # commensurate with the observed CD4 z-score scatter in Figure 2.
    addSd <- 1.43
    label("Additive residual SD on CD4 z-score (z-score units)")
    # Majekodunmi 2017 Table 2 'Residual error' = 1.43
  })

  model({
    # ====================================================================
    # Algebraic asymptotic-recovery model (Eq. 1):
    #   z(t) = asy + (int - asy) * exp(-krec * t)
    # where t (= the rxode2 'time' variable) is duration on ART in years.
    # The model takes no dose events; t = 0 is the time of ART initiation
    # and observations are CD4 z scores at subsequent visits.
    # ====================================================================

    # ----- 1. Individual asy (long-term z score) -----
    # Linear-additive age effect centred at 4.3 years (Table 2 footnote
    # reference age). Additive eta on z-score scale.
    asy_indiv <- asy +
                 e_age_asy * (AGE - 4.3) +
                 etaasy

    # ----- 2. Individual int (pre-ART z score) -----
    # Linear-additive age effect centred at 4.3 years, plus EPPICC-country
    # shifts with Ukraine as the implicit reference (all REGION_* = 0
    # reproduces the Ukrainian intercept). Additive eta on z-score scale.
    intercept_indiv <- intercept +
                 e_age_intercept * (AGE - 4.3) +
                 e_region_poland_intercept      * REGION_POLAND      +
                 e_region_russia_intercept      * REGION_RUSSIA      +
                 e_region_switzerland_intercept * REGION_SWITZERLAND +
                 e_region_uk_intercept          * REGION_UK          +
                 e_region_spain_intercept       * REGION_SPAIN       +
                 e_region_germany_intercept     * REGION_GERMANY     +
                 e_region_italy_intercept       * REGION_ITALY       +
                 etaintercept

    # ----- 3. Individual krec (recovery-rate constant, 1/year) -----
    # Multiplicative fractional reduction for HCV coinfection: when
    # HCV_POS = 0, krec = krec + etakrec (typical-value 1.55); when
    # HCV_POS = 1, krec = (krec + etakrec) * (1 - 0.77) =
    # (krec + etakrec) * 0.23 (typical-value 0.357). Additive normal
    # IIV; for stochastic simulation the user may wish to truncate krec
    # at a small positive value (see vignette).
    krec_indiv <- (krec + etakrec) * (1 + e_hcv_pos_krec * HCV_POS)

    # ----- 4. CD4 z-score trajectory (single declared observation Cc) -----
    # Asymptotic recovery: starts at intercept_indiv at t = 0, trends to
    # asy_indiv as t -> Inf, at exponential rate krec_indiv. ln(2)/krec
    # is the time to half the total recovery. The output Cc here is the
    # CD4 z-score (the model's only observation endpoint); the name Cc
    # is the nlmixr2lib convention for the single output variable and
    # does NOT refer to a drug concentration in this disease-progression
    # model. Additive residual error on the z-score scale.
    Cc <- asy_indiv + (intercept_indiv - asy_indiv) * exp(-krec_indiv * time)

    Cc ~ add(addSd)
  })
}
