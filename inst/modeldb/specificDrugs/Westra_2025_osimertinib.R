Westra_2025_osimertinib <- function() {
  description <- "Joint population PK model for osimertinib and its active metabolite AZ5104 in advanced EGFR-mutation-positive non-small cell lung cancer (NSCLC) patients selected for low osimertinib trough exposure, with concomitant cobicistat used as a CYP3A4-inhibiting pharmacokinetic booster (Westra 2025, OSIBOOST trial NCT03858491). First-order oral absorption feeds an osimertinib (parent) compartment, from which a serial AZ5104 metabolite compartment is formed at a rate fixed to 25 percent of the parent elimination rate constant. Body weight is scaled allometrically a priori on parent and metabolite apparent clearance (exponent 0.75), on both apparent volumes (exponent 1), and on the absorption rate constant (exponent -0.25). Concomitant cobicistat is retained as a multiplicative factor of 0.704 on osimertinib apparent clearance, a 29.6 percent reduction. Between-subject variability is a correlated block on parent and metabolite clearance; residual error is proportional."
  reference <- "Westra N, Kruithof PD, Croes S, van Geel RMJM, Hendriks LEL, Touw DJ, Kosterink JGW, Stevens J, Oude Munnink TH, Mian P. Osimertinib Cost Minimization in Non-Small Cell Lung Cancer (NSCLC) Treatment: Hypothesis Generation for a Population Pharmacokinetic Approach for Equivalent Dose Optimization of Osimertinib in Combination with Cobicistat. J Clin Pharmacol. 2025;65(12):1687-1698. doi:10.1002/jcph.70085"
  vignette <- "Westra_2025_osimertinib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; reported in kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power effects applied a priori (not estimated) with reference weight 70 kg: exponent 0.75 on parent CL/F and on AZ5104 CL/F, exponent 1 on parent V/F and on AZ5104 V/F, and exponent -0.25 on ka. Westra 2025 Methods 'PopPK Model Development' and Table 2 footnote b; the exponents appear literally in the supplementary Part SI NONMEM control stream ($PK block).",
      source_name        = "BW"
    ),
    CONMED_COBICISTAT = list(
      description        = "Concomitant cobicistat 150 mg once-daily coadministration indicator (1 = boosted with cobicistat, 0 = osimertinib monotherapy).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (osimertinib monotherapy; the reference state is the pre-boost osimertinib 80 mg QD phase of the OSIBOOST trial).",
      notes              = "Multiplicative power-form effect on osimertinib apparent clearance only: cl is multiplied by 0.704^CONMED_COBICISTAT, i.e. a 29.6 percent reduction in CL/F when cobicistat is coadministered (Westra 2025 Results 'PopPK Model' and Table 2 row 'Effect of cobicistat on CL/F'). Not confounded with dose level: every OSIBOOST patient contributed both an unboosted and a boosted steady-state profile in a within-subject sequential design, so the coefficient is a pure drug-drug-interaction effect rather than a dose-normalisation term. The paper reports that adding the same covariate to AZ5104 CL/F and/or to relative bioavailability did not further improve the fit (P > .05), so AZ5104 clearance carries no cobicistat effect. Time-varying within a subject: cobicistat was added after the monotherapy sampling occasion and steady state was re-established over 21 days before the boosted samples were drawn.",
      source_name        = "COBI"
    )
  )

  # Covariates screened by Westra 2025 but not retained in the final model.
  # Documentation only - these are deliberately absent from model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Westra 2025 Methods 'PopPK Model Development' states the categorical covariate gender was tested on CL/F and V/F; it is absent from the final model of Table 2 and from the supplementary Part SI control stream $PK block, so it was not retained. The OSIBOOST cohort was 63.6 percent female (Table 1). Carried in the analysis dataset as $INPUT column SEX."
    ),
    ALB = list(
      description = "Baseline serum albumin concentration.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Westra 2025 Methods 'PopPK Model Development' states albumin was tested, after normalisation on the population median, as a power model on CL/F and V/F; it is absent from the final model of Table 2 and from the supplementary Part SI control stream $PK block, so it was not retained. Albumin IS retained on parent V/F in the companion model Westra_2025_osimertinib_brownbase.R, where it is inherited fixed from Brown 2017. Carried in the analysis dataset as $INPUT column ALB."
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the supplementary Part SI NONMEM
  # $MODEL block, which names the three compartments ABSORB, PARENT and
  # METABOLITE, and against the $PK scaling statements S2 = V1/1000 and
  # S3 = V2/1000, which place both observed concentrations in ug/L.
  compartmentData <- list(
    depot          = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "osimertinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_az5104 = list(analyte = "AZ5104", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 11,
    n_studies      = 1,
    age_median     = "69.0 years",
    weight_median  = "78.5 kg",
    height_median  = "166.0 cm",
    bmi_median     = "23.6 kg/m^2",
    bsa_median     = "1.9 m^2",
    sex_female_pct = 63.6,
    race_ethnicity = c(Caucasian = 100.0),
    disease_state  = "Advanced non-small cell lung cancer, adenocarcinoma histology in 100 percent of the cohort, ECOG/WHO performance status 0-1 in 100 percent. All patients were on established osimertinib treatment and were selected for a relatively low steady-state trough concentration (Cmin,ss at or below 195 ug/L) under osimertinib 80 mg once daily. Smoking status: never 36.4 percent, current 9.1 percent, former 54.5 percent. No healthy volunteers.",
    dose_range     = "Oral osimertinib 80 mg once daily (one patient received an increased dose of 160 mg once daily), first as monotherapy and then with oral cobicistat 150 mg once daily added for at least 21 days to reach steady state.",
    regions        = "Netherlands (Maastricht University Medical Centre and the Antoni van Leeuwenhoek hospital).",
    n_observations = "88 osimertinib and AZ5104 plasma concentrations in total from the 11 patients. Steady-state sampling at pre-dose, 0.5-1.5 h, 2.5-3.5 h and 7-8 h on each of the monotherapy and cobicistat-boosted occasions. Seven additional t = 24 h therapeutic-drug-monitoring observations were excluded because they were not sampled per protocol (supplementary Part SI $DATA IGNORE(INACT=1) comment).",
    notes          = "Baseline characteristics reproduced from Westra 2025 Table 1, OSIBOOST cohort column; the table reports medians for continuous characteristics and percentages for categorical ones, with no ranges. Data are from the OSIBOOST trial (NCT03858491, ethics approval METC19-013). The osimertinib metabolite AZ7550 was measured in the trial but excluded from this analysis. Selection for low osimertinib exposure is a deliberate enrichment and, as the paper's own Limitations note, may bias the cobicistat effect upward relative to a general osimertinib population, because high CYP3A4/A5 activity is one cause of low osimertinib exposure."
  )

  ini({
    # Structural parameters. Reference values are typical values at the
    # allometric reference body weight of 70 kg (Westra 2025 Table 2
    # footnote b). All values are the population estimates of Table 2,
    # which are identical to the final estimates printed in the
    # supplementary Part SI NONMEM $THETA block.

    # ka could not be identified from the sparse OSIBOOST data and was
    # fixed at the literature value of Brown 2017. Westra 2025 Results
    # 'PopPK Model' reports a sensitivity analysis over ka in 0.05-0.45
    # that did not improve the visual diagnostics, so the value remained
    # fixed.
    lka        <- fixed(log(0.24)); label("First-order oral absorption rate constant at 70 kg (1/h)")               # Table 2 row 'Ka (/h)' = 0.24 (fixed); Part SI $THETA 2 '(0.24) FIX'
    lcl        <- log(19.0);        label("Apparent osimertinib clearance, CL/F, at 70 kg without cobicistat (L/h)")  # Table 2 row 'CL/F osimertinib (L/h)' = 19.0 (RSE 8.5%); Part SI $THETA 3
    lvc        <- log(990);         label("Apparent osimertinib central volume, V/F, at 70 kg (L)")                   # Table 2 row 'V/F osimertinib (L)' = 990 (RSE 29%); Part SI $THETA 5
    lcl_az5104 <- log(47.3);        label("Apparent AZ5104 clearance, CL/F, at 70 kg (L/h)")                          # Table 2 row 'CL/F AZ5104 (L/h)' = 47.3 (RSE 7.8%); Part SI $THETA 4
    lvc_az5104 <- log(184);         label("Apparent AZ5104 central volume, V/F, at 70 kg (L)")                        # Table 2 row 'V/F AZ5104 (L)' = 184 (RSE 48.5%); Part SI $THETA 6

    # Allometric exponents. Applied a priori (Westra 2025 Methods:
    # 'Allometrically scaling was applied a priori ... with power
    # exponents of 1, 0.75, and -0.25'), therefore fixed, not estimated.
    e_wt_cl        <- fixed(0.75);  label("Allometric exponent for body weight on osimertinib CL/F (unitless)")  # Methods 'PopPK Model Development'; Table 2 footnote b; Part SI $PK 'CL1 = THETA(3) * ((BW/70)**0.75)'
    e_wt_vc        <- fixed(1);     label("Allometric exponent for body weight on osimertinib V/F (unitless)")   # Methods 'PopPK Model Development'; Table 2 footnote b; Part SI $PK 'V1 = THETA(5)*(BW/70)'
    e_wt_ka        <- fixed(-0.25); label("Allometric exponent for body weight on ka (unitless)")                # Methods 'PopPK Model Development'; Table 2 footnote b; Part SI $PK 'KA = THETA(2)*((BW/70)**(-0.25))'
    e_wt_cl_az5104 <- fixed(0.75);  label("Allometric exponent for body weight on AZ5104 CL/F (unitless)")       # Table 2 footnote b; Part SI $PK 'CL2 = THETA(4)*((BW/70)**0.75)'
    e_wt_vc_az5104 <- fixed(1);     label("Allometric exponent for body weight on AZ5104 V/F (unitless)")        # Table 2 footnote b; Part SI $PK 'V2 = THETA(6)*(BW/70)'

    # Cobicistat drug-drug-interaction effect. Enters as a power form on
    # the binary indicator, so cl is multiplied by this factor when
    # cobicistat is coadministered and is unchanged otherwise.
    e_cobi_cl <- 0.704; label("Multiplicative factor of concomitant cobicistat on osimertinib CL/F (unitless)")  # Table 2 row 'Effect of cobicistat on CL/F' = 0.704 (RSE 3.5%); Part SI $THETA 7

    # Between-subject variability. Westra 2025 Table 2 reports the
    # variances 0.0691 and 0.0598 with the covariance 0.0292 alongside
    # them, and footnote a gives the correlation coefficient as 0.45.
    # 0.0292 / sqrt(0.0691 * 0.0598) = 0.454 confirms that the two
    # diagonal entries are VARIANCES on the log scale and that 0.0292 is
    # the covariance; the values are reproduced verbatim in the
    # supplementary Part SI $OMEGA BLOCK(2). Note that the Results text
    # 'PopPK model evaluation and validation' calls these percentages
    # 'coefficient of variation (CV%) of the BSV on CL/F' (quoting them as
    # 21.5% and 27.5%, the second differing from the 27.2% printed in
    # Table 2) whereas the Table 2 column header labels the same two
    # numbers 'RSE%'; the covariance arithmetic above shows the table
    # header is the correct reading and that these percentages are
    # relative standard errors, not CVs. The
    # implied log-scale CVs are sqrt(exp(0.0691)-1) = 26.7% for
    # osimertinib CL/F and sqrt(exp(0.0598)-1) = 24.8% for AZ5104 CL/F.
    etalcl + etalcl_az5104 ~ c(0.0691, 0.0292, 0.0598)                       # Table 2 rows 'BSV on CL/F osimertinib' 0.0691, 'Covariance' 0.0292, 'BSV on CL/F AZ5104' 0.0598; Part SI $OMEGA BLOCK

    # Residual error. The supplementary Part SI $ERROR block defines a
    # single weight W = IPRED*THETA(1) with $SIGMA 1 FIX, i.e. a purely
    # proportional model whose standard deviation is THETA(1) = 0.178,
    # applied to BOTH observed analytes through the one shared theta.
    # The two nlmixr2 endpoints therefore carry the same numeric value.
    propSd        <- 0.178; label("Proportional residual error on osimertinib (fraction)")                                     # Table 2 row 'Proportional error' = 0.178 (RSE 10.6%); Part SI $THETA 1
    propSd_az5104 <- 0.178; label("Proportional residual error on AZ5104 (fraction; the same shared theta as osimertinib)")    # Table 2 row 'Proportional error' = 0.178; Part SI $ERROR uses one W for both compartments
  })

  model({
    # Fraction of the osimertinib elimination flux that appears as
    # AZ5104, fixed at 0.25 per Westra 2025 Methods ('a fixed
    # metabolization rate of 25% as demonstrated previously by Brown et
    # al.'). In the supplementary Part SI $PK this is K23 = K20 * 0.25,
    # and the $DES block removes only K20*A(2) from the parent while
    # adding K23*A(2) to the metabolite. The metabolite formation flux is
    # therefore a fraction OF the parent elimination flux and does not
    # add to it; this is reproduced literally below.
    fmet <- 0.25

    # Allometric reference body weight (Westra 2025 Table 2 footnote b).
    ref_wt <- 70

    # Individual parameters.
    ka <- exp(lka) * (WT / ref_wt)^e_wt_ka

    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^e_wt_cl *
          e_cobi_cl^CONMED_COBICISTAT

    vc <- exp(lvc) * (WT / ref_wt)^e_wt_vc

    cl_az5104 <- exp(lcl_az5104 + etalcl_az5104) *
                 (WT / ref_wt)^e_wt_cl_az5104

    vc_az5104 <- exp(lvc_az5104) * (WT / ref_wt)^e_wt_vc_az5104

    # Elimination / transfer micro-constants, named as in the Part SI $PK
    # block (K20 parent elimination, K30 metabolite elimination,
    # K23 parent-to-metabolite formation).
    k20 <- cl / vc
    k30 <- cl_az5104 / vc_az5104
    k23 <- k20 * fmet

    # ODE system, amounts in mg. Reproduces the Part SI $DES block
    # exactly. No molar correction is applied between parent and
    # metabolite: Westra 2025 works on a mass scale throughout and its
    # control stream carries no molecular-weight factor, so the AZ5104
    # state and concentration are expressed in osimertinib mass
    # equivalents.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - k20 * central
    d/dt(central_az5104) <-  k23 * central - k30 * central_az5104

    # Observations. The Part SI $PK scaling S2 = V1/1000 and S3 = V2/1000
    # convert an amount in mg divided by a volume in L into ug/L, which
    # is the unit Westra 2025 reports Cmax and Cmin in (Table 3) and the
    # unit of the 125-259 ug/L provisional therapeutic window.
    Cc        <- 1000 * central / vc
    Cc_az5104 <- 1000 * central_az5104 / vc_az5104

    Cc        ~ prop(propSd)
    Cc_az5104 ~ prop(propSd_az5104)
  })
}
