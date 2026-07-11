Oishi_2014_fesoterodine <- function() {
  description <- "One-compartment population PK model with first-order absorption and a lag time for the active metabolite 5-hydroxymethyl tolterodine (5-HMT) following oral administration of the prodrug fesoterodine as a sustained-release tablet (Oishi 2014). Pooled analysis of 10 pharmacokinetic and 3 efficacy/safety studies in Western and East Asian healthy volunteers and adult patients with overactive bladder. Dose entered as mg of fesoterodine fumarate; CL/F and V/F are reported by the paper as molecular-weight-adjusted values for 5-HMT so that dose (mg fumarate) can be used directly with the reported parameters to yield 5-HMT concentrations. CL/F varies with Cockcroft-Gault creatinine clearance (power exponent 0.303 centred at 80 mL/min), moderate hepatic impairment (Child-Pugh Class B multiplier 0.422), CYP2D6 poor-metaboliser phenotype (multiplier 0.626), concomitant CYP3A inhibitor coadministration (multiplier 0.504, source population coadministered ketoconazole), concomitant CYP3A inducer coadministration (multiplier 3.90, source population coadministered rifampicin), female sex (multiplier 0.903), and Japanese ethnicity (multiplier 1.12). Inter-individual variability is a full block on CL/F, V/F, and ka; Tlag IIV was fixed to zero by the source authors. Residual variability is a combined additive + proportional error model on 5-HMT plasma concentrations (LOQ 0.02 ng/mL for most studies)."
  reference <- "Oishi M, Tomono Y, Yamagami H, Malhotra B. Population Pharmacokinetics of the 5-Hydroxymethyl Metabolite of Tolterodine After Administration of Fesoterodine Sustained Release Tablet in Western and East Asian Populations. Journal of Clinical Pharmacology. 2014;54(8):928-936. doi:10.1002/jcph.274"
  vignette <- "Oishi_2014_fesoterodine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalised)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts both BSA-normalised and raw Cockcroft-Gault values; per-model units are documented here). Precedent for using raw Cockcroft-Gault in mL/min under CRCL: Delattre 2010 amikacin. Oishi 2014 Methods 'Covariate Model' derived CRCL from the Cockcroft-Gault equation using serum creatinine, age, body weight, and sex; for study #3 CRCL was measured from serum and urine creatinine over a 24 h urine-collection interval. Reference value is 80 mL/min (Oishi 2014 final equation and typical-subject definition on p. 931). Continuous power effect on CL/F: (CRCL/80)^0.303. Cohort mean (SD) 87.0 (30.3) mL/min, range 19.8-273.1 (Table 2 (b)).",
      source_name        = "CLCR"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = Child-Pugh Class B, 0 = healthy or non-moderate)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy or non-moderate hepatic function)",
      notes              = "Classification scheme: Child-Pugh Class B (composite score 7-9). Only Child-Pugh B patients from Oishi 2014 study #4 were flagged as HEPIMP_MOD = 1; severe hepatic impairment (Child-Pugh C) subjects were excluded from all contributing studies (Oishi 2014 Methods 'Covariate Model'). Multiplicative effect on CL/F: exp(e_hepimp_mod_cl) = 0.422 (Table 3 theta_6). Cohort count: 8 subjects (0.5%) had Child-Pugh B, all from study #4 (Table 2 (a)).",
      source_name        = "HEP"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metaboliser phenotype indicator (1 = poor metaboliser, 0 = extensive metaboliser)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive metaboliser)",
      notes              = "Phenotype assigned from genotype: homozygous non-functional CYP2D6*3 or CYP2D6*4 -> poor metaboliser (PM); all other genotypes -> extensive metaboliser (EM). Study #10 (Japanese single-dose healthy male) and study #13 (Phase 2 Asian OAB) did not genotype CYP2D6; per Oishi 2014 Methods 'Covariate Model' those subjects were classified as EM based on the reported ~1-2% PM frequency in East Asian populations. Multiplicative effect on CL/F: exp(e_cyp2d6_pm_cl) = 0.626 (Table 3 theta_7). Cohort count: 71 PM (4.6%), 1475 EM (95.4%) (Table 2 (a)).",
      source_name        = "CYP2D6 genotype"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A4 inhibitor indicator (1 = coadministered with a CYP3A4 inhibitor, 0 = no inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor coadministration)",
      notes              = "Pooled indicator captures ketoconazole coadministration only in the Oishi 2014 cohort (Methods 'Study Design': studies #1 and #6 investigated 200 mg ketoconazole QD or BID for 6 days with a fesoterodine single dose). Multiplicative effect on CL/F: exp(e_cyp3a4_inh_cl) = 0.504 (Table 3 theta_8). Cohort count: 118 subjects (7.4%) coadministered a CYP3A4 inhibitor (Table 2 (a)).",
      source_name        = "CYP3A inhibitor"
    ),
    CONMED_CYP3A4_IND = list(
      description        = "Concomitant CYP3A4 inducer indicator (1 = coadministered with a CYP3A4 inducer, 0 = no inducer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inducer coadministration)",
      notes              = "Pooled indicator captures rifampicin coadministration only in the Oishi 2014 cohort (Methods 'Study Design': study #5 investigated 600 mg rifampicin QD for 8 days with a fesoterodine single dose). Multiplicative effect on CL/F: exp(e_cyp3a4_ind_cl) = 3.90 (Table 3 theta_9). Cohort count: 12 subjects (0.8%) coadministered a CYP3A4 inducer (Table 2 (a)).",
      source_name        = "CYP3A inducer"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effect on CL/F: exp(e_sexf_cl) = 0.903, i.e. females have ~10% lower typical CL/F than males (Table 3 theta_11; Discussion 'about 10% decrease in females'). The magnitude was considered not clinically significant by the authors (bootstrap 90% CI within 0.8-1.25 bioequivalence range). Cohort split: 1085 female (70.2%), 461 male (29.8%) (Table 2 (a)).",
      source_name        = "Sex"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese ethnicity indicator (1 = Japanese, 0 = non-Japanese)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese; the modelled cohort combines Westerners, Koreans, and other Asians into the reference)",
      notes              = "Multiplicative effect on CL/F: exp(e_race_japanese_cl) = 1.12, i.e. Japanese subjects have ~12% higher typical CL/F than non-Japanese (Table 3 theta_13; Discussion 'about 10% increase in Japanese'). The magnitude was considered not clinically significant by the authors (bootstrap 90% CI within 0.8-1.25 bioequivalence range). Cohort composition (Table 2 (a)): Westerners 878 (56.8%), Japanese 522 (33.8%), Korean 105 (6.8%), Other Asians 41 (2.7%).",
      source_name        = "Japanese"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry",
      units       = "years",
      type        = "continuous",
      notes       = "Oishi 2014 Methods 'Covariate Model' tested age on CL/F but it was not retained in the final model (not selected in forward inclusion). Cohort mean (SD) 55.5 (14.6) years, range 19-91 (Table 2 (b))."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Oishi 2014 Methods 'Covariate Model' tested body weight on both CL/F and V/F. Not retained in the final model on either parameter. Cohort mean (SD) 70.2 (18.8) kg, range 31.2-192.8 (Table 2 (b))."
    ),
    RACE_KOREAN = list(
      description = "Korean ethnicity indicator (1 = Korean, 0 = non-Korean)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Oishi 2014 Results 'Population Pharmacokinetic Modeling Results' selected Korean on V/F in forward analysis but excluded it in backward analysis based on Delta-OFV < 10.8. Cohort count 105 Korean subjects (6.8%) (Table 2 (a))."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1546L,
    n_studies      = 13L,
    n_observations = 10922L,
    age_range      = "19-91 years",
    age_median     = "56 years",
    weight_range   = "31.2-192.8 kg",
    weight_median  = "67.0 kg",
    sex_female_pct = 70.2,
    race_ethnicity = c(Westerner = 56.8, Japanese = 33.8, Korean = 6.8, `Other Asian` = 2.7),
    disease_state  = "Pooled cohort of adult healthy volunteers and patients with overactive bladder (OAB). Contributing studies included dedicated renal-impairment (study #3), hepatic-impairment (study #4), and CYP3A inhibitor/inducer drug-drug interaction (studies #1, #5, #6) pharmacokinetic sub-studies alongside Phase 2 and Phase 3 OAB efficacy/safety studies (see Table 1).",
    dose_range     = "4-28 mg fesoterodine (as fesoterodine fumarate; molecular weight 527.65 g/mol) once daily as sustained-release oral tablets; single-dose and multiple-dose regimens. Clinically recommended doses are 4 and 8 mg QD.",
    regions        = "Multi-regional pool: North America and Europe (7 Western PK studies + 2 Western OAB efficacy/safety studies); Japan (2 Japanese healthy-volunteer PK studies); Korea (1 Korean healthy-volunteer PK study); Japan / Korea / Taiwan / Hong Kong (1 Asian OAB Phase 2 study).",
    renal_function = "Cockcroft-Gault CRCL 19.8-273.1 mL/min (mean 87.0, SD 30.3); one dedicated renal-impairment study (#3) with 24 patients across the renal-function spectrum.",
    hepatic_function = "Only Child-Pugh Class B moderate hepatic impairment enrolled (n = 8, all from study #4); severe hepatic impairment excluded from all contributing studies.",
    cyp2d6_distribution = "Extensive metaboliser 1475 (95.4%); Poor metaboliser 71 (4.6%) (Table 2 (a)).",
    conmed_distribution = "CYP3A inhibitor coadministration 118 (7.4%); CYP3A inducer coadministration 12 (0.8%) (Table 2 (a)).",
    notes          = "NONMEM version V, level 1.1 with FOCEI. Structural model ADVAN2 / TRANS2 (1-cpt, first-order absorption). Preliminary 2-compartment model was unstable with the sparse Phase 2/3 sampling and abandoned in favour of the 1-cpt structural form. 5-HMT concentrations below the analytical LOQ (0.02 ng/mL for most studies) were excluded. Predictive check by VPC (1000 Monte Carlo simulations) and nonparametric bootstrap (1000 replicates, 56.7% converged) reported in Results."
  )

  ini({
    # Structural typical values for the reference subject: CYP2D6 EM Western
    # male with Cockcroft-Gault CRCL 80 mL/min and no concomitant CYP3A
    # inhibitor / inducer (Oishi 2014 p. 931 explicit reference definition).
    # Dose (fesoterodine fumarate mg) can be used directly with the reported
    # CL/F and V/F; the paper's Methods 'Data Analysis' states that CL/F and
    # V/F were adjusted by the molecular weights of fesoterodine fumarate
    # (527.65) and 5-HMT (341.5) so no additional MW conversion is required
    # in the model equation to reproduce 5-HMT plasma concentrations.
    lcl   <- log(93.8)   ; label("Apparent oral clearance CL/F at reference subject (L/h)")     # Oishi 2014 Table 3 final CL/F = 93.8 L/h (%SE 2.70)
    lvc   <- log(144)    ; label("Apparent volume of distribution V/F at reference subject (L)") # Oishi 2014 Table 3 final V/F = 144 L (%SE 3.64)
    lka   <- log(0.0935) ; label("Absorption rate constant ka (1/h)")                             # Oishi 2014 Table 3 final ka = 0.0935 1/h (%SE 1.39)
    ltlag <- log(0.447)  ; label("Absorption lag time Tlag (h)")                                  # Oishi 2014 Table 3 final Tlag = 0.447 h (%SE 6.38)

    # Covariate effects on CL/F. Oishi 2014 final-model equation (Results,
    # p. 931-932) encodes each binary covariate as
    #   (theta_k * I + 1 * (1 - I))
    # which equals theta_k when I = 1 and 1 when I = 0. Multiplying the
    # binary factors is mathematically equivalent to
    #   CL/F = CL/F_ref * (CRCL/80)^theta_5 * prod_k theta_k^{I_k} * exp(eta)
    # so the log-additive shift encoding
    #   log(CL/F) = lcl + e_crcl_cl * log(CRCL/80) + sum_k e_k_cl * I_k + eta
    # with e_k_cl = log(theta_k) reproduces the paper values exactly. The
    # power exponent on CRCL is a genuine exponent (not a log-additive shift)
    # and is stored linearly per the parameter-names convention.
    e_crcl_cl          <- 0.303       ; label("Power exponent on (CRCL / 80 mL/min) for CL/F (unitless)")                              # Oishi 2014 Table 3 final theta_5 = 0.303 (%SE 12.8)
    e_hepimp_mod_cl    <- log(0.422)  ; label("Log-additive shift on CL/F for Child-Pugh B (unitless; equivalent to a 0.422 multiplier)") # Oishi 2014 Table 3 final theta_6 = 0.422 (%SE 10.2)
    e_cyp2d6_pm_cl     <- log(0.626)  ; label("Log-additive shift on CL/F for CYP2D6 poor metaboliser (unitless; equivalent to a 0.626 multiplier)") # Oishi 2014 Table 3 final theta_7 = 0.626 (%SE 5.24)
    e_cyp3a4_inh_cl    <- log(0.504)  ; label("Log-additive shift on CL/F for concomitant CYP3A inhibitor (unitless; equivalent to a 0.504 multiplier)") # Oishi 2014 Table 3 final theta_8 = 0.504 (%SE 3.77)
    e_cyp3a4_ind_cl    <- log(3.90)   ; label("Log-additive shift on CL/F for concomitant CYP3A inducer (unitless; equivalent to a 3.90 multiplier)")   # Oishi 2014 Table 3 final theta_9 = 3.90 (%SE 5.23)
    e_sexf_cl          <- log(0.903)  ; label("Log-additive shift on CL/F for female sex (unitless; equivalent to a 0.903 multiplier)") # Oishi 2014 Table 3 final theta_11 = 0.903 (%SE 2.95)
    e_race_japanese_cl <- log(1.12)   ; label("Log-additive shift on CL/F for Japanese ethnicity (unitless; equivalent to a 1.12 multiplier)") # Oishi 2014 Table 3 final theta_13 = 1.12 (%SE 2.95)

    # Full block IIV on CL/F, V/F, ka. Oishi 2014 Table 3 reports the
    # variances (v^2) and covariances (COV) directly on the log-scale eta
    # basis (each eta is exponential IIV per Oishi 2014 Eq. 1). The reported
    # %CV values in Table 3 are computed by the source authors as
    # sqrt(v^2) * 100 (an approximation), so the internal variances are the
    # numbers reported in the v^2 column.
    #   Var(etalcl)         = 0.217   (%CV 46.6)
    #   Cov(etalcl, etalvc) = 0.148   (correlation 0.408)
    #   Var(etalvc)         = 0.606   (%CV 77.8)
    #   Cov(etalcl, etalka) = 0.0457  (correlation 0.381)
    #   Cov(etalvc, etalka) = 0.0953  (correlation 0.475)
    #   Var(etalka)         = 0.0663  (%CV 25.7)
    # Tlag IIV was fixed to zero by the source authors (Results,
    # 'Population Pharmacokinetic Modeling Results': "Interindividual variance
    # of lag time was fixed to zero because it was hard to estimate
    # interindividual variance of absorption rate (ka) and Tlag (lag time)
    # simultaneously." and "Interindividual variance of lag time was fixed to
    # zero because unstructured 4x4 variance-covariance matrix for the etas
    # was unstable."). No eta on ltlag is declared here.
    etalcl + etalvc + etalka ~ c(0.217,
                                 0.148,   0.606,
                                 0.0457,  0.0953, 0.0663)   # Oishi 2014 Table 3 Interindividual Variance block

    # Residual variability: combined additive + proportional error model on
    # 5-HMT plasma concentrations. Oishi 2014 Table 3 reports the residual
    # variances (sigma^2) directly; the paper's Methods 'Data Analysis'
    # states the values are "reported as square root sigma", i.e. the
    # tabulated 'proportional 0.0951 (30.8% CV)' entry gives sigma^2 = 0.0951
    # and %CV = sqrt(0.0951) * 100 = 30.8, and the additive entry
    # '0.00123 (SD = 0.0351)' gives sigma^2 = 0.00123 and SD = sqrt(0.00123)
    # = 0.0351 ng/mL. nlmixr2 propSd / addSd expect the SD directly, so
    # propSd = sqrt(0.0951) = 0.30838 and addSd = sqrt(0.00123) = 0.03507.
    propSd <- 0.30838 ; label("Proportional residual error SD (fraction; = sqrt(0.0951))")  # Oishi 2014 Table 3 sigma^2 proportional = 0.0951 (%SE 3.40)
    addSd  <- 0.03507 ; label("Additive residual error SD (ng/mL; = sqrt(0.00123))")        # Oishi 2014 Table 3 sigma^2 additive = 0.00123 (%SE 42.8)
  })

  model({
    # Individual PK parameters. All covariate effects live on CL/F only; the
    # paper's final model retained no covariates on V/F or ka. Reference
    # subject: CYP2D6 EM Western male with CRCL = 80 mL/min and no
    # concomitant CYP3A modulator (all binary indicators = 0, CRCL = 80
    # gives factor 1).
    cl <- exp(lcl
              + e_hepimp_mod_cl    * HEPIMP_MOD
              + e_cyp2d6_pm_cl     * CYP2D6_PM
              + e_cyp3a4_inh_cl    * CONMED_CYP3A4_INH
              + e_cyp3a4_ind_cl    * CONMED_CYP3A4_IND
              + e_sexf_cl          * SEXF
              + e_race_japanese_cl * RACE_JAPANESE
              + etalcl) * (CRCL / 80) ^ e_crcl_cl
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    # 5-HMT plasma concentrations reported in ng/mL. Dose is entered in mg
    # of fesoterodine fumarate, CL/F in L/h and V/F in L; central/vc gives
    # concentration in mg/L (equivalent to ug/mL), so multiply by 1000 to
    # convert to ng/mL. The paper's CL/F and V/F values are already adjusted
    # by the fesoterodine-fumarate/5-HMT molecular-weight ratio so this
    # unit-scaling factor is the only conversion needed.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
