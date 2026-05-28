Mazzocco_2015_temozolomide <- function() {
  description <- "Tumour growth inhibition (TGI) model for low-grade glioma (LGG) treated with first-line temozolomide chemotherapy (Mazzocco 2015): three tumour-tissue compartments (proliferative, non-damaged quiescent, damaged quiescent) coupled to a K-PD virtual drug compartment, with logistic proliferative growth (carrying capacity K fixed at 100 mm), treatment-induced damage of both proliferative and quiescent tissues, time-dependent acquired resistance of the proliferative tissue only, and tumour-genotype covariate effects of TP53 mutation status on TMZ efficacy and 1p/19q codeletion status on the damaged-quiescent-to-proliferative repair rate. Observation is mean tumour diameter (MTD = P + Q + Qp) in millimetres."
  reference <- paste(
    "Mazzocco P, Barthelemy C, Kaloshi G, Lavielle M, Ricard D,",
    "Idbaih A, Psimaras D, Renard MA, Alentorn A, Honnorat J,",
    "Delattre JY, Ducray F, Ribba B.",
    "Prediction of Response to Temozolomide in Low-Grade Glioma Patients",
    "Based on Tumor Size Dynamics and Genetic Characteristics.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(12):728-737.",
    "doi:10.1002/psp4.54.",
    "Carrying capacity K is fixed at 100 mm (= 10 cm), inherited from the",
    "upstream Ribba 2012 LGG model (Clin Cancer Res. 2012;18:5071-5080;",
    "doi:10.1158/1078-0432.CCR-12-0084). Provenance for the K = 10 cm value",
    "is the on-disk Ribba 2014 review (Table 1 footnote on the Ribba 2012",
    "row: 'theta was not estimated and was fixed to 10 cm');",
    "CPT Pharmacometrics Syst Pharmacol. 2014;3:e113; doi:10.1038/psp.2014.12.",
    sep = " "
  )
  vignette <- "Mazzocco_2015_temozolomide"
  units <- list(
    time          = "month",
    dosing        = "arbitrary unit per TMZ cycle (K-PD bolus to the kpdConc virtual drug compartment; the source paper represents each 5-day daily-dosing TMZ cycle as a single bolus of arbitrary magnitude, with the dose units absorbed into the typical-value gamma)",
    concentration = "mm (mean tumour diameter MTD = P + Q + Qp; not a drug concentration)"
  )

  covariateData <- list(
    TUM_TP53_MUT = list(
      description        = "Tumour TP53 mutation indicator (1 = TP53 mutant tumour, 0 = TP53 wild-type tumour). p53 protein overexpression by IHC is used by the source paper as a surrogate marker for TP53 missense mutations (Gillet et al. J Neurooncol 2014).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed per subject (somatic tumour genotype call at diagnosis). Acts on the TMZ tumour-cell death-rate constant gamma: gamma_mut / gamma_wt = 0.143 / 0.254 = 0.563 (44% lower TMZ efficacy in TP53-mutant LGG, per Mazzocco 2015 Table 2 and Results). Cohort frequency: 24 mutant / 59 with p53 status known = 41% mutant.",
      source_name        = "p53 mutation"
    ),
    TUM_1P19Q_CODEL = list(
      description        = "Tumour 1p/19q chromosomal codeletion indicator (1 = combined 1p and 19q loss, 0 = non-codeleted).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed per subject (somatic tumour genotype call at diagnosis). Acts on the damaged-quiescent-to-proliferative repair rate constant kQpP: kQpP_codel / kQpP_noncodel = 0.00807 / 0.00947 = 0.852 (15% lower DNA-damage-repair rate in 1p/19q-codeleted tumours, consistent with the longer duration of response reported in codeleted LGG patients; Ricard 2007). Mutually exclusive with TP53 mutation in this cohort. Cohort frequency: 23 codeleted / 70 with 1p/19q status known = 33% codeleted.",
      source_name        = "1p/19q codeletion"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 77L,
    n_studies           = 1L,
    age_range           = "25-71 years",
    age_median          = "40 years",
    sex_female_pct      = 45.5,
    disease_state       = "WHO grade II low-grade glioma (oligodendroglioma 73%, oligoastrocytoma 21%, astrocytoma 6%) at first-line chemotherapy onset",
    dose_range          = "Temozolomide 200 mg/m^2/day orally on days 1-5 of each 28-day cycle (median 18 cycles, range 2-24)",
    regions             = "France",
    notes               = "77 patients with at least one molecular characteristic known (1p/19q codeletion, TP53 / p53 mutation, or IDH mutation) from a 120-patient single-centre French cohort treated 1999-2007 (Ricard 2007, Ann Neurol). Used for population-parameter estimation. The remaining 43 patients lacked molecular characterisation and were held out as an external evaluation cohort (not used for the parameters extracted here). Tumour size was measured as mean tumour diameter (MTD) from printed MRI scans via MTD = (2V)^(1/3) where V = (D1 * D2 * D3) / 2 with D1, D2, D3 the three largest perpendicular tumour diameters. 952 MTD observations total (mean 12 per patient, range 4-28). Median post-treatment follow-up 21 months (range 5 months to 9.5 years). 1p/19q codeletion and TP53 missense mutation were mutually exclusive in this cohort."
  )

  ini({
    # Structural parameters -- Mazzocco 2015 Table 2 final-model column.
    # All rate constants are in 1/month (paper-native time unit). Initial
    # sizes are in mm.
    lp0       <- log(1.72);    label("Initial proliferative tissue P0 (mm)")                                       # Mazzocco 2015 Table 2: P0 = 1.72 mm (RSE 21%)
    lq0       <- log(32.1);    label("Initial non-damaged quiescent tissue Q0 (mm)")                                # Mazzocco 2015 Table 2: Q0 = 32.1 mm (RSE 7%)
    llambdap  <- log(0.143);   label("Proliferative-tissue growth rate lambda_P (1/month)")                         # Mazzocco 2015 Table 2: lambda_P = 0.143 (RSE 12%)
    lkpq      <- log(0.0429);  label("Proliferative->quiescent transition rate kPQ (1/month)")                      # Mazzocco 2015 Table 2: kPQ = 0.0429 (RSE 21%)
    lkqpp     <- log(0.00947); label("Damaged-quiescent->proliferative repair rate kQpP, 1p/19q non-codeleted reference (1/month)") # Mazzocco 2015 Table 2: kQpP non-codel = 0.00947 (RSE 42%)
    ldqp      <- log(0.0188);  label("Damaged-quiescent death rate dQp (1/month)")                                  # Mazzocco 2015 Table 2: dQp = 0.0188 (RSE 19%)
    lgamma    <- log(0.254);   label("TMZ tumour-cell death rate constant gamma, p53 wild-type reference (unitless)") # Mazzocco 2015 Table 2: gamma p53 wild = 0.254 (RSE 18%)
    lres      <- log(0.1);     label("Acquired-resistance rate res (1/month; multiplies exp(-res*t) on the proliferative-tissue drug effect)") # Mazzocco 2015 Table 2: res = 0.1 (RSE 22%)

    # KDE is FIXED at 8.3 /month per Mazzocco 2015 Table 2 (paper text:
    # "The population value of KDE parameter was fixed to 8.3 month^-1
    # corresponding to a half-life of 2.5 days, allowing for a residual
    # active concentration of TMZ after 5 days of treatment.").
    lkde      <- fixed(log(8.3)); label("Virtual-drug elimination rate constant KDE (1/month)")                     # Mazzocco 2015 Table 2: KDE = 8.3 FIXED

    # Carrying capacity K fixed at 100 mm (= 10 cm), inherited from the
    # upstream Ribba 2012 LGG model. Ribba 2014 review Table 1 footnote on
    # the Ribba 2012 row records: "theta was not estimated and was fixed
    # to 10 cm". The Ribba 2014 review PDF is the on-disk provenance
    # source for this value; the original Ribba 2012 publication was not
    # on disk for the present extraction.
    lK        <- fixed(log(100)); label("Carrying capacity K (mm)")                                                  # Ribba 2014 review Table 1 footnote: K = 10 cm = 100 mm

    # Covariate effects on log-scale typical values. Both coefficients are
    # back-derived from the per-stratum point estimates in Mazzocco 2015
    # Table 2:
    #   gamma = gamma_p0 * exp(e_p53_gamma   * TUM_TP53_MUT)
    #   kQpP  = kQpP_p0  * exp(e_codel_kqpp * TUM_1P19Q_CODEL)
    # e_p53_gamma   = log(0.143 / 0.254) = -0.5743
    # e_codel_kqpp  = log(0.00807 / 0.00947) = -0.1601
    e_p53_gamma   <- -0.5743; label("Effect of TP53 mutation on log(gamma) (log-fold change)")                      # log(0.143 / 0.254) per Mazzocco 2015 Table 2: gamma_mut = 0.143, gamma_wt = 0.254
    e_codel_kqpp  <- -0.1601; label("Effect of 1p/19q codeletion on log(kQpP) (log-fold change)")                    # log(0.00807 / 0.00947) per Mazzocco 2015 Table 2: kQpP_codel = 0.00807, kQpP_noncodel = 0.00947

    # Inter-individual variability (IIV). Mazzocco 2015 Table 2 reports
    # CV(%) for each random effect; for log-normal random effects the
    # corresponding NONMEM-style variance is omega^2 = log((CV/100)^2 + 1).
    etalp0       ~ 1.1136   # CV 143% -> omega^2 = log(1.43^2 + 1) per Mazzocco 2015 Table 2
    etalq0       ~ 0.2712   # CV 55.8% -> omega^2 = log(0.558^2 + 1) per Mazzocco 2015 Table 2
    etallambdap  ~ 0.3352   # CV 63.1% -> omega^2 = log(0.631^2 + 1) per Mazzocco 2015 Table 2
    etalkpq      ~ 0.5045   # CV 81% -> omega^2 = log(0.81^2 + 1) per Mazzocco 2015 Table 2
    etalkqpp     ~ 1.2876   # CV 162% -> omega^2 = log(1.62^2 + 1) per Mazzocco 2015 Table 2
    etaldqp      ~ 0.5556   # CV 86.2% -> omega^2 = log(0.862^2 + 1) per Mazzocco 2015 Table 2
    etalgamma    ~ 0.3857   # CV 68.6% -> omega^2 = log(0.686^2 + 1) per Mazzocco 2015 Table 2
    etalres      ~ 0.4996   # CV 80.5% -> omega^2 = log(0.805^2 + 1) per Mazzocco 2015 Table 2

    # KDE IIV is reported as "50 (FIXED)" in Table 2 -- the value of the
    # variance is held constant rather than estimated. The corresponding
    # log-normal omega^2 is log(0.50^2 + 1) = 0.2231.
    etalkde      ~ fixed(0.2231)  # CV 50% (FIXED) -> omega^2 = log(0.50^2 + 1) per Mazzocco 2015 Table 2

    # Constant additive residual error on MTD ("constant error model with
    # parameter value a", Mazzocco 2015 Methods).
    addSd <- 1.73; label("Additive residual SD on MTD (mm)")  # Mazzocco 2015 Table 2: a = 1.73 mm (RSE 3%)
  })

  model({
    # Individual structural parameters (lognormal IIV). Covariate-modified
    # parameters apply their log-fold-change coefficients here on the log
    # scale so the exponent gives a clean multiplicative covariate effect.
    p0       <- exp(lp0 + etalp0)
    q0       <- exp(lq0 + etalq0)
    lambdap  <- exp(llambdap + etallambdap)
    kpq      <- exp(lkpq + etalkpq)
    kqpp     <- exp(lkqpp + etalkqpp + e_codel_kqpp * TUM_1P19Q_CODEL)
    dqp      <- exp(ldqp + etaldqp)
    gamma    <- exp(lgamma + etalgamma + e_p53_gamma * TUM_TP53_MUT)
    res      <- exp(lres + etalres)
    kde      <- exp(lkde + etalkde)
    K        <- exp(lK)

    # K-PD virtual drug + three tumour-tissue ODEs (Mazzocco 2015 p.731
    # equation set; Figure 2 schematic):
    #   dC/dt   = -KDE * C
    #   dP/dt   = lambdaP * P * (1 - Pstar/K) + kQpP * Qp - kPQ * P
    #             - gamma * exp(-res*t) * C * KDE * P
    #   dQ/dt   = kPQ * P - gamma * C * KDE * Q
    #   dQp/dt  = gamma * C * KDE * Q - kQpP * Qp - dQp * Qp
    #   Pstar   = P + Q + Qp                                       (observed MTD)
    # The acquired-resistance factor exp(-res*t) decays the TMZ effect on
    # the proliferative tissue only (Mazzocco 2015 Methods: "the final
    # selected model incorporated a resistance term for the proliferative
    # tissue only"); the effect on the non-damaged quiescent tissue Q is
    # not subject to resistance. KDE multiplies inside the drug-effect
    # term as in Mazzocco 2015 Eq. 2.
    #
    # The time variable t in the resistance term is the rxode2 simulation
    # time. Convention for users: set t = 0 at the start of TMZ therapy;
    # supply pretreatment observations at t < 0. The K-PD compartment
    # initial value is 0 so the drug-effect term vanishes for t < 0
    # regardless of the sign of exp(-res*t).
    tumor_size <- prolif + quiesc + quiescDam

    drugEffect       <- gamma * kde * kpdConc
    drugEffectProlif <- drugEffect * exp(-res * t)

    d/dt(kpdConc)   <- -kde * kpdConc
    d/dt(prolif)    <-  lambdap * prolif * (1 - tumor_size / K) +
                        kqpp * quiescDam -
                        kpq * prolif -
                        drugEffectProlif * prolif
    d/dt(quiesc)    <-  kpq * prolif - drugEffect * quiesc
    d/dt(quiescDam) <-  drugEffect * quiesc -
                        kqpp * quiescDam -
                        dqp * quiescDam

    prolif(0)    <- p0
    quiesc(0)    <- q0
    quiescDam(0) <- 0
    kpdConc(0)   <- 0

    # Observation: mean tumour diameter (MTD) in mm, constant additive
    # residual error.
    tumor_size ~ add(addSd)
  })
}
