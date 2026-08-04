Walsh_2024_buprenorphine_desireToUse <- function() {
  description <- paste(
    "Maximum-inhibition (Imax) exposure-response model of the desire to use",
    "visual analog scale (VAS) score under subcutaneous long-acting",
    "buprenorphine (CAM2038 Q1W) in non-treatment-seeking adults with",
    "moderate-to-severe opioid use disorder. The endpoint is the",
    "pre-hydromorphone-challenge desire to use score on a 100 mm unipolar",
    "VAS; the driver is the time-matched buprenorphine plasma concentration",
    "CP_BPN_NGML. Structure (Supplementary Methods, 'Model development:",
    "Desire to use VAS score'):",
    "IPRED = BASE * (1 - Imax * Cp^gamma / (IC50^gamma + Cp^gamma)),",
    "with gamma fixed to 1 and the typical Imax fixed to 1. Unlike the",
    "sister drug-liking model this one is NOT direct: a delay in the onset",
    "of effect was observed in many participants and is modelled by making",
    "Imax an exponential function of time since treatment start,",
    "Imax = TVImax * (1 - exp(-kt * time)), parameterised through an onset",
    "half-life of 0.330 day. The baseline is logit-transformed onto the",
    "range -1 to 101 and its IIV carries a Box-Cox transformation (shape",
    "-12.9) to capture the heavily skewed baseline distribution.",
    "Predictions are logit-transformed onto the same -1 to 101 range with",
    "an additive residual error on that logit scale, encoded as",
    "logitNorm(addSd, -1, 101). The model has NO PK ODE and no",
    "compartments: it is a static algebraic exposure-response relationship",
    "in Cp with an explicit time term, and CP_BPN_NGML is a required",
    "time-varying input covariate supplied per observation record. The",
    "CAM2038 popPK model is not reported in this paper and is not packaged",
    "in nlmixr2lib at extraction time.",
    sep = " "
  )
  reference <- paste(
    "Walsh SL, Comer SD, Aguiar Zdovc J, Sarr C, Bjornsson M,",
    "Strandgarden K, Hjelmstrom P, Tiberg F.",
    "Pharmacokinetic-pharmacodynamic analysis of drug liking blockade by",
    "buprenorphine subcutaneous depot (CAM2038) in participants with opioid",
    "use disorder.",
    "Neuropsychopharmacology. 2024;49(7):1050-1057.",
    "doi:10.1038/s41386-023-01793-z.",
    "Parameter estimates from Supplementary Table S3; structural equation,",
    "onset-delay equation and the logit residual-error transformation from",
    "the Supplementary Appendix (41386_2023_1793_MOESM1_ESM.pdf).",
    "Underlying phase 2 study NCT02611752, reported in",
    "Walsh SL, Comer SD, Lofwall MR, et al. JAMA Psychiatry.",
    "2017;74(9):894-902. doi:10.1001/jamapsychiatry.2017.1874.",
    "Sister endpoint models from the same paper:",
    "modellib('Walsh_2024_buprenorphine_drugLiking').",
    sep = " "
  )
  vignette <- "Walsh_2024_buprenorphine_opioidBlockade"

  depends <- c("CP_BPN_NGML")

  units <- list(
    time          = "day",
    dosing        = "N/A (PD-only; buprenorphine plasma concentration is a required input covariate)",
    concentration = "VAS units (observation: desire to use unipolar VAS score, 0-100 mm); ng/mL (CP_BPN_NGML input covariate)"
  )

  covariateData <- list(
    CP_BPN_NGML = list(
      description        = "Time-matched buprenorphine plasma concentration driving the Imax inhibition of the desire to use VAS score. Time-varying; supplied per observation record in the event table rather than computed from a coupled PK model.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Supplementary Methods, 'Model development: Desire to use VAS score': Cp is 'the time-matched BPN concentration'. Unlike the drug-liking analysis, all hydromorphone challenge levels (0, 6 and 18 mg) contributed here because the desire to use score used was the PRE-challenge recording. Observed BPN plasma concentration range across the phase 2 study was 0.636-12.3 ng/mL (Supplementary Methods, 'Model application'). Member of the canonical CP_<DRUG>_<UNITS> plasma-PD-driver family (siblings CP_MORPH_NGML, CP_OXY_NGML, CP_FBX_NGML). The upstream CAM2038 popPK model is not reported in this paper and is not packaged in nlmixr2lib; concentrations must be supplied externally.",
      source_name        = "Cp (Walsh 2024 Supplementary Methods, desire to use Imax equation)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 1L,
    n_observations = 692L,
    age_range      = "18-54 years (median 36.0; mean 35.8, SD 9.1)",
    age_median     = "36.0 years",
    weight_range   = "53.0-110.0 kg (median 73.4; mean 75.9, SD 14.0)",
    weight_median  = "73.4 kg",
    sex_female_pct = 25.5,
    race_ethnicity = c(Black = 51, White = 47, Other = 2),
    disease_state  = "Non-treatment-seeking adults with moderate to severe opioid use disorder (DSM-5), physically dependent on opioids.",
    dose_range     = "CAM2038 (subcutaneous long-acting buprenorphine) 24 mg (n = 22) or 32 mg (n = 25) once weekly on days 0 and 7. Intramuscular hydromorphone challenges of 0, 6 or 18 mg; all challenge levels contributed to this model because the pre-challenge desire to use score was analysed.",
    regions        = "USA (multisite: University of Kentucky; Columbia University / New York State Psychiatric Institute)",
    biomarkers     = "Desire to use VAS score, collected on pen and paper on a 100 mm unipolar scale (0 = no effect); the pre-hydromorphone-challenge observation was used.",
    notes          = "Baseline demographics from Supplementary Table S1; observation counts from Supplementary Table S2 (desire to use: 692 PD records, 328 + 364 across the 24 mg and 32 mg arms). Study NCT02611752. No covariate analysis is reported for this endpoint (the Methods covariate screen applies to the drug liking IC50 only). Estimation used FOCEI in NONMEM 7.5."
  )

  ini({
    # =====================================================================
    # Structural parameters. Walsh 2024 Supplementary Table S3 ("Parameter
    # estimates of the final desire to use VAS score PK/PD model").
    #
    # Paper structural equation (Supplementary Methods, "Model development:
    # Desire to use VAS score"):
    #
    #                             Imax * Cp^gamma
    #   IPRED = BASE * ( 1 - ------------------------- )
    #                        IC50^gamma + Cp^gamma
    #
    # "Imax was fixed to 1 (i.e. complete inhibition). The sigmoidicity
    # parameter, gamma, was fixed to 1."
    #
    # Back-check on the fixed Hill coefficient: IC90 = IC50 * 9^(1/gamma);
    # with gamma = 1 and IC50 = 0.0129 ng/mL this gives 0.1161 ng/mL, which
    # reproduces the IC90 of 0.116 ng/mL reported independently in the main
    # paper (Results, "Final PK/PD models"). Confirms gamma = 1 as encoded.
    # =====================================================================

    # Baseline desire to use VAS score. Supplementary Methods: "The
    # distribution of baseline desire to use VAS score included values of 0
    # and 100, and therefore a logit transformation of the baseline,
    # similar to what is described for the predictions of desire to use,
    # was applied." The prediction transform pads the 0-100 scale to
    # (-1, 101) via TI = (y + 1) / 102, so the baseline is carried on the
    # logit scale over those same bounds.
    #   logit(92.5, -1, 101) = log(93.5 / 8.5) = log(11) = 2.397895
    logitbase <- logit(92.5, -1, 101)
    label("BASE: baseline desire to use VAS score, on the logit scale bounded (-1, 101) (unitless)")
    # Walsh 2024 Table S3: Baseline = 92.5 VAS units (RSE 0.296%).

    # Box-Cox shape parameter applied to the baseline IIV. Supplementary
    # Methods: "The baseline distribution was also heavily skewed, which
    # was modeled using a Box-Cox transformation of the baseline
    # variability." The standard NONMEM Box-Cox eta transform is
    #   eta_tr = (exp(eta)^lambda - 1) / lambda,
    # which is 0 at eta = 0 (so the typical value is unchanged) and, for a
    # large negative lambda, is bounded above by -1/lambda = 0.0775 while
    # remaining unbounded below -- producing the strong left skew the paper
    # describes for a baseline sitting near the top of the VAS range.
    lambda_base <- -12.9
    label("Box-Cox shape parameter for the baseline IIV distribution (unitless)")
    # Walsh 2024 Table S3: Baseline Box-Cox shape parameter = -12.9 (RSE 4.91%).

    # Buprenorphine plasma concentration producing half of the maximum
    # inhibition of desire to use.
    lic50 <- log(0.0129)
    label("IC50: BPN plasma concentration giving half-maximal desire-to-use inhibition (ng/mL)")
    # Walsh 2024 Table S3: IC50 = 0.0129 ng/mL (RSE 12.3%).

    # Typical maximum fractional inhibition, fixed to 1 by the authors.
    imax <- fixed(1)
    label("TVImax: typical maximum fractional inhibition of desire to use (unitless)")
    # Walsh 2024 Table S3: Imax = 1.00 (FIXED).

    # Sigmoidicity (Hill) coefficient, fixed to 1 by the authors.
    hill <- fixed(1)
    label("gamma: sigmoidicity (Hill) coefficient of the Imax relationship (unitless)")
    # Walsh 2024 Supplementary Methods: "The sigmoidicity parameter, gamma, was fixed to 1."

    # Onset half-life of the drug effect. Supplementary Methods: "In many
    # participants a delay in the onset of the effect was observed. This was
    # modeled by describing Imax as an exponential function of time,
    # Imax = TVImax * (1 - exp(-kt * time))". Table S3 reports the onset
    # HALF-LIFE (0.330 day) rather than kt itself; kt = ln(2) / t_half =
    # 2.100 /day is derived in model().
    lthalf_onset <- log(0.330)
    label("Onset half-life of the buprenorphine effect on desire to use (day)")
    # Walsh 2024 Table S3: Onset half-life = 0.330 day (RSE 18.7%).

    # =====================================================================
    # Inter-individual variability. Table S3 heads the IIV column "CV",
    # whereas the main-paper Table 1 heads the equivalent column "SD". The
    # tabulated quantities are omega (SD on the log scale) in both tables:
    #
    #  (a) The Supplementary Methods ("Parameter Estimation") describe IIV
    #      in exponential form P_i = TVP * exp(eta_i) with eta_i ~ N(0,
    #      omega^2), i.e. the natural reporting quantity is omega.
    #  (b) The same Table S3 column also carries the RUV entry (0.783).
    #      That RUV is ADDITIVE on the logit scale (Supplementary Methods
    #      equations), a quantity that has no coefficient-of-variation
    #      interpretation at all -- it can only be an SD. Its value is also
    #      near-identical to the drug-liking model's additive logit-scale
    #      RUV, which Table 1 explicitly labels "SD" (0.744).
    #
    # The "CV" column header in Table S3 is therefore read as a labelling
    # slip for the SD scale, and the values are squared here to give the
    # variances rxode2 expects. See the vignette "Assumptions and
    # deviations" section.
    # =====================================================================
    etalogitbase ~ 0.167^2
    # Walsh 2024 Table S3: IIV Baseline = 0.167 (RSE 14.1%, shrinkage 5.76%).
    etalic50 ~ 9.29^2
    # Walsh 2024 Table S3: IIV IC50 = 9.29 (RSE 21.7%, shrinkage 42.2%).
    etalthalf_onset ~ 1.91^2
    # Walsh 2024 Table S3: IIV onset half-life = 1.91 (RSE 17.6%, shrinkage 33.9%).

    # =====================================================================
    # Residual unexplained variability. Additive on the LOGIT-transformed
    # prediction scale. Supplementary Methods ("Model development: Desire
    # to use VAS score") give the transformation explicitly:
    #
    #   TI_yhat  = (yhat + 1) / 102
    #   PHI_yhat = ln( TI_yhat / (1 - TI_yhat) )
    #   y        = expit(PHI_yhat + eps_add) * 102 - 1
    #
    # i.e. a logit-normal residual bounded on (-1, 101) -- exactly
    # rxode2's logitNorm(sd, low = -1, hi = 101). The bounds pad the 0-100
    # unipolar VAS by one unit either side so that predictions "of exactly
    # 0 and 100" remain attainable.
    # =====================================================================
    addSd <- 0.783
    label("Additive residual SD on the logit-transformed desire to use VAS scale (unitless)")
    # Walsh 2024 Table S3: RUV = 0.783 (RSE 8.45%, shrinkage 5.68%).
  })

  model({
    # -------------------------------------------------------------------
    # Baseline. Box-Cox transform of the baseline eta, then back-transform
    # from the logit scale bounded (-1, 101). At eta = 0 the Box-Cox term
    # is 0 and base = expit(logit(92.5, -1, 101), -1, 101) = 92.5.
    # -------------------------------------------------------------------
    etabase_bc <- (exp(etalogitbase)^lambda_base - 1) / lambda_base
    base       <- expit(logitbase + etabase_bc, -1, 101)

    # -------------------------------------------------------------------
    # Remaining individual parameters (exponential IIV).
    # -------------------------------------------------------------------
    ic50        <- exp(lic50 + etalic50)
    thalf_onset <- exp(lthalf_onset + etalthalf_onset)

    # -------------------------------------------------------------------
    # Time-dependent onset of the maximum inhibitory effect:
    #   Imax(t) = TVImax * (1 - exp(-kt * time)),  kt = ln(2) / t_half
    # `time` is time since the start of CAM2038 treatment, in days. At
    # time = 0 there is no effect (Imax = 0) and the prediction equals the
    # baseline; Imax approaches TVImax with a half-life of 0.330 day.
    # -------------------------------------------------------------------
    kt     <- log(2) / thalf_onset
    imax_t <- imax * (1 - exp(-kt * time))

    # -------------------------------------------------------------------
    # Imax exposure-response on the desire to use VAS score.
    # -------------------------------------------------------------------
    dtu <- base * (1 - imax_t * CP_BPN_NGML^hill / (ic50^hill + CP_BPN_NGML^hill))

    # -------------------------------------------------------------------
    # Logit-normal residual error on the bounded range -1 to 101.
    # -------------------------------------------------------------------
    dtu ~ logitNorm(addSd, -1, 101)
  })
}
