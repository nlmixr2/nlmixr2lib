vanHasselt_2015_eribulin <- function() {
  description <- "Disease-progression (DP) model for prostate-specific antigen (PSA) dynamics in metastatic castration-resistant prostate cancer (CRPC) patients treated with eribulin mesilate (van Hasselt 2015). K-PD framework: the per-dose predicted eribulin AUC enters a single transient drug-effect compartment depot_kpd that decays with rate KP (fixed to 6000 /day so the effect is nearly instantaneous after each dose); PSA evolves under a first-order growth rate KG counteracted by an inhibition rate KD0 multiplied by the K-PD state depot_kpd and an exponentially decaying resistance factor exp(-k_res*t). PSA0, KD0, KG, k_res have correlated lognormal IIV; proportional residual error on PSA (log-transform-both-sides). Prior taxane treatment (binary PRIOR_TAXANE) multiplies PSA0; cumulative number of days of prior taxane treatment (continuous PRIOR_TAXANE_DAYS) enters KD0 as (1 + NTRT/720)^theta. The companion parametric Weibull survival sub-model fit in R survreg is documented in the vignette but not encoded here (not an ODE / nlmixr2 structure)."
  reference <- paste(
    "van Hasselt JGC, Gupta A, Hussein Z, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Disease Progression/Clinical Outcome Model for Castration-Resistant Prostate Cancer",
    "in Patients Treated With Eribulin.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(7):386-395.",
    "doi:10.1002/psp4.49.",
    sep = " "
  )
  vignette <- "vanHasselt_2015_eribulin"
  units <- list(
    time          = "day",
    dosing        = "ng*h/mL (predicted per-dose eribulin AUC, supplied as the amt for each dose event into compartment depot_kpd via the K-PD framework; the upstream eribulin popPK model -- Majid 2014 J Clin Pharmacol 54:1134; van Hasselt 2013 Br J Clin Pharmacol 76:412 -- is not encoded here)",
    concentration = "ng/mL (serum PSA)"
  )

  covariateData <- list(
    PRIOR_TAXANE = list(
      description        = "Binary indicator of prior taxane chemotherapy at study entry (1 = received any prior taxane regimen, 0 = taxane-naive).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (taxane-naive)",
      notes              = "van Hasselt 2015 Methods 'Clinical study': 108 CRPC patients, 50 prior-docetaxel and 58 taxane-naive (paper notation PTAX = 1/0). Multiplies baseline PSA0 as psa0 = exp(lpsa0 + etalpsa0) * e_prior_taxane_psa0^PRIOR_TAXANE, encoding the Table 2 footnote 'Individual parameters were defined as: ... PSA0 = hPSA0 * hPSA0-PTAX * exp(gPSA0)'. Time-invariant within a subject (treatment history at study entry).",
      source_name        = "PTAX"
    ),
    PRIOR_TAXANE_DAYS = list(
      description        = "Cumulative number of days of prior taxane treatment at study entry (0 for taxane-naive patients).",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "van Hasselt 2015 Methods 'Covariates model development': normalised by the population median 720 days. Enters KD0 as a power covariate kd0 = exp(lkd0 + etalkd0) * (1 + PRIOR_TAXANE_DAYS / 720)^e_prior_taxane_days_kd0 per the paper's Eq. 4 (Pg = hP * (1 + NTRT_i / 720)^h_NTRT). For patients with PRIOR_TAXANE = 0 (taxane-naive), PRIOR_TAXANE_DAYS = 0 by construction and the multiplier (1 + 0/720)^theta = 1 collapses to no covariate effect on KD0. Time-invariant within a subject.",
      source_name        = "NTRT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 108L,
    n_studies       = 1L,
    age_range       = "adult metastatic CRPC patients (median ages by stratum reported in Supporting Table S1 of van Hasselt 2015; specific range not on disk in the main paper)",
    weight_range    = "not reported in the main paper text on disk",
    sex_female_pct  = 0,
    race_ethnicity  = "not reported in the main paper text on disk",
    disease_state   = "Metastatic castration-resistant prostate cancer (CRPC). Phase II trial of eribulin mesilate (E7389), de Bono et al. Ann Oncol 2012;23:1241; van Hasselt 2015 reference 33. Subset stratification: 50 patients with prior docetaxel/taxane therapy, 58 taxane-naive.",
    dose_range      = "eribulin mesilate IV (per de Bono 2012 phase II protocol); per-dose AUC is provided as the amt input to the K-PD depot_kpd compartment in this model, not the mg dose itself. The upstream eribulin popPK model used to predict individual AUC (3-compartment linear elimination with albumin / alkaline phosphatase / total bilirubin on CL, Majid 2014 J Clin Pharmacol 54:1134-1143 and van Hasselt 2013 Br J Clin Pharmacol 76:412-424) is referenced but not encoded here.",
    regions         = "not reported in the main paper text on disk",
    notes           = "108 metastatic CRPC patients from a single Phase II trial of eribulin mesilate (de Bono et al. 2012). PSA-time profiles analysed with NONMEM 7.2 FOCE; survival analysed separately in R survreg. The PD model uses a K-PD approach where the predicted eribulin AUC per dose (from the external popPK model) drives a transient drug-effect compartment. PK data were not available for this study (van Hasselt 2015 Methods 'Pharmacokinetic model')."
  )

  ini({
    # Structural parameters from van Hasselt 2015 Table 2 ('Fixed effect
    # parameters' block, final-model column 'Estimate').
    lkel   <- fixed(log(6000));    label("K-PD elimination rate constant KP (1/day; FIXED per paper)") # van Hasselt 2015 Table 2: hKP = 6000 /day, marked '+' = FIXED in the table footnote. Paper Results: 'The optimal value of the drug exposure parameter (KP) was selected based on evaluation of different fixed values in the model. We selected a large value of 6,000 to allow for a nearly instantaneous dosing event.'
    lkd0   <- log(0.241);          label("Drug PSA inhibition rate KD0 (L/(ng*h*day); the paper's unit string 'ng*h/L /day /day' is dimensionally the same when KD0*D has units 1/day with D in ng*h/mL)") # van Hasselt 2015 Table 2: hKD0 = 0.241 (RSE 32.6%)
    lkres  <- log(0.0113);         label("Drug-resistance development rate k_res (1/day; multiplies exp(-k_res*t))") # van Hasselt 2015 Table 2: hk = 0.0113 (RSE 44.3%)
    lkg    <- log(0.00879);        label("PSA first-order growth rate KG (1/day)") # van Hasselt 2015 Table 2: hKG = 0.00879 (RSE 12.6%)
    lpsa0  <- log(23.2);           label("Baseline PSA at start of treatment (ng/mL; taxane-naive typical value)") # van Hasselt 2015 Table 2: hPSA0 = 23.2 ng/mL (RSE 16.5%)

    # Covariate-effect parameters from van Hasselt 2015 Table 2 (final-model
    # multivariate covariate column). Functional forms are recorded in the
    # Table 2 footnote: 'PSA0 = hPSA0 * hPSA0-PTAX * exp(gPSA0); KD0 = hKD0 *
    # (1 + (NTRT/720))^hKD-NTRT * exp(gKD0)' (the same equation appears as
    # Eq. 4 of the main text for the continuous-covariate power form).
    e_prior_taxane_psa0      <- 3.23;  label("Multiplicative effect of binary PRIOR_TAXANE on PSA0 (unitless; PSA0 ratio prior-taxane vs taxane-naive)") # van Hasselt 2015 Table 2: hPSA0-PTAX = 3.23 (RSE 27.6%)
    e_prior_taxane_days_kd0  <- -4.00; label("Power exponent of (1 + PRIOR_TAXANE_DAYS/720) on KD0 (unitless)") # van Hasselt 2015 Table 2: hKD0-NTRT = -4.00 (RSE 52.5%)

    # Inter-individual variability -- van Hasselt 2015 Table 2
    # 'Interindividual variability variances'. Paper reports CV%; converted
    # to log-normal variances via omega^2 = log(1 + CV^2).
    #   xKD0  = 127.3 % -> log(1 + 1.273^2)  = 0.96334
    #   xk    =  88.3 % -> log(1 + 0.883^2)  = 0.57665
    #   xKG   =  53.7 % -> log(1 + 0.537^2)  = 0.25342
    #   xPSA0 = 130.4 % -> log(1 + 1.304^2)  = 0.99340
    # Off-diagonal covariances from Table 2 footnote correlations:
    #   xk *xKD0  =  0.802  -> cov = 0.802 * sqrt(0.96334 * 0.57665) = 0.59776
    #   xKG*xKD0  = -0.293  -> cov = -0.14474
    #   xKG*xk    = -0.111  -> cov = -0.04243
    #   xPSA0*xk  = -0.094  -> cov = -0.07116
    #   xPSA0*xKG = -0.032  -> cov = -0.01605
    # The xPSA0 x xKD0 covariance is fixed at 0 per the paper Results
    # ('Statistical model' paragraph): 'Covariances between IIV random
    # effects could be estimated for all model parameters, except for the
    # correlation between the drug effect parameter KD and PSA0, which
    # approached zero in the final obtained estimate.'
    etalkd0 + etalkres + etalkg + etalpsa0 ~ c(
      0.96334,
      0.59776,    0.57665,
      -0.14474,  -0.04243,   0.25342,
      fixed(0),  -0.07116,  -0.01605,   0.99340
    )

    # Residual unexplained variability -- van Hasselt 2015 Table 2
    # 'Residual variability variance' row: proportional CV% = 34.2% (RSE
    # 27.5%, shrinkage 14%). The paper Methods 'Structural model' notes the
    # PSA observations were log-transformed before fitting (LTBS); a
    # proportional error on linear-scale PSA is the standard nlmixr2
    # equivalent of NONMEM additive-on-log-PSA (the SD of the log-additive
    # error equals the linear-space CV at small CV, and our CV = 0.342 is
    # in that regime to within ~5%).
    propSd <- 0.342; label("Proportional residual SD on serum PSA (fraction; LTBS-equivalent)") # van Hasselt 2015 Table 2: rprop = 34.2 CV%
  })

  model({
    # Individual structural parameters with lognormal IIV. Covariate forms
    # follow the Table 2 footnote ('Individual parameters were defined as:
    # KD0 = hKD0 * (1 + (NTRT/720))^hKD-NTRT * exp(gKD0); ...; PSA0 = hPSA0
    # * hPSA0-PTAX * exp(gPSA0)') and Eq. 4 of the main text.
    kel    <- exp(lkel)                                             # K-PD decay rate, fixed
    kd0    <- exp(lkd0   + etalkd0)  * (1 + PRIOR_TAXANE_DAYS / 720)^e_prior_taxane_days_kd0
    kres   <- exp(lkres  + etalkres)
    kg     <- exp(lkg    + etalkg)
    psa0   <- exp(lpsa0  + etalpsa0) * e_prior_taxane_psa0^PRIOR_TAXANE

    # K-PD drug-effect compartment depot_kpd (the paper's 'D'). Each dose event
    # delivers the per-dose predicted eribulin AUC (units ng*h/mL) as a
    # bolus into depot_kpd via the standard rxode2 dosing mechanism (amt column
    # on an evid=1 row with cmt = 'depot_kpd'). depot_kpd decays at the fixed
    # rate kel = 6000 /day, giving a near-instantaneous impulse whose integrated
    # area equals amt/kel. The upstream popPK model that produces the
    # per-dose AUC values is NOT encoded here -- the user assembles the
    # event table with AUC values supplied externally (see vignette).
    d/dt(depot_kpd) <- -kel * depot_kpd

    # PSA dynamics. Paper Eq. 2:
    #   dPSA/dt = KG * PSA - KD0 * exp(-k_res * t) * D(t) * PSA
    # with PSA(0) = PSA0_i (individual baseline). t is rxode2 simulation
    # time in days; treatment start is taken as t = 0 so exp(-k_res * t)
    # decays from 1 at study start, encoding progressive drug-resistance
    # development.
    d/dt(PSA)   <- PSA * (kg - kd0 * exp(-kres * t) * depot_kpd)
    PSA(0)      <- psa0

    PSA ~ prop(propSd)
  })
}
