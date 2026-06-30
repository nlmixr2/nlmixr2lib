Wright_2016_allopurinol <- function() {
  description <- "One-compartment population PK-PD model for allopurinol (via the active metabolite oxypurinol) and plasma urate in adults with gout (Wright 2016 BJCP). Oxypurinol disposition is a one-compartment first-order absorption / first-order elimination model with Ka fixed at 1.09 1/h; apparent oxypurinol clearance (CL/F_oxy) is allometrically scaled on fat-free mass (Janmahasatian formula, exponent fixed at 0.75) and power-scaled on Cockcroft-Gault creatinine clearance standardised to 70 kg, with a multiplicative reduction when a thiazide or loop diuretic is coadministered; apparent volume (V/F_oxy) is allometrically scaled on total body weight (exponent fixed at 1.0) and shares its IIV with CL/F_oxy via a fixed fractional scaler (Bonate 2006 fractional-effect parameterisation). Plasma urate is described by a direct-effect sigmoidal Emax inhibition of urate production on top of a baseline urate U0 that is power-scaled on renal function and multiplicatively higher with concomitant diuretic. The dose entered into the model is allopurinol oral mg; the implicit 1:1 molar conversion to oxypurinol is absorbed into the apparent CL/F_oxy and V/F_oxy. Cc is oxypurinol concentration (umol/L) and Eurate is plasma urate (mmol/L)."
  reference   <- "Wright DFB, Duffull SB, Merriman TR, Dalbeth N, Barclay ML, Stamp LK. Predicting allopurinol response in patients with gout. Br J Clin Pharmacol. 2016 Feb;81(2):277-289. doi:10.1111/bcp.12799"
  vignette    <- "Wright_2016_allopurinol"
  units       <- list(time = "hour", dosing = "mg", concentration = "umol/L (oxypurinol Cc); mmol/L (urate Eurate)")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass derived from total body weight, height, and sex by the Janmahasatian 2005 semi-mechanistic formula.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wright 2016 Methods (Covariate models): 'Fat-free mass (FFM) was calculated using the formula developed by Janmahasatian et al. [34]'. Enters apparent oxypurinol clearance CL/F_oxy via allometric power scaling with fixed exponent 0.75 and reference 70 kg (Wright 2016 final-model equation set; Table 3 footnote *: 'clearance expressed per 70 kg FFM/CLcr 6 l h^-1').",
      source_name        = "FFM"
    ),
    WT = list(
      description        = "Total body weight (TBW).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wright 2016 Methods (Covariate models). Used for allometric scaling of apparent oxypurinol volume V/F_oxy with fixed exponent 1.0 and reference 70 kg (Wright 2016 final-model equation; Table 3 footnote dagger: 'volume expressed per 70 kg body weight'). Not used on CL/F_oxy; the paper explicitly notes that TBW as a covariate on oxypurinol clearance did not provide a better fit than FFM (Wright 2016 Results paragraph 'Total body weight as a covariate on oxypurinol clearance did not provide a better fit to the data').",
      source_name        = "TBW"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault formula and standardised to 70 kg (NOT BSA-normalised). Expressed as L/h, not as mL/min.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wright 2016 Methods (Covariate models): 'Renal function was calculated using the Cockroft-Gault formula [38] and expressed as creatinine clearance (CLcr) standardized to 70 kg [39]. Renal function (RF) was then normalized to a standard creatinine clearance (CLcrSTD) of 6 l h^-1/70 kg (100 ml min^-1/70 kg)'. Values are stored in L/h (so a typical normal-renal-function adult is 6 L/h; the cohort median 68 mL/min / 70 kg ~= 4.08 L/h). Enters apparent oxypurinol clearance via power scaling (CRCL/6)^0.587 and enters baseline urate U0 via power scaling (CRCL/6)^(-0.119). Stored under the canonical CRCL register entry with the explicit Cockcroft-Gault-standardised-to-70 kg assay form documented here (mirrors the Stocker 2012 oxypurinol precedent which uses raw Cockcroft-Gault mL/min normalised on lean body weight; the two papers use different normalisations of the same biological quantity).",
      source_name        = "CLcr"
    ),
    CONMED_DIUR = list(
      description        = "Concomitant diuretic indicator. Wright 2016's CONMED_DIUR captures thiazide diuretics OR loop diuretics; potassium-sparing diuretics (spironolactone, amiloride) are NOT pooled in.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant thiazide or loop diuretic)",
      notes              = "Wright 2016 Methods (Covariate models): 'drugs associated with an increased or decreased risk of hyperuricaemia were tested in the PKPD model including thiazide or loop diuretics ...'. Wright 2016 Table 4 footer: 'Diuretics include thiazides and loop diuretics'. 33% of the cohort (44 of 133 patients) were on a thiazide or loop diuretic at study entry (Wright 2016 Table 2 totals). This is NARROWER than the Stocker_2012_oxypurinol.R definition (which pools thiazide + loop + spironolactone into the same CONMED_DIUR column): Wright excludes potassium-sparing diuretics on the documented clinical rationale that thiazide and loop diuretics raise serum urate (anti-uricosuric) whereas potassium-sparing diuretics tend to lower it (uricosuric). Users simulating across the Wright 2016 and Stocker 2012 models must populate the column accordingly per the paper definition. Multiplicative effects in Wright 2016: CL/F_oxy *= 0.740^CONMED_DIUR (-26% on diuretic), U0 *= 1.14^CONMED_DIUR (+14% on diuretic); both from Wright 2016 Table 3 final-model thetadiuretic and thetaE0_diuretic.",
      source_name        = "diuretic"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 133L,
    n_studies       = 5L,
    age_range       = "27 to 83 years",
    age_median      = "60 years",
    weight_range    = "51 to 171 kg total body weight",
    weight_median   = "94 kg total body weight",
    sex_female_pct  = 12.7,
    race_ethnicity  = "European 78.2%, Maori or Pacific Islander 20.3%, East Asian 0.8% (one Korean subject), South Asian 0.8% (one Indian subject) (Wright 2016 Table 2 totals).",
    disease_state   = "Adults with gout; 29 patients were allopurinol-naive (initiating therapy) and 104 were on chronic allopurinol for >=1 month at study entry. Baseline urate range 0.18 to 0.89 mmol/L (median 0.38 mmol/L) across the pooled five-study cohort.",
    dose_range      = "Allopurinol 50-700 mg/day oral (median 300 mg/day at study entry; Wright 2016 Table 2 totals).",
    regions         = "New Zealand (multi-centre; Christchurch and Auckland).",
    crcl_range      = "12 to 125 mL/min total cohort range (median 68 mL/min ~= 4.08 L/h standardised to 70 kg; Wright 2016 Table 2).",
    co_medication   = "Concomitant medications tested as PD covariates and retained: thiazide or loop diuretic (33%). Tested without retention: beta-adrenoceptor blockers (41%), ACE inhibitors (39%), ARBs (10%), calcium-channel blockers (20%), statins (41%), NSAIDs (13%), uricosurics (3%). On PK, frusemide and probenecid were tested in the upstream Wright 2013 oxypurinol popPK (reference [21]) which informed model selection.",
    genotypes_tested = "Renal-urate-transporter SNPs rs11942223 (SLC2A9), rs2231142 (ABCG2), rs1183201 (NPT1/SLC17A1), and rs3825018 (URAT1) were tested. None were retained in the final model. An ABCG2-T-allele homozygote effect on C50 (+40%) approached significance but did not meet the chi-squared P<0.01 backward-elimination criterion (Wright 2016 Discussion).",
    samples         = "1105 oxypurinol and 1162 urate plasma concentrations from 133 gout patients (Wright 2016 Results). Concentrations >LLOQ in this analysis.",
    studies_notes   = "Five-study pooled analysis: Study 1 (n=74; Stamp 2011 [19], the largest single cohort), Study 2 (n=10; Stamp 2012 [23] frusemide-interaction study), Study 3 (n=30; Stamp 2013 [24] dose-escalation study), Study 4 (n=19; the unpublished low-dose study used for graphical analysis of urate response time delays), Study 5 (n=8; the second unpublished Christchurch study). Per-study urate additive residual SDs differ from 0.021 to 0.054 mmol/L; see addSd_eurate in ini() for the representative-cohort choice.",
    notes           = "Allopurinol (parent) is NOT modelled because (a) the allopurinol-parent model was found to be unstable in prior work (Wright 2013 [21]), (b) allopurinol has a short half-life and most clinical samples are below the assay limit of quantitation, and (c) allopurinol does not contribute meaningfully to the urate-lowering effect compared with oxypurinol (Wright 2016 PK models section). Final model selection used IPP-framework PKPD; bootstrap n=1000 (5% nonconvergence retained); base/final/PPPD/simultaneous parameter estimates agreed within ~10% (Wright 2016 Supplementary Table S3)."
  )

  ini({
    # ----- Structural PK parameters (Wright 2016 Table 3 'Final model'). -----
    # All structural PK parameters are APPARENT (CL/F_oxy, V/F_oxy):
    # the dose is allopurinol mg and only oxypurinol concentrations were
    # measured, so the bioavailability F_oxy = F_abs * f_metabolism
    # (allopurinol -> oxypurinol, molar 1:1) is absorbed into the
    # apparent parameters. The model() block converts mg of allopurinol
    # in the central compartment to umol/L of oxypurinol using the
    # allopurinol molecular weight 136.11 g/mol (the dose's molar
    # form).
    lka <- fixed(log(1.09));  label("Apparent first-order absorption rate constant Ka (1/h, fixed)")        # Wright 2016 Table 3: Ka = 1.09 fixed (both base and final columns)
    lcl <- log(1.32);         label("Apparent oral CL/F_oxy at reference FFM = 70 kg, CRCL = 6 L/h, no diuretic (L/h)") # Wright 2016 Table 3: thetaCL = 1.32 (RSE 3.9%); bootstrap median 1.31 [1.22, 1.44]
    lvc <- log(41.6);         label("Apparent oral V/F_oxy at reference TBW = 70 kg (L)")                   # Wright 2016 Table 3: thetaV = 41.6 (RSE 3.0%); bootstrap median 41.5 [39.5, 43.5]

    # Allometric exponents -- fixed per Wright 2016 Methods ('Clearance
    # was allometrically scaled to an exponent of 0.75 and volume with
    # an exponent of 1') citing Anderson & Holford [37]. Not reported
    # with RSE in Table 3, consistent with the canonical fixed-exponent
    # convention.
    e_ffm_cl <- fixed(0.75);  label("Allometric exponent of FFM on CL/F_oxy (unitless, fixed)")             # Wright 2016 Methods: CL scaled to FFM with exponent 0.75 fixed (Anderson Holford 2008/2009)
    e_wt_vc  <- fixed(1.0);   label("Allometric exponent of TBW on V/F_oxy (unitless, fixed)")              # Wright 2016 Methods: V scaled to TBW with exponent 1 fixed (Anderson Holford 2008/2009)

    # Covariate effects on CL/F_oxy.
    e_crcl_cl        <- 0.587; label("Power exponent of (CRCL/6) on CL/F_oxy (unitless)")                    # Wright 2016 Table 3: thetaRFexp = 0.587 (RSE 11.7%); bootstrap median 0.588 [0.476, 0.742]
    e_conmed_diur_cl <- 0.740; label("Multiplicative factor on CL/F_oxy when CONMED_DIUR = 1 (unitless)")    # Wright 2016 Table 3: thetadiuretic = 0.740 (RSE 6.4%); bootstrap median 0.748 [0.64, 0.86]; CL reduced by 26% on diuretics

    # ----- PD parameters (Wright 2016 Table 3 'Final model'). -----
    # Direct-effect sigmoidal Emax inhibition of urate production on
    # top of a baseline urate U0. A turnover (indirect-response) model
    # for urate did not improve the fit and was unstable (Wright 2016
    # Results: 'A turnover model for urate did not provide a better
    # description of the data and was unstable.').
    lemax  <- log(0.409); label("Maximum fractional inhibition of urate by oxypurinol (Emax, unitless)")    # Wright 2016 Table 3: Emax = 0.409 (RSE 12%); bootstrap median 0.414 [0.323, 0.595]
    lec50  <- log(83.9);  label("Oxypurinol concentration at half-maximum urate inhibition (C50, umol/L)")  # Wright 2016 Table 3: C50 = 83.9 (RSE 17.4%); bootstrap median 87.9 [61.3, 173]
    lhill  <- log(1.30);  label("Empirical Hill coefficient lambda (unitless)")                              # Wright 2016 Table 3: lambda = 1.30 (RSE 11%); bootstrap median 1.26 [1.05, 1.59]
    lrbase <- log(0.511); label("Baseline plasma urate U0 at reference CRCL = 6 L/h, no diuretic (mmol/L)") # Wright 2016 Table 3: U0 = 0.511 (RSE 2.3%); bootstrap median 0.508 [0.487, 0.530]

    # Covariate effects on baseline urate U0.
    e_crcl_rbase        <- -0.119; label("Power exponent of (CRCL/6) on baseline urate U0 (unitless)")        # Wright 2016 Table 3: thetaE0_RFexp = -0.119 (RSE 21.6%); bootstrap median -0.121 [-0.18, -0.07]
    e_conmed_diur_rbase <-  1.14;  label("Multiplicative factor on baseline urate U0 when CONMED_DIUR = 1 (unitless)") # Wright 2016 Table 3: thetaE0_diuretic = 1.14 (RSE 1.8%); bootstrap median 1.14 [1.09, 1.19]; U0 14% higher on diuretics

    # ----- IIV (Wright 2016 Table 3 'Final model'). -----
    # Wright 2016 reports between-subject variability as 'omega (CV%)',
    # which is the log-scale standard deviation expressed as a
    # percentage (the back-transformed approximation CV ~= omega for
    # small omega, common in the NONMEM tradition). This convention is
    # required because the reported PD covariance block is positive-
    # definite only under this interpretation; the alternative
    # omega^2 = log(1 + CV^2) interpretation yields a non-positive-
    # definite block for (Emax, U0, C50). The encoding used here is:
    #   omega^2 = (CV% / 100)^2
    # which equates the reported '%CV' to the log-scale SD directly.

    # CL/F_oxy IIV. The paper's '24.2%' is treated as log-scale SD =
    # 0.242, giving omega^2 = 0.0586.
    etalcl ~ 0.0586                                                                                          # Wright 2016 Table 3: omegaCL = 24.2% (RSE 15.5%); bootstrap median 24.0 [19.3, 28.3]; encoded as omega^2 = 0.242^2

    # V/F_oxy fractional-effect scaler. Wright 2016 Results:
    # 'only a single random effects parameter was estimated for both CL
    # and V and the variance of V estimated using a fractional effect
    # parameter [50]' (Bonate 2006). F_v_oxy is fixed at 0.0355; V's
    # individual log-scale deviation is F_v_oxy * etalcl (perfect
    # correlation between CL and V on the log scale, magnitude scaled
    # by F).
    F_v_oxy <- fixed(0.0355); label("Fractional shared-eta scaler from etalcl onto V/F_oxy (unitless, fixed)") # Wright 2016 Table 3: F_omegaV_oxy = 0.0355 fixed (both base and final)

    # Ka IIV. Wright 2016 Table 3: omega_Ka = 58.9% (fixed in both
    # base and final). Encoded as omega^2 = 0.589^2 = 0.347.
    etalka ~ fixed(0.347)                                                                                    # Wright 2016 Table 3: omegaKa = 58.9% fixed; encoded as omega^2 = 0.589^2

    # PD block IIV (Wright 2016 Table 3) with reported off-diagonal
    # covariances on the log scale. Order in the block matrix: lrbase
    # (= log U0), lemax, lec50. Variances (with reported %CV as log-
    # scale SD):
    #   var_lrbase = 0.142^2 = 0.0202
    #   var_lemax  = 0.359^2 = 0.1289
    #   var_lec50  = 0.607^2 = 0.3684
    # Reported covariances (Wright 2016 Table 3 'Final model'):
    #   Cov(eta_Emax, eta_U0)   = 0.025
    #   Cov(eta_Emax, eta_C50)  = 0.193
    #   Cov(eta_U0,   eta_C50)  = 0.011
    # Resulting Pearson correlations: r(Emax,U0)=0.49, r(Emax,C50)=0.89,
    # r(U0,C50)=0.13. The block is positive-definite under the encoding
    # above.
    etalrbase + etalemax + etalec50 ~ c(0.0202,
                                        0.025,    0.1289,
                                        0.011,    0.193,    0.3684)

    # ----- Residual error (Wright 2016 Table 3). -----
    # Oxypurinol Cc: proportional 19.9% CV plus a tiny additive fixed
    # at 0.001 umol/L per the Table 3 footnote 'Oxypurinol sigma_add
    # fixed at 0.001'. The additive is effectively negligible.
    propSd <- 0.199;          label("Proportional residual error on oxypurinol Cc (fraction)")              # Wright 2016 Table 3: oxypurinol sigma_prop = 19.9% (final); bootstrap median 19.9% [17.7, 22.2]
    addSd  <- fixed(0.001);   label("Additive residual error on oxypurinol Cc (umol/L, fixed)")             # Wright 2016 Table 3 footer: 'Oxypurinol sigma_add fixed at 0.001'

    # Urate Eurate: additive only per the Table 3 footer 'Urate
    # sigma_prop fixed at 0.001'. The per-study additive SD varied
    # 0.021 to 0.054 mmol/L across the five pooled studies:
    #   Study 1 (n = 74, Stamp 2011 [19]) sigma_add = 0.037 mmol/L
    #   Study 2 (n = 10, Stamp 2012 [23]) sigma_add = 0.022 mmol/L
    #   Study 3 (n = 30, Stamp 2013 [24]) sigma_add = 0.037 mmol/L
    #   Study 4 (n = 19, unpublished)     sigma_add = 0.021 mmol/L
    #   Study 5 (n =  8, unpublished)     sigma_add = 0.054 mmol/L
    # The single representative value retained here is 0.037 mmol/L
    # (Studies 1 and 3), which is the cohort median and the value
    # associated with the largest sub-cohort (Study 1, n = 74 of 133
    # = 55.6%). Documented in this file's Errata-equivalent population$
    # studies_notes and re-stated in the validation vignette.
    addSd_Eurate  <- 0.037;         label("Additive residual error on plasma urate Eurate (mmol/L; Study 1 / Study 3 representative)") # Wright 2016 Table 3: Study 1 urate sigma_add = 0.037 (RSE 6.2%); bootstrap median 0.038 [0.033, 0.044]
    propSd_Eurate <- fixed(0.001);  label("Proportional residual error on plasma urate Eurate (fraction, fixed)") # Wright 2016 Table 3 footer: 'Urate sigma_prop fixed at 0.001'
  })

  model({
    # Reference values used by the final covariate model (Wright 2016
    # Table 3 footnotes: clearance expressed per 70 kg FFM / CLcr 6 L/h;
    # volume expressed per 70 kg TBW). CRCL units are L/h
    # (standardised to 70 kg per the paper's covariateData notes).
    ref_ffm  <- 70
    ref_wt   <- 70
    ref_crcl <- 6

    # Allopurinol molecular weight (g/mol = mg/mmol). Used to convert
    # mg of allopurinol (the dose unit) in the central compartment to
    # umol/L of oxypurinol (the observed concentration unit, after the
    # implicit F_oxy absorbed into the apparent CL/F_oxy and V/F_oxy).
    MW_allop <- 136.11

    # ----- Individual PK parameters -----
    # Wright 2016 final-model equations:
    #   CL/F_oxy = thetaCL * (CRCL/6)^thetaRFexp * (FFM/70)^0.75
    #                       * thetadiuretic^CONMED_DIUR
    #   V/F_oxy  = thetaV  * (TBW/70)^1
    #   Ka       = 1.09 (fixed)
    # Per Wright 2016 Results, V's eta is a fractional shared scaling
    # of CL's eta: V_i = V_pop * exp(F_v_oxy * eta_CL).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
          (CRCL / ref_crcl)^e_crcl_cl *
          (FFM / ref_ffm)^e_ffm_cl *
          e_conmed_diur_cl^CONMED_DIUR
    vc <- exp(lvc + F_v_oxy * etalcl) *
          (WT / ref_wt)^e_wt_vc
    kel <- cl / vc

    # ----- Individual PD parameters -----
    emax_pd <- exp(lemax + etalemax)
    ec50_pd <- exp(lec50 + etalec50)
    hill    <- exp(lhill)

    # Baseline urate U0 with renal-function power scaling and diuretic
    # multiplier:
    #   U0 = thetaU0 * (CRCL/6)^thetaE0_RFexp * thetaE0_diuretic^CONMED_DIUR
    rbase <- exp(lrbase + etalrbase) *
             (CRCL / ref_crcl)^e_crcl_rbase *
             e_conmed_diur_rbase^CONMED_DIUR

    # ----- ODE system -----
    # depot and central carry mg of allopurinol-equivalent mass (the
    # implicit F_oxy fraction and molar 1:1 conversion to oxypurinol
    # are absorbed into the apparent CL/F_oxy and V/F_oxy). Cc below
    # converts mg/L in central to umol/L of oxypurinol using the
    # allopurinol molecular weight (the dose's molar form).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- Observations -----
    # Oxypurinol concentration in umol/L. Conversion factor
    # 1000 / 136.11 ~= 7.347 umol/mg (umol/L per mg/L).
    Cc <- central / vc * 1000 / MW_allop

    # Direct-effect sigmoidal Emax inhibition of urate (Wright 2016
    # equation set):
    #   Eurate = U0 * (1 - Emax * Cc^lambda / (C50^lambda + Cc^lambda))
    Eurate <- rbase * (1 - emax_pd * Cc^hill / (ec50_pd^hill + Cc^hill))

    # Residual error: combined add + prop on Cc (Wright 2016 Table 3,
    # oxypurinol additive fixed); additive on Eurate (Wright 2016
    # Table 3, urate proportional fixed). The model has two algebraic
    # observables (Cc, Eurate) backed by the same single ODE state
    # (central). In event tables, route observation rows to one of
    # the two named observables via `cmt = "Cc"` or `cmt = "Eurate"`;
    # rxode2 returns BOTH algebraic columns at every observation row
    # regardless of which one the row was routed to. Only two ODE
    # states (depot, central) exist, so the auto-injected observable
    # cmt slots do not renumber any state reference (no slot-
    # renumbering bug).
    Cc     ~ add(addSd)        + prop(propSd)
    Eurate ~ add(addSd_Eurate) + prop(propSd_Eurate)
  })
}
