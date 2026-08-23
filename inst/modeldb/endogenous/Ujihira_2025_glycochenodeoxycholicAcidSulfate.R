Ujihira_2025_glycochenodeoxycholicAcidSulfate <- function() {
  description <- "Coupled three-drug turnover model for the endogenous OATP1B3 / OAT3 biomarker glycochenodeoxycholic acid 3-O-sulfate (GCDCA-S) in healthy adults (Ujihira 2025), fit simultaneously with the perpetrator PK of rifampicin and probenecid. GCDCA-S is produced at a zero-order synthesis rate ksyn and eliminated by hepatobiliary clearance CLh (the dominant route, ~98% of total systemic clearance, giving fe < 5%) and renal clearance CLR, with plasma concentration and cumulative urinary amount as the two biomarker outputs. Rifampicin competitively inhibits CLh through the unbound OATP1B3 constant Ki,u,OATP1B3 driven by its own one-compartment zero-order-absorption PK; probenecid competitively inhibits CLR through the unbound OAT3 constant Ki,u,OAT3 driven by its own one-compartment first-order-absorption PK, and additionally reduces CLh by a concentration-independent 1.7-fold factor X gated by the binary CONMED_PROBENECID indicator. All three drugs live in this one file because the paper fit them as a single coupled system (Figure 2); dosing rifampicin or probenecid alone drives the corresponding interaction, and dosing neither collapses the model to the inhibitor-free steady-state baseline of about 65 nmol/L (typical value). Perpetrator doses are given in umol, not mg, because the whole system runs in umol / umol per L."
  reference <- paste(
    "Ujihira Y, Georgiev V, Ogungbenro K, Galetin A.",
    "Population Pharmacokinetic Modeling of Glycochenodeoxycholic Acid",
    "3-O-Sulfate (GCDCA-S) as Endogenous Biomarker of OATP1B3 and OAT3",
    "Transporters.",
    "Clin Pharmacol Ther. 2025;118(6):1532-1543.",
    "doi:10.1002/cpt.70023.",
    "Structural equations from Eqs 1-4 of the main text; all parameter",
    "values from Table 2. Supplementary Material S1 (Tables S1-S7) was",
    "also used and is on disk. The unbound fractions that convert each",
    "perpetrator's modelled total plasma concentration to the unbound",
    "concentration driving inhibition are reported nowhere in the paper",
    "or its supplement; they are taken from the two upstream models the",
    "paper adopts -- rifampicin fu = 0.11 (Barnett 2018 Clin Pharmacol",
    "Ther 104:564-574 Table 1 footnote c) and probenecid fu = 0.062",
    "(Ahmad 2021 CPT Pharmacometrics Syst Pharmacol 10:467-477) -- and",
    "each is pinned by an on-disk answer key: the Table S4 AUCR of 13",
    "at 600 mg rifampicin and the Table S4 renal-clearance ratio of 0.1",
    "under 500 mg QID probenecid. See the vignette Errata.",
    sep = " "
  )
  vignette <- "Ujihira_2025_glycochenodeoxycholicAcidSulfate"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Checked against Ujihira 2025 Eqs 1-4 and Figure 2.
  compartmentData <- list(
    central      = list(analyte = "glycochenodeoxycholic acid 3-O-sulfate", units = "umol", specimen = "plasma", verified = TRUE),
    urine        = list(analyte = "glycochenodeoxycholic acid 3-O-sulfate", units = "umol", specimen = "urine",  verified = TRUE),
    central_rif  = list(analyte = "rifampicin",  units = "umol", specimen = "plasma", verified = TRUE),
    depot_prob   = list(analyte = "probenecid",  units = "umol", specimen = "administration site", verified = TRUE),
    central_prob = list(analyte = "probenecid",  units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_PROBENECID = list(
      description        = "Concomitant probenecid co-administration indicator (1 = subject is within a probenecid dosing period; 0 = no probenecid). Gates the concentration-independent fold reduction X in GCDCA-S hepatobiliary clearance that Ujihira 2025 estimated for probenecid's weak OATP1B3 inhibition.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no probenecid co-administration)",
      notes              = "Time-varying within subject. Probenecid regimen in the source studies (Ujihira 2025 Table 1, studies #2 and #3, adopted from Willemin et al. 2021): 500 mg orally at 6 pm and 11 pm on day 0, then 500 mg at 7 am, 1 pm, 6 pm and 11 pm on days 1 to 7; GCDCA-S plasma and urine sampling begins on day 1. The indicator modifies ONLY the hepatobiliary clearance arm cl_nonren, via the multiplicative factor (1 + e_prob_clh * CONMED_PROBENECID) = 1 / X. It deliberately does NOT gate the renal-clearance inhibition, which is concentration-driven through the probenecid PK model and Ki,u,OAT3 and therefore switches itself off whenever probenecid is not dosed. Set to 1 from the first probenecid dose through the end of the probenecid-phase sampling window; the paper reports no lag or washout term. Rifampicin needs no equivalent indicator because its only effect is the concentration-driven OATP1B3 term."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 3L,
    age_range      = "20-71 years across the three model-development studies (study #1 50-71; study #2 48-58; study #3 20-55).",
    weight_range   = "(not extracted; Ujihira 2025 Table 1 does not tabulate body weight, and no allometric or weight covariate appears in the final model.)",
    sex_female_pct = 37.5,
    race_ethnicity = "White in all three model-development studies (Ujihira 2025 Table 1).",
    disease_state  = "Healthy volunteers. GCDCA-S was monitored as an endogenous biomarker of hepatic OATP1B3 and renal OAT3 activity during transporter drug-drug-interaction studies.",
    dose_range     = "Endogenous biomarker (no exogenous GCDCA-S dose). Perpetrators: rifampicin single 600 mg oral dose (study #1); probenecid 500 mg orally four times daily (QID) on days 1 to 7, preceded by two loading doses on day 0 (studies #2 and #3).",
    regions        = "(not extracted; Ujihira 2025 Table 1 reports ethnicity but not study region for the three development studies.)",
    notes          = "Pooled from three healthy-volunteer crossover studies, each with a control occasion and an inhibitor occasion: study #1 Tatosian et al. 2021 (n = 6; 3 male, 3 female; microdose probe cocktail +/- rifampicin), study #2 Willemin et al. 2021 (n = 6; 6 female; +/- probenecid), study #3 unpublished data of the same design as study #2 (n = 12; 12 male). Sex percentage is 9 female of 24. 430 GCDCA-S plasma samples and 175 GCDCA-S urine samples were fit simultaneously with 153 perpetrator plasma samples in Monolix 2024R2. Observed GCDCA-S plasma baseline was approximately 80 nmol/L with high between-subject (CV 34-69%) and within-subject diurnal (CV 34-43%) variability. Neither a diurnal-fluctuation function (six forms tested, Table S2) nor a sex effect on ksyn (Table S3) improved the fit, so neither is in the final model. The model was verified against four independent studies (Table S1) not used in fitting."
  )

  # Covariates screened by Ujihira 2025 but NOT retained in the final model.
  # Documentation only: checkModelConventions() does not require these to be
  # referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on the GCDCA-S synthesis rate ksyn as ksyn,sex = ksyn,male * (1 - SEX * COVSEX) (Supplementary Material S1 Eq. 2, with SEX = 1 for women). The estimate COVSEX = 0.065 carried an RSE of 404% and moved the objective function by only 3 points (OFV -642 -> -639, p > 0.9; Table S3), so sex was excluded from the final model despite the paper's observation that baseline GCDCA-S is higher in men."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # GCDCA-S biomarker -- Ujihira 2025 Table 2, GCDCA-S block.
    # The turnover model is Eq. 1 (no inhibitor) / Eq. 2 (with
    # inhibitor) for plasma and Eq. 3 / Eq. 4 for cumulative urine.
    # ---------------------------------------------------------------
    lksyn <- log(1.0)
    label("GCDCA-S zero-order synthesis rate ksyn (umol/h)")
    # Table 2, GCDCA-S row 'k syn (umol/h)' = 1.0 (RSE 18%).

    lvc <- fixed(log(4.8))
    label("GCDCA-S volume of distribution Vc (L)")
    # Table 2, GCDCA-S row 'V c,GCDCA-S (L)' = 4.8 (FIXED). The paper
    # fixed Vc for model identifiability at the value initially
    # estimated by the POP-PK model.

    lcl_renal <- log(0.31)
    label("GCDCA-S renal clearance CLR (L/h)")
    # Table 2, GCDCA-S row 'CL R,GCDCA-S (L/h)' = 0.31 (RSE 8%).

    lcl_nonren <- log(15)
    label("GCDCA-S hepatobiliary clearance CLh (L/h)")
    # Table 2, GCDCA-S row 'CL h,GCDCA-S (L/h)' = 15 (RSE 13%). The
    # non-renal arm of the CL_total = CL_renal + CL_nonren
    # decomposition; the paper's CLh is entirely hepatobiliary, and
    # CLh / (CLh + CLR) = 15 / 15.31 = 98% reproduces the paper's
    # "biliary excretion is the primary route (~95%), fe < 5%".

    # ---------------------------------------------------------------
    # Rifampicin perpetrator PK -- Table 2, RIF block. One compartment,
    # zero-order absorption of duration Tk0 after a lag Tlag, linear
    # elimination (Figure 2, left panel). Structure adopted from
    # Takita et al. 2022 Clin Pharmacol Ther 112:615-626.
    # ---------------------------------------------------------------
    ltlag_rif <- log(0.56)
    label("Rifampicin absorption lag time Tlag,RIF (h)")
    # Table 2, RIF row 'Tlag RIF (h)' = 0.56 (RSE 20%).

    ld1_rif <- log(0.84)
    label("Rifampicin zero-order absorption duration Tk0,RIF (h)")
    # Table 2, RIF row 'Tk0 RIF (h)' = 0.84 (RSE 41%).

    lvc_rif <- log(35)
    label("Rifampicin volume of distribution V,RIF (L)")
    # Table 2, RIF row 'V RIF (L)' = 35 (RSE 5%).

    lcl_rif <- log(6.7)
    label("Rifampicin clearance CL,RIF (L/h)")
    # Table 2, RIF row 'CL RIF (L)' = 6.7 (RSE 8%). The unit is printed
    # as '(L)' in Table 2; it must be L/h for Eq. 1 to balance, and
    # 6.7 L/h matches the 6.35 L/h of the upstream Takita 2022 model.

    lkiu_oatp1b3 <- log(0.009)
    label("Rifampicin unbound OATP1B3 inhibition constant Ki,u,OATP1B3 (umol/L)")
    # Table 2, RIF row 'K i,unbound,OATP1B3 (uM)' = 0.009 (RSE 29%).

    fu_rif <- fixed(0.11)
    label("Rifampicin fraction unbound in plasma (unitless)")
    # NOT reported by Ujihira 2025. Taken from the lineage this model
    # is built on and confirmed three ways: (i) Barnett 2018 Table 1
    # footnote c states Ki,total 1.15 uM -> Ki,unbound 0.13 uM using
    # RIF fu = 0.11, and Ujihira Table S7 lists that same study at
    # 0.13; (ii) Barnett 2018 Table 1 rosuvastatin row, Ki,total
    # 2.23 -> Ki,unbound 0.25, same fu; (iii) Takita 2022 Table 1
    # reports Ki,total = 0.345 uM and Ujihira Table S7 lists that same
    # study at Ki,u = 0.038 uM, i.e. 0.038 / 0.345 = 0.110. See the
    # vignette Errata.

    # ---------------------------------------------------------------
    # Probenecid perpetrator PK -- Table 2, PROB block. One compartment
    # with first-order absorption and linear elimination (Figure 2,
    # right panel). Structure adopted from Ahmad et al. 2021 CPT
    # Pharmacometrics Syst Pharmacol 10:467-477.
    # ---------------------------------------------------------------
    lka_prob <- log(1.1)
    label("Probenecid absorption rate constant ka,PROB (1/h)")
    # Table 2, PROB row 'ka PROB (/h)' = 1.1 (RSE 37%).

    lvc_prob <- log(17)
    label("Probenecid volume of distribution V,PROB (L)")
    # Table 2, PROB row 'V PROB (L)' = 17 (RSE 10%).

    lcl_prob <- log(0.95)
    label("Probenecid clearance CL,PROB (L/h)")
    # Table 2, PROB row 'CL PROB (L/h)' = 0.95 (RSE 7%).

    lkiu_oat3 <- log(2.7)
    label("Probenecid unbound OAT3 inhibition constant Ki,u,OAT3 (umol/L)")
    # Table 2, PROB row 'K i,unbound,OAT3 (uM)' = 2.7 (RSE 12%).

    fu_prob <- fixed(0.062)
    label("Probenecid fraction unbound in plasma (unitless)")
    # NOT reported by Ujihira 2025 anywhere in the main text or the
    # supplement, and NOT independently verifiable on disk: it is
    # attributed to Ahmad et al. 2021, the paper Ujihira 2025 adopts
    # its probenecid structural model from (Methods, "Structural
    # models"; supplement reference 9), which is not among this
    # extraction's source files. It is instead PINNED by an on-disk
    # answer key: Ujihira 2025 Table S4 reports a GCDCA-S renal-
    # clearance ratio CL_R,+inhibitor / CL_R,control of 0.1 at steady
    # state under the 500 mg QID probenecid regimen. Reproducing 0.1
    # requires 1 + Cu,PROB/Ki,u,OAT3 ~ 10, i.e. Cu,PROB ~ 24 uM. With
    # V = 17 L and CL = 0.95 L/h the total average steady-state
    # concentration is 1752 umol / (0.95 L/h * 6 h) = 307 uM, so
    # fu = 0.062 gives Cu = 19 uM and a ratio of 0.12, consistent with
    # 0.1 at one significant figure. fu = 1 would give 0.009 and
    # fu = 0.03 would give 0.23; both are excluded. The vignette runs
    # this check as an assertion. See the vignette Errata.

    e_prob_clh <- -0.4117647
    label("Multiplicative change in GCDCA-S hepatobiliary clearance during probenecid co-administration (unitless)")
    # Table 2, PROB row 'X' = 1.7 (RSE 5%), defined by the paper as the
    # fold DECREASE in GCDCA-S CLh in the presence of probenecid, with
    # X = 1 in its absence (Eq. 2 carries CLh * 1/X). Encoded in the
    # standard binary-covariate multiplicative shape
    # (1 + e_prob_clh * CONMED_PROBENECID) so that the factor is 1 with
    # no probenecid and 1/X with probenecid: e = 1/1.7 - 1 = -0.4117647.
    # The paper declined to estimate a probenecid OATP1B3 Ki because the
    # inhibition was incomplete and no dose-range data were available.

    # ---------------------------------------------------------------
    # Between-subject variability. Ujihira 2025 Table 2 column 'IIV'
    # reports a percentage with its own RSE and states that random
    # effects are log-normally distributed; the values are read as
    # CV% and converted with omega^2 = log(1 + CV^2). Parameters whose
    # IIV cell is '-' were not given a random effect by the paper and
    # deliberately carry none here.
    # ---------------------------------------------------------------
    etalksyn      ~ 0.3163576   # Table 2, GCDCA-S IIV ksyn = 61% (RSE 15%); log(1 + 0.61^2)
    etalcl_renal  ~ 0.0560025   # Table 2, GCDCA-S IIV CLR  = 24% (RSE 25%); log(1 + 0.24^2)
    etaltlag_rif  ~ 0.1919877   # Table 2, RIF IIV Tlag     = 46% (RSE 32%); log(1 + 0.46^2)
    etald1_rif    ~ 0.4656585   # Table 2, RIF IIV Tk0      = 77% (RSE 37%); log(1 + 0.77^2)
    etalcl_rif    ~ 0.0318856   # Table 2, RIF IIV CL       = 18% (RSE 35%); log(1 + 0.18^2)
    etalvc_prob   ~ 0.0917476   # Table 2, PROB IIV V       = 31% (RSE 28%); log(1 + 0.31^2)
    etalcl_prob   ~ 0.0703679   # Table 2, PROB IIV CL      = 27% (RSE 20%); log(1 + 0.27^2)
    etalkiu_oat3  ~ 0.0974856   # Table 2, PROB IIV Ki,u,OAT3 = 32% (RSE 27%); log(1 + 0.32^2)

    # ---------------------------------------------------------------
    # Residual unexplained variability -- Table 2, 'Residual
    # unexplained variabilities' rows. Both GCDCA-S outputs are
    # proportional only; both perpetrators are combined
    # proportional + additive with the additive term FIXED.
    # ---------------------------------------------------------------
    propSd <- 0.45
    label("Proportional residual error, GCDCA-S plasma (fraction)")
    # Table 2, 'sigma prop - GCDCA-S plasma (%)' = 45 (RSE 3%).

    propSd_Ugcdcas <- 0.29
    label("Proportional residual error, GCDCA-S urine (fraction)")
    # Table 2, 'sigma prop - GCDCA-S urine (%)' = 29 (RSE 1%).

    propSd_rif <- 0.20
    label("Proportional residual error, rifampicin plasma (fraction)")
    # Table 2, 'sigma prop - RIF (%)' = 20 (RSE 14%).

    addSd_rif <- fixed(0.2)
    label("Additive residual error, rifampicin plasma (umol/L)")
    # Table 2, 'sigma add - RIF (uM)' = 0.2 (FIXED).

    propSd_prob <- 0.23
    label("Proportional residual error, probenecid plasma (fraction)")
    # Table 2, 'sigma prop - PROB (%)' = 23 (RSE 8%).

    addSd_prob <- fixed(0.001)
    label("Additive residual error, probenecid plasma (umol/L)")
    # Table 2, 'sigma add - PROB (uM)' = 0.001 (FIXED).
  })

  model({
    # --- Individual parameters -------------------------------------
    ksyn       <- exp(lksyn + etalksyn)
    vc         <- exp(lvc)
    cl_renal   <- exp(lcl_renal + etalcl_renal)
    cl_nonren  <- exp(lcl_nonren)

    tlag_rif   <- exp(ltlag_rif + etaltlag_rif)
    d1_rif     <- exp(ld1_rif + etald1_rif)
    vc_rif     <- exp(lvc_rif)
    cl_rif     <- exp(lcl_rif + etalcl_rif)
    kiu_oatp1b3 <- exp(lkiu_oatp1b3)

    ka_prob    <- exp(lka_prob)
    vc_prob    <- exp(lvc_prob + etalvc_prob)
    cl_prob    <- exp(lcl_prob + etalcl_prob)
    kiu_oat3   <- exp(lkiu_oat3 + etalkiu_oat3)

    # --- Perpetrator PK --------------------------------------------
    # Rifampicin: zero-order input of duration d1_rif into central_rif
    # after a lag. Dose records must carry rate = -2 so rxode2 uses the
    # modelled duration.
    dur(central_rif)  <- d1_rif
    alag(central_rif) <- tlag_rif

    Cc_rif  <- central_rif / vc_rif
    Cc_prob <- central_prob / vc_prob

    # Unbound perpetrator concentrations drive the inhibition terms.
    # Both collapse to 0 when the corresponding drug is not dosed.
    cu_rif  <- fu_rif * Cc_rif
    cu_prob <- fu_prob * Cc_prob

    # --- GCDCA-S turnover ------------------------------------------
    # Ujihira 2025 Eq. 2: CLh is divided by the competitive OATP1B3
    # term driven by unbound rifampicin AND by the probenecid fold
    # reduction X; Eq. 4: CLR is divided by the competitive OAT3 term
    # driven by unbound probenecid. With no perpetrator on board both
    # effective clearances reduce to their Eq. 1 / Eq. 3 baselines.
    cl_nonren_eff <- cl_nonren / (1 + cu_rif / kiu_oatp1b3) *
      (1 + e_prob_clh * CONMED_PROBENECID)
    cl_renal_eff  <- cl_renal / (1 + cu_prob / kiu_oat3)

    Cc <- central / vc

    # ODE declaration order sets compartment numbering: the GCDCA-S
    # biomarker states come first (central = 1, urine = 2) so the
    # subject of the model is compartment 1, followed by the
    # perpetrator states (central_rif = 3, depot_prob = 4,
    # central_prob = 5).
    d/dt(central)      <- ksyn - (cl_nonren_eff + cl_renal_eff) * Cc
    d/dt(urine)        <- cl_renal_eff * Cc
    d/dt(central_rif)  <- -cl_rif / vc_rif * central_rif
    d/dt(depot_prob)   <- -ka_prob * depot_prob
    d/dt(central_prob) <- ka_prob * depot_prob - cl_prob / vc_prob * central_prob

    # Steady-state baseline before any perpetrator is given, from
    # Eq. 1 with dC/dt = 0: Css = ksyn / (CLh + CLR). Typical value
    # 1.0 / (15 + 0.31) = 0.0653 umol/L = 65 nmol/L, which with the
    # 61% log-normal IIV on ksyn gives a population mean of about
    # 77 nmol/L against the roughly 80 nmol/L observed (Results).
    central(0) <- ksyn / (cl_nonren + cl_renal) * vc
    urine(0)   <- 0

    # --- Outputs ----------------------------------------------------
    Ugcdcas <- urine

    Cc      ~ prop(propSd)
    Ugcdcas ~ prop(propSd_Ugcdcas)
    Cc_rif  ~ add(addSd_rif)  + prop(propSd_rif)
    Cc_prob ~ add(addSd_prob) + prop(propSd_prob)
  })
}
