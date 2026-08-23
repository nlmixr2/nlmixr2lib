Said_2025_imatinib <- function() {
  description <- "Joint parent-metabolite population PK model for total imatinib (Cc), unbound imatinib (Cu), and total N-desmethyl imatinib (Cc_ndmima) in 335 pooled adults: COVID-19 ARDS patients from the InventCOVID (IV) and CounterCOVID (oral) trials plus a historical CML/GIST oncology cohort (Said 2025). Two-compartment parent disposition with first-order oral absorption and a one-compartment metabolite; saturable 1:1 molar binding of both analytes to alpha-1-acid glycoprotein solved in closed form, with the unbound fraction driving elimination, metabolite formation, and central-to-peripheral distribution. Retained covariates: COVID-19 on central volume, and IL-6R-inhibitor cotreatment on the AAG dissociation constant and on apparent metabolite clearance."
  reference <- "Said MM, Schippers JR, Bos LDJ, Atmowihardjo L, Mathot RAA, Li Y, van der Plas MS, Aman J, Bogaard HJ, Swart EL, Bartelink IH. Disease-Drug-Drug Interaction of Imatinib in COVID-19 ARDS: A Pooled Population Pharmacokinetic Analysis. CPT Pharmacometrics Syst Pharmacol. 2025;14(3):583-595. doi:10.1002/psp4.13299"
  vignette <- "Said_2025_imatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  dosing <- c("depot", "central")

  compartmentData <- list(
    depot          = list(analyte = "imatinib",              units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "imatinib",              units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1    = list(analyte = "imatinib",              units = "mg", specimen = "plasma",              verified = TRUE),
    central_ndmima = list(analyte = "N-desmethyl imatinib",  units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    AAG = list(
      description        = "Serum alpha-1-acid glycoprotein concentration; the saturable binding protein for both imatinib and N-desmethyl imatinib in this model",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING, and structural rather than a covariate effect: AAG enters the closed-form binding",
        "solution (Said 2025 Equations 1-2) rather than shifting a typical value, so it has no reference",
        "value and no `e_aag_*` coefficient. Measured by ELISA; a median of 4 (range 1-11) AAG levels were",
        "measured per patient during treatment, and missing values were imputed by last observation carried",
        "forward (Said 2025 Methods 2.2). Cohort medians (Table 1): 1.8 g/L ward C-ARDS, 2.1 g/L ICU C-ARDS,",
        "1.6 g/L ICU C-ARDS on IL-6R inhibitor, 0.8 g/L CML/GIST. The forest-plot typical patient uses",
        "0.9 g/L (Said 2025 Figure 2). The source control stream column is AAGI2, the LOCF-imputed series.",
        sep = " "
      ),
      source_name        = "AAGI2"
    ),
    DIS_COVID19 = list(
      description        = "1 = COVID-19 patient (pooled InventCOVID and CounterCOVID C-ARDS cohorts), 0 = CML/GIST oncology patient",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CML/GIST cancer patient from the historical oncology cohort)",
      notes              = paste(
        "Time-fixed per subject. The reference complement in this analysis is the CML/GIST oncology cohort,",
        "not healthy volunteers. Said 2025 calls the COVID-19 arm C-ARDS throughout, but the flag is not an",
        "ARDS indicator: InventCOVID required moderate-to-severe ARDS by Berlin criteria while CounterCOVID",
        "enrolled hospitalised COVID-19 patients requiring supplemental oxygen, and both are coded 1. The",
        "source control stream column is COVID. DIS_ARDS was therefore NOT used.",
        sep = " "
      ),
      source_name        = "COVID"
    ),
    CONMED_IL6RI = list(
      description        = "1 = preceding cotreatment with an interleukin-6 receptor inhibitor (8 mg/kg IV tocilizumab or 400 mg IV sarilumab, single dose on ICU admission), 0 = no IL-6R inhibitor",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no IL-6R inhibitor cotreatment)",
      notes              = paste(
        "Time-fixed per subject in the source analysis. All 29 IL-6R-inhibitor patients came from",
        "InventCOVID, where 90.6% of patients received one as standard care on ICU admission; no",
        "CounterCOVID or CML/GIST patient received one (Said 2025 Table 1). The source control stream",
        "column is IL6INHIB. Retained on Kd and on apparent metabolite clearance; it was also tested on",
        "unbound imatinib clearance and NOT retained (THETA(9) is 0 FIX in the Data S1 control stream).",
        sep = " "
      ),
      source_name        = "IL6INHIB"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Entered the full covariate model on CLu/F1 and V1/F1 but eliminated: the 95% CI fell inside the 80%-120% clinical-relevance bounds of the Said 2025 Figure 2 forest plot. No point estimate is reported."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Entered the full covariate model on CLu/F1 and V1/F1 but eliminated as not clinically significant (Said 2025 Figure 2 forest plot). No point estimate is reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not clinically significant in the Said 2025 Figure 2 forest plot. No point estimate is reported."
    ),
    WHO_ORDINAL = list(
      description = "WHO 9-point ordinal scale for clinical improvement, used as the COVID-19 disease-severity covariate",
      units       = "score (0-8)",
      type        = "continuous",
      notes       = "Tested both as a continuous WHO grade and as a binary WHO grade > 5 (ICU-admitted and invasively ventilated). Entered the full covariate model on CLu/F1 and V1/F1 but eliminated as not clinically significant (Said 2025 Figure 2 forest plot). No point estimate is reported."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Said 2025 Methods 2.3.2); not retained and no point estimate reported."
    ),
    AST = list(
      description = "Serum aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Said 2025 Methods 2.3.2); not retained and no point estimate reported."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Said 2025 Methods 2.3.2); not retained and no point estimate reported."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Said 2025 Methods 2.3.2); not retained and no point estimate reported."
    ),
    IL6 = list(
      description = "Serum interleukin-6 concentration",
      units       = "pg/mL",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Said 2025 Methods 2.3.2); not retained and no point estimate reported. Not measured at all in the CML/GIST cohort."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 335L,
    n_studies      = 3L,
    n_observations = paste(
      "InventCOVID contributed 160 total imatinib, 109 unbound imatinib, and 159 total DM-imatinib samples",
      "from 32 patients. 54 CounterCOVID samples were reanalysed (27 ICU, 27 matched hospitalised). The",
      "CML/GIST cohort contributed 475 total and 150 unbound imatinib concentrations at steady state,",
      "median 4 (range 1-10) samples per patient (Said 2025 Results 3.1 and Table S1).",
      sep = " "
    ),
    age_range      = "20-93 years",
    age_median     = "59 years (CML/GIST) to 64 years (ward and ICU/IL6RINH C-ARDS)",
    weight_range   = "40-167 kg",
    weight_median  = "70 kg (CML/GIST) to 88 kg (ICU C-ARDS); the forest-plot typical patient is 80.7 kg",
    sex_female_pct = 34.0,
    race_ethnicity = "Not reported in the source paper (Dutch multicentre trials plus a historical Dutch oncology cohort).",
    disease_state  = paste(
      "Pooled across four groups (Said 2025 Table 1): (1) COVID-19 ARDS patients who remained hospitalised",
      "on the ward, N = 158; (2) COVID-19 ARDS patients admitted to the ICU and invasively ventilated,",
      "N = 42; (3) COVID-19 ARDS patients admitted to the ICU, invasively ventilated, and given an IL-6R",
      "inhibitor on admission, N = 29; (4) chronic myelogenous leukemia / gastrointestinal stromal tumor",
      "outpatients, N = 106. Patients with pre-existing chronic pulmonary or cardiac disease, on strong",
      "CYP3A4 inducers, pregnant, breastfeeding, recently treated for malignancy, or on chronic home oxygen",
      "were excluded from InventCOVID.",
      sep = " "
    ),
    dose_range     = paste(
      "InventCOVID: 200 mg imatinib as a 2 h IV infusion twice daily for up to 7 days.",
      "CounterCOVID: 800 mg oral loading dose on day 0 then 400 mg orally once daily for 9 days.",
      "CML/GIST: 100-800 mg orally once daily at steady state.",
      sep = " "
    ),
    regions        = "Netherlands",
    notes          = paste(
      "Baseline characteristics are in Said 2025 Table 1 (grouped by ward / ICU / ICU+IL6RINH / CML/GIST)",
      "and Table S1 (grouped by study: InventCOVID N = 32, CounterCOVID N = 197, CML/GIST N = 106).",
      "n_subjects = 335 is the Table 1 / Table S1 column total (229 COVID-19 + 106 CML/GIST); Said 2025",
      "Results 3.1 describes the oncology cohort as 20 CML + 85 GIST = 105 patients, one fewer than both",
      "tables report. Trial registrations NCT04794088 / EudraCT 2020-005447-23 (InventCOVID).",
      "sex_female_pct = 34.0 is computed from the Table 1 counts (114 of 335).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # All fixed effects are the FINAL estimates from Said 2025 Table 2.
    # The Data S1 control stream's $THETA block holds INITIAL estimates in
    # NONMEM (lower, init, upper) syntax and is up to 5.1% away from Table 2
    # (e.g. exp(7.76) = 2345 L vs Table 2 V2/F1 = 2230.5 L), so it is not
    # used for fixed effects. See the vignette Errata for the full comparison.
    #
    # Everything is an apparent parameter: CL, V1, V2, and Q are divided by
    # F1, and the metabolite CL and V are divided by (Fm x F1). F1 and Fm are
    # both fixed, so the Table 2 numbers are used directly.
    # ---------------------------------------------------------------------

    # ---- Parent imatinib disposition (Said 2025 Table 2) ------------------
    lka <- log(0.17)    ; label("First-order oral absorption rate constant ka (1/h)")                     # Said 2025 Table 2 ka = 0.17 (RSE 3.5%)
    lcl <- log(298.9)   ; label("Apparent clearance of unbound imatinib CLu/F1 (L/h)")                    # Said 2025 Table 2 CLu/F1 = 298.9 (RSE 0.6%)
    lvc <- log(25.5)    ; label("Apparent central volume of distribution V1/F1 (L)")                      # Said 2025 Table 2 V1/F1 = 25.5 (RSE 3.3%)
    lvp <- log(2230.5)  ; label("Apparent peripheral volume of distribution V2/F1 (L)")                   # Said 2025 Table 2 V2/F1 = 2230.5 (RSE 0.5%)
    lq  <- log(626.4)   ; label("Apparent inter-compartmental clearance Q/F1 (L/h)")                      # Said 2025 Table 2 Q/F1 = 626.4 (RSE 1.4%)
    lkd <- log(368.7)   ; label("Imatinib-AAG equilibrium dissociation constant Kd (ug/L)")               # Said 2025 Table 2 Kd = 368.7 (RSE 0.7%)

    # F1 was fixed to the 98% absolute bioavailability determined in a prior
    # dedicated bioavailability study (Said 2025 Methods 2.3.1, reference 44).
    lfdepot <- fixed(log(0.98)) ; label("Oral bioavailability F1 (unitless fraction)")                    # Said 2025 Table 2 F1 = 0.98 (FIX)

    # ---- N-desmethyl imatinib (Said 2025 Table 2) ------------------------
    # Fm was fixed to 15% from literature data (Said 2025 Methods 2.3.1,
    # reference 47) to avoid an identifiability problem: only imatinib was
    # administered, and when equal parent and metabolite distribution volumes
    # were assumed Fm was estimated beyond 1. The model was therefore
    # reparameterised to estimate CLm/(Fm x F1) and Vm/(Fm x F1) instead.
    fm <- fixed(0.15)   ; label("Fraction of unbound imatinib clearance forming DM-imatinib Fm (unitless)")  # Said 2025 Table 2 Fm = 0.15 (FIX)
    lcl_ndmima <- log(190.6) ; label("Apparent clearance of unbound DM-imatinib CLm/(Fm x F1) (L/h)")     # Said 2025 Table 2 CLm/(Fm x F1) = 190.6 (RSE 1.4%)
    lvc_ndmima <- log(11.0)  ; label("Apparent volume of distribution of DM-imatinib Vm/(Fm x F1) (L)")   # Said 2025 Table 2 Vm/(Fm x F1) = 11.0 (RSE 5.8%)

    # ---- Retained covariate effects (Said 2025 Table 2) ------------------
    # The Table 2 footnote states that covariate coefficients are reported
    # already exponentially transformed, as exp(COVcat x theta_parameter).
    # They are therefore log-additive shifts here, written as log(<Table 2
    # value>) so each ini() value reads directly against the printed table.
    e_dis_covid19_vc         <- log(1.2)  ; label("Effect of COVID-19 (vs CML/GIST) on apparent central volume V1/F1 (log-scale shift)")            # Said 2025 Table 2 V1-COVID-19 = 1.2 (RSE 60%)
    e_conmed_il6ri_kd        <- log(1.7)  ; label("Effect of IL-6R-inhibitor cotreatment on the imatinib-AAG Kd (log-scale shift)")                 # Said 2025 Table 2 Kd-IL6R inhibitor = 1.7 (RSE 20.3%)
    e_conmed_il6ri_cl_ndmima <- log(0.46) ; label("Effect of IL-6R-inhibitor cotreatment on apparent DM-imatinib clearance (log-scale shift)")      # Said 2025 Table 2 CLm-IL6R inhibitor = 0.46 (RSE 13.1%)

    # ---- Inter-individual variability (Said 2025 Table 2) ----------------
    # IIV was modelled exponentially (Said 2025 Methods 2.3.1) and Table 2
    # reports it as %CV, so the log-scale variance is omega^2 = log(1 + CV^2):
    #   ka                77.3% -> log(1 + 0.773^2) = 0.468458
    #   CLu/F1            32.9% -> log(1 + 0.329^2) = 0.102774
    #   V1/F1             88.3% -> log(1 + 0.883^2) = 0.576439
    #   Kd                25.8% -> log(1 + 0.258^2) = 0.064442
    #   CLm/(Fm x F1)     41.7% -> log(1 + 0.417^2) = 0.160322
    #   Vm/(Fm x F1)      12.8% -> log(1 + 0.128^2) = 0.016251
    # The Data S1 $OMEGA block independently reproduces the first four of
    # these (0.103, 0.576, 0.0645, 0.468) but not the two metabolite terms;
    # Table 2 is used throughout. See vignette Errata.
    # Said 2025 reports no IIV on V2/F1 or Q/F1 and none is added here.
    etalka        ~ 0.468458   # Said 2025 Table 2 IIV ka = 77.3% (RSE 6.0%, shrinkage 22%)
    etalcl        ~ 0.102774   # Said 2025 Table 2 IIV CLu/F1 = 32.9% (RSE 8.0%, shrinkage 38%)
    etalvc        ~ 0.576439   # Said 2025 Table 2 IIV V1/F1 = 88.3% (RSE 20%, shrinkage 61%)
    etalkd        ~ 0.064442   # Said 2025 Table 2 IIV Kd = 25.8% (RSE 22%, shrinkage 68%)
    etalcl_ndmima ~ 0.160322   # Said 2025 Table 2 IIV CLm/(Fm x F1) = 41.7% (RSE 35%, shrinkage 13%)
    etalvc_ndmima ~ 0.016251   # Said 2025 Table 2 IIV Vm/(Fm x F1) = 12.8% (RSE 37%, shrinkage 72%)

    # ---- Residual error --------------------------------------------------
    # The model is log-transform-both-sides: Data S1 $ERROR sets
    # IPRED = LOG(TY) and Y = IPRED + ERR(n), so each $SIGMA is a variance on
    # the natural-log scale and maps to nlmixr2's lnorm() residual with
    # expSd = sqrt(sigma^2).
    #
    # The endpoint-to-variance assignment follows the Data S1 control stream,
    # per operator ruling on task sidecar oare_PMC11919263 request-001 q2.
    # Data S1 maps ERR(1) -> Q0 -> CTOT (total), ERR(2) -> Q1 -> CFREE
    # (unbound), ERR(3) -> Q2 -> CMTOT (metabolite) positionally, and its
    # inline comments (;CT, ;CF, ;CM) agree. Said 2025 Table 2 prints the
    # total and unbound rows the other way round. The two readings, and why
    # the control stream was preferred, are set out in the vignette Errata.
    #   total imatinib      sigma^2 = 0.103  -> sqrt = 0.320936
    #   unbound imatinib    sigma^2 = 0.176  -> sqrt = 0.419524
    #   total DM-imatinib   sigma^2 = 0.0501 -> sqrt = 0.223830
    expSd        <- 0.320936 ; label("Exponential (log-scale) residual SD for total imatinib Cc")            # Said 2025 Data S1 $SIGMA 0.103 ;CT
    expSd_Cu     <- 0.419524 ; label("Exponential (log-scale) residual SD for unbound imatinib Cu")          # Said 2025 Data S1 $SIGMA 0.176 ;CF
    expSd_ndmima <- 0.223830 ; label("Exponential (log-scale) residual SD for total DM-imatinib Cc_ndmima")  # Said 2025 Data S1 $SIGMA 0.0501 ;CM
  })

  model({
    # ---- Individual parameters -------------------------------------------
    # Covariate effects enter as log-additive shifts inside exp(), matching
    # Said 2025 Equation 3 and the Data S1 $PK block
    # (e.g. V = EXP(THETA(2) + THETA(8)*COVID + ETA(2))).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + e_dis_covid19_vc * DIS_COVID19 + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    kd <- exp(lkd + e_conmed_il6ri_kd * CONMED_IL6RI + etalkd)

    cl_ndmima <- exp(lcl_ndmima + e_conmed_il6ri_cl_ndmima * CONMED_IL6RI + etalcl_ndmima)
    vc_ndmima <- exp(lvc_ndmima + etalvc_ndmima)

    # ---- Micro-constants (Data S1 $PK) -----------------------------------
    kel        <- cl / vc               # K20
    k12        <- q / vc                # K23, central -> peripheral1
    k21        <- q / vp                # K32, peripheral1 -> central
    kel_ndmima <- cl_ndmima / vc_ndmima # K40

    # ---- Saturable AAG binding of imatinib (Said 2025 Equations 1-2) -----
    # Amounts are carried in mg and volumes in L, so central/vc is mg/L;
    # x 1000 gives the ug/L working unit of the binding equations (the Data S1
    # scaling S2 = V/1000 does the same thing).
    Cc <- 1000 * central / vc

    # bmax is the imatinib-equivalent concentration of AAG binding sites,
    # i.e. Bmax in the Langmuir form Cb = bmax * Cu / (Kd + Cu) of Equation 1.
    # The scale factor L converts AAG in g/L to ug/L of imatinib-equivalent
    # sites under a 1:1 molar binding ratio and is pure unit arithmetic on the
    # molecular weights, not a fitted parameter: 493.6 / 42000 x 1e6 = 11752
    # for imatinib and AAG respectively. Said 2025 Methods 2.3.1 fixed it to
    # 11700 (and Data S1 $PK to LMET = 11400 for DM-imatinib, from
    # 479.6 / 42000 x 1e6 = 11419), because L and Kd are correlated and cannot
    # be estimated independently. The rounded published constants are used.
    bmax <- 11700 * AAG

    # Closed-form positive root of the instantaneous-binding equilibrium
    # (Said 2025 Equation 2 / Data S1 $DES C2).
    Cu <- 0.5 * ((Cc - bmax - kd) + sqrt((Cc - bmax - kd)^2 + 4 * kd * Cc))

    # The 1e-6 guard reproduces Data S1 $DES `FU = C2/(C2T + 1E-6)` and keeps
    # fu finite at Cc = 0, where the numerator is also 0.
    fu <- Cu / (Cc + 1e-6)

    # ---- Saturable AAG binding of DM-imatinib ---------------------------
    # The metabolite binds AAG with the SAME Kd as the parent, including the
    # IL-6R-inhibitor effect on Kd: unbound DM-imatinib could not be measured,
    # so its binding affinity was fixed to the parent's (Said 2025 Discussion
    # limitations; Data S1 $DES uses KD in the CM2 expression).
    Cc_ndmima <- 1000 * central_ndmima / vc_ndmima
    bmax_ndmima <- 11400 * AAG
    Cu_ndmima <- 0.5 * ((Cc_ndmima - bmax_ndmima - kd) +
                        sqrt((Cc_ndmima - bmax_ndmima - kd)^2 + 4 * kd * Cc_ndmima))
    fu_ndmima <- Cu_ndmima / (Cc_ndmima + 1e-6)

    # ---- ODE system (Data S1 $DES) --------------------------------------
    # The unbound fraction drives elimination, metabolite formation, AND the
    # central -> peripheral distribution flux, because the bound drug acts as
    # a reservoir that continuously re-equilibrates with unbound drug
    # (Said 2025 Methods 2.3.1). The return flux from the peripheral
    # compartment is NOT scaled by fu: Data S1 $DES has `- K23*A(2)*FU` out of
    # central but plain `+ K32*A(3)` back in. That asymmetry is reproduced
    # verbatim here.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * fu * central -
                          k12 * fu * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * fu * central - k21 * peripheral1

    # Metabolite formation takes the Fm share of the unbound-driven parent
    # elimination flux, converted for the molecular-weight difference between
    # imatinib (493.6 g/mol) and DM-imatinib (479.6 g/mol). Data S1 $PK sets
    # CF = 0.97, the rounded ratio 479.6 / 493.6 = 0.9716, and applies it to
    # the formation flux only. Metabolite elimination is likewise driven by
    # the metabolite's own unbound fraction.
    d/dt(central_ndmima) <- 0.97 * fm * fu * kel * central -
                             fu_ndmima * kel_ndmima * central_ndmima

    # F1 applies to the oral depot only; IV doses go straight into central
    # (Data S1 $PK F1 = 0.98 on compartment 1 = ORALT).
    f(depot) <- exp(lfdepot)

    # ---- Observations ----------------------------------------------------
    # All three analytes are in ug/L. Log-transform-both-sides in NONMEM with
    # an additive error on the log scale is an exponential (log-normal)
    # residual in nlmixr2.
    Cc        ~ lnorm(expSd)
    Cu        ~ lnorm(expSd_Cu)
    Cc_ndmima ~ lnorm(expSd_ndmima)
  })
}
