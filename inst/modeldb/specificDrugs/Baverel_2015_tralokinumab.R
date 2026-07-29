Baverel_2015_tralokinumab <- function() {
  description <- "Two-compartment population PK model for tralokinumab in adolescent (12-17 y) and adult subjects with asthma or healthy volunteers (Baverel 2015), with parallel subcutaneous absorption (first-order with lag plus zero-order over a fixed duration), allometric body-weight scaling on disposition parameters, and an additional 15% lower clearance in adolescents."
  reference <- paste(
    "Baverel PG, Jain M, Stelmach I, She D, Agoram B, Sandbach S, Piper E, Kuna P.",
    "Pharmacokinetics of tralokinumab in adolescents with asthma: implications for future dosing.",
    "Br J Clin Pharmacol. 2015;80(6):1337-1349. doi:10.1111/bcp.12725."
  )
  vignette <- "Baverel_2015_tralokinumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline; allometric power scaling on disposition parameters with reference weight 73 kg (pooled-cohort median; Baverel 2015 Equation 5). Exponents fixed at 0.75 for CL and Q and at 1 for Vc and Vp.",
      source_name        = "WT"
    ),
    ADOLESCENT = list(
      description        = "Adolescent age-cohort indicator (1 = age 12-17 years; 0 = adult >= 18 years)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Categorical effect retained on clearance only; corresponds to Baverel 2015 Equation 7 with adolescents having a 15% lower CL than adults (Table 3 row 'CL decrease, % (adolescent)' = 15%).",
      source_name        = "ADOLESCENT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 578L,
    n_studies       = 8L,
    n_observations  = 5504L,
    age_range       = "12-75 years (adolescents 12-17; adults 18-75)",
    weight_range    = "36-115 kg overall; adolescents 40-94 kg",
    weight_median   = "73 kg overall (adolescents 60 kg)",
    sex_female_pct  = 53.5,
    race_ethnicity  = c(White = 69, Asian = 23, Black = 2, Other = 5),
    disease_state   = "Adolescents and adults with asthma plus healthy volunteers; pooled across 5 phase I studies and 3 phase II studies",
    dose_range      = "Single SC 300 mg in the adolescent study; pooled IV 0.1-30 mg/kg single dose and SC <= 600 mg multiple-dose Q2W or Q4W in adults",
    regions         = "Multinational; the adolescent phase I study was conducted in Poland (NCT01592396)",
    n_adolescents   = 20L,
    n_pk_samples_adolescent = 202L,
    notes           = "Demographics pooled across 8 studies (Baverel 2015 Table 1). Adolescent contribution is 20 subjects (3.5%) and 202 PK samples (3.7%) from NCT01592396; the remaining 558 subjects came from prior phase I (NCT01093040, NCT00638989, CAT-354-401, NCT00974675) and phase II (NCT00640016, NCT00873860, NCT01402986) studies in adults. Sex distribution computed from Table 1 row 'Male, n (%)': 269/578 male = 46.5% male, so 53.5% female. Japanese subtype contributed 64/578 (11.1%) but did not enter the final PK model as a covariate."
  )

  ini({
    # Structural population parameters from Baverel 2015 Table 3 (reference adult,
    # WT = 73 kg). Paper reports CL and Q in mL/day and Vc/Vp in mL; converted to
    # L/day and L by dividing by 1000 so the model writes out concentrations in
    # mg/L = ug/mL when doses are in mg.
    lcl  <- log(0.204);  label("Clearance for a 73 kg adult (CL, L/day)")                           # Table 3: CL adult 204 mL/day
    lvc  <- log(1.867);  label("Central volume of distribution for a 73 kg adult (Vc, L)")          # Table 3: Vc 1867 mL
    lvp  <- log(3.357);  label("Peripheral volume of distribution for a 73 kg adult (Vp, L)")       # Table 3: Vp 3357 mL
    lq   <- log(1.582);  label("Intercompartmental clearance for a 73 kg adult (Q, L/day)")         # Table 3: Q 1582 mL/day

    # SC absorption parameters. Absorption is split between a first-order pathway
    # (fraction Fr; rate Ka with absorption lag Tlag) and a zero-order pathway
    # (fraction 1-Fr; duration D0). Fsc is the SC bioavailability that scales the
    # entire SC dose; IV dosing implicitly carries F = 1. Fr is logit-transformed
    # in Baverel 2015 Equations 2 and 3 to constrain it in (0, 1).
    lka      <- log(0.34);             label("First-order absorption rate constant (Ka, 1/day)")   # Table 3: Ka 0.34 /day
    ld0      <- log(5.7);              label("Zero-order absorption duration (D0, day)")           # Table 3: D0 5.7 days
    ltlag    <- log(0.8);              label("First-order absorption lag time (Tlag, day)")        # Table 3: Tlag 0.8 day
    lfdepot  <- log(0.8);              label("Subcutaneous bioavailability (Fsc, fraction)")       # Table 3: Fsc 0.8
    logitfr  <- log(0.7 / (1 - 0.7));  label("Logit of fraction of SC dose absorbed first-order (Fr, unitless)")  # Table 3: Fr 0.7 (logit-transformed per Eqs. 2-3)

    # Allometric exponents on disposition parameters (Baverel 2015 Equation 5,
    # reference weight 73 kg). Fixed to canonical mAb values per the paper's
    # Methods and Results: "fixed to prior knowledge (0.75 for CL and Q, and 1
    # for Vc and Vp) ... the choice of final model was to keep fixed exponents".
    e_wt_cl_q   <- fixed(0.75);  label("Allometric exponent on CL and Q (unitless, fixed)")        # Table 3 row 'Effect of body weight' fixed per Methods
    e_wt_vc_vp  <- fixed(1.00);  label("Allometric exponent on Vc and Vp (unitless, fixed)")       # Table 3 row 'Effect of body weight' fixed per Methods

    # Categorical adolescent effect on CL. Baverel 2015 Equation 7 parameterises
    # this as an additive shift, with Table 3 reporting "CL decrease, %
    # (adolescent) = 15%". Encoded multiplicatively here (CL_adolescent =
    # CL_adult * (1 + e_adolescent_cl * ADOLESCENT)) which matches the table:
    # 204 * (1 - 0.15) = 173.4 mL/day, vs Table 3 adolescent CL of 173 mL/day.
    e_adolescent_cl <- -0.15;  label("Fractional decrease in CL for adolescents vs adults (unitless)")  # Table 3 'CL decrease, % (adolescent)' 15%

    # IIV. Baverel 2015 Table 3 reports bootstrap 95% CI for the between-subject
    # CV%; point estimates of the CV are not tabulated separately, so the
    # bootstrap-CI midpoints are used. The narrative on page 1340 confirms the
    # approximate magnitudes ("CV for CL and Vp was ~30%, although this was
    # slightly higher for Vc and Q (~65%)"). Log-normal IIV variance is
    # omega^2 = log(1 + CV^2):
    #   CL: CV = (27.0 + 31.2) / 2 = 29.1% -> omega^2 = log(1 + 0.291^2) = 0.08133
    #   Vc: CV = (59.4 + 70.0) / 2 = 64.7% -> omega^2 = log(1 + 0.647^2) = 0.34956
    #   Vp: CV = (20.2 + 33.3) / 2 = 26.75% -> omega^2 = log(1 + 0.2675^2) = 0.06915
    #   Q:  CV = (54.3 + 99.3) / 2 = 76.8% -> omega^2 = log(1 + 0.768^2) = 0.46362
    #
    # Block correlations are taken directly from Table 3
    # ("Correlations of interindividual variability estimates"; point estimates):
    #   Corr(CL, Vc) =  0.5  -> cov =  0.5  * sqrt(0.08133 * 0.34956) = 0.08431
    #   Corr(CL, Vp) =  0.2  -> cov =  0.2  * sqrt(0.08133 * 0.06915) = 0.01499
    #   Corr(CL, Q)  = -0.3  -> cov = -0.3  * sqrt(0.08133 * 0.46362) = -0.05832
    #   Corr(Vc, Vp) = -0.3  -> cov = -0.3  * sqrt(0.34956 * 0.06915) = -0.04663
    #   Corr(Vc, Q)  = -0.6  -> cov = -0.6  * sqrt(0.34956 * 0.46362) = -0.24161
    #   Corr(Vp, Q)  =  0.5  -> cov =  0.5  * sqrt(0.06915 * 0.46362) =  0.08964
    etalcl + etalvc + etalvp + etalq ~ c(
       0.08133,
       0.08431,  0.34956,
       0.01499, -0.04663,  0.06915,
      -0.05832, -0.24161,  0.08964,  0.46362
    )                                                                                              # Table 3 IIV block on {CL, Vc, Vp, Q}

    # Fr IIV on the logit-transformed scale (Equation 2-3, Table 3 row 'Fr'):
    # bootstrap CV% midpoint (8.6 + 32.8) / 2 = 20.7%, interpreted as the SD of
    # eta_logitfr on the logit scale -> omega^2 = 0.207^2 = 0.04285.
    etalogitfr ~ 0.04285                                                                           # Table 3 IIV row 'Fr', 8.6-32.8% on logit scale

    # Combined additive + proportional residual error (Baverel 2015 Equation 4,
    # Table 3 rows 'Additive residual error' and 'Proportional residual error').
    addSd  <- 0.2;    label("Additive residual error (ug/mL)")                                     # Table 3: epsilon_2 = 0.2 ug/mL
    propSd <- 0.198;  label("Proportional residual error (fraction)")                              # Table 3: epsilon_1 = 19.8 %
  })

  model({
    # Individual disposition parameters with allometric body-weight scaling
    # (Equation 5, reference weight 73 kg) and the adolescent categorical effect
    # on CL (Equation 7, multiplicatively reparameterised as documented in ini()).
    cl <- exp(lcl + etalcl) * (WT / 73)^e_wt_cl_q  * (1 + e_adolescent_cl * ADOLESCENT)
    vc <- exp(lvc + etalvc) * (WT / 73)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 73)^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * (WT / 73)^e_wt_cl_q

    # Absorption parameters. The logit transform of Fr keeps individual values in
    # (0, 1); 1-fr is the zero-order share. fdepot is the SC bioavailability that
    # scales the entire SC dose.
    ka     <- exp(lka)
    d0     <- exp(ld0)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot)
    fr     <- exp(logitfr + etalogitfr) / (1 + exp(logitfr + etalogitfr))

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Parallel SC absorption (Baverel 2015 Figure 2). The user supplies the SC
    # dose as two dose records at the same time and same amt: one to depot
    # (first-order with lag Tlag and rate Ka, bioavailability Fsc*Fr) and one
    # to depot2 (zero-order infusion over duration D0, bioavailability
    # Fsc*(1-Fr); use rate = -2 in the event record to invoke the model-defined
    # dur). depot2 holds the zero-order infusate transiently and drains into
    # central at a fast rate (100/day; t1/2 ~ 0.007 day) that is two orders of
    # magnitude faster than D0 so the input to central is effectively the
    # zero-order rate. Keeping the zero-order route on depot2 (rather than on
    # central directly) leaves f(central) and dur(central) free for IV
    # bolus / infusion dosing (the paper's pooled dataset includes adult IV
    # studies at 0.1-30 mg/kg, NCT00638989).
    kdepot2 <- 100  # 1/day; fast drainage of the SC zero-order transit pool into central
    d/dt(depot)       <- -ka * depot
    d/dt(depot2)      <- -kdepot2 * depot2
    d/dt(central)     <-  ka * depot + kdepot2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                                    k12 * central - k21 * peripheral1

    # First-order SC route: bioavailable fraction Fsc * Fr, with absorption lag.
    f(depot)    <- fdepot * fr
    alag(depot) <- tlag

    # Zero-order SC route: bioavailable fraction Fsc * (1-Fr), infused into
    # depot2 over duration D0.
    f(depot2)   <- fdepot * (1 - fr)
    dur(depot2) <- d0

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL (matches Baverel
    # 2015 reporting in ug/mL on the y-axis of Figures 3 and 5).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
