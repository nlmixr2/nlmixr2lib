Barnett_2018_rosuvastatin <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption for a single 5 mg dose of rosuvastatin in healthy adult males (Barnett 2018), refit from the Tzeng 2008 structural form with simultaneous plasma + urine fitting. The model includes separable biliary (CLb,RSV) and renal (CLr,RSV) clearance components from the central compartment, competitive rifampicin OATP1B inhibition of the biliary clearance via KiRSV driven by the instantaneous plasma rifampicin concentration, and a binary RIF-coadministration covariate that captures paper-reported reductions of V1, V2, and Q during the rifampicin phase (Barnett 2018 Table 1: V1 430 -> 2.98 L, V2 865 -> 128 L, Q 45.3 -> 5.03 L/h on RIF). Companion to modellib('Barnett_2018_coproporphyrin_I'); both share the rifampicin perpetrator parameterisation in modellib('Barnett_2018_rifampicin')."
  reference <- paste(
    "Barnett S, Ogungbenro K, Menochet K, Shen H, Lai Y, Humphreys WG, Galetin A.",
    "Gaining Mechanistic Insight Into Coproporphyrin I as Endogenous Biomarker",
    "for OATP1B-Mediated Drug-Drug Interactions Using Population Pharmacokinetic",
    "Modeling and Simulation.",
    "Clin Pharmacol Ther. 2018;104(3):564-574.",
    "doi:10.1002/cpt.983.",
    "Structural rosuvastatin popPK model adapted from",
    "Tzeng TB et al. Curr Med Res Opin. 2008;24(9):2575-2585.",
    "doi:10.1185/03007990802312807.",
    sep = " "
  )
  vignette <- "Barnett_2018_rosuvastatin"
  units <- list(time = "hour", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion / period indicator (Barnett 2018 study design: OCC1 = rifampicin-only period, OCC2 = rosuvastatin-only period, OCC3 = combined rifampicin + rosuvastatin period).",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Time-varying within subject; constant within an occasion. Rosuvastatin was dosed on OCC2 and OCC3 in the source clinical study (Lai et al. 2016 n=12 healthy-male SLCO1B1-wildtype cohort). Decomposed inside model() into binary indicators oc1, oc2, oc3 that multiplex the per-occasion IOV eta on log-Ka (only Ka carried IOV in the Barnett 2018 RSV fit, with a single shared variance across the two RSV occasions). OCC1 (RIF-only period) has no rosuvastatin data in the source fit; its eta is unused for in-paper simulations and is retained only for users who want to simulate at all three study occasions.",
      source_name        = "OCC"
    ),
    CONMED_RIF = list(
      description        = "Concomitant single-dose rifampicin co-administration indicator (1 = within a 600 mg rifampicin co-administration period in the Barnett 2018 study design; 0 = rosuvastatin alone).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (rosuvastatin alone)",
      notes              = "Time-varying within subject. Used here as the binary period-level covariate that captures Barnett 2018 Table 1's reductions of V1, V2, and Q during the rifampicin phase: V1_RSV 430 -> 2.98 L (~99% reduction), V2_RSV 865 -> 128 L (~85% reduction), Q_RSV 45.3 -> 5.03 L/h (~89% reduction). Distinct semantics from the multi-day CYP3A4-induction use of CONMED_RIF in Svensson_2014_bedaquiline: here rifampicin is acute (single 600 mg oral dose) and acts predominantly as a competitive OATP1B inhibitor, with the empirical distribution-parameter changes captured by this covariate likely reflecting OATP1B-inhibition-driven changes in hepatic distribution rather than enzyme induction.",
      source_name        = "CONMED_RIF"
    ),
    CP_RIF_UM = list(
      description        = "Instantaneous rifampicin plasma concentration as a time-varying perpetrator covariate driving competitive OATP1B inhibition of biliary rosuvastatin clearance (Barnett 2018 Methods, RSV model section).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the rifampicin co-administration window so the OATP1B inhibition term in the central-compartment ODE collapses to the baseline form. The PK trajectory of rifampicin is parameterised in modellib('Barnett_2018_rifampicin'); users typically simulate the rifampicin model first and then feed its central-compartment concentration (after MW conversion to umol/L; rifampicin MW 822.94 g/mol) as the CP_RIF_UM column on the RSV event table. Reference peak: ~29 umol/L (see modellib('Barnett_2018_rifampicin') typical-value simulation).",
      source_name        = "CRIF"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 12L,
    n_studies        = 1L,
    n_observations   = 264L,
    age_range        = "Healthy adult males; demographic distribution not tabulated by Barnett 2018 (source dataset Lai et al. 2016, 12 healthy male subjects, SLCO1B1 c.521 T>C wildtype only; no OATP1B1*5 / *15 carriers).",
    weight_range     = "(not extracted; Barnett 2018 Methods do not tabulate per-subject weights for the n=12 cohort.)",
    sex_female_pct   = 0,
    disease_state    = "Healthy adult male volunteers in a three-occasion (7-day washout between OCC1-2 and OCC2-3) drug-drug-interaction crossover study.",
    dose_range       = "Single 5 mg oral rosuvastatin dose at the start of OCC2 and at the start of OCC3 (co-administered with 600 mg rifampicin on OCC3).",
    regions          = "(not extracted; the underlying clinical study Lai et al. 2016 region was not explicitly stated in Barnett 2018 Methods.)",
    notes            = "Demographics inferred from Barnett 2018 Methods (Clinical data section) and the cited source clinical study Lai Y et al., Pharmacol Res Perspect 2016;4(3):e00207. The 264 RSV plasma samples (OCC2 = 132, OCC3 = 132) plus 44 RSV urine samples (cumulative excretion over 0-7 h and 7-24 h post-dose) were fit simultaneously in NONMEM using FOCE."
  )

  ini({
    # Structural parameters -- Barnett 2018 Table 1, rosuvastatin rows
    # (column 'Estimate (SE%) Population'). SE% in parentheses next to
    # each value. V1, V2, and Q carry a binary RIF-coadministration
    # covariate effect (separate rows in Table 1, suffixed '(RIF) b');
    # the multiplicative factor encoding is done in model(), not in ini().
    lka     <- log(0.287); label("Absorption rate constant Ka (1/h)")                                        # Table 1, RSV row 'ka RSV (1/h)' = 0.287 (SE 9%)
    lclr    <- log(8.48);  label("Apparent renal clearance CLr,RSV (L/h)")                                   # Table 1, RSV row 'CL R,RSV (L/h)' = 8.48 (SE 14%)
    lclb    <- log(124);   label("Apparent biliary clearance CLb,RSV (L/h, baseline phase)")                 # Table 1, RSV row 'CL b,RSV (L/h)' = 124 (SE 11%)
    lvc     <- log(430);   label("Apparent central volume of distribution V1,RSV (L, baseline phase)")       # Table 1, RSV row 'V1 RSV (L)' = 430 (SE 12%)
    lq      <- log(45.3);  label("Apparent intercompartmental clearance Q,RSV (L/h, baseline phase)")        # Table 1, RSV row 'Q RSV (L/h)' = 45.3 (SE 19%)
    lvp     <- log(865);   label("Apparent peripheral volume of distribution V2,RSV (L, baseline phase)")    # Table 1, RSV row 'V2 RSV (L)' = 865 (SE 37%)
    lki     <- log(2.23);  label("Total rifampicin OATP1B inhibition constant KiRSV (umol/L)")               # Table 1, RSV row 'Ki RSV (uM) c' = 2.23 (SE 15%); footnote c notes the unbound Ki = 0.25 umol/L after correction by RIF fu = 0.11.

    # Binary RIF-coadministration covariate effects on V1, Q, V2. The
    # paper reports separate '(RIF)' rows in Table 1; the encoding here
    # uses the multiplicative form param = param_baseline * (1 + e * CONMED_RIF)
    # with e = param_RIF / param_baseline - 1.
    e_rif_vc <- -0.99307;  label("Multiplicative reduction in V1,RSV during the rifampicin phase (unitless)")  # Table 1, RSV row 'V1 RSV (L) (RIF) b' = 2.98 (SE 50%); 2.98/430 - 1 = -0.99307
    e_rif_q  <- -0.88896;  label("Multiplicative reduction in Q,RSV during the rifampicin phase (unitless)")   # Table 1, RSV row 'Q RSV (L/h) (RIF) b' = 5.03 (SE 18%); 5.03/45.3 - 1 = -0.88896
    e_rif_vp <- -0.85202;  label("Multiplicative reduction in V2,RSV during the rifampicin phase (unitless)")  # Table 1, RSV row 'V2 RSV (L) (RIF) b' = 128 (SE 42%); 128/865 - 1 = -0.85202

    # IIV -- log-normal variances computed from Table 1 IIV CV% column
    # via the standard omega^2 = log(1 + CV^2) conversion. V2,RSV had
    # no IIV estimated in the source (Table 1 shows '-').
    etalka  ~ 0.02344  # Table 1, RSV IIV ka  = 15.4% CV (RSE 57%);   log(1 + 0.154^2) = 0.02344
    etalclr ~ 0.13427  # Table 1, RSV IIV CLr = 37.9% CV (RSE 28.8%); log(1 + 0.379^2) = 0.13427
    etalclb ~ 0.11187  # Table 1, RSV IIV CLb = 34.4% CV (RSE 22%);   log(1 + 0.344^2) = 0.11187
    etalvc  ~ 0.04643  # Table 1, RSV IIV V1  = 21.8% CV (RSE 48.4%); log(1 + 0.218^2) = 0.04643
    etalq   ~ 0.18581  # Table 1, RSV IIV Q   = 45.2% CV (RSE 34.4%); log(1 + 0.452^2) = 0.18581
    etalki  ~ 0.10577  # Table 1, RSV IIV Ki  = 33.4% CV (RSE 34.5%); log(1 + 0.334^2) = 0.10577

    # IOV on log-Ka across the two RSV occasions (OCC2 and OCC3 in the
    # source data; we keep three eta slots for OCC1..OCC3 to match the
    # study's overall occasion structure, with the same shared variance).
    etaiov_ka_1 ~ 0.03085          # Table 1, RSV IOV ka = 17.7% CV (RSE 42.7%); log(1 + 0.177^2) = 0.03085
    etaiov_ka_2 ~ fix(0.03085)
    etaiov_ka_3 ~ fix(0.03085)

    # Residual error -- combined additive + proportional for both
    # plasma and urine, with the additive components fixed (per Table 1
    # 'FIXED' annotation on the additive rows). Note the urine
    # additive row in Table 1 is labelled 'r prop (nM) - urine 0.1 FIXED'
    # which is interpreted here as a typo for 'r add (nM) - urine 0.1
    # FIXED' (a proportional residual with units of concentration
    # rather than fraction is not a well-defined NONMEM error structure;
    # the surrounding 'r prop (%) - urine 0.578' row already supplies
    # the proportional component). The proportional values are taken
    # verbatim from Table 1: 0.257 for plasma and 0.578 for urine.
    # These values are physically plausible only when read as fractions
    # (25.7% CV plasma, 57.8% CV urine), so they are encoded as the
    # propSd value directly rather than divided by 100; see the
    # validation vignette's 'Assumptions and deviations' section for
    # the unit-interpretation rationale (the column header '(%)' in
    # Table 1 conflicts with the magnitude under a strict percent
    # reading, especially for CPI and RIF where the same column shows
    # 13.9 / 34.2 / 31.3 as percents -- a within-table inconsistency
    # that the model encodes by reading each row in its physically
    # plausible scale).
    propSd      <- 0.257;        label("Proportional residual error, RSV plasma (fraction)")  # Table 1, RSV row 'r prop (%) - plasma' = 0.257 (SE 6%), read as fraction.
    addSd       <- fixed(0.05);  label("Additive residual error, RSV plasma (nmol/L, fixed)")  # Table 1, RSV row 'r add (nM) - plasma' = 0.05 FIXED
    propSd_Ursv <- 0.578;        label("Proportional residual error, RSV urine (fraction)")    # Table 1, RSV row 'r prop (%) - urine'  = 0.578 (SE 16%), read as fraction.
    addSd_Ursv  <- fixed(0.1);   label("Additive residual error, RSV urine (nmol/L, fixed)")   # Table 1, RSV row 'r prop (nM) - urine 0.1 FIXED' -- interpreted as the additive urine component (likely typo for 'r add (nM) - urine'); see in-file comment block above for rationale.
  })

  model({
    # Decompose the integer-valued OCC column into binary indicators
    # for IOV multiplexing on log-Ka.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3

    # Individual typical-value parameters with log-normal IIV / IOV.
    # The distribution parameters V1, V2, Q carry the binary RIF
    # covariate effect; CL components (CLr, CLb) and Ka, Ki do not
    # (paper only reports the RIF effect on the distribution-parameter
    # rows in Table 1).
    ka     <- exp(lka  + etalka + iov_ka)
    cl_r   <- exp(lclr + etalclr)
    cl_b   <- exp(lclb + etalclb)
    vc     <- exp(lvc  + etalvc)  * (1 + e_rif_vc * CONMED_RIF)
    q      <- exp(lq   + etalq)   * (1 + e_rif_q  * CONMED_RIF)
    vp     <- exp(lvp)            * (1 + e_rif_vp * CONMED_RIF)
    ki     <- exp(lki  + etalki)

    # Rifampicin competitive inhibition of biliary RSV clearance
    # (Barnett 2018 Methods, RSV model section: analogous Eq. to the
    # CPI form). cl_b_eff collapses to cl_b when CP_RIF_UM = 0.
    cl_b_eff <- cl_b / (1 + CP_RIF_UM / ki)

    # Micro-constants for the two-compartment disposition.
    cl_total <- cl_r + cl_b_eff
    kel <- cl_total / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system. depot -> central via first-order ka. central <->
    # peripheral1. urine accumulates the renal-cleared portion of the
    # central state (cl_r * Cc), in nmol.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                              k12 * central - k21 * peripheral1
    d/dt(urine)       <- cl_r * (central / vc)

    # Plasma concentration. Dose in mg, V1 in L -> central/vc has units
    # of mg/L; convert to nmol/L by multiplying by 1e6 / MW_rosuvastatin
    # (rosuvastatin MW = 481.54 g/mol).
    MW_rosuvastatin <- 481.54
    Cc   <- (central / vc) * 1e6 / MW_rosuvastatin
    Ursv <- urine

    Cc   ~ add(addSd)      + prop(propSd)
    Ursv ~ add(addSd_Ursv) + prop(propSd_Ursv)
  })
}
