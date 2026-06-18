PerezRuixo_2020_esketamine <- function() {
  description <- "Joint three-compartment esketamine + two-compartment apparent noresketamine population PK model with a hepato-portal first-pass compartment (well-stirred model) and three parallel absorption routes (intranasal direct, intranasal-swallowed via PO depot, and PO solution via PO depot) developed from 9784/9397 esketamine/noresketamine plasma observations in 820 healthy volunteers and patients with treatment-resistant depression receiving intranasal, intravenous, and oral esketamine (Perez-Ruixo 2020). Asian race decreases esketamine kel (x0.36) and noresketamine apparent CLn/F (x0.81); Japanese race increases the nasal-cavity-absorbed fraction FRn (x1.34); and hepatic blood flow Qh declines linearly by 2.19 L/h per year of age above 60."
  reference   <- paste(
    "Perez-Ruixo C, Rossenu S, Zannikos P, Nandy P, Singh J, Drevets WC,",
    "Perez-Ruixo JJ. Population Pharmacokinetics of Esketamine Nasal Spray",
    "and its Metabolite Noresketamine in Healthy Subjects and Patients with",
    "Treatment-Resistant Depression. Clin Pharmacokinet. 2021;60(4):501-516.",
    "doi:10.1007/s40262-020-00953-4",
    sep = " "
  )
  vignette    <- "PerezRuixo_2020_esketamine"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = non-Asian).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Multiplicative factor e_asian_kel = 0.36 on the esketamine hepatic elimination rate constant kel (Table 3: 'Asian on kel = 0.36'; 64.0% decrease in kel for Asian subjects), and a separate multiplicative factor e_asian_CLn = 0.81 on the noresketamine apparent clearance CLn/F (Table 3: 'Asian on CLn/F = 0.81'; 19.4% decrease). Both effects are multiplicative on the linear parameter scale. Per Table 2, 112 of 820 subjects (13.7%) were Asian (72 Japanese + 40 non-Japanese Asian).",
      source_name        = "RACE_ASIAN"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator (1 = Japanese, 0 = non-Japanese).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese, including Asian non-Japanese)",
      notes              = "Multiplicative factor e_japanese_FRn = 1.34 on the per-spray fraction of nasal dose absorbed via the nasal cavity (Table 3: 'Japanese on FRn = 1.34'; 34% increase). Applied on the linear FRn scale per the abstract text and Results Section. Per Table 2, 72 of 820 subjects (8.8%) were Japanese (out of 112 Asian subjects). Distinct from RACE_ASIAN: a Japanese subject has both RACE_ASIAN = 1 AND RACE_JAPANESE = 1.",
      source_name        = "RACE_JAPANESE"
    ),
    AGE = list(
      description        = "Subject age in years.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Piecewise-linear additive effect on hepatic blood flow Qh: Qh decreases by 2.19 L/h per year of age above 60 (Table 3: 'Age on Qh = -2.19 L/year'; equivalent to 21.9 L/h per decade). For AGE <= 60 the effect is zero. Encoded as Qh = exp(lqh + etalqh) - e_age_qh * max(0, AGE - 60). Per Table 2, AGE range 18-86 years, median 45.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 820L,
    n_studies      = 13L,
    age_range      = "18-86 years (median 45)",
    age_median     = "45 years",
    weight_range   = "39-170 kg (median 74)",
    weight_median  = "74 kg",
    sex_female_pct = 58.4,
    race_ethnicity = c(White = 72.4, Black = 6.82, Asian = 13.7, Other = 7.07),
    disease_state  = "256 (31.2%) healthy volunteers from Phase I single-dose / multi-dose studies; 564 (68.8%) patients with treatment-resistant depression (TRD) enrolled in Phase II and Phase III studies receiving twice-weekly intranasal esketamine 28-112 mg.",
    dose_range     = "Intranasal 14-112 mg single or multiple dose; intravenous 28 mg (study ESKETINTRD1009); oral 84 mg (study ESKETINTRD1009).",
    regions        = "Multi-regional (Phase I: USA, Belgium, Japan; Phase II/III: USA, Europe, Asia).",
    notes          = "Pooled analysis of 13 clinical studies (ESKETINTRD1001-1012, 54135419TRD1015, ESKETINTRD2003, ESKETINSUI2001, ESKETINTRD3001-3005). 9784 esketamine and 9397 noresketamine plasma concentrations. Demographics: Table 2; final parameter estimates: Table 3; covariate-stratified Cmax/AUC: Table 4. Hepatic function biomarkers (ALT, AST, ALP, GGT, LDH, total bilirubin), serum albumin, total protein, eGFR, sex, body weight, and disease state (healthy vs TRD) had no clinically relevant impact on PK (Results Section, r-squared <= 10%, p-value > 0.001) and are not in covariateData; see covariatesDataExcluded for the screened-but-not-retained set."
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained in the final model (Results Section: 'Sex, body weight, ALT, AST, ALP, gamma glutamyl transpeptidase, TB, albumin, TP, estimated glomerular filtration rate, and disease state had no discernable impact on the PK parameters of esketamine and noresketamine')."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Results Section)."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    GGT = list(
      description = "Gamma-glutamyl transpeptidase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Results Section)."
    )
  )

  ini({
    # ------------------- Esketamine absorption -------------------------------
    # Three parallel absorption routes:
    #   (1) Intranasal direct via nasal cavity (depot -> central, first-order
    #       ka_n with bioavailability FRn).
    #   (2) Intranasal swallowed portion via PO depot (depot2 -> liver,
    #       sequential zero-first-order with duration Dsw and rate ka_sw and
    #       bioavailability (1 - FRn) * Fgut).
    #   (3) PO solution via PO depot (depot3 -> liver, sequential zero-first-
    #       order with duration Dpo and rate ka_po and bioavailability Fgut).
    # Section 2.3 of Perez-Ruixo 2020 and Fig. 1 describe the structure.
    logitFRn       <- log(0.54 / (1 - 0.54));  label("Logit-transformed per-spray fraction of nasal dose absorbed via nasal cavity (FRn = 0.54, unitless)")  # Table 3: FRn = 0.54
    e_japanese_FRn <- 1.34;                    label("Multiplicative factor on FRn for Japanese race (unitless)")                                              # Table 3: Japanese on FRn = 1.34
    dose_FRn_effect <- 0.62;                   label("Per-spray FRn reduction factor for the second and subsequent 28-mg sprays of a single nasal dose (unitless)") # Table 3: Dose on FRn = 0.62 (38% decrease for subsequent 28-mg sprays). NOT applied dynamically in this model() block; consumers simulating multi-spray doses (56, 84, 112 mg) must compute the per-dose effective FRn externally and either set N sequential single-spray events or override logitFRn via params (see vignette).
    lka_n          <- log(2.93);               label("Esketamine nasal-cavity first-order absorption rate constant ka_n (1/h)")     # Table 3: ka_n = 2.93 L/h (table units typo; rate constants have units 1/h)
    lDsw           <- log(0.53);               label("Zero-order absorption duration of the swallowed nasal portion Dsw (h)")        # Table 3: Dsw = 0.53 h
    lka_sw         <- log(1.45);               label("Esketamine swallowed-nasal first-order absorption rate constant ka_sw (1/h)") # Table 3: ka_sw = 1.45 L/h (table units typo; rate constants have units 1/h)
    lDpo           <- log(0.32);               label("Zero-order absorption duration of the PO solution dose Dpo (h)")               # Table 3: Dpo = 0.32 (table h/L typo; duration has units h)
    lka_po         <- log(0.97);               label("Esketamine PO solution first-order absorption rate constant ka_po (1/h)")     # Table 3: ka_po = 0.97 L/h (table units typo; rate constants have units 1/h)
    lFgut          <- log(0.64 / (1 - 0.64));  label("Logit-transformed fraction of swallowed (oral or intranasal) dose that escapes intestinal metabolism and reaches the hepato-portal system (Fgut = 0.64, unitless)") # Table 3: Fgut = 0.64

    # ------------------- Esketamine disposition ------------------------------
    # Three-compartment esketamine with a hepato-portal (liver) compartment
    # representing the well-stirred first-pass system. Plasma clearance is
    # implicit: CL = Qh * E where E = Cl_int / (Qh + Cl_int) and
    # Cl_int = Vh * (kel + kmet); the ODE structure produces this
    # automatically through the bidirectional Qh exchange between central
    # and liver plus the kel + kmet first-order eliminations from liver.
    lvc            <- log(192);                label("Esketamine central volume of distribution Vc (L)")              # Table 3: Vc = 192 L
    lq1            <- log(84.3);               label("Esketamine inter-compartmental clearance Q1 to peripheral1 (L/h)")  # Table 3: Q1 = 84.3 L/h
    lvp            <- log(143);                label("Esketamine shallow peripheral volume of distribution Vp1 (L)")    # Table 3: Vp1 = 143 L
    lq2            <- log(37.6);               label("Esketamine inter-compartmental clearance Q2 to peripheral2 (L/h)")  # Table 3: Q2 = 37.6 L/h
    lvp2           <- log(417);                label("Esketamine deep peripheral volume of distribution Vp2 (L)")       # Table 3: Vp2 = 417 L
    lqh            <- log(151);                label("Esketamine hepatic blood flow Qh (L/h)")                          # Table 3: Qh = 151 L/h
    e_age_qh       <- 2.19;                    label("Linear decrement of Qh per year of age above 60 (L/h/year)")      # Table 3: Age on Qh = -2.19 L/year (applied only for AGE > 60)
    lvh            <- log(101);                label("Esketamine hepato-portal volume of distribution Vh (L)")          # Table 3: Vh = 101 L
    lkel           <- log(1.11);               label("Esketamine hepatic elimination rate constant kel to other metabolites (1/h)")  # Table 3: kel = 1.11 L/h (table units typo; rate constants have units 1/h)
    e_asian_kel    <- 0.36;                    label("Multiplicative factor on kel for Asian race (unitless)")                          # Table 3: Asian on kel = 0.36 (64.0% decrease)
    lkmet          <- log(2.77);               label("Esketamine hepatic rate constant of metabolism to noresketamine kmet (1/h)")     # Table 3: kmet = 2.77 L/h (table units typo; rate constants have units 1/h)

    # ------------------- Noresketamine (apparent) ----------------------------
    # Two-compartment apparent metabolite parameterization: because no IV
    # noresketamine data were available, the metabolite PK parameters are
    # 'apparent' in the sense that they are scaled by the unidentifiable
    # fraction of esketamine that survives elimination to become circulating
    # noresketamine. The kmet rate constant on the liver state captures
    # absolute mass flux into central_snk; the apparent CLn/F, Vcn/F,
    # Q3/F, Vpn/F values from Table 3 act on central_snk and peripheral1_snk
    # consistently because (CL_app)/(V_app) = CL_true/V_true (the apparent
    # scaling cancels in the elimination and inter-compartmental rate
    # constants). See Methods Section 2.3 of Perez-Ruixo 2020 for derivation.
    lvcn           <- log(70.0);               label("Noresketamine apparent central volume of distribution Vcn/F (L)")       # Table 3: Vcn/F = 70.0 L
    lcln           <- log(38.0);               label("Noresketamine apparent clearance CLn/F (L/h)")                          # Table 3: CLn/F = 38.0 L/h
    e_asian_CLn    <- 0.81;                    label("Multiplicative factor on CLn/F for Asian race (unitless)")              # Table 3: Asian on CLn/F = 0.81 (19.4% decrease)
    lvpn           <- log(115);                label("Noresketamine apparent peripheral volume of distribution Vpn/F (L)")    # Table 3: Vpn/F = 115 L
    lq3            <- log(26.1);               label("Noresketamine apparent inter-compartmental clearance Q3/F (L/h)")       # Table 3: Q3/F = 26.1 L/h

    # ------------------- Inter-individual variability ------------------------
    # IIV variances on the log scale; omega^2 = log(CV^2 + 1) for CV reported
    # in Table 3. The eta on FRn is additive in the logit domain (Section 2.4
    # of Perez-Ruixo 2020) and the 70.8% in Table 3 is interpreted as the SD
    # of the additive eta on the logit-FRn scale (so omega^2 ~ 0.708^2).
    # IIV on kmet is reported in Section 2.4 prose ("The IIV was quantified
    # using an exponential error model for k_a,n, k_a,po, V_c, V_p1, V_h,
    # Q_h, k_el, k_met, V_cn/F, and CL_n/F ...") but its numeric value
    # is OMITTED from Table 3 (the omega(k_met) Estimate column is blank;
    # publication omission, NOT a structural missing variance). We
    # therefore set etalkmet to a documented placeholder (CV = 30%;
    # omega^2 = log(1 + 0.30^2) = 0.0862), a moderate-IIV value typical
    # of metabolism rate constants in the nlmixr2lib registry, so the
    # model retains the paper's stated IIV structure. Downstream users
    # running stochastic VPCs should refit or override this entry once
    # a numeric value becomes available (e.g. via author correspondence
    # to cperezru@its.jnj.com or a future erratum). See the vignette
    # Errata section for the full discussion and the sidecar trail.
    etalogitFRn ~ 0.501             # Table 3 omega(FRn) = 70.8 (additive on logit); SD_logit ~ 0.708 -> variance 0.501
    etalka_n    ~ 0.319             # Table 3 omega(ka_n) = 61.5 CV%; omega^2 = log(1 + 0.615^2) = 0.319
    etalka_sw   ~ 1.009             # Table 3 omega(ka_sw) = 132 CV%; omega^2 = log(1 + 1.32^2) = 1.009
    etalka_po   ~ 1.009             # Table 3 omega(ka_po) = 132 CV%; same value as ka_sw (Table 3 row gives identical 132 CV%, 9.80 RSE)
    etalvc      ~ 0.0730            # Table 3 omega(Vc) = 27.5 CV%; omega^2 = log(1 + 0.275^2) = 0.0730
    etalvp      ~ 0.218             # Table 3 omega(Vp1) = 49.3 CV%; omega^2 = log(1 + 0.493^2) = 0.218
    etalqh      ~ 0.0524            # Table 3 omega(Qh) = 23.2 CV%; omega^2 = log(1 + 0.232^2) = 0.0524
    etalvh      ~ 0.113             # Table 3 omega(Vh) = 34.6 CV%; omega^2 = log(1 + 0.346^2) = 0.113
    etalkel     ~ 0.892             # Table 3 omega(kel) = 120 CV%; omega^2 = log(1 + 1.2^2) = 0.892
    etalkmet    ~ 0.0862            # Table 3 omega(kmet) reported (Section 2.4) but Estimate column blank; placeholder CV = 30% used (omega^2 = log(1 + 0.30^2) = 0.0862) - see vignette Errata
    etalvcn     ~ 0.0955            # Table 3 omega(Vcn/F) = 31.6 CV%; omega^2 = log(1 + 0.316^2) = 0.0955
    etalcln     ~ 0.0625            # Table 3 omega(CLn/F) = 25.4 CV%; omega^2 = log(1 + 0.254^2) = 0.0625

    # ------------------- Residual error -------------------------------------
    # Paper Section 2.4: 'additive error model after natural logarithmic
    # transformation of the observations and model predictions' (LTBS).
    # Equivalent to log-normal residual; the reported sigma is the SD on
    # the log-transformed scale. Phase I/II values are used (Phase III
    # values are reported separately but Phase I/II is the canonical
    # value for rich-PK simulation).
    expSd     <- 0.276;             label("Log-normal residual SD for esketamine plasma concentration (LTBS, log-units)")     # Table 3: sigma1(esketamine) = 27.6% (Phase I/II)
    expSd_snk <- 0.422;             label("Log-normal residual SD for noresketamine plasma concentration (LTBS, log-units)")  # Table 3: sigma2(noresketamine) = 42.2% (Phase I/II)
  })

  model({
    # -------- Population FRn after logit-additive IIV and Japanese effect ----
    # Section 2.4: 'an additive error model in the logit domain for FRn'.
    # The Japanese effect (Table 3: 1.34) is interpreted as a multiplicative
    # factor on the LINEAR per-spray FRn (matches abstract Results: '34%
    # increase in FRn'). The dose-effect reduction `dose_FRn_effect` (0.62)
    # for subsequent 28-mg sprays is exposed via a `dose_FRn_effect_export`
    # variable below; the model itself uses the per-spray FRn, and the
    # vignette demonstrates how to compute the per-dose effective FRn
    # externally for 56-/84-/112-mg doses by averaging across sprays:
    #   FRn_effective = FRn * (1 + (N_sprays - 1) * dose_FRn_effect) / N_sprays
    FRn_logit <- logitFRn + etalogitFRn
    FRn_per_spray <- 1 / (1 + exp(-FRn_logit))
    FRn <- FRn_per_spray * (e_japanese_FRn ^ RACE_JAPANESE)
    dose_FRn_effect_export <- dose_FRn_effect  # documented in ini(); not applied dynamically

    # -------- Individual structural parameters -------------------------------
    ka_n   <- exp(lka_n   + etalka_n)
    Dsw    <- exp(lDsw)
    ka_sw  <- exp(lka_sw  + etalka_sw)
    Dpo    <- exp(lDpo)
    ka_po  <- exp(lka_po  + etalka_po)
    Fgut   <- 1 / (1 + exp(-lFgut))
    Vc     <- exp(lvc     + etalvc)
    Q1     <- exp(lq1)
    Vp     <- exp(lvp     + etalvp)
    Q2     <- exp(lq2)
    Vp2    <- exp(lvp2)
    Qh     <- (exp(lqh) - e_age_qh * max(0, AGE - 60)) * exp(etalqh)
    Vh     <- exp(lvh     + etalvh)
    kel    <- exp(lkel    + etalkel)  * (e_asian_kel ^ RACE_ASIAN)
    kmet   <- exp(lkmet   + etalkmet)
    Vcn    <- exp(lvcn    + etalvcn)
    CLn    <- exp(lcln    + etalcln)  * (e_asian_CLn ^ RACE_ASIAN)
    Vpn    <- exp(lvpn)
    Q3     <- exp(lq3)

    # -------- Micro-constants ------------------------------------------------
    k_central_p1   <- Q1 / Vc        # central -> peripheral1
    k_p1_central   <- Q1 / Vp        # peripheral1 -> central
    k_central_p2   <- Q2 / Vc        # central -> peripheral2
    k_p2_central   <- Q2 / Vp2       # peripheral2 -> central
    k_central_h    <- Qh / Vc        # central -> liver (hepato-portal)
    k_h_central    <- Qh / Vh        # liver -> central
    ken            <- CLn / Vcn      # noresketamine apparent elimination rate constant
    k_csnk_p1snk   <- Q3 / Vcn       # noresketamine central -> peripheral1
    k_p1snk_csnk   <- Q3 / Vpn       # noresketamine peripheral1 -> central

    # -------- ODE system -----------------------------------------------------
    # Compartments:
    #   depot              : nasal cavity depot (intranasal dose with F = FRn,
    #                        first-order ka_n into central)
    #   depot2             : oral depot loaded by the swallowed portion of an
    #                        intranasal dose (F = (1 - FRn) * Fgut, zero-order
    #                        release over Dsw, first-order ka_sw into liver)
    #   depot3             : oral depot loaded by a PO solution dose
    #                        (F = Fgut, zero-order release over Dpo, first-
    #                        order ka_po into liver)
    #   central, peripheral1, peripheral2 : esketamine 3-compartment plasma
    #   liver              : hepato-portal compartment (Vh, Qh, kel, kmet)
    #   central_snk, peripheral1_snk : noresketamine apparent 2-compartment
    # Encoding note (see vignette Errata): the Section 2.3 printed ODEs
    # for dA1/dt and dA2/dt contain typos -- dA2/dt depends on A1 not A2,
    # so A2 would integrate to negative mass, and dA1/dt is missing the
    # symmetric (1 - FRn) drainage term that would balance the dA6/dt
    # input from the swallowed-nasal route. The standard pharmacometric
    # bioavailability-mechanism encoding below (f() and dur() per depot)
    # reproduces the published F_n = FRn + (1-FRn)*Fgut*(1-E) and
    # F_po = Fgut*(1-E) formulae exactly and matches the Fig. 1 schematic.
    d/dt(depot)            <- -ka_n  * depot
    d/dt(depot2)           <- -ka_sw * depot2
    d/dt(depot3)           <- -ka_po * depot3
    d/dt(central)          <-  ka_n  * depot +
                               k_p1_central * peripheral1 +
                               k_p2_central * peripheral2 +
                               k_h_central  * liver -
                               (k_central_p1 + k_central_p2 + k_central_h) * central
    d/dt(peripheral1)      <-  k_central_p1 * central - k_p1_central * peripheral1
    d/dt(peripheral2)      <-  k_central_p2 * central - k_p2_central * peripheral2
    d/dt(liver)            <-  ka_sw * depot2 +
                               ka_po * depot3 +
                               k_central_h * central -
                               (k_h_central + kel + kmet) * liver
    d/dt(central_snk)      <-  kmet * liver +
                               k_p1snk_csnk * peripheral1_snk -
                               (ken + k_csnk_p1snk) * central_snk
    d/dt(peripheral1_snk)  <-  k_csnk_p1snk * central_snk -
                               k_p1snk_csnk * peripheral1_snk

    # -------- Bioavailability and zero-order durations -----------------------
    # The user issues dosing events to specific compartments:
    #   - Intranasal dose: two events at the same time, both with amt =
    #     nasal-dose-mg. cmt = depot (instantaneous, F = FRn) and cmt =
    #     depot2 (zero-order over Dsw, F = (1 - FRn) * Fgut).
    #   - PO solution dose: cmt = depot3, amt = PO-dose-mg (zero-order over
    #     Dpo, F = Fgut).
    #   - IV dose: cmt = central, amt = IV-dose-mg, with rate / dur per the
    #     prescribed infusion.
    f(depot)    <- FRn
    f(depot2)   <- (1 - FRn) * Fgut
    dur(depot2) <- Dsw
    f(depot3)   <- Fgut
    dur(depot3) <- Dpo

    # -------- Observations ---------------------------------------------------
    # Dose units = mg; volumes in L => Cc = central / Vc has units mg/L
    # = 1000 ng/mL. Multiply by 1000 to express in ng/mL (matches Table 4).
    #
    # NORESKETAMINE AUC0-24 NOTE (see vignette Errata for the full discussion
    # and the sidecar-001 audit trail): the apparent-V denominator Vcn / 1
    # (no F_met scaling) is the most literal reading of Table 3 Vcn/F = 70 L
    # and the printed Section 2.3 equation dA7/dt = kmet*A6 - ken*A7 - ...
    # This faithful encoding reproduces Table 4 esketamine Cmax / AUC0-24
    # within 4% AND the noresketamine Cmax within 2%, but OVER-predicts the
    # noresketamine AUC0-24 by exactly 1.41 = 1/F_met = (k_el + k_met)/k_met
    # = 3.88/2.77 (a deterministic algebraic relationship, not noise). Alan
    # Maloney independently encountered the same discrepancy in his
    # implementation. Downstream users who specifically need to reproduce
    # the published noresketamine AUC0-24 can multiply Cc_snk by F_met
    # externally (a documented alternative; not the default because doing
    # so would under-predict Cmax by the same 1.41x factor). The committed
    # model preserves Table 3 verbatim and surfaces the discrepancy in
    # the vignette Errata so it is visible to every downstream reader.
    Cc     <- central / Vc * 1000
    Cc_snk <- central_snk / Vcn * 1000

    Cc     ~ lnorm(expSd)
    Cc_snk ~ lnorm(expSd_snk)
  })
}
