Seng_2015_isoniazid <- function() {
  description <- paste(
    "Parent + two-metabolite population pharmacokinetic model for oral",
    "isoniazid (INH), acetylisoniazid (AcINH), and isonicotinic acid",
    "(INA) in 33 healthy Asian adults (Seng 2015; Singapore single-dose",
    "300 mg oral INH study with crossover rifampin / efavirenz arms).",
    "Two-compartment INH disposition with first-order absorption,",
    "linked to a two-compartment AcINH disposition and a one-compartment",
    "INA disposition; metabolite formation splits via the fraction-of-",
    "clearance parameters F_AcINH (INH -> AcINH) and F_INA (AcINH ->",
    "INA), with the complementary (1 - F_AcINH) routing INH directly to",
    "INA. The NAT2-derived acetylator phenotype (rapid / intermediate /",
    "slow) selects between three typical-value INH clearances (65.2 /",
    "32.6 / 6.52 L/h at 63 kg). Creatinine clearance enters as a power-",
    "law covariate on AcINH clearance with exponent 0.4 referenced to",
    "the cohort median 113 mL/min. All clearance and volume terms are",
    "allometrically scaled by total body weight (0.75 exponent on",
    "clearance, 1.0 on volume) with reference weight 63 kg. AcINH and",
    "INA central volumes are fixed at 17 L (apparent central volume of",
    "AcINH from Boxenbaum & Riegelman 1976) to keep the metabolite",
    "model identifiable in the absence of intravenous data."
  )
  reference <- paste(
    "Seng K-Y, Hee K-H, Soon G-H, Chew N, Khoo SH, Lee LS-U (2015).",
    "Population pharmacokinetic analysis of isoniazid, acetylisoniazid,",
    "and isonicotinic acid in healthy volunteers.",
    "Antimicrob Agents Chemother 59(11):6791-6799.",
    "doi:10.1128/AAC.01244-15."
  )
  vignette <- "Seng_2015_isoniazid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on all clearance terms with exponent 0.75",
        "and on all volume terms with exponent 1.0, referenced to the",
        "rounded cohort median 63 kg (Seng 2015 Methods page 6792:",
        "'median body weight of 63 kg'; Table 1 reports median 62.5 kg",
        "with range 45.8-86.1 kg). Time-fixed per subject in the",
        "single-dose study window."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-law covariate on AcINH clearance referenced to the",
        "cohort median 113 mL/min (Seng 2015 Table 1 'Creatinine",
        "clearance' median 113.1, range 52.3-174.6 mL/min). Cockcroft",
        "& Gault equation using total body weight per Seng 2015",
        "Methods page 6792-6793. Time-fixed per subject."
      ),
      source_name        = "CRCL"
    ),
    NAT2_RAPID = list(
      description        = paste(
        "NAT2 rapid (fast) acetylator phenotype indicator: 1 = subject",
        "classified as a rapid (fast) NAT2 acetylator, 0 = otherwise",
        "(intermediate or slow). Paired with NAT2_SLOW to encode the",
        "three-level NAT2 phenotype; joint reference (both = 0) is the",
        "intermediate acetylator group."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or slow); joint reference with NAT2_SLOW = 0 is the intermediate acetylator group.",
      notes              = paste(
        "Seng 2015 Methods 'NAT2 genotyping and phenotyping' uses the",
        "Sequenom iPLEX ADME PGx panel to genotype NAT2 SNPs",
        "rs1805158, rs1801279, rs1041983, rs1801280, rs1799929,",
        "rs1799930, rs1799931, rs1208 and classifies subjects into",
        "rapid / intermediate / slow phenotypes using a published",
        "algorithm (ref 16). Cohort distribution (Table 1):",
        "7 rapid / 15 intermediate / 11 slow. Used together with",
        "NAT2_SLOW to select between three typical-value INH",
        "clearances (lcl_fast, lcl_int, lcl_slow) in model()."
      ),
      source_name        = "NAT2"
    ),
    NAT2_SLOW = list(
      description        = paste(
        "NAT2 slow acetylator phenotype indicator: 1 = subject",
        "classified as a slow NAT2 acetylator, 0 = otherwise",
        "(intermediate or rapid). Paired with NAT2_RAPID; joint",
        "reference (both = 0) is the intermediate acetylator group."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or rapid); joint reference with NAT2_RAPID = 0 is the intermediate acetylator group.",
      notes              = paste(
        "Same NAT2 phenotyping source as NAT2_RAPID (Seng 2015",
        "Methods 'NAT2 genotyping and phenotyping'). The Seng 2015",
        "three-level usage (fast / intermediate / slow each with its",
        "own typical CL/F) extends the binary slow-vs-nonslow",
        "convention from Horita 2018; see NAT2_SLOW canonical entry",
        "for the heritage. Cohort distribution (Table 1):",
        "7 rapid / 15 intermediate / 11 slow."
      ),
      source_name        = "NAT2"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "22-56 years",
    age_median     = "33 years",
    weight_range   = "45.8-86.1 kg",
    weight_median  = "62.5 kg",
    sex_female_pct = 30,
    race_ethnicity = c(Chinese = 64, Malay = 21, Indian = 15),
    disease_state  = paste(
      "Healthy Asian adults (Singapore), prior CYP2B6 516 GG (n = 23)",
      "or TT (n = 10) genotype stratification as part of an efavirenz",
      "/ rifampin interaction study. Subjects received a single 300 mg",
      "oral isoniazid dose on day 14 of a 14-day daily rifampin +",
      "isoniazid arm, plus single 600 mg oral efavirenz at the same",
      "timepoint; pharmacokinetic blood sampling at 0, 1, 2, 4, 6, 8,",
      "10, 12, 18, and 24 h after the final INH dose."
    ),
    dose_range     = paste(
      "Single 300 mg oral isoniazid; population pharmacokinetic",
      "simulations in the source paper additionally evaluated a 200 mg",
      "single oral dose for the 30-45 kg body-weight band per WHO",
      "guidance (Seng 2015 Methods page 6793)."
    ),
    regions        = "Singapore",
    nat2_phenotype_distribution = c(rapid = 7L, intermediate = 15L, slow = 11L),
    notes          = paste(
      "Demographics in Table 1 of Seng 2015; concentrations of INH,",
      "AcINH, and INA in plasma were assayed by LC-MS/MS with LLOQs",
      "of 5, 12.5, and 12.5 ng/mL respectively (Methods page 6791).",
      "Total INH / AcINH / INA observations available for modelling:",
      "298 / 309 / 310 (Results page 6794). Acetylator phenotype",
      "counts add to 33 subjects."
    )
  )

  ini({
    # ============================================================
    # All parameter values are taken from Seng 2015 Table 2, FINAL
    # MODEL ESTIMATES column (page 6794). RSE values from the
    # covariance step are quoted in inline comments. The model is
    # parameterised around the cohort-median reference body weight of
    # 63 kg (Methods page 6792) and the cohort-median creatinine
    # clearance of 113 mL/min (Table 1 + Results discussion page
    # 6796). All clearance and volume terms carry allometric scaling
    # (Results page 6794: 'All clearance and volume terms in our
    # population pharmacokinetic model were allometrically scaled by
    # body weight').
    # ============================================================

    # ---- INH absorption and disposition ------------------------------
    lka      <- log(0.6);   label("INH first-order absorption rate constant ka (1/h)")              # Table 2 final ka = 0.6 1/h (RSE 5.3%)

    # NAT2-phenotype-dependent typical INH apparent clearance. Three
    # explicit log-typical values follow Seng 2015's parameterisation
    # TVP_i = theta_1 * (1 + theta_2 * IND_i) with the rapid (fast)
    # group as the source-paper reference (theta_1 = 65.2 L/h) and the
    # tabulated rapid / intermediate / slow effects of 0 / -0.5 / -0.9
    # on CL/F (Table 2 final, last two rows of the structural-parameter
    # block; the magnitudes 0.5 and 0.9 are reported as positive
    # numbers, with the sign implicit in the abstract / Results page
    # 6794: '1.9- and 7.7-fold higher than in intermediate and slow
    # eliminators (65 versus 35 and 8 liters/h)').
    lcl_fast <- log(65.2);          label("INH apparent CL/F for NAT2 rapid acetylators at WT = 63 kg (L/h)")          # Table 2 final CL/F = 65.2 L/h (RSE 6.9%) -- rapid acetylator typical value
    lcl_int  <- log(65.2 * (1 - 0.5)); label("INH apparent CL/F for NAT2 intermediate acetylators at WT = 63 kg (L/h)") # Table 2 final CL/F = 65.2 L/h * (1 - 0.5) = 32.6 L/h; '0.5 (9.1)' = intermediate-on-CL/F reduction
    lcl_slow <- log(65.2 * (1 - 0.9)); label("INH apparent CL/F for NAT2 slow acetylators at WT = 63 kg (L/h)")         # Table 2 final CL/F = 65.2 L/h * (1 - 0.9) = 6.52 L/h; '0.9 (1.8)' = slow-on-CL/F reduction (matches Results: 8 L/h slow vs 65 L/h fast = 7.7-fold ratio)

    lvc      <- log(18);    label("INH apparent central volume V_2/F at WT = 63 kg (L)")             # Table 2 final V_2/F = 18 L (RSE 22.6%)
    lq       <- log(2.8);   label("INH apparent inter-compartmental clearance Q/F at WT = 63 kg (L/h)") # Table 2 final Q/F  = 2.8 L/h (RSE 14.5%)
    lvp      <- log(15.9);  label("INH apparent peripheral volume V_5/F at WT = 63 kg (L)")          # Table 2 final V_5/F = 15.9 L (RSE 11.3%)

    # Oral bioavailability F_INH. Fixed typical value = 1 (Methods page
    # 6791 / 6792: 'F_INH is not estimated due to a lack of pharmaco-
    # kinetic data from intravenous dosing') paired with an estimated
    # log-normal IIV (Table 2 final IIV F_INH 31.6%).
    lfdepot  <- fixed(log(1));  label("INH oral bioavailability F_INH (unitless, fixed at 1)")       # Table 2 final F_INH = 1 (fixed)

    # Logit-transformed fraction of INH clearance forming AcINH.
    # F_AcINH = 0.973 (Table 2 final, RSE 1.3%); the complementary
    # fraction 1 - F_AcINH = 0.027 routes directly to INA.
    logitfr_acinh <- logit(0.973); label("Fraction of INH clearance forming AcINH (proportion, logit scale)")  # Table 2 final F_AcINH = 0.973 (RSE 1.3%)

    # ---- AcINH disposition -----------------------------------------
    lcl_acinh <- log(21.3); label("AcINH clearance CL_A at WT = 63 kg, CRCL = 113 mL/min (L/h)")     # Table 2 final CL_A = 21.3 L/h (RSE 7.2%)
    lvc_acinh <- fixed(log(17)); label("AcINH central volume V_3 (L, fixed)")                       # Table 2 final V_3 = 17 L (fixed); Methods page 6792 / Fig 2 caption: 'fixed to ... 17 liters ... the calculated volume of distribution of the central compartment of AcINH in Boxenbaum and Riegelman'
    lq_acinh  <- log(69.2); label("AcINH inter-compartmental clearance Q_A at WT = 63 kg (L/h)")    # Table 2 final Q_A = 69.2 L/h (RSE 13.2%)
    lvp_acinh <- log(80.4); label("AcINH peripheral volume V_6 at WT = 63 kg (L)")                  # Table 2 final V_6 = 80.4 L (RSE 9.4%)

    # Logit-transformed fraction of AcINH clearance forming INA.
    # F_INA = 0.734 (Table 2 final, RSE 10.1%); the complementary
    # fraction 1 - F_INA = 0.266 is non-INA elimination of AcINH.
    logitfr_ina <- logit(0.734); label("Fraction of AcINH clearance forming INA (proportion, logit scale)")    # Table 2 final F_INA = 0.734 (RSE 10.1%)

    # ---- INA disposition -------------------------------------------
    lcl_ina <- log(44.6);   label("INA clearance CL_I at WT = 63 kg (L/h)")                          # Table 2 final CL_I = 44.6 L/h (RSE 5.6%)
    lvc_ina <- fixed(log(17)); label("INA volume V_4 (L, fixed)")                                  # Table 2 final V_4 = 17 L (fixed); Methods page 6792 / Fig 2 caption: 'fixed to a volume of 17 liters to ensure model identifiability'

    # ---- Allometric exponents --------------------------------------
    # Fixed at the canonical Anderson & Holford 2008 values; Seng 2015
    # Methods page 6792 reports the allometric form without estimating
    # the exponents.
    e_wt_cl  <- fixed(0.75); label("Allometric exponent of body weight on clearance (unitless, fixed)")        # Methods page 6792 allometric formula for CL/F with exponent 0.75
    e_wt_vc  <- fixed(1);    label("Allometric exponent of body weight on volume (unitless, fixed)")           # Methods page 6792 allometric formula for V_2/F with exponent 1

    # ---- Creatinine-clearance covariate on AcINH CL ----------------
    e_crcl_cl_acinh <- 0.4;  label("Power exponent of CRCL on AcINH clearance (unitless)")          # Table 2 final 'CrCL on CL_A' power-function exponent = 0.4 (RSE 21.5%)

    # ---- Inter-individual variability (omega^2 on the log scale) ---
    # Seng 2015 final-model IIVs from Table 2 column '%IIV (% RSE)'.
    # Each CV% maps to internal variance via omega^2 = log(1 + CV^2)
    # (log-normal). No IIV is reported on V_2/F, Q/F, V_3, Q_A, V_4,
    # F_AcINH, or F_INA; the single CL/F IIV (14.4%) is shared across
    # all three NAT2 phenotype groups (Seng 2015 Discussion page 6797:
    # 'the interindividual variability in CL/F cannot be estimated
    # separately for the fast, intermediate, and slow eliminator
    # subgroups').
    etalka         ~ 0.01651  # log(1 + 0.129^2) -- Table 2 final IIV ka = 12.9%
    etalcl         ~ 0.02053  # log(1 + 0.144^2) -- Table 2 final IIV CL/F = 14.4% (post-NAT2 covariate)
    etalvp         ~ 0.11183  # log(1 + 0.344^2) -- Table 2 final IIV V_5/F = 34.4%
    etalfdepot     ~ 0.09518  # log(1 + 0.316^2) -- Table 2 final IIV F_INH = 31.6%
    etalcl_acinh   ~ 0.01055  # log(1 + 0.103^2) -- Table 2 final IIV CL_A = 10.3% (post-CRCL covariate)
    etalvp_acinh   ~ 0.04077  # log(1 + 0.204^2) -- Table 2 final IIV V_6 = 20.4%
    etalcl_ina     ~ 0.03402  # log(1 + 0.186^2) -- Table 2 final IIV CL_I = 18.6%

    # ---- Residual variability --------------------------------------
    # The source paper used 'additive error model in log scale' for
    # all three outputs (Methods page 6791-6792; Table 2 footnote f).
    # In nlmixr2 this maps directly to lnorm() with expSd = the
    # reported SD on the log scale.
    expSd        <- 0.326;  label("INH residual SD on log(mg/L) scale")                              # Table 2 final residual INH = 0.326 (RSE 14.3%)
    expSd_acinh  <- 0.207;  label("AcINH residual SD on log(mg/L) scale")                            # Table 2 final residual AcINH = 0.207 (RSE 24.1%)
    expSd_ina    <- 0.269;  label("INA residual SD on log(mg/L) scale")                              # Table 2 final residual INA = 0.269 (RSE 13.8%)
  })

  model({
    # ---- NAT2-driven typical INH clearance -------------------------
    # Paired NAT2_RAPID + NAT2_SLOW indicators (joint reference =
    # intermediate; both = 0) select between three explicit
    # log-typical CL/F values. The selection reproduces Seng 2015's
    # tabulated typical CL/F per phenotype: rapid 65.2, intermediate
    # 32.6, slow 6.52 L/h at the reference 63 kg.
    lcl_typ <- lcl_fast * NAT2_RAPID +
               lcl_int  * (1 - NAT2_RAPID - NAT2_SLOW) +
               lcl_slow * NAT2_SLOW

    # ---- Individual structural parameters --------------------------
    # All clearance terms scale with (WT/63)^0.75 and all volume terms
    # with (WT/63)^1. AcINH CL also picks up the (CRCL/113)^0.4 power
    # covariate (Seng 2015 Methods page 6792-6793 + Table 2).
    ka         <- exp(lka       + etalka)
    cl         <- exp(lcl_typ   + etalcl)       * (WT / 63)^e_wt_cl
    vc         <- exp(lvc)                       * (WT / 63)^e_wt_vc
    q          <- exp(lq)                        * (WT / 63)^e_wt_cl
    vp         <- exp(lvp       + etalvp)        * (WT / 63)^e_wt_vc
    fdepot     <- exp(lfdepot   + etalfdepot)
    fr_acinh   <- expit(logitfr_acinh)
    cl_acinh   <- exp(lcl_acinh + etalcl_acinh)  * (WT / 63)^e_wt_cl * (CRCL / 113)^e_crcl_cl_acinh
    vc_acinh   <- exp(lvc_acinh)                  * (WT / 63)^e_wt_vc
    q_acinh    <- exp(lq_acinh)                   * (WT / 63)^e_wt_cl
    vp_acinh   <- exp(lvp_acinh + etalvp_acinh)  * (WT / 63)^e_wt_vc
    fr_ina     <- expit(logitfr_ina)
    cl_ina     <- exp(lcl_ina   + etalcl_ina)    * (WT / 63)^e_wt_cl
    vc_ina     <- exp(lvc_ina)                    * (WT / 63)^e_wt_vc

    # ---- Micro-constants (formation and elimination rates) ---------
    # Partial INH-clearance routes split via F_AcINH (forms AcINH) and
    # 1 - F_AcINH (forms INA directly from INH); analogously F_INA
    # splits AcINH-clearance to INA versus non-INA elimination. Q
    # terms set inter-compartmental rates within each parent /
    # metabolite sub-system.
    k_inh_acinh <-       fr_acinh  * cl / vc
    k_inh_ina   <- (1 - fr_acinh) * cl / vc
    k_inh_pc    <- q / vc
    k_inh_cp    <- q / vp
    k_acinh_ina <-       fr_ina   * cl_acinh / vc_acinh
    k_acinh_el  <- (1 - fr_ina)  * cl_acinh / vc_acinh
    k_acinh_pc  <- q_acinh / vc_acinh
    k_acinh_cp  <- q_acinh / vp_acinh
    k_ina_el    <- cl_ina / vc_ina

    # ---- ODE system (Fig 2 of Seng 2015) ---------------------------
    d/dt(depot)             <- -ka * depot
    d/dt(central)           <-  ka * depot -
                                (k_inh_acinh + k_inh_ina) * central -
                                k_inh_pc * central + k_inh_cp * peripheral1
    d/dt(peripheral1)       <-  k_inh_pc * central - k_inh_cp * peripheral1
    d/dt(central_acinh)     <-  k_inh_acinh * central -
                                (k_acinh_ina + k_acinh_el) * central_acinh -
                                k_acinh_pc * central_acinh + k_acinh_cp * peripheral1_acinh
    d/dt(peripheral1_acinh) <-  k_acinh_pc * central_acinh - k_acinh_cp * peripheral1_acinh
    d/dt(central_ina)       <-  k_inh_ina * central +
                                k_acinh_ina * central_acinh -
                                k_ina_el * central_ina

    # ---- Bioavailability -------------------------------------------
    f(depot) <- fdepot

    # ---- Observation variables and residual error ------------------
    # All concentrations in mg/L. Log-additive residual on each output
    # per Methods page 6791-6792 / Table 2 footnote f.
    Cc        <- central        / vc
    Cc_acinh  <- central_acinh  / vc_acinh
    Cc_ina    <- central_ina    / vc_ina

    Cc        ~ lnorm(expSd)
    Cc_acinh  ~ lnorm(expSd_acinh)
    Cc_ina    ~ lnorm(expSd_ina)
  })
}
