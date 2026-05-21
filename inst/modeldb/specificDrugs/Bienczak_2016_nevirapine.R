Bienczak_2016_nevirapine <- function() {
  description <- "One-compartment population PK model for oral nevirapine in African children (Bienczak 2016) with three-transit-compartment absorption, semi-mechanistic well-stirred hepatic extraction (Gordi-style) splitting oral bioavailability into a pre-hepatic component FpreH (age-driven exponential maturation toward an older-child reference fixed to 1) and a hepatic component FH derived algebraically from intrinsic clearance CLint via FH = QH / (QH + fu * CLint), allometric scaling of CLint and Vc to median weight 14.5 kg and of hepatic plasma flow QH to a 70 kg reference, CYP2B6 516G>T | 983T>C metabolizer-status (EM / IM / SM / USM) effects on CLint, and a 29% diurnal-variation cosine on CLint with zenith at noon."
  reference <- paste(
    "Bienczak A, Cook A, Wiesner L, Mulenga V, Kityo C, Kekitiinwa A,",
    "Walker AS, Owen A, Gibb DM, Burger D, McIlleron H, Denti P (2017).",
    "Effect of diurnal variation, CYP2B6 genotype and age on the",
    "pharmacokinetics of nevirapine in African children.",
    "Journal of Antimicrobial Chemotherapy 72(1):190-199.",
    "doi:10.1093/jac/dkw388.",
    "(Open Access; accepted 17 August 2016, published online 2016, print",
    "issue January 2017. The on-disk PDF filename is 'Bienczak_2017_*' for",
    "the print issue year while the paper text reports 'The Author 2016';",
    "this model file follows the print-issue-year convention only for the",
    "filename's task-dispatcher metadata and uses the 2017 doi-resolvable",
    "citation here.)",
    sep = " "
  )
  vignette <- "Bienczak_2016_nevirapine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying / baseline body weight. Drives the Anderson-Holford allometric scaling of intrinsic clearance CLint (exponent 0.75) and central volume Vc (exponent 1.0), with the structural typical-value parameters reported in Bienczak 2016 Table 3 corresponding to the cohort median 14.5 kg (paper Results 'Population pharmacokinetics' paragraph 4 and Table 3 footnote: 'All clearance and volume parameters scaled allometrically to the median weight of 14.5 kg.'). Hepatic plasma flow QH (50 L/h) and liver volume VH (1 L) are also allometrically scaled (CL-exponent 0.75 on QH, V-exponent 1.0 on VH) but use the 70-kg adult reference of Bienczak 2016 Methods 'Structural model' paragraph 1, not the 14.5-kg child median.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the exponential maturation of pre-hepatic bioavailability FpreH per Bienczak 2016 Results 'Population pharmacokinetics' paragraph 4 / Equation (7) in Appendix S1 (Appendix S1 was not on disk at extraction time; the equation is reconstructed from the paper's narrative as FpreH(AGE) = 1 - (1 - FpreH_birth) * exp(-ln(2) / t_half_FpreH * AGE), with FpreH(0) = 0.583, FpreH(infty) = 1, t_half = 1.54 years from Table 3, and FpreH(3.3 y) ~ 0.906 matching the paper's reported 90%). Treated as time-invariant per subject for steady-state PK simulation; in the source cohort age was the baseline value (Table 1 footnote b).",
      source_name        = "AGE"
    ),
    CYP2B6_IM = list(
      description        = "CYP2B6 intermediate-metabolizer indicator (1 = 516GT | 983TT or 516GG | 983TC)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-IM phenotype: EM, SM, or USM; EM is the structural reference when CYP2B6_IM = CYP2B6_SM = CYP2B6_USM = 0)",
      notes              = "First of three sibling binary indicators (CYP2B6_IM / CYP2B6_SM / CYP2B6_USM) encoding the four-level EM / IM / SM / USM CYP2B6 metabolizer phenotype defined by the combined 516G>T (rs3745274) | 983T>C (rs28399499) SNP vector (Bienczak 2016 Methods 'Covariate effects' paragraph 2 phenotype-assignment list). Multiplicative log-additive effect on CLint: -17% relative to EM (Table 3 absolute typical values 3.27 L/h EM vs 2.72 L/h IM; the encoded log-effect coefficient e_cyp2b6_im_cl = log(2.72/3.27) = -0.184). Cohort prevalence 44.6% IM (Bienczak 2016 Table 2; 141 of 319 genotyped CHAPAS-3 children).",
      source_name        = "metabolizer status (IM)"
    ),
    CYP2B6_SM = list(
      description        = "CYP2B6 slow-metabolizer indicator (1 = 516TT | 983TT or 516GT | 983TC)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-SM phenotype: EM, IM, or USM; EM is the structural reference when CYP2B6_IM = CYP2B6_SM = CYP2B6_USM = 0)",
      notes              = "Second of three sibling binary indicators (CYP2B6_IM / CYP2B6_SM / CYP2B6_USM) encoding the four-level EM / IM / SM / USM CYP2B6 metabolizer phenotype. Multiplicative log-additive effect on CLint: -50% relative to EM (Table 3 absolute typical values 3.27 L/h EM vs 1.65 L/h SM; the encoded log-effect coefficient e_cyp2b6_sm_cl = log(1.65/3.27) = -0.684). Cohort prevalence 21.7% SM (Bienczak 2016 Table 2; 70 of 319 genotyped CHAPAS-3 children).",
      source_name        = "metabolizer status (SM)"
    ),
    CYP2B6_USM = list(
      description        = "CYP2B6 ultra-slow-metabolizer indicator (1 = 516GG | 983CC, 983CC homozygote)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-USM phenotype: EM, IM, or SM; EM is the structural reference when CYP2B6_IM = CYP2B6_SM = CYP2B6_USM = 0)",
      notes              = "Third of three sibling binary indicators (CYP2B6_IM / CYP2B6_SM / CYP2B6_USM) encoding the four-level EM / IM / SM / USM CYP2B6 metabolizer phenotype. Multiplicative log-additive effect on CLint: -68% relative to EM (Table 3 absolute typical values 3.27 L/h EM vs 1.04 L/h USM; the encoded log-effect coefficient e_cyp2b6_usm_cl = log(1.04/3.27) = -1.146). Cohort prevalence 0.6% USM (Bienczak 2016 Table 2; 2 of 319 genotyped CHAPAS-3 children). Bienczak 2016 is the first study to quantify the 983CC homozygote (USM) effect on nevirapine clearance; the 983T>C loss-of-function allele is essentially absent from European-ancestry populations and reaches appreciable frequency only in sub-Saharan African cohorts.",
      source_name        = "metabolizer status (USM)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 414L,
    n_studies      = 2L,
    age_range      = "0.3-15.0 years (paediatric)",
    age_median     = "2.92 years",
    weight_range   = "3.5-29.6 kg",
    weight_median  = "12.2 kg (cohort) / 14.5 kg (allometric-scaling reference, Bienczak 2016 Table 3 footnote)",
    sex_female_pct = 47.5,
    race_ethnicity = "African (all 414 patients black African; recruited from Uganda and Zambia per Bienczak 2016 Table 1 footnote)",
    disease_state  = "HIV-1 infection on twice-daily nevirapine-based combination antiretroviral therapy (paediatric). Companion nucleoside reverse-transcriptase inhibitors (NRTI backbone) were abacavir (n = 115), stavudine (n = 191), or zidovudine (n = 114) per Bienczak 2016 Table 1; the NRTI backbone was tested as a covariate on nevirapine PK and was not retained in the final model (Bienczak 2016 Results 'Population pharmacokinetics' paragraph 5: 'No other covariates were identified as significant').",
    dose_range     = "Twice-daily oral nevirapine dosed by WHO weight-band guidelines (2006 guidelines in CHAPAS-1, 2010 guidelines in CHAPAS-3). Paediatric formulations used in CHAPAS-1: Triomune Baby (50 mg nevirapine FDC) and Triomune Junior (100 mg nevirapine FDC); paediatric and adult formulations used in CHAPAS-3 including Triomune Baby / Junior, Duovir-N Baby (50 mg nevirapine FDC), nevirapine 100 mg single-drug, Duovir-N (200 mg nevirapine FDC), and Triomune30 (200 mg nevirapine FDC). Unequal AM / PM splitting of the daily dose: CHAPAS-1 gave the larger dose at night; CHAPAS-3 gave the larger dose in the morning (Bienczak 2016 Methods 'CHAPAS-1' and 'CHAPAS-3' paragraphs).",
    regions        = "Uganda and Zambia (sub-Saharan Africa).",
    notes          = "Pooled cohort from two CHAPAS trials (CHAPAS-1: 84 children with intensive sampling at pre-dose and 1, 2, 4, 6, 8, 12 h after the morning dose; CHAPAS-3: 336 children with sparse sampling at clinic visits, two samples per visit at least 2 h apart). Six patients rolled over from CHAPAS-1 to CHAPAS-3 (414 = 84 + 336 - 6). Final analysis dataset had 3305 plasma nevirapine concentration measurements after exclusion of 246 samples (111 unclear dosage history, 87 visual outliers with |CWRESI| > 3, and 48 below-the-limit-of-quantification samples confirmed by undetectable companion-drug concentrations). Genotypes were available for 324 children (78.3%); the mixture-model imputation for the remaining 96 children (40.7% EM / 49.0% IM / 9.4% SM) is not encoded in this nlmixr2lib model -- this file assumes known genotype and the user supplies the metabolizer indicators directly. Baseline demographics from Bienczak 2016 Table 1; metabolizer-group prevalences from Bienczak 2016 Table 2 row 1."
  )

  ini({
    # =========================================================================
    # Structural disposition (Bienczak 2016 Table 3 'Final parameter estimates').
    # All clearance and volume parameters reported at the cohort median 14.5 kg
    # (Table 3 footnote: 'All clearance and volume parameters scaled
    # allometrically to the median weight of 14.5 kg.'). The reference here is
    # the EM phenotype with FpreH = 1 (older child, no diurnal modulation), and
    # the typical-value intrinsic clearance is 3.27 L/h.
    # =========================================================================
    lcl_em <- log(3.27);   label("Intrinsic clearance CLint at the EM reference (L/h), allometrically scaled to a 14.5 kg child")  # Bienczak 2016 Table 3 row 'CLint EM (L/h) = 3.27 (3.00-3.69)' (bootstrap median, 5th-95th)
    lvc    <- log(21.92);  label("Apparent central volume of distribution Vc (L), allometrically scaled to a 14.5 kg child")        # Bienczak 2016 Table 3 row 'VC (L) = 21.92 (20.24-26.23)'

    # =========================================================================
    # Transit-compartment absorption (Bienczak 2016 Methods 'Structural model'
    # and Table 3). The paper reports both a transit-chain mean transit time
    # MTT = 0.56 h with NTRANS = 3 fixed transit compartments AND a separate
    # absorption rate constant Ka = 0.84 1/h. The two are jointly identified
    # because the implementation places NN transit compartments in series with
    # a final first-order absorption step from the last transit to central:
    # depot -> transit_1 (rate ktr) -> ... -> transit_NN (rate ktr) -> central
    # (rate ka), with ktr = NN / MTT (i.e. MTT covers only the NN transit
    # compartments, not the depot-to-transit_1 step nor the last-transit-
    # to-central absorption). Appendix S1 of Bienczak 2016 (which would
    # disambiguate the parameterization) was not on disk at extraction time;
    # see vignette Assumptions and deviations for the chosen interpretation.
    # =========================================================================
    lmtt   <- log(0.56);   label("Mean transit time MTT through the NTRANS-compartment chain (h)")           # Bienczak 2016 Table 3 row 'MTT (h) = 0.56 (0.49-0.70)'
    lka    <- log(0.84);   label("First-order absorption rate from the last transit to central (1/h)")      # Bienczak 2016 Table 3 row 'Ka (1/h) = 0.84 (0.67-1.12)'
    nn_fix <- fixed(3);    label("Number of Savic-style transit compartments (integer, unitless)")           # Bienczak 2016 Table 3 row 'NTRANS (number) = 3 (fixed)'; paper Methods 'Structural model' paragraph 2 / Table 3 footnote: 'The number of transit compartments was first estimated and then fixed during the covariate analysis in order to improve model stability. The number was then re-estimated in the final model and proved not to be different from that previously fixed.'

    # =========================================================================
    # Well-stirred hepatic-extraction parameters (Bienczak 2016 Methods
    # 'Structural model' paragraph 1; semi-physiological model from Gordi et
    # al. 2003 reference 40). The three physiological constants are reported
    # for a typical 70 kg adult and allometrically scaled: fu is unitless
    # (no scaling), QH scales with allometric exponent 0.75, VH with 1.0.
    # =========================================================================
    fu_fix  <- fixed(0.40);  label("Nevirapine fraction unbound in plasma fu (unitless)")                # Bienczak 2016 Methods 'Structural model' paragraph 1: 'nevirapine fraction unbound in plasma (fu) 40%' (reference 41)
    lqh_70  <- fixed(log(50));  label("Hepatic plasma flow QH at the 70 kg adult reference (L/h)")       # Bienczak 2016 Methods 'Structural model' paragraph 1: 'hepatic plasma flow (QH) 50 L/h' (reference 42); allometrically scaled with WT
    lvh_70  <- fixed(log(1));   label("Liver volume VH at the 70 kg adult reference (L)")                # Bienczak 2016 Methods 'Structural model' paragraph 1: 'liver volume (VH) 1 L' (reference 40); allometrically scaled with WT (reported for transparency; VH does not enter the well-stirred extraction algebra below)

    # =========================================================================
    # Allometric exponents. Fixed by convention (Anderson and Holford,
    # reference 43 of the paper); the paper Methods 'Covariate effects'
    # paragraph 1: 'Allometric scaling was added to the model at an early
    # stage (before covariate testing), as suggested by Anderson and Holford,
    # and applied to all clearance and volume parameters.'
    # =========================================================================
    e_wt_cl <- fixed(0.75); label("Allometric exponent on intrinsic clearance CLint and hepatic plasma flow QH with body weight (unitless)") # Bienczak 2016 Methods 'Covariate effects' paragraph 1 (carried from Anderson and Holford 2008, ref 43)
    e_wt_vc  <- fixed(1.0);  label("Allometric exponent on central volume Vc and liver volume VH with body weight (unitless)")                # Bienczak 2016 Methods 'Covariate effects' paragraph 1 (carried from Anderson and Holford 2008, ref 43)

    # =========================================================================
    # FpreH parameters (Bienczak 2016 Results 'Population pharmacokinetics'
    # paragraph 4 + Table 3). The reference FpreH for an older child is fixed
    # to 1.0 (Table 3 row 'older children = 1 (fixed)'); FpreH at birth is
    # 58.30% of that, and the gap to the older-child plateau decays with a
    # half-life of 1.54 years (Table 3 row 't_1/2 (years) = 1.54'). The
    # functional form is reconstructed from the narrative as
    #     FpreH(AGE) = 1 - (1 - fpreh_birth) * exp(-ln(2) / t_half * AGE)
    # because Appendix S1 was not on disk at extraction time; the recon-
    # struction reproduces the paper's reported 90% FpreH at age 3.3 years
    # (FpreH(3.3) = 1 - (1 - 0.583) * exp(-0.450 * 3.3) = 0.906, paper says
    # 'approximately 90%'). See vignette Assumptions and deviations.
    # =========================================================================
    lfpreh_old <- fixed(log(1.0));  label("Pre-hepatic bioavailability FpreH in older children (unitless, fixed at 1)")  # Bienczak 2016 Table 3 row 'FpreH older children = 1 (fixed)'
    fpreh_birth <- 0.583;           label("Pre-hepatic bioavailability FpreH at birth, as a fraction of the older-child reference (unitless)") # Bienczak 2016 Table 3 row 'FpreH at birth (%) = 58.30 (50.48-68.24)'
    t_half_fpreh <- 1.54;           label("Half-life of the exponential FpreH approach to the older-child plateau (years)") # Bienczak 2016 Table 3 row 't_1/2 (years) = 1.54 (1.47-2.58)'

    # =========================================================================
    # Diurnal-variation cosine on CLint (Bienczak 2016 Results 'Population
    # pharmacokinetics' paragraph 2 + Table 3). Implementation:
    #     diurnal(t) = 1 + AMP * cos(2 * pi * (t - SHIFT_PEAK) / 24)
    # where SHIFT_PEAK = -12.30 h is the offset of the cosine zenith from
    # midnight; modulo 24 h that places the zenith at 11.70 h (~ noon),
    # matching the paper's 'cosine function with peak amplitude 29% at 12
    # noon' (Abstract Results) and 'cosine function with zenith around 12
    # noon and amplitude of /difference 29%' (Results 'Population
    # pharmacokinetics' paragraph 2). The model assumes simulation time t = 0
    # corresponds to midnight (clock time 00:00); see vignette Assumptions
    # and deviations.
    # =========================================================================
    amp_diurnal       <- 0.292;   label("Amplitude AMP of the diurnal cosine on CLint (fraction)")  # Bienczak 2016 Table 3 row 'AMP (%) = 29.2 (27.7-45.2)'
    shift_peak_diurnal <- -12.30; label("Offset SHIFT of the diurnal cosine zenith from midnight (h, in clock-time coordinates with t=0 at midnight)") # Bienczak 2016 Table 3 row 'SHIFT (h) = -12.30 (-13.32 to -10.38)'; modulo 24 h the zenith is at 11.70 h (~ noon)

    # =========================================================================
    # CYP2B6 metabolizer-status effects on CLint (Bienczak 2016 Results
    # 'Population pharmacokinetics' paragraph 3 + Table 3). Encoded as a
    # log-additive multiplicative effect with three binary indicators
    # (CYP2B6_IM / CYP2B6_SM / CYP2B6_USM); EM is the reference when all
    # three indicators are 0. Coefficients computed from the paper's
    # absolute CLint typical values: e.g. e_im = log(2.72/3.27) = -0.184
    # (paper's reported '17% lower' for IM relative to EM, since
    # 1 - exp(-0.184) = 0.168 ~ 17%).
    # =========================================================================
    e_cyp2b6_im_cl  <- -0.184;  label("Log-additive effect of CYP2B6 IM phenotype on CLint (unitless)")  # Bienczak 2016 Table 3 / Results 'Population pharmacokinetics' paragraph 3: IM = 17% lower than EM; log(2.72/3.27) = -0.184
    e_cyp2b6_sm_cl  <- -0.684;  label("Log-additive effect of CYP2B6 SM phenotype on CLint (unitless)")  # Bienczak 2016 Table 3 / Results 'Population pharmacokinetics' paragraph 3: SM = 50% lower than EM; log(1.65/3.27) = -0.684
    e_cyp2b6_usm_cl <- -1.146;  label("Log-additive effect of CYP2B6 USM phenotype on CLint (unitless)") # Bienczak 2016 Table 3 / Results 'Population pharmacokinetics' paragraph 3: USM = 68% lower than EM; log(1.04/3.27) = -1.146

    # =========================================================================
    # Inter-individual variability (Bienczak 2016 Table 3 'Variability (%)'
    # column). Source reports approximate %CV on the SD scale (footnote b);
    # nlmixr2lib uses variances on the log scale, converted via
    # omega^2 = log(1 + CV^2). The source distinguishes BSV (between-subject)
    # from BOV (between-occasion). nlmixr2lib has no idiomatic encoding for
    # BOV separate from BSV; per the convention used in Svensson_2018_
    # bedaquiline.R and Svensson_2016_rifampicin.R, BOV is dropped where a
    # BSV term is reported on the same parameter and folded in as a
    # BSV-equivalent where only BOV is reported. Specifically:
    #   - CLint: BSV 21.4% kept; no BOV reported
    #   - FpreH: BSV 18.7% kept; BOV 17.0% dropped (see vignette
    #     Assumptions and deviations)
    #   - MTT: only BOV 199.7% reported -> folded in as BSV-equivalent
    #   - Ka: only BOV 44.9% reported -> folded in as BSV-equivalent
    # =========================================================================
    etalcl_em ~ 0.04473  # Bienczak 2016 Table 3 BSV CLint = 21.40%; omega^2 = log(1 + 0.214^2) = 0.04473
    etalfpreh_old ~ 0.03439  # Bienczak 2016 Table 3 BSV FpreH = 18.72%; omega^2 = log(1 + 0.1872^2) = 0.03439 (BOV FpreH = 17.02% dropped per ini() comment above)
    etalmtt   ~ 1.39135  # Bienczak 2016 Table 3 BOV MTT  = 199.73% folded in as BSV-equivalent; omega^2 = log(1 + 1.9973^2) = 1.39135 (no separate BSV on MTT reported)
    etalka    ~ 0.18458  # Bienczak 2016 Table 3 BOV Ka   = 44.91% folded in as BSV-equivalent; omega^2 = log(1 + 0.4491^2) = 0.18458 (no separate BSV on Ka reported)

    # =========================================================================
    # Residual error (Bienczak 2016 Table 3 rows 'Additive error (mg/L)' and
    # 'Proportional error (%)'). The paper additionally reports a 1.56-fold
    # increase in the proportional error for sparse data (rows 'Increased
    # error for sparse data = 1.56 (1.49-1.81)'); this per-occasion / per-
    # observation residual-error scaling is dropped here, leaving the typical
    # (observed-intake) residual-error magnitude. See vignette Assumptions
    # and deviations.
    # =========================================================================
    addSd  <- 0.32;    label("Additive residual error (mg/L)")                       # Bienczak 2016 Table 3 row 'Additive error (mg/L) = 0.32 (0.21-0.38)'
    propSd <- 0.0526;  label("Proportional residual error (fraction)")               # Bienczak 2016 Table 3 row 'Proportional error (%) = 5.26 (4.26-6.18)'
  })

  model({
    # --- 1. Derived covariate terms ----------------------------------------
    # FpreH age-driven maturation (reconstructed exponential form; see ini()
    # comment block on FpreH parameters). Reaches 90% at age 3.3 y and
    # plateaus at 1.0 for older children.
    fpreh_age <- 1 - (1 - fpreh_birth) * exp(-log(2) / t_half_fpreh * AGE)

    # CYP2B6 metabolizer effect on CLint (log-additive composition; EM is the
    # reference when all three indicators are 0).
    cl_meta <- exp(e_cyp2b6_im_cl  * CYP2B6_IM  +
                   e_cyp2b6_sm_cl  * CYP2B6_SM  +
                   e_cyp2b6_usm_cl * CYP2B6_USM)

    # Diurnal cosine on CLint. t is in hours from midnight (simulation
    # convention: anchor t = 0 to clock midnight 00:00); modulo 24, the
    # zenith sits at clock 11.70 h ~ noon.
    diurnal_cl <- 1 + amp_diurnal * cos(2 * 3.141592653589793 * (t - shift_peak_diurnal) / 24)

    # --- 2. Individual PK parameters ---------------------------------------
    # Intrinsic clearance: typical EM CLint scaled allometrically to body
    # weight (reference 14.5 kg), modulated by metabolizer phenotype and the
    # diurnal cosine, with BSV on the log scale.
    clint <- exp(lcl_em + etalcl_em) *
             (WT / 14.5)^e_wt_cl *
             cl_meta *
             diurnal_cl

    # Central volume of distribution: allometric on body weight, no IIV
    # reported in Bienczak 2016 Table 3.
    vc <- exp(lvc) * (WT / 14.5)^e_wt_vc

    # Physiological constants (allometrically scaled from the 70-kg adult
    # reference of the well-stirred liver model).
    qh <- exp(lqh_70) * (WT / 70)^e_wt_cl
    vh <- exp(lvh_70) * (WT / 70)^e_wt_vc  # reported for transparency; not used in the algebraic well-stirred extraction below

    # Pre-hepatic bioavailability: older-child reference * age-driven
    # maturation, with BSV on the log scale (paper estimates an IIV term on
    # FpreH; treated here as a log-normal random effect on the maturation-
    # adjusted typical value).
    fpreh <- exp(lfpreh_old + etalfpreh_old) * fpreh_age

    # Hepatic bioavailability (well-stirred liver, Gordi 2003 form):
    #     FH = QH / (QH + fu * CLint)
    # and hepatic systemic clearance:
    #     CL_H = QH * fu * CLint / (QH + fu * CLint)
    # The apparent oral clearance is CL_H / (FpreH * FH) = fu * CLint / FpreH,
    # consistent with Bienczak 2016 Results 'Population pharmacokinetics'
    # paragraph 4: 'their values of oral clearance (CLoral) were 1.31 L/h EM
    # (reference)' for a 14.5 kg, 4.1 y, FpreH = 0.93 child, which matches
    # 0.40 * 3.27 / 0.93 ~ 1.41 L/h (~ 1.31 L/h after the diurnal-time-
    # average correction). The implementation feeds the dose into central
    # via the bioavailability fraction FpreH * FH and clears central at
    # rate CL_H / Vc, collapsing the explicit liver compartment into
    # algebraic adjustments (a steady-state-equivalent approximation valid
    # at the long ~25-30 h nevirapine half-life and the small VH ~ 1 L /
    # 70 kg liver volume; see vignette Assumptions and deviations).
    fh   <- qh / (qh + fu_fix * clint)
    cl_h <- qh * fu_fix * clint / (qh + fu_fix * clint)

    # --- 3. Absorption micro-constants -------------------------------------
    # NTRANS = 3 fixed transit compartments with shared inter-transit rate
    # ktr = NTRANS / MTT (the paper reports MTT and Ka as separate joint-
    # identified parameters; the chosen interpretation is that MTT covers
    # only the NN-compartment transit chain and Ka governs the additional
    # first-order absorption from the last transit into central; see ini()
    # comment block on transit-compartment absorption).
    mtt <- exp(lmtt + etalmtt)
    ka  <- exp(lka  + etalka)
    nn  <- nn_fix
    ktr <- nn / mtt

    # Micro-constant for systemic elimination.
    kel <- cl_h / vc

    # --- 4. ODE system -----------------------------------------------------
    # depot -> transit_1 -> transit_2 -> transit_3 -> central -> elimination
    # with transit-chain rate ktr (= NTRANS / MTT) on the depot-to-transit
    # and inter-transit steps, separate first-order rate ka from the last
    # transit to central, and systemic elimination via the well-stirred-
    # liver-derived CL_H / Vc rate constant.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ka  * transit3
    d/dt(central)  <-  ka  * transit3 - kel * central

    # --- 5. Bioavailability ------------------------------------------------
    # Total oral bioavailability = pre-hepatic * hepatic (well-stirred-liver
    # first-pass extraction). Applied at the entry to the absorption chain
    # via f(depot).
    f(depot) <- fpreh * fh

    # --- 6. Observation and residual error ---------------------------------
    # Dose in mg / volume in L -> mg/L (the paper's reported unit). Combined
    # additive (0.32 mg/L) and proportional (5.26%) residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
