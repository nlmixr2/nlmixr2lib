Thoueille_2026_salmeterol <- function() {
  description <- paste(
    "Joint plasma and urine population PK model for inhaled salmeterol and",
    "its major metabolite alpha-hydroxysalmeterol in healthy participants,",
    "chronic asthmatics and athletes / endurance-trained individuals",
    "(Thoueille 2026, WADA doping-control analysis). Intravenous-like bolus",
    "absorption into a two-compartment parent disposition; complete",
    "irreversible systemic conversion to alpha-hydroxysalmeterol in a",
    "metabolite plasma compartment whose volume is constrained equal to the",
    "parent central volume; separate cumulative urine compartments for",
    "parent and metabolite fed by first-order urinary excretion rate",
    "constants; and a cumulative urine-volume compartment driven by a",
    "constant estimated urine-production rate that approximates physiologic",
    "micturition. Athletes / endurance-trained individuals carry a 63%",
    "higher salmeterol plasma clearance and a 191% higher salmeterol",
    "urinary excretion rate constant."
  )
  reference <- paste(
    "Thoueille P, Danion A, Hostrup M, Petrou M, Deventer K, Buclin T,",
    "Girardin FR, Mazzoni I, Rabin O, Guidi M. Pharmacometric-Based",
    "Evaluation of Salmeterol and Its Metabolite alpha-Hydroxysalmeterol in",
    "Plasma and Urine: Practical Implications for Doping Control.",
    "CPT Pharmacometrics Syst Pharmacol. 2026. doi:10.1002/psp4.70187.",
    "Final NONMEM control stream and data-formatting description in Data S1."
  )
  vignette <- "Thoueille_2026_salmeterol"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Table 2 regroups the seven urine residual errors by magnitude tier
  # (low / mid / high) rather than by study, so the stratum suffixes below
  # mirror Table 2's own row labels. The study -> tier mapping is in the
  # STUDY_SALM covariate notes and in the model() comments.
  paper_specific_residual_sds <- c(
    "expSdUrineSalLow", "expSdUrineSalMid", "expSdUrineSalHigh",
    "expSdUrineOhsalMid", "expSdUrineOhsalHigh"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Data S1 $MODEL (NCOMP=6, compartment
  # names PARENT_PLASMA / PLASMA_PERIPH / METAB_PLASMA / PARENT_URINE /
  # METAB_URINE / UR_PROD) and Figure 2 of the paper.
  compartmentData <- list(
    central      = list(analyte = "salmeterol", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1  = list(analyte = "salmeterol", units = "nmol", specimen = "plasma", verified = TRUE),
    central_ohsal = list(analyte = "alpha-hydroxysalmeterol", units = "nmol", specimen = "plasma", verified = TRUE),
    urine        = list(analyte = "salmeterol", units = "nmol", specimen = "urine", verified = TRUE),
    urine_ohsal  = list(analyte = "alpha-hydroxysalmeterol", units = "nmol", specimen = "urine", verified = TRUE),
    urine_vol    = list(analyte = "not applicable", units = "L", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    ATHLETE = list(
      description        = "Athlete or healthy endurance-trained individual indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participants pooled with chronic asthmatics)",
      notes              = paste(
        "1 = athlete or healthy endurance-trained individual; 0 = the pooled",
        "reference group. The paper first fit a rich model with a separate PK",
        "parameter per individual type (healthy / chronic asthmatic /",
        "athlete-endurance-trained). Chronic asthmatics had the lowest CL_S/F",
        "and the slowest k14 but were not significantly different from healthy",
        "participants in either univariate or multivariate analysis, so the two",
        "groups were pooled into the reference category of a reduced model",
        "(Results 3.1.1 and 3.1.2). Drives CL_S/F x 1.63 (dOFV = -14, p < 0.001)",
        "and the salmeterol urinary excretion rate constant x 2.91",
        "(dOFV = -55, p < 0.0001). No effect of administration device (MDI vs",
        "DPI) was retained on either parameter. Data S1 $PK derives the",
        "indicator as ATHLETES = 1 when the source column TYPE == 3."
      ),
      source_name        = "TYPE (category 3 of healthy / asthmatic / athlete-endurance-trained)"
    ),
    USG_CORRECTED = list(
      description        = "Urine-specific-gravity-corrected urine concentration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (urine concentration not corrected for specific gravity)",
      notes              = paste(
        "1 = the urine concentration was corrected for urine specific gravity",
        "(USG) per the WADA technical document, using",
        "corrected = (1.020 - 1) / (USG + 0.002 - 1) x observed (Methods 2.1);",
        "0 = no USG value was recorded for the sample so no correction was",
        "applied. Gates the inter-individual variability on the urine-production",
        "rate only: Data S1 $PK branches UR_PROD onto ETA(9) when USG == 0 and",
        "onto ETA(10) when USG == 1, and $OMEGA fixes ETA(10) to 0. Correcting",
        "for USG removes the hydration-status component of urine-production",
        "variability, so no residual IIV was supported on the corrected branch",
        "(Results 3.1.2). USG correction was possible for 42% of salmeterol and",
        "50% of alpha-hydroxysalmeterol urine concentrations (Results 3).",
        "The paper's own Figure 3 / Table 3 / Table S1 simulations report",
        "USG-corrected concentrations, i.e. USG_CORRECTED = 1."
      ),
      source_name        = "USG (0/1 flag in the NONMEM $INPUT record)"
    ),
    STUDY_SALM = list(
      description        = "Source-study identifier in the Thoueille 2026 pooled salmeterol analysis",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "n/a -- selects a residual-error magnitude tier only",
      notes              = paste(
        "Data S1 $INPUT column SID. Values used by the final model:",
        "1 = Jacobson & Hostrup 2022 (urine, DPI, healthy + athletes);",
        "2 = Jessen 2021 (plasma + urine with recorded urine volumes,",
        "endurance-trained); 4 = Jacobson 2017 enantioselective (urine);",
        "5 = Hostrup 2012 (plasma + urine, healthy + asthmatics);",
        "6 = Deventer 2011 (urine with recorded urine volumes, healthy);",
        "7 = Petrou et al. unpublished (plasma + urine, healthy).",
        "SID 3 (Bozzolino 2019) was excluded from the analysis for unreliable",
        "dose / collection timing (Methods 2.1). The indicator selects only the",
        "urine residual-error magnitude tier of Table 2: parent urine uses the",
        "low tier for SID 1 / 5 / 6, the mid tier for SID 2 / 7 and the high",
        "tier for SID 4; metabolite urine uses the mid tier for SID 1 / 5 / 7",
        "and the high tier for SID 2. Plasma residual errors are common to all",
        "studies. Set STUDY_SALM to any of 1 / 5 / 6 to simulate with the",
        "lowest reported urine residual magnitudes."
      ),
      source_name        = "SID"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 85,
    n_studies      = 6,
    age_range      = "Adult (per-study demographics not reported in the pooled analysis)",
    weight_range   = "Not reported",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Pooled healthy participants, chronic asthmatics, and healthy",
      "endurance-trained individuals / athletes. All assumed to have normal",
      "renal function (Methods 2.2.1)."
    ),
    dose_range     = paste(
      "Inhaled salmeterol 50-400 ug: 50 ug and 100 ug MDI, 100 ug / 200 ug /",
      "400 ug DPI, single dose, plus a 200 ug 7-day self-administration arm",
      "(Table 1)."
    ),
    regions        = "Australia, Belgium, Cyprus, Denmark (study locations of the pooled sources)",
    notes          = paste(
      "Six studies pooled (Table 1): Jacobson & Hostrup 2022 (7 healthy +",
      "14 athletes), Jessen 2021 (11 endurance-trained), Jacobson 2017",
      "(10 healthy), Hostrup 2012 (10 healthy + 10 asthmatics), Deventer 2011",
      "(6 healthy), and Petrou et al. unpublished (24 healthy). Bozzolino 2019",
      "(10 asthmatics) was removed for unreliable dose / sampling times.",
      "92 individuals contributed 1175 concentrations; after M1 handling of",
      "209 below-quantification-limit values, 966 concentrations from 85",
      "individuals were analysed (Results 3). Median 6 observations per",
      "individual (range 1-58). Demographic covariates (weight, age, sex) were",
      "not available across the pooled studies, so only individual type and",
      "inhalation device could be screened (Methods 2.2.1)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All values are the final estimates of the
    # Data S1 $THETA block, which carries more significant digits than the
    # rounded Table 2 column; both are noted per line. Parameters are
    # apparent (/F) throughout: the analysis had no intravenous data
    # (Methods 2.2.1).
    # ---------------------------------------------------------------------
    lvc  <- log(446);   label("Apparent central volume of salmeterol, V1/F (L)")                 # $THETA 1 = 446; Table 2 V1/F 446 (RSE 17%)
    lq   <- log(1490);  label("Apparent intercompartmental clearance, Q/F (L/h)")                # $THETA 2 = 1490; Table 2 Q/F 1490 (RSE 22%)
    lvp  <- log(871);   label("Apparent peripheral volume of salmeterol, V2/F (L)")              # $THETA 3 = 871; Table 2 V2/F 871 (RSE 12%)
    lcl  <- log(193);   label("Apparent plasma clearance of salmeterol, CL_S/F (L/h)")           # $THETA 4 = 193; Table 2 CL_S/F 193 (RSE 10%)

    lcl_ohsal <- log(233); label("Apparent plasma clearance of alpha-hydroxysalmeterol, CL_alpha/F (L/h)") # $THETA 5 = 233; Table 2 CL_alpha/F 233 (RSE 22%)

    lkmet <- log(0.3); label("Plasma conversion rate constant of salmeterol to alpha-hydroxysalmeterol, k13 (1/h)") # $THETA 6 = 0.3; Table 2 k13 0.30 (RSE 16%)

    lkurine       <- log(0.000943); label("Urinary excretion rate constant of salmeterol, k14 (1/h)")                   # $THETA 7 = 0.000943; Table 2 k14 0.00094 (RSE 22%)
    lkurine_ohsal <- log(0.0147);   label("Urinary excretion rate constant of alpha-hydroxysalmeterol, k35 (1/h)")      # $THETA 8 = 0.0147; Table 2 k35 0.015 (RSE 22%)

    lurprod <- log(0.079); label("Urine production rate, UR_PROD (L/h)")                          # $THETA 9 = 0.079; Table 2 UR_PROD 0.079 (RSE 18%); Discussion: 79 mL/h, consistent with a typical 1-2 L daily urine volume

    # ---------------------------------------------------------------------
    # Covariate effects. Data S1 $PK writes both as a multiplicative power
    # of a binary indicator, TVCLP = THETA(4) * THETA(10)**ATHLETES and
    # TVK14 = THETA(7) * THETA(11)**ATHLETES, i.e. the theta is the
    # athlete-to-reference ratio itself.
    # ---------------------------------------------------------------------
    e_athlete_cl     <- 1.63; label("Ratio of salmeterol plasma clearance in athletes to the reference group (unitless)")            # $THETA 10 = 1.63; Table 2 theta_Athletes on CL_S 1.63 (RSE 21%); Results 3.1.1: 63% higher CL_S/F
    e_athlete_kurine <- 2.91; label("Ratio of the salmeterol urinary excretion rate constant in athletes to the reference group (unitless)") # $THETA 11 = 2.91; Table 2 theta_Athletes on k14 2.91 (RSE 15%); Results 3.1.2: 191% higher k14

    # ---------------------------------------------------------------------
    # Inter-individual variability. Data S1 $OMEGA reports variances on the
    # log scale; Table 2 reports them as CV% = sqrt(exp(omega^2) - 1) x 100
    # (Table 2 footnote). Each line below carries that back-calculation so
    # the Table 2 column can be checked against the stream.
    # ---------------------------------------------------------------------
    etalvc      ~ 0.0262   # $OMEGA 1 IIV V1;      sqrt(exp(0.0262) - 1) = 16.3% = Table 2 omega_V1 16 (RSE 54%)
    etalq       ~ 0.554    # $OMEGA 2 IIV Q;       sqrt(exp(0.554)  - 1) = 86.0% = Table 2 omega_Q 86 (RSE 37%)
    etalvp      ~ 0.171    # $OMEGA 3 IIV V2;      sqrt(exp(0.171)  - 1) = 43.2% = Table 2 omega_V2 43 (RSE 29%)
    etalcl      ~ 0.103    # $OMEGA 4 IIV CLP;     sqrt(exp(0.103)  - 1) = 32.9% = Table 2 omega_CL_S 33 (RSE 27%)
    etalkmet    ~ 0.157    # $OMEGA 6 IIV K13;     sqrt(exp(0.157)  - 1) = 41.2% = Table 2 omega_k13 41 (RSE 23%)
    etalkurine  ~ 0.081    # $OMEGA 7 IIV K14;     sqrt(exp(0.081)  - 1) = 29.0% = Table 2 omega_k14 29 (RSE 41%)
    etalurprod  ~ 0.429    # $OMEGA 9 IIV UR_PROD not corrected; sqrt(exp(0.429) - 1) = 73.2% = Table 2 omega_UR_PROD USG-uncorrected 73 (RSE 19%)

    # $OMEGA 5 (IIV on CL_alpha) and $OMEGA 8 (IIV on k35) are both '0 FIX',
    # and Table 2 reports no omega row for either: the metabolite plasma data
    # came from a single study and lacked the granularity to estimate them
    # (Discussion). $OMEGA 10 (IIV on UR_PROD for USG-corrected records) is
    # likewise '0 FIX' and is carried structurally by the USG_CORRECTED gate
    # in model() rather than as a second eta. All three are omitted rather
    # than written as `~ fixed(0)` because a zero-variance diagonal makes
    # OMEGA singular and breaks the Cholesky sampler used by rxSolve.

    # ---------------------------------------------------------------------
    # Residual unexplained variability. Data S1 $ERROR applies
    # Y = IPRED * EXP(ERR(n)) to every endpoint, i.e. log-normal residual
    # error, so each SD below is sqrt of the corresponding $SIGMA variance
    # and maps to `~ lnorm(...)`. Table 2 reports the same quantity as
    # CV% = sqrt(sigma^2) x 100 (Table 2 note), so SD x 100 reproduces the
    # published CV% directly.
    # ---------------------------------------------------------------------
    expSd       <- 0.21471; label("Log-normal residual SD for plasma salmeterol (fraction)")                # $SIGMA BLOCK(2) element 1,1 = 0.0461; sqrt = 0.21471 = Table 2 sigma plasma salmeterol 22% (RSE 13%)
    expSd_ohsal <- 0.38079; label("Log-normal residual SD for plasma alpha-hydroxysalmeterol (fraction)")   # $SIGMA BLOCK(2) element 2,2 = 0.145;  sqrt = 0.38079 = Table 2 sigma plasma alpha-hydroxysalmeterol 38% (RSE 13%)

    expSdUrineSalLow    <- 0.29103; label("Log-normal residual SD for urine salmeterol, low-variability studies (fraction)")            # $SIGMA 3 = 0.0847; sqrt = 0.29103 = Table 2 sigma urine salmeterol (low) 29% (RSE 19%); applies to SID 1, 5, 6
    expSdUrineSalMid    <- 0.40743; label("Log-normal residual SD for urine salmeterol, medium-variability studies (fraction)")         # $SIGMA BLOCK(2) element 1,1 = 0.166; sqrt = 0.40743 = Table 2 sigma urine salmeterol (mid) 41% (RSE 21%); applies to SID 2, 7
    expSdUrineSalHigh   <- 0.56745; label("Log-normal residual SD for urine salmeterol, high-variability study (fraction)")             # $SIGMA 7 = 0.322; sqrt = 0.56745 = Table 2 sigma urine salmeterol (high) 57% (RSE 26%); applies to SID 4
    expSdUrineOhsalMid  <- 0.37947; label("Log-normal residual SD for urine alpha-hydroxysalmeterol, medium-variability studies (fraction)") # $SIGMA 4 = 0.144; sqrt = 0.37947 = Table 2 sigma urine alpha-hydroxysalmeterol (mid) 38% (RSE 15%); applies to SID 1, 5, 7
    expSdUrineOhsalHigh <- 0.51381; label("Log-normal residual SD for urine alpha-hydroxysalmeterol, high-variability study (fraction)")     # $SIGMA BLOCK(2) element 2,2 = 0.264; sqrt = 0.51381 = Table 2 sigma urine alpha-hydroxysalmeterol (high) 51% (RSE 13%); applies to SID 2

    # Two residual-error correlations reported by the paper are not encoded:
    # the $SIGMA BLOCK(2) off-diagonals 0.0422 (plasma, Jessen 2021,
    # 0.0422 / sqrt(0.0461 x 0.145) = 52%, Table 2 reports 53%) and 0.127
    # (urine, Jessen 2021, 0.127 / sqrt(0.166 x 0.264) = 61%, Table 2 reports
    # 61%). nlmixr2 has no analogue of the NONMEM L2 record for correlated
    # residuals across endpoints. They affect the joint likelihood during
    # estimation but not the marginal distribution of either endpoint, so
    # simulation of the paper's per-endpoint percentiles is unaffected. See
    # the vignette Assumptions and deviations section.
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual parameters. Data S1 $PK. The athlete effect enters as
    #    a power of the binary indicator, exactly as the control stream
    #    writes it, so ATHLETE = 0 recovers the reference typical value.
    # -------------------------------------------------------------------
    vc        <- exp(lvc  + etalvc)
    q         <- exp(lq   + etalq)
    vp        <- exp(lvp  + etalvp)
    cl        <- exp(lcl  + etalcl) * e_athlete_cl^ATHLETE
    cl_ohsal  <- exp(lcl_ohsal)                                 # $OMEGA 5 is 0 FIX: no IIV
    kmet      <- exp(lkmet + etalkmet)
    kurine    <- exp(lkurine + etalkurine) * e_athlete_kurine^ATHLETE
    kurine_ohsal <- exp(lkurine_ohsal)                          # $OMEGA 8 is 0 FIX: no IIV

    # V3/F was constrained equal to V1/F for identifiability
    # (Data S1 $PK: V3 = V1; Methods 2.2.1 and Discussion).
    vc_ohsal <- vc

    # The urine-production IIV applies only to records whose concentration
    # was NOT USG-corrected. Data S1 $PK branches UR_PROD onto ETA(9) when
    # USG == 0 and onto ETA(10) when USG == 1, with ETA(10) fixed to 0.
    urprod <- exp(lurprod + etalurprod * (1 - USG_CORRECTED))

    # -------------------------------------------------------------------
    # 2. Micro-constants (Data S1 $PK). k10 and k30 are the residual
    #    elimination rates left after the explicitly modelled conversion
    #    and urinary-excretion routes are subtracted from the total
    #    apparent clearances, so total loss from each plasma compartment
    #    is exactly CL/V.
    # -------------------------------------------------------------------
    k12 <- q  / vc
    k21 <- q  / vp
    k10 <- cl / vc       - kmet - kurine        # $PK: K10 = CLP/V1 - K13 - K14
    k30 <- cl_ohsal / vc_ohsal - kurine_ohsal   # $PK: K30 = CLM/V3 - K35

    # -------------------------------------------------------------------
    # 3. ODE system (Data S1 $DES). Absorption is intravenous-like bolus
    #    into `central`: maximum plasma concentrations occur within 20 min
    #    of inhalation, so no absorption compartment was retained
    #    (Methods 2.2.1, Results 3.1.1). The dose therefore goes directly
    #    to `central` (NONMEM DEFDOSE on compartment 1).
    # -------------------------------------------------------------------
    d/dt(central)       <- k21 * peripheral1 - k12 * central - kmet * central - kurine * central - k10 * central
    d/dt(peripheral1)   <- k12 * central - k21 * peripheral1
    d/dt(central_ohsal) <- kmet * central - kurine_ohsal * central_ohsal - k30 * central_ohsal
    d/dt(urine)         <- kurine * central
    d/dt(urine_ohsal)   <- kurine_ohsal * central_ohsal
    d/dt(urine_vol)     <- urprod

    # -------------------------------------------------------------------
    # 4. Observations. Plasma concentrations use S1 = V1 and S3 = S1
    #    (Data S1 $PK). Urine concentrations divide the amount excreted
    #    since the last bladder voiding by the urine volume produced over
    #    the same period (Data S1 $ERROR: IPRED = A(4)/A(6)); a bladder
    #    voiding is an event-table record that resets `urine`,
    #    `urine_ohsal` and `urine_vol` to zero. For the two studies with
    #    recorded urine volumes the source model divided by the observed
    #    UVOL instead; that is a data property, not a structural
    #    difference, and both branches share the same ODEs.
    # -------------------------------------------------------------------
    Cc         <- central       / vc
    Cc_ohsal   <- central_ohsal / vc_ohsal
    urineSal   <- urine         / urine_vol
    urineOhsal <- urine_ohsal   / urine_vol

    # Study-stratified urine residual magnitudes, regrouped by the paper
    # into low / mid / high tiers (Results 3.1.2, Table 2). Written as a
    # sum of indicator products, following Friberg_2012_voriconazole.R.
    sdUrineSal <-
      expSdUrineSalLow  * ((STUDY_SALM == 1) + (STUDY_SALM == 5) + (STUDY_SALM == 6)) +
      expSdUrineSalMid  * ((STUDY_SALM == 2) + (STUDY_SALM == 7)) +
      expSdUrineSalHigh *  (STUDY_SALM == 4)
    sdUrineOhsal <-
      expSdUrineOhsalMid  * ((STUDY_SALM == 1) + (STUDY_SALM == 5) + (STUDY_SALM == 7)) +
      expSdUrineOhsalHigh *  (STUDY_SALM == 2)

    Cc         ~ lnorm(expSd)
    Cc_ohsal   ~ lnorm(expSd_ohsal)
    urineSal   ~ lnorm(sdUrineSal)
    urineOhsal ~ lnorm(sdUrineOhsal)
  })
}
