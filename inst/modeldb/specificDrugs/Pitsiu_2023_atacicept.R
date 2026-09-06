Pitsiu_2023_atacicept <- function() {
  description <- "Two-compartment quasi-steady-state (QSS) approximation of the target-mediated drug disposition (TMDD) model with first-order subcutaneous absorption for total atacicept (free plus BLyS/APRIL-bound) in healthy volunteers and patients with systemic lupus erythematosus (Pitsiu 2023, Table 1). Apparent clearance and central volume scale allometrically with body weight (exponents fixed at 0.75 and 1); the baseline total-target concentration Rmax scales as a power of baseline serum BLyS. Residual error is proportional and stratified by SLE status."
  reference <- "Pitsiu M, Yalkinoglu O, Farrell C, Girard P, Vazquez-Mateo C, Papasouliotis O. Population pharmacokinetics of atacicept in systemic lupus erythematosus: An analysis of three clinical trials. CPT Pharmacometrics Syst Pharmacol. 2023;12(8):1157-1169. doi:10.1002/psp4.12982"
  vignette <- "Pitsiu_2023_atacicept"

  # Unit system is inherited verbatim from the published NONMEM control stream
  # (supplement PSP4-12-1157-s001.txt, $PK): "S2=V2; assuming AMT in ug, DV in
  # ng/mL, V in L". Dose amounts are therefore MICROGRAMS -- a 150 mg dose is
  # amt = 150000. ug/L is identical to ng/mL, so dose/V lands directly in the
  # serum-concentration unit with no conversion factor.
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight. Allometric power effect on CL/F and Vc/F with reference weight 70 kg; both exponents were fixed (0.75 and 1.00), not estimated. Cohort median 65.0 kg, range 37.0-135 kg (Table S1).",
      source_name        = "WEIGHT"
    ),
    SBLYS = list(
      description        = "Baseline serum B lymphocyte stimulator (BLyS / BAFF / TNFSF13B) concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value. Power effect on the baseline total-target concentration Rmax,",
        "normalised to 2.56 ng/mL -- the median of the QUANTIFIABLE values only, not the overall",
        "cohort median (which is 1.79 ng/mL including BLQ subjects; Table S1). The control stream",
        "comments this reference as 'median excluding BLQ'. Values below the assay LLOQ of",
        "1.56 ng/mL are coded 0 in the source dataset and imputed to LLOQ/2 = 0.78 ng/mL; that",
        "imputation is reproduced inside model() so a 0 or negative input behaves as the paper",
        "intends. 235 of 540 subjects (44%) were BLQ (Table S1 footnotes d, e, f)."
      ),
      source_name        = "BLYS (imputed to BLYSB in $PK)"
    ),
    DIS_SLE = list(
      description        = "Systemic lupus erythematosus disease-state indicator (1 = SLE patient, 0 = healthy volunteer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer; the 37 phase I subjects)",
      notes              = paste(
        "Time-fixed per subject. Affects the RESIDUAL ERROR ONLY -- it selects between the two",
        "proportional error terms in $ERROR (IF (SLE.EQ.1) W=SQRT(IPRED**2*PROPSLE**2)). SLE was",
        "also tested as a structural covariate and on relative bioavailability, and was NOT",
        "retained: the paper reports no significant PK difference between healthy volunteers and",
        "patients with SLE. 503 of 540 subjects (93.1%) were SLE patients."
      ),
      source_name        = "SLE"
    )
  )

  # Covariates the paper SCREENED but did not retain in the final model. Kept
  # here for provenance; checkModelConventions() treats this list as
  # documentation only and does not require these to appear in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "The previously published model carried an age effect on Ka. Adding IIV on Ka in the present analysis made that effect non-significant, so it was dropped (Results, 'Final population PK model'). Median 37 y, range 16-75 y (Table S1)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the covariate search, not retained. Reported as raw mL/min in Table S1 (median 110, range 39.0-270), i.e. NOT BSA-normalised as the CRCL register entry's canonical definition assumes."
    ),
    SAPRIL = list(
      description = "Baseline serum a proliferation-inducing ligand (APRIL / TNFSF13) concentration",
      units       = "pg/mL",
      type        = "continuous",
      notes       = paste(
        "Atacicept's second soluble target. Screened in the covariate search alongside BLyS but",
        "NOT retained. Median 2011 pg/mL, range 0.00-5699 (Table S1); note the unit is pg/mL,",
        "unlike BLyS in ng/mL. DOCUMENTATION-ONLY KEY: `SAPRIL` is deliberately NOT registered in",
        "inst/references/covariate-columns.md, because registering a canonical no model uses is",
        "what a prior register audit had to undo. A future extraction that RETAINS an APRIL effect",
        "should ratify the name at that point."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened ('gender'), not retained. 89.6% female (Table S1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened, not retained: 'no significant differences in PKs were detected ... between different racial groups'. 20.0% Asian (Table S1)."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened, not retained (see RACE_ASIAN). 4.1% African/African American (Table S1)."
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the published control stream ($MODEL /
  # $DES, supplement s001.txt) and the Figure S1 schematic.
  compartmentData <- list(
    depot = list(
      analyte = "atacicept", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "atacicept (total: unbound plus BLyS/APRIL-bound)", units = "ug",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "atacicept (unbound)", units = "ug",
      specimen = "serum", verified = TRUE
    ),
    # Unlike the three drug states, this state holds a CONCENTRATION, not an
    # amount: the control stream sets A_0(4)=RMAX with RMAX in ng/mL and
    # KSYN = RMAX*KDEG, and $DES divides no volume into it.
    total_target = list(
      analyte = "BLyS + APRIL (total target: unbound plus atacicept-bound)", units = "ng/mL",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 540,
    n_studies      = 3,
    n_observations = 3640,
    age_range      = "16-75 years",
    age_median     = "37 years",
    weight_range   = "37.0-135 kg",
    weight_median  = "65.0 kg",
    sex_female_pct = 89.6,
    race_ethnicity = c(White = 69.1, `African/African American` = 4.1, Asian = 20.0, Other = 6.9),
    disease_state  = "503 patients (93.1%) with moderate-to-severe or active autoantibody-positive systemic lupus erythematosus; 37 healthy volunteers (6.9%).",
    dose_range     = "25, 75, or 150 mg subcutaneous: single dose (phase I), or bi-weekly for 4 weeks then once weekly to week 52 (APRIL-SLE), or once weekly for 24 weeks (ADDRESS II).",
    regions        = "Multinational; the phase I study enrolled Japanese and White healthy volunteers, ADDRESS II 46.8% Hispanic/Latino.",
    renal_function = "CrCL median 110 mL/min, range 39.0-270 mL/min.",
    notes          = "Pooled analysis of EMR700461-022 (phase I, n = 37, 533 observations), APRIL-SLE / NCT00624338 (phase II, n = 298, 1728 observations), and ADDRESS II / NCT01972568 (phase IIb, n = 205, 1379 observations). Baseline demographics from supplementary Table S1. Total atacicept was measured by acid-dissociation ELISA (LLOQ 100 ng/mL); post-dose BLQ records were under 2% of the data set and were excluded from the fit."
  )

  # The paper's two proportional residual SDs are cohort-stratified rather than
  # per-output, so they cannot use the bare canonical `propSd` name. Declared
  # here following the Shoji_2011_pregabalin / Ozdin_2025_dexamethasone
  # precedent; the effective SD is assembled inside model().
  paper_specific_residual_sds <- c("propSdNonSle", "propSdSle")

  ini({
    # --- Structural PK, typical 70 kg subject (Table 1, "NONMEM estimates") ---
    lcl    <- log(0.324);    label("Apparent nonspecific clearance for a 70 kg subject (CL/F, L/h)")           # Table 1 CL/F 0.324 (95% CI 0.298-0.350, RSE 4.10%)
    lvc    <- log(36.3);     label("Apparent central volume of distribution for a 70 kg subject (Vc/F, L)")    # Table 1 Vc/F 36.3 (95% CI 31.9-40.7, RSE 6.14%)
    lq     <- log(0.149);    label("Apparent inter-compartmental clearance (Q/F, L/h)")                        # Table 1 Q/F 0.149 (95% CI 0.114-0.184, RSE 11.9%)
    lvp    <- log(38.5);     label("Apparent peripheral volume of distribution (Vp/F, L)")                     # Table 1 Vp/F 38.5 (95% CI 31.0-46.0, RSE 9.90%)
    lka    <- log(0.0705);   label("First-order subcutaneous absorption rate constant (Ka, 1/h)")              # Table 1 Ka 0.0705 (95% CI 0.0595-0.0815, RSE 7.94%)

    # --- Target binding and turnover (QSS TMDD) ---
    lkss   <- log(19.9);     label("Quasi-steady-state constant for atacicept-target binding (Kss, ng/mL)")    # Table 1 Kss 19.9 (95% CI 14.4-25.4, RSE 14.1%)
    lkint  <- log(0.000618); label("Drug-target complex elimination rate constant (Kint, 1/h)")                # Table 1 Kint 0.000618 (95% CI 0.000572-0.000664, RSE 3.83%)
    lkdeg  <- log(0.00362);  label("Free-target elimination rate constant (Kdeg, 1/h)")                        # Table 1 Kdeg 0.00362 (95% CI 0.00307-0.00417, RSE 7.82%)
    lrbase <- log(715);      label("Baseline total-target concentration at the reference BLyS of 2.56 ng/mL (Rmax, ng/mL)")  # Table 1 Rmax 715 (95% CI 613-817, RSE 7.27%)

    # --- Covariate effects ---
    # Allometric exponents held constant during estimation ("0.75 fixed" /
    # "1.00 fixed" in Table 1, with no CI and no RSE reported).
    e_wt_cl       <- fixed(0.75); label("Allometric exponent on (WT/70) for CL/F (unitless)")   # Table 1 "Weight on CL" 0.75 fixed
    e_wt_vc       <- fixed(1.00); label("Allometric exponent on (WT/70) for Vc/F (unitless)")   # Table 1 "Weight on Vc" 1.00 fixed
    e_sblys_rbase <- 0.176;       label("Power exponent on (SBLYS/2.56) for Rmax (unitless)")   # Table 1 "BLyS on Rmax" 0.176 (95% CI 0.120-0.232, RSE 16.2%)

    # --- IIV: diagonal Omega, exponential on each parameter ---
    # Table 1 reports the VARIANCES omega^2 directly; the CV% column is
    # 100*sqrt(omega^2), which reproduces each printed CV% exactly
    # (e.g. 100*sqrt(0.233) = 48.3%), confirming the variance scale.
    etalcl    ~ 0.233  # Table 1 omega^2 CL   = 0.233 (CV 48.3%, shrinkage 23.0%)
    etalvc    ~ 0.284  # Table 1 omega^2 Vc   = 0.284 (CV 53.3%, shrinkage 48.8%)
    etalrbase ~ 0.102  # Table 1 omega^2 Rmax = 0.102 (CV 31.9%, shrinkage 32.8%)
    etalvp    ~ 0.532  # Table 1 omega^2 Vp   = 0.532 (CV 72.9%, shrinkage 43.3%)
    etalka    ~ 0.182  # Table 1 omega^2 Ka   = 0.182 (CV 42.7%, shrinkage 71.7%)

    # --- Residual error: two proportional terms, selected by SLE status ---
    # $ERROR sets ADD=0.0, so the model is purely proportional; W is the SD
    # fraction directly (the CV% column of Table 1 is 100 * the point estimate).
    propSdNonSle <- 0.188; label("Proportional residual SD, healthy volunteers (fraction)")        # Table 1 "Proportional error no SLE" 0.188 (95% CI 0.181-0.195, CV 18.8%)
    propSdSle    <- 0.251; label("Proportional residual SD, patients with SLE (fraction)")         # Table 1 "Proportional error SLE" 0.251 (95% CI 0.247-0.255, CV 25.1%)
  })

  model({
    # ---- 1. Derived covariate terms -------------------------------------
    # BLQ imputation, verbatim from $PK: BLYS is coded 0 when the baseline
    # sample was below the 1.56 ng/mL assay LLOQ, and is replaced by
    # LLOQ/2 = 0.78 ng/mL before the power term is formed.
    #   BLYSB = BLYS
    #   IF(BLYS.EQ.0) BLYSB = 0.78
    # Reproducing this matters: simulating the paper's "BLOQ" row of Table 3
    # with 1.56 instead of 0.78 overpredicts steady-state Cmax by up to 9%.
    sblysb <- SBLYS * (SBLYS > 0) + 0.78 * (SBLYS <= 0)

    # ---- 2. Individual parameters ---------------------------------------
    cl    <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc    <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q     <- exp(lq)                                    # no IIV on Q in the final model
    vp    <- exp(lvp + etalvp)
    ka    <- exp(lka + etalka)
    kss   <- exp(lkss)
    kint  <- exp(lkint)
    kdeg  <- exp(lkdeg)
    rbase <- exp(lrbase + etalrbase) * (sblysb / 2.56)^e_sblys_rbase

    # ---- 3. Micro-constants ---------------------------------------------
    kel <- cl / vc
    k12 <- q / vc     # K23 in the control stream (central -> peripheral)
    k21 <- q / vp     # K32 in the control stream (peripheral -> central)

    # Zero-order target synthesis, set so the target sits at its baseline in
    # the absence of drug: KSYN = RMAX * KDEG ($PK).
    ksyn <- rbase * kdeg

    # ---- 4. QSS algebra (Gibiansky 2008; Figure S1) ----------------------
    # `central` carries TOTAL atacicept amount (unbound + target-bound) and
    # `total_target` carries TOTAL target concentration Rtot. The unbound drug
    # concentration C is the positive root of the QSS binding quadratic:
    #   C = 0.5 * [(Ctot - Rtot - Kss) + sqrt((Ctot - Rtot - Kss)^2 + 4*Kss*Ctot)]
    ctot  <- central / vc
    rtot  <- total_target
    cfree <- 0.5 * ((ctot - rtot - kss) +
                      sqrt((ctot - rtot - kss)^2 + 4 * kss * ctot))

    # ---- 5. ODE system ($DES; identical to Figure S1) --------------------
    # Distribution and nonspecific elimination act on UNBOUND drug (hence the
    # cfree * vc products); the target-mediated term removes complex at Kint.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot -
                           (kel + k12) * cfree * vc -
                           kint * cfree * vc * rtot / (kss + cfree) +
                           k21 * peripheral1
    d/dt(peripheral1) <- k12 * cfree * vc - k21 * peripheral1
    d/dt(total_target) <- ksyn - kdeg * rtot -
                            (kint - kdeg) * rtot * cfree / (kss + cfree)

    # Target starts at its drug-free baseline (A_0(4) = RMAX in $PK).
    total_target(0) <- rbase

    # ---- 6. Observation and error ---------------------------------------
    # The assay measures TOTAL atacicept after an acid-dissociation step, so
    # the observation is Ctot, not the unbound concentration (IPRED=A(2)/V2).
    Cc <- ctot

    # $ERROR selects one of the two proportional SDs on the SLE flag.
    propSd <- propSdNonSle * (1 - DIS_SLE) + propSdSle * DIS_SLE
    Cc ~ prop(propSd)
  })
}
