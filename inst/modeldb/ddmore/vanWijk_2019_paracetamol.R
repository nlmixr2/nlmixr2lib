vanWijk_2019_paracetamol <- function() {
  description <- "PRECLINICAL (zebrafish): two-compartment paracetamol PK model fit to zebrafish (Danio rerio) larvae continuously exposed to a 1 mM paracetamol bath at 3, 4, or 5 days post-fertilization (van Wijk 2019, DDMODEL00000294). The medium reservoir (compartment 1) is held at constant amount, so K12 acts as a zero-order absorption rate from the bath into the larva; elimination from the larva (compartment 2) is first-order with rate K25. Larval age in dpf enters as a step factor on K12 (~2.06x at >= 4 dpf vs 3 dpf) and a per-day power factor on K25 (+17.4% per day post-fertilization), consistent with maturation of paracetamol absorption and elimination capacity across the 3-5 dpf window."
  reference <- paste(
    "van Wijk RC, Krekels EHJ, Hankemeier T, Spaink HP, van der Graaf PH (2019).",
    "Impact of post-hatching maturation on the pharmacokinetics of paracetamol in zebrafish larvae.",
    "Sci Rep 9(1):2149.",
    "doi:10.1038/s41598-019-38530-w.",
    "DDMORE Foundation Model Repository: DDMODEL00000294.",
    sep = " "
  )
  vignette <- "vanWijk_2019_paracetamol"
  units <- list(time = "minute", dosing = "pmol", concentration = "pmol/larva")

  ddmore_id    <- "DDMODEL00000294"
  replicate_of <- NULL

  covariateData <- list(
    AGE_DPF = list(
      description        = "Zebrafish-larval age in days post-fertilization at the start of paracetamol exposure (integer 3, 4, or 5 in van Wijk 2019). Time-fixed per subject under the destructive-sampling design (each larva is harvested at exactly one observation time).",
      units              = "days post-fertilization (dpf)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference age 3 dpf (the youngest cohort). Step-form covariate effect on K12 absorption",
        "(.mod $PK line: IF(AGE.GT.3) TVK12 = THETA(2)*(1+THETA(3))) -- so K12 doubles (~2.06x) at",
        ">= 4 dpf vs 3 dpf. Per-day power-form effect on K25 elimination",
        "(.mod $PK line: K25 = TVK25 * ((1 + THETA(4)) ** (AGE - 3))) -- so K25 increases ~17.4% per day post-fertilization.",
        "Source data column AGE (integer 3..5) renamed to canonical AGE_DPF on input to avoid",
        "collision with the human-PK canonical AGE (years)."
      ),
      source_name        = "AGE"
    )
  )

  population <- list(
    n_subjects     = 242L,
    n_studies      = 1L,
    age_range      = "3-5 days post-fertilization (dpf) at the start of the 1 mM paracetamol bath exposure (van Wijk 2019 zebrafish-larvae study).",
    weight_range   = "Not extractable from DDMORE bundle; the model parameterises elimination as a per-larva first-order rate constant rather than via an explicit larval volume.",
    sex_female_pct = "Not applicable (zebrafish larvae 3-5 dpf are sexually undifferentiated; sex is not a covariate in the model).",
    disease_state  = "Healthy zebrafish (Danio rerio) larvae at 3, 4, or 5 days post-fertilization.",
    dose_range     = paste(
      "Continuous environmental exposure to 1 mM paracetamol in the surrounding medium (E3 zebrafish water).",
      "The .mod represents this as a single dose AMT = 1 (arbitrary unit) into a depot compartment whose",
      "DADT is held at zero throughout the simulation, so the depot amount stays constant at AMT and the",
      "absorption rate K12 (pmol/min) sets the per-minute pmol load into each larva."
    ),
    regions        = "In vitro / aquatic study; no clinical regions.",
    species        = "Zebrafish (Danio rerio) larvae 3-5 dpf -- preclinical entry.",
    notes          = paste(
      "n_subjects (242) and n_obs (177) come from the Output_real_Paracetamol_Zebrafish_345dpf.lst",
      "run summary lines 'TOT. NO. OF INDIVIDUALS:' and 'TOT. NO. OF OBS RECS:' respectively.",
      "Destructive-sampling design: each larva contributes one DV observation, then is sacrificed for",
      "the assay (which is why the .mod $OMEGA on K25 was held at FIX 0 -- per-larva IIV is",
      "indistinguishable from residual error when each animal is sampled exactly once).",
      "The van Wijk 2019 PDF is not on disk under the literature tree, so the original n / cohort",
      "structure / inclusion criteria could not be cross-checked against the DDMORE bundle's run summary."
    )
  )

  ini({
    # Final parameter estimates from
    # /home/bill/github/mab_human_consensus/literature/from_people/ddmore/ddmore_scraping/294/Output_real_Paracetamol_Zebrafish_345dpf.lst
    # FINAL PARAMETER ESTIMATE block (.lst lines 313-352), captured after
    # `MINIMIZATION SUCCESSFUL` (.lst line 260, OBJV 466.583). NONMEM THETAs
    # are log-back-transformed values, so each `lX <- log(value)` wraps the
    # .lst value to keep the internal scale log. NONMEM SIGMA entries are
    # variances on the relevant linear scale; the nlmixr2 `prop()` / `add()`
    # SDs are the sqrt of those variances and the conversion is documented
    # per parameter.

    # Naming deviation: the model is parameterised in micro-constants
    # (K12 = pmol/min absorption rate from the constant-amount medium
    # depot; K25 = 1/min first-order elimination rate constant from the
    # larva). Neither maps cleanly onto the canonical (lka, lcl, lvc) PK
    # parameter set because there is no explicit larval volume in the
    # source -- elimination is parameterised on amount, not on
    # concentration. We keep the source's K12 / K25 names (with the `l`
    # log-transform prefix) so the source trace is unambiguous; the
    # vignette's Assumptions and deviations section records this as a
    # checkModelConventions() deviation.
    lk25 <- log(0.0193) ; label("First-order elimination rate constant K25 at AGE_DPF = 3 dpf (1/min)")          # DDMODEL00000294 .lst FINAL THETA(1) = 1.93E-02
    lk12 <- log(0.289)  ; label("Zero-order absorption rate K12 at AGE_DPF = 3 dpf (pmol/min per unit depot)")    # DDMODEL00000294 .lst FINAL THETA(2) = 2.89E-01

    e_age_dpf_k12 <- 1.06  ; label("Fractional change in K12 at AGE_DPF > 3 dpf (unitless step; K12 at >= 4 dpf is K12_3 * (1 + e_age_dpf_k12))")  # DDMODEL00000294 .lst FINAL THETA(3) = 1.06E+00
    e_age_dpf_k25 <- 0.174 ; label("Per-day fractional increase in K25 (unitless; K25(age) = K25_3 * (1 + e_age_dpf_k25)^(AGE_DPF - 3))")          # DDMODEL00000294 .lst FINAL THETA(4) = 1.74E-01

    # Inter-individual variability. The .mod declared $OMEGA 0 FIX
    # (.mod line 47, .lst INITIAL OMEGA 0.00E+00 with 'OMEGA CONSTRAINED
    # TO BE THIS INITIAL ESTIMATE'); per the .mod's own comment the IIV
    # on K25 is 'undistinguishable from residual variability due to
    # destructive sampling'. No eta is declared on either rate constant.

    # Residual error. The .mod $ERROR uses a combined linear-scale model:
    # `Y = IPRED * (1 + EPS(1)) + EPS(2)` (.mod line 40). NONMEM SIGMAs
    # are variances; nlmixr2 `prop()` / `add()` SDs are sqrt of those
    # variances. Note the unusual additive-error units: DV in van Wijk
    # 2019 is amount-per-larva (pmol/larva), so `addSd` is pmol/larva,
    # not the conventional concentration unit.
    propSd <- sqrt(0.10906)  ; label("Proportional residual error on per-larva paracetamol amount Cc (fraction)")            # DDMODEL00000294 .lst FINAL SIGMA EPS1 (variance) = 1.09E-01 -> SD 3.30E-01
    addSd  <- sqrt(0.0084383); label("Additive residual error on per-larva paracetamol amount Cc (pmol/larva)")              # DDMODEL00000294 .lst FINAL SIGMA EPS2 (variance) = 8.44E-03 -> SD 9.18E-02
  })
  model({
    # Typical-value PK parameters. The .mod $PK block carries:
    #   TVK12 = THETA(2)
    #   IF (AGE > 3) TVK12 = THETA(2) * (1 + THETA(3))
    #   TVK25 = THETA(1) * EXP(ETA(1))                       (ETA(1) is FIX 0 here)
    #   K12   = TVK12
    #   K25   = TVK25 * ((1 + THETA(4)) ** (AGE - 3))
    # which translates verbatim below. The (AGE_DPF > 3) indicator is
    # used as a (0/1) multiplier on the e_age_dpf_k12 step.
    age_step  <- (AGE_DPF > 3)
    k12 <- exp(lk12) * (1 + e_age_dpf_k12 * age_step)
    k25 <- exp(lk25) * (1 + e_age_dpf_k25)^(AGE_DPF - 3)

    # ODE system. The medium reservoir is held at constant amount via
    # `d/dt(depot) <- 0`, exactly reproducing the .mod's `DADT(1) = 0
    # ;constant infusion` (.mod line 35). Doses to the depot at t = 0
    # therefore set A(depot) once and that mass persists, so the K12
    # term acts as a zero-order influx of pmol/min into the larva.
    # Dosing in van Wijk 2019 uses AMT = 1 as the dose amount on every
    # dosing record; users simulating from this model should preserve
    # AMT = 1 unless they intend to rescale the bath exposure.
    d/dt(depot)   <- 0
    d/dt(central) <- k12 * depot - k25 * central

    # Observation. In van Wijk 2019 DV is reported as pmol/larva
    # (amount per larva), with no explicit larval volume in the model
    # (the .mod's `; V = total larval volume` comment notes that V is
    # held fixed and absorbed into the rate constants). We carry this
    # forward by setting Cc = central directly; Cc therefore reads as
    # pmol/larva, not the conventional mass/volume concentration.
    Cc <- central

    Cc ~ prop(propSd) + add(addSd)
  })
}
