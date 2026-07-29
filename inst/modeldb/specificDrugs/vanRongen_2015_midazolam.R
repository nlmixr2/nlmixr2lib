vanRongen_2015_midazolam <- function() {
  description <- paste(
    "Joint parent-and-sequential-metabolites population PK model for",
    "intravenous midazolam, its primary CYP3A oxidative metabolite",
    "1'-hydroxymidazolam (1-OH-midazolam), and the downstream phase-II",
    "1'-hydroxymidazolam glucuronide (1-OH-midazolam glucuronide) in 19",
    "overweight and obese adolescents (12.5-18.9 years, body weight",
    "62-149.8 kg) undergoing surgery (van Rongen 2015). Two-compartment",
    "disposition for midazolam (central + peripheral) routes the entire",
    "elimination clearance CL1 to 1-OH-midazolam formation. 1-OH-midazolam",
    "is described by a one-compartment model with apparent volume of",
    "distribution fixed at 0.9 times the midazolam central volume",
    "(Mandema 1992); the entire 1-OH-midazolam clearance CL3 is routed to",
    "1-OH-midazolam glucuronide formation. 1-OH-midazolam glucuronide is",
    "described by a two-compartment model with renal elimination",
    "clearance CL4. Total body weight (TBW) enters the peripheral volume",
    "of distribution of midazolam as a power function with reference 104.7",
    "kg (cohort median) and estimated exponent X = 1.68; no other",
    "covariate effect was retained. Concentrations are modeled in umol/L",
    "throughout (paper Methods), so dosing is in umol; dose_umol = dose_mg",
    "* 1000 / MW_midazolam where MW_midazolam = 325.77 g/mol."
  )
  reference <- paste(
    "van Rongen A, Vaughns JD, Moorthy GS, Barrett JS, Knibbe CAJ,",
    "van den Anker JN (2015). Population pharmacokinetics of midazolam",
    "and its metabolites in overweight and obese adolescents.",
    "Br J Clin Pharmacol 80(5):1185-1196.",
    "doi:10.1111/bcp.12693.",
    sep = " "
  )
  vignette <- "vanRongen_2015_midazolam"
  units <- list(time = "min", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-function scaling on the midazolam peripheral volume of",
        "distribution V_mdz_peripheral with reference TBW = 104.7 kg",
        "(cohort median) and estimated exponent X = 1.68 (paper Table 2",
        "final model: V_mdz_peripheral = V_{104.7 kg} * (TBW / 104.7)^X).",
        "No other PK parameter retained TBW as a covariate in the final",
        "model. Time-fixed at baseline."
      ),
      source_name        = "TBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 19L,
    n_studies      = 1L,
    age_range      = "12.5-18.9 years",
    age_median     = "15.9 years (mean +/- 1.6 SD)",
    weight_range   = "62-149.8 kg",
    weight_median  = "102.7 kg (mean +/- 24.9 SD); 104.7 kg (model reference)",
    sex_female_pct = 68.4,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Overweight and obese adolescents (BMI for age >=85th percentile",
      "overweight; >=95th percentile obese) undergoing general surgery",
      "(orthopaedics, tonsillectomy, bariatric surgery) with ASA",
      "physical status I, II, or III. Three patients were overweight and",
      "16 obese; BMI 24.8-55 kg/m^2, BMI z-score 1.5-2.7. Patients with",
      "prior benzodiazepine exposure, liver or renal disease, or CYP3A-",
      "modulating co-medication were excluded."
    ),
    dose_range     = paste(
      "Single intravenous bolus dose of 2 mg (n = 16) or 3 mg (n = 3)",
      "midazolam administered a few minutes before transfer to the",
      "operating room. Blood samples at t = 0, (5), 15, 30 min and 1, 2,",
      "4, 6, occasionally 8 h post-dose; 129 midazolam, 118",
      "1-OH-midazolam, and 128 1-OH-midazolam glucuronide samples."
    ),
    regions        = "United States (Children's National Medical Center, Washington DC)",
    notes          = paste(
      "Demographics from Table 1 of van Rongen 2015. Cohort median TBW",
      "104.7 kg is the reference used in the power-function TBW effect",
      "on V_mdz_peripheral (Table 2 final model). 1-OH-midazolam apparent",
      "volume of distribution was fixed at 0.9 x V_mdz_central per",
      "Mandema 1992 (paper Methods, citing reference [20])."
    )
  )

  ini({
    # Final population PK parameters from Table 2 of van Rongen 2015 (page
    # 1188), "Final model (RSE%)" column. Bootstrap 95% confidence intervals
    # in the right-most column confirm the point estimates. The paper
    # modeled concentrations on the umol/L molar scale (Methods page 1187:
    # "Concentrations were expressed as umol/l using the molecular weights
    # of midazolam, 1-OH-midazolam and 1-OH-midazolam glucuronide (325.77,
    # 341.77 and 517.9, respectively)"), and reported three independent
    # proportional residual error models, one per analyte.

    # Structural midazolam disposition. CL1 in the paper is the
    # metabolic clearance from midazolam to 1-OH-midazolam and -- because
    # the model assumes 1-OH-midazolam formation is the sole elimination
    # route for midazolam -- equals the total midazolam clearance. The
    # peripheral volume V_mdz_peripheral varies with TBW; the reported
    # value 154 L is the typical-individual estimate at the reference
    # median TBW of 104.7 kg.
    lcl    <- log(0.66)    ; label("Midazolam total clearance (= CL1, L/min)")                     # Table 2 final CL1 = 0.66 L/min (RSE 8.3%)
    lvc    <- log(39.8)    ; label("Midazolam central volume of distribution V_mdz,central (L)")   # Table 2 final V_mdz,central = 39.8 L (RSE 8.3%)
    lvp    <- log(154)     ; label("Midazolam peripheral volume at TBW = 104.7 kg (L)")            # Table 2 final V_104.7kg = 154 L (RSE 11.2%) -- typical V_mdz,peripheral at the cohort median TBW; the TBW power covariate uses this as reference
    lq     <- log(1.18)    ; label("Midazolam inter-compartmental clearance Q (L/min)")            # Table 2 final Q = 1.18 L/min (RSE 15.6%)

    # 1-OH-midazolam disposition. CL3 in the paper is the metabolic
    # clearance from 1-OH-midazolam to 1-OH-midazolam glucuronide and
    # equals the total 1-OH-midazolam clearance because 1-OH-midazolam
    # glucuronide formation is the sole elimination route for 1-OH-
    # midazolam in this model. The apparent volume of distribution of
    # 1-OH-midazolam is fixed at 0.9 x V_mdz,central per Mandema 1992,
    # so vc_1ohm is computed inline in model() and there is no
    # corresponding lvc_1ohm in ini().
    lcl_1ohm  <- log(1.85) ; label("1-OH-midazolam total clearance (= CL3, L/min)")                # Table 2 final CL3 = 1.85 L/min (RSE 9.3%)

    # 1-OH-midazolam glucuronide disposition. CL4 in the paper is the
    # renal elimination clearance of the glucuronide conjugate.
    lvc_1ohmg <- log(4.05) ; label("1-OH-midazolam glucuronide central volume V_1OHgluc,central (L)")    # Table 2 final V_1OHgluc,central = 4.05 L (RSE 17.5%)
    lvp_1ohmg <- log(15.9) ; label("1-OH-midazolam glucuronide peripheral volume V_1OHgluc,peripheral (L)") # Table 2 final V_1OHgluc,peripheral = 15.9 L (RSE 9.5%)
    lq_1ohmg  <- log(0.49) ; label("1-OH-midazolam glucuronide inter-compartmental clearance Q2 (L/min)")   # Table 2 final Q2 = 0.49 L/min (RSE 23.9%)
    lcl_1ohmg <- log(0.42) ; label("1-OH-midazolam glucuronide elimination clearance (= CL4, L/min)")       # Table 2 final CL4 = 0.42 L/min (RSE 6.4%)

    # Covariate effect: TBW on V_mdz,peripheral via a power function.
    # V_mdz,peripheral = V_{104.7 kg} * (TBW / 104.7)^X. Estimated with
    # RSE 12.1%, bootstrap 95% CI 0.9-2.63 (Table 2).
    e_wt_vp <- 1.68        ; label("Power exponent of TBW on V_mdz,peripheral (unitless)")         # Table 2 final X = 1.68 (RSE 12.1%)

    # Inter-individual variability. Table 2 reports IIV on the log-normal
    # CV% scale; the internal variance is omega^2 = log(1 + CV^2). The
    # paper notes "V_mdz,central = V_1-OH" in the IIV section of Table 2,
    # meaning the IIV on the midazolam central volume of distribution is
    # shared with that on the 1-OH-midazolam apparent volume of
    # distribution. Because vc_1ohm = 0.9 * vc inherits the etalvc random
    # effect, a single etalvc covers both volumes and no separate
    # etalvc_1ohm is needed. No IIV is reported for the 1-OH-midazolam
    # glucuronide central / peripheral volumes, Q2, or CL4.
    etalcl     ~ log(1 + 0.237^2)   # Table 2 final IIV CL1 = 23.7% CV (RSE 25)
    etalvc     ~ log(1 + 0.305^2)   # Table 2 final IIV V_mdz,central = V_1-OH = 30.5% CV (RSE 14.4)
    etalvp     ~ log(1 + 0.302^2)   # Table 2 final IIV V_mdz,peripheral = 30.2% CV (RSE 32.6)
    etalq      ~ log(1 + 0.395^2)   # Table 2 final IIV Q = 39.5% CV (RSE 18.8)
    etalcl_1ohm ~ log(1 + 0.267^2)  # Table 2 final IIV CL3 = 26.7% CV (RSE 20)

    # Residual variability. Table 2 reports three independent proportional
    # residual error models, one per analyte (Methods page 1187: "Residual
    # variability was tested using proportional, additive or combined
    # proportional and additive error models for midazolam and metabolites
    # ... best described by three proportional error models for the
    # midazolam, 1-OH-midazolam and 1-OH-midazolam glucuronide
    # concentrations"). The reported CV% maps directly to nlmixr2's propSd
    # in linear (umol/L) space.
    propSd       <- 0.265  ; label("Midazolam proportional residual SD (fraction)")                 # Table 2 final proportional error MDZ = 26.5% (RSE 14)
    propSd_1ohm  <- 0.260  ; label("1-OH-midazolam proportional residual SD (fraction)")            # Table 2 final proportional error 1-OH = 26.0% (RSE 13)
    propSd_1ohmg <- 0.234  ; label("1-OH-midazolam glucuronide proportional residual SD (fraction)") # Table 2 final proportional error 1-OH-gluc = 23.4% (RSE 7)
  })

  model({
    # Individual structural parameters. Only the midazolam peripheral
    # volume of distribution carries a covariate (TBW power function,
    # reference 104.7 kg, exponent e_wt_vp). The 1-OH-midazolam apparent
    # volume of distribution is structurally tied to the midazolam
    # central volume of distribution via the fixed 0.9 ratio (paper
    # Methods page 1187 citing Mandema 1992); the etalvc random effect
    # therefore propagates through vc_1ohm with identical magnitude, as
    # Table 2 reports ("V_mdz,central = V_1-OH" under IIV).
    cl        <- exp(lcl       + etalcl)
    vc        <- exp(lvc       + etalvc)
    vp        <- exp(lvp       + etalvp) * (WT / 104.7)^e_wt_vp
    q         <- exp(lq        + etalq)
    vc_1ohm   <- 0.9 * vc
    cl_1ohm   <- exp(lcl_1ohm  + etalcl_1ohm)
    vc_1ohmg  <- exp(lvc_1ohmg)
    vp_1ohmg  <- exp(lvp_1ohmg)
    q_1ohmg   <- exp(lq_1ohmg)
    cl_1ohmg  <- exp(lcl_1ohmg)

    # Micro-constants. Each metabolic clearance divided by the upstream
    # volume gives the first-order rate constant of parent disappearance
    # / daughter formation; metabolite elimination rates are computed in
    # the same way. Mole-for-mole conservation across each metabolic step
    # is preserved because the model is written in molar amounts (umol)
    # and the metabolic flux from the parent compartment is identical to
    # the formation flux into the daughter compartment.
    k_form_1ohm  <- cl       / vc        # MDZ central -> 1-OH-midazolam (CL1 / V_mdz,central)
    k_mdz_cp     <- q        / vc        # MDZ central -> peripheral1
    k_mdz_pc     <- q        / vp        # MDZ peripheral1 -> central
    k_form_1ohmg <- cl_1ohm  / vc_1ohm   # 1-OH-midazolam -> 1-OH-midazolam glucuronide (CL3 / V_1-OH-midazolam)
    k_el_1ohmg   <- cl_1ohmg / vc_1ohmg  # 1-OH-midazolam glucuronide central -> elimination
    k_1ohmg_cp   <- q_1ohmg  / vc_1ohmg  # 1-OH-midazolam glucuronide central -> peripheral1
    k_1ohmg_pc   <- q_1ohmg  / vp_1ohmg  # 1-OH-midazolam glucuronide peripheral1 -> central

    # ODE system matching Figure 1 of van Rongen 2015 (paper page 1188).
    # All amounts are in umol; concentrations are amount divided by the
    # corresponding L volume, yielding umol/L (matching the units the
    # paper reports throughout). Mole conservation across the metabolic
    # sequence is automatic: 1 umol of midazolam metabolized produces 1
    # umol of 1-OH-midazolam, which on further metabolism produces 1 umol
    # of 1-OH-midazolam glucuronide.
    d/dt(central)            <- -k_form_1ohm  * central          -
                                 k_mdz_cp     * central          +
                                 k_mdz_pc     * peripheral1
    d/dt(peripheral1)        <-  k_mdz_cp     * central          -
                                 k_mdz_pc     * peripheral1
    d/dt(central_1ohm)       <-  k_form_1ohm  * central          -
                                 k_form_1ohmg * central_1ohm
    d/dt(central_1ohmg)      <-  k_form_1ohmg * central_1ohm     -
                                 k_el_1ohmg   * central_1ohmg    -
                                 k_1ohmg_cp   * central_1ohmg    +
                                 k_1ohmg_pc   * peripheral1_1ohmg
    d/dt(peripheral1_1ohmg)  <-  k_1ohmg_cp   * central_1ohmg    -
                                 k_1ohmg_pc   * peripheral1_1ohmg

    # Plasma concentrations in umol/L. Users dosing in mg should convert
    # to umol via dose_umol = dose_mg * 1000 / MW_midazolam, where
    # MW_midazolam = 325.77 g/mol (Methods page 1187). The metabolite
    # outputs Cc_1ohm and Cc_1ohmg are directly in umol/L of the
    # respective metabolite species and may be converted to ng/mL via
    # the metabolite molecular weights 341.77 g/mol (1-OH-midazolam) and
    # 517.9 g/mol (1-OH-midazolam glucuronide).
    Cc       <- central       / vc
    Cc_1ohm  <- central_1ohm  / vc_1ohm
    Cc_1ohmg <- central_1ohmg / vc_1ohmg

    Cc       ~ prop(propSd)
    Cc_1ohm  ~ prop(propSd_1ohm)
    Cc_1ohmg ~ prop(propSd_1ohmg)
  })
}
