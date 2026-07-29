Brill_2014_midazolam <- function() {
  description <- paste(
    "Three-compartment population PK model for midazolam with two equalized",
    "peripheral volumes and a three-transit-compartment first-order oral",
    "absorption chain (Ka = Ktr), supporting oral and intravenous dosing, in",
    "20 morbidly obese patients (mean total body weight 144 kg, range 112-186;",
    "mean BMI 47, range 40-68) and 12 non-obese healthy volunteers (mean total",
    "body weight 76 kg, mean BMI 22). Total body weight enters as a linear",
    "covariate on central volume (reference 127 kg) and a power covariate on",
    "peripheral volume (reference 127 kg); morbid-obesity status (BMI > 40)",
    "shifts oral bioavailability up and the transit absorption rate down."
  )
  reference <- paste(
    "Brill MJE, van Rongen A, Houwink API, Burggraaf J, van Ramshorst B,",
    "Wiezer RJ, van Dongen EPA, Knibbe CAJ (2014). Midazolam Pharmacokinetics",
    "in Morbidly Obese Patients Following Semi-Simultaneous Oral and",
    "Intravenous Administration: A Comparison with Healthy Volunteers.",
    "Clin Pharmacokinet 53(10):931-941. doi:10.1007/s40262-014-0166-x.",
    sep = " "
  )
  vignette <- "Brill_2014_midazolam"
  units    <- list(time = "minute", dosing = "microgram", concentration = "microgram/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'TBW' (total body weight, kg) maps to the canonical WT.",
        "Time-fixed at baseline. Enters the model with two distinct functional",
        "forms, both centred at 127 kg (the published reference): linear-",
        "deviation on the central volume Vc, Vc_TV = 44.1 * (1 + 0.0105 * (WT - 127));",
        "power on the peripheral volume Vp, Vp_TV = 139 * (WT / 127)^3.06.",
        "127 kg is the median total body weight across the pooled 32-subject",
        "cohort and is hardcoded in the .docx supplementary control stream",
        "($PK lines `TVV2= THETA(3)*(1+THETA(8)*(TBW-127))` and",
        "`TVV3= THETA(5)*(TBW/127)**THETA(9)`)."
      ),
      source_name        = "TBW"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = "BMI <= 40 (non-morbidly-obese reference cohort).",
      notes              = paste(
        "Used as a binary stratifier `OBES = as.integer(BMI > 40)` per Brill",
        "2014 Methods 2.1 (morbid obesity inclusion threshold) and the .docx",
        "supplementary control stream's OBES indicator in the $INPUT block",
        "(values 0 / 1 for healthy volunteers and morbidly obese patients).",
        "OBES enters oral bioavailability and the transit absorption rate",
        "as a multiplicative shift: F = F_HV * (F_MO/F_HV)^OBES and",
        "Ka = Ka_HV * (Ka_MO/Ka_HV)^OBES per supplementary $PK lines",
        "`TVF1= THETA(2)*THETA(10)**OBES` and `TVKA= THETA(6)*THETA(11)**OBES`.",
        "The lidocaine NA_NA_lidocaine.R model (DDMODEL00000281) sets the",
        "registered-precedent for using BMI with a binary threshold inside",
        "model() rather than registering a separate canonical obesity-status",
        "column."
      ),
      source_name        = "BMI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 2L,
    age_range      = "Morbidly obese 26-57 years (mean 43.6, SD 7.6); healthy volunteers 18-27 years (mean 22.0, SD 3.1)",
    weight_range   = "Morbidly obese 112-186 kg (mean 144.4, SD 21.7); healthy volunteers 63-93 kg (mean 76.0, SD 8.7)",
    weight_median  = "127 kg (pooled-cohort reference for covariate centring)",
    sex_female_pct = 37.5,
    disease_state  = paste(
      "Morbidly obese cohort: BMI > 40 kg/m^2 (mean 47.1, range 40-68),",
      "undergoing laparoscopic gastric bypass or sleeve surgery. Healthy-",
      "volunteer control cohort: non-obese (mean BMI 22.3, range 19-26).",
      "Both cohorts had no CYP3A-inducing or CYP3A-inhibiting co-medications",
      "and no renal insufficiency (eGFR >= 60 mL/min). Liver function markers",
      "in morbidly obese subjects were within three times the upper limit of",
      "normal."
    ),
    dose_range     = paste(
      "Morbidly obese: 7.5 mg midazolam oral tablet (Dormicum, Roche), then",
      "~159 min later 5 mg IV bolus (Midazolam Actavis 5 mg/mL).",
      "Healthy volunteers: 2 mg midazolam oral solution (Synthon), then",
      "150 min later 1 mg IV (Midazolam, Synthon)."
    ),
    regions        = "Netherlands",
    notes          = paste(
      "Pooled analysis of 20 morbidly obese surgical patients (NTC01519726 /",
      "EudraCT 2011-003293-93, St. Antonius Hospital, Nieuwegein) and 12",
      "non-obese male healthy volunteers (EudraCT 2009-010331-40) receiving",
      "midazolam in a semi-simultaneous oral - intravenous dosing design.",
      "Per Table 1 of Brill 2014: 20 + 12 = 32 subjects total. Sex split:",
      "12 F / 8 M morbidly obese, 0 F / 12 M healthy volunteers (12/32 = 37.5%",
      "female overall). 22 samples per patient over 11 h in the morbidly",
      "obese cohort and 19 samples per healthy volunteer; 42 of 434",
      "morbidly-obese samples were below the LLOQ (0.8 ng/mL) and were",
      "retained in the analysis dataset, with no LLOQ data from healthy",
      "volunteers."
    )
  )

  ini({
    # Structural parameters from Brill 2014 Table 2 (Final model column, RSE%).
    # The .docx supplementary control stream (ESM 4) names these THETA(1)-THETA(11)
    # with the parameterisation lines `TVF1= THETA(2)*THETA(10)**OBES`,
    # `TVV2= THETA(3)*(1+THETA(8)*(TBW-127))`, `TVV3= THETA(5)*(TBW/127)**THETA(9)`,
    # and `TVKA= THETA(6)*THETA(11)**OBES`.

    # Clearance (no covariate effect; Table 2 final-model column).
    lcl     <- log(0.359);   label("Clearance (L/min)")                                              # Brill 2014 Table 2 (CL = 0.359 L/min, final)

    # Oral bioavailability for the non-obese reference cohort (THETA(2)).
    # The OBES-conditional shift is captured separately as e_obes_fdepot
    # (= THETA(10)) and combined inside model() as fdepot * Z_F^OBES.
    lfdepot <- log(0.284);   label("Oral bioavailability in non-obese healthy volunteers (unitless)") # Brill 2014 Table 2 (F healthy volunteers = 0.284, final)

    # Central volume Vc for the median (127 kg) subject (THETA(3)). The
    # additive WT-slope coefficient is captured separately as e_wt_vc
    # (= THETA(8)) and combined inside model() as vc * (1 + e_wt_vc * (WT - 127)).
    lvc     <- log(44.1);    label("Central volume of distribution at WT = 127 kg (L)")              # Brill 2014 Table 2 (V_central,127kg = 44.1 L, final)

    # Inter-compartmental clearance to peripheral1 (THETA(4)).
    lq      <- log(1.33);    label("Inter-compartmental clearance to peripheral1 (L/min)")           # Brill 2014 Table 2 (Q = 1.33 L/min, final)

    # Peripheral volume Vp for the median (127 kg) subject (THETA(5)). The
    # power WT-exponent is captured separately as e_wt_vp (= THETA(9)) and
    # combined inside model() as vp * (WT/127)^e_wt_vp. Peripheral2 volume
    # was constrained equal to peripheral1 (V4 = V3 in the .docx ESM 4
    # control stream); a single lvp parameter serves both.
    lvp     <- log(139);     label("Peripheral volume of distribution at WT = 127 kg (L; shared by peripheral1 and peripheral2)") # Brill 2014 Table 2 (V_peripheral,127kg = 139 L, final)

    # Inter-compartmental clearance to peripheral2 (THETA(7) in the .docx
    # ESM 4; called Q2 in the paper, Q4 in the control stream).
    lq2     <- log(0.15);    label("Inter-compartmental clearance to peripheral2 (L/min)")           # Brill 2014 Table 2 (Q2 = 0.15 L/min, final)

    # Absorption / transit-compartment rate constant for non-obese subjects
    # (THETA(6)). Ka = Ktr per Brill 2014 Methods 2.5; the .docx ESM 4
    # `KTR = KA` line confirms.
    lka     <- log(0.13);    label("First-order absorption / transit rate constant Ka = Ktr in non-obese subjects (1/min)") # Brill 2014 Table 2 (Ka = Ktr healthy volunteers = 0.130 1/min, final)

    # Covariate effects (paper-named THETA(8)-THETA(11)). Linear-deviation
    # coefficient on Vc, power exponent on Vp, multiplicative shift on F
    # and Ka for the morbidly obese (OBES = 1) cohort.
    e_wt_vc      <- 0.0105;             label("Linear-deviation coefficient of WT on central volume per kg (1/kg)")             # Brill 2014 Table 2 (Z = 0.0105, final)
    e_wt_vp      <- 3.06;               label("Power exponent of (WT / 127) on peripheral volume (unitless)")                   # Brill 2014 Table 2 (W = 3.06, final)
    e_obes_fdepot <- log(0.603 / 0.284); label("Log-additive ADULT OBES (BMI > 40) effect on oral bioavailability (unitless)")  # Brill 2014 Table 2 (F morbid 0.603 / F HV 0.284 = THETA(10) = 2.124)
    e_obes_ka    <- log(0.057 / 0.13);  label("Log-additive ADULT OBES (BMI > 40) effect on transit-rate constant (unitless)")  # Brill 2014 Table 2 (Ka morbid 0.057 / Ka HV 0.130 = THETA(11) = 0.4385)

    # IIV (final model). Brill 2014 Table 2 reports the IIVs as CV%.
    # Conversion to omega^2 on the log-eta scale: omega^2 = log(CV^2 + 1).
    #   CL  : 18.1% -> log(0.181^2 + 1) = 0.03224
    #   F   : 26.4% -> log(0.264^2 + 1) = 0.06738
    #   Ka  : 41.4% -> log(0.414^2 + 1) = 0.15833
    #   Vc  : 55.2% -> log(0.552^2 + 1) = 0.26601
    #   Vp  : 34.4% -> log(0.344^2 + 1) = 0.11189
    # The off-diagonal "Correlation between eta Vcentral and Vperipheral"
    # 0.12 is the correlation coefficient rho; the .docx ESM 4 implements
    # this via $OMEGA BLOCK(2). The block-covariance in linear units is
    # rho * sqrt(omega^2_Vc * omega^2_Vp) = 0.12 * sqrt(0.26601 * 0.11189) = 0.02070.
    # Per the .docx ESM 4, IIV on Q was held fixed at 0 (the $OMEGA line
    # `0 FIX ;Q` re-uses the ETA(3) slot in $PK as a structural placeholder
    # for the Q row); reproduced here as a fixed-zero eta to preserve
    # provenance.
    etalcl     ~ 0.03224                                                  # Brill 2014 Table 2 (CL IIV 18.1% CV)
    etalfdepot ~ 0.06738                                                  # Brill 2014 Table 2 (F IIV 26.4% CV)
    etalq      ~ fixed(0)                                                 # Brill 2014 .docx ESM 4 $OMEGA `0 FIX ;Q`
    etalka     ~ 0.15833                                                  # Brill 2014 Table 2 (Ka IIV 41.4% CV)
    etalvc + etalvp ~ c(0.26601, 0.02070, 0.11189)                        # Brill 2014 Table 2 ($OMEGA BLOCK(2); rho = 0.12)

    # Residual error. The .docx ESM 4 $ERROR block defines two distinct
    # proportional errors via OBES indicator: Y1 = IPRED*(1+ERR(1)) for
    # healthy volunteers (OBES==0) and Y2 = IPRED*(1+ERR(2)) for morbidly
    # obese (OBES==1). Both are proportional only; the additive term shown
    # for the simple-model column of Table 2 (3.1%) was dropped in the
    # final model. Final-column CV% are 10.0% (HV) and 46.7% (MO);
    # propSd is the linear-scale standard deviation.
    propSd_hv <- 0.100;   label("Proportional residual error in non-obese subjects (SD, fraction)") # Brill 2014 Table 2 (HV final, 10.0% CV)
    propSd_mo <- 0.467;   label("Proportional residual error in morbidly obese subjects (SD, fraction)") # Brill 2014 Table 2 (MO final, 46.7% CV)
  })
  model({
    # Binary morbid-obesity stratifier derived from BMI per the .docx
    # ESM 4 OBES column convention (Brill 2014 Methods 2.1: BMI > 40 kg/m^2
    # inclusion threshold; the published study had no subjects in the gap
    # between 26 and 40, so OBES = 1 collapses to the morbidly-obese
    # cohort).
    is_obes <- 1.0 * (BMI > 40)

    # Individual PK parameters. The .docx ESM 4 stream applies covariate
    # effects to typical values (TV*) and exponentiates eta on top:
    #   TVF1 = THETA(2) * THETA(10)^OBES;    F1 = TVF1 * exp(ETA(2))
    #   TVV2 = THETA(3) * (1 + THETA(8) * (TBW - 127)); V2 = TVV2 * exp(ETA(5))
    #   TVV3 = THETA(5) * (TBW / 127)^THETA(9);          V3 = TVV3 * exp(ETA(6))
    #   TVKA = THETA(6) * THETA(11)^OBES;    KA = TVKA * exp(ETA(4))
    # Translated to nlmixr2 log-additive form, the OBES effects on F and Ka
    # become e_obes_fdepot = log(THETA(10)) and e_obes_ka = log(THETA(11))
    # added inside the log expression.
    cl     <- exp(lcl + etalcl)
    fdepot <- exp(lfdepot + etalfdepot + e_obes_fdepot * is_obes)
    vc     <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 127))
    q      <- exp(lq  + etalq)
    vp     <- exp(lvp + etalvp) * (WT / 127)^e_wt_vp
    q2     <- exp(lq2)
    ka     <- exp(lka + etalka + e_obes_ka * is_obes)
    ktr    <- ka                                       # Brill 2014 Methods 2.5 / .docx ESM 4 `KTR = KA`

    # Peripheral2 shares its volume with peripheral1 (.docx ESM 4 `V4 = V3`
    # line; Brill 2014 Methods 3.2 "two equalized volumes of distribution
    # of the peripheral compartments").
    vp2 <- vp

    # Micro-constants for the three-compartment central-disposition model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Compartment ordering: depot (oral input), then three transit
    # compartments preceding the central pool. The .docx ESM 4 $MODEL block
    # has the same 3 transits between PODOSE and CENTRAL with rate constant
    # KA on the depot transition and KTR ( = KA) on the transit-to-transit
    # and transit-to-central transitions. Oral doses target the depot; IV
    # doses (bolus or infusion) load central directly via cmt = central in
    # the event table.
    d/dt(depot)       <- -ka  * depot
    d/dt(transit1)    <-  ka  * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2  - ktr * transit3
    d/dt(central)     <-  ktr * transit3 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Oral bioavailability applied to the depot (IV doses bypass).
    f(depot) <- fdepot

    # Observed plasma midazolam concentration: dose in microgram, volume in
    # L -> microgram/L which is numerically equal to ng/mL.
    Cc <- central / vc

    # OBES-conditional proportional residual error (per .docx ESM 4 $ERROR
    # `Y = Y1*COM1 + Y2*COM2` with COM1 = (OBES == 0) and COM2 = (OBES == 1)).
    propSdEff <- propSd_hv * (1 - is_obes) + propSd_mo * is_obes
    Cc ~ prop(propSdEff)
  })
}
