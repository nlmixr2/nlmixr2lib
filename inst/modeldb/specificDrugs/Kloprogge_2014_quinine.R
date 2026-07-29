# Population pharmacokinetic model for oral quinine in pregnant Ugandan
# women with uncomplicated Plasmodium falciparum malaria
# (Kloprogge 2014, J Antimicrob Chemother 69(11):3033-3040;
# doi:10.1093/jac/dku228).

Kloprogge_2014_quinine <- function() {
  description <- paste(
    "Population PK model for oral quinine in pregnant women with",
    "uncomplicated Plasmodium falciparum malaria in Uganda",
    "(Kloprogge 2014). First-order absorption into a two-compartment",
    "disposition model with allometric body-weight scaling on clearance",
    "and intercompartmental clearance (power 2/3) and on apparent volumes",
    "(power 1), centered at the cohort typical weight of 56 kg. Relative",
    "bioavailability F is fixed at 1 with log-normal IIV; a linear",
    "covariate effect of time-varying parasitaemia (per log10",
    "parasites/uL, last-observation-carried-forward) increases F by 38.9%",
    "per log10 parasitaemia, and an exponential effect of admission body",
    "temperature decreases elimination clearance by ~21.6% per degC",
    "(centered at the cohort median 37.2 degC).",
    sep = " "
  )
  reference <- paste(
    "Kloprogge F, Jullien V, Piola P, Dhorda M, Muwanga S, Nosten F,",
    "Day NPJ, White NJ, Guerin PJ, Tarning J (2014).",
    "Population pharmacokinetics of quinine in pregnant women with",
    "uncomplicated Plasmodium falciparum malaria in Uganda.",
    "Journal of Antimicrobial Chemotherapy 69(11):3033-3040.",
    "doi:10.1093/jac/dku228.",
    sep = " "
  )
  vignette <- "Kloprogge_2014_quinine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Kloprogge 2014 enrolled 22",
        "pregnant Ugandan women with median (range) body weight 56.5",
        "(44.0-71.0) kg (Table 1). Body weight is applied as an",
        "allometric scaler with power 2/3 on CL/F and Q/F and power 1 on",
        "Vc/F and Vp/F:",
        "(WT/56)^(2/3) on clearance parameters,",
        "(WT/56)^1 on volume parameters.",
        "The reference weight 56 kg is the typical-patient value used in",
        "the paper's Monte Carlo simulations (Figures 4 and 5) and is",
        "consistent with the median observed cohort weight. The paper",
        "does not state the reference weight explicitly; this",
        "implementation uses 56 kg so that the reported population",
        "estimates (CL/F = 10.4 L/h, Vc/F = 174 L, Q/F = 10.7 L/h,",
        "Vp/F = 54.3 L) reproduce at the typical patient.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    BODYTEMP = list(
      description        = "Body temperature at admission",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission body temperature, time-fixed per subject. Kloprogge",
        "2014 cohort median 37.2 degC (Table 1, range 36.0-38.9 degC).",
        "Applied as an exponential effect on elimination clearance:",
        "CL_indiv = TVCL * exp(e_bodytemp_cl * (BODYTEMP - 37.2)) with",
        "e_bodytemp_cl = -0.243 per degC, centered at the cohort median.",
        "Higher admission temperature is associated with lower apparent",
        "elimination clearance, reflecting acute-phase reductions in",
        "CYP3A4 metabolic activity (Discussion). The effect was retained",
        "across the entire 7-day treatment course in the source paper",
        "(time-restricted versions over the first 24 / 48 h gave",
        "significantly worse fits). Effect should not be extrapolated",
        "outside the observed 36.0-38.9 degC range.",
        sep = " "
      ),
      source_name        = "TEMP"
    ),
    PARA = list(
      description        = "Plasmodium falciparum parasitaemia (asexual parasites/uL)",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying covariate, propagated by last-observation",
        "-carried-forward (Methods, Discussion). Parasite slides were",
        "taken once daily during the 7-day treatment course; the most",
        "recent observed parasitaemia is used at each time point in the",
        "model. Kloprogge 2014 admission median 2240 (range 39-44500)",
        "parasites/uL (Table 1). Applied as a linear effect on relative",
        "bioavailability with the log10 transform applied inside the",
        "model:",
        "F_indiv = (1 + e_para_f * log10(max(PARA, 1))) *",
        "exp(lfdepot + etalfdepot)",
        "with e_para_f = +0.389 per log10 parasitaemia (38.9% increase",
        "per log10 parasites/uL). The max(PARA, 1) gate enforces the",
        "paper's wording that the covariate effect applies only during",
        "the acute phase when parasitaemia is above detection: values of",
        "PARA <= 1 (effectively zero / below LOQ) collapse to no",
        "covariate effect (F = F_typ). Distinct from the LNPC canonical",
        "(natural-log admission-only parasitaemia, used by Birgersson",
        "2019); PARA is the raw count and is time-varying here.",
        sep = " "
      ),
      source_name        = "PARA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 22L,
    n_studies       = 1L,
    age_range       = "18.0-37.0 years (Table 1)",
    age_median      = "21.0 years (Table 1)",
    weight_range    = "44.0-71.0 kg (Table 1)",
    weight_median   = "56.5 kg (Table 1)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria; second or third",
      "trimester of pregnancy (estimated gestational age 13.0-37.0",
      "weeks, median 26.0; trimester split 12/22 second, 10/22 third).",
      "All patients were enrolled at the Mbarara National Referral",
      "Hospital antenatal clinic in Uganda."
    ),
    dose_range      = paste(
      "Oral quinine sulphate (Remedica, Limassol, Cyprus): 300 mg",
      "salt/tablet; 10 mg salt/kg per dose, three times daily at 0, 8,",
      "and 16 h for 7 days. Salt-to-base conversion: quinine sulphate",
      "molecular weight 782.96 g/mol, quinine base 324.42 g/mol; the",
      "model dose should be supplied in quinine base equivalents",
      "(multiply quinine sulphate mg by 324.42/782.96 = 0.4144). Full /",
      "half replacement dose given if vomiting occurred within 30 min /",
      "30-60 min of intake."
    ),
    regions         = "Uganda (Mbarara National Referral Hospital antenatal clinic)",
    trial_registration = "ClinicalTrials.gov NCT00495508",
    notes           = paste(
      "Demographics from Kloprogge 2014 Table 1. Twenty-three women were",
      "enrolled in the pharmacokinetic study; one was excluded from the",
      "population analysis due to an unexplainable mismatch between",
      "dosing history and plasma concentration profile, leaving 22.",
      "Co-medications were ferrous sulphate + folic acid (n=1), unknown",
      "(n=2), and amoxicillin (n=1); none was expected to affect",
      "quinine pharmacokinetics. Sampling: 0, 1, 2, 3, 4, 8, 16, 24, 48,",
      "72, 96, 120, 144, 160, 161, 162, 163, 164, 168, 170, 172, 176",
      "and 184 h after the first dose (dense around first and last",
      "doses; daily during steady-state portion).",
      "Bioanalytical assay: LC with fluorimetric detection, LLOQ = 1",
      "ug/mL, accuracy < 5%, precision (RSD) < 9.9%."
    )
  )

  ini({
    # Structural population parameters come from Kloprogge 2014 Table 2,
    # "Population estimate" column. Parameters are reported on the
    # linear (back-transformed) scale; log() is applied here for the
    # nlmixr2 internal scale. Population estimates correspond to a
    # typical 56-kg patient with admission body temperature 37.2 degC
    # (cohort median) and parasitaemia centered at log10(PARA) = 0
    # (i.e. F_typ = 1 when PARA = 1 parasite/uL).
    lka  <- log(0.817)  ; label("Absorption rate constant ka (1/h)")                          # Kloprogge 2014 Table 2: ka     = 0.817 (RSE 18.8%; 95% CI 0.479-1.03)
    lcl  <- log(10.4)   ; label("Apparent elimination clearance CL/F (L/h, typical 56-kg)")    # Kloprogge 2014 Table 2: CL/F   = 10.4 (RSE 4.36%; 95% CI 9.51-11.4)
    lvc  <- log(174)    ; label("Apparent central volume of distribution Vc/F (L, typical 56-kg)") # Kloprogge 2014 Table 2: Vc/F = 174 (RSE 14.0%; 95% CI 112-195)
    lq   <- log(10.7)   ; label("Apparent intercompartmental clearance Q/F (L/h, typical 56-kg)") # Kloprogge 2014 Table 2: Q/F  = 10.7 (RSE 44.6%; 95% CI 7.06-36.9)
    lvp  <- log(54.3)   ; label("Apparent peripheral volume of distribution Vp/F (L, typical 56-kg)") # Kloprogge 2014 Table 2: Vp/F = 54.3 (RSE 29.1%; 95% CI 33.6-112)

    # Relative bioavailability anchored at 1 (structural, fixed by the
    # source paper, with all variability around F captured by IIV).
    # The paper additionally reports between-dose IOV on F (CV 21.4%);
    # IOV is omitted here because it requires per-dose occasion indices
    # that are not part of the typical simulation workflow (see vignette
    # Errata). Modelling F_typ = 1 with IIV CV 12.3% reproduces the
    # between-subject component of the source variability structure.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless, fixed)")          # Kloprogge 2014 Table 2: F = 100 (fixed)

    # Allometric exponents. Fixed at the canonical Mahidol-Oxford
    # malaria-popPK values used by the source paper: 2/3 on clearance
    # (CL, Q) and 1 on volume (Vc, Vp). Kloprogge 2014 explicitly
    # selected the 2/3 exponent on clearance over the alternative 3/4
    # ("a power coefficient of 2/3 on clearance parameters produced a
    # better fit of the model compared with a coefficient of 3/4 ...
    # in good agreement with the observed physiology since clearance
    # does not normally scale linearly with body weight").
    allo_cl  <- fixed(2/3) ; label("Allometric exponent on CL/F and Q/F (unitless, fixed)")    # Kloprogge 2014 Methods + Results paragraph 2
    allo_vc  <- fixed(1)   ; label("Allometric exponent on Vc/F and Vp/F (unitless, fixed)")   # Kloprogge 2014 Methods + Results paragraph 2

    # Covariate effects.
    # Body temperature is encoded as an exponential effect on CL/F,
    # centered at the cohort median 37.2 degC (Table 1). The Kloprogge
    # 2014 form gives a 51.8% lower CL/F at 39 degC vs 36 degC
    # (exp(-0.243 * 3) = 0.482).
    e_bodytemp_cl <- -0.243 ; label("Exponential effect of admission body temperature on CL/F (per degC, centered at 37.2 degC)") # Kloprogge 2014 Table 2: Temperature on CL/F = -0.243 (RSE 21.1%; 95% CI -0.427 to -0.180)
    # Parasitaemia is encoded as a linear effect on F with log10
    # transform applied inside model(). 38.9% increase per log10
    # parasites/uL. Reference parasitaemia is log10(PARA) = 0,
    # i.e. F_typ = 1 when PARA = 1 parasite/uL; below-detection
    # PARA values are gated to PARA = 1 to collapse the effect to
    # F = F_typ (Discussion: "effect only during the acute phase when
    # parasitaemia was above the limit of detection").
    e_para_f      <-  0.389 ; label("Linear effect of parasitaemia on F (per log10 parasites/uL)")    # Kloprogge 2014 Table 2: Parasitaemia (log10) on F = 38.9 (RSE 9.33%; 95% CI 32.4-47.2)

    # Inter-individual variability. Kloprogge 2014 Table 2 reports CV%
    # (footnote: 100 * sqrt(exp(omega^2) - 1)). The internal log-normal
    # variances are recovered by omega^2 = log((CV/100)^2 + 1):
    #   ka  CV 58.7% -> log(0.587^2 + 1) = 0.296221
    #   CL  CV  7.69% -> log(0.0769^2 + 1) = 0.005893
    #   Vp  CV 70.8% -> log(0.708^2 + 1)  = 0.406119
    #   F   CV 12.3% (IIV only) -> log(0.123^2 + 1) = 0.015015
    # Vc/F and Q/F carry no IIV in the source model (Table 2 dashes).
    etalka     ~ 0.296221  # Kloprogge 2014 Table 2: BSV on ka  = 58.7% CV (RSE 32.7%; 95% CI 40.5-107%)
    etalcl     ~ 0.005893  # Kloprogge 2014 Table 2: BSV on CL  =  7.69% CV (RSE 65.4%; 95% CI 1.16-47.4%)
    etalvp     ~ 0.406119  # Kloprogge 2014 Table 2: BSV on Vp  = 70.8% CV (RSE 65.3%; 95% CI 8.00-128%)
    etalfdepot ~ 0.015015  # Kloprogge 2014 Table 2: BSV on F   = 12.3% CV (RSE 77.0%; 95% CI 0.170-48.8%). IOV component (21.4% CV) omitted; see model description and vignette Errata.

    # Residual error. Kloprogge 2014 modelled the natural logarithm of
    # the quinine plasma concentration with an additive error on log
    # scale, which is equivalent to proportional error in nlmixr2's
    # linear-concentration space (see references/parameter-names.md
    # 'Residual error' and the Table 2 footnote "The additive error
    # variance will essentially be exponential on normal scale data").
    # Table 2 reports the additive variance on the log scale as 0.0158;
    # the corresponding SD is sqrt(0.0158) = 0.1257, encoded here as
    # propSd.
    propSd <- sqrt(0.0158) ; label("Proportional residual SD on linear concentration scale (= SD on log scale)") # Kloprogge 2014 Table 2: additive residual variance = 0.0158 (RSE 41.6%; 95% CI 0.0129-0.156)
  })

  model({
    # Individual structural parameters with allometric WT scaling and
    # the body-temperature exponential effect on CL/F.
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) * (WT / 56)^allo_cl * exp(e_bodytemp_cl * (BODYTEMP - 37.2))
    vc  <- exp(lvc)          * (WT / 56)^allo_vc
    q   <- exp(lq)           * (WT / 56)^allo_cl
    vp  <- exp(lvp + etalvp) * (WT / 56)^allo_vc

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order absorption from depot into a
    # two-compartment disposition model.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    # Relative bioavailability applied to the depot compartment. The
    # linear parasitaemia covariate is applied to F via the log10
    # transform inside the model; max(PARA, 1) gates PARA values below
    # detection (effectively zero) to no covariate effect (F = F_typ).
    f(depot) <- (1 + e_para_f * log10(max(PARA, 1))) * exp(lfdepot + etalfdepot)

    # Quinine plasma concentration in mg/L (== ug/mL). Dose units are
    # quinine-base mg, Vc is L, so central/vc has units mg/L. Confirm
    # that input doses have been converted from quinine sulphate to
    # quinine base by the factor 324.42/782.96 = 0.4144 (Methods).
    Cc <- central / vc

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale maps to proportional in nlmixr2;
    # see ini() comment on propSd).
    Cc ~ prop(propSd)
  })
}
