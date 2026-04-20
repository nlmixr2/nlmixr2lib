Tiraboschi_2025_amlitelimab <- function() {
  description <- "Two-compartment population PK model for amlitelimab (anti-OX40L mAb) in adults, with parallel first-order and Michaelis-Menten (TMDD) elimination, SC absorption with lag time, allometric body-weight scaling, and EASI / albumin covariate effects (Tiraboschi 2025)"
  reference <- "Tiraboschi JM, Zohar S, Quartino AL, Monnier R, Coulette V, Bizot JL, Jamois C. Population Pharmacokinetic and Pharmacodynamic Modeling for the Prediction of the Extended Amlitelimab Phase 3 Dosing Regimen in Atopic Dermatitis. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2161-2173. doi:10.1002/psp4.70121"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (not time-varying) per the Tiraboschi 2025 NONMEM control stream; allometric effects on CL, V1 (central), and V2 (peripheral) with reference weight 75 kg (population median).",
      source_name        = "BWT"
    ),
    EASI = list(
      description        = "Baseline Eczema Area and Severity Index score",
      units              = "(score, 0-72)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline EASI (BEASI in the source NONMEM code); enters linear clearance as an additive term 0.00111 * EASI (L/day). Healthy volunteers have EASI = 0. Renamed from source column BEASI to the canonical EASI per covariate-columns.md.",
      source_name        = "BEASI"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline albumin (BALB in the source NONMEM code); covariate on subcutaneous bioavailability via the additive term 0.598 * ((ALB/47) - 1) on the linear F1 scale before logit transformation. Reference value 47 g/L is the population median. Renamed from source column BALB to the canonical ALB per covariate-columns.md.",
      source_name        = "BALB"
    )
  )

  population <- list(
    n_subjects     = 439,
    n_studies      = 5,
    age_range      = "18-72 years",
    age_median     = "33 years",
    weight_range   = "40.5-148 kg",
    weight_median  = "74.5 kg",
    sex_female_pct = 38,
    race_ethnicity = c(White = 82.7, Black = 4.33, Asian = 11.2, Other = 1.82),
    disease_state  = "Pooled across 5 clinical studies: 78 healthy volunteers (phase 1) and 361 adults with moderate-to-severe atopic dermatitis (phase 2a n=59 and STREAM-AD phase 2b n=302); among AD subjects, baseline EASI mean 29.7 (SD 11.3) and 72.9% classified as severe (EASI > 21)",
    dose_range     = "Single or repeated IV and SC doses; labelled STREAM-AD regimens include 250 mg SC Q4W / Q12W with a 500 mg SC loading dose and 62.5 mg SC Q4W",
    regions        = "Multi-regional (STREAM-AD phase 2b primary driver; region breakdown not reported in the source)",
    notes          = "Demographics from Tiraboschi 2025 Table S1 (pooled PopPK analysis population, n=439). Three phase 1 studies in healthy volunteers (n=48 + 24 + 6) and two phase 2 studies in AD (n=59 phase 2a + n=302 STREAM-AD phase 2b). 80 of 439 subjects (18.2%) had at least one positive ADA sample; ADA was not a significant covariate in the final PopPK model."
  )

  ini({
    # Structural parameters — reference values for a 75 kg participant with albumin = 47 g/L; values from Tiraboschi 2025 Table S2
    lka      <- log(0.233);                      label("Absorption rate (Ka, 1/day)")                                                       # Table S2 TVKa
    lcl      <- log(0.115);                      label("Linear clearance for a 75 kg participant with EASI = 0 (CL, L/day)")                # Table S2 TVCLL
    lvc      <- log(3.46);                       label("Central volume of distribution for a 75 kg participant (V1, L)")                    # Table S2 TVV1
    lvp      <- log(2.48);                       label("Peripheral volume of distribution for a 75 kg participant (V2, L)")                 # Table S2 TVV2
    lq       <- log(0.569);                      label("Intercompartmental clearance (Q, L/day)")                                           # Table S2 TVQ2
    lvmax    <- log(0.0362);                     label("Maximum velocity of nonlinear (TMDD) elimination (Vmax, mg/day)")                   # Table S2 TVVM (labeled 'ug/day' in Table S2, verified mg/day via the paper's 66% TMDD-fraction at LLOQ 0.0469 ug/mL and 20% at 1 ug/mL)
    lkm      <- log(0.0783);                     label("Michaelis-Menten constant (Km, ug/mL)")                                             # Table S2 TVKM
    lalag    <- log(0.0351);                     label("Absorption lag time (ALAG, day)")                                                   # Table S2 TVALAG
    logitf1  <- log(0.888 / (1 - 0.888));        label("Typical subcutaneous bioavailability on the logit scale (linear F = 0.888 at population-median albumin 47 g/L)")  # Table S2 TVFsc = 0.888 on linear scale

    # Allometric exponents on body weight (reference 75 kg)
    allo_cl <- 1.20;  label("Allometric exponent on linear clearance (unitless)")   # Table S2 BWT effect on CL
    allo_v1 <- 0.901; label("Allometric exponent on central volume (unitless)")     # Table S2 BWT effect on V1
    allo_v2 <- 0.350; label("Allometric exponent on peripheral volume (unitless)")  # Table S2 BWT effect on V2

    # Covariate effects — both additive on the linear scale per the source NONMEM control stream
    e_easi_cl <- 0.00111; label("Additive EASI effect on linear CL (L/day per EASI unit)")                                                  # Table S2 BEASI effect on CL
    e_alb_f1  <- 0.598;   label("Additive albumin effect on Fsc on the linear scale, per unit of (ALB/47 - 1) (fraction)")                  # Table S2 BALB effect on Fsc

    # Inter-individual variability — V1 and CL are a correlated BLOCK(2); V2, Fsc (logit), ALAG, Ka are diagonal
    etalvc + etalcl ~ c(0.0491,
                        0.0240, 0.0482)  # Table S2 omega^2 V1, omega(V1,CL), omega^2 CL
    etalvp     ~ 0.0693  # Table S2 omega^2 V2
    etalogitf1 ~ 1.18    # Table S2 omega^2 Fsc (IIV applied on logit scale per the source NONMEM control stream)
    etalalag   ~ 0.151   # Table S2 omega^2 ALAG
    etalka     ~ 0.135   # Table S2 omega^2 Ka

    # Residual error (proportional; variance 0.0248 -> CV = sqrt(0.0248) = 0.1575 = 15.75%)
    propSd <- 0.1575; label("Proportional residual error (fraction)")  # Table S2 sigma^2 proportional = 0.0248
  })
  model({
    # Individual PK parameters
    ka   <- exp(lka + etalka)
    # Linear CL: allometric on the structural term + additive EASI effect (Tiraboschi 2025 Table S2 footnote e)
    cl   <- exp(lcl + etalcl) * (WT / 75)^allo_cl + e_easi_cl * EASI
    vc   <- exp(lvc + etalvc) * (WT / 75)^allo_v1
    vp   <- exp(lvp + etalvp) * (WT / 75)^allo_v2
    q    <- exp(lq)
    vmax <- exp(lvmax)
    km   <- exp(lkm)
    alag <- exp(lalag + etalalag)

    # Fsc: population on logit scale -> invert to linear, add additive albumin term, re-logit, then apply the logit-scale eta
    # (Tiraboschi 2025 Table S2 footnote f: Fsc = TVFsc + 0.598 * ((BALB/47) - 1); the source NONMEM code then uses
    # F1 = exp(PHI + ETAF1) / (1 + exp(PHI + ETAF1)) with PHI = log(Fsc / (1 - Fsc)).)
    f1_typ_lin <- exp(logitf1) / (1 + exp(logitf1))
    f1_cov_lin <- f1_typ_lin + e_alb_f1 * ((ALB / 47) - 1)
    phi_f1     <- log(f1_cov_lin / (1 - f1_cov_lin))
    f1         <- exp(phi_f1 + etalogitf1) / (1 + exp(phi_f1 + etalogitf1))

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL (defined before the ODE so the MM term's Cc is the pre-evaluated concentration)
    Cc <- central / vc

    # 2-compartment ODE system with parallel linear + Michaelis-Menten elimination from central
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - (vmax * Cc) / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- f1
    alag(depot) <- alag

    Cc ~ prop(propSd)
  })
}
