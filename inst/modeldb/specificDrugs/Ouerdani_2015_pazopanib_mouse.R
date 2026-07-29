Ouerdani_2015_pazopanib_mouse <- function() {
  description <- "Preclinical (mouse, CB-17 SCID with CAKI-2 renal-cell carcinoma xenografts). Semi-mechanistic tumour growth and angiogenesis-inhibition (TGI) model for pazopanib (Ouerdani 2015): logistic tumour growth (state tumor_size) limited by a separately tracked vasculature-determined carrying capacity (state carrying_capacity), with an antiangiogenic drug effect on carrying-capacity loss (power form in AUC_PAZO) and a putative cytotoxic drug effect on tumour decay (exponentially declining resistance, no AUC effect after the cytotoxic exponent was fixed to 0)."
  reference <- paste(
    "Ouerdani A, Struemper H, Suttle AB, Ouellet D, Ribba B.",
    "Preclinical modeling of tumor growth and angiogenesis inhibition",
    "to describe pazopanib clinical effects in renal cell carcinoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(11):660-668.",
    "doi:10.1002/psp4.12001.",
    "Erratum/revised version published online 2015-11-12 (Table 1 was replaced;",
    "the on-disk PDF is the corrected version).",
    sep = " "
  )
  vignette <- "Ouerdani_2015_pazopanib_mouse"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; pazopanib exposure enters as the per-period AUC_PAZO covariate, not via a PK ODE)",
    concentration = "mm^3 (tumour volume measured by handheld caliper, computed as (length * width^2) / 2; not a drug concentration)"
  )

  covariateData <- list(
    AUC_PAZO = list(
      description        = "Per-period mean AUC of pazopanib driving the antiangiogenic and (formally) the cytotoxic drug-effect rates in the Ouerdani 2015 TGI model.",
      units              = "ug*h/mL (= mg*h/L)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per dose group in the preclinical CAKI-2 xenograft experiment. Ouerdani 2015 cites the FDA Pharmacology Review for pazopanib NDA 022465 for the per-dose-group AUC values: 220.2 ug*h/mL at 10 mg/kg/day, 656.8 ug*h/mL at 30 mg/kg/day, and 1140.8 ug*h/mL at 100 mg/kg/day. Set to 0 for vehicle-control mice; the model() block gates both drug-effect rates on AUC_PAZO > 0 so that vehicle animals follow the pure-growth dynamics (a = 0, c = 0 at AUC_PAZO = 0).",
      source_name        = "AUC"
    ),
    TUM_VOL = list(
      description        = "Per-animal observed baseline tumour volume at randomisation; used as the per-subject initial condition for the tumor_size ODE state.",
      units              = "mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ouerdani 2015 sets P0 to the observed baseline tumour volume per animal so only six structural parameters are estimated. Range at randomisation: 100-250 mm^3 (paper Methods, preclinical data section).",
      source_name        = "P0"
    )
  )

  population <- list(
    species        = "mouse (female CB-17 SCID with subcutaneous CAKI-2 renal-cell carcinoma xenograft)",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "8-10 weeks at randomisation",
    weight_range   = "not reported in the modelling paper",
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state  = "renal-cell carcinoma CAKI-2 subcutaneous xenograft; tumour volumes 100-250 mm^3 at randomisation",
    dose_range     = "vehicle, 10, 30, or 100 mg/kg pazopanib once daily by oral gavage for 24 days",
    regions        = "preclinical (in-vivo xenograft); housing in specific-pathogen-free environments",
    notes          = "8 mice per dose group (4 dose groups, including vehicle) = 32 mice. Tumour volumes measured twice weekly by handheld caliper; calculated as (length * width^2) / 2. The paper's preclinical Methods state that 8 observations of each type (tumour volume and body weight) were collected per mouse; the model fits the tumour-volume time course only. Drug exposure (AUC_PAZO) is reported per dose group (not per animal) from a separate preclinical PK study cited as FDA NDA 022465 Pharmacology Review."
  )

  ini({
    # Structural parameters -- preclinical column of Ouerdani 2015 Table 1.
    # Note that the on-disk PDF is the November 2015 corrected version of
    # Table 1 (the paper's first-page footnote records that the original
    # Table 1 was replaced on 2015-11-12).
    lk_tumor  <- log(0.166);   label("Tumour growth rate k (1/day)")                              # Ouerdani 2015 Table 1 preclinical k = 0.166 (RSE 24%)
    lk_cap    <- log(0.0183);  label("Carrying-capacity rate constant b (1/day)")                  # Ouerdani 2015 Table 1 preclinical b = 0.0183 (RSE 58%); IIV fixed to 0 per Table 1
    lk_aa0    <- log(0.007);   label("Baseline antiangiogenic effect rate c0 (1/day)")             # Ouerdani 2015 Table 1 preclinical c = 0.007 (RSE 29%)
    lk_cyto0  <- log(0.251);   label("Baseline cytotoxic effect rate a0 (1/day)")                  # Ouerdani 2015 Table 1 preclinical a = 0.251 (RSE 13%)
    lk_res    <- log(0.196);   label("Cytotoxic resistance rate d (1/day; multiplies exp(-d*t))")  # Ouerdani 2015 Table 1 preclinical d = 0.196 (RSE 26%)
    lK0       <- log(543);     label("Initial carrying capacity K0 (mm^3)")                        # Ouerdani 2015 Table 1 preclinical K0 = 543 (RSE 15%)

    # Drug-exposure covariate effects on the two drug-effect rates. The
    # antiangiogenic exponent was estimated; the cytotoxic exponent was
    # initially estimated at 0.0002 and then fixed to 0 (paper text:
    # "the level of pazopanib exposure had no impact on the cytotoxic effect").
    e_auc_pazo_k_aa   <- 0.332;        label("Power exponent of AUC_PAZO on antiangiogenic rate c") # Ouerdani 2015 Table 1 preclinical b_c = 0.332 (RSE 15%)
    e_auc_pazo_k_cyto <- fixed(0);     label("Power exponent of AUC_PAZO on cytotoxic rate a")      # Ouerdani 2015 Table 1 preclinical b_a fixed to 0 (initially estimated at 0.0002, then fixed)

    # Empirical structural exponent on tumour-driven capacity growth dK/dt = b * P^n - c * K.
    # Fixed at 1 for the mouse fit per the paper (likelihood-profile-tested
    # values: 0.5, 2/3, 1, 1.5, 2, 2.5, 3; best fit at n = 1).
    n <- fixed(1); label("Empirical exponent on tumour volume in capacity ODE")                    # Ouerdani 2015 Methods (Equation 1) and Results: best preclinical fit at n = 1

    # Inter-individual variability (IIV) -- Ouerdani 2015 Table 1 footnote
    # reports IIV as the omega (square root of NONMEM variance) expressed
    # as a percentage. Variance entered as omega^2 = (IIV / 100)^2.
    etalk_tumor  ~ 0.2809   # (0.53)^2; preclinical IIV on k = 53% (RSE 108%)
    etalk_aa0    ~ 0.0361   # (0.19)^2; preclinical IIV on c0 = 19% (RSE 127%)
    etalk_cyto0  ~ 0.0576   # (0.24)^2; preclinical IIV on a0 = 24% (RSE 96%)
    etalk_res    ~ 0.1764   # (0.42)^2; preclinical IIV on d  = 42% (RSE 103%)
    etalK0       ~ 0.1296   # (0.36)^2; preclinical IIV on K0 = 36% (RSE 77%)
    # No etalk_cap: Ouerdani 2015 Table 1 fixes IIV on b to 0 (paper text:
    # "Random effects were assumed for all model parameters except for
    # parameter b, as its estimation was associated with numerical instabilities").

    # Residual error -- combined additive + proportional on tumour volume.
    # Paper footnote on Table 1: "e1 is presented as a percentage, whereas
    # e2 is the standard deviation in the unit of the observed variable
    # (mm^3 ... for preclinical ...)". So 14% maps to a proportional SD of
    # 0.14 and 3 maps to an additive SD of 3 mm^3.
    propSd <- 0.14;  label("Proportional residual SD on tumour volume (fraction)")               # Ouerdani 2015 Table 1 preclinical e1 = 14% (RSE 23%)
    addSd  <- 3;     label("Additive residual SD on tumour volume (mm^3)")                       # Ouerdani 2015 Table 1 preclinical e2 = 3 (RSE 17%)
  })

  model({
    # ----- Individual structural parameters (lognormal IIV where present) -----
    k_tumor <- exp(lk_tumor + etalk_tumor)
    k_cap   <- exp(lk_cap)                                  # no IIV (Table 1 footnote)
    a0      <- exp(lk_cyto0 + etalk_cyto0)
    c0      <- exp(lk_aa0   + etalk_aa0)
    k_res   <- exp(lk_res   + etalk_res)
    K0      <- exp(lK0      + etalK0)

    # ----- AUC-driven drug-effect rates -----
    # Ouerdani 2015 Methods Equations 2 and 3:
    #     a = a0 * AUC^b_a       (cytotoxic effect rate, 1/day)
    #     c = c0 * AUC^b_c       (antiangiogenic effect rate, 1/day)
    # In the mouse fit, b_a = 0 (fixed) so the power form collapses to
    # a = a0 mathematically for AUC > 0. To keep the model physically
    # sensible at AUC = 0 (vehicle control mice), both rates are gated on
    # AUC_PAZO > 0 -- a and c both vanish in the absence of drug, so
    # vehicle animals follow pure logistic tumour growth and unimpeded
    # capacity expansion. Without this gate, AUC^0 = 1 (R / NONMEM
    # convention) would keep a non-zero cytotoxic rate active in vehicle
    # animals, which contradicts the paper's vehicle-arm tumour trajectory.
    if (AUC_PAZO > 0) {
      cyto_rate <- a0 * AUC_PAZO^e_auc_pazo_k_cyto
      aa_rate   <- c0 * AUC_PAZO^e_auc_pazo_k_aa
    } else {
      cyto_rate <- 0
      aa_rate   <- 0
    }

    # ----- ODE system (Ouerdani 2015 Equation 1, mouse fit with n = 1) -----
    #   dP/dt = k * P * (1 - P/K) - a * exp(-d*t) * P
    #   dK/dt = b * P^n           - c * K
    # P  = tumor_size         (mm^3)
    # K  = carrying_capacity  (mm^3); biologically the maximum tumour volume
    #      sustainable by current vasculature.
    # t  is rxode2 simulation time in days; treatment start is taken as t = 0
    # so the exp(-d*t) resistance term decays from 1 at study start.
    d/dt(tumor_size)        <- k_tumor * tumor_size * (1 - tumor_size / carrying_capacity) -
                              cyto_rate * exp(-k_res * t) * tumor_size
    d/dt(carrying_capacity) <- k_cap * tumor_size^n - aa_rate * carrying_capacity

    # Initial state. P0 is per-animal (covariate TUM_VOL); K0 is a
    # population-typical-value parameter with IIV (etalK0).
    tumor_size(0)        <- TUM_VOL
    carrying_capacity(0) <- K0

    # ----- Observation and combined residual error -----
    # Tumour volume is observed in mm^3 by caliper. Combined additive +
    # proportional residual error per Table 1 footnote.
    tumor_size ~ add(addSd) + prop(propSd)
  })
}
