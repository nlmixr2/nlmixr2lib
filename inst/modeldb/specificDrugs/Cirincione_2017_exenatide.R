Cirincione_2017_exenatide <- function() {
  description <- "Population PK model for exenatide immediate-release (Cirincione 2017): two-compartment, parallel linear and Michaelis-Menten elimination, sequential zero-order then saturable first-order absorption after SC dosing."
  reference <- "Cirincione B, Mager DE. Population pharmacokinetics of exenatide. British Journal of Clinical Pharmacology. 2017;83(3):517-526. doi:10.1111/bcp.13135"
  units <- list(time = "hour", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric-like effect on central volume; reference weight 84.8 kg (population median).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Modification of Diet in Renal Disease (MDRD) estimated glomerular filtration rate (creatinine-based, BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear clearance; reference CRCL 80 mL/min/1.73 m^2. Source column 'eGFR' (MDRD eGFR) maps to the canonical general-scope CRCL covariate (which also accepts measured CrCl BSA-normalized in the same units); the MDRD estimation method is documented here in the description.",
      source_name        = "eGFR"
    ),
    STUDY1 = list(
      description        = "Indicator for Study 1 of the Cirincione 2017 pooled analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other studies)",
      notes              = "Selects the Study 1 log-scale residual error magnitude. Source used a character-valued DVID column with values 'study1' / 'study5' / 'otherStudy'; split into STUDY1 and STUDY5 binary indicators per covariate-columns.md register.",
      source_name        = "DVID"
    ),
    STUDY5 = list(
      description        = "Indicator for Study 5 of the Cirincione 2017 pooled analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other studies)",
      notes              = "Selects the Study 5 log-scale residual error magnitude. Paired with STUDY1; STUDY1 = STUDY5 = 0 selects the pooled 'other studies' residual error.",
      source_name        = "DVID"
    )
  )

  population <- list(
    n_subjects     = 195,
    n_studies      = 8,
    age_range      = "18-76 years",
    age_mean       = "48 years (SD 12.9)",
    weight_range   = "52.6-162 kg",
    weight_mean    = "86 kg (SD 16.2)",
    sex_female_pct = 22,
    race_ethnicity = c(Caucasian = 67, Other = 33),
    egfr_range     = "4.1 to >=90 mL/min/1.73 m^2",
    egfr_mean      = "88.7 mL/min/1.73 m^2 (SD 31.4)",
    renal_function = c(normal_ge_90 = 129, mild_60_89 = 47, moderate_30_59 = 10,
                       severe_15_29 = 1, ESRD_lt_15 = 8),
    disease_state  = "Pooled analysis: healthy volunteers and patients with type 2 diabetes mellitus, including normal renal function and mild/moderate/severe/ESRD renal impairment.",
    dose_range     = "2.5-20 ug single or repeated SC dose",
    regions        = "Not reported in detail",
    notes          = "Demographics from Cirincione 2017 Table 1. 195 subjects across 8 clinical trials; concentrations quantified by 2 immunoassays (Studies 1 and 5 had distinct residual error magnitudes)."
  )

  ini({
    # Structural parameters — Cirincione 2017 Table 2 (reference CRCL = 80 mL/min/1.73 m^2, reference WT = 84.8 kg)
    lcl       <- log(4.58);  label("Linear clearance at reference CRCL (L/hr)")                      # Table 2: Cl_int = 4.58 L/hr
    e_cl_crcl <- 0.838;      label("Power effect of CRCL (MDRD eGFR) on linear clearance (unitless)") # Table 2: Cl_eGFR = 0.838
    lq       <- log(3.72);   label("Intercompartmental clearance (L/hr)")                           # Table 2: Cld = 3.72 L/hr
    lvc      <- log(7.03);   label("Central volume at reference body weight (L)")                   # Table 2: Vc_int = 7.03 L
    e_vc_wt  <- 2.67;        label("Power effect of body weight on central volume (unitless)")       # Table 2: Vc_wtkg = 2.67
    lvp      <- log(7.04);   label("Peripheral volume (L)")                                          # Table 2: Vp = 7.04 L
    lvmax    <- log(1.55);   label("Maximum Michaelis-Menten clearance (ug/hr)")                     # Table 2: Vmax = 1.55 ug/hr
    lkm      <- log(0.567);  label("Michaelis-Menten constant for saturable clearance (ng/mL)")      # Table 2: Km = 567 pg/mL = 0.567 ng/mL (expressed in concentration units matching Cc = central(ug)/vc(L))
    lkamax   <- log(12.8);   label("Maximum first-order absorption rate (1/hr)")                     # Table 2: ka_max = 12.8 /hr (RSE 42.5%)
    lkmka    <- log(16.9);   label("Michaelis-Menten constant for saturable absorption (ug)")        # Table 2: Km_ka = 16.9 ug
    ttau     <- fixed(1.35); label("Duration of zero-order absorption (hr)")                         # Results p. 520: tau fixed at 1.35 hr
    fdepot   <- fixed(1);    label("Bioavailability fraction (fixed)")                               # Results p. 520: F fixed to 1 (absolute F not identifiable)
    logitfr  <- logit(0.628); label("Fraction of dose entering first-order absorption (fraction)")   # Table 2: fr = 0.628

    # IIV — Table 2 reports CV% (Cl_int = 33.9%, Vc_int = 80.5%, Km = 95.7%); log-normal variance = log(1 + CV^2)
    etalcl ~ log(1 + 0.339^2)                                                                        # Table 2: CV Cl_int = 33.9%
    etalvc ~ log(1 + 0.805^2)                                                                        # Table 2: CV Vc_int = 80.5%
    etalkm ~ log(1 + 0.957^2)                                                                        # Table 2: CV Km = 95.7%

    # Residual error — Table 2 / NONMEM $TABLE reports log-scale SDs per study
    expSdOther  <- 0.373; label("Log-scale residual SD, all studies other than 1 and 5")             # Table 2: residual sigma, other studies = 0.373
    expSdStudy1 <- 0.39;  label("Log-scale residual SD, Study 1 assay")                              # Table 2: residual sigma, Study 1 = 0.39
    expSdStudy5 <- 0.08;  label("Log-scale residual SD, Study 5 assay")                              # Table 2: residual sigma, Study 5 = 0.08
  })

  model({
    # Individual PK parameters
    cl    <- exp(lcl + etalcl) * (CRCL / 80)^e_cl_crcl
    q     <- exp(lq)
    vc    <- exp(lvc + etalvc) * (WT / 84.8)^e_vc_wt
    vp    <- exp(lvp)
    vmax  <- exp(lvmax)
    km    <- exp(lkm + etalkm)
    kamax <- exp(lkamax)
    kmka  <- exp(lkmka)
    fr    <- expit(logitfr)

    # Micro-constants
    k12 <- q / vc
    k21 <- q / vp
    # kel has a linear and a saturable component; the saturable component uses the central amount (ug) divided by
    # (km * vc + central) so that km is carried in concentration units consistent with Cc = central/vc.
    kel <- cl / vc + vmax / (km * vc + central)

    # Sequential absorption: zero-order over [0, tau], then saturable first-order
    mtime(tau) <- ttau

    kzero <- (1 - fr) * podo(depot) / tau
    if (tad(depot) > tau) kzero <- 0.0

    ka <- fr * kamax / (kmka + depot)
    if (tad(depot) <= tau) ka <- 0.0

    d/dt(depot)       <- -ka * depot - kzero
    d/dt(central)     <-  ka * depot + kzero - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                       k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # central is in ug, vc is in L, so central/vc is in ug/L = ng/mL. Multiply by 1000 to report Cc in pg/mL
    # for comparison with Table 2 / Figure 5.
    Cc <- (central / vc) * 1000

    # Residual error: a single model switched by binary study indicators (STUDY1, STUDY5)
    expSd <- expSdStudy1 * STUDY1 + expSdStudy5 * STUDY5 + expSdOther * (1 - STUDY1 - STUDY5)
    Cc ~ lnorm(expSd)
  })
}
