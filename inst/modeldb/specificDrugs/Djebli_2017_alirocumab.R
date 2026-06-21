Djebli_2017_alirocumab <- function() {
  description <- "Quasi-steady-state target-mediated drug disposition (TMDD-QSS) population PK model for alirocumab and total PCSK9 in healthy adults and adults with hypercholesterolemia (Djebli 2017, final model on expanded data set n=2870). Two-compartment disposition with first-order SC absorption (lag time and bioavailability), linear catabolic clearance from central, and PCSK9 binding / complex internalization described by QSS algebra; allometric weight scaling on CLL, Q, and Vc plus a statin-coadministration effect on CLL."
  reference   <- "Djebli N, Martinez JM, Lohan L, Khier S, Brunet A, Hurbin F, Fabre D. Target-Mediated Drug Disposition Population Pharmacokinetics Model of Alirocumab in Healthy Volunteers and Patients: Pooled Analysis of Randomized Phase I/II/III Studies. Clin Pharmacokinet. 2017;56(10):1155-1171. doi:10.1007/s40262-016-0505-1"
  vignette    <- "Djebli_2017_alirocumab"
  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling P_i = TVP * (WT/WT_med)^EXP with EXP = 0.75 for CLL and Q and EXP = 1 for Vc (Djebli 2017 Sect. 3.6, n=2870 expanded model). Reference body weight 85 kg is the mean baseline weight reported in Table 3 for the expanded data set; the paper labels the reference as WT_med (median) but only reports the mean, so the mean is used here as the closest reported summary.",
      source_name        = "WT"
    ),
    CONMED_STATIN = list(
      description        = "Concomitant statin (HMG-CoA reductase inhibitor) coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no statin coadministration)",
      notes              = "Multiplicative effect on linear clearance CLL: CLL = TVCLL * COV1^STATIN with COV1 = 1.27 (Djebli 2017 Table 4 expanded-model column and Sect. 3.6 equation). STATIN = 1 captures coadministration of any statin (rosuvastatin, atorvastatin, or simvastatin at any reported dose; low- and high-dose pooled). The covariate replaces the disease-state (DISST) effect on Vc retained in the smaller n=527 model because of strong collinearity between disease state and statin therapy in the n=527 cohort (Djebli 2017 Sect. 4).",
      source_name        = "CONMED_STATIN"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 2870L,
    n_observations   = 19595L,
    n_studies        = 13L,
    phases           = "Phase I, II, and III",
    age_mean         = "58.2 years (SD 11.7)",
    weight_mean      = "85.0 kg (SD 18.4)",
    bmi_mean         = "29.5 kg/m^2 (SD 5.40)",
    sex_female_pct   = 37.9,
    race_ethnicity   = c(Caucasian = 87.2, Black = 4.77, Asian = 5.02, Other = 3.03),
    disease_state    = "Pooled healthy adults (5.23%) and adults with hypercholesterolemia (94.77%); 28.5% familial hypercholesterolemia and 71.5% non-FH. Most patients on background lipid-lowering therapy (any statin 92.8%, ezetimibe 16.2%, fibrate 4.81%).",
    dose_range       = "Subcutaneous regimens 50-300 mg Q2W or Q4W (phase II/III) and 75 or 150 mg Q2W (ODYSSEY phase III). One phase-I study (NCT01026597) administered single intravenous doses of 0.3, 1, 3, 6, or 12 mg/kg.",
    regions          = "Multi-regional pool of 13 Sanofi/Regeneron trials including Japanese cohorts (NCT01448317 phase I, NCT01812707 phase II) and the ODYSSEY phase III programme (MONO, COMBO II, FH I, LONG TERM).",
    pcsk9_baseline   = "Mean baseline total PCSK9 = 9.14 nM (SD 6.84) in the expanded data set (Djebli 2017 Table 3).",
    notes            = "Baseline characteristics from Djebli 2017 Table 3 expanded-data-set column. The pooled dependent variable in the original NONMEM model combined total alirocumab and total PCSK9 concentrations on a common nM scale; the same propSd and addSd are applied to both outputs in this implementation."
  )

  ini({
    # ---- Structural PK parameters (Djebli 2017 Table 4, final model of expanded data set n=2870) ----
    # Reference body weight 85 kg (mean baseline weight in the expanded data set; Table 3).
    # All rates are reported per day in the paper; no unit conversion needed.
    lcl   <- log(0.176); label("Linear catabolic clearance CLL at reference covariates (L/day)")           # Table 4 expanded model column (CLL 0.176 L/day, 95% CI 0.152-0.200; CONMED_STATIN = 0, WT = 85 kg)
    lvc   <- log(4.67);  label("Central volume of distribution Vc at reference WT (L)")                    # Table 4 expanded model column (Vc 4.67 L, 95% CI 4.40-4.94; WT = 85 kg)
    lq    <- log(0.343); label("Inter-compartmental clearance Q at reference WT (L/day)")                  # Table 4 expanded model column (Q 0.343 L/day, 95% CI 0.321-0.365; WT = 85 kg)
    lvp   <- fixed(log(2.61)); label("Peripheral volume of distribution Vp (L; FIXED)")                    # Table 4 (Vp FIXED to 2.61 L in both the n=527 and expanded models; no allometric scaling)
    lka   <- log(0.307); label("First-order SC absorption rate constant Ka (1/day)")                       # Table 4 expanded model column (Ka 0.307 day^-1, 95% CI 0.286-0.329)
    lfdepot <- log(0.556); label("Subcutaneous bioavailability F1 (fraction)")                             # Table 4 expanded model column (F1 0.556, 95% CI 0.531-0.582)
    ltlag   <- log(0.0535); label("SC absorption lag time LAG (days)")                                    # Table 4 expanded model column (LAG 0.0535 days, 95% CI 0.0533-0.0537)

    # ---- TMDD-QSS parameters (Djebli 2017 Table 4 expanded model) ----
    lkdeg <- log(1.10);    label("First-order PCSK9 degradation rate constant Kdeg (1/day)")              # Table 4 expanded model column (Kdeg 1.10 day^-1, 95% CI 1.06-1.14)
    lkint <- log(0.112);   label("First-order alirocumab-PCSK9 complex internalization rate Kint (1/day)") # Table 4 expanded model column (Kint 0.112 day^-1, 95% CI 0.109-0.115)
    # The paper parameterizes the QSS approximation by Kon (FIXED to 559 nM^-1 day^-1; Table 4) together
    # with Kdeg, Kint, and a measured in-vitro dissociation constant K_D ~ 0.58 nM (Methods Sect. 2.3.2,
    # Fig. 3 caption and Sect. 3.2). The standard QSS algebra (Gibiansky et al. 2008) requires
    # K_SS = (K_off + K_int) / K_on. With K_off = K_on * K_D = 559 * 0.58 = 324 day^-1 and K_int = 0.112 day^-1,
    # K_off >> K_int so K_SS ~ K_D = 0.58 nM (the K_int / K_on contribution is ~0.0002 nM, negligible).
    # K_SS is therefore fixed at the published in-vitro K_D; K_on is preserved as a documented constant
    # in the model body for reproducibility of the binding-rate parameterization.
    lkss  <- fixed(log(0.58)); label("Quasi-steady-state dissociation constant K_SS for alirocumab-PCSK9 (nM; FIXED)") # Methods Sect. 2.3.2 (in-vitro K_D 0.58 nM) + Fig. 3 caption (K_SS = (K_off + K_int)/K_on)
    lkon  <- fixed(log(559));  label("Binding association rate constant K_on (nM^-1 day^-1; FIXED; informational)")    # Table 4 (K_on FIXED to 559 nM^-1 day^-1 across all models)

    # Baseline total PCSK9 concentration. The paper does not estimate baseline PCSK9 as a population
    # parameter; the initial condition for the target compartment is set to the reported population
    # mean baseline total PCSK9 in the expanded data set (Table 3). No IIV is encoded because the paper
    # does not report one for this quantity.
    lrbase <- fixed(log(9.14)); label("Baseline total PCSK9 R_base (nM; FIXED at population mean)")        # Table 3 expanded data set (baseline total PCSK9 mean 9.14 nM, SD 6.84)

    # ---- Allometric exponents (Djebli 2017 Sect. 3.6 equation; theory-based, FIXED) ----
    # P_i = TVP * (WT / WT_med)^EXP with EXP = 0.75 for CLL and Q and EXP = 1 for Vc. The paper
    # describes these as "theory-based allometric scaling" (i.e., FIXED at the canonical values).
    allo_cl  <- fixed(0.75); label("Allometric exponent of (WT/85)^EXP on CLL (FIXED)")                    # Sect. 3.6 (theory-based allometric scaling on CLL)
    allo_q   <- fixed(0.75); label("Allometric exponent of (WT/85)^EXP on Q (FIXED)")                      # Sect. 3.6 (theory-based allometric scaling on Q)
    allo_vc  <- fixed(1.0);  label("Allometric exponent of (WT/85)^EXP on Vc (FIXED)")                     # Sect. 3.6 (theory-based allometric scaling on Vc)

    # ---- Covariate effect: statin coadministration on CLL ----
    # CLL = TVCLL * COV1^STATIN; COV1 = 1.27 with statin coadministration (Table 4 final expanded
    # model). Implemented as a power of CONMED_STATIN: CLL = TVCLL * 1.27^CONMED_STATIN.
    e_conmed_statin_cl <- 1.27; label("Multiplicative effect of CONMED_STATIN on CLL (unitless; CLL fold-change with statin)") # Table 4 expanded model column (COV1 STATIN on CLL 1.27, 95% CI 1.10-1.45)

    # ---- Inter-individual variability (Djebli 2017 Table 4 expanded model) ----
    # Paper reports exponential IIV with CV%; estimate column gives the omega^2 value with
    # CV% = sqrt(omega^2) * 100 (e.g., 0.193 -> sqrt(0.193) ~ 43.9% matches Table 4). LAG, K_on,
    # and V_p had no IIV (Table 4 NA entries). No covariance terms were retained (Methods Sect. 2.3.2
    # and Sect. 3.2: "No combination (among 56 tested) of estimated covariance between the seven
    # estimated inter-individual variabilities ... improved the quality of the model").
    etalcl    ~ 0.193  # Table 4 expanded model (omega^2 = 0.193, CV 43.9%)
    etalkint  ~ 0.0620 # Table 4 expanded model (omega^2 = 0.0620, CV 24.9%)
    etalkdeg  ~ 0.219  # Table 4 expanded model (omega^2 = 0.219, CV 46.8%)
    etalq     ~ 0.426  # Table 4 expanded model (omega^2 = 0.426, CV 65.3%)
    etalvc    ~ 0.0515 # Table 4 expanded model (omega^2 = 0.0515, CV 22.7%)
    etalka    ~ 0.723  # Table 4 expanded model (omega^2 = 0.723, CV 85.0%)
    etalfdepot ~ 0.433 # Table 4 expanded model (omega^2 = 0.433 for F1, CV 65.8%)

    # ---- Residual error (Djebli 2017 Table 4 expanded model; pooled DV on the nM scale) ----
    # Paper combines total alirocumab and total PCSK9 into a single NONMEM dependent variable on
    # the nM scale, with a combined additive + proportional error model. nlmixr2 requires distinct
    # residual-parameter names per endpoint, so the same paper-reported values are applied to both
    # observation symbols Cc (total alirocumab) and Ct (total PCSK9). The additive nM-scale error
    # is converted to the mg/L-scale observation per output using the analyte molecular weight:
    #   addSd_mgL = 10.4 nM * MW (g/mol) / 1e6 = mg/L
    # Cc:  10.4 * 144100 / 1e6 = 1.499 mg/L   (alirocumab MW 144100 g/mol)
    # Ct:  10.4 *  74000 / 1e6 = 0.7696 mg/L  (mature secreted PCSK9 MW ~74 kDa; see Assumptions)
    propSd    <- 0.152;  label("Proportional residual error for total alirocumab Cc (fraction)")           # Table 4 expanded model (pooled DV proportional 0.152, CV 15.2%, 95% CI 0.150-0.154)
    addSd     <- 1.499;  label("Additive residual error for total alirocumab Cc (mg/L; = 10.4 nM * MW_alirocumab)") # Table 4 expanded model (10.4 nM additive converted via MW 144100 g/mol)
    propSd_Ct <- 0.152;  label("Proportional residual error for total PCSK9 Ct (fraction)")                # Table 4 expanded model (pooled DV proportional, same value as for Cc)
    addSd_Ct  <- 0.7696; label("Additive residual error for total PCSK9 Ct (mg/L; = 10.4 nM * MW_PCSK9)")  # Table 4 expanded model (10.4 nM additive converted via PCSK9 mature-form MW 74000 g/mol)
  })

  model({
    # ---- Physical constants (top of model() per template convention) ----
    # Alirocumab is a fully human IgG1 monoclonal antibody. The reference molecular weight
    # 144,100 g/mol is the value reported in the alirocumab USPI / Praluent (Sanofi/Regeneron).
    # It is used to convert mg doses into nmol so the central compartment carries amounts on
    # the same molar scale as the QSS binding algebra (nM = nmol/L).
    mw_alirocumab   <- 144100               # g/mol (alirocumab; Praluent prescribing information)
    nmol_per_mg     <- 1e6 / mw_alirocumab  # 1 mg = 1e6 / MW nmol ~ 6.939 nmol per mg

    # ---- Individual parameters (covariate model: allometric on CLL/Q/Vc, STATIN on CLL) ----
    cl   <- exp(lcl + etalcl) * (WT / 85)^allo_cl * e_conmed_statin_cl^CONMED_STATIN
    vc   <- exp(lvc + etalvc) * (WT / 85)^allo_vc
    q    <- exp(lq + etalq)   * (WT / 85)^allo_q
    vp   <- exp(lvp)
    ka   <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag)
    kdeg <- exp(lkdeg + etalkdeg)
    kint <- exp(lkint + etalkint)
    kss  <- exp(lkss)
    rbase <- exp(lrbase)
    kon  <- exp(lkon)  # documented constant; not used in QSS algebra

    # Zero-order PCSK9 synthesis (steady-state balance with no drug -> tfree = rbase, complex = 0)
    ksyn <- kdeg * rbase

    # ---- QSS algebra (Gibiansky et al. 2008) ----
    # central tracks the TOTAL alirocumab amount (free + complex) in nmol.
    # total_target tracks TOTAL PCSK9 concentration (free + complex) in central (nM).
    # ctot is total drug concentration in nM (nmol/L); ttot is total target concentration in nM.
    ctot <- central / vc
    ttot <- total_target

    # Numerically stable discriminant: (Ctot + Ttot + Kss)^2 - 4 Ctot Ttot
    #                                = (Ctot - Ttot)^2 + 2 (Ctot + Ttot) Kss + Kss^2
    qss_disc <- (ctot - ttot)^2 + 2 * (ctot + ttot) * kss + kss^2
    complex  <- ((ctot + ttot + kss) - sqrt(qss_disc)) / 2
    cfree    <- ctot - complex
    tfree    <- ttot - complex

    # ---- ODE system (Djebli 2017 Sect. 2.3.2 / Sect. 3.2 with QSS) ----
    # depot: SC depot with first-order Ka and bioavailability F1 + LAG (subcutaneous administration).
    # central: total alirocumab amount in central compartment (nmol). Free drug eliminates linearly
    #          (CL * cfree), distributes to peripheral (Q * cfree out, Q * Cp in), and the complex
    #          eliminates via Kint (Kint * complex * Vc, converting nM/day back to nmol/day).
    # peripheral1: non-binding tissue distribution of free alirocumab; carries amount in nmol.
    # total_target: total PCSK9 concentration in central (nM). Synthesised at constant rate ksyn,
    #               free PCSK9 degrades at kdeg * tfree, and complex internalizes at kint * complex.
    d/dt(depot)         <- -ka * depot
    d/dt(central)       <-  ka * depot - cl * cfree - q * cfree + q * (peripheral1 / vp) - kint * complex * vc
    d/dt(peripheral1)   <-                            q * cfree - q * (peripheral1 / vp)
    d/dt(total_target)  <-  ksyn - kdeg * tfree - kint * complex
    total_target(0)     <-  rbase

    # ---- SC absorption: bioavailability multiplier and lag time on the depot ----
    # SC doses (the marketed regimen) deposit mg of drug in `depot`; convert to nmol via the
    # nmol_per_mg factor and apply F1 + LAG. IV doses (phase-I NCT01026597) bypass the depot via
    # the event table (cmt = "central"); to put them on the same nmol scale, the same nmol_per_mg
    # factor is also applied to central (bioavailability is 1 for IV but the unit conversion is the
    # same). The vignette uses `cmt = "depot"` (SC) or `cmt = "central"` (IV) on dose events.
    f(depot)    <- fdepot * nmol_per_mg
    alag(depot) <- tlag
    f(central)  <- nmol_per_mg

    # ---- Observations (paper applies same residual to total alirocumab and total PCSK9 on nM scale) ----
    # Internal nM quantities ctot (total alirocumab) and ttot (total PCSK9) are converted to mg/L
    # using the analyte molecular weight for the user-facing observations Cc and Ct. The additive
    # residual is converted with the same scale (see ini()); the proportional residual is unitless.
    # Observation rows in event tables use cmt = "central" for Cc and cmt = "total_target" for Ct;
    # rxode2 maps each endpoint to the right output column based on the cmt of the observation row.
    mw_pcsk9 <- 74000  # PCSK9 mature secreted form, g/mol (UniProt Q8NBP7; see Assumptions)
    Cc <- ctot * mw_alirocumab / 1e6   # total alirocumab concentration (mg/L)
    Ct <- ttot * mw_pcsk9     / 1e6   # total PCSK9 concentration       (mg/L)
    Cc ~ add(addSd) + prop(propSd)
    Ct ~ add(addSd_Ct) + prop(propSd_Ct)
  })
}
