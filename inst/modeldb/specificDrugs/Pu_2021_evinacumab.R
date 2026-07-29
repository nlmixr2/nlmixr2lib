Pu_2021_evinacumab <- function() {
  description <- "Population PK/PD model for evinacumab in healthy volunteers and adults / pediatric patients with homozygous familial hypercholesterolemia (Pu 2021): two-compartment PK with first-order SC absorption (with lag time) and parallel linear plus Michaelis-Menten elimination from the central compartment, linked to a Type 1 indirect-response model for low-density lipoprotein cholesterol (LDL-C) where evinacumab inhibits LDL-C production."
  reference   <- "Pu X, Sale M, Yang F, Zhang Y, Davis JD, Al-Huniti N. Population pharmacokinetics and exposure-response modeling for evinacumab in homozygous familial hypercholesterolemia. CPT Pharmacometrics Syst Pharmacol. 2021;10(11):1412-1421. doi:10.1002/psp4.12711"
  vignette    <- "Pu_2021_evinacumab"
  paper_specific_compartments <- c("LDL")

  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Allometric power scaling on linear CL with reference 74.1 kg (PK reference; Table 1) and on Vc, Vp, and Q with the same 74.1 kg reference (the Pu 2021 K23, K32 rate constants are weight-invariant, so Vp = Q / K32 inherits the Vc allometric exponent). Power effect on Imax with reference 71 kg (PD reference; Table 2 footnote, equal to the median weight of the patients with HoFH in the PK/PD analysis set).",
      source_name        = "WGTBL"
    ),
    DIS_HOFH = list(
      description        = "Homozygous familial hypercholesterolemia disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer reference)",
      notes              = "Time-fixed per subject. Multiplicative log-effect on Vmax: HoFH patients exhibit ~25% lower target-mediated Vmax than the HV reference (Pu 2021 Table 1: theta = -0.289, exp(-0.289) ~= 0.749). Renamed from source column DISTYPN to canonical DIS_HOFH per inst/references/covariate-columns.md.",
      source_name        = "DISTYPN"
    ),
    ANGPTL3 = list(
      description        = "Baseline total serum angiopoietin-like protein 3 concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline only (per-subject time-fixed). Power-form effect on Vmax with reference 0.08 mg/L (Pu 2021 typical-patient median, paper text). Higher baseline target predicts a faster saturable elimination; biologically consistent with evinacumab being co-cleared along with bound ANGPTL3. Bioanalytical assay detects free + target-bound ANGPTL3 after acid pretreatment of serum (Pu 2021 Methods, LLOQ 0.0195 mg/L). Renamed from source column ANGBL to canonical ANGPTL3 per inst/references/covariate-columns.md.",
      source_name        = "ANGBL"
    ),
    LDLC = list(
      description        = "Baseline serum low-density lipoprotein cholesterol concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline only (per-subject time-fixed). Used in two roles in this model: (1) initialises the indirect-response LDL state (LDL(0) <- LDLC), and (2) drives a power-form effect on IC50 with reference 211 mg/dL (Pu 2021 Table 2 footnote). Patients with higher baseline LDL-C are more sensitive to evinacumab (smaller IC50). Renamed from source column LDLBL to canonical LDLC per inst/references/covariate-columns.md.",
      source_name        = "LDLBL"
    ),
    RACE_WHITE = list(
      description        = "White (Caucasian) race indicator (1 = White, 0 = non-White / Other)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (White; reference for the typical-value PD parameter Imax = 0.7435)",
      notes              = "Multiplicative effect on Imax: non-White subjects have a lower maximum drug-induced inhibitory effect on LDL-C production than the White reference (Pu 2021 Table 2: theta = -0.191, exp(-0.191) ~= 0.83 multiplier when non-White). Source paper (and supplementary control stream) used a binary White-vs-other dichotomy; the typical-value reference is the White (RACE_WHITE = 1) subgroup, which is the inverse of the Lin 2024 use of the same canonical column. Renamed from source column RAC1 (where RAC1 = 1 means non-White / others) to canonical RACE_WHITE per inst/references/covariate-columns.md; the model uses (1 - RACE_WHITE) inline so the parameter values reproduce the published equations as-is.",
      source_name        = "RAC1"
    )
  )

  population <- list(
    n_subjects                 = 278L,
    n_subjects_pk              = 278L,
    n_subjects_pkpd            = 95L,
    n_studies                  = 6L,
    n_studies_pkpd             = 3L,
    phases                     = "Pooled phase I, II, and III (PK); pooled phase II + two phase III (PK/PD)",
    age_range                  = "12 years and older (FDA-approved population includes adolescents and adults; the analysis set spans the phase I HV cohorts and the HoFH phase II/III studies)",
    weight_range               = "42.4-152 kg (Pu 2021 Figure 2 covariate ranges)",
    weight_median              = "74.1 kg (PK analysis set median; Pu 2021 Table 1 typical-patient definition); 71 kg (PK/PD analysis set median; Pu 2021 Table 2 footnote)",
    sex_female_pct             = NA_real_,
    race_ethnicity             = "Pooled across HoFH phase II/III + HV phase I cohorts; race coded as White (RAC1 = 0) vs Other (RAC1 = 1) per the source NM-TRAN dataset.",
    disease_state              = "Healthy volunteers (n = 183) and patients with homozygous familial hypercholesterolemia (n = 95).",
    dose_range                 = "IV: 5-20 mg/kg single dose or repeated weekly, every 4 weeks, or every 12 weeks. SC: 75-450 mg single dose or repeated weekly or every 2 weeks. Approved adult/adolescent (>= 12 y) HoFH dose is 15 mg/kg IV every 4 weeks.",
    regions                    = "Multi-regional pool of six Regeneron-sponsored clinical trials (NCT01749878, NCT03146416, NCT02107872, NCT02265952, NCT03399786) and one phase I HV study referenced in the paper.",
    angptl3_baseline_median    = "0.08 mg/L (paper text; comparable between HoFH and HV cohorts)",
    ldlc_baseline_typical      = "211 mg/dL (PK/PD typical patient)",
    samples_pk                 = "Total evinacumab serum concentrations measured by validated ELISA (LLOQ 0.078 mg/L; assay detects free evinacumab + evinacumab bound to one or two ANGPTL3 molecules).",
    samples_pkpd               = "1194 LDL-C concentration timepoints from 95 HoFH patients (after removing 9 outliers with CWRES > 5) used for the PK/PD model.",
    notes                      = "PK and PK/PD reference covariates differ slightly: PK reference is 74.1 kg with 0.08 mg/L baseline ANGPTL3 (Table 1 caption: 'CL, L/day for 74.1 kg subject'); PK/PD reference is 71 kg with 211 mg/dL baseline LDL-C (Table 2 footnote). The model file preserves both reference values to reproduce the published parameter equations exactly. Disease-state effect on Vmax was estimated from the pooled HV+HoFH PK dataset; the PK/PD layer was estimated only from HoFH patients (no DIS_HOFH covariate appears at the PD level). Race effect on Imax was tested only at the PK/PD level. The bioanalytical assay for ANGPTL3 included acid pretreatment of serum and detected both free and bound ANGPTL3."
  )

  ini({
    # ----- PK structural parameters (Pu 2021 Table 1 final-model typical values) -----
    # Reference covariates: WT = 74.1 kg, ANGPTL3 = 0.08 mg/L, DIS_HOFH = 0 (HV reference).
    # Units: time = day, dose = mg, concentration = mg/L (= ug/mL).
    lka     <- log(0.181);    label("First-order SC absorption rate Ka (1/day)")                                          # Pu 2021 Table 1 final-model Ka
    lcl     <- log(0.0955);   label("Linear clearance CL at 74.1 kg HV reference (L/day)")                               # Pu 2021 Table 1 final-model CL
    lvc     <- log(2.56);     label("Central volume of distribution Vc (V2 in source) at 74.1 kg reference (L)")          # Pu 2021 Table 1 final-model V2
    lq      <- log(0.109 * 2.56); label("Inter-compartmental clearance Q at 74.1 kg reference (L/day; = K23 * V2)")       # Pu 2021 Table 1 final-model K23 * V2 = 0.109 * 2.56
    lvp     <- log(0.109 * 2.56 / 0.124); label("Peripheral volume of distribution Vp (V3 in source) at 74.1 kg reference (L; = K23 * V2 / K32)")  # Pu 2021 Table 1 final-model K23 * V2 / K32 = 0.109 * 2.56 / 0.124
    lvmax   <- log(3.16);     label("Maximum target-mediated Michaelis-Menten elimination rate Vmax at HV / typical ANGPTL3 reference (mg/day)")   # Pu 2021 Table 1 final-model Vmax
    lkm     <- log(1.02);     label("Michaelis-Menten constant Km (mg/L)")                                                # Pu 2021 Table 1 final-model Km
    ltlag    <- log(0.168);    label("SC absorption lag time (day)")                                                       # Pu 2021 Table 1 final-model Alag
    lfdepot <- log(0.714);    label("Log SC bioavailability (F = 0.714, fraction)")                                        # Pu 2021 Table 1 final-model F = 0.714

    # ----- PK covariate effects (Pu 2021 Table 1 covariate block) -----
    e_wt_cl       <- fixed(0.75);  label("Allometric power exponent: (WT/74.1)^e_wt_cl on linear CL (unitless)")            # Pu 2021 Table 1 weight ~ linear clearance, fixed
    e_wt_vc_vp    <- 0.875;        label("Allometric power exponent: (WT/74.1)^e_wt_vc_vp on Vc, shared with Q and Vp via the weight-invariance of K23/K32 (unitless)")  # Pu 2021 Table 1 weight ~ central V (estimated)
    e_angptl3_vmax <- 0.405;       label("Power exponent of (ANGPTL3 / 0.08 mg/L) on Vmax (unitless)")                      # Pu 2021 Table 1 ANGPTL3 ~ Vmax
    e_dis_hofh_vmax <- -0.289;     label("Exponential effect of DIS_HOFH on Vmax: Vmax * exp(theta * DIS_HOFH); HoFH gives ~25% reduction (unitless)")  # Pu 2021 Table 1 disease state ~ Vmax

    # ----- PK inter-individual variability (Pu 2021 Table 1 OMEGA block) -----
    # Paper Table 1 reports the lower-triangular OMEGA in the SD/correlation form:
    #   sigma(eta_CL) = 0.355, sigma(eta_V) = 0.213, corr(eta_CL, eta_V) = 0.213.
    # Convert to NONMEM-style variance/covariance:
    #   var(eta_CL) = 0.355^2  = 0.126025
    #   var(eta_V)  = 0.213^2  = 0.045369
    #   cov(CL, V)  = 0.213 * 0.355 * 0.213 = 0.016108
    # Approximate CV%: CV = sqrt(exp(sigma^2) - 1) -> 36% (CL) and 21% (Vc), matching the paper text.
    etalcl + etalvc ~ c(0.126025, 0.016108, 0.045369)  # Pu 2021 Table 1 sigma(CL)=0.355, sigma(V)=0.213, corr(CL,V)=0.213
    etalka          ~ 0.470596                          # Pu 2021 Table 1 sigma(Ka) = 0.686 -> var = 0.686^2
    etaltlag         ~ 1.4161                            # Pu 2021 Table 1 sigma(Alag1) = 1.19 -> var = 1.19^2

    # ----- PK residual error (Pu 2021 Table 1 final-model RUV) -----
    addSd  <- 0.303; label("Additive PK residual error on Cc (mg/L)")                            # Pu 2021 Table 1 sigma additive
    propSd <- 0.189; label("Proportional PK residual error on Cc (fraction)")                     # Pu 2021 Table 1 sigma proportional, CV%

    # ----- PD structural parameters (Pu 2021 Table 2 final-model typical values) -----
    # Reference covariates: WT = 71 kg, baseline LDL-C = 211 mg/dL, RACE_WHITE = 1 (White reference).
    # LDL-C concentration units: mg/dL (paper Table 2: sigma additive = 17.97 mg/dL).
    lkin   <- log(38.99); label("Zero-order LDL-C production rate constant Kin (mg/dL/day)")    # Pu 2021 Table 2 final-model Kin
    limax  <- log(0.7435); label("Log of maximum drug-induced inhibitory effect Imax at White / 71 kg reference (unitless, < 1)")  # Pu 2021 Table 2 final-model Imax = 0.743459
    lic50  <- log(57.4);  label("Evinacumab concentration giving 50% of Imax at 211 mg/dL baseline-LDL reference (IC50, mg/L)")     # Pu 2021 Table 2 final-model IC50

    # ----- PD covariate effects (Pu 2021 Table 2 covariate block; supplement final NONMEM expressions) -----
    e_wt_imax       <- -0.27;  label("Power exponent of (WT/71 kg) on Imax (unitless)")                                       # Pu 2021 Table 2 Imax ~ weight (power)
    e_nonwhite_imax <- -0.191; label("Exponential effect of (1 - RACE_WHITE) on Imax: exp(theta * (1 - RACE_WHITE)) (unitless)")  # Pu 2021 Table 2 Imax ~ race (footnote: theta = -0.191248); model uses (1 - RACE_WHITE) so non-White subjects (RACE_WHITE = 0) get the exp(-0.191) ~ 0.83 multiplier and the White-reference TVImax = 0.7435 is preserved unchanged
    e_ldlc_ic50     <- -1.17;  label("Power exponent of (baseline LDLC / 211 mg/dL) on IC50 (unitless)")                       # Pu 2021 Table 2 IC50 ~ baseline LDL (power)

    # ----- PD inter-individual variability (Pu 2021 Table 2 OMEGA block) -----
    # Paper Table 2 reports the diagonal OMEGA values directly as variances in the
    # NONMEM internal (log) scale (matched to the supplement $OMEGA initials of 0.602 / 0.064);
    # see the Pu 2021 final $OMEGA block for IC50 (3.11) and Kin (0.47).
    etalic50 ~ 3.11   # Pu 2021 Table 2 sigma, eta(IC50)
    etalkin  ~ 0.47   # Pu 2021 Table 2 sigma, eta(Kin)

    # ----- PD residual error (Pu 2021 Table 2 final-model RUV) -----
    addSd_LDL  <- 17.97; label("Additive PD residual error on LDL-C (mg/dL)")              # Pu 2021 Table 2 sigma additive
    propSd_LDL <- 0.18;  label("Proportional PD residual error on LDL-C (fraction)")        # Pu 2021 Table 2 sigma proportional
  })
  model({
    # ----- Individual PK parameters (Pu 2021 Table 1 reference WT = 74.1 kg) -----
    # K23 and K32 are weight-invariant rate constants in the source paper, so Q = K23 * Vc
    # and Vp = Q / K32 inherit the Vc allometric exponent. The paper does not report
    # IIV on Q, Vp, Vmax, Km, or F, so these parameters are typical-value-only in this layer.
    ka    <- exp(lka + etalka)
    cl    <- exp(lcl + etalcl) * (WT / 74.1)^e_wt_cl
    vc    <- exp(lvc + etalvc) * (WT / 74.1)^e_wt_vc_vp
    vp    <- exp(lvp)          * (WT / 74.1)^e_wt_vc_vp
    q     <- exp(lq)           * (WT / 74.1)^e_wt_vc_vp
    km    <- exp(lkm)
    vmax  <- exp(lvmax) * (ANGPTL3 / 0.08)^e_angptl3_vmax * exp(e_dis_hofh_vmax * DIS_HOFH)
    alag1 <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot)

    # ----- PK ODE system (Pu 2021 Figure 1 schematic) -----
    # Two-compartment evinacumab disposition with parallel linear (CL) and Michaelis-Menten
    # elimination from the central compartment. The MM term is Vmax * Cc / (Km + Cc) so that
    # Vmax has units of mg/day (mass per time) and matches Pu 2021 Table 1.
    Cc <- central / vc
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q  / vc) * central +
                          (q  / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # SC absorption lag and bioavailability apply to the depot (oral / SC).
    # IV doses bypass the depot via cmt = central on the dose event.
    f(depot)    <- fdepot
    alag(depot) <- alag1

    # ----- Individual PD parameters (Pu 2021 Table 2 reference WT = 71 kg, LDLC = 211 mg/dL,
    #       RACE_WHITE = 1 (White)) -----
    kin   <- exp(lkin + etalkin)
    imax  <- exp(limax) * (WT / 71)^e_wt_imax * exp(e_nonwhite_imax * (1 - RACE_WHITE))
    ic50  <- exp(lic50 + etalic50) * (LDLC / 211)^e_ldlc_ic50

    # Kout enforces baseline LDL = LDLC at steady state (no drug): Kin = Kout * LDLC.
    # Per-subject baseline LDLC is the time-fixed baseline value carried in the data.
    kout <- kin / LDLC

    # ----- Indirect-response PD ODE for LDL-C (Pu 2021 Figure 1, Type 1 indirect-response) -----
    # Drug inhibits LDL-C production via Imax * Cc / (IC50 + Cc); LDL-C eliminates with rate Kout.
    LDL(0) <- LDLC
    d/dt(LDL) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * LDL

    # ----- Observation and error models -----
    # Cc in mg/L (= ug/mL); LDL in mg/dL.
    Cc  ~ add(addSd)     + prop(propSd)
    LDL ~ add(addSd_LDL) + prop(propSd_LDL)
  })
}
