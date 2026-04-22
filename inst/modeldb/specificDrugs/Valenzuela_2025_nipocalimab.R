Valenzuela_2025_nipocalimab <- function() {
  description <- "Integrated PK/RO/IgG/MG-ADL QSS TMDD model for nipocalimab in healthy adults and generalized myasthenia gravis (Valenzuela 2025)"
  reference <- "Valenzuela B, Neyens M, Zhu Y, Ramchandren S, Dosne AG, Leu JH, Faelens R, Ling LE, Perez-Ruixo JJ. Nipocalimab Dose Selection in Generalized Myasthenia Gravis. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2074-2085. doi:10.1002/psp4.70109"
  vignette <- "Valenzuela_2025_nipocalimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline); allometric scaling (WT/75)^exponent on CL, Q, Vc, Vp with exponents fixed at 0.75 (CL, Q) and 1 (Vc, Vp).",
      source_name        = "WT"
    ),
    ELISA = list(
      description        = "Bioanalytical assay indicator: 1 = ELISA (LLOQ 0.150 ug/mL), 0 = ECLIA (LLOQ 0.010 ug/mL)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECLIA)",
      notes              = "Switches the additive PK residual magnitude. Assay is study-fixed: ELISA for studies MOM-M281-001, MOM-M281-007, MOM-M281-010; ECLIA for studies EDI1001, EDI1002, and MOM-M281-004 (Vivacity-MG).",
      source_name        = "ELISA"
    ),
    PHASE1 = list(
      description        = "Study-phase indicator: 1 = Phase 1 study (healthy participants), 0 = Phase 2 Vivacity-MG study in gMG",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 2 MOM-M281-004 / Vivacity-MG)",
      notes              = "Switches the proportional PK residual magnitude (0.0834 for Phase 1 vs 0.367 for Phase 2).",
      source_name        = "PHASE1"
    ),
    STUDY_M281_004 = list(
      description        = "Vivacity-MG study indicator: 1 = subject enrolled in MOM-M281-004 (NCT03772587), 0 = other study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Vivacity-MG studies)",
      notes              = "Scales IgG baseline by FRIgG0_M281_004 = 0.777 for MOM-M281-004 participants; IgG baseline is equal to IgG0 (typical 11.4 g/L) for other studies.",
      source_name        = "M281_004"
    ),
    MGADL = list(
      description        = "Baseline Myasthenia Gravis Activities of Daily Living score (0-24; higher = more severe)",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline); power-form effect (MGADL/7)^E_IgG on SIgG and (MGADL/7)^E_ADL on IDec_placebo. Set MGADL = 0 for healthy participants so that MG-ADL response predictions collapse to 0. Vivacity-MG participants ranged 6-24 points (mean 9.8).",
      source_name        = "MGADL"
    )
  )

  population <- list(
    n_subjects     = 228L,
    n_studies      = 6L,
    age_range      = "18-83 years",
    age_median     = "41.3 years (mean; SD 16.2)",
    weight_range   = "45-188 kg",
    weight_median  = "75.9 kg (mean; SD 19.2)",
    sex_female_pct = 52.6,
    race_ethnicity = c(White = 68.0, Black = 3.51, Asian = 25.0, Other = 3.51),
    disease_state  = "Pooled healthy adult participants (160) and adult participants with generalized myasthenia gravis (68, from Vivacity-MG / MOM-M281-004)",
    dose_range     = "0.3-60 mg/kg single IV dose; 5-60 mg/kg IV Q4W or Q2W (gMG cohort)",
    regions        = "Pooled phase 1 (healthy) and phase 2 (gMG) studies; populations include Japanese (study MOM-M281-010) and Chinese (study EDI1002) healthy participants",
    n_subjects_healthy = 160L,
    n_subjects_gmg     = 68L,
    mgadl_range_gmg    = "6-24 points (mean 9.8; Vivacity-MG baseline)",
    notes              = "Demographics from Valenzuela 2025 Table 2. Trials: NCT02828046 (MOM-M281-001), sequential-dose phase 1 (MOM-M281-007), dose-ranging Japanese phase 1 (MOM-M281-010), NCT04848558 (EDI1001), NCT05151692 (EDI1002), NCT03772587 (MOM-M281-004 / Vivacity-MG)."
  )

  ini({
    # ---- PK disposition (Table 3; reference body weight 75 kg) ----
    lcl     <- log(0.655);  label("Linear serum clearance for a 75 kg adult (L/day)")                     # Table 3: CL = 0.655 L/d
    lvc     <- log(3.23);   label("Central volume of distribution for a 75 kg adult (L)")                 # Table 3: Vc = 3.23 L
    lq      <- log(0.250);  label("Intercompartmental clearance for a 75 kg adult (L/day)")               # Table 3: Q = 0.250 L/d
    lvp     <- log(0.622);  label("Peripheral volume of distribution for a 75 kg adult (L)")              # Table 3: Vp = 0.622 L
    allo_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")                           # Eq. 3; fixed per paper footnote
    allo_v  <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)")                          # Eq. 3; fixed per paper footnote

    # ---- FcRn target turnover (Table 3) ----
    # FcRn0 reported in nmol/L by the paper; converted to ug/mL inside model()
    # using nipocalimab MW 142 kDa so all drug-FcRn QSS arithmetic is carried
    # out in a single ug/mL scale consistent with dosing in mg and volumes in L.
    lFcRn0  <- log(143);    label("Total FcRn concentration at baseline (nmol/L; converted to ug/mL inside model)") # Table 3: FcRn0 = 143 nmol/L (= 20.3 ug/mL)
    FRmax   <- 0.947;       label("Maximal fraction of FcRn available for nipocalimab binding (fraction)") # Table 3: FRmax = 0.947
    lkss    <- log(6.05);   label("Quasi-steady-state dissociation constant (ug/mL)")                      # Table 3: Kss = 6.05 ug/mL
    lkint   <- log(62.4);   label("Internalization rate of nipocalimab-FcRn complex (1/day)")             # Table 3: kint = 62.4 1/d
    lkdeg   <- log(1.3);    label("Degradation rate of free FcRn (1/day)")                                # Table 3: kdeg = 1.3 1/d (fixed per paper)

    # ---- IgG indirect-response (Table 3) ----
    lIgG0          <- log(11.4);   label("IgG baseline concentration in non-Vivacity-MG studies (g/L)")    # Table 3: IgG0 = 11.4 g/L
    FRIgG0_M281_004 <- 0.777;      label("Fractional scaling of IgG baseline for study MOM-M281-004 (fraction)") # Table 3: FR_IgG0,M281-004 = 0.777
    lkdeg_IgG      <- log(0.217);  label("IgG degradation rate in absence of FcRn recycling (1/day)")     # Table 3: kdeg,IgG = 0.217 1/d
    lIgK           <- log(5.08);   label("FcRn-mediated IgG recycling parameter IgK (unitless)")          # Table 3: IgK = 5.08 (treated on log scale for exponential IIV)

    # ---- MG-ADL PD (Table 3) ----
    lke0        <- log(0.414); label("IgG effect-compartment rate constant for MG-ADL (1/day)")           # Table 3: ke0 = 0.414 1/d
    SIgG        <- -0.216;     label("Slope between MG-ADL change and IgG reduction (points per 10% IgG reduction)") # Table 3: S_IgG = -0.216 pts/10%reduction; note the 10x unit conversion in model()
    EIgG        <- 0.871;      label("Exponent of MG-ADL baseline effect on S_IgG (unitless)")            # Table 3: E_IgG = 0.871
    IDecplacebo <- -1.08;      label("Typical initial placebo decrease in MG-ADL after treatment start (points)") # Table 3: IDec_placebo = -1.08 pts
    EADL        <- 1.23;       label("Exponent of MG-ADL baseline effect on IDec_placebo (unitless)")      # Table 3: E_ADL = 1.23
    Splacebo    <- -0.0594;    label("Slope of placebo effect over time on MG-ADL (points/week)")          # Table 3: S_placebo = -0.0594 pts/week

    # ---- IIV (exponential on PK / FcRn / IgG params; additive on MG-ADL params) ----
    # omega^2 = log(CV^2 + 1) for log-normal parameters:
    #   CV 24.1% -> 0.05644   (CL)
    #   CV 15.0% -> 0.02225   (Vc)
    #   CV 25.0% -> 0.06062   (FcRn0)
    #   CV 21.9% -> 0.04686   (IgG0)
    #   CV 16.9% -> 0.02816   (kdeg_IgG)
    #   CV 26.3% -> 0.06690   (IgK)
    etalcl       ~ 0.05644    # Table 3 IIV CV 24.1%
    etalvc       ~ 0.02225    # Table 3 IIV CV 15.0%
    etalFcRn0    ~ 0.06062    # Table 3 IIV CV 25.0%
    etalIgG0     ~ 0.04686    # Table 3 IIV CV 21.9%
    etalkdeg_IgG ~ 0.02816    # Table 3 IIV CV 16.9%
    etalIgK      ~ 0.06690    # Table 3 IIV CV 26.3%

    # MG-ADL additive IIV values reported in Table 3 (1.89 and 3.76 under "IIV, CV%" column
    # but MG-ADL model uses additive IIV per paper Methods; values are interpreted as NONMEM
    # OMEGA variances in squared parameter units). Correlation between etaIDecplacebo and
    # etaSIgG is -0.733; covariance = -0.733 * sqrt(1.89) * sqrt(3.76) = -1.954.
    etaIDecplacebo + etaSIgG ~ c(1.89,
                                 -1.954, 3.76)
    etaSplacebo ~ 0.125     # Table 3 IIV (additive; variance in (pts/week)^2)

    # ---- Residual error (Table 3) ----
    # Paper reports PK additive RUVs in nmol/L; converted to ug/mL here using
    # nipocalimab MW 142 kDa so the residual applies directly to Cc (ug/mL).
    CcaddSdELISA   <- 0.0632; label("Additive PK residual SD for ELISA-assay observations (ug/mL)")        # Table 3: 0.445 nmol/L x 142/1000 = 0.0632 ug/mL
    CcaddSdECLIA   <- 0.00486; label("Additive PK residual SD for ECLIA-assay observations (ug/mL)")       # Table 3: 0.0342 nmol/L x 142/1000 = 0.00486 ug/mL
    CcpropSdPhase1 <- 0.0834; label("Proportional PK residual for Phase 1 observations (fraction)")        # Table 3
    CcpropSdPhase2 <- 0.367;  label("Proportional PK residual for Phase 2 observations (fraction)")        # Table 3

    pctUnoccupiedFcRnaddSd  <- 2.98;  label("Additive residual SD for % unoccupied FcRn (percent)")        # Table 3
    pctUnoccupiedFcRnpropSd <- 0.227; label("Proportional residual for % unoccupied FcRn (fraction)")       # Table 3

    IgG_obspropSd <- 0.0858; label("Proportional residual for total serum IgG (fraction)")                  # Table 3
    dMGADLaddSd   <- 1.50;   label("Additive residual SD for MG-ADL change from baseline (points)")          # Table 3
  })

  model({
    # ---- Physical constants ----
    # Nipocalimab MW derived from paper: 15 mg/kg x 75 kg = 1125 mg; paper states this delivers 7.92 umol.
    # MW = 1125 mg / 7.92e-3 mmol = 142045 g/mol; rounded to 142000 g/mol (~142 kDa).
    MW <- 142000             # Nipocalimab molecular weight, g/mol (derived from paper results section)

    # ---- Individual PK parameters with allometric weight scaling (reference 75 kg) ----
    cl <- exp(lcl + etalcl) * (WT / 75)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / 75)^allo_v
    q  <- exp(lq)           * (WT / 75)^allo_cl
    vp <- exp(lvp)          * (WT / 75)^allo_v

    # ---- Individual FcRn target parameters ----
    # FcRn0 stored in nmol/L (paper reporting unit); converted to ug/mL here so
    # all drug-FcRn QSS arithmetic runs on one consistent mass-concentration
    # scale (drug Cc is in ug/mL with dosing in mg and volumes in L).
    # 143 nmol/L * MW[g/mol] / 1e6 = 143 * 142000 / 1e6 = 20.3 ug/mL.
    FcRn0 <- exp(lFcRn0 + etalFcRn0) * MW / 1e6        # ug/mL (total FcRn at baseline)
    kint  <- exp(lkint)                                # 1/day
    kdeg  <- exp(lkdeg)                                # 1/day
    kss   <- exp(lkss)                                 # ug/mL (paper reports Kss in ug/mL directly)

    # Ksyn from steady-state constraint on the accessible FcRn pool:
    #   dtotal_target/dt|t=0 = 0 -> Ksyn = kdeg * total_target(0) = kdeg * FcRn0 * FRmax
    ksyn <- kdeg * FcRn0 * FRmax

    # ---- Individual IgG parameters ----
    # IgG baseline: IgG0 for non-Vivacity-MG subjects; IgG0 * FRIgG0_M281_004 for MOM-M281-004 subjects.
    IgG_BL <- exp(lIgG0 + etalIgG0) * (1 - (1 - FRIgG0_M281_004) * STUDY_M281_004)
    kdeg_IgG <- exp(lkdeg_IgG + etalkdeg_IgG)                 # 1/day
    IgK_i <- exp(lIgK + etalIgK)                              # unitless
    # FcRn-mediated IgG recycling rate (Eq. 5): krec = kdeg * IgK / (1 + IgK).
    krec_IgG <- kdeg_IgG * IgK_i / (1 + IgK_i)
    # Steady-state IgG synthesis (Eq. 4): ksyn_IgG = (kdeg_IgG - krec_IgG) * IgG_baseline.
    ksyn_IgG <- (kdeg_IgG - krec_IgG) * IgG_BL

    # ---- Individual MG-ADL parameters ----
    ke0         <- exp(lke0)
    IDec_i      <- (IDecplacebo + etaIDecplacebo) * (MGADL / 7)^EADL
    SIgG_i      <- (SIgG + etaSIgG)
    Splacebo_i  <- (Splacebo + etaSplacebo)

    # Dose in mg (rxode2 convention) goes directly into the central compartment
    # as an amount in mg; with vc in L, Cc = central/vc is in mg/L = ug/mL.

    # ---- Initial conditions for endogenous states ----
    # Accessible FcRn pool at baseline (Eq. 2: Rtotal(time 0) = Rtotal,0 * FRmax):
    total_target(0) <- FcRn0 * FRmax
    # Total IgG at baseline:
    total_IgG(0) <- IgG_BL
    # effect (IgG effect compartment for MG-ADL) starts at 0 (no IgG reduction yet).
    effect(0) <- 0

    # ---- QSS TMDD relations in the central compartment ----
    # Ctotal = central / vc (ug/mL; total drug concentration, free + FcRn-bound).
    ctot  <- central / vc
    disc  <- ctot - total_target - kss
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    complex <- total_target * cfree / (kss + cfree)

    # ---- Drug, FcRn, and IgG ODEs (Supplement Equations 1, 4) ----
    # Central total drug (mg): elimination via linear CL on Cfree, TMDD internalization on complex,
    # and distribution to peripheral (Q on Cfree; return as Q/Vp on peripheral1).
    d/dt(central)      <- -cl * cfree - kint * vc * complex - q * cfree + (q / vp) * peripheral1
    d/dt(peripheral1)  <-  q * cfree - (q / vp) * peripheral1
    # Accessible FcRn pool (ug/mL). Eq. 1 rewritten with complex = total_target*cfree/(kss+cfree).
    d/dt(total_target) <-  ksyn - kdeg * total_target - (kint - kdeg) * complex
    # Total serum IgG (g/L). Eq. 4; F_free is the fraction of accessible FcRn still free vs baseline.
    f_free <- (total_target - complex) / (FcRn0 * FRmax)
    d/dt(total_IgG)    <-  ksyn_IgG - (kdeg_IgG - krec_IgG * f_free) * total_IgG
    # IgG effect compartment (dimensionless fraction, 0 -> 1 at full blockade). Supplement Eq. 6.
    d/dt(effect)       <-  ke0 * ((1 - total_IgG / IgG_BL) - effect)

    # ---- Observation variables ----
    # Total nipocalimab concentration (ug/mL) -- ELISA/ECLIA assays measure total drug.
    Cc <- ctot

    # Percent unoccupied FcRn (Eq. 2). Includes the inaccessible (1-FRmax) fraction in the
    # numerator so the observation matches the paper's "% unoccupied FcRn" scale (100% at baseline,
    # near 0% at full RO). Unit-agnostic ratio, so ug/mL vs nmol/L makes no difference here.
    Rfree <- total_target - complex + FcRn0 * (1 - FRmax)
    pctUnoccupiedFcRn <- Rfree / FcRn0 * 100

    # Total serum IgG (g/L).
    IgG_obs <- total_IgG

    # MG-ADL absolute change from baseline (points). Placebo effect turns on after first dose
    # (t > 0 gate). IgG-driven effect: S_IgG is reported per 10% IgG reduction; multiply by 10
    # to convert the fractional effect-compartment output to "number of 10% reductions".
    placebo_cfb <- (t > 0) * (IDec_i + Splacebo_i * t / 7)
    igg_cfb     <- 10 * SIgG_i * effect * (MGADL / 7)^EIgG
    dMGADL      <- placebo_cfb + igg_cfb

    # ---- Residual-error models (per-assay / per-phase switches on PK; single form on others) ----
    CcaddSd   <- CcaddSdELISA   * ELISA  + CcaddSdECLIA   * (1 - ELISA)
    CcpropSd  <- CcpropSdPhase1 * PHASE1 + CcpropSdPhase2 * (1 - PHASE1)

    Cc                ~ add(CcaddSd) + prop(CcpropSd)
    pctUnoccupiedFcRn ~ add(pctUnoccupiedFcRnaddSd) + prop(pctUnoccupiedFcRnpropSd)
    IgG_obs           ~ prop(IgG_obspropSd)
    dMGADL            ~ add(dMGADLaddSd)
  })
}
