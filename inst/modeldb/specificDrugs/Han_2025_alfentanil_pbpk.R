Han_2025_alfentanil_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, dynamic age-dependent; custom Phoenix WinNonlin 8.4",
    "implementation). Alfentanil disposition across the human age range,",
    "from preterm neonates through adults to the elderly. Fifteen",
    "perfusion-limited tissues plus a five-segment gut lumen / gut wall",
    "cascade; CYP3A4 elimination in liver and in the duodenal, jejunal and",
    "ileal wall. Structurally identical to the midazolam, fentanyl,",
    "alfentanil and sufentanil siblings from this paper; only the",
    "drug-specific tissue:plasma partition coefficients, blood:plasma ratio,",
    "unbound fraction and intrinsic clearance differ. Organ volumes, organ",
    "blood flows, hepatic microsomal protein, CYP3A4 activity and plasma",
    "protein binding are all driven by AGE, so one model covers paediatric,",
    "adult and geriatric populations. INTRAVENOUS ONLY: Han 2025 report no",
    "effective permeability for alfentanil (Table S1 gives '/'), and every one",
    "of the 15 clinical data sets in Table S10 is an intravenous bolus or",
    "infusion, so peff is fixed at 0 and the luminal transit chain is inert.",
    "Han 2025 build their 1000 virtual individuals by varying Peff, fu,b,",
    "CLint,l, CLint,i, Rb and PBSF over 0.67-1.5x the ideal value; that sweep",
    "is carried here as exponential IIV of variance 0.05415 on the seven",
    "corresponding control-stream parameters that are non-zero for this drug",
    "(see the vignette Errata for the derivation)."
  )
  reference <- paste(
    "Han J, Zhang Z, Liu X, Yang H, Liu L (2025). Prediction of",
    "Pharmacokinetics for CYP3A4-Metabolized Drugs in Pediatrics and",
    "Geriatrics Using Dynamic Age-Dependent Physiologically Based",
    "Pharmacokinetic Models. Pharmaceutics 17(2):214.",
    "doi:10.3390/pharmaceutics17020214.",
    "ODE system from Supplementary Information Eq S1-S10; drug parameters",
    "from Table S1; adult system physiology from Table S2; paediatric",
    "scaling from Tables S3-S5 and Eq 3-17; geriatric scaling from Tables",
    "S6-S7 and Eq 18-22; clinical data sets in Table S10; validation",
    "targets from Table S14. The deposited",
    "Supplementary File S2 (Sobol global-sensitivity R script, rxode2) is",
    "the authors' own executable encoding of the adult midazolam ODEs and",
    "is used here as a second source for the equation forms.",
    sep = " "
  )
  vignette <- "Han_2025_cyp3a4_pediatric_geriatric_pbpk"
  units    <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # Perfused gut-WALL tissue states. The gut LUMEN segments use the canonical
  # stomach / duodenum / jejunum / ileum / cecum / colon names; the tissue
  # (wall) states between lumen and portal drainage have no canonical, so they
  # are declared paper-specific following the `wall_<segment>` precedent set by
  # Luo_2024_remimazolam_pbpk.R (same research group, same model lineage).
  paper_specific_compartments <- c(
    "wall_stomach", "wall_duodenum", "wall_jejunum",
    "wall_ileum", "wall_cecum", "wall_colon"
  )

  compartmentData <- list(
    lung          = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    kidney        = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    heart         = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    liver         = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    muscle        = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    skin          = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    adipose       = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    brain         = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    venous        = list(analyte = "alfentanil", units = "mg", specimen = "whole blood",  verified = TRUE),
    spleen        = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    arterial      = list(analyte = "alfentanil", units = "mg", specimen = "whole blood",  verified = TRUE),
    other         = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    stomach       = list(analyte = "alfentanil", units = "mg", specimen = "administration site",  verified = TRUE),
    duodenum      = list(analyte = "alfentanil", units = "mg", specimen = "administration site",  verified = TRUE),
    jejunum       = list(analyte = "alfentanil", units = "mg", specimen = "administration site",  verified = TRUE),
    ileum         = list(analyte = "alfentanil", units = "mg", specimen = "administration site",  verified = TRUE),
    cecum         = list(analyte = "alfentanil", units = "mg", specimen = "administration site",  verified = TRUE),
    colon         = list(analyte = "alfentanil", units = "mg", specimen = "administration site",  verified = TRUE),
    wall_stomach  = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    wall_duodenum = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum  = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum    = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    wall_cecum    = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE),
    wall_colon    = list(analyte = "alfentanil", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Chronological age in years; drives every age-dependent system parameter.",
      units              = "years",
      type               = "continuous",
      notes              = paste(
        "The model switches between three physiology regimes on AGE:",
        "paediatric (AGE < 18, Tables S3-S5 and Eq 3-17), adult",
        "(18 <= AGE <= 60, the fixed Table S2 reference values) and geriatric",
        "(AGE > 60, Tables S6-S7 and Eq 18-22). Han 2025 Section 2.7 bins",
        "the paediatric range as 0-0.083, 0.083-0.5, 0.5-1, 1-2, 2-5, 5-9,",
        "9-12, 12-15 and 15-18 years and the geriatric range as 60-65,",
        "65-75, 75-85 and >= 85 years, but the model itself is continuous in",
        "AGE. Because the three regimes come from three different literature",
        "sources, the physiology is NOT continuous across the 18-year and",
        "60-year boundaries; see the vignette's Errata section."
      ),
      source_name        = "Age (year)"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Read only by the geriatric branch (AGE > 60). Table S6 / S7",
        "footnote states 'Male sex = 0. Female sex = 1', which matches the",
        "canonical SEXF orientation exactly, so no sign inversion is needed.",
        "The paediatric equations of Tables S3-S5 and the adult Table S2",
        "reference values carry no sex term, so SEXF has no effect below",
        "60 years."
      ),
      source_name        = "Sex"
    ),
    PAGE = list(
      description        = "Postmenstrual age in months (gestational age at birth + postnatal age).",
      units              = "months",
      type               = "continuous",
      notes              = paste(
        "Han 2025 writes postmenstrual age (PMA) in YEARS in Eq 8 and",
        "Eq 13; the canonical PAGE carries months, so model() converts with",
        "PMA_years = PAGE / 12. Load-bearing only for neonates and preterm",
        "infants (PNA <= 1 month), where it drives the hepatic CYP3A4",
        "activity ratio (Eq 8) and the preterm plasma albumin curve",
        "(Eq 13). For any older subject set PAGE = 12 * AGE + GA / 4.348,",
        "where GA is gestational age at birth in weeks; the value is then",
        "not read."
      ),
      source_name        = "PMA (postmenstrual age)"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth) in months.",
      units              = "months",
      type               = "continuous",
      notes              = paste(
        "Two roles. (1) PNA <= 1 month selects the neonatal / preterm",
        "branch, which uses Eq 8 for hepatic CYP3A4 activity and Eq 13 for",
        "plasma albumin instead of Eq 7 and Eq 14. (2) Eq 16 is written in",
        "DAYS after birth, so model() converts with days = 30.4375 * PNA.",
        "For a non-neonate set PNA = 12 * AGE."
      ),
      source_name        = "PNA (postnatal age)"
    ),
    ROUTE_IV = list(
      description        = "Route indicator for the dose being simulated, 1 = intravenous, 0 = oral.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral administration)",
      notes              = paste(
        "Gates the gut-wall unbound fraction only. Supplementary Eq S5 and",
        "its surrounding text state that fu,gut 'was defaulted to 1 for oral",
        "administration and equaled to free fraction of drug (fu,b) in blood",
        "for intravenous administration', so ROUTE_IV = 0 gives fu_gut = 1",
        "(full first-pass gut extraction of the absorbed flux) and",
        "ROUTE_IV = 1 gives fu_gut = fu_b. This is the same gating role the",
        "ROUTE_IV canonical plays in vandenBerg_2021_uprifosbuvir_pbpk.R.",
        "Note the reference category here is ORAL, not the SC reference of",
        "the canonical's founding paper. Han 2025 dose this drug",
        "intravenously only, so ROUTE_IV = 1 on every record the paper",
        "simulates; the oral setting is retained for structural parity with",
        "the midazolam sibling but cannot absorb, because peff is 0."
      ),
      source_name        = "Route of administration"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1000L,
    n_studies      = 15L,
    age_range      = "gestational age 29.5 weeks (preterm neonates) to 91 years",
    weight_range   = "0.96 kg (preterm neonates) to 85 kg (adults)",
    disease_state  = paste(
      "Preterm and term neonates undergoing stressful intensive-care procedures, children in elective orthopaedic and general surgery, healthy adult volunteers, cardiac-surgery patients, and elderly patients in intra-abdominal surgery.",
      "Han 2025 Table S10 lists 15 alfentanil data sets digitised from the",
      "published literature (5 paediatric, 9 adult, 1 geriatric);",
      "Table S14 gives the observed and predicted AUC0-t for each."
    ),
    dose_range     = paste(
      "11.8-125 ug/kg intravenous bolus over 0.5-2 min, plus continuous infusions of 3-6 ug/kg/min"
    ),
    regions        = "Multinational literature compilation",
    notes          = paste(
      "The 1000 'subjects' are virtual individuals, not real patients: Han",
      "2025 Section 2.1 constructs them by varying Peff, fu,b, CLint,l,",
      "CLint,i, Rb and PBSF uniformly between 0.67 and 1.5 times the ideal",
      "value for the age stage. That sweep is carried here as exponential IIV",
      "of variance 0.05415 on the corresponding control-stream parameters;",
      "rxode2::zeroRe() recovers the typical-value predictor.",
      "Table S14 reports observed and predicted AUC0-t only; no Cmax column is given for alfentanil.",
      "Model performance over all four drugs and all three age bands is",
      "reported in Table S16: AFE 1.1 / 0.9 / 1.1 and PE 24.0 / 15.1 /",
      "24.2 % for paediatrics / adults / geriatrics."
    )
  )

  ini({
    # ================= ALFENTANIL drug-specific parameters =================
    # All from Han 2025 Table S1, alfentanil column. The Kt:p values of fentanyl,
    # alfentanil and sufentanil were calculated by the method of Table S1
    # reference [21] from pKa 6.5 and logP 2.1.
    kp_lung    <- fixed(0.32);        label("alfentanil lung:plasma partition coefficient Kt:p (unitless)")      # Table S1, Lung
    kp_kidney  <- fixed(0.26);        label("alfentanil kidney:plasma partition coefficient Kt:p (unitless)")    # Table S1, Kidney
    kp_heart   <- fixed(0.24);        label("alfentanil heart:plasma partition coefficient Kt:p (unitless)")     # Table S1, Heart
    kp_liver   <- fixed(0.28);        label("alfentanil liver:plasma partition coefficient Kt:p (unitless)")     # Table S1, Liver
    kp_muscle  <- fixed(0.19);        label("alfentanil muscle:plasma partition coefficient Kt:p (unitless)")    # Table S1, Muscle
    kp_skin    <- fixed(0.65);        label("alfentanil skin:plasma partition coefficient Kt:p (unitless)")      # Table S1, Skin
    kp_adipose <- fixed(8.22);        label("alfentanil adipose:plasma partition coefficient Kt:p (unitless)")   # Table S1, Adipose
    kp_brain   <- fixed(0.46);        label("alfentanil brain:plasma partition coefficient Kt:p (unitless)")     # Table S1, Brain
    kp_spleen  <- fixed(0.18);        label("alfentanil spleen:plasma partition coefficient Kt:p (unitless)")    # Table S1, Spleen
    kp_stomach <- fixed(0.47);        label("alfentanil stomach:plasma partition coefficient Kt:p (unitless)")   # Table S1, Stomach
    kp_gut     <- fixed(0.47);        label("alfentanil intestine:plasma partition coefficient Kt:p (unitless)") # Table S1, Intestine (all five gut-wall segments)
    kp_other   <- fixed(0.001);   label("alfentanil rest-of-body:plasma partition coefficient Kt:p (unitless)")  # Table S1, ROB (footnote a: assumed, not calculated)
    lbpr       <- fixed(log(0.63));   label("alfentanil blood:plasma concentration ratio Rb (log, unitless)")         # Table S1, Rb [ref 16]
    lfu_p      <- fixed(log(0.086));  label("alfentanil fraction unbound in adult plasma (log, unitless)")            # Table S1, fu,p [ref 13]
    lcl_int_h  <- fixed(log(0.096));  label("alfentanil CYP3A4 hepatic intrinsic clearance (log mL/min/mg microsomal protein)")  # Table S1, CLint [ref 19]
    # Table S1 reports "/" for the alfentanil effective permeability: Han 2025
    # never simulate an oral alfentanil dose, and every data set in Table S10 is
    # intravenous. Encoded as an explicit fixed(0) rather than invented, so
    # the luminal states carry transit but never absorb. Do NOT dose into
    # `stomach` with this model.
    peff       <- fixed(0);       label("alfentanil effective intestinal permeability Peff (cm/min)")           # Table S1: not reported (intravenous administration only)

    # ================= System parameters carrying between-subject variability =================
    # PBSF and the three intestinal-wall intrinsic clearances are system /
    # drug-system hybrids rather than Table S1 drug properties, but Han 2025
    # give them random effects, so they are lifted out of model() into ini().
    lpbsf      <- fixed(log(64220));  label("Adult hepatic microsomal protein scaling factor PBSF (log mg)")  # Table S2 footnote: 1690 g liver x 38 mg microsomal protein/g; Phoenix control stream tvPBSF = 64220
    # Supplementary Eq S5 gives duodenal / jejunal / ileal CYP3A contents of
    # 9.7 / 38.4 / 22.4 nmol and states that CLint,gwi was scaled on CYP3A
    # content in intestine and liver, but never prints the resulting mL/min.
    # The Phoenix control stream does, for midazolam only: tvClintg1/2/3 =
    # 34 / 134.6 / 78.51 at CLint,l = 0.389. Those are exactly proportional to
    # the three nmol contents, giving 9.0107 mg microsomal protein per nmol
    # CYP3A (111 pmol CYP3A4 per mg hepatic microsomal protein), so the
    # mg-equivalents are 87.4036 / 346.0154 / 201.8252. The alfentanil values
    # below are that constant times the Table S1 alfentanil CLint,l of 0.096:
    # DERIVED, not tabulated. See the vignette Errata.
    lcl_gw1    <- fixed(log(8.390746));  label("alfentanil duodenal wall CYP3A4 intrinsic clearance (log mL/min)")  # 87.4036 x 0.096
    lcl_gw2    <- fixed(log(33.217478));  label("alfentanil jejunal wall CYP3A4 intrinsic clearance (log mL/min)")   # 346.0154 x 0.096
    lcl_gw3    <- fixed(log(19.375219));  label("alfentanil ileal wall CYP3A4 intrinsic clearance (log mL/min)")     # 201.8252 x 0.096

    # ================= Between-subject variability =================
    # Han 2025 Section 2.1 sets Peff, fu,b, CLint,l, CLint,i, Rb and PBSF as
    # the random-effect parameters of the virtual individuals and states they
    # vary between 0.67 and 1.5 times the ideal value for the age stage. The
    # deposited Phoenix control stream implements them as exp(eta) on exactly
    # those parameters, splitting CLint,i into its three absorbing segments:
    # Clint, Clintg1, Clintg2, Clintg3, fub, PBSF, Rb, Peff1.
    #
    # The control stream prints ranef(diag(n)) = 0.405465108 for all eight.
    # That number is ln(1.5) to nine decimals -- the UPPER BOUND of the paper's
    # own 0.67-1.5x sweep pasted into a variance slot, not a variance. Read
    # literally it gives about 2.7x the between-subject spread Han 2025
    # themselves publish. The value encoded here is instead the log-scale
    # variance of Uniform[ln 0.67, ln 1.5], (log(1.5) - log(0.67))^2 / 12 =
    # 0.05413, which reproduces both published CV% values in Table S12 for
    # midazolam (42.6 simulated vs 41.8 published, oral 7.5 mg; 27.3 vs 26.4,
    # intravenous 5 mg). Derived, not printed: see the vignette Errata.
    #
    # Peff carries no random effect in this file: Table S1 reports no alfentanil
    # permeability, peff is fixed at 0, and log(0) has no lognormal.
    etalbpr       ~ fixed(0.05415);  label("IIV on blood:plasma ratio Rb (variance, log scale)")
    etalfu_p      ~ fixed(0.05415);  label("IIV on fraction unbound in plasma (variance, log scale)")
    etalcl_int_h  ~ fixed(0.05415);  label("IIV on hepatic intrinsic clearance (variance, log scale)")
    etalpbsf      ~ fixed(0.05415);  label("IIV on hepatic microsomal protein scaling factor (variance, log scale)")
    etalcl_gw1    ~ fixed(0.05415);  label("IIV on duodenal wall intrinsic clearance (variance, log scale)")
    etalcl_gw2    ~ fixed(0.05415);  label("IIV on jejunal wall intrinsic clearance (variance, log scale)")
    etalcl_gw3    ~ fixed(0.05415);  label("IIV on ileal wall intrinsic clearance (variance, log scale)")

    # ================= Residual error =================
    # Han 2025 report no residual-error model: the PBPK model is a forward
    # predictor assessed by 5th-95th percentile coverage and 0.5-2-fold
    # AUC/Cmax agreement (Section 2.3), not by a fitted sigma. Non-fitted
    # placeholder per the in-repo PBPK convention (Luo_2024_remimazolam_pbpk.R,
    # An_2012_mitoxantrone_human_pbpk.R).
    propSd     <- fixed(0.10);    label("Proportional residual error placeholder, alfentanil (fraction)")    # not reported by Han 2025
  })

  model({
    # =====================================================================
    # 0. Drug-independent switches
    #    bind_alb selects which plasma protein drives the age-dependent
    #    unbound fraction of Eq 17. Eq 17 admits either protein ("P is the
    #    level of plasma ALB (or AAG)"). Alfentanil is a lipophilic base
    #    (Table S1 pKa 6.5, logP 2.1) and binds predominantly to
    #    alpha1-acid glycoprotein, so the AAG ontogeny of Eq 15-16 drives
    #    its unbound fraction; the midazolam sibling uses albumin and sets
    #    bind_alb <- 1. Han 2025 do not state the assignment per drug; see
    #    the vignette Errata.
    # =====================================================================
    bind_alb <- 0

    # Individual values of the seven random-effect parameters. At eta = 0
    # every one collapses to its Table S1 / Table S2 / control-stream typical
    # value, so the model as shipped predicts the paper's typical profile and
    # rxode2::zeroRe() recovers it exactly.
    bpr      <- exp(lbpr + etalbpr)
    fu_p     <- exp(lfu_p + etalfu_p)
    cl_int_l <- exp(lcl_int_h + etalcl_int_h)
    pbsf     <- exp(lpbsf + etalpbsf)

    # =====================================================================
    # 1. Age regime. Han 2025 develops the adult model on the fixed Table S2
    #    physiology (Section 2.3), extrapolates DOWN to paediatrics with
    #    Tables S3-S5 (Section 2.4) and UP to geriatrics with Tables S6-S7
    #    (Section 2.5). Section 2.3 defines adults as 20-59 years and
    #    Section 2.7 starts the geriatric bins at 60, so the adult plateau
    #    is taken as 18 < AGE <= 60.
    #    a_ped / a_ger are CLAMPED copies of AGE: the Table S4 muscle
    #    blood-flow rational function has a pole at A = 21.9 years, so the
    #    unclamped polynomial would return Inf on an adult record and
    #    0 * Inf = NaN would poison the sum below.
    # =====================================================================
    is_ped <- (AGE < 18)
    is_ger <- (AGE > 60)
    is_adl <- 1 - is_ped - is_ger
    a_ped  <- AGE * is_ped + 18 * (1 - is_ped)
    a_ger  <- AGE * is_ger + 65 * (1 - is_ger)

    # =====================================================================
    # 2. Body size.
    #    Paediatric BW (kg) and HT (cm): Eq 3 with the Table S5 coefficients,
    #    which switch at 2 years. Geriatric HT / BW: Table S6 rows 1-2.
    #    Adult anchor: BW 70 kg is PROVEN by Table S2 -- the thirteen adult
    #    tissue volumes sum to exactly 70,000 mL (1170 + 310 + 1450 + 35000
    #    + 10000 + 7800 + 280 + 190 + 160 + 1690 + 1650 + 5200 + 5100).
    #    Adult HT 176 cm is the intercept of the Table S6 male height
    #    equation at A = 0.
    #    BSA: Eq 5 (BW > 15 kg) and Eq 6 (BW < 15 kg).
    # =====================================================================
    u_lt2  <- (a_ped <= 2)
    cw_a   <- 3.41 * u_lt2 + 9.86 * (1 - u_lt2)          # Table S5 weight a
    cw_b   <- 20.1 * u_lt2 + 0.370 * (1 - u_lt2)         # Table S5 weight b
    cw_c   <- 1.46 * u_lt2 + (-0.0789) * (1 - u_lt2)     # Table S5 weight c
    cw_d   <- (-0.107) * u_lt2 + 0.00205 * (1 - u_lt2)   # Table S5 weight d
    ch_a   <- 50.5 * u_lt2 + 74.8 * (1 - u_lt2)          # Table S5 height a
    ch_b   <- 117 * u_lt2 + 3.35 * (1 - u_lt2)           # Table S5 height b
    ch_c   <- 1.28 * u_lt2 + (-0.0354) * (1 - u_lt2)     # Table S5 height c
    ch_d   <- (-0.0696) * u_lt2 + 0.00125 * (1 - u_lt2)  # Table S5 height d
    bw_ped <- (cw_a + cw_b * a_ped) / (1 + cw_c * a_ped + cw_d * a_ped * a_ped)  # Eq 3
    ht_ped <- (ch_a + ch_b * a_ped) / (1 + ch_c * a_ped + ch_d * a_ped * a_ped)  # Eq 3
    ht_ger <- -0.0039 * a_ger * a_ger + 0.238 * a_ger - 12.5 * SEXF + 176        # Table S6 Height
    bw_ger <- -0.0039 * a_ger * a_ger + 1.12 * ht_ger + 0.611 * a_ger - 0.424 * SEXF - 137  # Table S6 Body weight
    bw     <- is_ped * bw_ped + is_adl * 70 + is_ger * bw_ger
    ht     <- is_ped * ht_ped + is_adl * 176 + is_ger * ht_ger
    bsa    <- (bw > 15) * (0.007184 * bw^0.425 * ht^0.725) +
              (bw <= 15) * (0.024265 * bw^0.5378 * ht^0.3964)   # Eq 5 / Eq 6

    # =====================================================================
    # 3. Organ masses (g; Table S3 footnote "assuming an average specific
    #    gravity of 1 g/mL", so g == mL throughout).
    #    Paediatric: Table S3 polynomials in A (years).
    #    Geriatric:  Table S6 (kg and, for stomach, g -> converted to g).
    #    Adult:      Table S2 volume column.
    # =====================================================================
    wp_lung    <- -0.0346 * a_ped^4 + 1.5069 * a_ped^3 - 20.31 * a_ped^2 + 123.99 * a_ped + 59.231
    wp_heart   <- -0.0132 * a_ped^4 + 0.5051 * a_ped^3 - 5.7113 * a_ped^2 + 32.213 * a_ped + 20.364
    wp_brain   <- 10000 * (a_ped + 0.213) / (6.030 + 6.895 * a_ped)
    wp_muscle  <- 0.535 * a_ped^3 + 56.937 * a_ped^2 - 124.25 * a_ped + 1051.3
    wp_adipose <- 0.0165 * a_ped^5 - 1.9784 * a_ped^4 + 51.963 * a_ped^3 - 459.38 * a_ped^2 + 1566.8 * a_ped + 1004.2
    wp_skin    <- -0.0992 * a_ped^4 + 4.2762 * a_ped^3 - 62.165 * a_ped^2 + 432.78 * a_ped + 203.2
    wp_kidney  <- 0.0009737 * a_ped^5 - 0.0561 * a_ped^4 + 1.1729 * a_ped^3 - 10.34 * a_ped^2 + 44.604 * a_ped + 28.291
    wp_spleen  <- -0.0091 * a_ped^4 + 0.3457 * a_ped^3 - 4.0754 * a_ped^2 + 22.269 * a_ped + 11.05
    wp_stomach <- 0.0008 * a_ped^5 - 0.0356 * a_ped^4 + 0.5823 * a_ped^3 - 4.0437 * a_ped^2 + 17.888 * a_ped + 7.54
    wp_liver   <- 0.0072 * a_ped^5 - 0.3975 * a_ped^4 + 7.9052 * a_ped^3 - 65.624 * a_ped^2 + 262.02 * a_ped + 157.52
    wp_gut     <- -0.047817 * a_ped^4 + 1.925 * a_ped^3 - 22.382 * a_ped^2 + 107.09 * a_ped + 51.125
    wp_blood   <- -0.0623 * a_ped^5 + 2.4425 * a_ped^4 - 31.37 * a_ped^3 + 149.98 * a_ped^2 + 31.305 * a_ped + 393.7

    wg_lung    <- 1000 * exp(0.028 * ht_ger + 0.0077 * a_ger - 5.6)          # Table S6 Lung (kg)
    wg_adipose <- 1000 * (0.68 * bw_ger - 0.56 * ht_ger + 6.1 * SEXF + 65)   # Table S6 Adipose (kg)
    wg_brain   <- 1000 * exp(-0.0075 * a_ger + 0.0078 * ht_ger - 0.97)       # Table S6 Brain (kg)
    wg_heart   <- 1000 * (0.34 * bsa + 0.0018 * a_ger - 0.36)                # Table S6 Heart (kg)
    wg_kidney  <- 1000 * (-0.00038 * a_ger - 0.056 * SEXF + 0.33)            # Table S6 Kidney (kg)
    wg_muscle  <- 1000 * (17.9 * bsa - 0.0667 * a_ger - 5.68 * SEXF - 1.22)  # Table S6 Muscle (kg)
    wg_skin    <- 1000 * exp(-0.0058 * a_ger - 0.37 * SEXF + 1.13)           # Table S6 Skin (kg)
    wg_stomach <- 160                                                        # Table S6 Stomach (g), footnote a: assumed equal to adults
    wg_gut     <- 1000 * (3e-6 * ht_ger^2.49)                                # Table S6 Intestine (kg)
    wg_spleen  <- 1000 * exp(1.13 * bsa - 3.93)                              # Table S6 Spleen (kg)
    wg_liver   <- 1000 * exp(0.87 * bsa - 0.0014 * a_ger - 1)                # Table S6 Liver (kg)
    wg_blood   <- 1000 * exp(0.067 * bsa - 0.0025 * a_ger - 0.38 * SEXF + 1.7)  # Table S6 Blood (kg)

    v_lung    <- is_ped * wp_lung    + is_adl * 1170  + is_ger * wg_lung     # Table S2 Lung 1170 mL
    v_heart   <- is_ped * wp_heart   + is_adl * 310   + is_ger * wg_heart    # Table S2 Heart 310 mL
    v_brain   <- is_ped * wp_brain   + is_adl * 1450  + is_ger * wg_brain    # Table S2 Brain 1450 mL
    v_muscle  <- is_ped * wp_muscle  + is_adl * 35000 + is_ger * wg_muscle   # Table S2 Muscle 35000 mL
    v_adipose <- is_ped * wp_adipose + is_adl * 10000 + is_ger * wg_adipose  # Table S2 Adipose 10000 mL
    v_skin    <- is_ped * wp_skin    + is_adl * 7800  + is_ger * wg_skin     # Table S2 Skin 7800 mL
    v_kidney  <- is_ped * wp_kidney  + is_adl * 280   + is_ger * wg_kidney   # Table S2 Kidneys 280 mL
    v_spleen  <- is_ped * wp_spleen  + is_adl * 190   + is_ger * wg_spleen   # Table S2 Spleen 190 mL
    v_stomach <- is_ped * wp_stomach + is_adl * 160   + is_ger * wg_stomach  # Table S2 Stomach 160 mL
    v_liver   <- is_ped * wp_liver   + is_adl * 1690  + is_ger * wg_liver    # Table S2 Liver 1690 mL
    v_gut     <- is_ped * wp_gut     + is_adl * 1650  + is_ger * wg_gut      # Table S2 duodenum+jejunum+ileum+cecum+colon = 70+209+139+116+1116
    v_blood   <- is_ped * wp_blood   + is_adl * 5200  + is_ger * wg_blood    # Table S2 Vein 3470 + Artery 1730

    # Venous / arterial split and gut-segment split held at the Table S2 adult
    # proportions: Tables S3 and S6 give only a single lumped blood volume and
    # a single lumped intestine mass, so the five segments and the two blood
    # pools are apportioned by their adult fractions. Mass is conserved by
    # construction.
    v_venous   <- 0.667307692 * v_blood   # 3470 / 5200
    v_arterial <- 0.332692308 * v_blood   # 1730 / 5200
    v_gw1 <- 0.042424242 * v_gut          # duodenum 70 / 1650
    v_gw2 <- 0.126666667 * v_gut          # jejunum 209 / 1650
    v_gw3 <- 0.084242424 * v_gut          # ileum 139 / 1650
    v_gw4 <- 0.070303030 * v_gut          # cecum 116 / 1650
    v_gw5 <- 0.676363636 * v_gut          # colon 1116 / 1650

    # Rest of body closes the whole-body volume balance. In the adult this is
    # exact: 70,000 mL body volume minus the 64,900 mL of tabulated tissue and
    # blood leaves precisely the 5100 mL Table S2 ROB entry.
    v_other0 <- 1000 * bw - (v_lung + v_heart + v_brain + v_muscle + v_adipose +
                             v_skin + v_kidney + v_spleen + v_stomach + v_liver +
                             v_gut + v_blood)
    v_other  <- (v_other0 > 100) * v_other0 + (v_other0 <= 100) * 100

    # =====================================================================
    # 4. Blood flows (mL/min).
    #    Adult: Table S2 blood-flow column.
    #    Paediatric: Table S4. Its rows carry THREE different units, which
    #      must be reconciled before use (see the vignette Errata):
    #        - cardiac output and brain are stated in L/h  -> x 1000/60;
    #        - heart and adipose are explicitly per gram of organ
    #          (footnote a: 0.79 * W_heart, 0.024 * W_adipose);
    #        - the seven remaining ref-[29] rational functions are tissue
    #          perfusion per 100 g of organ -> x W_organ / 100. Read this
    #          way skin, kidney, muscle and intestine land within 20 % of
    #          their adult Table S2 values at 18 years, and the paediatric
    #          predictions reproduce the paper's own Table S12 Cmax to
    #          within 1 %; read as flat mL/min they do not (muscle would be
    #          2.5 mL/min at age 18).
    #    Geriatric: Table S7, where every tissue except stomach is a percent
    #      of cardiac output.
    # =====================================================================
    co_ped  <- (0.012 * a_ped^3 - 1.2144 * a_ped^2 + 40.324 * a_ped + 44.414) * 1000 / 60
    qp_heart   <- 0.79 * wp_heart
    qp_brain   <- (-0.0024 * a_ped^4 + 0.1305 * a_ped^3 - 2.4822 * a_ped^2 + 18.025 * a_ped + 15.197) * 1000 / 60
    qp_muscle  <- (4.73 - 0.189 * a_ped) / (1 + 0.0672 * a_ped - 0.00515 * a_ped * a_ped) * wp_muscle / 100
    qp_adipose <- 0.024 * wp_adipose
    qp_skin    <- (9.65 + 25.7 * a_ped) / (1 + 1.04 * a_ped + 0.0585 * a_ped * a_ped) * wp_skin / 100
    qp_kidney  <- (229 + 700 * a_ped) / (1 + 1.49 * a_ped + 0.0108 * a_ped * a_ped) * wp_kidney / 100
    qp_spleen  <- (223 + 102 * a_ped) / (1 + 1.03 * a_ped + 0.004 * a_ped * a_ped) * wp_spleen / 100
    qp_stomach <- (128 + 63.9 * a_ped) / (1 + 1.11 * a_ped + 0.00568 * a_ped * a_ped) * wp_stomach / 100
    qp_ha      <- (25.6 + 2.33 * a_ped) / (1 - 0.0227 * a_ped + 0.00508 * a_ped * a_ped) * wp_liver / 100
    qp_gut     <- (216 + 106 * a_ped) / (1 + 1.09 * a_ped + 0.00562 * a_ped * a_ped) * wp_gut / 100

    co_ger  <- (159 * bsa - 1.56 * a_ger + 114) * 1000 / 60                                # Table S7 CO (L/h)
    qg_adipose <- ((0.044 + 0.027 * SEXF) * a_ger + 2.4 * SEXF + 3.9) / 100 * co_ger        # Table S7 Adipose (CO%)
    qg_brain   <- exp(-0.48 * bsa + 0.044 * SEXF + 3.5) / 100 * co_ger                      # Table S7 Brain (CO%)
    qg_heart   <- (-0.72 * ht_ger - 10 * SEXF + 134) / 100 * co_ger                         # Table S7 Heart (CO%)
    qg_kidney  <- (-8.7 * bsa + 0.29 * ht_ger - 0.081 * a_ger - 13) / 100 * co_ger          # Table S7 Kidney (CO%)
    qg_muscle  <- (-6.4 * SEXF + 17.5) / 100 * co_ger                                       # Table S7 Muscle (CO%)
    qg_skin    <- 5 / 100 * co_ger                                                          # Table S7 Skin (CO%)
    qg_stomach <- 38.33                                                                     # Table S7 Stomach, footnote a: assumed equal to adults
    qg_gut     <- (2 * SEXF + 14) / 100 * co_ger                                            # Table S7 Intestine (CO%)
    qg_spleen  <- 3 / 100 * co_ger                                                          # Table S7 Spleen (CO%)
    qg_liver   <- (-0.108 * a_ger + 1.04 * SEXF + 27.9) / 100 * co_ger                      # Table S7 Liver (CO%) = TOTAL hepatic flow
    qg_ha      <- qg_liver - qg_gut - qg_spleen - qg_stomach                                # hepatic artery by difference (Eq S6 defines Ql = Qsp + Qst + Qliv + sum Qgwi)

    co       <- is_ped * co_ped     + is_adl * 5600    + is_ger * co_ger      # Table S2 Lung flow 5600 mL/min = cardiac output
    q_heart  <- is_ped * qp_heart   + is_adl * 240     + is_ger * qg_heart    # Table S2 Heart 240
    q_brain  <- is_ped * qp_brain   + is_adl * 700     + is_ger * qg_brain    # Table S2 Brain 700
    q_muscle <- is_ped * qp_muscle  + is_adl * 750     + is_ger * qg_muscle   # Table S2 Muscle 750
    q_adipos <- is_ped * qp_adipose + is_adl * 260     + is_ger * qg_adipose  # Table S2 Adipose 260
    q_skin   <- is_ped * qp_skin    + is_adl * 300     + is_ger * qg_skin     # Table S2 Skin 300
    q_kidney <- is_ped * qp_kidney  + is_adl * 1240    + is_ger * qg_kidney   # Table S2 Kidneys 1240
    q_spleen <- is_ped * qp_spleen  + is_adl * 80      + is_ger * qg_spleen   # Table S2 Spleen 80
    q_stomac <- is_ped * qp_stomach + is_adl * 38.33   + is_ger * qg_stomach  # Table S2 Stomach 38.33
    q_ha     <- is_ped * qp_ha      + is_adl * 300     + is_ger * qg_ha       # Table S2: Ql 1518.33 - (80 + 38.33 + 1100) = 300 exactly
    q_gut    <- is_ped * qp_gut     + is_adl * 1100    + is_ger * qg_gut      # Table S2 duodenum+jejunum+ileum+cecum+colon = 118+413+244+44+281

    q_gw1 <- 0.107272727 * q_gut   # duodenum 118 / 1100
    q_gw2 <- 0.375454545 * q_gut   # jejunum 413 / 1100
    q_gw3 <- 0.221818182 * q_gut   # ileum 244 / 1100
    q_gw4 <- 0.040000000 * q_gut   # cecum 44 / 1100
    q_gw5 <- 0.255454545 * q_gut   # colon 281 / 1100
    q_liver <- q_ha + q_spleen + q_stomac + q_gut   # Eq S6: Ql = Qsp + Qst + Qliv + sum(Qgwi); adult 300+80+38.33+1100 = 1518.33 = Table S2

    # Rest-of-body flow closes the cardiac-output balance. Adult check:
    # 5600 - (240 + 700 + 750 + 260 + 300 + 1240 + 1518.33) = 591.67, i.e. the
    # Table S2 ROB flow of 592 mL/min.
    q_other0 <- co - (q_heart + q_brain + q_muscle + q_adipos + q_skin + q_kidney + q_liver)
    q_other  <- (q_other0 > 1) * q_other0 + (q_other0 <= 1) * 1

    # =====================================================================
    # 5. CYP3A4 activity, hepatic microsomal protein and intrinsic clearance.
    #    Eq 7  paediatric hepatic CYP3A ratio (infants onward)
    #    Eq 8  neonatal / preterm hepatic CYP3A ratio, in PMA years
    #    Eq 9  paediatric MPPGL
    #    Eq 10 CLint,pediatric = CLint,adult,l * R_CYP3A * MPPGL * LW
    #    Eq 11 paediatric INTESTINAL CYP3A fraction
    #    Eq 18/19 geriatric MPPGL below / above 80 years
    #    Eq 20 geriatric CYP3A ratio, 8 % lost per decade from age 40
    #    Eq 21 CLint,elderly = CLint,adult,l * R_CYP3A * MPPGL * LW
    #    Adult PBSF = 64220 mg = 1690 g liver x 38 mg protein/g (Table S2
    #    footnote); Eq 9 returns 38.2 mg/g at 18 years, so Eq 10 and the
    #    adult anchor agree at the paediatric boundary.
    # =====================================================================
    pma_y   <- PAGE / 12
    is_neo  <- (PNA <= 1)
    r_cyp_neo   <- 0.003 + 0.1331 * pma_y                    # Eq 8
    r_cyp_child <- a_ped^0.83 / (0.31 + a_ped^0.83)          # Eq 7
    r_cyp_ped   <- is_neo * r_cyp_neo + (1 - is_neo) * r_cyp_child
    mppgl_ped   <- 10^(1.407 + 0.0158 * a_ped - 0.00038 * a_ped * a_ped + 0.0000024 * a_ped^3)  # Eq 9
    is_lt80     <- (a_ger < 80)
    mppgl_ger   <- is_lt80 * (0.0001653 * a_ger^3 - 0.02739 * a_ger * a_ger + 1.143 * a_ger + 25.52) +
                   (1 - is_lt80) * (-0.08155 * a_ger * a_ger + 14.48 * a_ger - 612.7)           # Eq 18 / Eq 19
    r_cyp_ger   <- 0.92^((a_ger - 40) / 10)                  # Eq 20: 8 % decline per decade from 40 y

    # PBSF carries its own random effect, so the paediatric and geriatric
    # branches -- which build the microsomal-protein amount from MPPGL x liver
    # weight rather than from the adult PBSF -- are scaled by the same draw.
    # pbsf_scl is 1 at eta = 0, so the typical-value model is unchanged.
    pbsf_scl <- pbsf / 64220
    cl_int_h <- is_ped * (cl_int_l * r_cyp_ped * mppgl_ped * wp_liver * pbsf_scl) +
                is_adl * (cl_int_l * pbsf) +
                is_ger * (cl_int_l * r_cyp_ger * mppgl_ger * wg_liver * pbsf_scl)

    # Intestinal-wall CYP3A4. The adult typical values are derived in ini();
    # Eq 11 scales them down in paediatrics.
    r_gut_cyp <- is_ped * (0.639 * a_ped / (2.36 + a_ped) + 0.42) + is_adl * 1 + is_ger * 1   # Eq 11; Section 2.5 keeps the elderly intestinal CLint at the adult value
    cl_gw1 <- exp(lcl_gw1 + etalcl_gw1) * r_gut_cyp     # duodenum, 9.7 nmol CYP3A
    cl_gw2 <- exp(lcl_gw2 + etalcl_gw2) * r_gut_cyp     # jejunum, 38.4 nmol CYP3A
    cl_gw3 <- exp(lcl_gw3 + etalcl_gw3) * r_gut_cyp     # ileum, 22.4 nmol CYP3A

    # =====================================================================
    # 6. Plasma protein binding.
    #    Eq 13 preterm albumin (PMA years), Eq 14 paediatric albumin,
    #    Eq 22 adult / elderly albumin (stated valid from 20 years),
    #    Eq 15-16 paediatric alpha1-acid glycoprotein, Eq 17 unbound
    #    fraction. Eq 16 as printed returns 53.4 at birth, which cannot be a
    #    ratio; it is read as a PERCENT of the adult level, so
    #    AAG_ped = 0.61 * (0.01137 * days + 53.4) / 100 -- 0.33 g/L in a
    #    neonate rising to the 0.61 g/L adult value, which is the accepted
    #    AAG ontogeny. The adult albumin anchor is Eq 22 at 40 years
    #    (47.12 g/L), 40 years being the reference age Eq 20 uses.
    #    a_log is a_ped floored at 0.01 y so that log() cannot return -Inf
    #    on a term that is then multiplied by zero.
    # =====================================================================
    a_log    <- a_ped * (a_ped > 0.01) + 0.01 * (a_ped <= 0.01)
    alb_ped  <- is_neo * (41.3 * pma_y^2.70 / (0.383^2.70 + pma_y^2.70)) +
                (1 - is_neo) * (1.1287 * log(a_log) + 33.746)              # Eq 13 / Eq 14
    alb_age  <- 51.4 - 0.107 * (is_ger * a_ger + (1 - is_ger) * 40)        # Eq 22
    aag_ped  <- 0.61 * (0.01137 * (PNA * 30.4375) + 53.4) / 100            # Eq 15 + Eq 16
    p_ped    <- bind_alb * alb_ped + (1 - bind_alb) * aag_ped
    p_adl    <- bind_alb * 47.12 + (1 - bind_alb) * 0.61                   # Eq 22 at 40 y; adult AAG 0.61 g/L [ref 45]
    p_ger    <- bind_alb * alb_age + (1 - bind_alb) * 0.58                 # Section 2.5: elderly AAG set to 0.58 g/L
    p_ratio  <- is_ped * (p_ped / p_adl) + is_adl * 1 + is_ger * (p_ger / p_adl)
    fu_p_i   <- 1 / (1 + p_ratio * (1 - fu_p) / fu_p)                      # Eq 17
    # Eq S6 defines fu,b = fu,p * (1 - Hct) / Rb with Hct 42 %, but the
    # deposited Supplementary File S2 script parameterises fub = 0.044, i.e.
    # fu,p itself, and that is the value which reproduces the Table S12
    # predicted AUC (ratio 1.01 vs 1.03). fu,b = fu,p is used here; see the
    # vignette Errata.
    fu_b     <- fu_p_i

    # =====================================================================
    # 7. Absorption. Eq S4: ka,i = 2 * Peff / r_i. Adult radii are the
    #    Table S2 values (2, 1.63, 1.45 cm). Paediatric / geriatric radii
    #    come from Eq 4, r = (0.016 * BSA + 0.0159) / 2, which is written in
    #    METRES and is multiplied by 100 here to match the centimetre radii
    #    of Table S2 and the cm/min units of Peff (Eq 4 returns 2.29 cm at
    #    the adult BSA of 1.87 m2, i.e. the duodenal radius). Eq 4 gives one
    #    radius only, so the jejunal and ileal radii keep their adult ratios
    #    to the duodenum (1.63/2 and 1.45/2).
    #    Gastric emptying and intestinal transit are held at the adult
    #    Table S2 values at every age per Section 2.4 ("the gastric emptying
    #    time and intestinal transit time in pediatrics were assumed to be
    #    similar to those for adults").
    # =====================================================================
    r_duo0 <- (1.6 * bsa + 1.59) / 2                       # Eq 4, m -> cm
    r_duo  <- is_adl * 2 + (1 - is_adl) * r_duo0
    r_jej  <- is_adl * 1.63 + (1 - is_adl) * r_duo0 * 0.815
    r_ile  <- is_adl * 1.45 + (1 - is_adl) * r_duo0 * 0.725
    ka_duo <- 2 * peff / r_duo                             # Eq S4
    ka_jej <- 2 * peff / r_jej                             # Eq S4
    ka_ile <- 2 * peff / r_ile                             # Eq S4
    kt_sto <- 0.08    # Table S2 stomach transit rate constant, 1/min
    kt_duo <- 0.07    # Table S2 duodenum (printed "0..07")
    kt_jej <- 0.03    # Table S2 jejunum
    kt_ile <- 0.04    # Table S2 ileum
    kt_cec <- 0.003   # Table S2 cecum
    kt_col <- 0.001   # Table S2 colon
    fu_gut <- ROUTE_IV * fu_b + (1 - ROUTE_IV) * 1         # Eq S5 surrounding text

    # =====================================================================
    # 8. Tissue concentrations. Every state is an amount in mg; dividing by
    #    the organ volume in mL gives mg/mL, and dividing again by
    #    (Kt:p / Rb) gives the BLOOD concentration leaving that tissue
    #    (Eq S1). Cc converts the venous blood concentration to plasma
    #    (Eq S9 mixes blood) and rescales mg/mL to ng/mL.
    # =====================================================================
    c_lung    <- lung     / v_lung
    c_kidney  <- kidney   / v_kidney
    c_heart   <- heart    / v_heart
    c_liver   <- liver    / v_liver
    c_muscle  <- muscle   / v_muscle
    c_skin    <- skin     / v_skin
    c_adipose <- adipose  / v_adipose
    c_brain   <- brain    / v_brain
    c_venous  <- venous   / v_venous
    c_spleen  <- spleen   / v_spleen
    c_arteria <- arterial / v_arterial
    c_other   <- other    / v_other
    c_gws     <- wall_stomach  / v_stomach
    c_gw1     <- wall_duodenum / v_gw1
    c_gw2     <- wall_jejunum  / v_gw2
    c_gw3     <- wall_ileum    / v_gw3
    c_gw4     <- wall_cecum    / v_gw4
    c_gw5     <- wall_colon    / v_gw5

    kb_lung    <- kp_lung    / bpr
    kb_kidney  <- kp_kidney  / bpr
    kb_heart   <- kp_heart   / bpr
    kb_liver   <- kp_liver   / bpr
    kb_muscle  <- kp_muscle  / bpr
    kb_skin    <- kp_skin    / bpr
    kb_adipose <- kp_adipose / bpr
    kb_brain   <- kp_brain   / bpr
    kb_spleen  <- kp_spleen  / bpr
    kb_stomach <- kp_stomach / bpr
    kb_gut     <- kp_gut     / bpr
    kb_other   <- kp_other   / bpr

    # =====================================================================
    # 9. Mass balance.
    #    Eq S1  perfusion-limited tissue        Eq S2  stomach lumen
    #    Eq S3  intestinal lumen segment        Eq S5  intestinal wall
    #    Eq S6  liver                           Eq S8  arterial blood
    #    Eq S9  mixed venous blood              Eq S10 lung
    # =====================================================================
    d/dt(lung)     = co * c_venous - co * c_lung / kb_lung                              # Eq S10
    d/dt(kidney)   = q_kidney * c_arteria - q_kidney * c_kidney / kb_kidney             # Eq S1
    d/dt(heart)    = q_heart * c_arteria - q_heart * c_heart / kb_heart                 # Eq S1
    d/dt(liver)    = q_stomac * c_gws / kb_stomach +
                     q_gw1 * c_gw1 / kb_gut + q_gw2 * c_gw2 / kb_gut +
                     q_gw3 * c_gw3 / kb_gut + q_gw4 * c_gw4 / kb_gut +
                     q_gw5 * c_gw5 / kb_gut +
                     q_spleen * c_spleen / kb_spleen +
                     q_ha * c_arteria -
                     q_liver * c_liver / kb_liver -
                     cl_int_h * fu_b * c_liver / kb_liver                               # Eq S6
    d/dt(muscle)   = q_muscle * c_arteria - q_muscle * c_muscle / kb_muscle             # Eq S1
    d/dt(skin)     = q_skin * c_arteria - q_skin * c_skin / kb_skin                     # Eq S1
    d/dt(adipose)  = q_adipos * c_arteria - q_adipos * c_adipose / kb_adipose           # Eq S1
    d/dt(brain)    = q_brain * c_arteria - q_brain * c_brain / kb_brain                 # Eq S1
    d/dt(venous)   = q_kidney * c_kidney / kb_kidney + q_heart * c_heart / kb_heart +
                     q_brain * c_brain / kb_brain + q_liver * c_liver / kb_liver +
                     q_muscle * c_muscle / kb_muscle + q_skin * c_skin / kb_skin +
                     q_adipos * c_adipose / kb_adipose + q_other * c_other / kb_other -
                     co * c_venous                                                      # Eq S9
    d/dt(spleen)   = q_spleen * c_arteria - q_spleen * c_spleen / kb_spleen             # Eq S1
    d/dt(arterial) = co * c_lung / kb_lung - co * c_arteria                             # Eq S8
    d/dt(other)    = q_other * c_arteria - q_other * c_other / kb_other                 # Eq S1

    d/dt(stomach)  = -kt_sto * stomach                                                  # Eq S2
    d/dt(duodenum) = kt_sto * stomach  - kt_duo * duodenum - ka_duo * duodenum          # Eq S3
    d/dt(jejunum)  = kt_duo * duodenum - kt_jej * jejunum  - ka_jej * jejunum           # Eq S3
    d/dt(ileum)    = kt_jej * jejunum  - kt_ile * ileum    - ka_ile * ileum             # Eq S3
    d/dt(cecum)    = kt_ile * ileum    - kt_cec * cecum                                 # Eq S3, no absorption below the ileum
    d/dt(colon)    = kt_cec * cecum    - kt_col * colon                                 # Eq S3, no absorption below the ileum

    d/dt(wall_stomach)  = q_stomac * c_arteria - q_stomac * c_gws / kb_stomach          # Eq S1; Eq S2 text: no absorption or metabolism in the stomach
    d/dt(wall_duodenum) = q_gw1 * c_arteria + ka_duo * duodenum -
                          fu_gut * cl_gw1 * c_gw1 / kb_gut - q_gw1 * c_gw1 / kb_gut     # Eq S5
    d/dt(wall_jejunum)  = q_gw2 * c_arteria + ka_jej * jejunum -
                          fu_gut * cl_gw2 * c_gw2 / kb_gut - q_gw2 * c_gw2 / kb_gut     # Eq S5
    d/dt(wall_ileum)    = q_gw3 * c_arteria + ka_ile * ileum -
                          fu_gut * cl_gw3 * c_gw3 / kb_gut - q_gw3 * c_gw3 / kb_gut     # Eq S5
    d/dt(wall_cecum)    = q_gw4 * c_arteria - q_gw4 * c_gw4 / kb_gut                    # Eq S5 without absorption or CYP3A
    d/dt(wall_colon)    = q_gw5 * c_arteria - q_gw5 * c_gw5 / kb_gut                    # Eq S5 without absorption or CYP3A

    # Plasma concentration in ng/mL: blood -> plasma via Rb, mg/mL -> ng/mL
    # via 1e6, matching the ng/mL*min AUC units of Table S14.
    Cc <- c_venous / bpr * 1e6
    Cc ~ prop(propSd)
  })
}
