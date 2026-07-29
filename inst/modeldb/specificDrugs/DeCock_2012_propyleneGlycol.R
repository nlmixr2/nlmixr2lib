DeCock_2012_propyleneGlycol <- function() {
  description <- "One-compartment population PK model for intravenous propylene glycol (PG) excipient exposure in preterm and term neonates receiving paracetamol-PG or phenobarbital-PG (De Cock 2012)."
  reference <- paste(
    "De Cock RFW, Knibbe CAJ, Kulo A, de Hoon J, Verbesselt R, Danhof M,",
    "Allegaert K (2013). Developmental pharmacokinetics of propylene glycol",
    "in preterm and term neonates.",
    "British Journal of Clinical Pharmacology 75(1):162-171.",
    "doi:10.1111/j.1365-2125.2012.04312.x.",
    sep = " "
  )
  vignette <- "DeCock_2012_propyleneGlycol"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT_BIRTH = list(
      description        = "Birth weight (time-fixed per subject)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at birth. Drives an allometric power effect on CL relative to a",
        "2.72 kg (= 2720 g) reference: (WT_BIRTH / 2.72)^1.69. The source paper",
        "(Table 2 and abstract) reports the formula in grams: (bBW/2720)^1.69; the",
        "model converts the canonical kg unit at use-site by retaining the kg",
        "reference 2.72 so the dimensionless ratio is identical to the source.",
        "Distinct from WT (current bodyweight, time-varying)."
      ),
      source_name        = "bBW (g; multiply by 1/1000 to obtain canonical WT_BIRTH in kg)"
    ),
    WT = list(
      description        = "Current bodyweight (time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying (weight at day of sampling). Drives an allometric power effect",
        "on V relative to a 2.72 kg reference: (WT / 2.72)^1.45. The source paper",
        "(Table 2 and abstract) reports the formula in grams: (cBW/2720)^1.45; the",
        "model retains the kg reference 2.72 so the dimensionless ratio is identical",
        "to the source. Distinct from WT_BIRTH (birth weight, time-fixed)."
      ),
      source_name        = "cBW (g; multiply by 1/1000 to obtain canonical WT in kg)"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives an allometric power effect on CL relative to a 3-day",
        "reference: (PNA_days / 3)^0.201. The canonical PNA is in months (per",
        "inst/references/covariate-columns.md); the source paper reports the formula",
        "in days, so the model converts internally as PNA_days = PNA_months *",
        "30.4375 and keeps the 3-day reference on the source-paper days scale for",
        "traceability against the published equation. Source range 1-30 days; the",
        "model is not calibrated for PNA < 1 day and supplying PNA = 0 will collapse",
        "the (PNA/3)^0.201 term to 0."
      ),
      source_name        = "PNA (days; multiply by 1/30.4375 to obtain canonical PNA in months)"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital coadministration indicator (1 = subject is receiving phenobarbital-PG; 0 = paracetamol-PG only)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (paracetamol-PG, no phenobarbital coadministration)",
      notes              = paste(
        "Drives a multiplicative effect on V: V is 1.77x higher when PG is",
        "co-administered with phenobarbital compared with paracetamol (De Cock 2012,",
        "Table 2 parameter p and abstract). Encoded log-additively on lvc as",
        "e_conmedpb_vc * CONMED_PB with e_conmedpb_vc = log(1.77). In the source",
        "cohort 25 subjects received phenobarbital-PG only, 34 received",
        "paracetamol-PG only, and 3 received both (Table 1); the dual-treatment",
        "subjects carry CONMED_PB = 1 whenever any phenobarbital-PG is in the",
        "regimen, consistent with the paper's binary stratification."
      ),
      source_name        = "p (Table 2 multiplicative factor on V for phenobarbital-PG vs paracetamol-PG)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 62,
    n_studies      = 1,
    age_range      = "Postnatal age 1-30 days; gestational age 24-41 weeks (Table 1)",
    age_median     = "PNA 3 days (paracetamol-PG) / 2 days (phenobarbital-PG); gestational age 38 weeks (paracetamol) / 34 weeks (phenobarbital) (Table 1)",
    weight_range   = "Birth weight 630-3980 g; current bodyweight 700-4100 g (Table 1)",
    weight_median  = "Birth weight 2990 g (paracetamol) / 1965 g (phenobarbital); pooled cohort median 2720 g (used as the allometric reference in Table 2)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (single-centre Belgian cohort, University Hospitals Leuven)",
    disease_state  = "Preterm and term neonates receiving intravenous paracetamol-PG or phenobarbital-PG for analgesia / sedation; some after perinatal asphyxia",
    dose_range     = paste(
      "IV paracetamol 10 mg/mL containing 800 mg PG per 1000 mg paracetamol (loading",
      "20 mg/kg paracetamol then 5-10 mg/kg q6h; PG load ~16 mg/kg, PG maintenance",
      "~4-8 mg/kg q6h). IV phenobarbital 200 mg/mL containing 700 mg PG per 200 mg",
      "phenobarbital (loading 20 mg/kg phenobarbital then 5 mg/kg/day; PG load",
      "~70 mg/kg, PG maintenance ~17.5 mg/kg/day). See Table 3."
    ),
    regions        = "Belgium (University Hospitals Leuven NICU)",
    notes          = paste(
      "372 PG plasma concentrations from 62 (pre)term neonates; samples 1-11 per",
      "neonate collected 20 min to 20.5 h after dose. Six outliers excluded after",
      "visual inspection of individual chromatographies. Methods 'Patients' section",
      "and Table 1 give the cohort demographics; the 2720 g pooled-median birth",
      "weight is the allometric reference for both CL and V."
    )
  )

  ini({
    # Structural typical values from the De Cock 2012 final pharmacokinetic
    # covariate model (Table 2, also reproduced in the abstract). Reference
    # subject: WT_BIRTH = 2.72 kg (2720 g pooled-cohort median bBW),
    # WT = 2.72 kg, PNA = 3 days, CONMED_PB = 0 (paracetamol-PG only).
    lcl <- log(0.0849); label("Typical CL at reference WT_BIRTH = 2.72 kg, PNA = 3 d (L/h)")  # abstract; Table 2 CLp = 0.085 rounded
    lvc <- log(0.967);  label("Typical V at reference WT = 2.72 kg, CONMED_PB = 0 (L)")        # abstract; Table 2 Vp = 0.97 rounded

    # Allometric exponents (Table 2 parameters m, n, o). All three exponents
    # were estimated, not fixed, per the Table 2 bootstrap CV%.
    e_wtbirth_cl <- 1.69;  label("Allometric exponent on WT_BIRTH for CL (unitless)")  # Table 2 m = 1.69 (CV 10.2%)
    e_pna_cl     <- 0.201; label("Allometric exponent on PNA for CL (unitless)")       # Table 2 n = 0.20 / abstract 0.201 (CV 31.9%)
    e_wt_vc      <- 1.45;  label("Allometric exponent on WT for V (unitless)")          # Table 2 o = 1.45 (CV 10.4%)

    # Phenobarbital-coadministration effect on V (Table 2 parameter p): V is
    # 1.77x higher with phenobarbital-PG vs paracetamol-PG. Encoded log-
    # additively on lvc; e^log(1.77) = 1.77 so CONMED_PB = 1 multiplies V by
    # 1.77 and CONMED_PB = 0 leaves V unchanged.
    e_conmedpb_vc <- log(1.77); label("Log multiplicative effect of phenobarbital coadministration on V (unitless)")  # Table 2 p = 1.77 (CV 12.1%; 95% CI 1.35-2.19)

    # Inter-individual variability. Table 2 reports omega^2 directly on the
    # log-normal scale; final-model values 0.12 (CL) and 0.18 (V) are
    # equivalent to ~36% and ~44% apparent CV. NONMEM $OMEGA was diagonal
    # (no block correlation reported), so etalcl and etalvc are independent.
    etalcl ~ 0.12  # Table 2 omega^2(CL) = 0.12 (CV of estimate 26.3%)
    etalvc ~ 0.18  # Table 2 omega^2(V)  = 0.18 (CV of estimate 25.6%)

    # Residual error: Methods states "proportional error model" with reported
    # sigma^2 = 0.036 (Table 2). In nlmixr2 the proportional residual takes
    # the SD directly, so propSd = sqrt(0.036) ~= 0.190.
    propSd <- sqrt(0.036); label("Proportional residual error (fraction; SD scale)")  # Table 2 sigma^2 = 0.036 (CV 11.8%)
  })

  model({
    # 1. Convert canonical PNA (months) to source-paper days so the 3-day
    #    reference and 0.201 exponent apply on the same scale as Table 2.
    pna_days <- PNA * 30.4375

    # 2. Allometric / covariate factors (Table 2 final-model equations).
    #    The kg reference 2.72 is numerically equivalent to the source's
    #    2720 g median because both ratios cancel units.
    cl_size <- (WT_BIRTH / 2.72)^e_wtbirth_cl * (pna_days / 3)^e_pna_cl
    vc_size <- (WT       / 2.72)^e_wt_vc

    # 3. Individual PK parameters with covariate effects and IIV
    cl <- exp(lcl + etalcl) * cl_size
    vc <- exp(lvc + etalvc + e_conmedpb_vc * CONMED_PB) * vc_size

    # 4. Micro-constant for the 1-compartment ODE
    kel <- cl / vc

    # 5. ODE system: IV doses go directly to central. PG enters as the
    #    excipient mass within paracetamol-PG or phenobarbital-PG infusions.
    d/dt(central) <- -kel * central

    # 6. Observation (PG plasma concentration; dose in mg, V in L)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
