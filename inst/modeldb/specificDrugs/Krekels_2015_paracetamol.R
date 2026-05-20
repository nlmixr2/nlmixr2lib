Krekels_2015_paracetamol <- function() {
  description <- paste(
    "Parent-and-metabolites population PK model for intravenous paracetamol",
    "(administered as the prodrug propacetamol; doses expressed as paracetamol",
    "equivalents) and its glucuronide and sulphate phase-II conjugates in 54",
    "preterm and term neonates and infants (Krekels 2015). One-compartment",
    "plasma disposition for paracetamol with three parallel elimination",
    "pathways from the central compartment: glucuronide formation (CL_gluc),",
    "sulphate formation (CL_sulf), and unchanged renal excretion (CL_renal).",
    "Each metabolite distributes into a one-compartment plasma space whose",
    "volume is fixed at 18% of the parent volume (Vc_gluc = Vc_sulf =",
    "0.18 * Vc, based on the previously reported adult paracetamol model in",
    "Allegaert et al. and adult literature). The two metabolites share a",
    "common urinary excretion rate constant kE_met = mf * kE_renal, where",
    "kE_renal = CL_renal / Vc is the parent unchanged-renal rate constant",
    "and mf (multiplication factor) is estimated to be 11.3. Cumulative",
    "urinary amounts of parent paracetamol and the two metabolites are",
    "tracked as elimination-amount compartments and exposed as additive-error",
    "observations. Bodyweight enters linearly on Vc (and so on the",
    "metabolite volumes by inheritance), on the glucuronide formation",
    "clearance (CL_gluc), and on the unchanged renal clearance (CL_renal);",
    "the sulphate formation clearance (CL_sulf) scales with bodyweight as",
    "a power with an estimated exponent of 1.40. No postnatal age,",
    "postmenstrual age, sex, term-vs-preterm, or study-protocol covariate",
    "was retained in the final model, and no time-varying (up-regulation)",
    "component was detected on the glucuronidation pathway. Parameter",
    "values reported throughout (mL/min/kg, L/kg) are per-kg quantities;",
    "individual structural parameters are obtained in model() by",
    "multiplying by body weight in kg (linear) or body weight in kg raised",
    "to n (power)."
  )
  reference <- paste(
    "Krekels EHJ, van Ham S, Allegaert K, de Hoon J, Tibboel D,",
    "Danhof M, Knibbe CAJ (2015). Developmental changes rather than",
    "repeated administration drive paracetamol glucuronidation in",
    "neonates and infants.",
    "European Journal of Clinical Pharmacology 71(9):1075-1082.",
    "doi:10.1007/s00228-015-1887-y.",
    sep = " "
  )
  vignette <- "Krekels_2015_paracetamol"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at start of study",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at study start in Krekels 2015 (median bodyweight",
        "2.5 kg, range 0.5-6.3 kg). The per-kg parameterization in Table 1",
        "implies a reference weight of 1 kg, with linear scaling Vc_i =",
        "Vc * WT, CL_gluc_i = CL_gluc * WT, CL_renal_i = CL_renal * WT,",
        "and power scaling CL_sulf_i = CL_sulf * WT^n (n = 1.40). The",
        "metabolite volumes inherit the linear scaling because",
        "Vc_gluc = Vc_sulf = 0.18 * Vc. Source column 'WT' in the",
        "original NONMEM dataset; same orientation, no value",
        "transformation."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 2L,
    age_range      = "PNA 1-140 days; PMA 27-60 weeks",
    age_median     = "PNA 1 day; PMA 36 weeks",
    weight_range   = "0.5-6.3 kg",
    weight_median  = "2.5 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Preterm and term neonates and young infants admitted to the",
      "neonatal intensive care unit (NICU). Paracetamol was administered",
      "as the water-soluble prodrug propacetamol for minor painful",
      "procedures or as additional treatment in neonates on opioid",
      "treatment. Urine collections were available for 22 of the 54",
      "subjects."
    ),
    dose_range     = paste(
      "Single-dose study: 20 or 40 mg/kg propacetamol (paracetamol",
      "equivalent 10 or 20 mg/kg) on the first day of postnatal life,",
      "as a 15-minute IV infusion. Repeated-dose study: 30 mg/kg",
      "propacetamol loading (paracetamol equivalent 15 mg/kg),",
      "followed by 1 to 11 maintenance doses of 20 mg/kg propacetamol",
      "(paracetamol equivalent 10 mg/kg) at PMA-dependent intervals of",
      "q6h to q12h."
    ),
    regions        = "Belgium (Leuven), Netherlands (Rotterdam)",
    notes          = paste(
      "Two studies pooled: a single-dose study (plasma sampling up to",
      "10 h post-dose) and a repeated-dose study with urine collection.",
      "353 paracetamol plasma concentrations and 435 urine observations",
      "of parent and metabolite amounts in 6-, 8-, or 12-hour aliquots",
      "were available. Demographics in Krekels 2015 Methods and Results;",
      "model parameter values from Table 1 'Model fit' column."
    )
  )

  ini({
    # Structural parameter values are the 'Model fit' column of Krekels
    # 2015 Table 1. Per-kg parameterization: each value is the typical
    # parameter at a reference body weight of 1 kg; individual values
    # are derived in model() by multiplying by WT (linear) or WT^n
    # (power). RSE percentages from Table 1 are noted alongside each
    # value and confirm that every structural parameter was estimated
    # with good precision (all < 30%).

    lvc        <- log(1.06)  ; label("Paracetamol central volume at WT = 1 kg (L/kg)")                    # Table 1: V1 = 1.06 L/kg (RSE 4.34%); linear bodyweight covariate
    lcl_gluc   <- log(0.266) ; label("Glucuronide formation clearance at WT = 1 kg (mL/min/kg)")          # Table 1: CL1 = 0.266 mL/min/kg (RSE 17.4%); linear bodyweight covariate
    lcl_sulf   <- log(1.46)  ; label("Sulphate formation clearance at WT = 1 kg (mL/min/kg^n)")           # Table 1: CL2 = 1.46 mL/min/kg^n (RSE 14.8%); power bodyweight covariate with exponent n = 1.40
    lcl_renal  <- log(0.285) ; label("Unchanged paracetamol renal clearance at WT = 1 kg (mL/min/kg)")    # Table 1: CL4 = 0.285 mL/min/kg (RSE 6.98%); linear bodyweight covariate
    e_wt_cl_sulf <- 1.40     ; label("Power exponent for bodyweight on CL_sulf (unitless)")               # Table 1: n = 1.40 (RSE 9.29%)
    lmf        <- log(11.3)  ; label("Log multiplication factor for metabolite urinary excretion rate (unitless)") # Table 1: mf = 11.3 (RSE 21.5%); kE_met = mf * kE_renal where kE_renal = CL_renal / Vc
    frac_vmet  <- fixed(0.18); label("Fraction of Vc used as metabolite volume (unitless)")               # Methods: Vc_gluc = Vc_sulf = 0.18 * Vc, fixed from prior adult paracetamol model (Allegaert 2011 ref. 13; Critchley 2005 ref. 14)

    # Inter-individual variability: P_i = P_p * exp(eta_i) (Methods Eq. 1).
    # Variance values are the omega^2 of the log-normal eta as reported
    # in Table 1 (no transformation needed). The final model retains IIV
    # on Vc, CL_gluc, CL_sulf, and CL_renal only; mf, the WT exponent n,
    # and the residual-error parameters carry no IIV.
    etalvc       ~ 0.0925  # Table 1: IIV variance for V1  = 0.0925 (RSE 29.5%)
    etalcl_gluc  ~ 0.599   # Table 1: IIV variance for CL1 = 0.599  (RSE 41.6%)
    etalcl_sulf  ~ 0.312   # Table 1: IIV variance for CL2 = 0.312  (RSE 33.3%)
    etalcl_renal ~ 0.0879  # Table 1: IIV variance for CL4 = 0.0879 (RSE 52.2%)

    # Residual error: combination (additive + proportional) for the
    # paracetamol plasma observations; additive only for the cumulative
    # urinary amounts of parent paracetamol, paracetamol-glucuronide,
    # and paracetamol-sulphate. Table 1 reports the variance of each
    # error component; nlmixr2 ini() takes the SD, so we convert by
    # taking sqrt() of the published variance.
    addSd          <- sqrt(0.354)  ; label("Paracetamol plasma additive residual SD (mg/L)")                      # Table 1: P plasma additive variance = 0.354 (mg/L)^2 (RSE 29.4%)
    propSd         <- sqrt(0.0198) ; label("Paracetamol plasma proportional residual SD (fraction)")              # Table 1: P plasma proportional variance = 0.0198 (RSE 16.5%)
    addSd_urineP   <- sqrt(0.188)  ; label("Cumulative unchanged paracetamol in urine additive residual SD (mg)") # Table 1: P urine additive variance = 0.188 mg^2 (RSE 27.5%)
    addSd_urinePG  <- sqrt(0.223)  ; label("Cumulative paracetamol-glucuronide in urine additive residual SD (mg)") # Table 1: PG urine additive variance = 0.223 mg^2 (RSE 14.5%)
    addSd_urinePS  <- sqrt(0.332)  ; label("Cumulative paracetamol-sulphate in urine additive residual SD (mg)") # Table 1: PS urine additive variance = 0.332 mg^2 (RSE 35.9%)
  })

  model({
    # Reference bodyweight is 1 kg because Table 1 parameters are
    # reported per kg (linear) or per kg^n (power).
    wt_ref <- 1

    # Individual structural parameters. Bodyweight WT enters linearly
    # on Vc, CL_gluc, CL_renal and as a power with estimated exponent
    # on CL_sulf. The metabolite volumes Vc_gluc = Vc_sulf =
    # frac_vmet * Vc carry the Vc bodyweight scaling by construction.
    # Internal clearance units are L/min (Table 1 reports mL/min/kg, so
    # we divide by 1000 to convert mL -> L while the WT term carries
    # the kg).
    vc       <- exp(lvc        + etalvc)       * (WT / wt_ref)
    cl_gluc  <- exp(lcl_gluc   + etalcl_gluc)  * (WT / wt_ref)                    / 1000
    cl_sulf  <- exp(lcl_sulf   + etalcl_sulf)  * (WT / wt_ref)^e_wt_cl_sulf       / 1000
    cl_renal <- exp(lcl_renal  + etalcl_renal) * (WT / wt_ref)                    / 1000
    mf       <- exp(lmf)

    # Metabolite distribution volumes - fixed fraction of Vc (Methods).
    vc_gluc <- frac_vmet * vc
    vc_sulf <- frac_vmet * vc

    # Micro-constants. Formation rate constants k_form_gluc (P -> PG),
    # k_form_sulf (P -> PS), and the unchanged-paracetamol urinary
    # excretion rate kE_renal (P -> urine) are derived from the parent
    # central compartment. Metabolite urinary excretion rate constants
    # kE_gluc = kE_sulf = mf * kE_renal (Methods text under the
    # structural model: 'the same excretion rate constant (k) was
    # assumed for both metabolites, and this value was estimated to be
    # a multiple of the excretion rate constant of the parent
    # compound').
    k_form_gluc <- cl_gluc  / vc
    k_form_sulf <- cl_sulf  / vc
    kE_renal    <- cl_renal / vc
    kE_gluc     <- mf * kE_renal
    kE_sulf     <- mf * kE_renal

    # ODE system (Krekels 2015 Figure 1 schematic). Three parallel
    # losses from the paracetamol central compartment (k_form_gluc,
    # k_form_sulf, kE_renal) feed the glucuronide central, the sulphate
    # central, and the cumulative urine-of-unchanged-paracetamol
    # compartments. Each metabolite central is then drained at rate
    # mf * kE_renal into its respective cumulative-urine compartment.
    d/dt(central)      <- -(k_form_gluc + k_form_sulf + kE_renal) * central
    d/dt(central_gluc) <-  k_form_gluc * central - kE_gluc * central_gluc
    d/dt(central_sulf) <-  k_form_sulf * central - kE_sulf * central_sulf
    d/dt(urine)        <-  kE_renal    * central
    d/dt(urine_gluc)   <-  kE_gluc     * central_gluc
    d/dt(urine_sulf)   <-  kE_sulf     * central_sulf

    # Observations. Cc is the paracetamol plasma concentration; urineP,
    # urinePG, urinePS are the cumulative amounts (mg paracetamol
    # equivalents) recovered in urine. Krekels 2015 fit per-interval
    # urine amounts; users wanting per-interval amounts difference
    # successive values of urineP / urinePG / urinePS across
    # collection times.
    Cc      <- central / vc
    urineP  <- urine
    urinePG <- urine_gluc
    urinePS <- urine_sulf

    Cc      ~ add(addSd) + prop(propSd)
    urineP  ~ add(addSd_urineP)
    urinePG ~ add(addSd_urinePG)
    urinePS ~ add(addSd_urinePS)
  })
}
