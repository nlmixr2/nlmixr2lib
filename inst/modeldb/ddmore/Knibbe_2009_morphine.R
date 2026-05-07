# Joint parent-metabolite population PK model for morphine and its two
# glucuronide metabolites (M3G, M6G) in preterm neonates, infants, and
# children younger than 3 years, extracted from DDMORE Foundation Model
# Repository entry DDMODEL00000248 (Knibbe et al. 2009, Clin
# Pharmacokinet 48(6):371-385, doi:10.2165/00003088-200948060-00003).

Knibbe_2009_morphine <- function() {
  description <- "Joint parent-metabolite population PK model for morphine and its glucuronide metabolites M3G and M6G in preterm neonates, infants and toddlers <3 years, with bodyweight allometric scaling and a postnatal-age-stratified glucuronidation step at PNA = 10 days"
  reference <- paste(
    "Knibbe CAJ, Krekels EHJ, van den Anker JN, DeJongh J,",
    "Santen GWE, van Dijk M, Simons SHP, van Lingen RA,",
    "Jacqz-Aigrain EM, Danhof M, Tibboel D (2009).",
    "Morphine glucuronidation in preterm neonates, infants and",
    "children younger than 3 years.",
    "Clin Pharmacokinet 48(6):371-385.",
    "doi:10.2165/00003088-200948060-00003.",
    "DDMORE Foundation Model Repository: DDMODEL00000248.",
    sep = " "
  )
  vignette <- "Knibbe_2009_morphine"
  ddmore_id <- "DDMODEL00000248"
  replicate_of <- NULL
  units <- list(
    time = "min",
    dosing = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight, time-varying",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying body weight at the time of each dose / observation.",
        "Source NONMEM column BWS reports body weight in grams; convert to canonical WT (kg) via WT = BWS / 1000.",
        "The source model uses the un-normalised expression WT^allo_cl on each clearance and WT^allo_v on each volume,",
        "i.e. an implicit reference weight of 1 kg (the typical-value parameters are the per-kg^exponent constants)."
      ),
      source_name        = "BWS"
    ),
    PNA = list(
      description        = "Postnatal age, time-varying",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying postnatal age. Source NONMEM column PNA reports age in days; canonical unit is months.",
        "The Knibbe 2009 maturation step is at PNA = 10 days = 10/30.4375 months ~= 0.32855 months:",
        "morphine -> M3G and morphine -> M6G formation clearances increase from THETA(1)/THETA(6)",
        "(PNA <= 10 days) to THETA(10)/THETA(11) (PNA > 10 days)."
      ),
      source_name        = "PNA"
    )
  )

  population <- list(
    n_subjects     = 248,
    n_studies      = "pooled (multi-study) cohort; exact study count not stated in the available abstract",
    age_range      = "preterm newborns to <3 years (postnatal age range covers neonatal day 0 through ~36 months)",
    age_median     = "TODO: full Knibbe 2009 publication not on disk during DDMORE extraction",
    weight_range   = "TODO: see age_median note",
    weight_median  = "TODO: see age_median note",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state  = "Postoperative neonates / infants / toddlers (preterm + term) receiving IV morphine for analgesia",
    dose_range     = "IV bolus + continuous infusion; doses and infusion rates not captured (full text not on disk)",
    regions        = "Pooled European paediatric cohorts (ICU and postoperative settings; original studies span the Netherlands and France per author affiliations)",
    notes          = paste(
      "Per the publication abstract (PMID 19650676): 248 infants contributing 2,159 morphine concentrations.",
      "The DDMORE-shipped Output_real_run4.lst was fit to a pooled `Combined_InternalExternalData.csv` covering",
      "338 individuals / 2,809 observations / 5,302 records (likely an extended-cohort post-publication re-run).",
      "Demographic detail fields are TODO because the full Knibbe 2009 publication PDF was not on disk during",
      "extraction; the abstract reports only headline numbers (n = 248, BW exponent on CL = 1.44, PNA cutoff = 10 d)."
    )
  )

  ini({
    # All structural-parameter values come from the FINAL PARAMETER ESTIMATE
    # block of the DDMORE-shipped Output_real_run4.lst (status:
    # MINIMIZATION SUCCESSFUL, NSIG=3.0). The `.mod` $THETA / $OMEGA /
    # $SIGMA blocks list the NONMEM initial values, not the final
    # estimates; per ddmore-source.md final values come from the `.lst`.

    # Morphine central volume V1 (NONMEM CMT 1, the dosed compartment).
    lvc <- log(1.99)
    label("Morphine central volume V1 per kg^allo_v (L)")
    # Output_real_run4.lst, FINAL PARAMETER ESTIMATE: TH 2 = 1.99E+00.

    # Morphine inter-compartmental clearance Q (CMT1 <-> CMT4).
    # Source `.mod`: Q = THETA(3) -- note the absence of WT scaling on Q.
    lq <- log(0.0289)
    label("Morphine inter-compartmental clearance Q (L/min, no WT scaling)")
    # Output_real_run4.lst FINAL: TH 3 = 2.89E-02.

    # Volume of M3G and M6G compartments expressed as a fraction of
    # morphine V1 (V2 = V3 = THETA(8) * V1; V4 = V1).
    lfvm <- log(0.119)
    label("Metabolite distribution volume as fraction of morphine V1 (unitless)")
    # Output_real_run4.lst FINAL: TH 8 = 1.19E-01.

    # Morphine-to-M3G formation clearance (NONMEM CL2). PNA-stratified:
    # `_le10` for postnatal age <= 10 days, `_gt10` for PNA > 10 days.
    # Units of the bare typical-value: L/min per kg^allo_cl.
    lcl_form_m3g_le10 <- log(0.00309)
    label("Morphine->M3G formation clearance per kg^allo_cl, PNA <= 10 days (L/min)")
    # Output_real_run4.lst FINAL: TH 1 = 3.09E-03.
    lcl_form_m3g_gt10 <- log(0.00825)
    label("Morphine->M3G formation clearance per kg^allo_cl, PNA > 10 days (L/min)")
    # Output_real_run4.lst FINAL: TH10 = 8.25E-03.

    # Morphine-to-M6G formation clearance (NONMEM CL3). PNA-stratified.
    lcl_form_m6g_le10 <- log(0.000408)
    label("Morphine->M6G formation clearance per kg^allo_cl, PNA <= 10 days (L/min)")
    # Output_real_run4.lst FINAL: TH 6 = 4.08E-04.
    lcl_form_m6g_gt10 <- log(0.000699)
    label("Morphine->M6G formation clearance per kg^allo_cl, PNA > 10 days (L/min)")
    # Output_real_run4.lst FINAL: TH11 = 6.99E-04.

    # M3G and M6G excretion clearances (NONMEM CL4 and CL5).
    lcl_m3g <- log(0.00219)
    label("M3G elimination (excretion) clearance per kg^allo_cl (L/min)")
    # Output_real_run4.lst FINAL: TH 7 = 2.19E-03.
    lcl_m6g <- log(0.00111)
    label("M6G elimination (excretion) clearance per kg^allo_cl (L/min)")
    # Output_real_run4.lst FINAL: TH 9 = 1.11E-03.

    # Allometric exponents on body weight.
    allo_cl <- 1.44
    label("Allometric exponent of body weight on clearances (unitless)")
    # Output_real_run4.lst FINAL: TH 4 = 1.44E+00.
    allo_v <- fix(1.00)
    label("Allometric exponent of body weight on volumes (unitless, FIXED in source)")
    # Output_real_run4.lst FINAL: TH 5 = 1.00E+00; source `.mod` $THETA flag `1 FIX`.

    # IIV. The source `.mod` $OMEGA structure is:
    #   $OMEGA          0.104                         ; OM1 = ETA(1) on CL2
    #                   0.230                         ; OM2 = ETA(2) on V1
    #   $OMEGA BLOCK(2) 0.185                         ; OM3 = ETA(3) on CL5 (M6G excretion)
    #                   0.178   0.258                 ; cov, OM4 = ETA(4) on CL4 (M3G excretion)
    # Note that the BLOCK pairs ETA(3) and ETA(4); .mod ordering is
    # row-by-row lower triangle (var3, cov34, var4), preserved here.
    # No IIV on Q (THETA(3)) or on the morphine->M6G formation pathway.
    etalcl_form_m3g ~ 0.104
    # Output_real_run4.lst FINAL OMEGA: ETA1|ETA1 = 1.04E-01.
    etalvc ~ 0.230
    # Output_real_run4.lst FINAL OMEGA: ETA2|ETA2 = 2.30E-01.
    etalcl_m6g + etalcl_m3g ~ c(0.185, 0.178, 0.258)
    # Output_real_run4.lst FINAL OMEGA block:
    #   ETA3|ETA3 = 1.85E-01, ETA4|ETA3 = 1.78E-01, ETA4|ETA4 = 2.58E-01.

    # Residual error. The source .mod $ERROR block applies additive errors
    # on the natural-log of the typical prediction:
    #   IPRE = LOG(F)
    #   Y1 = IPRE + ERR(1)        # morphine
    #   Y2 = IPRE + ERR(2)        # M3G
    #   Y3 = IPRE + ERR(3)        # M6G
    # Per naming-conventions.md `$ERROR block patterns`,
    # NONMEM "additive on log-scale" maps to proportional residual error
    # in nlmixr2's linear concentration space. For small SIGMA, the
    # proportional SD equals sqrt(SIGMA). The fourth $SIGMA slot in the
    # .mod (NKOD-conditional 24h additive error) is run-summary-only; the
    # `Output_real_run4.lst` $SIGMA SIG4 = 9.36E+00 corresponds to the
    # `Y_i = IPRE + ERR(i) + TEH * ERR(4)` term that is gated by an
    # NKOD/TIME>1900 indicator (see source `.mod`); it is dropped here
    # because the gate column NKOD is study-specific to the `Combined`
    # dataset and not part of the canonical-public model.
    propSd <- sqrt(0.371)
    label("Proportional residual SD for morphine plasma concentration (fraction)")
    # Output_real_run4.lst FINAL SIGMA: EPS1|EPS1 = 3.71E-01; propSd = sqrt(SIGMA).
    propSd_m3g <- sqrt(0.206)
    label("Proportional residual SD for M3G plasma concentration (fraction)")
    # Output_real_run4.lst FINAL SIGMA: EPS2|EPS2 = 2.06E-01; propSd = sqrt(SIGMA).
    propSd_m6g <- sqrt(0.0967)
    label("Proportional residual SD for M6G plasma concentration (fraction)")
    # Output_real_run4.lst FINAL SIGMA: EPS3|EPS3 = 9.67E-02; propSd = sqrt(SIGMA).
  })

  model({
    # Postnatal-age threshold expressed in canonical PNA months.
    # Source `.mod` uses days (PNA > 10); 10 days / (365.25/12) months = 10/30.4375.
    pna_threshold <- 10 / 30.4375

    # Body-weight scaling factors. allo_v is fixed at 1.0 in the source
    # (source NONMEM `$THETA 1 FIX`); allo_cl is estimated.
    bw_cl_factor <- WT^allo_cl
    bw_v_factor  <- WT^allo_v

    # Individual structural parameters.
    # Morphine central volume V1.
    vc <- exp(lvc + etalvc) * bw_v_factor
    # Morphine peripheral volume V4 = V1 per source (`V4 = V1`).
    vp <- vc
    # M3G and M6G compartment volumes V2 = V3 = TH8 * V1 per source.
    fvm <- exp(lfvm)
    vc_m3g <- fvm * vc
    vc_m6g <- fvm * vc

    # Morphine inter-compartmental clearance Q (no BW scaling per source).
    q <- exp(lq)

    # PNA-stratified typical-value formation clearances. The source `.mod`
    # selects between the two THETAs via `IF (PNA.GT.10) CL = THETA(...)`,
    # i.e. a hard step at PNA = 10 days. The single ETA on CL2 is shared
    # across PNA strata and is applied multiplicatively after the typical
    # value is selected (source: `CL2 = ... * EXP(ETA(1))` in both branches).
    typ_cl_form_m3g <-
      ifelse(PNA <= pna_threshold,
             exp(lcl_form_m3g_le10),
             exp(lcl_form_m3g_gt10))
    typ_cl_form_m6g <-
      ifelse(PNA <= pna_threshold,
             exp(lcl_form_m6g_le10),
             exp(lcl_form_m6g_gt10))

    cl_form_m3g <- typ_cl_form_m3g * exp(etalcl_form_m3g) * bw_cl_factor
    cl_form_m6g <- typ_cl_form_m6g * bw_cl_factor   # source has no IIV on CL3
    cl_m3g <- exp(lcl_m3g + etalcl_m3g) * bw_cl_factor
    cl_m6g <- exp(lcl_m6g + etalcl_m6g) * bw_cl_factor

    # ODE system. Compartments mapped from source `.mod` $MODEL:
    #   `central`     = NONMEM CMT 1 (CENTRAL morphine, DEFDOSE)
    #   `peripheral1` = NONMEM CMT 4 (PERIPH morphine)
    #   `central_m3g` = NONMEM CMT 2 (M3G compartment)
    #   `central_m6g` = NONMEM CMT 3 (M6G compartment)
    # Doses arrive on `central` (the compartment named first in the d/dt
    # ordering). The source NONMEM `.mod` does not include a unit-of-mass
    # conversion between morphine and its glucuronides; the formation
    # flux moves morphine-mass-equivalents into the metabolite
    # compartments and the clearance/volume estimates absorb the
    # underlying molecular-weight ratio implicitly. The same convention
    # is preserved here so the parameter values match the source `.lst`.
    d/dt(central) <- q * peripheral1 / vp -
                     (q + cl_form_m3g + cl_form_m6g) * central / vc
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(central_m3g) <- cl_form_m3g * central / vc -
                         cl_m3g * central_m3g / vc_m3g
    d/dt(central_m6g) <- cl_form_m6g * central / vc -
                         cl_m6g * central_m6g / vc_m6g

    # Plasma concentrations. With AMT in micrograms and V in litres, the
    # raw concentration is in micrograms per litre, which equals
    # nanograms per millilitre (1 ug/L = 1 ng/mL). No further unit
    # conversion is applied here; `units$concentration = "ng/mL"`.
    Cc <- central / vc
    Cc_m3g <- central_m3g / vc_m3g
    Cc_m6g <- central_m6g / vc_m6g

    Cc ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
    Cc_m6g ~ prop(propSd_m6g)
  })
}
