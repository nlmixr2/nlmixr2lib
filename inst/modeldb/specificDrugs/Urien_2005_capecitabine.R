Urien_2005_capecitabine <- function() {
  description <- paste(
    "Population PK model for oral capecitabine and its three sequential",
    "metabolites 5'-DFCR (5'-deoxy-5-fluorocytidine), 5'-DFUR",
    "(5'-deoxy-5-fluorouridine), and 5-FU (5-fluorouracil) in 40 adult",
    "patients with metastatic cancer (Urien 2005). Capecitabine PK is a",
    "one-compartment apparent V1/F model with first-order absorption (Ka)",
    "and a lag time; non-transformation elimination CL10/F runs in parallel",
    "with the formation clearance CL12/F to 5'-DFCR. Each metabolite has",
    "its own central compartment with apparent volume fixed to 1 L (only",
    "output rate constants are identifiable in the source NONMEM ADVAN6",
    "fit), so the chain 5'-DFCR -> 5'-DFUR -> 5-FU -> output is described",
    "by first-order rate constants K23, K34, K40 (paper's notation).",
    "Total bilirubin (canonical TBILI; source column BILT, umol/L) is the",
    "only retained covariate: power exponent +0.32 on CL10/F and -0.36 on",
    "K34, both centred on the median bilirubin 8.8 umol/L. Inter-individual",
    "variability is reported on TLAG, V1, CL10, K23, K34, and K40; ISV on",
    "CL12 was fixed to 0 and ISV on Ka was deleted in favour of a large",
    "inter-occasion variability on Ka that is not represented in this",
    "static model file (see vignette Errata)."
  )
  reference <- "Urien S, Rezai K, Lokiec F. Pharmacokinetic modelling of 5-FU production from capecitabine--a population study in 40 adult patients with metastatic cancer. J Pharmacokinet Pharmacodyn. 2005;32(5-6):817-833. doi:10.1007/s10928-005-0018-2"
  vignette  <- "Urien_2005_capecitabine"
  units     <- list(
    time          = "hour",
    dosing        = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    TBILI = list(
      description        = "Total serum bilirubin concentration (baseline; constant within an individual in the source dataset).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling centred on the median bilirubin of 8.8 umol/L: (TBILI/8.8)^+0.32 on capecitabine non-transformation CL10/F and (TBILI/8.8)^-0.36 on the 5'-DFUR -> 5-FU rate constant K34. SI units in the source paper (Table I lists median 8.8 umol/L, range 3-22 umol/L); convert mg/dL to umol/L by multiplying by 17.1 if the simulation dataset reports US units.",
      source_name        = "BILT"
    )
  )

  population <- list(
    n_subjects     = 40L,
    n_studies      = 2L,
    age_range      = "30-73 years",
    age_median     = "54.5 years",
    weight_range   = "41.5-95 kg",
    weight_median  = "68 kg",
    height_range   = "150-178 cm",
    height_median  = "169 cm",
    bsa_range      = "1.40-2.10 m^2",
    bsa_median     = "1.80 m^2",
    sex_female_pct = 37.5,
    disease_state  = "Adult patients with metastatic cancer receiving second- or third-line chemotherapy. Capecitabine was combined with either irinotecan 200-250 mg/m^2 (90-min infusion on day 1) or irofulven 0.4 mg/kg (30-min infusion on day 1), depending on the phase I study.",
    dose_range     = "Oral capecitabine 1400-2300 mg/m^2/day (median 2000 mg/m^2/day; total daily dose 1900 mg/m^2 on average), divided into twice-daily doses every 12 hours within 30 minutes of breakfast or dinner.",
    regions        = "France (Centre Rene Huguenin, Saint-Cloud).",
    n_observations = "1426 plasma concentrations across 75 PK courses (40 patients): 373 capecitabine, 354 5'-DFCR, 363 5'-DFUR, and 336 5-FU. Limits of quantification 0.03 / 0.05 / 0.05 / 0.05 umol/L respectively. Most subjects had two PK evaluations on days 1 and 15.",
    baseline_chem  = "Serum albumin median 37 g/L (range 30-46), serum creatinine median 85 umol/L (range 58-113), total bilirubin median 8.8 umol/L (range 3-22).",
    notes          = "Demographics from Urien 2005 Table I. Modelling software NONMEM V level 1.1, FO method (FOCE produced abnormal terminations on this complex four-compartment chain). Final estimates are the bootstrap means in Table II from 1922 successful runs of 2000 programmed runs."
  )

  ini({
    # Capecitabine absorption (first-order with lag time). Ka and TLAG
    # estimates are bootstrap means from Urien 2005 Table II. ISV on Ka
    # was deleted in favour of a large IOV (167%); IOV is not represented
    # in this static file -- see vignette Errata.
    lka  <- log(2.07)
    label("Capecitabine first-order absorption rate constant (1/h)")  # Urien 2005 Table II: Ka mean 2.07 1/h (median 2.29; 95% CI 0.91-2.66)
    ltlag <- log(0.28)
    label("Capecitabine absorption lag time TLAG (h)")  # Urien 2005 Table II: TLAG mean 0.28 h (median 0.16; 95% CI 0.01-0.36)

    # Capecitabine apparent central volume V1/F. Bootstrap mean 338 L
    # (median 290; 95% CI 206-391).
    lvc <- log(338)
    label("Capecitabine apparent central volume V1/F (L)")  # Urien 2005 Table II: V1 mean 338 L

    # Capecitabine non-transformation apparent clearance CL10/F (the
    # arm that does NOT proceed to 5'-DFCR). Typical value 218 L/h
    # holds at the median bilirubin 8.8 umol/L; covariate effect is
    # multiplicative power.
    lcl <- log(218)
    label("Capecitabine non-transformation apparent clearance CL10/F at TBILI = 8.8 umol/L (L/h)")  # Urien 2005 Table II: TV.CL10 mean 218 L/h; covariate equation page 826: CL10 = 218 * (BILT/8.8)^+0.32
    e_tbili_cl <- 0.32
    label("Power exponent of (TBILI/8.8) on capecitabine non-transformation CL10/F (unitless)")  # Urien 2005 Table II: BILT effect on CL10 mean +0.32 (median +0.24; 95% CI +0.10 to +0.44); page 826 equation

    # Capecitabine -> 5'-DFCR formation apparent clearance CL12/F. ISV
    # on CL12 was fixed to 0 in the source (no detrimental effect on the
    # OFV), so no IIV on this parameter.
    lcl_dfcr <- log(12.9)
    label("Capecitabine -> 5'-DFCR formation apparent clearance CL12/F (L/h)")  # Urien 2005 Table II: CL12 mean 12.9 L/h (median 10.4; 95% CI 7.25-20)

    # Sequential metabolite chain. Apparent metabolite volumes are not
    # identifiable and are FIXED to 1 L in the paper, so the published
    # rate constants K23, K34, K40 are output rate constants in the
    # ADVAN6 sense and act directly on the metabolite-compartment
    # amount. The naming `lk_<destination>_form` matches the parallel
    # NA_NA_lidocaine.R sequential-metabolite model (k_megx_form etc.).
    lvc_dfcr <- fixed(log(1))
    label("5'-DFCR apparent central volume (L) -- FIXED to 1 per Urien 2005 (not identifiable)")  # Urien 2005 page 820: "metabolite distribution volumes are not identifiable and fixed to 1"
    lvc_dfur <- fixed(log(1))
    label("5'-DFUR apparent central volume (L) -- FIXED to 1 per Urien 2005 (not identifiable)")  # Urien 2005 page 820
    lvc_5fu  <- fixed(log(1))
    label("5-FU apparent central volume (L) -- FIXED to 1 per Urien 2005 (not identifiable)")    # Urien 2005 page 820

    lk_dfur_form <- log(10.7)
    label("5'-DFCR -> 5'-DFUR first-order rate constant K23 (1/h)")  # Urien 2005 Table II: K23 mean 10.7 1/h (median 7.9; 95% CI 5.8-14.3)

    # K34 typical value 5.30 1/h at the median bilirubin 8.8 umol/L.
    # Note: the paper text on page 826 reports the typical value in the
    # covariate equation as 5.70 1/h, while Table II reports the bootstrap
    # mean as 5.30 1/h with bootstrap median 4.30 (95% CI 3.0-7.60). This
    # file uses 5.30 (Table II) per the standing rule that the official
    # bootstrap-summary table is authoritative; see vignette Errata for
    # the discrepancy.
    lk_5fu_form <- log(5.30)
    label("5'-DFUR -> 5-FU first-order rate constant K34 at TBILI = 8.8 umol/L (1/h)")  # Urien 2005 Table II: K34 mean 5.30 1/h (page 826 covariate equation reports 5.70 -- see vignette Errata)
    e_tbili_k_5fu_form <- -0.36
    label("Power exponent of (TBILI/8.8) on K34 (unitless)")  # Urien 2005 Table II: BILT effect on K34 mean -0.36 (median -0.42; 95% CI -0.68 to -0.20); page 826 equation

    lkel_5fu <- log(66)
    label("5-FU first-order elimination rate constant K40 (1/h)")  # Urien 2005 Table II: K40 mean 66 1/h (median 59; 95% CI 38-117)

    # Inter-individual variability. Source paper reports ISV as %CV
    # (Table II); convert to log-normal omega^2 via omega^2 = log(1 +
    # CV^2). ISV on Ka was deleted in favour of IOV; ISV on CL12 was
    # fixed to 0 -- no eta on lka or lcl_dfcr.
    etaltlag         ~ 0.79299  # Urien 2005 Table II: ISV TLAG 110% -> log(1 + 1.10^2) = 0.79299
    etalvc          ~ 1.04706  # Urien 2005 Table II: ISV V1 136% -> log(1 + 1.36^2) = 1.04706
    etalcl          ~ 0.03187  # Urien 2005 Table II: ISV CL10 18% -> log(1 + 0.18^2) = 0.03187
    etalk_dfur_form ~ 0.22314  # Urien 2005 Table II: ISV K23 50% -> log(1 + 0.50^2) = 0.22314
    etalk_5fu_form  ~ 0.06062  # Urien 2005 Table II: ISV K34 25% -> log(1 + 0.25^2) = 0.06062
    etalkel_5fu     ~ 0.10947  # Urien 2005 Table II: ISV K40 34% -> log(1 + 0.34^2) = 0.10947

    # Residual error. Source reports residual variability in umol/L
    # (Table II "Res. variability" rows); the values are additive SDs
    # in the linear-concentration space. Bootstrap covariances between
    # residual-error pairs (cov(eps_cap, eps_dfcr) = 3.1; cov(eps_dfur,
    # eps_5fu) = 2.6) are NOT propagated into this independent-residual
    # parameterisation -- see vignette Errata.
    addSd      <- 3.83
    label("Capecitabine additive residual SD (umol/L)")  # Urien 2005 Table II: Res. variability CAP mean 3.83 umol/L (median 4.10; 95% CI 3.1-4.9)
    addSd_dfcr <- 3.72
    label("5'-DFCR additive residual SD (umol/L)")        # Urien 2005 Table II: Res. variability 5'-DFCR mean 3.72 umol/L (median 3.7; 95% CI 2.7-4.8)
    addSd_dfur <- 5.81
    label("5'-DFUR additive residual SD (umol/L)")        # Urien 2005 Table II: Res. variability 5'-DFUR mean 5.81 umol/L (median 6.6; 95% CI 5.2-7.8)
    addSd_5fu  <- 0.64
    label("5-FU additive residual SD (umol/L)")           # Urien 2005 Table II: Res. variability 5-FU mean 0.64 umol/L (median 0.67; 95% CI 0.40-0.92)
  })

  model({
    # 1. Individual capecitabine PK parameters with bilirubin power
    # scaling on CL10/F (page 826 equation).
    ka       <- exp(lka)
    tlag     <- exp(ltlag + etaltlag)
    vc       <- exp(lvc + etalvc)
    cl       <- exp(lcl + etalcl) * (TBILI / 8.8)^e_tbili_cl
    cl_dfcr  <- exp(lcl_dfcr)

    # 2. Metabolite apparent volumes (FIXED to 1 L in the source).
    vc_dfcr <- exp(lvc_dfcr)
    vc_dfur <- exp(lvc_dfur)
    vc_5fu  <- exp(lvc_5fu)

    # 3. Metabolite first-order rate constants with bilirubin power
    # scaling on K34 (page 826 equation).
    k_dfur_form <- exp(lk_dfur_form + etalk_dfur_form)
    k_5fu_form  <- exp(lk_5fu_form  + etalk_5fu_form) *
                     (TBILI / 8.8)^e_tbili_k_5fu_form
    kel_5fu     <- exp(lkel_5fu + etalkel_5fu)

    # 4. Capecitabine micro-constants. Total apparent capecitabine
    # elimination is the sum of the non-transformation (CL10) and the
    # 5'-DFCR-formation (CL12) pathways; each is a separate flux.
    kel       <- cl       / vc
    k_dfcr_form <- cl_dfcr / vc

    # 5. ODE system (Urien 2005 Appendix A, equations A.1-A.4, with the
    # explicit depot equivalent of the analytical "Ka * D" input).
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * central - k_dfcr_form * central
    d/dt(central_dfcr) <-  k_dfcr_form * central      - k_dfur_form * central_dfcr
    d/dt(central_dfur) <-  k_dfur_form * central_dfcr - k_5fu_form  * central_dfur
    d/dt(central_5fu)  <-  k_5fu_form  * central_dfur - kel_5fu     * central_5fu

    # 6. Absorption lag time on the depot (Urien 2005 Appendix A
    # implicitly via TLAG; explicit lag() is the rxode2 equivalent).
    lag(depot) <- tlag

    # 7. Observation variables. Concentrations are amount divided by
    # apparent volume; metabolite volumes equal 1 L by fixed parameter
    # so the metabolite "concentration" is numerically equal to the
    # compartment amount in umol.
    Cc       <- central      / vc
    Cc_dfcr  <- central_dfcr / vc_dfcr
    Cc_dfur  <- central_dfur / vc_dfur
    Cc_5fu   <- central_5fu  / vc_5fu

    Cc      ~ add(addSd)
    Cc_dfcr ~ add(addSd_dfcr)
    Cc_dfur ~ add(addSd_dfur)
    Cc_5fu  ~ add(addSd_5fu)
  })
}
