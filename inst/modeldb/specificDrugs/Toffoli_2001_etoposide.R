Toffoli_2001_etoposide <- function() {
  description <- "Two-compartment population PK model for oral and IV etoposide in adult patients with solid tumours (Toffoli 2001). Additive-linear creatinine-clearance covariate on CL from Equation 2."
  reference <- "Toffoli G, Corona G, Sorio R, Robieux I, Basso B, Colussi AM, Boiocchi M. Population pharmacokinetics and pharmacodynamics of oral etoposide. Br J Clin Pharmacol. 2001;52(5):511-519. doi:10.1046/j.0306-5251.2001.01468.x"
  vignette <- "Toffoli_2001_etoposide"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Jelliffe (1973) method, reported as absolute mL/min (NOT BSA-normalized). Precedent for a raw-mL/min value under the CRCL canonical is Delattre_2010_amikacin (Cockcroft-Gault).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (baseline). Median 70 mL/min, range 22-121 in the Toffoli 2001 cohort (Table 1). Enters CL via the additive linear equation CL (L/h) = 0.74 + 0.0057 * CRCL (Toffoli 2001 Equation 2, page 4). The intercept-slope form is preserved rather than converted to a proportional form to remain faithful to the source. The Eq. 2 slope 0.0057 is calibrated to raw mL/min values (not BSA-normalized); do NOT feed BSA-normalized eGFR (mL/min/1.73 m^2) into this model without rescaling.",
      source_name        = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported as significant on Vc and k12 in Table 2 (final population model) but no covariate-effect coefficient is given in the paper text or in Equations 2-5. Cannot be encoded in the model without the coefficient; documented here for provenance.",
      source_name        = "BSA"
    ),
    FU = list(
      description        = "Unbound fraction of etoposide in plasma (100% - percent protein binding). Per-subject value determined by ultrafiltration of the Cmax sample.",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported as significant on Vc and k12 (via %PB) in Table 2 but no covariate-effect coefficient is given. Population mean protein binding 91.5% (fu = 0.085), CV 5%, range 79-98% (i.e., fu range 0.02-0.21). Documented here for provenance.",
      source_name        = "%PB"
    ),
    TBILI = list(
      description        = "Serum total bilirubin",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a candidate covariate on CL (Table 2 'Covariables' column) but did not further reduce residual variability once CRCL was in the model; not retained in the final CL relationship (page 4, 'Influence of hepatic function on clearance').",
      source_name        = "bili"
    ),
    AGE = list(
      description        = "Age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Correlates with CL (via correlation with CRCL, Eq. 3) and with free AUCp.o. (Eq. 4) but not retained as an independent covariate after CRCL is accounted for (stepwise regression; page 5).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50,
    n_studies      = 1,
    age_range      = "50-83 years",
    age_median     = "65 years (pooled across tumour types)",
    weight_range   = NULL,
    weight_median  = NULL,
    sex_female_pct = 26,
    race_ethnicity = NULL,
    disease_state  = "Advanced solid tumours (hepatocellular carcinoma n=17, non-small-cell lung cancer n=19, gastric cancer n=5, breast cancer n=6, other n=3)",
    dose_range     = "100 mg oral soft-gelatin capsule (Vepesid or Lastet) daily for 14 days every 3 weeks. One oral dose replaced by 50 mg 1-h IV infusion on day 1 or day 7 (crossover; day randomized).",
    regions        = "Italy (Aviano)",
    notes          = "Prospective open-label bioavailability study, recruitment 1994-1998, single-center (Centro di Riferimento Oncologico, Aviano). Baseline creatinine clearance median 70 mL/min (per Table 1 medians across tumour groups, Jelliffe method). 728 plasma concentration samples across 50 patients (page 4)."
  )

  ini({
    # Structural parameters (Toffoli 2001 Table 2, 'Mean with covariables' column, page 4).
    # ka is NOT reported in Toffoli 2001 and is fixed here at 0.6/h, a literature-typical
    # value for oral etoposide (Slevin et al. 1989 [Toffoli ref 35] / Arbuck et al. 1986
    # [Toffoli ref 26], which the paper cites as the source of initial PK parameter values).
    # See vignette Assumptions and deviations for the full rationale.
    lka       <- fixed(log(0.6));   label("Absorption rate ka (1/h) - ASSUMED, not reported in Toffoli 2001")
    lcl       <- log(1.14);         label("Typical CL at reference CRCL=70 mL/min (L/h) - Toffoli 2001 Table 2 mean-with-covariables (page 4)")
    lvc       <- log(6.1);          label("Typical central volume of distribution (L) - Toffoli 2001 Table 2 mean-with-covariables (page 4)")
    lvp       <- log(12.174);       label("Peripheral volume of distribution (L) - derived Vp = (k12 * Vc) / k21 = (0.14 * 6.1) / 0.07 from Toffoli 2001 Table 2 microconstants (page 4)")
    lq        <- log(0.854);        label("Intercompartmental clearance Q (L/h) - derived Q = k12 * Vc = 0.14 * 6.1 from Toffoli 2001 Table 2 microconstants (page 4)")
    lfdepot   <- log(0.44);         label("Oral bioavailability F - Toffoli 2001 Table 2 (page 4)")

    # Additive-linear CL slope on CRCL (Toffoli 2001 Equation 2, page 4).
    # Equation 2 (least-squares regression, r=0.57, P<0.001):
    #   CL (L/h) = 0.74 + 0.0057 * CRCL (mL/min)
    # Recentred on median CRCL = 70 mL/min: CL = 1.14 + 0.0057 * (CRCL - 70).
    e_crcl_cl <- fixed(0.0057);     label("Additive CL slope on CRCL (L/h per mL/min) - Toffoli 2001 Equation 2 (page 4)")

    # IIV -- residual CV% after covariates from Toffoli 2001 Table 2 (page 4),
    # converted to log-scale variance via omega^2 = log(CV^2 + 1).
    # No IIV on ka (ka is fixed).
    # For Vp and Q the paper reports IIV on the microconstants k12 (13% CV) and
    # k21 (20% CV); these are mapped to the derived macroconstants Q and Vp
    # respectively as a pragmatic approximation (see vignette Errata).
    etalcl     ~ 0.02840  # 17% CV on CL (Toffoli 2001 Table 2)
    etalvc     ~ 0.02226  # 15% CV on Vc (Toffoli 2001 Table 2)
    etalq      ~ 0.01677  # 13% CV on k12 mapped to Q (Toffoli 2001 Table 2)
    etalvp     ~ 0.03922  # 20% CV on k21 mapped to Vp (Toffoli 2001 Table 2)
    etalfdepot ~ 0.04725  # 22% CV on F (Toffoli 2001 Table 2)

    # Residual error is not explicitly reported in Toffoli 2001. Using a proportional
    # error consistent with the assay intra-assay CV (5% at 0.4 mg/mL, 11% at 0.06
    # mg/mL; page 3, 'Drug assay'); prominently flagged in vignette Errata.
    propSd <- 0.10; label("Proportional residual error (fraction) - ASSUMED from assay validation; not reported in Toffoli 2001")
  })

  model({
    # Additive-linear typical CL on CRCL (Toffoli 2001 Equation 2, recentred on median).
    #   Eq. 2:  CL = 0.74 + 0.0057 * CRCL
    #   At median CRCL = 70:  CL_typical = 0.74 + 0.0057 * 70 = 1.139 ~ 1.14 (matches Table 2)
    #   Recentred:  CL_typical(CRCL) = 1.14 + 0.0057 * (CRCL - 70)
    cl_typical <- exp(lcl) + e_crcl_cl * (CRCL - 70)

    ka  <- exp(lka)
    cl  <- cl_typical * exp(etalcl)
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq  + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot + etalfdepot)

    # Total etoposide plasma concentration (mg/L); dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
