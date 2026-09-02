Cortez_2015_nevirapine_rat_la1 <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model with first-order absorption and",
    "linear elimination for the NVP LA-1 long-acting injectable nevirapine",
    "microparticle suspension (nevirapine cores, 411 um median diameter,",
    "coated with poly(vinyl alcohol) and poly(D,L-lactide) to a 13 um",
    "membrane) after a single subcutaneous dose to Sprague-Dawley rats",
    "(study 1, JNN00006). Fitted by naive pooling in Phoenix NLME; no",
    "between-subject variability was estimated because only n = 4 animals",
    "contributed. Vc, Vp, CL and CLp are apparent (/F) formulation-dependent",
    "quantities distorted by flip-flop kinetics, not true systemic",
    "disposition parameters -- the authors state this explicitly because no",
    "intravenous rat reference arm was available. Body weight enters de",
    "facto through fixed allometric exponents (-0.25 on ka, 1 on volumes,",
    "0.75 on clearances) referenced to the 0.33 kg mean study weight, which",
    "is what makes the companion human projection",
    "(Cortez_2015_nevirapine_human_la1) possible.",
    sep = " ")
  reference <- paste(
    "Cortez JM Jr, Quintero R, Moss JA, Beliveau M, Smith TJ, Baum MM.",
    "Pharmacokinetics of injectable, long-acting nevirapine for HIV",
    "prophylaxis in breastfeeding infants.",
    "Antimicrob Agents Chemother. 2015;59(1):59-66.",
    "doi:10.1128/AAC.03906-14.",
    sep = " ")
  vignette <- "Cortez_2015_nevirapine"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot       = list(analyte = "nevirapine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "nevirapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nevirapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size covariate, reference 0.33 kg (the mean weight of",
        "the n = 4 study-1 rats, Cortez 2015 Methods 'Animals' and Table 2",
        "column header). Cortez 2015 Methods 'Structural PK model buildup'",
        "states the covariate model is theta_i = theta_typical *",
        "(Cov_i / Cov_reference)^theta_eff with nominal exponents -0.25 on",
        "absorption rate, 1 on volumes and 0.75 on clearances; the exponents",
        "were imposed, not estimated, so they are fixed() here."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 4L,
    n_studies      = 1L,
    age_range      = "Adult (age not reported)",
    weight_range   = "Mean 0.33 kg (individual weights not reported)",
    weight_median  = "0.33 kg (reported as the mean)",
    disease_state  = "Healthy",
    dose_range     = paste(
      "Single 60 mg subcutaneous injection of NVP LA-1 suspended in 1 mL",
      "saline through a 14-gauge needle. Cortez 2015 Methods 'Population PK",
      "data set construction' assumes 75% of the nominal dose was actually",
      "delivered because of injection inefficiencies, so the modelled dose",
      "is 45 mg (135 mg/kg)."
    ),
    regions        = "United States (Charles River Laboratories, Shrewsbury, MA)",
    notes          = paste(
      "Study 1 (JNN00006, March 2006). Serial jugular-vein sampling on days",
      "1 (1 and 6 h post-dose), 2, 3, 4, 6, 10, 14 and 28. Plasma nevirapine",
      "by LC-MS/MS (API 5000, APCI+, m/z 267.1 -> 226.1; glyburide internal",
      "standard) with an LLOQ of 1 ng/mL; concentrations below the LLOQ were",
      "set to missing rather than imputed. Subcutaneous bioavailability was",
      "assumed to be 100% (Cortez 2015 Discussion, citing reference 38).",
      "See Cortez 2015 Methods ('Animals', 'NVP LA administration and",
      "sampling', 'Analysis of NVP in rat plasma', 'Population PK data set",
      "construction', 'Population PK model development'), Table 2 (study 1",
      "column) and Figure 3A."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Cortez 2015 Table 2, "Study 1, NVP LA-1
    # (n = 4, mean wt 0.33 kg)" column. Values are reported in litres and
    # litres/day for a 0.33 kg rat and are apparent (/F) quantities; the
    # paper's Discussion is explicit that they do not correspond to true
    # systemic parameters. The paper's CL_p is the inter-compartmental
    # ("peripheral") clearance and maps onto the canonical lq.
    # ------------------------------------------------------------------
    lka <- log(1.89)  ; label("Absorption rate constant ka (1/day) for a 0.33 kg rat")                       # Cortez 2015 Table 2: Ka = 1.89 /day (RSE 26.1%); consistent with the tabulated t1/2,a = 0.367 day
    lvc <- log(8.66)  ; label("Apparent central volume Vc/F (L) for a 0.33 kg rat")                          # Cortez 2015 Table 2: Vc/F = 8.66 L (RSE 17.9%)
    lvp <- log(14.3)  ; label("Apparent peripheral volume Vp/F (L) for a 0.33 kg rat")                       # Cortez 2015 Table 2: Vp/F = 14.3 L (RSE not available)
    lcl <- log(3.43)  ; label("Apparent clearance CL/F (L/day) for a 0.33 kg rat")                           # Cortez 2015 Table 2: CL/F = 3.43 L/day (RSE 11.3%)
    lq  <- log(0.541) ; label("Apparent inter-compartmental clearance CLp/F (L/day) for a 0.33 kg rat")      # Cortez 2015 Table 2: CLp/F = 0.541 L/day (RSE 35.8%)

    # Subcutaneous bioavailability was ASSUMED to be 100%, not estimated.
    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (fraction)")                     # Cortez 2015 Discussion 'Pharmacokinetics of sustained-release NVP in rats': "The subcutaneous bioavailability, F, was assumed to be 100% for both formulations based on previous studies in rats (38)."

    # ------------------------------------------------------------------
    # Allometric exponents. Cortez 2015 Methods 'Structural PK model
    # buildup': "The nominal covariate effect was -0.25 on absorption
    # rate, 1 on volumes, and 0.75 on clearance (19)." These are nominal
    # (imposed) values with no reported uncertainty, hence fixed().
    # ------------------------------------------------------------------
    e_wt_ka <- fixed(-0.25) ; label("Allometric exponent on ka (unitless)")                           # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_vc <- fixed(1)     ; label("Allometric exponent on Vc (unitless)")                           # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_vp <- fixed(1)     ; label("Allometric exponent on Vp (unitless)")                           # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_cl <- fixed(0.75)  ; label("Allometric exponent on CL (unitless)")                           # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_q  <- fixed(0.75)  ; label("Allometric exponent on CLp (unitless)")                          # Cortez 2015 Methods 'Structural PK model buildup'

    # ------------------------------------------------------------------
    # No IIV. Cortez 2015 Methods 'Population PK model development':
    # "Due to the low number of rats (n = 4 or n = 10) available for each
    # formulation, no between-subject variability (BSV) was added in the
    # PK model". The naive pooled approach was used. No eta terms are
    # therefore declared -- see the vignette Assumptions and deviations.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Residual error -- Cortez 2015 equation 2,
    #   y_ij = yhat_ij * (1 + eps1_ij) + eps2_ij,
    # i.e. combined proportional + additive, which is nlmixr2's
    # add() + prop(). Table 2 (study 1) reports +/-28.6 (RSE 40.4%)
    # ug/liter additive and +/-50.1 (RSE 32.4%) % proportional.
    # 28.6 ug/L = 0.0286 mg/L = 0.0286 ug/mL, the model's concentration
    # unit (dose in mg / volume in L).
    # ------------------------------------------------------------------
    addSd  <- 0.0286 ; label("Additive residual error SD (ug/mL)")           # Cortez 2015 Table 2 study 1 'Error' row: +/-28.6 (40.4) ug/liter
    propSd <- 0.501  ; label("Proportional residual error SD (fraction)")    # Cortez 2015 Table 2 study 1 'Error' row: +/-50.1 (32.4) %
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. Allometric scaling on body weight with
    # the nominal exponents, referenced to the 0.33 kg study mean.
    # ------------------------------------------------------------------
    ka <- exp(lka) * (WT / 0.33)^e_wt_ka
    vc <- exp(lvc) * (WT / 0.33)^e_wt_vc
    vp <- exp(lvp) * (WT / 0.33)^e_wt_vp
    cl <- exp(lcl) * (WT / 0.33)^e_wt_cl
    q  <- exp(lq)  * (WT / 0.33)^e_wt_q

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 3. ODE system -- depot (subcutaneous injection site) feeding a
    # two-compartment linear disposition model.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 4. Bioavailability (assumed 100%).
    # ------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 5. Observation and error.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
