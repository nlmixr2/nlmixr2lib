Solana_2014_omeprazole <- function() {
  description <- "Two-compartment intravenous-infusion population PK model for omeprazole in 40 critically ill children (Solana 2014), with fixed Anderson-Holford allometric body-weight scaling on all four disposition parameters (exponents 0.75 on CL and Q, 1.00 on Vc and Vp; reference 70 kg). Between-patient variability was retained on CL only; residual error is proportional."
  reference   <- "Solana MJ, Colom H, Lopez-Herce J, Urbano J, Gonzalez R, Lopez J, Manzanares C, Carrillo A. Population pharmacokinetics of omeprazole in critically ill pediatric patients. Ther Drug Monit. 2014;36(4):519-527. doi:10.1097/FTD.0000000000000033"
  vignette    <- "Solana_2014_omeprazole"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Fixed Anderson-Holford allometric scaling on all four disposition parameters with reference 70 kg: exponent 0.75 on CL and Q, 1.00 on Vc and Vp (Solana 2014 abstract Results paragraph: 'Allometric size models seemed to predict changes adequately in all the pharmacokinetic parameters'; the source abstract reports the typical values per 70 kg and does not estimate the allometric exponents, so the exponents are taken as the theoretical Anderson-Holford values held fixed).",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 40L,
    n_studies       = 1L,
    n_observations  = 186L,
    age_range       = "Critically ill paediatric patients (full numeric age / weight ranges not reported in the abstract).",
    weight_range    = "Critically ill paediatric patients (per-patient weight range not reported in the abstract; allometric scaling implies the cohort spans the typical paediatric ICU body-weight range relative to the 70 kg adult reference).",
    sex_female_pct  = NA_real_,
    disease_state   = "Critically ill children (paediatric intensive care).",
    dose_range      = "Intravenous omeprazole 0.5 or 1 mg/kg twice daily, randomized between the two dose levels.",
    sampling        = "Plasma samples drawn at 0.5, 2, 6, 12, 24, and 48 hours after the first infusion.",
    regions         = "Spain (Hospital General Universitario Gregorio Maranon, Madrid).",
    notes           = "ABSTRACT-ONLY EXTRACTION: the full-text Therapeutic Drug Monitoring article (DOI 10.1097/FTD.0000000000000033) was not available on disk for this extraction; all structural and parameter information was taken from the PubMed abstract (PMID 24365987). The numeric covariate ranges (per-patient age, weight, sex distribution) are not enumerated in the abstract and are recorded here as NA / narrative; they are not load-bearing for simulation because the only covariate used by the model is body weight via allometric scaling. The infusion duration is not stated in the abstract and is assumed to be a short IV bolus-like infusion in the vignette simulation; users with the full text can adjust the infusion rate accordingly."
  )

  ini({
    # Final-model typical values from the Solana 2014 abstract (PMID 24365987,
    # Results paragraph). The abstract reports per-70-kg typical values with the
    # %RSE in parentheses; the values below mirror that table exactly.
    lcl <- log(24.9);  label("Plasma clearance CL at reference 70 kg (L/h)")              # Solana 2014 abstract Results: CL = 24.9 L/h per 70 kg (RSE 10.08%)
    lq  <- log(53.9);  label("Distributional clearance Q at reference 70 kg (L/h)")       # Solana 2014 abstract Results: Q  = 53.9 L/h per 70 kg (RSE 11.00%)
    lvc <- log(4.23);  label("Central volume of distribution Vc at reference 70 kg (L)")  # Solana 2014 abstract Results: Vc = 4.23 L per 70 kg (RSE 19.62%)
    lvp <- log(674);   label("Peripheral volume of distribution Vp at reference 70 kg (L)") # Solana 2014 abstract Results: Vp = 674 L per 70 kg (RSE 0.89%); see vignette Assumptions and deviations for the source-faithful note about the unusually large Vp:Vc ratio.

    # Fixed Anderson-Holford allometric exponents (reference 70 kg). The
    # abstract reports that allometric size models adequately predicted changes
    # in all PK parameters and does not estimate the exponents; they are taken
    # as the theoretical Anderson-Holford values held fixed (0.75 on
    # clearances, 1.00 on volumes) as is standard practice in paediatric popPK.
    e_wt_cl <- fixed(0.75);  label("Allometric exponent of WT on CL (unitless)")   # Solana 2014 abstract Results paragraph; standard Anderson-Holford clearance exponent held fixed
    e_wt_q  <- fixed(0.75);  label("Allometric exponent of WT on Q (unitless)")    # Solana 2014 abstract Results paragraph; standard Anderson-Holford clearance exponent held fixed
    e_wt_vc <- fixed(1.00);  label("Allometric exponent of WT on Vc (unitless)")   # Solana 2014 abstract Results paragraph; standard Anderson-Holford volume exponent held fixed
    e_wt_vp <- fixed(1.00);  label("Allometric exponent of WT on Vp (unitless)")   # Solana 2014 abstract Results paragraph; standard Anderson-Holford volume exponent held fixed

    # Between-patient variability (IIV) -- Solana 2014 abstract Results paragraph:
    # "Between-patient variability could only be associated with plasma clearance
    # (CL). ... High values of between-patient variability of CL [75.50% (2.60%)]
    # ... were still found in the final model."
    # The reported 75.5% is the CV% of the log-normal CL distribution; conversion
    # to the log-scale variance for nlmixr2's `eta ~ var` syntax follows the
    # standard log-normal identity omega^2 = log(1 + CV^2):
    #   omega^2 = log(1 + 0.755^2) = log(1 + 0.570025) = log(1.570025) = 0.45106.
    # No IIV is reported on Q, Vc, or Vp in the abstract.
    etalcl ~ 0.45106  # Solana 2014 abstract: IIV(CL) = 75.5% CV -> log(1 + 0.755^2) = 0.45106

    # Residual variability -- Solana 2014 abstract Results paragraph:
    # "residual variability [130.0% (5.26%)] were still found in the final
    # model." The single reported residual variability is treated as
    # proportional (the standard interpretation when a high %CV residual is
    # reported without an additive component). The proportional SD is the
    # fractional CV directly, so propSd = 1.30.
    propSd <- 1.30;  label("Proportional residual error (fraction)")  # Solana 2014 abstract Results: residual variability = 130.0% (RSE 5.26%)
  })

  model({
    # Individual disposition parameters with fixed Anderson-Holford allometric
    # body-weight scaling on all four disposition parameters (reference 70 kg).
    cl <- exp(lcl + etalcl) * (WT / 70) ^ e_wt_cl
    q  <- exp(lq)           * (WT / 70) ^ e_wt_q
    vc <- exp(lvc)          * (WT / 70) ^ e_wt_vc
    vp <- exp(lvp)          * (WT / 70) ^ e_wt_vp

    # Micro-constants for the linear two-compartment ODE
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Linear two-compartment intravenous-infusion model. The dose is
    # administered directly into the central compartment; the event table
    # supplies the infusion rate (or duration) on the dosing rows. There is
    # no first-order absorption depot in this model.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration in the central compartment; dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
