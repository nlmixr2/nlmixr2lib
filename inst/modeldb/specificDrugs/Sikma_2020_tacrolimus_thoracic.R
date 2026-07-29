Sikma_2020_tacrolimus_thoracic <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for oral whole-blood ",
    "tacrolimus in 30 adult thoracic organ transplant recipients (10 heart, ",
    "20 lung) during the first 6 postoperative days at the University Medical ",
    "Center Utrecht intensive care unit (Sikma 2020 EJDMP). Apparent ",
    "clearance CL/F, apparent volumes V1/F and V2/F, inter-compartmental ",
    "clearance Q/F, and first-order absorption rate ka are estimated; ",
    "bioavailability F is fixed at 1. Only the inter-individual variability ",
    "of CL/F was identifiable in the source dataset; all other IIV elements ",
    "were not estimated. Inter-occasion (dose-to-dose) variability dominated ",
    "the variance structure but is not encoded structurally in this ",
    "extraction. No covariates were retained in the final model."
  )
  reference <- paste0(
    "Sikma MA, Hunault CC, Van Maarseveen EM, Huitema ADR, Van de Graaf EA, ",
    "Kirkels JH, Verhaar MC, Grutters JC, Kesecioglu J, De Lange DW. ",
    "High Variability of Whole-Blood Tacrolimus Pharmacokinetics Early After ",
    "Thoracic Organ Transplantation. Eur J Drug Metab Pharmacokinet. ",
    "2020;45(1):123-134. doi:10.1007/s13318-019-00591-7."
  )
  vignette <- "Sikma_2020_tacrolimus_thoracic"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 30L,
    n_studies       = 1L,
    n_observations  = "1180 whole-blood tacrolimus concentrations across 119 twelve-hour profiles (median 5 profiles per patient, range 1-6); 46 of 1180 (3.9%) observations below the lower limit of quantification (0.5 ng/mL) were discarded.",
    age_range       = "34-60 years (Q1-Q3 of 30 patients)",
    age_median      = "43 years",
    weight_range    = "61-86 kg (Q1-Q3)",
    weight_median   = "73.5 kg",
    height_median   = "173.5 cm",
    sex_female_pct  = 50.0,
    race_ethnicity  = "Not reported in source paper (single-centre Utrecht, Netherlands cohort).",
    disease_state   = paste0(
      "Adult thoracic organ transplant recipients admitted to the intensive ",
      "care unit during the first 6 postoperative days. Heart transplant ",
      "indications (n=10): ischaemic cardiomyopathy 5, non-ischaemic ",
      "cardiomyopathy 5. Lung transplant indications (n=20): cystic fibrosis ",
      "10, chronic obstructive pulmonary disease 4, interstitial lung disease ",
      "6. Eighteen of 20 lung recipients underwent double-lung ",
      "transplantation. All patients met at least one Systemic Inflammatory ",
      "Response Syndrome criterion during the observation window; 93% (28 of ",
      "30) had at least one shock episode; gut dysmotility was observed in ",
      "97% (29 of 30). Postoperative ECMO was used in 27% (8 of 30) for a ",
      "median of 4 days (Q1-Q3 2-6)."
    ),
    dose_range      = paste0(
      "Oral tacrolimus (Prograf, Astellas Pharma Europe) twice daily, ",
      "starting at 0.1 mg/kg/dose for lung recipients and 2 mg/dose for ",
      "heart recipients on the day after transplantation. Dose adjustments ",
      "during the first 6 postoperative days were guided by the 12-hour ",
      "post-dose whole-blood trough concentration (target 9-15 ng/mL) and ",
      "the attending physician's assessment of drug-drug interactions, gut ",
      "dysmotility, and liver injury; steady state was not necessarily ",
      "reached at the time of dose adjustment."
    ),
    regions         = "Netherlands (University Medical Center Utrecht, June 2013 - March 2015).",
    co_medication   = paste0(
      "Triple immunosuppression with corticosteroids (prednisolone) plus ",
      "mycophenolate mofetil. Lung recipients additionally received ",
      "basiliximab induction on postoperative days 1 and 4. Drug-drug ",
      "interactions with potential CYP3A4/5 or ABCB1 inhibitors / inducers ",
      "occurred in 100% of patients; 0-6 interacting drugs per patient that ",
      "would increase tacrolimus exposure and 0-2 that would decrease ",
      "exposure. No covariate, including these co-medications, was retained ",
      "in the final model."
    ),
    notes           = paste0(
      "Single-centre prospective study (NTR 3912 / EudraCT 2012-001909-24). ",
      "Whole-blood samples drawn pre-dose and at 1, 1.5, 2, 2.5, 3, 4, 6, 8, ",
      "and 12 hours after a daily index dose during ICU admission. Tacrolimus ",
      "quantified by HPLC-MS/MS (lower limit of quantification 0.5 ng/mL, ",
      "linear range 1-50 ng/mL; intra-day imprecision < 5%). NCA summary ",
      "across all 119 profiles (Table 1): median Cmax 18.5 ng/mL (range ",
      "2.1-74.7), median Tmax 1.6 h (range 0.4-8.0), median AUC 151.2 ng*h/mL ",
      "(range 31.2-2525), median terminal half-life 9.4 h (range 6.0-31.4)."
    )
  )

  ini({
    # Final population parameter estimates from Sikma 2020 Table 2 (NONMEM
    # 7.3.0; SAEM estimation; importance-sampling likelihood; SIR-based 95%
    # confidence intervals). The published model parameterises apparent
    # clearances and apparent volumes (CL/F, V1/F, Q/F, V2/F); the
    # absolute bioavailability F is not identifiable in oral-only data and
    # was fixed at 1.

    # ---- Structural disposition (Sikma 2020 Table 2) -----------------------
    lcl <- log(19.6) ; label("Apparent clearance CL/F (L/h)")                          # Sikma 2020 Table 2 CL/F = 19.6 (95% CI 16.2-22.9)
    lvc <- log(231)  ; label("Apparent central volume of distribution V1/F (L)")        # Sikma 2020 Table 2 V1/F = 231 (95% CI 199-267)
    lq  <- log(58.2) ; label("Apparent inter-compartmental clearance Q/F (L/h)")        # Sikma 2020 Table 2 Q/F  = 58.2 (95% CI 49.7-69.3)
    lvp <- log(521)  ; label("Apparent peripheral volume of distribution V2/F (L)")     # Sikma 2020 Table 2 V2/F = 521 (95% CI 441-634)
    lka <- log(0.579); label("First-order absorption rate constant ka (1/h)")           # Sikma 2020 Table 2 ka   = 0.579 (95% CI 0.456-0.778)

    # Bioavailability fixed at 1 (Sikma 2020 Section 2.8: "Although the
    # absolute bioavailability was not identifiable, the variability in the
    # relative bioavailability was estimated similarly with theta_k fixed at
    # 1").
    lfdepot <- fixed(log(1)) ; label("Oral bioavailability F (FIXED at 1)")             # Sikma 2020 Table 2 F = Fixed at 1

    # ---- Inter-individual variability (Sikma 2020 Table 2) ------------------
    # Reported as %CV; converted to log-scale variance via
    #   omega^2 = log(1 + CV^2).
    # Only the IIV on CL/F was estimated; "the estimation of all other IIV
    # elements yielded estimates close to 0 and/or unsuccessful runs, most
    # likely due to the fact that variability was dominated by the
    # corresponding IOV" (Sikma 2020 Section 3.4). Following the source
    # paper's final model, only the etalcl term is included here.
    #   CL/F  CV 34.6% -> log(1 + 0.346^2) = 0.11307
    etalcl ~ 0.11307                                                                    # Sikma 2020 Table 2 IIV CL/F = 34.6% (95% CI 24.2-48.6)

    # Inter-occasion (dose-to-dose) variability was reported in Sikma 2020
    # Table 2 and dominates the variance structure (CL/F 29.5%, V1/F 35.1%,
    # ka 98.3%, F 55.0%; Q/F and V2/F not estimated). nlmixr2lib popPK
    # extractions do not standardise IOV encoding (it requires an OCC
    # variable in the dataset), so the source paper's IOV terms are
    # documented in the vignette's Assumptions and deviations section rather
    # than coded structurally in ini(). Downstream users can add a per-
    # occasion eta in rxode2 if required.

    # ---- Residual error (Sikma 2020 Table 2) --------------------------------
    # Proportional residual error 14.0% (Sikma 2020 Section 2.8 Eq. 3:
    # residual error proportional to the predicted concentration).
    propSd <- 0.140 ; label("Proportional residual error (fraction)")                   # Sikma 2020 Table 2 RUV = 14.0% (95% CI 13.3-14.6)
  })

  model({
    # ---- Individual structural parameters ----------------------------------
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    ka     <- exp(lka)
    fdepot <- exp(lfdepot)

    # Micro-constants for the two-compartment disposition (Sikma 2020 Fig. 1
    # caption: k23 = Q/V1, k32 = Q/V2; here k12 = q/vc and k21 = q/vp using
    # nlmixr2 compartment naming).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: two-compartment first-order oral absorption -----------
    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot   - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability fixed at 1 (Sikma 2020 Table 2).
    f(depot) <- fdepot

    # ---- Observation and error ---------------------------------------------
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000 to
    # report whole-blood tacrolimus concentration in ng/mL (Sikma 2020
    # reported tacrolimus concentrations in ng/mL throughout).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
