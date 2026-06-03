Boer_2015_cisplatin <- function() {
  description <- "Two-compartment population PK model for long-term circulating platinum (Pt) decay after cisplatin-based chemotherapy in adult testicular cancer survivors followed 1-13 years post-treatment (Boer 2015). Dose is the cumulative cisplatin dose expressed as elemental Pt in mg (multiply cumulative cisplatin in mg by 0.6502, the Pt/cisplatin mass ratio 195.08/300.05). An apparent bioavailability F1 (fdepot) accounts for the fraction of the administered Pt remaining in the body after the rapid pre-measurement urinary-excretion phase. Pt is assumed to be cleared solely via urine."
  reference <- "Boer H, Proost JH, Nuver J, Bunskoek S, Gietema JQ, Geubels BM, Altena R, Zwart N, Oosting SF, Vonk JM, Lefrandt JD, Uges DRA, Meijer C, de Vries EGE, Gietema JA. Long-term exposure to circulating platinum is associated with late effects of treatment in testicular cancer survivors. Ann Oncol. 2015;26(11):2305-2310. doi:10.1093/annonc/mdv369"
  vignette <- "Boer_2015_cisplatin"
  units <- list(time = "day", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at start of chemotherapy",
      units       = "years",
      type        = "continuous",
      notes       = "Tested in covariate screening (Boer 2015 Results); did not significantly improve OFV and was not retained in the final population PK model. Univariate correlation between age and Pt AUC(1-3 yr) disappeared after adjustment for renal function."
    ),
    WT = list(
      description = "Body weight at start of chemotherapy",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested in covariate screening (Boer 2015 Results); did not significantly improve OFV and was not retained in the final population PK model. No correlation with Pt AUC(1-3 yr)."
    ),
    HT = list(
      description = "Body height at start of chemotherapy",
      units       = "m",
      type        = "continuous",
      notes       = "Tested in covariate screening (Boer 2015 Results); did not significantly improve OFV and was not retained in the final population PK model."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 99L,
    n_studies       = 1L,
    n_observations  = 240L,
    age_range       = "17-53 years (at start of chemotherapy)",
    age_median      = "29 years (at start of chemotherapy); 39 years (range 23-64) at follow-up assessment",
    weight_range    = "not reported as range; demographic Table 1 covers cohort baseline",
    sex_female_pct  = 0,
    race_ethnicity  = "not reported; single-centre Dutch cohort (University Medical Centre Groningen)",
    disease_state   = "Non-seminomatous testicular cancer survivors treated with cisplatin-based chemotherapy (BEP, EP, 3xBEP/1xEP, or other cisplatin-based regimen) between 1988 and 2000. Stage II 53%, stage III 5%, stage IV 42%; IGCCCG risk 'good' 56%, 'intermediate' 34%, 'poor' 9%.",
    dose_range      = "Cumulative cisplatin 554-1713 mg (median 809 mg); 275-800 mg/m^2 (median 400 mg/m^2); standard testicular-cancer cisplatin regimens delivered over 9-12 weeks. The model dose is the elemental-Pt equivalent (multiply cumulative cisplatin in mg by 195.08/300.05 = 0.6502).",
    regions         = "Netherlands (single-centre, University Medical Centre Groningen)",
    notes           = "240 serum Pt measurements from 98 patients (one to three samples per patient) collected 0.9-13.2 years post-start-of-chemotherapy (median 5.0 years); one 24-h urine sample from 91 patients (median 6.6 years post-chemotherapy). Serum and urinary excretion rate were fit simultaneously. Pt assayed by adsorptive voltammetry after high-pressure decomposition; LLOQ 6 pg/g serum; intra/inter-assay CVs 6%/5%. Cohort follow-up assessment at median 9 years (range 3-15) post-chemotherapy. Baseline demographics per Table 1 of Boer 2015."
  )

  ini({
    # Structural parameters from Boer 2015 Supplementary Table S1 (final
    # population PK model). The dose-input column is amount of elemental Pt
    # in mg (computed in the source dataset as cumulative cisplatin mg/m^2
    # times body surface area, then times the Pt/cisplatin mass ratio
    # 195.08/300.05 = 0.6502).
    lcl     <- log(0.0220);    label("Clearance (L/day)")                      # Boer 2015 Suppl Table S1: CL = 0.0220 L/day (95% CI 0.0197-0.0246)
    lvc     <- log(6.61);      label("Central volume of distribution V1 (L)")  # Boer 2015 Suppl Table S1: V1 = 6.61 L (95% CI 5.20-8.03)
    lvp     <- log(8.05);      label("Peripheral volume of distribution V2 (L)") # Boer 2015 Suppl Table S1: V2 = 8.05 L (95% CI 6.55-10.17)
    lq      <- log(0.00531);   label("Intercompartmental clearance Q (L/day)") # Boer 2015 Suppl Table S1: Q = 0.00531 L/day (95% CI 0.00461-0.00611)
    lfdepot <- log(0.000158);  label("Apparent bioavailability F1 (fraction of administered Pt remaining after the pre-measurement phase, unitless)") # Boer 2015 Suppl Table S1: F1 = 0.000158 (95% CI 0.000137-0.000188)

    # Inter-individual variability (omega^2 = log(CV^2 + 1) for lognormal IIV)
    # Boer 2015 Suppl Table S1 reports IIV (%) only for CL (13%), V2 (22%),
    # and F1 (17%). IIV on V1 and Q was tested but did not improve OFV and
    # was not retained (Boer 2015 Results, "population PK model" section).
    etalcl     ~ 0.016759   # Boer 2015 Suppl Table S1: IIV CL = 13% -> log(1 + 0.13^2) = 0.016759
    etalvp     ~ 0.047284   # Boer 2015 Suppl Table S1: IIV V2 = 22% -> log(1 + 0.22^2) = 0.047284
    etalfdepot ~ 0.028471   # Boer 2015 Suppl Table S1: IIV F1 = 17% -> log(1 + 0.17^2) = 0.028471

    # Residual error. Boer 2015 Results, "population PK model" section:
    # "Proportional residual error in the final model was 34%". Suppl Table
    # S1 does not list a numeric SD; the 34% in the main text is the CV.
    propSd <- 0.34; label("Proportional residual error (fraction)")  # Boer 2015 Results: 34% proportional residual error
  })

  model({
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc)
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Apparent bioavailability applied to IV-bolus dosing into central.
    # Dose enters via event-table rows with cmt = central. The model is the
    # "post-treatment residual" representation of Boer 2015 (the cisplatin
    # treatment period of 9-12 weeks is treated as negligible relative to
    # the 1-13 year follow-up); F1 captures the fraction of the cumulative
    # Pt dose surviving the rapid urinary-excretion phase that precedes the
    # first serum measurement.
    f(central) <- fdepot

    # Concentration scaling: dose in mg, vc in L -> central/vc has units
    # mg/L. Multiply by 1000 to express in ug/L, matching the y-axis of
    # Boer 2015 Figure 1.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
