Mondick_2006_dactinomycin <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for actinomycin-D",
    "(dactinomycin) in 33 pediatric and young-adult patients (1.58-20.3",
    "years) with Wilms' tumor or rhabdomyosarcoma. All disposition",
    "parameters are allometrically scaled by total body weight,",
    "normalized to a reference weight of 70 kg, with theory-based fixed",
    "exponents (0.75 on clearances, 1.0 on volumes; not explicitly stated",
    "in the abstract). Inter-individual variability was reported only for",
    "V1 (54.4% CV) and CL (57.2% CV); residual error and the remaining",
    "IIV terms (V2, V3, Q2, Q3) were not reported in the source",
    "conference abstract and are encoded as fixed(0). Mondick 2006",
    "(PAGE 15 Abstr 938)."
  )
  reference <- paste(
    "Mondick JT, Gibiansky L, Gastonguay MR, Veal GJ, Barrett JS.",
    "Acknowledging Parameter Uncertainty in the Simulation-Based Design",
    "of an Actinomycin-D Pharmacokinetic Study in Pediatric Patients",
    "with Wilms' Tumor or Rhabdomyosarcoma.",
    "PAGE 15 (2006) Abstr 938.",
    "https://www.page-meeting.org/?abstract=938"
  )
  vignette <- "Mondick_2006_dactinomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on all disposition parameters with reference",
        "weight 70 kg (Mondick 2006 Results). Clearances (CL, Q2, Q3)",
        "scale as (WT/70)^0.75 and volumes (V1, V2, V3) as (WT/70)^1.0;",
        "the abstract states 'all parameters allometrically scaled by",
        "total body weight, normalized to a weight of 70 kg' but does",
        "not print the exponents, so the theory-based Holford/Anderson",
        "values are assumed and held fixed. Cohort age range 1.58-20.3",
        "years; extrapolation to infants under one year is the principal",
        "scenario the source paper investigates by simulation."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human (pediatric and young adult)",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "1.58-20.3 years",
    weight_range   = "not reported in abstract",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Pediatric / young-adult patients with Wilms' tumor or",
      "rhabdomyosarcoma receiving actinomycin-D (AMD) as part of",
      "standard chemotherapy. Age, gender, and body size were screened",
      "as potential covariates; only body-weight allometric scaling was",
      "retained in the final model. Infants under 12 months are",
      "underrepresented (the cohort starts at 1.58 years); the abstract",
      "motivates a prospective study to fill that gap."
    ),
    dose_range     = paste(
      "Not reported in the abstract. Sampling schemes designed for the",
      "prospective Children's Oncology Group Phase I Consortium trial",
      "are (1) 5 min, 10 min, 2-3 h, 24-28 h, 48-96 h post-dose and (2)",
      "5 min, 0.75-1.5 h, 5-6 h, 24-28 h, 48-96 h post-dose, consistent",
      "with the typical short IV bolus / brief IV infusion of",
      "actinomycin-D in pediatric oncology."
    ),
    regions        = paste(
      "Children's Hospital of Philadelphia (USA);",
      "Northern Institute for Cancer Research, Newcastle upon Tyne (UK)."
    ),
    notes          = paste(
      "Source is a one-page PAGE 2006 conference abstract; no",
      "supplement or subsequent full-text publication of this PK model",
      "is available on disk. NIH Award #CA098543-0251.",
      "Final NONMEM fit reports %CV next to V1 and CL only; the",
      "remaining IIVs and the residual error are not reported and are",
      "encoded here as fixed(0)."
    )
  )

  ini({
    # Structural disposition parameters at reference weight 70 kg.
    # All values from the Results section of Mondick 2006 (PAGE 15 Abstr 938):
    # "Population mean (%CV) parameter estimates for V1, V2, V3, CL, Q2, and Q3
    # were 3.88 (54.4%) L, 443 L, 23.7 L, 8.27 (57.2%) L/h, 252 L/h, and 12.7 L/h."
    lvc  <- log(3.88);  label("Central volume V1 at reference 70 kg (L)")               # Mondick 2006 Results: V1 = 3.88 L
    lvp  <- log(443);   label("Peripheral volume V2 at reference 70 kg (L)")            # Mondick 2006 Results: V2 = 443 L
    lvp2 <- log(23.7);  label("Second peripheral volume V3 at reference 70 kg (L)")     # Mondick 2006 Results: V3 = 23.7 L
    lcl  <- log(8.27);  label("Clearance CL at reference 70 kg (L/h)")                  # Mondick 2006 Results: CL = 8.27 L/h
    lq   <- log(252);   label("Intercompartmental clearance Q2 at reference 70 kg (L/h)")  # Mondick 2006 Results: Q2 = 252 L/h
    lq2  <- log(12.7);  label("Intercompartmental clearance Q3 at reference 70 kg (L/h)")  # Mondick 2006 Results: Q3 = 12.7 L/h

    # Allometric exponents (theory-based, Holford/Anderson):
    # 0.75 on clearance terms, 1.0 on volume terms. The Mondick 2006 abstract
    # states only "all parameters allometrically scaled by total body weight,
    # normalized to a weight of 70 kg" and does not print the exponents; the
    # canonical theory-based values are assumed and held fixed. The assumption
    # is documented in the vignette Assumptions and deviations section.
    e_wt_cl  <- fixed(0.75);  label("Allometric exponent on clearance (unitless, fixed)")   # standard theory-based; not printed in Mondick 2006 abstract
    e_wt_vc  <- fixed(1.0);   label("Allometric exponent on volume (unitless, fixed)")      # standard theory-based; not printed in Mondick 2006 abstract

    # IIV (log-normal) for V1 and CL. omega^2 = log(CV^2 + 1).
    #   V1: CV = 54.4%; omega^2 = log(1 + 0.544^2) = 0.25923
    #   CL: CV = 57.2%; omega^2 = log(1 + 0.572^2) = 0.28306
    # The abstract reports %CV alongside V1 and CL only; no eta correlation is
    # printed. The remaining disposition parameters (V2, V3, Q2, Q3) have no
    # %CV in the abstract and are encoded with IIV fixed at 0.
    etalvc ~ 0.25923                                                                       # Mondick 2006 Results: V1 %CV = 54.4
    etalcl ~ 0.28306                                                                       # Mondick 2006 Results: CL %CV = 57.2
    etalvp  ~ fixed(0)                                                                     # not reported in Mondick 2006 abstract
    etalvp2 ~ fixed(0)                                                                     # not reported in Mondick 2006 abstract
    etalq   ~ fixed(0)                                                                     # not reported in Mondick 2006 abstract
    etalq2  ~ fixed(0)                                                                     # not reported in Mondick 2006 abstract

    # Residual error is not reported in the Mondick 2006 abstract. Encoded as
    # fixed(0) so the typical-value structural simulation is reproducible;
    # the assumption is documented in the vignette Assumptions and deviations
    # section.
    propSd <- fixed(0); label("Proportional residual error (fraction, fixed)")             # not reported in Mondick 2006 abstract
  })

  model({
    # Individual disposition parameters at the subject's body weight, with the
    # theory-based allometric scaling described in the Mondick 2006 abstract.
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl
    q   <- exp(lq   + etalq)   * (WT / 70)^e_wt_cl
    q2  <- exp(lq2  + etalq2)  * (WT / 70)^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc
    vp  <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vc
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central)     <- -(kel + k12 + k13) * central +
                          k21 * peripheral1 +
                          k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Plasma concentration in ng/mL (dose in mg / volume in L gives mg/L = ug/mL;
    # multiply by 1000 to convert to ng/mL, the conventional reporting unit for
    # actinomycin-D plasma concentrations).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
