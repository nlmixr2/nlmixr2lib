Azechi_2024_tamibarotene_pediatric <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral tamibarotene (synthetic retinoid",
    "RAR-alpha/beta agonist) in pediatric and young-adult patients (4-23",
    "years) with recurrent or refractory solid tumors (Azechi 2024).",
    "Apparent oral clearance CL/F, apparent central volume V1/F, and apparent",
    "peripheral volume V2/F scale linearly with body surface area (BSA)",
    "referenced to the cohort mean of 0.995 m^2 (Table 1); inter-compartmental",
    "clearance Q/F, the absorption rate constant ka, and the absorption lag",
    "time tlag have no covariate effects. tlag was held fixed at 0.95 h in",
    "the published final model (the authors judged the post-covariate Tlag",
    "estimate of ~1.8 h to have low physiological validity and fell back to",
    "the pre-covariate value). Residual error is proportional with a 42.4%",
    "magnitude. The Methods section specifies an exponential IIV model on",
    "all five PK parameters but the paper reports no per-parameter omega",
    "magnitudes and no supplement exists; per operator decision (sidecar",
    "request-001 q1 = A, 2026-06-21) the five eta terms are encoded as",
    "fixed(0) so the published structural IIV declaration is preserved while",
    "remaining faithful to the absence of reported variance values. See the",
    "vignette Assumptions and deviations section for the resulting",
    "limitations on VPC-style validation.",
    sep = " "
  )
  reference <- paste(
    "Azechi T, Fukaya Y, Nitani C, Hara J, Kawamoto H, Taguchi T,",
    "Yoshimura K, Sato A, Hattori N, Ushijima T, Kimura T.",
    "Population Pharmacokinetics of Tamibarotene in Pediatric and Young",
    "Adult Patients with Recurrent or Refractory Solid Tumors.",
    "Curr Oncol. 2024;31(11):7155-7164.",
    "doi:10.3390/curroncol31110527.",
    "PMID 39590158.",
    sep = " "
  )
  vignette <- "Azechi_2024_tamibarotene_pediatric"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline). Azechi 2024 Table 1 reports",
        "cohort BSA mean = 0.995 m^2 (SD 0.389; median 0.805; range",
        "0.63-1.97). BSA computation formula (DuBois / Mosteller /",
        "Haycock) is not stated in the paper; recorded as unspecified.",
        "Enters the final model as linear scaling on CL/F, V1/F, and",
        "V2/F via the form parameter = typical_value * (BSA / 0.995)",
        "per Azechi 2024 Table 4 Final Model column",
        "(CL/F = tvCL/F * BSA/mean, V1/F = tvV1/F * BSA/mean,",
        "V2/F = tvV2/F * BSA/mean). The 'mean' centering value is the",
        "cohort BSA mean 0.995 m^2 from Table 1. Q/F has no BSA",
        "covariate. The Methods Section 2.6 calls the screened form an",
        "'allometric function' but the final-model column in Table 4",
        "shows linear BSA/mean scaling without an explicit power",
        "exponent; this model file encodes the exponents as fixed(1) to",
        "make the linear-scaling provenance explicit in ini().",
        sep = " "
      ),
      source_name        = "BSA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "4-23 years",
    age_median     = "8 years",
    weight_range   = "14.5-76.6 kg",
    weight_median  = "19.9 kg",
    sex_female_pct = 32.8,
    disease_state  = paste(
      "Pediatric and young adult patients with histologically confirmed",
      "advanced or recurrent solid tumors (sarcomas, blastomas, germ",
      "cell tumors, or CNS tumors) qualifying for the phase I",
      "dose-escalation study. Malignant lymphoma was excluded by",
      "protocol.",
      sep = " "
    ),
    dose_range     = paste(
      "Oral tamibarotene 4, 6, 8, 10, or 12 mg/m^2/day divided twice",
      "daily; 1-mg soft capsule formulation suitable for pediatric",
      "administration. The dosing schedule began 2-days-on / 2-days-off",
      "and was modified to 3-days-on / 1-day-off; 12 mg/m^2/day was",
      "established as the maximum dose after no adverse events were",
      "observed at the lower doses.",
      sep = " "
    ),
    regions        = paste(
      "Japan (4 medical institutions). Trial registration:",
      "UMIN-CTR identifier UMIN000017053.",
      sep = " "
    ),
    notes          = paste(
      "22 patients (15 male = 68.2%, 7 female = 32.8%) contributed 109",
      "concentration samples to the popPK dataset; one patient with an",
      "abnormally high terminal half-life was excluded from the",
      "single-dose NCA. PK sampling on day 1: predose and 2, 4, 8, 10 h",
      "post-dose; on day 14: predose only. Baseline demographics from",
      "Azechi 2024 Table 1 (Mean +/- SD): age 9.8 +/- 5.7 yr; body",
      "weight 28.4 +/- 17.2 kg; height 129.1 +/- 26.9 cm; body surface",
      "area 0.995 +/- 0.389 m^2.",
      sep = " "
    )
  )

  ini({
    # Final-model fixed-effect parameter estimates from Azechi 2024 Table 4
    # (Original Estimate column). Reference subject: BSA = 0.995 m^2 (cohort
    # mean, used as the centering value in the linear BSA-scaling form).
    # Apparent clearances in L/h, apparent volumes in L, ka in 1/h, tlag in h.
    lka   <- log(0.415);       label("Absorption rate constant ka (1/h)")                                # Azechi 2024 Table 4 final ka = 0.415 /h (95% CI 0.270-0.560; bootstrap median 0.429, 95% CI 0.350-0.546)
    ltlag <- fixed(log(0.95)); label("Absorption lag time tlag (h) -- FIXED per source")                 # Azechi 2024 Table 4 Tlag = 0.95 h FIXED (Section 2.6 / Discussion: post-covariate estimate ~1.8 h was judged of low physiological validity, so Tlag was held at the pre-covariate value of 0.95 h)
    lcl   <- log(8.73);        label("Apparent oral clearance CL/F at BSA = 0.995 m^2 (L/h)")            # Azechi 2024 Table 4 final tvCL/F = 8.73 L/h (95% CI 7.12-10.35; bootstrap median 9.1, 95% CI 7.61-10.81)
    lvc   <- log(9.17);        label("Apparent central volume V1/F at BSA = 0.995 m^2 (L)")              # Azechi 2024 Table 4 final tvV1/F = 9.17 L (95% CI 1.84-16.50; bootstrap median 10.13, 95% CI 4.47-15.40)
    lq    <- log(3.45);        label("Apparent inter-compartmental clearance Q/F (L/h)")                 # Azechi 2024 Table 4 final tvQ/F = 3.45 L/h (95% CI 1.25-5.65; bootstrap median 3.39, 95% CI 2.89-4.70). No covariate on Q/F in the final model.
    lvp   <- log(60.28);       label("Apparent peripheral volume V2/F at BSA = 0.995 m^2 (L)")           # Azechi 2024 Table 4 final tvV2/F = 60.28 L (95% CI 11.10-109.47; bootstrap median 48.64, 95% CI 28.67-68.94)

    # BSA scaling exponents -- Azechi 2024 Table 4 Final Model column states
    # CL/F = tvCL/F * BSA/mean, V1/F = tvV1/F * BSA/mean, V2/F = tvV2/F *
    # BSA/mean (no explicit power exponent, i.e. linear scaling with
    # exponent = 1). The Methods section calls the screened form an
    # "allometric function"; the final-model form retained by the authors
    # is linear in (BSA / 0.995). Encoded as fixed(1) so the linear-scaling
    # provenance is explicit in ini() (matches the Andrews 2017 tacrolimus
    # convention of encoding theory-fixed allometric exponents as
    # fixed(<exp>)).
    e_bsa_cl <- fixed(1); label("Power exponent of (BSA/0.995) on CL/F (unitless; fixed at 1 per Table 4 form)") # Azechi 2024 Table 4 Final Model: CL/F = tvCL/F * BSA/mean
    e_bsa_vc <- fixed(1); label("Power exponent of (BSA/0.995) on V1/F (unitless; fixed at 1 per Table 4 form)") # Azechi 2024 Table 4 Final Model: V1/F = tvV1/F * BSA/mean
    e_bsa_vp <- fixed(1); label("Power exponent of (BSA/0.995) on V2/F (unitless; fixed at 1 per Table 4 form)") # Azechi 2024 Table 4 Final Model: V2/F = tvV2/F * BSA/mean

    # Inter-individual variability. Azechi 2024 Methods Section 2.6
    # specifies an exponential IIV model on all five PK parameters
    # (ka, CL/F, V1/F, V2/F, Q/F) as theta_i = tvtheta * exp(eta_i).
    # No per-parameter omega magnitudes are reported in Table 4, any
    # other table, or any supplement (Europe PMC: hasSuppl = N for
    # PMC11592880). Per operator decision (sidecar request-001
    # q1 = A, 2026-06-21), the five eta terms are encoded as fixed(0)
    # so the published structural IIV declaration is preserved while
    # remaining faithful to the absence of reported variance values.
    # Stochastic VPCs built from this model show no between-subject
    # variability around the typical-value predictions; see vignette
    # Assumptions and deviations. Precedent: Chi_2018_propofol.R uses
    # the same ~ fixed(0) pattern for a paper that reported only
    # final-model THETAs without OMEGAs.
    etalka ~ fixed(0)
    etalcl ~ fixed(0)
    etalvc ~ fixed(0)
    etalvp ~ fixed(0)
    etalq  ~ fixed(0)

    # Residual error. Azechi 2024 Section 2.6 specifies a proportional
    # (exponential) residual-error model: Cobs = Cpred * (1 + eps),
    # eps ~ N(0, sigma^2). Table 4 reports "Residual variability (%) =
    # 42.4" as the single numeric, encoded here as 0.424 on the
    # proportional SD scale.
    propSd <- 0.424; label("Proportional residual SD (fraction)")                                        # Azechi 2024 Table 4 final Residual variability = 42.4%
  })

  model({
    # BSA centering value: cohort mean from Azechi 2024 Table 1 (0.995 m^2).
    bsa_norm <- BSA / 0.995

    # Individual PK parameters with the Azechi 2024 covariate equations
    # (Table 4 Final Model column). With every eta fixed at 0, exp(eta) = 1
    # and the typical-value structural form is reproduced exactly.
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl) * bsa_norm ^ e_bsa_cl
    vc   <- exp(lvc + etalvc) * bsa_norm ^ e_bsa_vc
    q    <- exp(lq  + etalq)
    vp   <- exp(lvp + etalvp) * bsa_norm ^ e_bsa_vp
    tlag <- exp(ltlag)

    # Two-compartment oral disposition. Dose lands in `depot`; the
    # bioavailability F is absorbed into the apparent CL/F and V/F
    # parameterisation (no separate F1 term reported).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tamibarotene serum concentrations reported in ng/mL. With dose in mg
    # and vc in L, central/vc is in mg/L = ug/mL; multiply by 1000 to
    # convert to ng/mL (matches the units in Table 2 NCA results, e.g.,
    # Cmax range 52.7-259 ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
