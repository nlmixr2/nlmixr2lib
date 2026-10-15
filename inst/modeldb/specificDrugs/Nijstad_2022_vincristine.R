Nijstad_2022_vincristine <- function() {
  description <- paste(
    "Semi-mechanistic population PK model for intravenous vincristine in",
    "children, adolescents and young adults (Nijstad 2022; n = 206, 0.04-33.9",
    "years, 1297 plasma concentrations). Two-compartment linear disposition",
    "(CL, Vc, Q, Vp; allometric on body weight, exponents fixed at 0.75 and 1,",
    "70 kg reference) plus a third, saturable compartment for vincristine",
    "bound to beta-tubulin, filled from the central amount at",
    "kon x (1 - bound/Bmax) and emptied at koff. The binding capacity Bmax",
    "scales with (WT/70)^1 x (AGE/18)^-0.199, so younger children carry more",
    "binding capacity per kg. IIV on CL, Q, Vc, Vp, kon and koff;",
    "inter-occasion variability on Bmax (one occasion per dose);",
    "proportional residual error."
  )
  reference <- paste(
    "Nijstad AL, Chu WY, de Vos-Kerkhof E, Enters-Weijnen CF,",
    "van de Velde ME, Kaspers GJL, Barnett S, Veal GJ, Lalmohamed A,",
    "Zwaan CM, Huitema ADR. A Population Pharmacokinetic Modelling Approach",
    "to Unravel the Complex Pharmacokinetics of Vincristine in Children.",
    "Pharm Res. 2022;39(10):2487-2495. doi:10.1007/s11095-022-03364-1"
  )
  vignette <- "Nijstad_2022_vincristine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "vincristine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vincristine", units = "mg", specimen = "tissue", verified = TRUE),
    complex = list(analyte = "vincristine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source term BW. Allometric size descriptor on CL and Q (exponent",
        "0.75) and on Vc, Vp and Bmax (exponent 1), all fixed and referenced",
        "to 70 kg (Nijstad 2022 Methods 'Covariate Analysis' and Results",
        "'Covariate Analysis'). Cohort median 27.1 kg (range 2.9-126.0) per",
        "Table I. Missing BW (2 UK patients) was imputed by the authors from",
        "age and BSA with UK growth charts and the Du Bois equation."
      ),
      source_name = "BW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on Bmax only, normalised to 18 years:",
        "(AGE/18)^-0.199 (Nijstad 2022 Results 'Covariate Analysis' and",
        "Table II). Cohort median 8.3 years (range 0.04-33.9) per Table I;",
        "25 patients were younger than 1 year. The age term grows without",
        "bound as AGE approaches 0 (3.4-fold at 2 weeks), so supply a",
        "positive postnatal age; the cohort's youngest patient was about",
        "2 weeks old. Missing age (8 UK patients) was imputed by the authors",
        "from BW and height with UK growth charts."
      ),
      source_name = "age"
    ),
    OCC = list(
      description = "Integer-valued occasion index (one occasion per vincristine dose, 1-5)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Nijstad 2022 Methods 'Model Development': inter-occasion",
        "variability 'was implemented similarly as IIV, with each dose and",
        "subsequent sampling defined as a separate occasion'. Patients",
        "contributed 1-5 occasions (Table I), so OCC takes values 1-5 and is",
        "decomposed in model() into binary indicators that multiplex the",
        "per-occasion IOV etas on Bmax. Observation records carry the",
        "occasion of the dose that preceded them. OCC = 0 or any value",
        "outside 1-5 zeros every indicator and gives the IIV-only Bmax."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    PLT = list(
      description = "Platelet (thrombocyte) count",
      units = "10^9/L",
      type = "continuous",
      notes = paste(
        "Tested as a power covariate on Bmax normalised to 300 x 10^9/L, with",
        "300 x 10^9/L imputed for the 46% of occasions without a count. Not",
        "retained: the models were unstable, the IOV on Bmax did not fall",
        "and the IIV on other parameters rose (Nijstad 2022 Results",
        "'Covariate Analysis')."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 206L,
    n_studies = 5L,
    n_pk_samples = 1297L,
    n_occasions = 253L,
    age_range = "0.04-33.9 years",
    age_median = "8.3 years",
    weight_range = "2.9-126.0 kg",
    weight_median = "27.1 kg",
    sex_female_pct = 48,
    race_ethnicity = "Not reported in source",
    disease_state = paste(
      "Children and young adults with cancer receiving vincristine as",
      "standard of care: an unrestricted paediatric oncology cohort",
      "(Princess Maxima Center, the Netherlands; Down syndrome excluded),",
      "Ewing sarcoma patients up to 24 years (20 UK centres; GFR < 60",
      "mL/min/1.73 m^2 excluded), and three historical paediatric cohorts",
      "(Lee et al., UK patients only, n = 24; van de Velde et al., n = 37;",
      "Barnett et al., n = 26). 25 patients were younger than 1 year",
      "(7 aged 0-3 months, 8 aged 3-6 months, 4 aged 6-9 months, 6 aged",
      "9-12 months)."
    ),
    dose_range = paste(
      "1-2 mg/m^2 capped at 2 mg, with infant reductions per local protocol;",
      "median 1.6 mg (0.1-2.0), 1.4 mg/m^2 (0.4-2.5), 0.05 mg/kg",
      "(0.02-0.09). Given as an IV bolus (214 occasions) or a 15-113 min",
      "infusion (39 occasions)."
    ),
    regions = "the Netherlands and the United Kingdom",
    sampling_window = "4-8 samples per patient (median 5 per occasion, range 1-8), up to about 75 h after the dose; 30 of 1297 samples below the LLOQ (0.10, 0.25 or 0.50 ng/mL by assay), the first of which were included at half the LLOQ.",
    notes = "NONMEM 7.3.0 with FOCE-I; parameter precision by sampling importance resampling. Platelet counts were available for 137 of 253 occasions (median 224 x 10^9/L, range 5-1063)."
  )

  ini({
    # Structural parameters: Nijstad 2022 Table II, for a 70 kg subject (Bmax
    # also for an 18-year-old). Allometric exponents 0.75 (CL, Q) and 1 (Vc,
    # Vp, Bmax) are fixed a priori (Methods 'Covariate Analysis').
    lcl <- log(30.6) ; label("Clearance CL at 70 kg (L/h)")                                     # Table II 'CL 70kg (L/h)' = 30.6 (95% CI 27.6-33.0)
    lq <- log(63.2) ; label("Intercompartmental clearance Q at 70 kg (L/h)")                     # Table II 'Q 70kg (L/h)' = 63.2 (95% CI 57.2-70.1)
    lvc <- log(5.39) ; label("Central volume of distribution Vc at 70 kg (L)")                   # Table II 'Vc 70kg (L)' = 5.39 (95% CI 4.23-6.46)
    lvp <- log(400) ; label("Peripheral volume of distribution Vp at 70 kg (L)")                 # Table II 'Vp 70kg (L)' = 400 (95% CI 357-463)
    lbmax <- log(0.525) ; label("Maximal beta-tubulin binding capacity Bmax at 70 kg and 18 years (mg)") # Table II 'Bmax 18yrs,70 kg (mg)' = 0.525 (95% CI 0.479-0.602)
    lkon <- fixed(log(1300)) ; label("First-order association rate constant kon (1/h)")         # Table II 'k on (/h)' = 1300 fixed; Results 'Model Development' fixed at the value with the lowest OFV
    lkoff <- log(11.5) ; label("First-order dissociation rate constant koff (1/h)")              # Table II 'k off (/h)' = 11.5 (95% CI 9.2-14.5)
    e_age_bmax <- -0.199 ; label("Power exponent of (AGE/18) on Bmax (unitless)")              # Table II 'Age on Bmax' = -0.199 (95% CI -0.304 to -0.090)

    # Inter-individual variability, exponential (Methods Eq. 2:
    # Pi = Ppop x exp(eta_i)). Table II prints percentages; they are read as
    # CV% and converted with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.205003    # Table II 'IIV CL (%)' = 47.7 (95% CI 41.0-54.3)
    etalq ~ 0.135545     # Table II 'IIV Q (%)' = 38.1 (95% CI 26.2-49.0)
    etalvc ~ 0.916541    # Table II 'IIV Vc (%)' = 122.5 (95% CI 98.7-158.3)
    etalvp ~ 0.282198    # Table II 'IIV Vp (%)' = 57.1 (95% CI 48.8-69.7)
    etalkon ~ 0.955598   # Table II 'IIV k on (%)' = 126.5 (95% CI 108.7-147.8)
    etalkoff ~ 0.0564569 # Table II 'IIV k off (%)' = 24.1 (95% CI 11.1-33.8)

    # Inter-occasion variability on Bmax, one occasion per dose (Methods
    # 'Model Development'). One IOV variance shared by the five occasions:
    # occasion 1 carries the estimate and occasions 2-5 are fixed equal to it
    # (the equivalent of NONMEM $OMEGA BLOCK(1) SAME).
    etaiov_bmax_1 ~ 0.299572        # Table II 'IOV Bmax (%)' = 59.1 (95% CI 50.7-66.1)
    etaiov_bmax_2 ~ fixed(0.299572) # equal to the occasion-1 IOV variance on Bmax
    etaiov_bmax_3 ~ fixed(0.299572) # equal to the occasion-1 IOV variance on Bmax
    etaiov_bmax_4 ~ fixed(0.299572) # equal to the occasion-1 IOV variance on Bmax
    etaiov_bmax_5 ~ fixed(0.299572) # equal to the occasion-1 IOV variance on Bmax

    # Residual error: proportional (Table II).
    propSd <- 0.301 ; label("Proportional residual error (fraction)") # Table II 'Proportional residual error (%)' = 30.1 (95% CI 28.9-31.4)
  })

  model({
    # Occasion indicators for the IOV multiplexing (one occasion per dose).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_bmax <- oc1 * etaiov_bmax_1 + oc2 * etaiov_bmax_2 + oc3 * etaiov_bmax_3 +
      oc4 * etaiov_bmax_4 + oc5 * etaiov_bmax_5

    # Individual parameters. Bmax scales with body weight like a volume and
    # with age by a power function normalised to 18 years (Results 'Covariate
    # Analysis'). kon and koff are rate constants, which the paper's
    # allometric rule (exponents for clearances and volumes only) leaves
    # unscaled.
    cl <- exp(lcl + etalcl) * (WT / 70)^0.75
    q <- exp(lq + etalq) * (WT / 70)^0.75
    vc <- exp(lvc + etalvc) * (WT / 70)
    vp <- exp(lvp + etalvp) * (WT / 70)
    bmax <- exp(lbmax + iov_bmax) * (WT / 70) * (AGE / 18)^e_age_bmax
    kon <- exp(lkon + etalkon)
    koff <- exp(lkoff + etalkoff)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Supplementary Table S1. Binding is written in amounts: the central
    # amount fills the beta-tubulin pool at kon x (1 - bound/Bmax) and the
    # bound amount returns at koff (Methods Eq. 1).
    binding <- kon * central * (1 - complex / bmax) - koff * complex
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - binding
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(complex) <- binding

    # Dose in mg over volume in L gives mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
