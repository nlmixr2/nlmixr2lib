Pohl_2022_linzagolix <- function() {
  description <- "Two-compartment population pharmacokinetic model for linzagolix, an oral gonadotropin-releasing hormone (GnRH) receptor antagonist developed for endometriosis, in healthy women and women with endometriosis (Pohl 2022). Sequential zero-order then first-order oral absorption; fixed allometric body-weight scaling on apparent clearance (0.75) and apparent central volume (1) with a 58 kg reference; a categorical non-Caucasian effect on apparent clearance; and a proportional residual error whose variance differs between the EDELWEISS phase 2b study and the phase 1 studies and additionally carries a subject-level random effect."
  reference <- "Pohl O, Baron K, Riggs M, French J, Garcia R, Gotteland J-P. A model-based analysis to guide gonadotropin-releasing hormone receptor antagonist use for management of endometriosis. British Journal of Clinical Pharmacology. 2022;88(5):2359-2371. doi:10.1111/bcp.15171"
  vignette <- "Pohl_2022_linzagolix"

  # Non-canonical residual-SD and eta names carried by this paper: the residual
  # error variance is study-specific AND carries its own subject-level random
  # effect (Table 4 "IIV-sigma2"). Same shape as Friberg_2012_voriconazole.
  paper_specific_residual_sds <- c("propSdEdelweiss", "propSdOther")
  paper_specific_etas <- c("eta_re")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F (exponent 0.75) and V2/F (exponent 1), both fixed to their allometric values rather than estimated, normalised to a 58 kg reference subject (Table 4 covariate-effect rows 'CL/F ~ (weight 58 kg)' and 'V2/F ~ (weight 58 kg)'; Supporting Information section 4.1). Q/F and V3/F carry no weight effect in the final model.",
      source_name        = "WT"
    ),
    RACE_WHITE = list(
      description        = "Caucasian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (Caucasian; the typical-value reference in this model)",
      notes              = "Source paper dichotomises race as Caucasian vs non-Caucasian (Table 3 'Percent Caucasian'). The Caucasian subgroup is the typical-value reference, so the effect is implemented on (1 - RACE_WHITE): non-Caucasian subjects have 8% higher CL/F (Table 4 'CL/F ~ non-Caucasian' = 1.08). Same reference-category orientation as Hu_2014_bapineuzumab.R.",
      source_name        = "race group (Caucasian / non-Caucasian)"
    ),
    STUDY_EDELWEISS = list(
      description        = "EDELWEISS phase 2b study indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the phase 1 studies KLH1101, 16-OBE2109-011, 17-OBE2109-008 and the phase 1 study KLH1204)",
      notes              = "1 = observation from the EDELWEISS phase 2b dose-ranging study (NCT02778399), which contributed sparse PK; 0 = observation from any other pooled study, all of which contributed rich PK. Selects between the two proportional residual error magnitudes of Table 4 ('EDELWEISS data' vs 'All other studies'). Per-record study-fixed indicator.",
      source_name        = "study"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "linzagolix", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "linzagolix", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "linzagolix", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 756,
    n_studies      = 5,
    age_range      = "18-48 years",
    age_median     = "32-40 years across the contributing studies (Table 3)",
    weight_range   = "median 53.9-65.5 kg across the contributing studies (Table 3)",
    weight_median  = "58 kg (model reference weight)",
    sex_female_pct = 100,
    race_ethnicity = "0-100% Caucasian by study (Table 3); pooled analysis deliberately retained both Caucasian and non-Caucasian subjects to reduce parameter uncertainty",
    disease_state  = "endometriosis, or endometriosis with co-existing uterine fibroids, plus healthy pre- and postmenopausal volunteers (approximately 24% of subjects and 55% of observations were from healthy volunteers)",
    dose_range     = "12.5-400 mg single dose; 25-200 mg once daily for 42 days to 24 weeks",
    regions        = "Europe and Japan",
    notes          = "4250 linzagolix concentration observations from 756 subjects pooled across five clinical trials (KLH1101 phase 1 SAD/MAD, EDELWEISS phase 2b, 16-OBE2109-011, 17-OBE2109-008, KLH1204). Subjects receiving placebo contributed no records to the population PK analysis set. About 5% of concentrations were below the quantitation limit and were dropped (Supporting Information section 1.1). Baseline demographics are in Table 3."
  )

  ini({
    # Structural parameters -- typical values for a 58 kg Caucasian subject
    lcl <- log(0.422); label("Apparent clearance CL/F (L/h)")                            # Table 4: CL/F 0.422 (95% CI 0.393-0.455)
    lvc <- log(5.13);  label("Apparent central volume of distribution V2/F (L)")         # Table 4: V2/F 5.13 (95% CI 4.19-6.18)
    lq  <- log(0.168); label("Apparent intercompartmental clearance Q/F (L/h)")          # Table 4: Q/F 0.168 (95% CI 0.130-0.225)
    lvp <- log(3.12);  label("Apparent peripheral volume of distribution V3/F (L)")      # Table 4: V3/F 3.12 (95% CI 2.83-3.41)
    lka <- log(2.49);  label("First-order absorption rate constant (1/h)")               # Table 4: KA 2.49 (95% CI 2.04-3.08); the Table 4 unit "L/h" is a typographical slip, KA is a first-order rate constant
    ld1 <- log(0.644); label("Duration of the zero-order input into the depot (h)")      # Table 4: D1 0.644 (95% CI 0.314-1.24)

    # Allometric exponents -- fixed to their allometric values, not estimated
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")              # Table 4: "CL/F ~ (weight 58 kg)" 0.75 FIXED; Supporting Information 4.1
    e_wt_vc <- fixed(1.00); label("Allometric exponent on V2/F (unitless)")              # Table 4: "V2/F ~ (weight 58 kg)" 1.00 FIXED; Supporting Information 4.1

    # Covariate effects
    e_nonwhite_cl <- 0.08; label("Fractional increase in CL/F for non-Caucasian vs Caucasian reference (unitless)")  # Table 4: "CL/F ~ non-Caucasian" ratio 1.08 (95% CI 1.05-1.12), i.e. a fractional increase of 0.08

    # Interindividual variability -- log-normal; Table 4 reports variances
    etalcl ~ 0.0354   # Table 4: IIV-CL/F 0.0354 (95% CI 0.0271-0.0498), shrinkage 16.5%
    etalvc ~ 0.0444   # Table 4: IIV-V2/F 0.0444 (95% CI 0.0203-0.115), shrinkage 62.0%
    etald1 ~ 0.510    # Table 4: IIV-D1 0.510 (95% CI 0.230-1.41), shrinkage 46.2%
    eta_re ~ 0.764    # Table 4: IIV-sigma2 0.764 (95% CI 0.505-1.11), shrinkage 24.8%; subject-level random effect on the residual error VARIANCE

    # Residual error -- proportional, study-specific. Table 4 reports variances,
    # so each SD below is the square root of the tabulated value.
    propSdEdelweiss <- 0.3435; label("Proportional residual SD for EDELWEISS phase 2b data (fraction)")  # Table 4: "EDELWEISS data" variance 0.118 (95% CI 0.0698-0.206); sqrt(0.118) = 0.3435
    propSdOther     <- 0.1972; label("Proportional residual SD for all other studies (fraction)")        # Table 4: "All other studies" variance 0.0389 (95% CI 0.0309-0.0502); sqrt(0.0389) = 0.1972
  })

  model({
    # 1. Individual parameters. Weight is normalised to the 58 kg reference
    #    subject; Caucasian is the typical-value reference so the race effect
    #    enters on (1 - RACE_WHITE).
    cl <- exp(lcl + etalcl) * (WT / 58)^e_wt_cl * (1 + e_nonwhite_cl * (1 - RACE_WHITE))
    vc <- exp(lvc + etalvc) * (WT / 58)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Absorption is sequential zero-order then first-order:
    #    the dose enters the depot as a zero-order input of duration D1 and
    #    then transfers to the central compartment with rate constant ka.
    #    Dose the depot with rate = -2 so rxode2 applies the modelled dur().
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot) <- d1

    # 4. Observation and error. Volumes are apparent (V/F), so with mg doses
    #    and litre volumes Cc is in mg/L = ug/mL.
    Cc <- central / vc

    # The residual error variance is study-specific and additionally carries a
    # subject-level random effect (Table 4 "IIV-sigma2", described in Supporting
    # Information 3.1 as "a subject-level random effect on the residual error
    # variance"). eta_re is therefore on the log-VARIANCE scale, so it enters
    # the SD as exp(eta_re / 2).
    propSdInd <-
      (propSdEdelweiss * STUDY_EDELWEISS + propSdOther * (1 - STUDY_EDELWEISS)) *
      exp(eta_re / 2)

    Cc ~ prop(propSdInd)
  })
}
