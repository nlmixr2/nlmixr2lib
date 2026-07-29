Bienczak_2025_ligelizumab <- function() {
  description <- "Two-compartment population PK model for ligelizumab in adolescent and adult patients with chronic spontaneous urticaria and healthy adult volunteers (Bienczak 2025)"
  reference <- "Bienczak A, Gautier A, Hua E, Ji Y, Scosyrev E, Smeets S, Severin T, Drollmann A, Patekar M, Savelieva M. Model-Informed Drug Development for Ligelizumab in Patients With Chronic Spontaneous Urticaria. CPT Pharmacometrics Syst Pharmacol. 2025. doi:10.1002/psp4.70090"
  vignette <- "Bienczak_2025_ligelizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL/F, Q/F (exponent 0.993, fixed) and Vc/F, Vp/F (exponent 0.597, fixed), centered at 70 kg.",
      source_name        = "WT"
    ),
    IGE = list(
      description        = "Baseline serum total immunoglobulin E concentration",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL/F (exponent 0.106) and Vp/F (exponent -0.0816), centered on 90 IU/mL. Source paper reports IgE in IU/mL (Table S6 footnote a); the per-model units are kept as IU/mL to preserve the source reference value rather than converting to the canonical ng/mL convention.",
      source_name        = "IgE"
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody status (1 if subject had at least one ADA-positive sample at any time in the study, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative throughout)",
      notes              = "Time-fixed (ever-positive); multiplicative effect on CL/F (+27.5%, exp(0.243)) and Vp/F (-40.9%, exp(-0.526)). Renamed from source column ADA to the canonical ADA_POS per covariate-columns.md.",
      source_name        = "ADA"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with chronic spontaneous urticaria)",
      notes              = "1 = pooled healthy-volunteer cohort from studies A2103 and C2101; 0 = CSU patient. Multiplicative effect on CL/F (exp(-0.087); ~8.3% lower).",
      source_name        = "healthy_volunteer"
    ),
    STUDY_C2201 = list(
      description        = "Study C2201 cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other study in the Bienczak 2025 pooled PopPK analysis: A2103, C2101, C2202, C2302, or C2303)",
      notes              = "1 = subject enrolled in study C2201 (NCT02477332; Phase 2b ligelizumab in adult CSU patients). Multiplicative effect on CL/F (exp(0.176); ~19.2% higher). Time-fixed (subject-level).",
      source_name        = "C2201"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1907,
    n_studies      = 6,
    age_range      = "12-80 years (adolescent and adult)",
    weight_range   = "31.0-181.3 kg",
    sex_female_pct = 65.6,
    race_ethnicity = c(White = 73, Asian = 21, Black = 2, Other = 4),
    disease_state  = "Adolescent and adult patients with chronic spontaneous urticaria (n=1706: 113 adolescent + 1593 adult); pooled with adult healthy volunteers (n=201) for absorption and disposition.",
    dose_range     = "0.2-4.0 mg/kg q2w sc, 24-240 mg q4w sc, and 120 mg single dose sc across studies A2103, C2101, C2201, C2202, C2302, and C2303.",
    regions        = "Global (multi-region Phase 1, 2b, and 3 program)",
    notes          = "Baseline demographics from Supplementary Information S1, Tables S4 and S5 of Bienczak 2025."
  )

  ini({
    # Structural parameters - typical values for a 70 kg, ADA-negative, CSU patient
    # with IgE = 90 IU/mL, not enrolled in study C2201 (Table S6 of Supplementary Information S1)
    lka  <- log(0.218);   label("Absorption rate constant (ka, 1/day)")                                # Table S6, structural ka
    lcl  <- log(0.602);   label("Apparent clearance (CL/F, L/day) for a 70 kg CSU patient")            # Table S6, structural CL/F
    lvc  <- log(5.47);    label("Apparent central volume of distribution (Vc/F, L) for a 70 kg CSU patient")  # Table S6, structural Vc/F
    lq   <- log(0.969);   label("Apparent inter-compartmental clearance (Q/F, L/day) for a 70 kg CSU patient")  # Table S6, structural Q/F
    lvp  <- log(10.2);    label("Apparent peripheral volume of distribution (Vp/F, L) for a 70 kg CSU patient")  # Table S6, structural Vp/F

    # Allometric exponents (shared between CL+Q and Vc+Vp)
    # Footnote d of Table S6: "Fixed to the values estimated in the final structural model."
    e_wt_cl_q  <- fixed(0.993); label("Shared allometric exponent on CL/F and Q/F (unitless)")  # Table S6 weight on CL/F and Q/F, footnote d (fixed)
    e_wt_vc_vp <- fixed(0.597); label("Shared allometric exponent on Vc/F and Vp/F (unitless)") # Table S6 weight on Vc/F and Vp/F, footnote d (fixed)

    # Power covariate effects on continuous IgE (centered at 90 IU/mL per Table S6 footnote a)
    e_ige_cl <-  0.106;   label("Power exponent of baseline IgE on CL/F (unitless)")  # Table S6 IgE on CL/F
    e_ige_vp <- -0.0816;  label("Power exponent of baseline IgE on Vp/F (unitless)")  # Table S6 IgE on Vp/F

    # Categorical / binary covariate effects (log-additive on the structural parameter)
    e_ada_pos_cl     <-  0.243;  label("Effect of ever-positive ADA on CL/F (log-additive)")          # Table S6 ADA-positive on CL/F
    e_ada_pos_vp     <- -0.526;  label("Effect of ever-positive ADA on Vp/F (log-additive)")          # Table S6 ADA-positive on Vp/F
    e_dis_healthy_cl <- -0.087;  label("Effect of healthy-volunteer status on CL/F (log-additive)")   # Table S6 healthy-volunteer status on CL/F
    e_study_c2201_cl <-  0.176;  label("Effect of study C2201 enrollment on CL/F (log-additive)")     # Table S6 study C2201 on CL/F

    # IIV: variances on the log-scale eta of each structural parameter
    # CL/F and Vp/F are correlated (r = 0.484) per Table S6
    # var_cl    = 0.335^2  = 0.112225
    # var_vp    = 0.39^2   = 0.1521
    # cov_cl_vp = 0.484 * 0.335 * 0.39 = 0.06323286
    etalcl + etalvp ~ c(0.112225,
                        0.06323286, 0.1521)              # Table S6 IIV SD for CL/F and Vp/F and their correlation
    etalka          ~ 2.0164e-6                          # Table S6 IIV SD for ka = 0.00142 (100% shrinkage); var = 0.00142^2 ~ 0
    etalvc          ~ 0.837225                           # Table S6 IIV SD for Vc/F = 0.915 -> var = 0.837225
    etalq           ~ 0.067081                           # Table S6 IIV SD for Q/F = 0.259 -> var = 0.067081

    # Combined residual error: additive (a) = 142 ng/mL = 0.142 ug/mL, proportional (b) = 0.177
    addSd  <- 0.142;  label("Additive residual error (ug/mL)")           # Table S6 additive error (a) = 142 ng/mL
    propSd <- 0.177;  label("Proportional residual error (fraction)")    # Table S6 proportional error (b)
  })
  model({
    # Covariate multipliers
    # Power form for body weight (centered 70 kg) and baseline IgE (centered 90 IU/mL)
    wt_cl_q  <- (WT  / 70)^e_wt_cl_q
    wt_vc_vp <- (WT  / 70)^e_wt_vc_vp
    ige_cl   <- (IGE / 90)^e_ige_cl
    ige_vp   <- (IGE / 90)^e_ige_vp

    # Log-additive form for ADA-positive, healthy-volunteer, and study C2201 indicators
    cov_cl <- exp(e_ada_pos_cl * ADA_POS + e_dis_healthy_cl * DIS_HEALTHY + e_study_c2201_cl * STUDY_C2201)
    cov_vp <- exp(e_ada_pos_vp * ADA_POS)

    # Individual PK parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * wt_cl_q  * ige_cl * cov_cl
    vc <- exp(lvc + etalvc) * wt_vc_vp
    q  <- exp(lq  + etalq)  * wt_cl_q
    vp <- exp(lvp + etalvp) * wt_vc_vp * ige_vp * cov_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
