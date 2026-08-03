Watson_2010_liraglutide <- function() {
  description <- paste(
    "One-compartment population PK model for subcutaneous liraglutide",
    "(once-daily human GLP-1 analog) with sequential zero-order-then-",
    "first-order (sequential dual) absorption, pooled across four Phase 1",
    "studies (Studies A-D) in healthy volunteers and adults with type 2",
    "diabetes (Watson 2010). The combined-data final model estimates",
    "apparent clearance and volume (CL/F, V/F) on a per-kg body-weight",
    "basis. Absolute bioavailability F=0.51 was identifiable only in the",
    "Study A submodel (with IV data) and is documented in the vignette."
  )
  reference <- paste(
    "Watson E, Jonker DM, Jacobsen LV, Ingwersen SH. Population",
    "pharmacokinetics of liraglutide, a once-daily human glucagon-like",
    "peptide-1 analog, in healthy volunteers and subjects with type 2",
    "diabetes, and comparison to twice-daily exenatide. J Clin Pharmacol.",
    "2010;50(8):886-894. doi:10.1177/0091270009354996"
  )
  vignette <- "Watson_2010_liraglutide"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L (nM)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Watson 2010 reports apparent clearance and volume on a per-kg",
        "basis (CL/F in L/h/kg, V/F in L/kg). WT enters as a linear",
        "(allometric-exponent-1) multiplier on CL/F and V/F relative to",
        "the paper's per-kg (1 kg) reference. The paper does not tabulate",
        "cohort baseline weight; typical adult 70-90 kg is a reasonable",
        "simulation range. Time-fixed at baseline."
      ),
      source_name        = "WT"
    ),
    DIS_DIAB = list(
      description        = "Type 2 diabetes status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = paste(
        "Watson 2010 (Results, Single Dose in HV and Subjects With T2D)",
        "tested the statistical significance of the HV-vs-T2D disease-",
        "status covariate on all model parameters one at a time; the",
        "only significant differences (P<.001, chi-square test) were on",
        "the duration of the zero-order absorption phase (T0) and the",
        "first-order absorption rate constant (kA). Encoded as additive",
        "shifts on log(T0) and log(kA) so DIS_DIAB=0 yields the HV",
        "typical values (T0=6.0 h, kA=0.104 /h) and DIS_DIAB=1 yields",
        "the T2D typical values (T0=8.0 h, kA=0.154 /h) per Table II",
        "(All Studies columns). CL/F and V/F did not differ significantly",
        "between HV and T2D so no DIS_DIAB effect is applied there.",
        "Source column 'T2D' in the paper's dataset maps to this",
        "canonical DIS_DIAB indicator (Type 1 vs Type 2 diabetes are",
        "pooled at the column level per the canonical register; the",
        "Watson 2010 cohort is Type 2 only)."
      ),
      source_name        = "T2D"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 98L,
    n_studies      = 4L,
    age_range      = "Adults (Watson 2010 does not tabulate a pooled age range; individual ranges are in the source studies A-D references)",
    weight_range   = "Adults (Watson 2010 does not tabulate a pooled weight range; per-kg parameterisation scales linearly with WT)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not tabulated in Watson 2010",
    disease_state  = paste(
      "Healthy volunteers (Study A n=54; Study B n=20 HV) pooled with",
      "adults with type 2 diabetes (Study B n=2 T2D; Study C n=11;",
      "Study D n=11). Studies A and B included fasted single-dose SC",
      "arms; Studies C and D were single-dose SC crossover trials."
    ),
    dose_range     = paste(
      "SC single dose 1.25-20 ug/kg (Study A: 1.25, 2.5, 5.0, 10.0, 12.5,",
      "15.0, 17.5, or 20.0 ug/kg, n=8 per level); SC single or once-daily",
      "for 5 days 1.25-12.5 ug/kg (Study B); SC single dose 10 ug/kg",
      "(Study C) or 7.5 ug/kg (Study D). Study A also included a 5",
      "ug/kg SC IV arm (n=8) used only in the initial Study-A-only",
      "submodel. Molar equivalents via MW 3751.2 g/mol."
    ),
    regions        = "Not tabulated in Watson 2010; source studies were Novo Nordisk single-country Phase 1 trials (Denmark).",
    trials         = c("Elbrond 2002 (Study A)",
                       "Agerso 2002 (Study B, HV + 2 T2D unpublished)",
                       "Juhl 2002 (Study C)",
                       "Nauck 2003 (Study D)"),
    notes          = paste(
      "The final combined model was fit to Studies A-D SC-only data per",
      "Watson 2010 Methods ('reestimated on the combined pharmacokinetic",
      "data from studies A through D ... exposure after subcutaneous",
      "injection only'). Eight outliers were excluded during model",
      "building: 5 from Study A (low SC dosing volume; 3 in 1.25 ug/kg,",
      "2 in 5 ug/kg) and 3 from Study D (7.5 ug/kg with dose-normalised",
      "Cmax about 2x other subjects). Study E (Vilsboll 2007, n=168 T2D,",
      "SC OD 0.65-1.9 mg for up to 10 weeks) provided sparse steady-",
      "state samples that were used only to verify CL/F after fixing",
      "kA, T0, and V/F to the Studies A-D estimates; it is not part of",
      "the parameters encoded here. Study F (Calara 2005, n=28) was an",
      "exenatide study used only for the peak-to-trough comparison and",
      "is not encoded (see Cirincione_2017_exenatide.R for a proper",
      "popPK exenatide model)."
    )
  )

  ini({
    # Structural parameters from Watson 2010 Table II ('All Studies'
    # columns; final combined model on Studies A-D SC-only data). CL/F
    # and V/F are apparent because absolute F was only estimable in
    # Study A (with IV data). Doses are entered in nmol so that
    # Cc = central / vc directly returns nmol/L (= nM), matching the
    # residual-error scale (Table II additive = 0.042 nM).
    lcl <- log(0.013);  label("Apparent clearance per kg body weight (CL/F, L/h/kg)")     # Table II, All Studies: 0.013 L/h/kg (rSE 5%)
    lvc <- log(0.16);   label("Apparent central volume per kg body weight (V/F, L/kg)")   # Table II, All Studies: 0.16 L/kg (rSE 42%)
    lka <- log(0.104);  label("First-order absorption rate constant in HV (kA, 1/h)")     # Table II, All Studies HV: 0.104 /h (rSE 71%)
    lt0 <- log(6.0);    label("Duration of the zero-order absorption phase in HV (T0, h)") # Table II, All Studies HV: 6.0 h (rSE 0.1%)

    # T2D disease-status effects on absorption (Watson 2010 Results:
    # 'The only significant differences were found on the absorption
    # parameters (Table II) yielding a slightly longer duration of the
    # zero-order phase in subjects with T2D (P<.001, chi-square test),
    # as well as slightly higher (P<.001, chi-square test) absorption
    # rate.'). Encoded as additive log-scale shifts so DIS_DIAB=0 gives
    # the HV typical values and DIS_DIAB=1 gives the T2D typical values.
    # Magnitudes derived from the reported HV and T2D estimates.
    e_dis_diab_lt0 <- log(8.0 / 6.0);     label("Additive shift in log(T0) for T2D vs HV (T0_T2D=8.0 h vs T0_HV=6.0 h)")   # Table II, All Studies T2D: 8.0 h (rSE 0.2%)
    e_dis_diab_lka <- log(0.154 / 0.104); label("Additive shift in log(kA) for T2D vs HV (kA_T2D=0.154 /h vs kA_HV=0.104 /h)") # Table II, All Studies T2D: 0.154 /h (rSE 34%)

    # IIV -- Watson 2010 Methods ('individual parameter estimate P_i
    # equals P_pop * exp(eta), where eta is normal with mean zero and
    # variance omega^2'; log-normal). Table II reports BSV as CV%; the
    # log-normal variance is omega^2 = log(1 + CV^2):
    #   BSV CL/F = 37% -> omega^2 = log(1 + 0.37^2) = 0.12833
    #   BSV V/F  = 47% -> omega^2 = log(1 + 0.47^2) = 0.19960
    # No IIV was reported on kA or T0 in the combined model.
    etalcl ~ 0.12833  # Table II, All Studies: BSV CL/F = 37% (rSE 16%)
    etalvc ~ 0.19960  # Table II, All Studies: BSV V/F = 47% (rSE 33%)

    # Residual error (Watson 2010 Methods 'A proportional residual error
    # model as well as the combination of an additive and proportional
    # residual error model was tested'; final combined model uses both).
    propSd <- 0.154;  label("Proportional residual error (fraction)")           # Table II, All Studies: 15.4% (rSE 5%)
    addSd  <- 0.042;  label("Additive residual error (nmol/L)")                 # Table II, All Studies: 0.042 nM (rSE 4%)
  })

  model({
    # Individual PK parameters. The paper's per-kg parameterisation is
    # preserved by multiplying the per-kg estimates (from ini() lcl and
    # lvc) by body weight WT (kg) to obtain absolute L/h and L. IIV
    # (log-normal) is applied on the per-kg scale so the fractional CV
    # is invariant under WT scaling. The DIS_DIAB effect on kA and T0
    # is a categorical log-scale shift; DIS_DIAB=0 (HV) is the reference.
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT
    ka <- exp(lka + e_dis_diab_lka * DIS_DIAB)
    t0 <- exp(lt0 + e_dis_diab_lt0 * DIS_DIAB)

    kel <- cl / vc

    # Sequential zero-order-then-first-order absorption (Watson 2010
    # Results, 'Single Dose in HV (Study A)'). Equations 1 and 2:
    #   Equation 1 (t <= T0):  dA1/dT = -D * F * kA / (1 + kA * T0)   (constant rate)
    #   Equation 2 (t  > T0):  dA1/dT = -kA * A1
    # This parameterisation ensures dA1/dT is continuous at T0: at t=T0,
    # A1 = D * F * (1 - kA*T0/(1+kA*T0)) = D * F / (1 + kA * T0), so
    # kA * A1(T0) = D*F*kA/(1+kA*T0) = the zero-order rate. In the
    # combined model, D and F are not separately identifiable (only D*F
    # via CL/F and V/F); the dose amount placed into the depot is
    # D*F, i.e., podo(depot). mtime() inserts an ODE integration step
    # exactly at T0 so the piecewise absorption dynamics are captured;
    # tad(depot) is time-after-dose relative to the last depot dose.
    mtime(tswitch) <- t0

    absorption_rate <- ka * depot
    if (tad(depot) <= tswitch) absorption_rate <- podo(depot) * ka / (1 + ka * tswitch)

    d/dt(depot)   <- -absorption_rate
    d/dt(central) <-  absorption_rate - kel * central

    # Observation. Dose enters the depot in nmol; vc is in L; central/vc
    # is nmol/L = nM, matching Watson 2010 Table II residual-error units
    # (proportional 15.4%, additive 0.042 nM).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
