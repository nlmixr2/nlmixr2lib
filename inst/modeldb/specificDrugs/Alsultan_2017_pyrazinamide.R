Alsultan_2017_pyrazinamide <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and first-order elimination for oral pyrazinamide in adults with drug-susceptible pulmonary tuberculosis (Alsultan 2017); body weight is an allometric covariate on CL/F and V/F (fixed exponents 0.75 and 1) and biological sex is an exponential covariate on V/F"
  reference <- "Alsultan A, Savic R, Dooley KE, Weiner M, Whitworth W, Mac Kenzie WR, Peloquin CA, Tuberculosis Trials Consortium. Population pharmacokinetics of pyrazinamide in patients with tuberculosis. Antimicrob Agents Chemother. 2017;61(6):e02625-16. doi:10.1128/AAC.02625-16"
  vignette <- "Alsultan_2017_pyrazinamide"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight, used in allometric scaling on both CL/F and V/F with fixed exponents 0.75 and 1 respectively, normalised to a 70 kg reference subject (Alsultan 2017 Table 3 footnote a).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- the source paper uses female as the reference category for V/F (V/F = 46.5 L for a 70 kg female).",
      notes              = "Alsultan 2017 encodes sex as a male-indicator (1 = male, 0 = female) with female as the reference category (Table 3 footnote a: 'female was considered the reference'). To store under the canonical SEXF (1 = female, 0 = male) while preserving Alsultan's female-reference V/F = 46.5 L, the effect is applied in model() as exp(e_sex_vc * (1 - SEXF)), so SEXF = 1 (female) yields factor 1 and SEXF = 0 (male) yields the paper's male-vs-female exp-coefficient (V/F = 46.5 * exp(0.148) ~= 53.9 L for a 70 kg male, matching the paper's quoted 54.2 L within rounding).",
      source_name        = "SEX"
    )
  )

  population <- list(
    n_subjects     = 72L,
    n_studies      = 2L,
    n_observations = 499L,
    age_range      = "19-76 years",
    age_mean       = "36.7 years",
    weight_range   = "40-101.9 kg",
    weight_mean    = "59.8 kg",
    sex_female_pct = 16.7,
    africa_site_pct = 51.4,
    serum_creatinine_range = "0.5-1.2 mg/dL",
    serum_creatinine_mean  = "0.75 mg/dL",
    disease_state  = "Adults with drug-susceptible, smear-positive pulmonary tuberculosis enrolled in the PK substudies of Tuberculosis Trials Consortium (TBTC) studies 27 and 28.",
    dose_range     = "Oral pyrazinamide given as weight-banded daily doses: 1,000 mg (40-55 kg), 1,500 mg (56-75 kg), 2,000 mg (76-90 kg). Mean dose 1,351 mg (1,000-2,000 mg) corresponding to 22.7 mg/kg on average. PK sampling after the fourth or fifth daily dose (steady state).",
    regions        = "Uganda, South Africa, and the United States.",
    notes          = "Phase 2 prospective, placebo-controlled, randomised clinical trials. PZA was given alongside rifampin, isoniazid (or moxifloxacin), and ethambutol (or moxifloxacin). Plasma sampling at predose and 1, 2, 6, 8, 12, and 24 h postdose. PZA quantified by validated GC-MS (LLOQ 0.5 ug/mL, range 0.5-100 ug/mL). Baseline demographics from Alsultan 2017 Table 1; final population PK parameter estimates from Table 3."
  )

  ini({
    # Structural parameters -- typical values reference a 70 kg female
    # subject. Alsultan 2017 fits PZA with a one-compartment open model
    # with first-order absorption and first-order elimination
    # (Population pharmacokinetics paragraph, page 3). Parameters were
    # log-normally distributed: pi = p * exp(eta_i), var(eta) = omega^2.
    lka <- log(3.63);  label("Absorption rate constant (1/h)")                                # Alsultan 2017 Table 3 final-model ka
    lcl <- log(5.06);  label("Apparent oral clearance for a 70 kg subject (L/h)")             # Alsultan 2017 Table 3 final-model CL/F (70 kg)
    lvc <- log(46.5);  label("Apparent volume of distribution for a 70 kg female (L)")        # Alsultan 2017 Table 3 final-model V/F (70 kg female reference)

    # Covariate effects.
    # Allometric WT exponents are FIXED in the source paper: "Allometric
    # scaling was used to describe the effect of body weight on both V/F
    # and CL/F, using fixed exponents of 1 and 0.75, respectively"
    # (Population pharmacokinetics paragraph, page 4).
    e_wt_cl  <- fixed(0.75);  label("Allometric exponent on CL/F (fixed, unitless)")          # Alsultan 2017 Methods/Results: fixed allometric exponent
    e_wt_vc  <- fixed(1.0);   label("Allometric exponent on V/F (fixed, unitless)")           # Alsultan 2017 Methods/Results: fixed allometric exponent
    # Sex on V/F: source uses sex(male)=1 with female as reference.
    # Stored exponential coefficient is +0.148; applied in model() as
    # exp(e_sex_vc * (1 - SEXF)) to preserve the female-reference V/F.
    e_sex_vc <- 0.148;        label("Exponential coefficient of male sex on V/F (unitless; applied as (1 - SEXF))")  # Alsultan 2017 Table 3 footnote a

    # IIV. Source reports CV% on the natural-parameter scale (log-normal
    # IIV in Monolix). omega^2 = log(1 + CV^2), with CV expressed as a
    # fraction.
    # ka:   CV = 220%  -> omega^2 = log(1 + 2.20^2)  = 1.7647
    # V/F:  CV = 10.9% -> omega^2 = log(1 + 0.109^2) = 0.01181
    # CL/F: CV = 23%   -> omega^2 = log(1 + 0.23^2)  = 0.05155
    etalka ~ 1.7647   # Alsultan 2017 Table 3 final-model IIV ka  = 220% CV
    etalvc ~ 0.01181  # Alsultan 2017 Table 3 final-model IIV V/F = 10.9% CV
    etalcl ~ 0.05155  # Alsultan 2017 Table 3 final-model IIV CL/F = 23% CV

    # Combined residual error: y = f + (a + b*f) * eps in Monolix.
    # 'a' is the additive SD in observation units and 'b' is the
    # proportional SD as a fraction.
    addSd <- 0.94;  label("Additive residual SD (ug/mL)")                                     # Alsultan 2017 Table 3 final-model 'a' = 0.94
    propSd <- 0.10; label("Proportional residual SD (fraction)")                              # Alsultan 2017 Table 3 final-model 'b' = 10%
  })

  model({
    # Sex effect: source paper encodes sex(male) = 1 with female as the
    # reference category, so (1 - SEXF) reproduces the paper's male = 1
    # column while keeping SEXF (1 = female) as the canonical storage
    # convention.
    sex_male <- 1 - SEXF

    # Individual PK parameters. Allometric WT scaling on CL/F and V/F
    # with fixed exponents (0.75 and 1) and exponential male-sex effect
    # on V/F per Alsultan 2017 Table 3 footnote a:
    #   log(CL/F) = log(5.06) + 0.75 * (log(WT) - log(70))
    #   log(V/F)  = log(46.5) + 0.148 * sex(male) + (log(WT) - log(70))
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * exp(e_sex_vc * sex_male) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # ODEs: one-compartment open model with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation: oral CL/F and V/F absorb F into the apparent terms,
    # so dose/F is implicit in the typical-value parameterisation.
    # Concentration units: dose mg / volume L = mg/L = ug/mL, matching
    # the source paper.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
