Momper_2016_fluconazole <- function() {
  description <- "One-compartment population PK model for fluconazole with first-order oral absorption and IV administration in extremely premature infants with birth weights < 750 g (Momper 2016)"
  reference <- "Momper JD, Capparelli EV, Wade KC, Kantak A, Dhanireddy R, Cummings JJ, Nedrelow JH, Hudak ML, Mundakel GT, Natarajan G, Gao J, Laughon M, Smith PB, Benjamin DK Jr; Fluconazole Prophylaxis Study Team. Population pharmacokinetics of fluconazole in premature infants with birth weights less than 750 grams. Antimicrob Agents Chemother. 2016;60(9):5539-5545. doi:10.1128/AAC.00963-16"
  vignette <- "Momper_2016_fluconazole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Allometrically scaled: CL ~ WT^0.75 and V ~ WT^1.0 with both exponents fixed per the Methods section. Paper Table 1 reports weight in grams (median 710 g, range 345-2,680 g); convert to kg before use.",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Multiplicative power covariate on CL with reference 0.8 mg/dL: (CREAT / 0.8)^-0.41. Paper source column is SCR. Subjects with SCR > 2 mg/dL were excluded from the original trial.",
      source_name        = "SCR"
    ),
    PAGE = list(
      description        = "Postmenstrual age = gestational age (weeks) / 4.35 + postnatal age (months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Multiplicative power covariate on CL with reference 28 weeks PMA (i.e., 28/4.35 = 6.4368 months): (PAGE / (28/4.35))^2.05. Momper 2016 expresses PMA in weeks with reference 28; the canonical PAGE column stores months, so the reference is rescaled to months here. Paper source column is PMA.",
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 141,
    n_samples      = 604,
    n_studies      = 1,
    age_range      = "PNA 3-47 days (median 23); GA 22.6-28.7 weeks (median 24.7); PMA 23.7-35.1 weeks (median 28.3)",
    weight_range   = "0.345-2.680 kg at first PK evaluation (median 0.710 kg)",
    sex_female_pct = 60,
    race_ethnicity = c(`Black or African American` = 53, White = 40, `American Indian or Alaska Native` = 5, Asian = 1),
    disease_state  = "Extremely premature infants with birth weights < 750 g receiving fluconazole prophylaxis for invasive candidiasis",
    dose_range     = "6 mg/kg IV (60-min infusion) or oral, twice weekly (Tuesdays and Fridays) for up to 42 days",
    regions        = "United States (multicenter randomized placebo-controlled trial)",
    notes          = "Table 1 baseline demographics. Exclusion criteria: AST/ALT > 250 U/L, SCR > 2 mg/dL, invasive candidiasis or congenital Candida infection at randomization, azole-antifungal hypersensitivity. Cesarean delivery in 67%; intubated in 81%. Final analysis dataset 604 plasma samples (61% from scavenged residual laboratory samples)."
  )

  ini({
    # Structural parameters: reference values per Table 3 (Momper 2016).
    # Time is in hours; allometric scaling references are absolute (WT in kg).
    lka     <- log(0.96);   label("Absorption rate constant (ka, 1/h)")                                   # Table 3 theta_KA point estimate
    lcl     <- log(0.0127); label("Allometric-scaled clearance (CL, L/h/kg^0.75)")                        # Table 3 theta_CL point estimate
    lvc     <- log(1.00);   label("Allometric-scaled volume of distribution (V, L/kg)")                   # Table 3 theta_V point estimate
    lfdepot <- log(1.00);   label("Oral bioavailability (F1, fraction)")                                  # Table 3 theta_F1 point estimate (estimated, not fixed)

    # Allometric exponents (held fixed per Methods: "Population PK parameters were
    # scaled by body size prior to evaluation of potential covariates.
    # Clearance was scaled by allometric weight (WT^0.75), and volume of
    # distribution was scaled by weight (WT^1.0).").
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")                                  # Methods, paragraph on covariate scaling
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on V  (unitless)")                                  # Methods, paragraph on covariate scaling

    # Covariate-effect power exponents on CL (multiplicative power form).
    # Final CL equation: CL = theta_CL * WT^0.75 * (SCR/0.8)^-0.41 * (PMA/28)^2.05.
    e_creat_cl <- -0.41; label("Serum creatinine power exponent on CL (unitless)")                         # Table 3 theta_SCR
    e_page_cl  <-  2.05; label("Postmenstrual age power exponent on CL (unitless)")                        # Table 3 theta_PMA

    # IIV (log-normal). NONMEM CV% from Table 3 mapped via omega^2 = log(CV^2 + 1):
    #   V  CV 13% -> omega^2 = log(0.13^2 + 1) = 0.01676
    #   CL CV 23% -> omega^2 = log(0.23^2 + 1) = 0.05153
    #   F1 CV 31% -> omega^2 = log(0.31^2 + 1) = 0.09175
    # ka has no BSV (Results: "Due to limited numbers of early samples after
    # oral administration in the data set, between-subject variability (BSV)
    # was not estimated for the absorption rate constant (ka).").
    etalcl     ~ 0.05153                                                                                    # Table 3 omega^2 CL (CV 23%)
    etalvc     ~ 0.01676                                                                                    # Table 3 omega^2 V  (CV 13%)
    etalfdepot ~ 0.09175                                                                                    # Table 3 omega^2 F1 (CV 31%)

    # Combined proportional + additive residual error (Methods:
    # "The initial model used a combined proportional error and additive
    # residual error value.").
    propSd <- 0.46; label("Proportional residual error (fraction)")                                         # Table 3 sigma proportional (46% CV)
    addSd  <- 505;  label("Additive residual error (ng/mL)")                                                # Table 3 sigma additive (505 ng/mL)
  })
  model({
    # Covariate ratios. CL covariate references: SCR 0.8 mg/dL, PMA 28 weeks
    # (28/4.35 months) per Table 2 / final CL equation.
    creat_ratio <- CREAT / 0.8
    page_ratio  <- PAGE / (28 / 4.35)

    # Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * WT^e_wt_cl * creat_ratio^e_creat_cl * page_ratio^e_page_cl
    vc <- exp(lvc + etalvc) * WT^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- exp(lfdepot + etalfdepot)

    # Concentration: dose in mg, vc in L gives mg/L; multiply by 1000 to get ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
