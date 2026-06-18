Fanta_2007_ciclosporin <- function() {
  description <- "Three-compartment population PK model with first-order absorption for ciclosporin in paediatric renal transplant candidates (Fanta 2007)"
  reference <- "Fanta S, Jonsson S, Backman JT, Karlsson MO, Hoppu K. Developmental pharmacokinetics of ciclosporin -- a population pharmacokinetic study in paediatric renal transplant candidates. Br J Clin Pharmacol. 2007;64(6):772-784. doi:10.1111/j.1365-2125.2007.03003.x"
  paper_specific_residual_sds <- c("propSdPo", "propSdIv", "addSdIv")
  vignette <- "Fanta_2007_ciclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; allometric scaling with reference weight 13 kg (population median per Fanta 2007 Table 2 footnote). Fixed exponent 3/4 on clearance parameters (CL, Q3, Q4) and 1 on volume parameters (V2, V3, V4) per Fanta 2007 Results 'Covariate model'.",
      source_name        = "BW"
    ),
    TCHOL = list(
      description        = "Total plasma cholesterol (fasting)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect (1 - 0.0542 * (TCHOL - 5.4)) applied identically to CL, Q3, Q4, V2, V3, V4 per Fanta 2007 Table 2 footnote and Results 'Covariate model'. Reference 5.4 mmol/L is the population median fasting plasma cholesterol.",
      source_name        = "CHOL"
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect (1 - 0.00732 * (HCT - 31)) applied identically to CL, Q3, Q4, V2, V3, V4 per Fanta 2007 Table 2 footnote and Results 'Covariate model'. Reference 31% is the population median haematocrit.",
      source_name        = "haematocrit"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect (1 + 0.000214 * (CREAT - 524)) applied identically to CL, Q3, Q4, V2, V3, V4 per Fanta 2007 Table 2 footnote and Results 'Covariate model'. Reference 524 umol/L is the population median serum creatinine.",
      source_name        = "serum creatinine"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous (IV) administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral, microemulsion formulation)",
      notes              = "Per-occasion dosing-route indicator (1 = IV infusion, 0 = oral microemulsion). Fanta 2007 carried a 23% lower typical CL on oral occasions than on IV occasions to correct an assay artifact (the monoclonal-RIA cross-reacts with ciclosporin metabolites that are formed more abundantly after oral than after IV administration); the effect is encoded as exp(e_route_iv_cl * ROUTE_IV), so the lcl parameter represents the apparent PO clearance and ROUTE_IV = 1 multiplies it by exp(0.2614) = 1/0.77 to recover the IV typical CL of 6.1 L/h. Distinct from the rxode2 cmt event column (cmt = central for IV doses, cmt = depot for oral doses).",
      source_name        = "route"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 162L,
    n_observations   = 2437L,
    n_studies        = 1L,
    study            = "Single-centre pretransplantation pharmacokinetic study at Hospital for Children and Adolescents, University of Helsinki (1988-2005); 166 patients enrolled, 162 eligible for the population analysis.",
    age_range        = "0.36-17.5 years",
    age_median       = "3.8 years",
    weight_range     = "6.9-64 kg",
    weight_mean      = "22.2 +/- 15.7 kg (mean +/- SD; Fanta 2007 Table 1)",
    sex_female_pct   = 34.0,
    race_ethnicity   = c(White_Finnish = 161/162 * 100, EastAfrican = 1/162 * 100),
    disease_state    = "Paediatric renal transplant candidates with renal disease (congenital nephrosis of the Finnish type, urethral valve, polycystic renal disease, nephronophtisis, other) waiting for renal transplantation; all on continuous ambulatory or continuous cycling peritoneal dialysis.",
    dose_range       = "IV 3 mg/kg as a 4-hour infusion (Sandimmun, Novartis); oral 10 mg/kg single dose as the microemulsion formulation (Sandimmun Neoral, Novartis). One patient was accidentally given 1 mg/kg orally. IV given to all 162 patients; oral microemulsion data contributed by 89 patients (the 73 patients who received the conventional oral formulation are excluded from the population PK model since that formulation is no longer used).",
    regions          = "Finland (Helsinki)",
    sampling         = "Whole blood ciclosporin samples (1 mL EDTA) drawn at 0, 1, 2, 3, 4, 6, 9, 12, 16, 24 h after oral dose; and at 0 (pre-dose), 2 (mid-infusion), 4 (end of infusion), 1, 2, 3, 4, 6, 9, 12, 16, 24 h after end of infusion. Assayed by specific monoclonal radioimmunoassay (Sandoz Sandimmune Kit until May 1994; Incstar/DiaSorin CycloTrac thereafter; detection limit ~5 ug/L; within- and between-run CV < 7% for concentrations > 30 ug/L).",
    reference_subject = "Body weight 13 kg, total plasma cholesterol 5.4 mmol/L, haematocrit 31%, serum creatinine 524 umol/L (population medians per Fanta 2007 Table 2 footnote).",
    notes            = "Comedications: 8 patients on potentially interacting drugs (1 carbamazepine, 2 oxcarbazepine, 5 phenobarbital); removing them did not change the final parameter estimates and they were retained in the dataset. 23% of patients had 1-4 ciclosporin concentration data points missing. Patients with polycystic renal disease and severe hepatic involvement were excluded."
  )

  ini({
    # Structural parameters at the reference subject (BW = 13 kg, TCHOL = 5.4
    # mmol/L, HCT = 31%, CREAT = 524 umol/L) per Fanta 2007 Table 2 and the
    # equation in the Table 2 footnote / Figure 5 caption.
    lka     <- log(0.68);    label("Absorption rate constant ka (1/h)")                                                                                                                                # Fanta 2007 Table 2
    lcl     <- log(0.77 * 6.1); label("Apparent clearance after oral administration at reference covariates (L/h); IV typical CL of 6.1 L/h is recovered via exp(lcl + e_route_iv_cl) when ROUTE_IV=1") # Fanta 2007 Table 2 (IV CL = 6.1) and Results 'Structural and stochastic models' paragraph 4 (PO CL is 23% lower than IV CL)
    lvc     <- log(5.3);     label("Central volume of distribution V2 at reference covariates (L)")                                                                                                    # Fanta 2007 Table 2
    lq      <- log(1.5);     label("First inter-compartmental clearance Q3 at reference covariates (L/h)")                                                                                             # Fanta 2007 Table 2
    lvp     <- log(19.6);    label("First peripheral volume of distribution V3 at reference covariates (L)")                                                                                            # Fanta 2007 Table 2
    lq2     <- log(3.0);     label("Second inter-compartmental clearance Q4 at reference covariates (L/h)")                                                                                            # Fanta 2007 Table 2
    lvp2    <- log(4.4);     label("Second peripheral volume of distribution V4 at reference covariates (L)")                                                                                          # Fanta 2007 Table 2
    lfdepot <- log(0.36);    label("Oral bioavailability F (fraction)")                                                                                                                                # Fanta 2007 Table 2

    # Fixed allometric exponents (3/4 on clearances, 1 on volumes). Fanta 2007
    # Results 'Covariate model' paragraph 1: "Because the fixed, allometrically
    # sound exponents (3/4 for clearance parameters and 1 for volume parameters)
    # were very close to the BW exponents estimated using NONMEM (0.787 for
    # clearance parameters and 0.951 for volume parameters), fixed values from
    # the literature were chosen."
    e_wt_cl <- fixed(0.75); label("Allometric exponent on clearance parameters CL, Q3, Q4 (unitless; fixed at the theoretical 3/4)")  # Fanta 2007 Results, fixed
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on volume parameters V2, V3, V4 (unitless; fixed at the theoretical 1)")       # Fanta 2007 Results, fixed

    # Linear-deviation covariate effects (same coefficient applied identically
    # to CL, Q3, Q4, V2, V3, V4 per Fanta 2007 Results 'Covariate model'
    # paragraph 2: "the same covariate coefficients were assigned for all
    # clearance and volume parameters").
    e_tchol_cl <- -0.0542;  label("Linear cholesterol effect on CL/V parameters (1/(mmol/L))")    # Fanta 2007 Table 2 footnote / Figure 5 caption
    e_hct_cl   <- -0.00732; label("Linear haematocrit effect on CL/V parameters (1/%)")           # Fanta 2007 Table 2 footnote / Figure 5 caption
    e_creat_cl <-  0.000214; label("Linear serum creatinine effect on CL/V parameters (1/(umol/L))") # Fanta 2007 Table 2 footnote / Figure 5 caption

    # Route-of-administration effect on CL: oral CL is 23% lower than IV CL
    # (Fanta 2007 Results 'Structural and stochastic models' paragraph 4). The
    # ROUTE_IV = 0 (oral) baseline encoded in lcl is multiplied by
    # exp(0.2614) = 1/0.77 = 1.299 when ROUTE_IV = 1 (IV) to give the Table 2
    # typical CL of 6.1 L/h.
    e_route_iv_cl <- 0.26136; label("Log-shift in CL for IV vs oral administration (unitless; IV typical CL = exp(lcl) * exp(e_route_iv_cl) = 6.1 L/h)") # = log(1/0.77); Fanta 2007 Results

    # Inter-individual variability. Fanta 2007 Table 2 reports CV% on the
    # log-normal scale; convert to variance via omega^2 = log(CV^2 + 1).
    # V2, Q4, V4 share a single eta with different magnitudes per Fanta 2007
    # Results 'Structural and stochastic models' paragraph 1; encoded here as a
    # 3x3 fully positively-correlated block (the published prose says
    # "complete positive OR negative correlation" without specifying signs).
    # See vignette Assumptions and deviations.
    # Shared standardized random effect for the V2 / Q4 / V4 block.
    # Fanta 2007 Results 'Structural and stochastic models' paragraph 1
    # reports a perfectly-correlated 3-eta block for V2, Q4, V4 (the
    # published prose: "complete positive OR negative correlation"). A
    # full 3x3 covariance with r = +1 between every pair is rank-1 and
    # singular, so rxode2's Cholesky-based simulator cannot decompose
    # it. Encoded here as one shared standardized eta scaled by the
    # per-parameter sqrt(variance) in model() — mathematically identical
    # to the published rank-1 block but numerically well-conditioned.
    eta_v2q4v4 ~ 1.0  # standardized shared eta for the V2/Q4/V4 perfect-correlation block
    etalcl     ~ 0.028501  # Fanta 2007 Table 2: IIV CL  CV 17% -> log(1 + 0.17^2)
    etalq      ~ 0.091705  # Fanta 2007 Table 2: IIV Q3  CV 31% -> log(1 + 0.31^2)
    etalvp     ~ 0.162519  # Fanta 2007 Table 2: IIV V3  CV 42% -> log(1 + 0.42^2)
    etalka     ~ 0.103358  # Fanta 2007 Table 2: IIV Ka  CV 33% -> log(1 + 0.33^2)
    etalfdepot ~ 0.012027  # Fanta 2007 Table 2: IIV F   CV 11% -> log(1 + 0.11^2)

    # Residual error. Fanta 2007 Table 2: IV arm uses combined proportional
    # (CV 8.9%) + additive (SD 1.5 ug/L); oral arm uses proportional-only
    # (CV 20%). The route-conditional residual is assembled in model().
    propSdPo <- 0.20;  label("Proportional residual error, oral arm (fraction)") # Fanta 2007 Table 2
    propSdIv <- 0.089; label("Proportional residual error, IV arm (fraction)")   # Fanta 2007 Table 2
    addSdIv  <- 1.5;   label("Additive residual error, IV arm (ug/L)")           # Fanta 2007 Table 2
  })
  model({
    # Shared linear-deviation covariate factor applied identically to every
    # clearance and every volume per Fanta 2007 Table 2 footnote.
    cov_factor <- (1 + e_tchol_cl * (TCHOL - 5.4)) *
                  (1 + e_hct_cl   * (HCT   - 31))  *
                  (1 + e_creat_cl * (CREAT - 524))

    # Allometric weight scaling: 3/4 on clearances, 1 on volumes (reference
    # 13 kg per the Table 2 footnote).
    allom_cl <- (WT / 13)^e_wt_cl
    allom_v  <- (WT / 13)^e_wt_vc

    # Per-parameter scaling of the shared standardized V2/Q4/V4 eta.
    # Variances and CVs are unchanged vs. the published rank-1 block:
    # var(etalvc) = (sqrt(0.014297))^2 = 0.014297 -> CV 12% on V2
    # var(etalq2) = (sqrt(0.141562))^2 = 0.141562 -> CV 39% on Q4
    # var(etalvp2) = (sqrt(0.060625))^2 = 0.060625 -> CV 25% on V4
    etalvc  <- 0.119571 * eta_v2q4v4  # sqrt(0.014297)
    etalq2  <- 0.376247 * eta_v2q4v4  # sqrt(0.141562)
    etalvp2 <- 0.246222 * eta_v2q4v4  # sqrt(0.060625)

    # Individual parameters
    ka      <- exp(lka  + etalka)
    cl      <- exp(lcl  + etalcl)  * allom_cl * cov_factor * exp(e_route_iv_cl * ROUTE_IV)
    vc      <- exp(lvc  + etalvc)  * allom_v  * cov_factor
    q       <- exp(lq   + etalq)   * allom_cl * cov_factor
    vp      <- exp(lvp  + etalvp)  * allom_v  * cov_factor
    q2      <- exp(lq2  + etalq2)  * allom_cl * cov_factor
    vp2     <- exp(lvp2 + etalvp2) * allom_v  * cov_factor
    fdepot  <- exp(lfdepot + etalfdepot)

    # Micro-constants for the three-compartment system
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODEs
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Oral bioavailability applies only to depot doses; IV doses go directly to
    # central via cmt = central in the event table.
    f(depot) <- fdepot

    # Concentration: dose in mg, volume in L -> central / vc is mg/L; multiply
    # by 1000 to express in ug/L (the unit Fanta 2007 uses throughout).
    Cc <- central / vc * 1000

    # Route-conditional residual error per Fanta 2007 Table 2 (IV: prop + add;
    # oral: prop only). When ROUTE_IV = 0, the additive term collapses to 0.
    propSd <- propSdPo + (propSdIv - propSdPo) * ROUTE_IV
    addSd  <- addSdIv  * ROUTE_IV
    Cc ~ add(addSd) + prop(propSd)
  })
}
