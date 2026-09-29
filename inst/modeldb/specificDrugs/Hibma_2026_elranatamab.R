Hibma_2026_elranatamab <- function() {
  description <- "Two-compartment semi-mechanistic target-binding population PK model for the BCMA x CD3 bispecific antibody elranatamab in adults with relapsed or refractory multiple myeloma (Hibma 2026; MagnetisMM-1, -2, -3 and -9, N = 321). Elranatamab binds soluble BCMA (sBCMA) in the central compartment under rapid (quasi-)equilibrium with dissociation constant Kd; free drug, free sBCMA and the elranatamab-sBCMA complex each carry their own clearance and volume. States are total elranatamab (free + complex) in the central compartment, free elranatamab in the peripheral compartment, and total sBCMA with zero-order synthesis; the complex amount is the closed-form root of the binding quadratic. SC absorption is first order with bioavailability F; IV doses are 1-h infusions. Covariates: sex on elranatamab CL (linear), baseline body weight on elranatamab Vc and age on ka (power). Outputs: free (Cc) and total (Cc_total) elranatamab in ng/mL, free (Ctarget) and total (Ctotal_target) sBCMA in nM, each with proportional residual error. The companion landmark exposure-response model for cytokine release syndrome is Hibma_2026_elranatamab_crs."
  reference <- paste(
    "Hibma JE, Irby D, Liu A, Elmeliegy M, King LE, Gifondorwa D, Jiang S,",
    "Poels KE, Soltantabar P, Lon H-K, Shtylla B, Wang D, Williams JH,",
    "Nicholas T. Elranatamab population pharmacokinetics and",
    "exposure-response for cytokine release syndrome in patients with",
    "relapsed or refractory multiple myeloma.",
    "Clin Pharmacokinet. 2026;65:1173-1192.",
    "doi:10.1007/s40262-026-01663-z.",
    sep = " "
  )
  vignette <- "Hibma_2026_elranatamab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ng/mL (free and total elranatamab: Cc, Cc_total); nM (free and total soluble BCMA: Ctarget, Ctotal_target)"
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Hibma 2026 Section 3.3 final-model equation: TVCL = 0.324 L/d *",
        "(1 - 0.492 * CLcov,SEXF), 'where CLcov SEXF is 1 if applicable to",
        "each participant and 0 otherwise'. Already female-coded, so no",
        "transformation is needed. The sign is confirmed by the Results:",
        "females have the lower clearance and the higher exposure (median",
        "free AUCtau,ss 1.33 vs 0.579 in males, a ~49% difference in",
        "typical clearance). Analysis population 154/321 (48%) female",
        "(Table 1)."
      ),
      source_name = "SEXF"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value, time-fixed. Power effect on elranatamab central",
        "volume centred at 71.45 kg, the value printed in the Section 3.3",
        "equation (Table 1 median is 71.50 kg, range 36.5-159.6)."
      ),
      source_name = "BWT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value, time-fixed. Power effect on the SC absorption rate",
        "constant centred at 66 years (Section 3.3 equation; Table 1 median",
        "66, range 36-89)."
      ),
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator; 1 = Asian, 0 = other",
      units = "(binary)",
      type = "binary",
      notes = "Race was tested on CL and Vc in the stepwise covariate search (Section 2.5) and not retained. Analysis population Asian 49/321 (15%), Black 29 (9%), White 193 (60%), missing 50 (16%) (Table 1)."
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Tested on CL (Section 2.5) and not retained; no statistically significant effect of baseline eGFR or renal-impairment category (Discussion). Table 1 median 81.21, range 22.50-124.58."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Tested on CL as a hepatic-function marker (Section 2.5) and not retained. Table 1 median 3.60, range 1.90-4.80."
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested on CL as a hepatic-function marker (Section 2.5) and not retained. Table 1 median 0.40, range 0.10-3.80."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL as a hepatic-function marker (Section 2.5) and not retained. Table 1 median 23.00, range 8.00-247.00."
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status indicator; 1 = ADA-positive, 0 = ADA-negative",
      units = "(binary)",
      type = "binary",
      notes = "Baseline ADA, treatment-induced ADA and ADA over time were each tested on CL and none was retained (Discussion). Baseline ADA-positive 40/321 (12%), on-treatment ADA-positive 28/321 (9%) (Table 1)."
    ),
    SBCMA = list(
      description = "Baseline soluble BCMA concentration",
      units = "nM",
      type = "continuous",
      notes = "Tested on CL and Vc (Section 2.5) and not retained as a covariate. In this model baseline sBCMA is instead a structural parameter (rbase_target, 6.914 nM, with IIV) that sets the initial condition of the total_target state. Table 1 median 8.29 nM, range 0.00-266.67, 9 missing. 1 nM = 5.4 ng/mL (sBCMA 5.4 kDa, Section 2.2)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "elranatamab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "elranatamab (total: free + sBCMA-bound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "elranatamab (free)", units = "mg", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "soluble BCMA (total: free + elranatamab-bound)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 321L,
    n_studies = 4L,
    n_observations = "13,233 non-BLQ observations: 3739 total elranatamab, 2947 free elranatamab, 3812 total sBCMA, 2735 free sBCMA (data cutoff 2022)",
    age_range = "36-89 years (median 66)",
    weight_range = "36.5-159.6 kg (median 71.5)",
    sex_female_pct = 48,
    race_ethnicity = c(White = 60, Asian = 15, Black = 9, Missing = 16),
    disease_state = "Relapsed or refractory multiple myeloma",
    dose_range = "IV (1-h infusion) 0.1-50 ug/kg QW and SC 80-1000 ug/kg QW (MagnetisMM-1 Part 1); SC 600 then 1000 ug/kg QW or Q2W (Part 1.1, MagnetisMM-2); SC 44 then 76 mg QW (Part 2A); SC 12/32/76 mg on C1D1/C1D4/C1D8 then 76 mg QW, Q2W after 24 weeks in responders (MagnetisMM-3); SC 4/20/76 mg step-up (MagnetisMM-9)",
    regions = "Multinational (MagnetisMM-2 enrolled Japanese patients only)",
    baseline_sbcma = "median 8.29 nM, range 0.00-266.67 (Table 1)",
    notes = paste(
      "MagnetisMM-1 (NCT03269136) Parts 1, 1.1 and 2A (n = 53 + 19 + 15),",
      "MagnetisMM-2 (NCT04798586, n = 4), MagnetisMM-3 (NCT04649359)",
      "Cohorts A and B (n = 123 + 64) and MagnetisMM-9 (NCT05014412)",
      "Parts 1 and 2A (n = 33 + 10); Hibma 2026 Table 1 and ESM Table S1.",
      "Temporal validation used a later MagnetisMM-3 cutoff (March 2024,",
      "14,002 non-BLQ observations including Q4W dosing) without refitting."
    )
  )

  ini({
    # All values: Hibma 2026 Table 2 ('Final population PK model parameter
    # estimates'), final covariate model. The paper works on a molar scale:
    # amounts in nmol, concentrations in nM, Kd in nM (Section 2.2: MW 148 kDa
    # for elranatamab, 5.4 kDa for sBCMA). Time unit is day.

    # --- Free elranatamab disposition ---
    lcl <- log(0.324); label("Clearance of free elranatamab, typical male (L/day)") # Table 2 CLelranatamab 0.324 (RSE 9.114%); Section 3.3 equation
    lvc <- log(4.777); label("Central volume of free elranatamab at 71.45 kg (L)") # Table 2 Vc,elranatamab 4.777 (RSE 5.745%); Section 3.3 equation
    lvp <- log(2.83); label("Peripheral volume of free elranatamab (L)") # Table 2 Vp,elranatamab 2.83 (RSE 1.766%)
    lq <- log(0.225); label("Inter-compartmental clearance of free elranatamab (L/day)") # Table 2 Q 0.225 (RSE 1.933%)
    lka <- log(0.287); label("First-order SC absorption rate constant at age 66 years (1/day)") # Table 2 ka 0.287 (RSE 4.432%); Section 3.3 equation
    lfdepot <- log(0.562); label("Absolute SC bioavailability (fraction)") # Table 2 F 0.562 (RSE 0.564%)

    # --- Soluble BCMA and the elranatamab-sBCMA complex ---
    lcl_target <- log(0.273); label("Clearance of free soluble BCMA (L/day)") # Table 2 CLsBCMA 0.273 (RSE 21.03%)
    lvc_target <- log(15.418); label("Volume of distribution of soluble BCMA (L)") # Table 2 Vc,sBCMA 15.418 (RSE 10.904%)
    lrbase_target <- log(6.914); label("Baseline (drug-free) soluble BCMA concentration (nM)") # Table 2 BLsBCMA 6.914 (RSE 7.793%)
    lcl_complex <- log(0.164); label("Clearance of the elranatamab-sBCMA complex (L/day)") # Table 2 CLcomplex 0.164 (RSE 9.201%)
    lvc_complex <- log(3.802); label("Volume of distribution of the elranatamab-sBCMA complex (L)") # Table 2 Vc,complex 3.802 (RSE 4.848%)
    lkd <- log(3.138); label("Equilibrium dissociation constant of elranatamab:sBCMA (nM)") # Table 2 Kd 3.138 (RSE 3.312%)

    # --- Covariate effects (Section 3.3 equations; Table 2) ---
    e_sexf_cl <- -0.492; label("Fractional change in elranatamab CL for females, linear (unitless)") # Table 2 'Sex on CLelranatamab' -0.492 (RSE 12.533%)
    e_wt_vc <- 1.017; label("Power exponent of baseline body weight on elranatamab Vc (unitless)") # Table 2 'BWT on Vc,elranatamab' 1.017 (RSE 22.351%)
    e_age_ka <- -1.459; label("Power exponent of baseline age on ka (unitless)") # Table 2 'Age on ka' -1.459 (RSE 20.025%); printed CI '(2.031; -0.886)' drops the minus sign on the lower bound

    # --- Inter-individual variability ---
    # Table 2 reports IIV as 'CV (%)'. That CV is 100 * sqrt(omega^2) (the
    # first-order approximation), NOT the exact log-normal sqrt(exp(omega^2) - 1):
    # the printed 95% CIs are reproduced exactly by a symmetric Wald interval
    # on omega^2 (RSE on the variance scale) mapped back through sqrt(), e.g.
    # Vc: omega^2 = 0.6856^2 = 0.47005, SE = 0.13346 * 0.47005 = 0.06273,
    # CI (0.34709, 0.59301) -> sqrt -> (58.91%; 77.01%), exactly as printed;
    # CLsBCMA: 4.4807^2 = 20.0767, RSE 15.702% -> (372.8%; 512.4%), as printed.
    # Hence omega^2 = (CV / 100)^2 throughout. The four IIVs fixed 'to a
    # small positive value (~15% CV)' print as 14.83%, i.e. omega^2 = 0.022;
    # the Methods prose quotes the fixed variance as 0.025 (~15.8% CV), and
    # the table value is used here (see vignette Assumptions and deviations).
    # No off-diagonal elements are reported.
    etalcl ~ 1.0020 # Table 2 'IIV on CLelranatamab' 100.1% -> 1.001^2
    etalvc ~ 0.47005 # Table 2 'IIV on Vc,elranatamab' 68.56% -> 0.6856^2
    etalvp ~ fixed(0.022) # Table 2 'IIV on Vp elranatamab' 14.83%, held constant per the Table 2 footnote -> 0.1483^2
    etalq ~ fixed(0.022) # Table 2 'IIV on Q' 14.83%, held constant per the Table 2 footnote -> 0.1483^2
    etalka ~ 0.46799 # Table 2 'IIV on ka' 68.41% -> 0.6841^2
    etalfdepot ~ fixed(0.022) # Table 2 'IIV on F' 14.83%, held constant per the Table 2 footnote -> 0.1483^2
    etalcl_target ~ 20.0767 # Table 2 'IIV on CLsBCMA' 448.07% -> 4.4807^2
    etalvc_target ~ 1.84906 # Table 2 'IIV on Vc,sBCMA' 135.98% -> 1.3598^2
    etalrbase_target ~ 1.81199 # Table 2 'IIV on theta BLsBCMA' 134.61% -> 1.3461^2
    etalcl_complex ~ 0.62695 # Table 2 'IIV on CLcomplex' 79.18% -> 0.7918^2
    etalvc_complex ~ 0.49098 # Table 2 'IIV on Vc,complex' 70.07% -> 0.7007^2
    etalkd ~ fixed(0.022) # Table 2 'IIV on Kd' 14.83%, held constant per the Table 2 footnote -> 0.1483^2

    # --- Residual error: proportional on each analyte (Section 2.4) ---
    # Table 2 lists the four residual terms among the fixed effects with an
    # RSE, i.e. estimated as SD-scale THETAs (the usual SAEM Mu-referencing
    # set-up with SIGMA fixed to 1), and they are encoded as proportional SDs.
    propSd <- 0.347; label("Proportional residual error, free elranatamab (fraction)") # Table 2 'Residual error free elranatamab' 0.347 (RSE 2.61%)
    propSd_Cc_total <- 0.422; label("Proportional residual error, total elranatamab (fraction)") # Table 2 'Residual error total elranatamab' 0.422 (RSE 3.038%)
    propSd_Ctarget <- 0.518; label("Proportional residual error, free sBCMA (fraction)") # Table 2 'Residual error free sBCMA' 0.518 (RSE 3.503%)
    propSd_Ctotal_target <- 0.347; label("Proportional residual error, total sBCMA (fraction)") # Table 2 'Residual error total sBCMA' 0.347 (RSE 4.028%)
  })

  model({
    # Molecular weight of elranatamab (Section 2.2: 148 kDa). 1 mg = 1000/148
    # nmol and 1 nM = 148 ng/mL. Drug states are carried in mg so that SC and
    # IV doses are both given in mg; the paper's molar amounts A1, A2 and A3
    # are the mg states times mg_to_nmol.
    mw_elra <- 148
    mg_to_nmol <- 1000 / mw_elra

    # --- Individual parameters (Section 3.3 equations) ---
    cl <- exp(lcl + etalcl) * (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc + etalvc) * (WT / 71.45)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka) * (AGE / 66)^e_age_ka
    fdepot <- exp(lfdepot + etalfdepot)

    cl_target <- exp(lcl_target + etalcl_target)
    vc_target <- exp(lvc_target + etalvc_target)
    rbase_target <- exp(lrbase_target + etalrbase_target)
    cl_complex <- exp(lcl_complex + etalcl_complex)
    vc_complex <- exp(lvc_complex + etalvc_complex)
    kd <- exp(lkd + etalkd)

    k12 <- q / vc
    k21 <- q / vp

    # Zero-order sBCMA synthesis. Not printed in the paper: with BLsBCMA
    # estimated as a parameter and the system drug-free at baseline, the
    # steady state of dA4/dt = ksyn - CLsBCMA * Cf,sBCMA gives
    # ksyn = CLsBCMA * BLsBCMA (nmol/day), and A4(0) = BLsBCMA * Vc,sBCMA.
    ksyn <- cl_target * rbase_target

    # --- Rapid-binding equilibrium (Section 2.4 equations) ---
    # a1 = total elranatamab amount in the central compartment (A1, nmol);
    # total_target = total sBCMA amount (A4, nmol). Kd = Cf,elra * Cf,sBCMA /
    # Ccomplex with Cf,elra = (A1 - X)/Vc, Cf,sBCMA = (A4 - X)/Vc,sBCMA and
    # Ccomplex = X/Vc,complex gives the quadratic in the complex amount X
    # whose smaller root is the printed Xcomplex equation.
    a1 <- central * mg_to_nmol
    kdv <- kd * vc * vc_target / vc_complex
    sb <- kdv + a1 + total_target
    x_complex <- 0.5 * (sb - sqrt(sb * sb - 4 * a1 * total_target))

    cf_elra <- (a1 - x_complex) / vc # free elranatamab (nM)
    cf_sbcma <- (total_target - x_complex) / vc_target # free sBCMA (nM)
    c_complex <- x_complex / vc_complex # elranatamab-sBCMA complex (nM)

    # --- ODEs (Section 2.4 equations dA1/dt-dA4/dt) ---
    # The molar fluxes of the dA1/dt and dA2/dt equations are converted to
    # mg/day by dividing by mg_to_nmol. Only free drug distributes to the
    # peripheral compartment: k12 * (A1 - Xcomplex).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot +
      (-k12 * (a1 - x_complex) + k21 * peripheral1 * mg_to_nmol -
        cl * cf_elra - cl_complex * c_complex) / mg_to_nmol
    d/dt(peripheral1) <- (k12 * (a1 - x_complex)) / mg_to_nmol - k21 * peripheral1
    d/dt(total_target) <- ksyn - cl_target * cf_sbcma - cl_complex * c_complex

    f(depot) <- fdepot
    total_target(0) <- rbase_target * vc_target

    # --- Observations ---
    # Ct,elra = Cf,elra + Ccomplex and Ct,sBCMA = Cf,sBCMA + Ccomplex
    # (Section 2.4). Elranatamab is reported in ng/mL (nM * 148), sBCMA in nM.
    Cc <- cf_elra * mw_elra
    Cc_total <- (cf_elra + c_complex) * mw_elra
    Ctarget <- cf_sbcma
    Ctotal_target <- cf_sbcma + c_complex

    Cc ~ prop(propSd)
    Cc_total ~ prop(propSd_Cc_total)
    Ctarget ~ prop(propSd_Ctarget)
    Ctotal_target ~ prop(propSd_Ctotal_target)
  })
}
