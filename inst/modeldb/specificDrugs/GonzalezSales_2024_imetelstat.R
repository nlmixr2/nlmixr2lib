GonzalezSales_2024_imetelstat <- function() {
  description <- "Three-compartment population PK model for imetelstat (GRN163L), a 13-mer N3'-P5' thio-phosphoramidate oligonucleotide telomerase inhibitor, fit to 4375 plasma concentrations from 424 adults with hematologic malignancies (lower-risk MDS, myelofibrosis, multiple myeloma, ET/PV, CLD) or solid tumors who received IV imetelstat 0.4-11.7 mg/kg weekly to every-4-weeks (Gonzalez-Sales 2024). Imetelstat is described by a two-compartment nonlinear disposition model with saturable binding/distribution to a peripheral binding (BIND) compartment (Snoeck 1999 / Peletier 2017 parameterisation): free drug binds reversibly to a target pool with capacity Bmax (Kon, Koff); bound drug is internalised to a deep peripheral tissue (Kint) and returns to central as free drug (Kback); free drug also undergoes linear elimination from central (CL). Theory-based allometric exponents for body weight (1 on Vc, 0.75 on CL, -0.25 on Kback) are fixed. Final covariates: sex, dose, time, and MF / MM malignancy on CL; sex and MM malignancy on Vc; MF malignancy and baseline spleen volume on Bmax. The time effect encodes a hyperbolic decay of baseline CL: CL(t) = CL * t50_cl_time / (t + t50_cl_time)."
  reference <- "Gonzalez-Sales M, Lennox AL, Huang F, Pamulapati C, Wan Y, Sun L, Berry T, Kelly Behrs M, Feller F, Morcos PN. (2024). Population pharmacokinetics of imetelstat, a first-in-class oligonucleotide telomerase inhibitor. CPT Pharmacometrics Syst Pharmacol 13(7):1264-1277. doi:10.1002/psp4.13160."
  vignette <- "GonzalezSales_2024_imetelstat"
  units <- list(
    time          = "hour",
    dosing        = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Theory-based allometric scaling on Vc (exponent 1), CL (exponent 0.75), and Kback (exponent -0.25) with reference body weight 70 kg (Gonzalez-Sales 2024 Methods, 'Base structural model', Equation 1). The mean baseline body weight in the analysis dataset was 77.2 kg (range 44.0-161 kg).",
      source_name        = "WEIGHT"
    ),
    SEXF = list(
      description        = "Sex indicator: 1 if female, 0 if male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Categorical effect on CL (multiplicative, exp scale) and on Vc (multiplicative, exp scale). The source column SEX uses the same 0=male / 1=female encoding as the canonical SEXF (Gonzalez-Sales 2024 supplement S1, $INPUT block: 'SEX ;Sex 0 = Male; 1 = Female; -99 = Missing'). 41.7% of the analysis dataset was female (Results, 247 of 424 patients were male).",
      source_name        = "SEX"
    ),
    DIS_MF = list(
      description        = "Myelofibrosis disease-state indicator: 1 if myelofibrosis (MF), 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-MF subject (the analysis-dataset complement is solid tumors, MDS, MM, ET/PV, or CLD).",
      notes              = "Derived from the source NONMEM column MTYPE (0=solid tumors, 1=MF, 2=MDS, 3=MM, 4=ET/PV, 5=CLD) as DIS_MF = as.integer(MTYPE == 1) (Gonzalez-Sales 2024 supplement S1, MTYPE definition). Categorical effect on CL and on Bmax (multiplicative, exp scale). For Bmax, an additional power-form spleen-volume covariate is applied only to MF subjects (gating via DIS_MF in the model).",
      source_name        = "MTYPE == 1"
    ),
    DIS_MM = list(
      description        = "Active multiple myeloma disease-state indicator: 1 if active (non-smoldering) multiple myeloma, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-MM subject (the analysis-dataset complement is solid tumors, MF, MDS, ET/PV, or CLD).",
      notes              = "Renamed from the working covariate name MM to the canonical disease-indicator DIS_MM per the 2026-06-19 canonical-register standardization audit (parallels DIS_MF). Derived from the source NONMEM column MTYPE (0=solid tumors, 1=MF, 2=MDS, 3=MM, 4=ET/PV, 5=CLD) as DIS_MM = as.integer(MTYPE == 3) (Gonzalez-Sales 2024 supplement S1, MTYPE definition). Categorical effect on Vc only (multiplicative, exp scale).",
      source_alias       = "MM (working column name prior to the 2026-06-19 canonical rename)",
      source_name        = "MTYPE == 3"
    ),
    SPLV = list(
      description        = "Baseline spleen volume (cm^3).",
      units              = "cm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on Bmax centred at the MF-cohort median 3010 cm^3 (Gonzalez-Sales 2024 Figure 4 caption, 'baseline spleen volume 3010 cm3'; supplement S1 SPLEFF_BMAX = (SPLV0/3010)^THETA(9)). The effect is gated to MF subjects only (STDY.EQ.6 in the supplement); for non-MF subjects the spleen-volume term collapses to 1 in the model. Baseline spleen volume was only available from patients with MF (Study MYF2001) per Results, 'Patient characteristics at baseline'.",
      source_name        = "SPLV0"
    ),
    DOSE = list(
      description        = "Dose level for the current administration occasion.",
      units              = "umol",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL centred at the reference dose 108 umol (Gonzalez-Sales 2024 supplement S1: DOSE_CL = (DOSE/108)^THETA(12)). 108 umol of imetelstat sodium corresponds to approximately 7.5 mg/kg in a 70 kg patient (525 mg / 4896 g/mol = 107.2 umol). Carried as a per-record covariate that takes the current dose level; constant across observation records between dose changes.",
      source_name        = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age (years).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL and Vc but not retained in the final model (Gonzalez-Sales 2024 Results, 'Covariate analysis': 'Effects of other baseline covariates including ... age ... did not have a statistically significant or meaningful impact on imetelstat PK')."
    ),
    RACE_BLACK = list(
      description = "Race indicator: 1 if Black/African-American, 0 otherwise.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained in the final model (Gonzalez-Sales 2024 Results, no race effect identified)."
    ),
    RACE_ASIAN = list(
      description = "Race indicator: 1 if Asian, 0 otherwise.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained in the final model."
    ),
    ADA_POS = list(
      description = "Antidrug-antibody status: 1 if ADA-positive, 0 if ADA-negative.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the forward addition (time-variant ADA on CL, dOFV = -48.6) but eliminated in backward elimination (dOFV = +7.7). ADA status data were available from two studies (MYF2001 and MDS3001), contributing 65.6% of the analysis dataset (Gonzalez-Sales 2024 Results, 'Patient characteristics'; 'Covariate analysis')."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault, mL/min).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL (forward addition dOFV = -11.4) but eliminated in backward elimination (dOFV = +8.7). Mild-to-moderate renal impairment had no effect on imetelstat PK in the final model."
    ),
    HEPF = list(
      description = "Hepatic-function class (NCI ODWG; 0 = normal, 1 = mild, 2 = moderate, 3 = severe).",
      units       = "(category)",
      type        = "categorical",
      notes       = "Screened but not retained in the final model. No effect of mild-to-moderate hepatic impairment on imetelstat PK was identified (Gonzalez-Sales 2024 Results, Discussion)."
    ),
    ALB = list(
      description = "Serum albumin.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained in the final model. Canonical units standardized to SI g/L per the 2026-06-19 canonical-register audit; the source paper reports serum albumin in g/dL (1 g/dL = 10 g/L). No inline conversion is needed because ALB is an excluded covariate and is not referenced in model()/ini()."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained in the final model. Renamed from the working covariate name BILI to the canonical TBILI and units standardized to SI umol/L per the 2026-06-19 canonical-register audit; the source paper reports total bilirubin in mg/dL (1 mg/dL = 17.1 umol/L). No inline conversion is needed because TBILI is an excluded covariate and is not referenced in model()/ini().",
      source_alias = "BILI (working column name prior to the 2026-06-19 canonical rename)"
    ),
    AST = list(
      description = "Aspartate aminotransferase (U/L).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained in the final model."
    ),
    ALT = list(
      description = "Alanine aminotransferase (U/L).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 424L,
    n_studies      = 7L,
    age_range      = "adults; full age range not enumerated in the main paper (no age effect on PK identified)",
    age_median     = "not reported in main-paper Table S1 quotation",
    weight_range   = "44.0-161 kg",
    weight_median  = "77.2 kg (mean)",
    sex_female_pct = 41.7,
    race_ethnicity = "not enumerated in main paper text (race / ethnicity available in supplement Table S1; no race effect on PK identified)",
    disease_state  = "Pooled adults with hematologic malignancies (myelofibrosis, lower-risk MDS, multiple myeloma, essential thrombocythemia / polycythemia vera, chronic lymphoproliferative disease) or solid tumors. Studies included: CP04-151 (chronic lymphoproliferative disease, phase I), CP05-101 (refractory solid tumors, phase I), CP14A004 (refractory multiple myeloma, phase I), CP14B013 (multiple myeloma + lenalidomide maintenance, phase II), CP14B015 (essential thrombocythemia / polycythemia vera, phase II), MYF2001 (intermediate-2 / high-risk myelofibrosis refractory to JAK inhibitors, phase II), MDS3001 (transfusion-dependent low- or intermediate-1-risk MDS refractory to erythropoiesis-stimulating agents, phase II/III).",
    dose_range     = "0.4-11.7 mg/kg of imetelstat sodium (MW 4896 g/mol) administered as a 2-hour or 6-hour IV infusion (per-study schedules ranged from weekly to once every 4 weeks); the registered MDS3001 dose was 7.5 mg/kg Q4W and the registered MYF2001 / MYF3001 myelofibrosis dose was 9.4 mg/kg Q3W.",
    regions        = "Multinational (specific regions not enumerated in the main paper text on disk).",
    notes          = "Baseline-covariate summary statistics quoted from Gonzalez-Sales 2024 main-paper Results (Table S1 referenced but not on disk). Mean baseline body weight 77.2 kg, range 44.0-161 kg; 247 of 424 patients (58.3%) were male, so 41.7% female. Baseline spleen volume was available only from MF patients (Study MYF2001); MF-cohort median spleen volume 3010 cm^3 (Figure 4 caption). ADA status was categorised as negative / positive / missing (ADA data available from MYF2001 and MDS3001, contributing 65.6% of the analysis dataset). 374 of 4375 plasma observations (8.55%) were below the limit of quantification and modelled by the M3 method (Beal 2001)."
  )

  ini({
    # Structural parameters -- typical values for a 70 kg reference patient with male sex,
    # solid-tumor reference disease, reference dose 108 umol, time = 0.
    # All values from Gonzalez-Sales 2024 Table 2 (final population PK parameter estimates).
    lcl     <- log(1.00);    label("Clearance from central compartment (CL, L/h per 70 kg)")                  # Table 2: CL = 1.00 L/h/70 kg (RSE 3.50%)
    lvc     <- log(4.08);    label("Central volume of distribution (Vc, L per 70 kg)")                       # Table 2: Vc = 4.08 L/70 kg (RSE 2.55%)
    lkback  <- log(0.0253);  label("Transfer rate constant from peripheral binding pool to central (Kback, 1/h per 70 kg)") # Table 2: Kback = 0.0253 1/h/70 kg (RSE 7.58%)
    lbmax   <- log(15.0);    label("Total concentration of target binding sites (Bmax, umol/L)")             # Table 2: Bmax = 15.0 umol/L (RSE 7.08%)
    lkint   <- log(0.103);   label("Internalisation rate constant from bound complex to peripheral pool (Kint, 1/h)") # Table 2: Kint = 0.103 (RSE 9.08%); supplement S1 units 'L/h/70 kg' interpreted as 1/h with implicit unit volume per the Snoeck saturable-distribution convention -- see vignette Assumptions and deviations
    lkon    <- log(0.159);   label("Binding rate constant (Kon, L/(umol*h))")                                # Table 2: Kon = 0.159 (RSE 8.52%); supplement S1 units 'L^2/(umol/L * h)' interpreted as L/(umol*h) per dimensional analysis -- see vignette Assumptions and deviations
    lkoff   <- log(0.609);   label("Dissociation rate constant from bound complex (Koff, 1/h)")              # Table 2: Koff = 0.609 (RSE 10.7%); supplement S1 units 'L/h' interpreted as 1/h per dimensional analysis

    # Allometric exponents -- theory-based, fixed (Gonzalez-Sales 2024 Methods, 'Base structural model',
    # Equation 1: 'theory-based allometric exponents for body weight of 1, 0.75, and -0.25 were
    # incorporated on Vc, CL, and Kback'). Reference weight 70 kg.
    e_wt_cl    <- fixed(0.75);  label("Allometric exponent of (WT/70) on CL (unitless, fixed)")              # Methods Eq 1; supplement S1 line 70 FSIZE_CL = (WEIGHT/70)**0.75
    e_wt_vc    <- fixed(1);     label("Allometric exponent of (WT/70) on Vc (unitless, fixed)")              # Methods Eq 1; supplement S1 line 72 FSIZE_VC = (WEIGHT/70)**1
    e_wt_kback <- fixed(-0.25); label("Allometric exponent of (WT/70) on Kback (unitless, fixed)")           # Methods Eq 1; supplement S1 line 74 FSIZE_KBACK = (WEIGHT/70)**(-0.25)

    # Covariate effects on CL (multiplicative; categorical effects via exp(coef * indicator),
    # continuous effects via power(coef) at the reference value).
    e_dose_cl  <- -0.401; label("Power exponent of (DOSE/108 umol) on CL (unitless)")                        # Table 2: 'Effect of dose on CL' = -0.401 (RSE 9.69%)
    e_dis_mf_cl <- 0.511; label("Exponential coefficient of DIS_MF on CL (unitless)")                        # Table 2: 'Effect of MF on CL' = 0.511 (RSE 10.9%); MF patients have exp(0.511) = 1.67-fold higher CL
    e_sexf_cl  <- -0.299; label("Exponential coefficient of SEXF on CL (unitless)")                          # Table 2: 'Effect of sex on CL' = -0.299 (RSE 17.1%); females have exp(-0.299) = 0.74-fold lower CL

    # Time-on-CL effect: hyperbolic decay of baseline CL with characteristic half-time t50_cl_time.
    # Supplement S1 line 107-108 / 122: TIME_CL = TIME/(TIME+T50); CL = TVCL * ... * (1 - TIME_CL)
    #                                         = TVCL * ... * T50 / (T50 + TIME).
    # At t = 0, CL = TVCL (baseline). At t = t50_cl_time, CL = TVCL/2 (half of baseline). As t -> Inf,
    # CL -> 0 in the limit; in practice the function decays slowly given t50_cl_time ~ 245 days.
    lt50_cl_time <- log(5880); label("Half-time of the time effect on CL (t50_cl_time, hours)")              # Table 2: 'Effect of time on CL' = 5880 h (RSE 6.37%)

    # Covariate effects on Vc (multiplicative; categorical effects via exp(coef * indicator)).
    e_dis_mm_vc <- -0.233; label("Exponential coefficient of DIS_MM on Vc (unitless)")                       # Table 2: 'Effect of MM malignancy on Vc' = -0.233 (RSE 28.3%); MM patients have exp(-0.233) = 0.79-fold Vc
    e_sexf_vc  <- -0.122; label("Exponential coefficient of SEXF on Vc (unitless)")                          # Table 2: 'Effect of sex on Vc' = -0.122 (RSE 27.3%); females have exp(-0.122) = 0.89-fold Vc

    # Covariate effects on Bmax (continuous power form for SPLV, categorical exp form for DIS_MF).
    # The SPLV power effect is gated to DIS_MF = 1 subjects only (supplement S1 line 76-80:
    # 'IF(STDY.EQ.6) ... SPLEFF_BMAX = (SPLV0/3010)**THETA(9); ELSE SPLEFF_BMAX = 1'). Reference
    # SPLV = 3010 cm^3 (median in MF cohort; Figure 4 caption).
    e_splv_bmax    <- 0.772; label("Power exponent of (SPLV/3010 cm^3) on Bmax for MF subjects (unitless)")  # Table 2: 'Effect of spleen volume on Bmax' = 0.772 (RSE 27.3%)
    e_dis_mf_bmax  <- 1.44;  label("Exponential coefficient of DIS_MF on Bmax (unitless)")                   # Table 2: 'Effect of MF malignancy on Bmax' = 1.44 (RSE 7.25%); MF patients have exp(1.44) = 4.22-fold higher Bmax

    # Inter-individual variability (Gonzalez-Sales 2024 Table 2).
    # CL and Vc share a 2x2 BLOCK (correlation r = 0.545; supplement S1 $OMEGA BLOCK(2)).
    # Published CV% on CL was 43.7% and on Vc was 25.7%; the supplement OMEGA values 0.189 and
    # 0.0655 are reported as omega^2 (variance on the log-eta scale). With sqrt(omega^2) ~ omega
    # interpreted as the CV approximation, sqrt(0.189) ~= 0.435 (43.5%) and sqrt(0.0655) ~= 0.256
    # (25.6%) -- consistent with the published CV%. Off-diagonal covariance 0.0602 implies
    # correlation cov / (sd1 * sd2) = 0.0602 / (sqrt(0.189) * sqrt(0.0655)) = 0.541 ~ 0.545.
    etalcl + etalvc ~ c(0.189,
                        0.0602, 0.0655)                              # Table 2: IIV CL 43.7% CV, IIV Vc 25.7% CV, r(CL,Vc) = 0.545
    etalbmax ~ 0.181                                                  # Table 2: IIV Bmax 45.8% CV -> omega^2 = log(1+0.458^2) = 0.190, or omega = 0.425 if reported as SD-on-log-scale; supplement S1 OMEGA = 0.181

    # Residual error -- log-additive in NONMEM (Eq. 3): Ln(Yobs) = Ln(Ypred) + epsilon,
    # W = THETA(8) * exp(ETA(8)). Linear-space equivalent is proportional with propSd = sqrt(THETA(8)^2)
    # at the typical-value level. Published 'Residual variability' CV% = 21.8%; THETA(8) = 0.224
    # in the supplement gives CV = sqrt(exp(0.224^2) - 1) = 22.7% (close to 21.8%). The published
    # 'IIV residual variability' (omega^2 = 0.302 on log(W)) inflates the marginal residual variance
    # but is not encoded here -- see vignette Assumptions and deviations.
    propSd <- 0.218; label("Proportional residual error (fraction; log-additive in source NONMEM)") # Table 2: 'Residual variability' CV% = 21.8%
  })
  model({
    # ---- Disease-state indicators derived from MTYPE (covariate column) ----
    # The source NONMEM column MTYPE codes 0=solid tumors, 1=MF, 2=MDS, 3=MM, 4=ET/PV, 5=CLD.
    # The covariate dataset is expected to carry canonical DIS_MF (1 = MF) and DIS_MM (1 = active MM)
    # binary indicators derived from MTYPE: DIS_MF = (MTYPE == 1) and DIS_MM = (MTYPE == 3). The
    # spleen-volume effect on Bmax is gated to MF subjects via DIS_MF; non-MF subjects collapse
    # the spleen term to 1 regardless of their SPLV value.

    # ---- Time-on-CL hyperbolic decay (Gonzalez-Sales 2024 supplement S1) ----
    # f_time_cl(t) = T50 / (T50 + t)  =  1 - t / (T50 + t).
    # Equivalent forms; the second matches the supplement's TIME_CL = TIME/(TIME+T50) variable.
    t50_cl_time  <- exp(lt50_cl_time)
    f_time_cl    <- t50_cl_time / (t + t50_cl_time)

    # ---- Spleen-volume effect on Bmax (gated to MF subjects) ----
    # For MF subjects (DIS_MF = 1): f_splv_bmax = (SPLV / 3010)^e_splv_bmax.
    # For non-MF subjects: f_splv_bmax = 1 (the gating multiplier in the exponent collapses to 0).
    # Implemented as a power expression on (SPLV/3010) raised to (DIS_MF * e_splv_bmax) so that
    # non-MF subjects (whose SPLV may be NA or zero in the dataset) do not affect the model.
    f_splv_bmax  <- (SPLV / 3010)^(DIS_MF * e_splv_bmax)

    # ---- Individual-level PK parameters ----
    # CL: log-typical * size scaling * dose effect (continuous, power) * MF effect (binary, exp)
    #     * sex effect (binary, exp) * time effect (hyperbolic).
    cl    <- exp(lcl    + etalcl)  * (WT / 70)^e_wt_cl                                      *
             (DOSE / 108)^e_dose_cl                                                        *
             exp(e_dis_mf_cl * DIS_MF)                                                     *
             exp(e_sexf_cl   * SEXF)                                                       *
             f_time_cl

    # Vc: log-typical * size scaling * DIS_MM effect (binary, exp) * sex effect (binary, exp).
    vc    <- exp(lvc    + etalvc)  * (WT / 70)^e_wt_vc                                      *
             exp(e_dis_mm_vc * DIS_MM)                                                      *
             exp(e_sexf_vc * SEXF)

    # Kback: log-typical * size scaling (-0.25 exponent).
    kback <- exp(lkback)           * (WT / 70)^e_wt_kback

    # Bmax: log-typical * MF effect (binary, exp) * SPLV effect (MF-gated, power).
    bmax  <- exp(lbmax  + etalbmax) *
             exp(e_dis_mf_bmax * DIS_MF) *
             f_splv_bmax

    # Internalisation, binding, dissociation rate constants -- no IIV in the source model.
    kint  <- exp(lkint)
    kon   <- exp(lkon)
    koff  <- exp(lkoff)

    # ---- Compartment-level concentrations driving the binding kinetics ----
    # central holds drug amount in umol (the dosing unit); Cc = central / vc has units umol/L = uM.
    # peripheral1 (deep tissue) and complex (drug bound to target sites) follow the Snoeck 1999 /
    # Peletier 2017 saturable-distribution convention in which both compartments are treated as
    # concentration-like states with an implicit unit volume (1 L); their numerical values equal
    # uM concentrations directly (BMAX in uM is subtracted from the complex state). See vignette
    # Assumptions and deviations for the unit-convention rationale.
    Cc          <- central / vc
    free_target <- bmax - complex

    # ---- ODE system (direct translation of Gonzalez-Sales 2024 supplement S1 $DES) ----
    # NONMEM source equations:
    #   DADT(1) = -CL/VC*A(1) - KON*C1*(BMAX-A(3)) + KOFF*A(3) + KBACK*A(2)
    #   DADT(2) =                                              - KBACK*A(2) + KINT*A(3)
    #   DADT(3) =               KON*C1*(BMAX-A(3)) - KOFF*A(3)              - KINT*A(3)
    # Mass conservation: linear elimination CL*Cc out of central is the only true sink;
    # binding / internalisation / return cycles drug between central <-> complex <-> peripheral1.
    d/dt(central)     <- -(cl / vc) * central - kon * Cc * free_target + koff * complex + kback * peripheral1
    d/dt(peripheral1) <-                                                                  -kback * peripheral1 + kint * complex
    d/dt(complex)     <-                         kon * Cc * free_target - koff * complex                       - kint * complex

    # ---- Observation and error ----
    # Proportional residual on the linear scale, which approximates the source NONMEM log-additive
    # form for small CV (here ~22%).
    Cc ~ prop(propSd)
  })
}
