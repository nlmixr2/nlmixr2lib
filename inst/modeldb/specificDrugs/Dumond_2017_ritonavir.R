Dumond_2017_ritonavir <- function() {
  description <- "One-compartment population PK model with a first-order absorption lag time for low-dose (100 mg) boosting ritonavir in 31 HIV-infected adults aged 24-61 years, with total and unbound (protein-free) plasma concentrations described simultaneously. The ODE system is written on the unbound side: cl and vc are the apparent UNBOUND parameters CLu/F and Vu/F, the unbound concentration is Cu = central / vc, and the total concentration is Cc = Cu / fu. Ritonavir is the only one of the paper's three drugs to carry covariates: body weight scales CLu/F and Vu/F allometrically with exponents fixed at 0.75 and 1 about a 70 kg reference, and the fraction unbound is 1.52-fold higher (52% higher) in participants with a body mass index below 30 kg/m^2 than in the BMI at or above 30 kg/m^2 reference group. Between-occasion variability on CLu/F over two sampling occasions. Chronologic age, the Fried frailty phenotype and p16INK4a expression had no significant effect. The companion atazanavir model from the same paper is Dumond_2017_atazanavir."
  reference <- paste(
    "Dumond JB, Chen J, Cottrell M, Trezza CR, Prince HMA, Sykes C, Torrice C,",
    "White N, Malone S, Wang R, Patterson KB, Sharpless NE, Forrest A.",
    "Population pharmacokinetics modeling of unbound efavirenz, atazanavir,",
    "and ritonavir in HIV-infected subjects with aging biomarkers.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(2):128-135.",
    "doi:10.1002/psp4.12151.",
    sep = " "
  )
  vignette <- "Dumond_2017_unbound_antiretrovirals"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Drives allometric scaling of CLu/F",
        "(exponent fixed at 0.75) and Vu/F (exponent fixed at 1) about a",
        "70 kg reference, per Methods ('using a power model and the",
        "empirical scaling factor, 0.75 and 1, for CLu/F and Vu/F,",
        "respectively. Body weight was centered on 70 kg') and the",
        "control stream's LWT = LOG(WT/70) with MU_1 = THETA(1) +",
        "0.75*LWT and MU_2 = THETA(2) + LWT. Table 2 footnote *:",
        "'Population parameter estimates of subjects with body weight of",
        "70 kg.' The paper tabulates BMI rather than weight, so no",
        "weight range is published for this cohort.",
        sep = " "
      ),
      source_name = "WT"
    ),
    BMI = list(
      description = "Body mass index at baseline, used as a binary stratification at 30 kg/m^2",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Enters the model only through the binary",
        "indicator BMI < 30 kg/m^2 (the control stream's BMIC data item),",
        "which multiplies the fraction unbound by 1.52. Methods: 'Besides",
        "the original continuous value, BMI was also converted to 0 or 1",
        "split by 30 kg/m^2', and Results: 'Using BMI as a categorical",
        "variable gave a more stable model in terms of parameter estimates",
        "and was therefore retained.' The reference group is BMI at or",
        "above 30 kg/m^2, per Table 2 footnote e. Cohort median BMI 30.3",
        "kg/m^2 (range 20.2-40.4).",
        sep = " "
      ),
      source_name = "BMI"
    ),
    OCC = list(
      description = "Sampling-occasion index used for the between-occasion random effect on unbound clearance",
      units = "(index)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions, matching the supplementary control stream's",
        "$ABBR REPLACE ETA(OCC_CL)=ETA(6,7) multiplexer over the two",
        "$OMEGA BLOCK(1) / BLOCK(1) SAME interoccasion etas on CLu/F.",
        "Sparse-sampling participants contributed one to three occasions in",
        "the study; the final control stream allocates two IOV etas. For a",
        "single-occasion simulation set OCC = 1 throughout.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ritonavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "24-61 years",
    age_median = "49 years",
    sex_female_pct = 39,
    race_ethnicity = c(`African American` = 61, White = 32, Other = 6),
    disease_state = "HIV-1 infection on stable antiretroviral therapy for at least 2 weeks; adherence at least 27 of the previous 30 doses. Median HIV duration 10 years (range 1-24); median CD4 count 692 cells/mm^3 (range 375-1,501). Participants with hemoglobin < 10 mg/dL, estimated creatinine clearance < 30 mL/min (Cockcroft-Gault, total body weight) or DAIDS grade 2 or higher laboratory abnormalities were excluded, except that drug-attributed total bilirubin elevations were permitted.",
    dose_range = "Ritonavir 100 mg with atazanavir 300 mg by mouth once daily at steady state, co-administered with tenofovir disoproxil fumarate 300 mg and emtricitabine 200 mg once daily.",
    regions = "United States (UNC HealthCare Infectious Diseases Clinic, Chapel Hill NC; Cone Health Regional Center for Infectious Diseases, Greensboro NC). ClinicalTrials.gov NCT01180075.",
    bmi_range = "20.2-40.4 kg/m^2 (median 30.3)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) median 100 mL/min (range 67-227)",
    aging_markers = "Fried frailty phenotype: 21 participants (68%) with no positive components, 9 (29%) prefrail with 1-2 components, 1 (3%) frail with 3 or more. Log2(p16INK4a) median 2.2 (range 0.16-3.9). All three frail participants had BMI < 30 kg/m^2, but no relationship between frailty and BMI was observed.",
    notes = "Table 1 (ATV/RTV arm, n = 31). 25 participants were sparsely sampled (4 samples per occasion: predose, 2 h, 4-6 h, 10-14 h postdose, over one to three occasions) and 6 provided intensive sampling (11 samples around a single observed dose). Only nonfrail participants aged 55 or older were enrolled in the intensive group. Total and unbound concentrations were measured at every sampling time (unbound by rapid equilibrium dialysis with LC-MS/MS). 9% of RTV samples were below the limit of quantitation and were imputed as half the lower limit; five RTV samples had insufficient volume for the unbound assay and were treated as missing. NONMEM 7.3, ADVAN13, SAEM followed by importance sampling. Because of sparse sampling and substantial absorption-phase variability, ka and the lag time were estimated with RSE > 30%."
  )

  ini({
    # ==================================================================
    # Structural parameters - Dumond 2017 Table 2, RTV column
    # ("Estimates (RSE%)"). The paper parameterises the model on the
    # UNBOUND side: the ODE central compartment is the unbound plasma
    # compartment, so cl and vc below are CLu/F and Vu/F. The
    # total-drug parameters quoted in Table 2 are the products with fu
    # (Table 2 footnotes a-b): CL/F = CLu/F * fu = 4.85 L/h and
    # V/F = Vu/F * fu = 52.9 L. The micro-constant kel = cl/vc is
    # identical either way, so this is a reparameterisation and not a
    # different model.
    # Supplementary control stream PSP4-6-128-s003.txt
    # ('RTV protein binding model with free as the core') confirms the
    # structure: COMP(DEPOT) / COMP(UNBOUND) only (no peripheral),
    # ALAG1 on the depot, CU = A(2)/V, CTOTAL = CU/FU.
    # The starred Table 2 values are reported at a body weight of 70 kg
    # (Table 2 footnote *), and the fu value is the BMI at or above
    # 30 kg/m^2 reference group (Table 2 footnote e).
    # ==================================================================
    lka <- log(1.59); label("Absorption rate constant ka (1/h)") # Table 2 RTV: ka = 1.59 1/h (RSE 32%); bootstrap 1.61 (1.19-2.47)
    ltlag <- log(0.604); label("Absorption lag time (h)") # Table 2 RTV: lag time = 0.604 h (RSE 51%); bootstrap 0.678 (0.323-0.873)
    lcl <- log(772); label("Apparent unbound oral clearance CLu/F at 70 kg (L/h)") # Table 2 RTV: CLu/F = 772 L/h (RSE 3%) at 70 kg body weight (footnote *); bootstrap 749 (566-962). Total CL/F = 772 * 0.00628 = 4.85 L/h (Table 2 row 'CL/F', footnote a)
    lvc <- log(8430); label("Apparent unbound central volume Vu/F at 70 kg (L)") # Table 2 RTV: Vu/F = 8,430 L (RSE 2%) at 70 kg body weight (footnote *); bootstrap 8,100 (5,380-10,900). Total V/F = 8,430 * 0.00628 = 52.9 L (Table 2 row 'V/F', footnote b)
    lfu <- log(0.00628); label("Fraction of ritonavir unbound in plasma fu in the BMI at or above 30 kg/m^2 reference group (unitless)") # Table 2 RTV: fu = 0.628% (RSE 1%); bootstrap 0.667% (0.551-0.782%). Footnote e: fu,pop is the estimate for BMI > 30 kg/m^2. Entered as the fraction 0.00628

    # ==================================================================
    # Allometric exponents. Methods: the weight effects on CLu/F and
    # Vu/F were investigated 'using a power model and the empirical
    # scaling factor, 0.75 and 1, for CLu/F and Vu/F, respectively',
    # i.e. the exponents were imposed rather than estimated. The
    # control stream hardcodes them: MU_1 = THETA(1) + 0.75*LWT and
    # MU_2 = THETA(2) + LWT with LWT = LOG(WT/70). Neither exponent
    # appears as a row in Table 2, consistent with both being fixed.
    # Results: 'Empirical inclusion of the effects of body weight on
    # apparent unbound clearance (CLu/F) and apparent volume of
    # distribution (Vu/F) was significant for RTV, but not for ATV or
    # EFV.'
    # ==================================================================
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CLu/F (unitless)") # Methods 'Population pharmacokinetic modeling'; control stream MU_1 = THETA(1) + 0.75*LOG(WT/70)
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vu/F (unitless)") # Methods 'Population pharmacokinetic modeling'; control stream MU_2 = THETA(2) + LOG(WT/70)

    # ==================================================================
    # Covariate effect of BMI on the fraction unbound. Table 2 row
    # 'Influence of BMI < 30 on fu' = 1.52 (RSE 1%), bootstrap 1.51
    # (1.22-2.00). Table 2 footnote e defines the form as
    # fu = fu,pop * COEF^BMICAT, where fu,pop is the estimate for
    # BMI > 30 kg/m^2 and BMICAT is 1 when BMI < 30 kg/m^2 and 0
    # otherwise. Discussion: 'the unbound fraction of RTV is 52% higher
    # in participants with BMI < 30 kg/m^2, or 34% lower with
    # BMI > 30 kg/m^2' (1/1.52 = 0.658, a 34% reduction). The control
    # stream writes the same effect on the log scale as
    # MU_4 = THETA(4) + BMIC*THETA(6), so the log-scale coefficient
    # stored here is log(1.52) = 0.419 (the control stream's initial
    # estimate for THETA(6) was 0.4).
    # ==================================================================
    e_bmi_fu <- log(1.52); label("Log fold-change in fu for body mass index below 30 kg/m^2 relative to 30 kg/m^2 or above (unitless)") # Table 2 RTV row 'Influence of BMI < 30 on fu' = 1.52 (RSE 1%); bootstrap 1.51 (1.22-2.00)

    # ==================================================================
    # Inter-individual variability - Dumond 2017 Table 2, first
    # 'IIV (CV%)' block of the RTV column. Methods: IIV was 'assumed to
    # be normally distributed and exponentially related to the
    # population parameters', i.e. P_i = TVP * exp(eta_i), so the
    # log-scale variance is omega^2 = log(1 + CV^2). Shrinkage as
    # printed in brackets by Table 2.
    #
    # DEVIATION: the supplementary control stream fits a full
    # $OMEGA BLOCK(5) across all five etas. Table 2 prints only the
    # diagonal CV% values, so the off-diagonal covariances are not
    # recoverable and the IIV is encoded here as diagonal. See the
    # vignette's Assumptions and deviations section.
    # ==================================================================
    etalcl ~ 0.125712 # Table 2 RTV IIV CLu/F = 36.6% CV (RSE 3%) [shrinkage 4%]; bootstrap 37.1 (20.4-54.4) -> log(1 + 0.366^2)
    etalvc ~ 0.45013 # Table 2 RTV IIV Vu/F = 75.4% CV (RSE 25%) [shrinkage 2%]; bootstrap 80.8 (50.8-108) -> log(1 + 0.754^2)
    etalka ~ 0.343306 # Table 2 RTV IIV ka = 64% CV (RSE 17%) [shrinkage 27%]; bootstrap 68.8 (26.7-137) -> log(1 + 0.64^2)
    etaltlag ~ 0.891998 # Table 2 RTV IIV lag time = 120% CV (RSE 10%) [shrinkage 24%]; bootstrap 128 (68.6-165) -> log(1 + 1.20^2)
    etalfu ~ 0.0663906 # Table 2 RTV IIV fu = 26.2% CV (RSE 17%) [shrinkage 14%]; bootstrap 23.3 (11.2-34.6) -> log(1 + 0.262^2)

    # ==================================================================
    # Between-occasion variability on CLu/F. Results: 'Interoccasional
    # variability was incorporated on CLu/F for the three drugs.'
    # Table 2 prints this as a SECOND block also headed 'IIV (CV%)'
    # containing only a CLu/F row; that heading is a typesetting error
    # in the published table - the block is the interoccasion variance,
    # matching the control stream's $OMEGA BLOCK(1) + BLOCK(1) SAME
    # pair multiplexed by $ABBR REPLACE ETA(OCC_CL)=ETA(6,7).
    # The SAME keyword means both occasions share one variance, so the
    # second occasion is fixed to the first estimate.
    # ==================================================================
    etaiov_cl_1 ~ 0.311905 # Table 2 RTV second 'IIV (CV%)' block, CLu/F = 60.5% CV (RSE 55%); bootstrap 58.0 (32.8-87.1) -> log(1 + 0.605^2)
    etaiov_cl_2 ~ fix(0.311905) # control stream $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance

    # ==================================================================
    # Residual error - Dumond 2017 Table 2 'Residual error (CV%)'.
    # Control stream $ERROR: W = IPRED and
    # Y = IPRED + W*EPS(1)*(1-TYPE) + W*EPS(2)*TYPE with TYPE = 0 for
    # DVID 1 (unbound) and TYPE = 1 for DVID 2 (total), i.e. a
    # proportional error with a separate magnitude per analyte stream.
    #
    # DEVIATION: the control stream's $SIGMA BLOCK(2) estimates the
    # correlation between the two residual streams (the paper's L2 data
    # item). Only the two diagonal CV% values are published, and
    # nlmixr2 has no cross-endpoint residual correlation, so the two
    # streams are encoded as independent here.
    # ==================================================================
    propSd <- 0.414; label("Proportional residual error on total plasma ritonavir Cc (fraction)") # Table 2 RTV residual error, Total = 41.4% CV (RSE 20%); bootstrap 45.5 (35.3-56.0)
    propSd_Cu <- 0.344; label("Proportional residual error on unbound plasma ritonavir Cu (fraction)") # Table 2 RTV residual error, Unbound = 34.4% CV (RSE 33%); bootstrap 36.7 (28.7-45.7)
  })

  model({
    # ---- Derived covariate terms ------------------------------------
    # Control stream: LWT = LOG(WT/70); BMIC is the BMI < 30 kg/m^2
    # indicator (1 below 30, 0 at or above 30).
    bmi_lt30 <- (BMI < 30)

    # ---- Between-occasion variability multiplexer -------------------
    # Control stream: ETA(OCC_CL) selects one of two interoccasion etas
    # by the OCC data item; the selected eta is added to the CLu/F
    # log-scale typical value alongside the subject-level eta.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # ---- Individual parameters (unbound-side parameterisation) ------
    # Control stream $PK:
    #   MU_1 = THETA(1) + 0.75*LWT ; CL = EXP(MU_1 + ETA(1) + ETA(OCC_CL))
    #   MU_2 = THETA(2) + LWT      ; V  = EXP(MU_2 + ETA(2))
    #   MU_3 = THETA(3)            ; KA = EXP(MU_3 + ETA(3))
    #   MU_4 = THETA(4) + BMIC*THETA(6) ; FU = EXP(MU_4 + ETA(4))
    #   MU_5 = THETA(5)            ; ALAG1 = EXP(MU_5 + ETA(5))
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    fu <- exp(lfu + e_bmi_fu * bmi_lt30 + etalfu)
    tlag <- exp(ltlag + etaltlag)

    # ---- Concentrations --------------------------------------------
    # Control stream $ERROR: CU = A(2)/V and CTOTAL = CU/FU. The
    # central state holds the amount of drug; dividing by the apparent
    # unbound volume gives the unbound concentration and dividing that
    # by fu gives the total concentration. Equivalently
    # Cc = central / (vc * fu) = central / (V/F).
    Cu <- central / vc
    Cc <- Cu / fu

    # ---- ODE system -------------------------------------------------
    # Control stream $DES:
    #   DADT(1) = -KA*A(1)
    #   DADT(2) =  KA*A(1) - CL/V*A(2)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cu

    # ---- Absorption lag ---------------------------------------------
    # Control stream ALAG1 targets compartment 1 = COMP(DEPOT).
    alag(depot) <- tlag

    # ---- Residual error ---------------------------------------------
    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
