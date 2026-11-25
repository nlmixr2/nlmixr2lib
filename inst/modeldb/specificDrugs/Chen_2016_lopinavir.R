Chen_2016_lopinavir <- function() {
  description <- "One-compartment population PK model with first-order absorption and an absorption lag time for oral lopinavir (co-administered with ritonavir) in 12 HIV-infected women studied at 20-24, 30 and 32 weeks of gestation and about 8 weeks postpartum, with total and unbound (protein-free) plasma concentrations described simultaneously. The ODE system is written on the unbound side: cl and vc are the apparent UNBOUND parameters CLu/F and Vu/F, the unbound concentration is Cu = central / vc, and the total concentration is Cc = Cu / fu with a linear (concentration-independent) unbound fraction. CLu/F rises exponentially with gestational week (exp(0.0148 * EGA); EGA = 0 postpartum), a 61% increase at 32 weeks; fu does not change with pregnancy. Interindividual variability on CLu/F and fu and interoccasion variability over four study visits on CLu/F, ka and fu. The companion ritonavir model from the same paper is Chen_2016_ritonavir."
  reference <- paste(
    "Chen J, Malone S, Prince HMA, Patterson KB, Dumond JB.",
    "Model-based analysis of unbound lopinavir pharmacokinetics in",
    "HIV-infected pregnant women supports standard dosing in the third",
    "trimester. CPT Pharmacometrics Syst Pharmacol. 2016;5(3):147-157.",
    "doi:10.1002/psp4.12065.",
    sep = " "
  )
  vignette <- "Chen_2016_lopinavir_ritonavir_pregnancy"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    EGA = list(
      description = "Maternal gestational age at the PK visit; 0 at the postpartum visit",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The control stream's WEEK data item. Methods ('Covariate",
        "evaluation'): 'Week of gestation was coded as zero for the",
        "postpartum visit', so EGA = 0 is the postpartum (non-pregnant)",
        "reference at which the typical CLu/F of 167 L/h applies.",
        "Study visits were at 20-24 weeks (visit 1), 30 weeks (visit 2),",
        "32 weeks (visit 3) and about 8 weeks postpartum (visit 4).",
        "Enters CLu/F as exp(e_ega_cl * EGA) (paper Eq 14; control",
        "stream CL = THETA(1)*EXP(ETA(1))*EXP(COEF*WEEK)*EXP(ETA(OCC_CL))).",
        "Time-varying across visits, constant within a visit.",
        sep = " "
      ),
      source_name = "WEEK"
    ),
    OCC = list(
      description = "Study-visit (occasion) index for the interoccasion random effects on CLu/F, ka and fu",
      units = "(index)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Four occasions = the four intensive PK visits (1: 20-24 weeks,",
        "2: 30 weeks, 3: 32 weeks, 4: about 8 weeks postpartum).",
        "Control stream $ABBR REPLACE ETA(OCC_CL)=ETA(3,4,5,6),",
        "ETA(OCC_KA)=ETA(7,8,9,10), ETA(OCC_FU)=ETA(11,12,13,14), each a",
        "$OMEGA BLOCK(1) followed by three BLOCK(1) SAME. Decomposed",
        "into oc1..oc4 inside model(); OCC = 0 (or any value outside",
        "1..4) gives the occasion-free typical value.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric weight scaling of CLu/F (exponent 0.75) and Vu/F (exponent 1) was tested for lopinavir and rejected by the likelihood ratio test (Results, 'LPV individual model'; Supplementary Material S2). Lean body mass was likewise not retained.",
      source_name = "WT"
    ),
    RACE_BLACK = list(
      description = "African-American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (African-American was the reference population)",
      notes = "Race was tested on CLu/F, Vu/F, the absorption parameters and fu (Methods 'Covariate evaluation') and was not retained (Results; Supplementary Material S2).",
      source_name = "RACE"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on fu and not significant (Results 'LPV individual model'; Supplementary Material S2). The paper reports albumin in g/dL (x 10 for the canonical g/L); it fell significantly during pregnancy (postpartum median 4.2 g/dL = 42 g/L).",
      source_name = "ALBUMIN"
    ),
    AAG = list(
      description = "Alpha-1 acid glycoprotein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on fu and not significant (Results 'LPV individual model'; Supplementary Material S2). The paper reports AAG in mg/dL (x 0.01 for the canonical g/L); it fell significantly during pregnancy (postpartum median 82 mg/dL = 0.82 g/L).",
      source_name = "AAG"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on fu and not significant (Results 'LPV individual model').",
      source_name = "BMI"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lopinavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lopinavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    n_observations = 330L,
    age_range = "18-35 years",
    age_median = "28 years",
    weight_range = "45.3-142.6 kg (at screening)",
    weight_median = "87.4 kg (at screening)",
    sex_female_pct = 100,
    race_ethnicity = c(`African American` = 75, Other = 25),
    disease_state = "HIV-1 infection in pregnancy; 67% antiretroviral-treatment-naive; CD4 > 400 cells/uL throughout; viral load suppressed (< 48 copies/mL) in all subjects in the latter part of the study.",
    dose_range = "Lopinavir/ritonavir tablets 400/100 mg twice daily from enrolment to 30 weeks, 500/125 mg twice daily from 30 weeks to 2 weeks postpartum, then 400/100 mg twice daily; co-administered with nucleoside agents.",
    ga_range = "20-24, 30 and 32 weeks of gestation, and about 8 weeks postpartum (EGA = 0)",
    regions = "United States (University of North Carolina at Chapel Hill; Northwestern University). ClinicalTrials.gov NCT00766818.",
    notes = "Secondary analysis of the Patterson 2012 study (12 women, four intensive PK visits each). 330 plasma samples, each assayed for total (HPLC/UV) and unbound (rapid equilibrium dialysis + LC/MS) lopinavir; no BLQ data. Observations with |CWRES| > 6 were excluded. Albumin 4.2 (4.1-4.4) g/dL and AAG 82 (72-116) mg/dL at 8 weeks postpartum (median, IQR), both significantly lower during pregnancy. NONMEM 7.3, ADVAN13, FOCE-I. Demographics from the Results first paragraph."
  )

  ini({
    # ==================================================================
    # Structural parameters - Chen 2016 Table 2, Lopinavir 'Estimate'.
    # The model is parameterised on the UNBOUND side (Methods 'LPV and
    # RTV structural model': 'Unbound and total concentrations were
    # comodeled with unbound drug as the central compartment'), so cl and
    # vc are CLu/F and Vu/F. Table 2 footnotes a-b give the total-drug
    # products CL/F = fu * CLu/F = 0.0258 * 167 = 4.31 L/h and
    # V/F = fu * Vu/F = 0.0258 * 3430 = 88.5 L; kel = cl/vc is the same
    # either way. CLu/F is the value at EGA = 0 (postpartum; paper Eq 14).
    # Supplementary Material S5 Part I (the authors' control stream)
    # confirms the structure: COMP(DEPOT), COMP(FREE), ALAG1 on the depot,
    # CU = A(2)/V, CTOTAL = CU/FU.
    # ==================================================================
    lka <- log(1.07); label("Absorption rate constant ka (1/h)") # Table 2 LPV ka = 1.07 1/h; bootstrap 0.906 (0.582, 1.78)
    ltlag <- log(1.77); label("Absorption lag time (h)") # Table 2 LPV 'Lag time' = 1.77 h; bootstrap 1.75 (1.25, 1.93)
    lcl <- log(167); label("Apparent unbound (intrinsic) clearance CLu/F at EGA = 0, postpartum (L/h)") # Table 2 LPV CLu/F = 167 L/h; bootstrap 167 (136, 202)
    lvc <- log(3430); label("Apparent unbound volume of distribution Vu/F (L)") # Table 2 LPV Vu/F = 3430 L; bootstrap 3325 (2655, 4256)
    lfu <- log(0.0258); label("Fraction of lopinavir unbound in plasma fu (unitless)") # Table 2 LPV fu = 2.58 percent; bootstrap 2.58 (2.48, 2.71); entered as the fraction 0.0258

    # Gestational-week effect on CLu/F, paper Eq 14:
    # CLu/F = theta * exp(coef * WEEK); control stream COEF = THETA(6).
    e_ega_cl <- 0.0148; label("Exponential coefficient of gestational week on CLu/F (1/week)") # Table 2 LPV 'Coefficient of gestational week impact on CLu/F' = 0.0148; bootstrap 0.0148 (0.00835, 0.0223)

    # ==================================================================
    # Interindividual variability - Table 2 LPV 'Interindividual
    # variability (%CV)'. The paper's %CV is 100 * sqrt(omega^2): the
    # ritonavir omegas the authors carried into the interaction control
    # streams (Supplementary Material S5 Parts III-V, $OMEGA BLOCK(2) FIX
    # 0.0073, -0.0199, 0.0597 and the IOV values 2.01, 0.0551, 0.131)
    # reproduce every ritonavir %CV and the -0.95 correlation in Table 2
    # exactly under that convention, so variance = (CV/100)^2 here.
    # Diagonal: the LPV control stream ($OMEGA 0.1 / 0.1) has no
    # off-diagonal term.
    # ==================================================================
    etalcl ~ 0.000729 # Table 2 LPV IIV CLu/F = 2.70 %CV -> 0.0270^2
    etalfu ~ 0.002116 # Table 2 LPV IIV fu = 4.60 %CV -> 0.0460^2

    # ==================================================================
    # Interoccasion variability over the four study visits - Table 2 LPV
    # 'Interoccasion variability (%CV)'. Control stream $OMEGA BLOCK(1)
    # followed by three BLOCK(1) SAME per parameter, so occasions 2-4
    # share the occasion-1 variance and are entered with fix().
    # ==================================================================
    etaiov_cl_1 ~ 0.061009 # Table 2 LPV IOV CLu/F = 24.7 %CV -> 0.247^2
    etaiov_cl_2 ~ fix(0.061009) # control stream $OMEGA BLOCK(1) SAME
    etaiov_cl_3 ~ fix(0.061009) # control stream $OMEGA BLOCK(1) SAME
    etaiov_cl_4 ~ fix(0.061009) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 1.9044 # Table 2 LPV IOV ka = 138 %CV -> 1.38^2
    etaiov_ka_2 ~ fix(1.9044) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fix(1.9044) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fix(1.9044) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fu_1 ~ 0.01 # Table 2 LPV IOV fu = 10.0 %CV -> 0.100^2
    etaiov_fu_2 ~ fix(0.01) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fu_3 ~ fix(0.01) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fu_4 ~ fix(0.01) # control stream $OMEGA BLOCK(1) SAME

    # ==================================================================
    # Residual error - Table 2 LPV 'Residual error (%CV)'. Control stream
    # $ERROR: W = IPRED; Y = IPRED + W*EPS(1)*(1-TYPE) + W*EPS(2)*TYPE,
    # TYPE = 0 for DVID 1 (unbound), 1 for DVID 2 (total): a proportional
    # error per analyte. The $SIGMA BLOCK(2) correlation (Table 2
    # 'sigma unbound - sigma total' r = 0.614) cannot be expressed across
    # two nlmixr2 endpoints and is dropped (see the vignette).
    # ==================================================================
    propSd <- 0.184; label("Proportional residual error on total plasma lopinavir Cc (fraction)") # Table 2 LPV residual error Total = 18.4 %CV; bootstrap 18.1 (14.0, 23.2)
    propSd_Cu <- 0.244; label("Proportional residual error on unbound plasma lopinavir Cu (fraction)") # Table 2 LPV residual error Unbound = 24.4 %CV; bootstrap 24.0 (21.1, 27.8)
  })

  model({
    # ---- Interoccasion multiplexer (control stream $ABBR REPLACE) ---
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_fu <- oc1 * etaiov_fu_1 + oc2 * etaiov_fu_2 + oc3 * etaiov_fu_3 + oc4 * etaiov_fu_4

    # ---- Individual parameters (control stream $PK) ------------------
    #   CL = THETA(1)*EXP(ETA(1))*EXP(COEF*WEEK)*EXP(ETA(OCC_CL))
    #   V  = THETA(2)
    #   KA = THETA(3)*EXP(ETA(OCC_KA))
    #   FU = THETA(4)*EXP(ETA(2))*EXP(ETA(OCC_FU))
    #   ALAG1 = THETA(5)
    cl <- exp(lcl + etalcl + iov_cl + e_ega_cl * EGA)
    vc <- exp(lvc)
    ka <- exp(lka + iov_ka)
    fu <- exp(lfu + etalfu + iov_fu)
    tlag <- exp(ltlag)
    kel <- cl / vc

    # ---- ODE system (control stream $DES) ----------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # ---- Concentrations (control stream $ERROR) ----------------------
    # CU = A(2)/V; CTOTAL = CU/FU. mg / L = mg/L (= ug/mL).
    Cu <- central / vc
    Cc <- Cu / fu

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
