Chen_2016_ritonavir <- function() {
  description <- "One-compartment population PK model with first-order absorption and an absorption lag time for oral low-dose (100-125 mg) boosting ritonavir co-administered with lopinavir in 12 HIV-infected women studied at 20-24, 30 and 32 weeks of gestation and about 8 weeks postpartum, with total and unbound (protein-free) plasma concentrations described simultaneously. The ODE system is written on the unbound side: cl and vc are the apparent UNBOUND parameters CLu/F and Vu/F, the unbound concentration is Cu = central / vc, and the total concentration is Cc = Cu / fu with a linear (concentration-independent) unbound fraction. CLu/F and Vu/F scale allometrically with body weight (exponents fixed at 0.75 and 1, 70 kg reference) and CLu/F rises exponentially with gestational week (exp(0.0184 * EGA); EGA = 0 postpartum); fu does not change with pregnancy. Correlated interindividual variability on CLu/F and fu (r = -0.95) and interoccasion variability over four study visits on ka, fu and relative bioavailability. The companion lopinavir model from the same paper is Chen_2016_lopinavir."
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
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CLu/F (exponent fixed at 0.75) and Vu/F",
        "(exponent fixed at 1) about 70 kg (Results 'RTV individual",
        "model': 'Both CLu/F and Vu/F were adjusted for body weight using",
        "empirical allometric rule'; control stream CL = THETA(1) *",
        "(WT/70)**0.75 ... and V = THETA(2)*(WT/70)). Table 2 labels the",
        "ritonavir CLu/F, CL/F, Vu/F and V/F estimates 'per 70kg'. The",
        "$INPUT also carries WTCHG, so weight was recorded per visit;",
        "median 87.4 kg (range 45.3-142.6) at screening.",
        sep = " "
      ),
      source_name = "WT"
    ),
    EGA = list(
      description = "Maternal gestational age at the PK visit; 0 at the postpartum visit",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The control stream's WEEK data item. Methods ('Covariate",
        "evaluation'): 'Week of gestation was coded as zero for the",
        "postpartum visit', so EGA = 0 is the postpartum (non-pregnant)",
        "reference. Study visits were at 20-24 weeks (visit 1), 30 weeks",
        "(visit 2), 32 weeks (visit 3) and about 8 weeks postpartum",
        "(visit 4). Enters CLu/F as exp(e_ega_cl * EGA) (control stream",
        "EXP(WEEK*THETA(6))). Time-varying across visits, constant within",
        "a visit.",
        sep = " "
      ),
      source_name = "WEEK"
    ),
    OCC = list(
      description = "Study-visit (occasion) index for the interoccasion random effects on ka, fu and relative bioavailability",
      units = "(index)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Four occasions = the four intensive PK visits (1: 20-24 weeks,",
        "2: 30 weeks, 3: 32 weeks, 4: about 8 weeks postpartum).",
        "Control stream $ABBR REPLACE ETA(OCC_KA)=ETA(3,4,5,6),",
        "ETA(OCC_FU)=ETA(7,8,9,10), ETA(OCC_BIO)=ETA(11,12,13,14), each",
        "a $OMEGA BLOCK(1) followed by three BLOCK(1) SAME. Decomposed",
        "into oc1..oc4 inside model(); OCC = 0 (or any value outside",
        "1..4) gives the occasion-free typical value.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
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
      notes = "Tested on fu and not significant (Results 'RTV individual model'; Supplementary Material S2). The paper reports albumin in g/dL (x 10 for the canonical g/L); postpartum median 4.2 g/dL = 42 g/L.",
      source_name = "ALBUMIN"
    ),
    AAG = list(
      description = "Alpha-1 acid glycoprotein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on fu and not significant (Results 'RTV individual model'; Supplementary Material S2). The paper reports AAG in mg/dL (x 0.01 for the canonical g/L); postpartum median 82 mg/dL = 0.82 g/L.",
      source_name = "AAG"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on fu and not significant (Results 'RTV individual model').",
      source_name = "BMI"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ritonavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE)
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
    notes = "Secondary analysis of the Patterson 2012 study (12 women, four intensive PK visits each). Total (HPLC/UV) and unbound (rapid equilibrium dialysis + LC/MS) ritonavir measured in each plasma sample; 330 plasma samples in total, no BLQ data. Observations with |CWRES| > 6 were excluded. NONMEM 7.3, ADVAN13, FOCE-I. Demographics from the Results first paragraph."
  )

  ini({
    # ==================================================================
    # Structural parameters - Chen 2016 Table 2, Ritonavir 'Estimate'.
    # Parameterised on the UNBOUND side (cl = CLu/F, vc = Vu/F); CLu/F
    # and Vu/F are 'per 70kg' and CLu/F is the value at EGA = 0
    # (postpartum). Table 2 footnotes a-b: CL/F = fu * CLu/F =
    # 0.149 * 160 = 23.8 L/h and V/F = fu * Vu/F = 0.149 * 1560 = 232 L.
    # The authors carried the same values, fixed, into their joint
    # interaction control streams (Supplementary Material S5 Parts III-V:
    # 160, 1560, 1.01, 0.149, 1.86 FIX and WEEK 0.0184 fix).
    # ==================================================================
    lka <- log(1.01); label("Absorption rate constant ka (1/h)") # Table 2 RTV ka = 1.01 1/h; bootstrap 1.01 (0.340, 2.00)
    ltlag <- log(1.86); label("Absorption lag time (h)") # Table 2 RTV 'Lag time' = 1.86 h; bootstrap 1.86 (1.67, 1.94)
    lcl <- log(160); label("Apparent unbound (intrinsic) clearance CLu/F at 70 kg and EGA = 0, postpartum (L/h)") # Table 2 RTV CLu/F = 160 L/h per 70 kg; bootstrap 157 (120, 202)
    lvc <- log(1560); label("Apparent unbound volume of distribution Vu/F at 70 kg (L)") # Table 2 RTV Vu/F = 1560 L per 70 kg; bootstrap 1630 (1040, 2100)
    lfu <- log(0.149); label("Fraction of ritonavir unbound in plasma fu (unitless)") # Table 2 RTV fu = 14.9 percent; bootstrap 14.8 (12.1, 18.1); S5 Part III THETA(9) 0.149

    # Allometric exponents fixed (Methods: 'empirical allometric scaling
    # factors, 0.75 and 1'; control stream (WT/70)**0.75 and (WT/70)).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CLu/F (unitless)") # Methods 'LPV and RTV structural model'; S5 Part II CL = THETA(1) * (WT/70)**0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vu/F (unitless)") # Methods 'LPV and RTV structural model'; S5 Part II V = THETA(2)*(WT/70)

    # Gestational-week effect on CLu/F: control stream EXP(WEEK*THETA(6)).
    e_ega_cl <- 0.0184; label("Exponential coefficient of gestational week on CLu/F (1/week)") # Table 2 RTV 'Coefficient of gestational week impact on CLu/F' = 0.0184; bootstrap 0.0188 (0.0129, 0.0262)

    # ==================================================================
    # Interindividual variability - correlated CLu/F and fu (control
    # stream $OMEGA BLOCK(2)). Entered as the final estimates the authors
    # fixed into the interaction control streams (S5 Parts III-V,
    # $OMEGA BLOCK (2) FIX 0.0073 / -0.0199 0.0597), which reproduce
    # Table 2 under the %CV = 100 * sqrt(omega^2) convention:
    # sqrt(0.0073) = 8.5 %CV (Table 2 'CLu/F' 8.5), sqrt(0.0597) =
    # 24.4 %CV (Table 2 'fu' 24.4), -0.0199 / sqrt(0.0073 * 0.0597) =
    # -0.95 (Table 2 'eta CL - eta fu' r = -0.95).
    # ==================================================================
    etalcl + etalfu ~ c(0.0073, -0.0199, 0.0597) # S5 Part III $OMEGA BLOCK(2) FIX; Table 2 RTV IIV CLu/F 8.5 %CV, fu 24.4 %CV, r = -0.95

    # ==================================================================
    # Interoccasion variability over the four study visits (control
    # stream $OMEGA BLOCK(1) + three BLOCK(1) SAME per parameter; the
    # SAME occasions are entered with fix()). Values from S5 Parts III-V
    # ($OMEGA BLOCK(1) 2.01E+00 FIX, 5.51E-02 FIX, 1.31E-01 FIX on
    # ETA(OCC_KAR), ETA(OCC_FUR), ETA(OCC_BIOR); the control-stream
    # comments label all three 'IOV KAR', but the $ABBR REPLACE order
    # assigns them to ka, fu and F), matching Table 2 RTV
    # 'Interoccasion variability (%CV)': ka 141.8, fu 23.5, F 36.2.
    # Bioavailability IOV: control stream F1 = EXP(ETA(OCC_BIO)), i.e.
    # a log-normal occasion effect about F = 1 on the depot.
    # ==================================================================
    etaiov_ka_1 ~ 2.01 # S5 $OMEGA 2.01E+00 FIX (ETA(OCC_KAR)); Table 2 RTV IOV ka = 141.8 %CV = 100 * sqrt(2.01)
    etaiov_ka_2 ~ fix(2.01) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fix(2.01) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fix(2.01) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fu_1 ~ 0.0551 # S5 $OMEGA 5.51E-02 FIX (ETA(OCC_FUR)); Table 2 RTV IOV fu = 23.5 %CV = 100 * sqrt(0.0551)
    etaiov_fu_2 ~ fix(0.0551) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fu_3 ~ fix(0.0551) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fu_4 ~ fix(0.0551) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_1 ~ 0.131 # S5 $OMEGA 1.31E-01 FIX (ETA(OCC_BIOR)); Table 2 RTV IOV F = 36.2 %CV = 100 * sqrt(0.131)
    etaiov_fdepot_2 ~ fix(0.131) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fix(0.131) # control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fix(0.131) # control stream $OMEGA BLOCK(1) SAME

    # ==================================================================
    # Residual error - Table 2 RTV 'Residual error (%CV)'; proportional
    # per analyte (control stream W = IPRED, EPS(1) unbound, EPS(2)
    # total). The $SIGMA BLOCK(2) correlation r = 0.611 is dropped (see
    # the vignette).
    # ==================================================================
    propSd <- 0.269; label("Proportional residual error on total plasma ritonavir Cc (fraction)") # Table 2 RTV residual error Total = 26.9 %CV; bootstrap 27.8 (23.1, 33.4)
    propSd_Cu <- 0.427; label("Proportional residual error on unbound plasma ritonavir Cu (fraction)") # Table 2 RTV residual error Unbound = 42.7 %CV; bootstrap 43.0 (38.2, 50.5)
  })

  model({
    # ---- Interoccasion multiplexer (control stream $ABBR REPLACE) ---
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_fu <- oc1 * etaiov_fu_1 + oc2 * etaiov_fu_2 + oc3 * etaiov_fu_3 + oc4 * etaiov_fu_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # ---- Individual parameters (control stream $PK) ------------------
    #   CL = THETA(1) * (WT/70)**0.75 * EXP(ETA(1)) * EXP(WEEK*THETA(6))
    #   V  = THETA(2)*(WT/70)
    #   KA = THETA(3)*EXP(ETA(OCC_KA))
    #   FU = THETA(4)*EXP(ETA(OCC_FU))*EXP(ETA(2))
    #   ALAG1 = THETA(5); F1 = EXP(ETA(OCC_BIO))
    cl <- exp(lcl + etalcl + e_ega_cl * EGA) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + iov_ka)
    fu <- exp(lfu + etalfu + iov_fu)
    tlag <- exp(ltlag)
    fdepot <- exp(iov_fdepot)
    kel <- cl / vc

    # ---- ODE system (control stream $DES) ----------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag
    f(depot) <- fdepot

    # ---- Concentrations (control stream $ERROR) ----------------------
    # CU = A(2)/V; CTOTAL = CU/FU. mg / L = mg/L (= ug/mL).
    Cu <- central / vc
    Cc <- Cu / fu

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
