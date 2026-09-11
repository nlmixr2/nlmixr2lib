Wei_2025_methotrexate <- function() {
  description <- "Three-compartment population PK model for high-dose intravenous methotrexate in Chinese adults with primary central nervous system lymphoma (Wei 2025), NONGENE variant. Clearance carries power effects of BSA-normalized eGFR, blood urea nitrogen and alanine aminotransferase; the inter-compartmental clearance to the first peripheral compartment carries a power effect of total serum protein. Inter-individual variability on clearance, both peripheral volumes, the central volume and the first inter-compartmental clearance, with proportional residual error. Parameter values are taken from the publication's Table 5 ('Final nongene-model' column) and the covariate equations 12, 13 and 16 to 19. A companion file Wei_2025_methotrexate_genotype.R holds the authors' second final model, which adds a composite ABCC4-ABCG2-ADORA2A genotype effect on clearance."
  reference <- paste(
    "Wei S, Zhang S, Wang D, Zhang D, Lu Q, Mo J, Yang Z, Guan L, He Y,",
    "Zhao Z, Mei S. (2025). Population pharmacokinetics of high-dose",
    "methotrexate in patients with primary central nervous system lymphoma.",
    "Front Pharmacol 16:1578033.",
    "doi:10.3389/fphar.2025.1578033.",
    sep = " "
  )
  vignette <- "Wei_2025_methotrexate"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the 2021 race-free CKD-EPI creatinine equation (Methods 'Study design', Equation 3, citing Inker 2021), NOT by Cockcroft-Gault -- the paper computes Cockcroft-Gault CLcr as well (Equation 2) but reports that 'eGFR emerged as a superior predictor of drug clearance than CLcr and Scr' (Discussion) and retains only eGFR. Normalized to 101.8 mL/min/1.73 m^2 in Equations 12 and 14, which matches the Table 3 cohort median of 101.8 exactly. Cohort range 5.4-162.9 mL/min/1.73 m^2 (Table 3). Time-varying: hepatic and renal function were 'routinely evaluated prior to MTX administration and monitored daily for at least three consecutive days thereafter', with the nearest value within 1-3 days substituted when a given day was missing and the record dropped when nothing fell inside a 7-day window (Methods 'Study design'). Enters clearance only, as a power term with the estimated exponent `e_crcl_cl` = 0.67. The temporal distribution of eGFR is plotted in Supplementary Appendix SA3.",
      source_name        = "eGFR"
    ),
    BUN = list(
      description        = "Blood urea nitrogen concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "SI units (mmol/L), NOT mg/dL -- Table 3 reports a median of 4.6 with range 0.5-19, which is the mmol/L scale (the mg/dL equivalent would be about 12.9). Normalized to 4.6 mmol/L in Equations 12 and 14, matching the Table 3 median exactly. Time-varying on the same daily-monitoring schedule as eGFR. Enters clearance only, as a power term with the estimated exponent `e_bun_cl` = -0.08, so a higher urea burden lowers methotrexate clearance (Discussion: 'higher BUN are significantly linked to decreased MTX clearance'). The paper argues this is an independent renal marker rather than collinearity with eGFR, reporting a between-covariate R^2 of only 0.22 and noting that BUN does not enter the 2021 CKD-EPI equation it used.",
      source_name        = "BUN"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalized to 25 U/L in Equations 12 and 14, matching the Table 3 cohort median of 25 exactly. Cohort range 2.2-1141.7 U/L (Table 3). Time-varying on the same daily-monitoring schedule as eGFR. Enters clearance only, as a power term with the estimated exponent `e_alt_cl` = +0.03. The sign is POSITIVE, i.e. liver injury RAISES methotrexate clearance in this cohort; the authors flag this as counterintuitive and offer a release-from-damaged-hepatocytes mechanism while conceding 'the exact mechanisms remain uncertain' (Discussion). The same Discussion paragraph supplies a check on this term: it states clearance rises 'approximately 5%-9%' at 5 to 20 times the upper limit of normal, and (5)^0.03 = 1.050 with (20)^0.03 = 1.094 reproduce those two figures when the multiplier of the upper limit of normal is applied to the 25 U/L centering value itself.",
      source_name        = "ALT"
    ),
    TPRO = list(
      description        = "Total serum protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the inter-compartmental clearance Q1 only (never clearance), as a power term with the estimated exponent `e_tpro_q` = -1.68 in this nongene model. NOTE a paper-internal mismatch, transcribed as printed: Equations 13 and 15 and the Abstract all normalize to 58 g/L, but the Table 3 cohort median is 61.8 g/L, and Methods 'Covariate model' states that 'all continuous covariates were standardized to their median values'. The printed equation constant 58 is used here, following the register's standing rule that the printed equation rather than the demographics table is the authority for a centering value; at the exponent -1.68 the difference moves Q1 by about 11%. Cohort range 27.4-95.7 g/L (Table 3). Time-varying on the same daily-monitoring schedule as eGFR. Methotrexate is about 50% protein bound, and the authors read the negative exponent as more bound (hence less transportable) drug at higher total protein (Discussion).",
      source_name        = "TP"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 752L,
    n_studies      = 1L,
    age_range      = "18.12-86.65 years (median 57.445)",
    age_median     = "57.445 years",
    weight_range   = "30-115 kg (median 68)",
    weight_median  = "68 kg",
    bsa_range      = "1.16-2.39 m^2 (median 1.73)",
    bsa_median     = "1.73 m^2",
    sex_female_pct = 44.1,
    renal_function = "eGFR 5.4-162.9 mL/min/1.73 m^2 (median 101.8) by the 2021 CKD-EPI equation; serum creatinine 24.8-641.7 umol/L (median 64.6); Cockcroft-Gault CLcr 5.9-361.8 mL/min (median 100.1). At least 17.4% of the cohort met the label definition of delayed elimination.",
    hepatic_function = "ALT 2.2-1141.7 U/L (median 25); AST 5-1915.2 U/L (median 20.4); total protein 27.4-95.7 g/L (median 61.8); albumin 19.9-51.8 g/L (median 37.5).",
    disease_state  = "Adults with primary central nervous system lymphoma (PCNSL) receiving high-dose methotrexate, most commonly combined with rituximab or cytarabine.",
    dose_range     = "Intravenous methotrexate 3.5 g/m^2, median infusion duration 3.1 h. Median of four infusions per patient (range 1-34). Given either as a single infusion or as a divided regimen in which 2 g was infused over 0.5 h and the remainder over the following 2.5 h. Inclusion required a dose of at least 0.5 g/m^2.",
    co_medication  = "Omeprazole in 41.03% and levetiracetam in 35.65% of concentration records; also ilaprazole (6.27%), furosemide (15.82%), torasemide (24.91%), bumetanide (16.10%) and NSAIDs (3.51%). None reached statistical significance as a covariate. Leucovorin rescue began 6 h post-infusion in every patient and is not represented in the model.",
    regions        = "China (single center: Beijing Tiantan Hospital, Capital Medical University), September 2016 through August 2023.",
    notes          = "Retrospective therapeutic-drug-monitoring cohort of 752 adults contributing 6074 methotrexate plasma concentrations. Each methotrexate administration was treated as an INDEPENDENT event in the dataset because dosing intervals exceeded five elimination half-lives (Methods 'Base model'), so the inter-individual variance terms are estimated across administrations rather than across patients. Concentrations were total (protein-bound plus free) drug by UHPLC-MS/MS with a lower limit of quantification of 0.002 umol/L; records below that limit were excluded. Estimation was by first-order conditional estimation extended least squares in Phoenix NLME 8.3. Model evaluation used 200 bootstrap replicates (not 1000, for run-time reasons) and a 1000-replicate visual predictive check. The authors caution that predictive accuracy beyond 120 h post-dose is limited because only 5% (304/6074) of samples fall there, and that the sparse distribution-phase sampling may bias the peripheral volumes and inter-compartmental clearances. Demographics from Table 3; parameter estimates from Table 5 ('Final nongene-model' column)."
  )

  ini({
    # Structural PK parameters -- Wei 2025 Table 5, 'Final nongene-model'
    # column, cross-checked against the closed forms printed as Equations 12,
    # 13 and 16 to 19. Each value is the typical value at the reference
    # covariates eGFR = 101.8 mL/min/1.73 m^2, BUN = 4.6 mmol/L, ALT = 25 U/L
    # and TP = 58 g/L. The paper's Q1/Vp1 pair maps onto the nlmixr2lib
    # canonical q/vp and its Q2/Vp2 pair onto q2/vp2; Equations 5 and 8 pair
    # Q1 with Vp1, and Equations 6 and 9 pair Q2 with Vp2.
    lcl  <- log(8.2)   ; label("Clearance CL at the reference covariates (L/h)")                                  # Table 5 nongene CL = 8.2 (%RSE 2.83) [7.75, 8.66]; also Equation 12
    lvc  <- log(33.39) ; label("Central volume of distribution Vc (L)")                                           # Table 5 nongene Vc = 33.39 (%RSE 3.51) [31.09, 35.69]; Equation 17 rounds this to 33.3
    lq   <- log(0.04)  ; label("Inter-compartmental clearance Q1 to peripheral1 at TP = 58 g/L (L/h)")            # Table 5 nongene Q1 = 0.04 (%RSE 8.13) [0.03, 0.05]; also Equation 13
    lvp  <- log(17.9)  ; label("First peripheral volume of distribution Vp1 (L)")                                 # Table 5 nongene Vp1 = 17.9 (%RSE 11.52) [13.86, 21.94]; also Equation 18
    lq2  <- log(0.09)  ; label("Inter-compartmental clearance Q2 to peripheral2 (L/h)")                           # Table 5 nongene Q2 = 0.09 (%RSE 5.18) [0.08, 0.10]; also Equation 16
    lvp2 <- log(1.14)  ; label("Second peripheral volume of distribution Vp2 (L)")                                # Table 5 nongene Vp2 = 1.14 (%RSE 4.40) [1.04, 1.23]; also Equation 19

    # Covariate effects. The Table 5 estimates are restated as the closed forms
    #   Equation 12: CL (L/h) = 8.2 * (eGFR/101.8)^0.67 * (BUN/4.6)^-0.08 * (ALT/25)^0.03
    #   Equation 13: Q1 (L/h) = 0.04 * (TP/58)^-1.68
    # Every centering constant equals the Table 3 cohort median EXCEPT the
    # 58 g/L for total protein, whose median is tabulated as 61.8; see the
    # `TPRO` covariateData note and the vignette Errata.
    e_crcl_cl <-  0.67 ; label("Power exponent on (CRCL / 101.8 mL/min/1.73 m^2) for CL (unitless)")              # Table 5 nongene theta_eGFR = 0.67 (%RSE 2.10) [0.65, 0.70]
    e_bun_cl  <- -0.08 ; label("Power exponent on (BUN / 4.6 mmol/L) for CL (unitless)")                          # Table 5 nongene theta_BUN = -0.08 (%RSE 10.23) [-0.09, -0.06]
    e_alt_cl  <-  0.03 ; label("Power exponent on (ALT / 25 U/L) for CL (unitless)")                              # Table 5 nongene theta_ALT = 0.03 (%RSE 14.47) [0.02, 0.03]
    e_tpro_q  <- -1.68 ; label("Power exponent on (TPRO / 58 g/L) for Q1 (unitless)")                             # Table 5 nongene theta_TP = -1.68 (%RSE 8.3) [-1.96, -1.41]

    # Inter-individual variability. Equation 10 is the exponential model
    # theta_i = theta_TV * exp(eta) with eta ~ N(0, omega^2) (Methods 'Base
    # model'), so each eta below is on the log scale. Table 5 heads these rows
    # 'IIV<param> (CV%)'; the tabulated numbers are read as 100 * omega, the
    # log-scale SD, so the variance is (CV%/100)^2. Two checks support the SD
    # reading over a variance reading: the reported %RSE on IIV_CL is 1.86%,
    # which is the order of the 1/sqrt(2N) precision an SD estimate carries on
    # the roughly 2000 administrations in the dataset, whereas a variance would
    # carry the much looser sqrt(2/N); and reading 27.3 as a variance would
    # imply an implausible log-scale SD of 5.2. The residual-error row in the
    # same table is likewise 100 * SD. See the vignette Errata for the one
    # residual ambiguity -- whether the authors tabulated 100*omega directly or
    # the exact log-normal CV sqrt(exp(omega^2)-1)*100 -- which is immaterial
    # for CL, Vc and Vp2 but moves omega_Q1 and omega_Vp1 by 16% and 12%.
    # There is NO inter-individual variability on Q2: Table 5 has no IIV_Q2
    # row, so q2 is a typical value only.
    etalcl  ~ 0.074529  # Table 5 nongene IIV_CL  = 27.3 CV%  (%RSE 1.86);  variance = 0.273^2
    etalvc  ~ 0.043890  # Table 5 nongene IIV_Vc  = 20.95 CV% (%RSE 5.77);  variance = 0.2095^2
    etalq   ~ 0.978319  # Table 5 nongene IIV_Q1  = 98.91 CV% (%RSE 9.43);  variance = 0.9891^2
    etalvp  ~ 0.613872  # Table 5 nongene IIV_Vp1 = 78.35 CV% (%RSE 25.38); variance = 0.7835^2
    etalvp2 ~ 0.072630  # Table 5 nongene IIV_Vp2 = 26.95 CV% (%RSE 5.22);  variance = 0.2695^2

    # Residual error. Equation 11 is Cobs = Cpred * (1 + epsilon) with
    # epsilon ~ N(0, sigma^2), i.e. proportional on the linear scale. Methods
    # 'Base model' records that proportional, exponential, additive and
    # combined additive-proportional models were all evaluated and 'the
    # proportional error model provided the best fit'. Table 5 reports the
    # same sigma for both final models.
    propSd <- 0.7379 ; label("Proportional residual error (fraction)")                                            # Table 5 sigma (proportional) = 73.79 (%RSE 1.07) [72.24, 75.35]
  })

  model({
    # 1. Individual PK parameters. Clearance carries the three power-scaled
    #    covariate terms of Equation 12; the first inter-compartmental
    #    clearance carries the total-protein term of Equation 13. Reference
    #    covariates: eGFR 101.8 mL/min/1.73 m^2, BUN 4.6 mmol/L, ALT 25 U/L,
    #    TP 58 g/L.
    cl  <- exp(lcl + etalcl) * (CRCL / 101.8)^e_crcl_cl *
      (BUN / 4.6)^e_bun_cl * (ALT / 25)^e_alt_cl
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + etalq) * (TPRO / 58)^e_tpro_q
    vp  <- exp(lvp + etalvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. Three-compartment intravenous disposition, Equations 4 to 9. Written
    #    there in concentration form as dA1/dt = -CL*Cc - Q1*(Cc - Cp1) -
    #    Q2*(Cc - Cp2) with Cc = A1/Vc, Cp1 = A2/Vp1 and Cp2 = A3/Vp2; the
    #    amount form below is the identical system. Methotrexate is given as
    #    an intravenous infusion, so the dose enters `central` directly and
    #    there is no absorption compartment.
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 4. Observation. Dose units umol and vc units L give Cc in umol/L, the
    #    unit the paper reports every concentration and every delayed-
    #    elimination threshold in. Methotrexate has a molar mass of
    #    454.44 g/mol, so 1 mg = 2.2005 umol.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
