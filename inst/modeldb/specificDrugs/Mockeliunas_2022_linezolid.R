Mockeliunas_2022_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model for oral linezolid in adults with",
    "multidrug- and extensively drug-resistant tuberculosis (Mockeliunas",
    "2022), with transit absorption (the dose passes through five transit",
    "compartments at ktr = 6/MTT into an absorption compartment emptied at",
    "first-order ka) and concentration- and time-dependent auto-inhibition",
    "of elimination after Plock et al.: an empirical inhibition compartment",
    "equilibrates with the plasma concentration at rate kIC and scales the",
    "uninhibited apparent clearance by RCLF + (1 - RCLF) * IC50 / (IC50 +",
    "Ci), with kIC (0.0005 /h) and IC50 (0.38 mg/L) fixed to literature",
    "values from Keel et al. Body weight enters allometrically (exponents",
    "0.75 on CL/F and 1 on V/F, reference 70 kg). HIV co-infection raises",
    "CL/F, female sex raises ka, and concomitant P-glycoprotein inhibitors",
    "lengthen MTT. Inter-individual variability on CL/F and MTT;",
    "inter-occasion variability over seven sampling occasions on CL/F, V/F,",
    "ka and MTT; combined additive and proportional residual error."
  )
  reference <- paste(
    "Mockeliunas L, Keutzer L, Sturkenboom MGG, Bolhuis MS, Hulskotte LMG,",
    "Akkerman OW, Simonsson USH (2022).",
    "Model-Informed Precision Dosing of Linezolid in Patients with",
    "Drug-Resistant Tuberculosis.",
    "Pharmaceutics 14(4):753.",
    "doi:10.3390/pharmaceutics14040753.",
    sep = " "
  )
  vignette <- "Mockeliunas_2022_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte = "linezolid",
      units = "mg/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at admission",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed: the control stream (Text S1) uses BWAD, the body",
        "weight registered on the day of admission (Table 1 footnote a).",
        "Allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1),",
        "each scaled to 70 kg (Methods 2.2.3). Cohort mean 61.2 kg, range",
        "35.3-88.9 kg (Table 1)."
      ),
      source_name = "BWAD"
    ),
    HIV_POS = list(
      description = "HIV co-infection indicator (1 = HIV-positive, 0 = HIV-negative)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Time-fixed. Text S1: IF(HIV.EQ.1) CLHIV = (1 + THETA(10)), a",
        "multiplicative shift on CL/F; typical CL/F is 43% higher in",
        "HIV-positive patients (Discussion). Only 5 of 70 patients (7.1%,",
        "Table 1) were HIV-positive, and the 90% CI of the effect is wide",
        "(0.07-0.90); the authors ask that it be interpreted with caution."
      ),
      source_name = "HIV"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Time-fixed. Text S1: IF(SEX.EQ.1) KASEX = (1 + THETA(11)), and the",
        "Discussion states that typical ka was 95% higher in females, so",
        "the source SEX = 1 is female and maps to SEXF with no",
        "transformation. 38 of 70 patients were male (Table 1)."
      ),
      source_name = "SEX"
    ),
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor indicator (1 = on a P-gp inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Time-varying (Table S2 lists P-gp inhibitors among the",
        "time-varying categorical covariates). Text S1: IF(PGP_INH.EQ.1)",
        "MTTPGP_INH = (1 + THETA(12)), a multiplicative shift on MTT;",
        "typical MTT is 96% higher with a P-gp inhibitor (Discussion). The",
        "paper does not list which agents were classified as P-gp",
        "inhibitors or how many patients received one."
      ),
      source_name = "PGP_INH"
    ),
    OCC = list(
      description = "Integer-valued pharmacokinetic sampling-occasion index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Linezolid concentrations were collected at up to seven",
        "independent sampling occasions per patient (Methods 2.1), and Text",
        "S1 multiplexes the IOV etas with IF(OCC.EQ.1) ... IF(OCC.EQ.7)",
        "blocks. The index is decomposed inside model() into the",
        "mutually-exclusive indicators oc1..oc7; a record with OCC outside",
        "1..7 carries no IOV, matching the IOV_CL = 0 default of Text S1.",
        "For a single-occasion simulation pass OCC = 1."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL/F in the stepwise covariate search (Table S2) but not retained. Cohort mean 32 years, range 15-70 (Table 1)."
    ),
    ALCOHOL_ABUSE = list(
      description = "Alcohol abuse indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F, V/F, ka and MTT (Table S2) but not retained. Defined as more than 1 or 2 glasses of alcohol per day and fewer than 2 alcohol-free days per week before TB treatment (Text S1 $INPUT comment); 6 of 70 patients (Table 1)."
    ),
    DIS_DIAB = list(
      description = "Diabetes indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F, V/F, ka and MTT (Table S2) but not retained; 9 of 70 patients (Table 1)."
    ),
    SMOKE = list(
      description = "Smoking indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F, V/F, ka and MTT (Table S2) but not retained; 26 of 70 patients (Table 1)."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A4 inhibitor indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F as a time-varying covariate (Table S2) but not retained."
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F as a time-varying covariate (Table S2) but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_studies = 1L,
    age_range = "15-70 years",
    age_mean = "32 years",
    weight_range = "35.3-88.9 kg",
    weight_mean = "61.2 kg",
    height_mean = "1.70 m (range 1.50-1.93)",
    bmi_mean = "21.2 kg/m^2 (range 15.5-32.6)",
    sex_female_pct = 45.7,
    race_ethnicity = paste(
      "Origin of birth by WHO region (Table 1): African 10 (14.3%),",
      "Americas 2 (2.9%), South-East Asia 6 (8.6%), European 27 (38.6%),",
      "Eastern Mediterranean 15 (21.4%), Western Pacific 10 (14.3%)."
    ),
    disease_state = paste(
      "Multidrug- or extensively drug-resistant tuberculosis. HIV",
      "co-infection 5 (7.1%), diabetes 9 (12.9%), smoking 26 (37.1%),",
      "alcohol abuse 6 (8.6%), pregnancy 3 (4.3%)."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance mean 116.1 mL/min (range",
      "40.7-150.0), using lean body weight when BMI > 25 and truncated at",
      "150 mL/min (Table 1 footnote a)."
    ),
    dose_range = paste(
      "Oral linezolid 150-1200 mg per day, once or twice daily, for up to",
      "542 days in combination with other anti-TB drugs. Regimens by share",
      "of the 155 PK sampling occasions (Table S1): 150 mg QD 4.5%, 200 mg",
      "QD 3.2%, 200 mg BID 5.2%, 300 mg QD 21.3%, 300 mg BID 34.2%, 400 mg",
      "QD 0.6%, 600 mg QD 21.9%, 600 mg BID 9.0%."
    ),
    regions = "The Netherlands (Tuberculosis Center Beatrixoord, University Medical Center Groningen).",
    notes = paste(
      "Retrospective routine therapeutic-drug-monitoring data collected",
      "2007-2019: 811 total plasma linezolid concentrations from 70",
      "patients at up to seven sampling occasions each, mostly at steady",
      "state and usually including a pre-dose sample. LC-MS/MS assay with",
      "LLOQ 0.05 mg/L; the two BLQ observations were set to LLOQ/2. Fit in",
      "NONMEM 7.4.3 with FOCE-I. Besides the covariates in",
      "covariatesDataExcluded, the stepwise search also tested creatinine",
      "clearance (Cockcroft-Gault, mL/min, on CL/F), WHO region of birth,",
      "pre-emptive erythropoietin, and concomitant P-gp inducers (Table",
      "S2); none was retained. Pregnancy, antiretroviral therapy and",
      "therapeutic erythropoietin were not evaluated (3, 4 and 0 patients)."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters (Mockeliunas 2022 Table 2). CL/F and V/F are
    # the typical values for a 70-kg patient (Methods 2.2.3). Text S1 gives
    # the model structure; its $THETA values are initial estimates, not the
    # final estimates, which come from Table 2.
    # =========================================================================
    lcl <- log(6.3)
    label("Apparent uninhibited clearance at 70 kg (L/h)")
    # Table 2 row 'CL/F (L/h/70 kg) = 6.3 (90% CI 5.6-7.0; RSE 6.4%)'; THETA(1)

    lvc <- log(50.6)
    label("Apparent volume of distribution at 70 kg (L)")
    # Table 2 row 'Vd/F (L/70 kg) = 50.6 (48.5-53.1; 3.1%)'; THETA(2)

    lka <- log(1.8)
    label("First-order absorption rate constant from the absorption compartment, males (1/h)")
    # Table 2 row 'ka (h-1) = 1.8 (1.5-2.1; 13.8%)'; THETA(3). Male is the
    # reference category (Text S1 IF(SEX.EQ.0) KASEX = 1).

    lmtt <- log(0.53)
    label("Mean transit time through the absorption transit chain, no P-gp inhibitor (h)")
    # Table 2 row 'MTT (h) = 0.53 (0.44-0.61; 10.6%)'; THETA(6). Text S1:
    # NN = 5 transit compartments hard-coded, KTR = (NN + 1) / MTT.

    lke0 <- fixed(log(0.0005))
    label("Rate constant into the inhibition compartment (1/h)")
    # Table 2 row 'kIC (h-1) = 0.0005 FIX' with footnote e (value from Keel
    # et al., paper reference 25); Results 3.1 'kIC was fixed to the best
    # fitting literature value of 0.0005 h-1'. THETA(7) '0.0005 FIX' in Text S1.

    lic50 <- fixed(log(0.38))
    label("Inhibition-compartment concentration giving half of the maximum clearance inhibition (mg/L)")
    # Table 2 row 'IC50 (mg/L) = 0.38 FIX' with footnote e (Keel et al.);
    # THETA(9) '0.38 FIX' in Text S1.

    lfcl_noinh <- log(0.798)
    label("Remaining fraction of clearance that cannot be inhibited (fraction)")
    # Table 2 row 'RCLF = 0.798 (0.69-0.92; 11.3%)'; THETA(8), bounded (0, 1)
    # in Text S1.

    # =========================================================================
    # Allometric exponents (Methods 2.2.3: 'fixed to 0.75 and 1 for CL/F and
    # V/F, respectively ... scaled to 70 kg'; Text S1 TVCL and TVV lines).
    # =========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")
    # Methods 2.2.3; Text S1 TVCL = THETA(1)*(BWAD/70)**0.75

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V/F (unitless)")
    # Methods 2.2.3; Text S1 TVV = THETA(2)*(BWAD/70)

    # =========================================================================
    # Covariate effects (Table 2 'Covariates'). Each is a fractional change
    # applied as (1 + theta) when the indicator is 1 (Text S1).
    # =========================================================================
    e_hiv_pos_cl <- 0.43
    label("Fractional change in CL/F with HIV co-infection (unitless)")
    # Table 2 row 'HIV co-infection on CL/F = 0.43 (0.07-0.90; 122.0%)'; THETA(10)

    e_sexf_ka <- 0.95
    label("Fractional change in ka for females (unitless)")
    # Table 2 row 'Sex on ka = 0.95 (0.78-1.10; 14.0%)'; THETA(11)

    e_conmed_pgp_inh_mtt <- 0.96
    label("Fractional change in MTT with a concomitant P-gp inhibitor (unitless)")
    # Table 2 row 'P-gp inhibitor on MTT = 0.96 (0.84-1.09; 9.0%)'; THETA(12)

    # =========================================================================
    # Inter-individual variability (Table 2). Footnote a: IIV is 'expressed
    # as the standard deviation'; Text S1 applies the etas exponentially, so
    # the log-scale variance is the printed SD squared.
    # =========================================================================
    etalcl ~ 0.0676
    # Table 2 row 'IIV CL/F = 0.26 (0.21-0.31; 13.0%)' as an SD; 0.26^2

    etalmtt ~ 0.3844
    # Table 2 row 'IIV MTT = 0.62 (0.40-0.80; 19.9%)' as an SD; 0.62^2

    # =========================================================================
    # Inter-occasion variability over seven sampling occasions (Table 2,
    # footnote b: 'expressed as the standard deviation'). Text S1 codes one
    # $OMEGA BLOCK(1) per parameter for occasion 1 followed by six BLOCK(1)
    # SAME copies, so occasions 2-7 share occasion 1's variance.
    # =========================================================================
    etaiov_cl_1 ~ 0.0729
    # Table 2 row 'IOV CL/F = 0.27 (0.23-0.30; 9.0%)' as an SD; 0.27^2
    etaiov_cl_2 ~ fixed(0.0729)
    etaiov_cl_3 ~ fixed(0.0729)
    etaiov_cl_4 ~ fixed(0.0729)
    etaiov_cl_5 ~ fixed(0.0729)
    etaiov_cl_6 ~ fixed(0.0729)
    etaiov_cl_7 ~ fixed(0.0729)

    etaiov_vc_1 ~ 0.0676
    # Table 2 row 'IOV V/F = 0.26 (0.23-0.30; 8.5%)' as an SD; 0.26^2
    etaiov_vc_2 ~ fixed(0.0676)
    etaiov_vc_3 ~ fixed(0.0676)
    etaiov_vc_4 ~ fixed(0.0676)
    etaiov_vc_5 ~ fixed(0.0676)
    etaiov_vc_6 ~ fixed(0.0676)
    etaiov_vc_7 ~ fixed(0.0676)

    etaiov_ka_1 ~ 0.8649
    # Table 2 row 'IOV ka = 0.93 (0.71-1.16; 15.2%)' as an SD; 0.93^2
    etaiov_ka_2 ~ fixed(0.8649)
    etaiov_ka_3 ~ fixed(0.8649)
    etaiov_ka_4 ~ fixed(0.8649)
    etaiov_ka_5 ~ fixed(0.8649)
    etaiov_ka_6 ~ fixed(0.8649)
    etaiov_ka_7 ~ fixed(0.8649)

    etaiov_mtt_1 ~ 0.4761
    # Table 2 row 'IOV MTT = 0.69 (0.53-0.85; 13.4%)' as an SD; 0.69^2
    etaiov_mtt_2 ~ fixed(0.4761)
    etaiov_mtt_3 ~ fixed(0.4761)
    etaiov_mtt_4 ~ fixed(0.4761)
    etaiov_mtt_5 ~ fixed(0.4761)
    etaiov_mtt_6 ~ fixed(0.4761)
    etaiov_mtt_7 ~ fixed(0.4761)

    # =========================================================================
    # Residual error: combined additive and proportional on the normal scale
    # (Results 3.1; Text S1 W = SQRT((THETA(4)*IPRED)**2 + THETA(5)**2) with
    # $SIGMA 1 held constant, so THETA(4) and THETA(5) are the SDs).
    # =========================================================================
    propSd <- 0.054
    label("Proportional residual error (fraction)")
    # Table 2 row 'Proportional error (%) = 0.054 (0.045-0.065; 12.5%)'; THETA(4)

    addSd <- 0.53
    label("Additive residual error (mg/L)")
    # Table 2 row 'Additive error (mg/L) = 0.53 (0.483-0.570; 7.0%)'; THETA(5)
  })

  model({
    # --- 1. Occasion indicators and inter-occasion variability ------------
    # Mutually-exclusive indicators over the seven occasions of Text S1; a
    # record with OCC outside 1..7 zeroes every term and carries no IOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5 + oc6 * etaiov_vc_6 +
      oc7 * etaiov_vc_7
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6 +
      oc7 * etaiov_ka_7
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6 +
      oc7 * etaiov_mtt_7

    # --- 2. Individual parameters -----------------------------------------
    # Text S1 $PK: TVCL = THETA(1) * (BWAD/70)^0.75 * CLHIV;
    # TVV = THETA(2) * (BWAD/70); TVKA = THETA(3) * KASEX;
    # TVMTT = THETA(6) * MTTPGP_INH. CL = TVCL * EXP(ETA(1) + IOV_CL),
    # V = TVV * EXP(IOV_V), KA = TVKA * EXP(IOV_KA),
    # MTT = TVMTT * EXP(ETA(2) + IOV_MTT).
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl *
      (1 + e_hiv_pos_cl * HIV_POS)
    vc <- exp(lvc + iov_vc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + iov_ka) * (1 + e_sexf_ka * SEXF)
    mtt <- exp(lmtt + etalmtt + iov_mtt) *
      (1 + e_conmed_pgp_inh_mtt * CONMED_PGP_INH)
    ke0 <- exp(lke0)
    ic50 <- exp(lic50)
    fcl_noinh <- exp(lfcl_noinh)

    # --- 3. Micro-constants -----------------------------------------------
    # Five transit compartments hard-coded (Results 3.1, Text S1 NN = 5);
    # KTR = (NN + 1) / MTT (Methods 2.2.1).
    ntr <- 5
    ktr <- (ntr + 1) / mtt
    kel <- cl / vc

    # --- 4. Auto-inhibition of clearance ----------------------------------
    # Equation (6) / Text S1 DADT(2): elimination is scaled by
    # RCLF + (1 - RCLF) * (1 - Ci / (IC50 + Ci)), written here in the
    # algebraically identical IC50 / (IC50 + Ci) form. The factor is 1 when
    # the inhibition compartment is empty (first dose) and falls toward
    # RCLF as Ci accumulates.
    Cc <- central / vc
    inh <- fcl_noinh + (1 - fcl_noinh) * ic50 / (ic50 + effect)

    # --- 5. ODE system ----------------------------------------------------
    # Text S1 $DES. The dose enters the first transit compartment (A(1),
    # DEFDOSE, here depot); A(3)-A(6) are transit1-transit4; A(7), the
    # absorption compartment emptied at ka, is transit5; A(8), the
    # inhibition compartment, is effect and holds Ci in mg/L. Equation (2)
    # of the paper prints the inflow term of the transit chain with a minus
    # sign; the control stream's DADT(3)-DADT(6) carry the physically
    # correct plus sign, which is used here.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ka * transit5
    d/dt(central) <- ka * transit5 - kel * central * inh
    # Equation (7) / Text S1 DADT(8) = KIC * (CP - A(8))
    d/dt(effect) <- ke0 * (Cc - effect)

    # --- 6. Observation and residual error --------------------------------
    # Dose in mg / volume in L gives mg/L, the unit of the total plasma
    # linezolid concentrations the model was fit to.
    Cc ~ add(addSd) + prop(propSd)
  })
}
