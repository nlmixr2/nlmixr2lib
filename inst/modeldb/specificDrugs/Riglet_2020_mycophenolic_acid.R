Riglet_2020_mycophenolic_acid <- function() {
  description <- "Population PK model for plasma total, plasma unbound and intracellular (peripheral blood mononuclear cell, PBMC) mycophenolic acid (MPA) after oral mycophenolate mofetil in adult kidney transplant recipients of the CIMTRE study (Riglet 2020). Unbound MPA follows a two-compartment disposition with zero-order absorption into the central compartment and first-order unbound clearance; total plasma MPA is the unbound concentration times (1 + kns), a linear (non-saturable) protein-binding ratio; a third compartment holds PBMC MPA, exchanging with the unbound central concentration through an influx and an efflux clearance. Covariates: Cockcroft-Gault creatinine clearance on unbound clearance and peripheral volume, baseline serum albumin on the binding ratio, and ABCB1 3435C>T (rs1045642) TT homozygosity on the PBMC efflux clearance. Inter-individual and inter-occasion (visit) variability, proportional residual error on each of the three concentrations. All clearances and volumes are apparent (divided by the unestimated bioavailability F)."
  reference <- paste(
    "Riglet F, Bertrand J, Barrail-Tran A, Verstuyft C, Michelon H, Benech H,",
    "Durrbach A, Furlan V, Barau C. Population Pharmacokinetic Model of Plasma",
    "and Cellular Mycophenolic Acid in Kidney Transplant Patients from the CIMTRE",
    "Study. Drugs R D. 2020;20(4):331-342. doi:10.1007/s40268-020-00319-y"
  )
  vignette <- "Riglet_2020_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault formula (raw, NOT BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw Cockcroft-Gault creatinine clearance in mL/min (Riglet 2020 Methods section 2.1),",
        "NOT normalised to 1.73 m^2. Time-varying: collected at every visit and applied per",
        "occasion, because both parameters it acts on (CL/F and Vp/F) carry inter-occasion",
        "variability (Methods section 2.6). Power effects centred on the study median 54.81",
        "mL/min (Results section 3.1). Observed range 7.3-133.4 mL/min across the four PK visits",
        "(Table 2).",
        sep = " "
      ),
      source_name = "CrCL"
    ),
    ALB = list(
      description = "Serum albumin (human serum albumin, HSA) at baseline (day of transplant)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Riglet 2020 Methods section 2.6: when no inter-occasion variability could be estimated",
        "for the parameter carrying the covariate, the baseline (D0, day of transplant) value was",
        "used. The binding ratio kns has no IOV, so supply the per-subject D0 albumin as a",
        "time-fixed value. The power effect is centred on the 'study median', which the paper",
        "does not print; 37 g/L is back-solved from the two worked examples in Results section",
        "3.2.2 (fu = 1.3% at 45.8 g/L and 3.1% at 24.7 g/L are both reproduced to the printed",
        "digit only for a centre between 36.6 and 37.5 g/L). Post-transplant albumin by visit is",
        "in Table 2 (medians 30.3-36.4 g/L).",
        sep = " "
      ),
      source_name = "HSA"
    ),
    SNP_ABCB1_RS1045642_HOM = list(
      description = "ABCB1 3435C>T (rs1045642) homozygous-variant (TT) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC or CT genotype, i.e. carriers of the C3435 allele)",
      notes = paste(
        "1 = ABCB1 3435TT, 0 = CC or CT. Riglet 2020 retained the recessive model (Results",
        "section 3.2.2), in which heterozygotes pool with wild-type homozygotes. Genotype",
        "frequencies CC/CT/TT 30/28/11 among the 69 genotyped patients (Table 1); nine patients",
        "with a missing genotype were imputed to the most common genotype (CC).",
        sep = " "
      ),
      source_name = "ABCB1 3435C>T"
    ),
    OCC = list(
      description = "Occasion (visit) index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Occasion = study visit (Methods section 2.5: 'Interoccasion (or visit) variability').",
        "PK data were modelled at four visits: OCC = 1 at day 15 (D15), 2 at month 1 (M1), 3 at",
        "month 2 (M2) and 4 at month 6 (M6) after transplantation (Results section 3.1).",
        "Decomposed inside model() into the indicators occ1..occ4 that select the per-occasion",
        "IOV etas; a record with OCC outside 1-4 receives no IOV.",
        sep = " "
      ),
      source_name = "visit"
    )
  )

  compartmentData <- list(
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    pbmc = list(analyte = "mycophenolic acid", units = "mg", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    age_range = "21-78 years",
    age_median = "50 years",
    weight_range = "36-125 kg",
    weight_median = "66.5 kg",
    sex_female_pct = 42.3,
    disease_state = "Adult kidney transplant recipients followed for 6 months after transplantation (CIMTRE study)",
    dose_range = "Oral mycophenolate mofetil 1000 mg twice daily starting dose (adjusted for adverse effects), with tacrolimus and prednisone",
    regions = "France (Bicetre Hospital, Paris)",
    renal_function = "Cockcroft-Gault CrCL median 54.81 mL/min overall; 7.3-133.4 mL/min across visits (Table 2)",
    co_medication = "Tacrolimus (target trough 5-15 ng/mL) and prednisone in all patients",
    notes = paste(
      "82 enrolled 2005-2008, 78 analysed (45 men). 1931 concentrations: 925 total plasma, 560",
      "unbound plasma and 446 PBMC MPA, none below the limit of quantitation. PK data at D15,",
      "M1, M2 and M6 for 71, 73, 70 and 57 patients. Demographics Table 1; CrCL and albumin by",
      "visit Table 2.",
      sep = " "
    )
  )

  ini({
    # Structural parameters (Riglet 2020 Table 3, 'Covariate model' columns).
    # The reference subject has CrCL = 54.81 mL/min, baseline albumin 37 g/L
    # and an ABCB1 3435 CC or CT genotype.
    ld1 <- log(1.29); label("Zero-order absorption duration Tk0 (h)") # Table 3 Tk0 = 1.29 h (RSE 8%)
    lvc <- log(1620); label("Apparent central volume of unbound MPA Vcu/F (L)") # Table 3 Vcu/F = 1620 L (RSE 9%)
    lcl <- log(900); label("Apparent clearance of unbound MPA CLu/F (L/h)") # Table 3 CLu/F = 900 L/h (RSE 4%)
    lq <- log(2040); label("Apparent intercompartmental clearance of unbound MPA Qu/F (L/h)") # Table 3 Qu/F = 2040 L/h (RSE 15%)
    lvp <- log(19400); label("Apparent peripheral volume of unbound MPA Vpu/F (L)") # Table 3 Vpu/F = 19,400 L (RSE 29%)
    lkns <- log(56.5); label("Linear protein-binding ratio theta_pb, bound/unbound plasma MPA (unitless)") # Table 3 theta_pb = 56.5 (RSE 3%)
    lclin <- log(1200); label("Apparent influx clearance from unbound plasma into PBMC CLin/F (L/h)") # Table 3 CLin/F = 1200 L/h (RSE 12%)
    lclef <- log(43.8); label("Apparent efflux clearance from PBMC to unbound plasma CLout/F (L/h)") # Table 3 CLout/F = 43.8 L/h (RSE 16%)
    lvpbmc <- log(1980); label("Apparent PBMC volume Vcell/F (L)") # Table 3 Vcell/F = 1980 L (RSE 18%)

    # Covariate effects (Table 3 'Covariate model'; forms Eq. 5 and Eq. 6)
    e_crcl_cl <- 0.38; label("Power exponent of CrCL on CLu/F, reference 54.81 mL/min (unitless)") # Table 3 beta CLu/F,CrCL = 0.38 (RSE 19%)
    e_crcl_vp <- -1.03; label("Power exponent of CrCL on Vpu/F, reference 54.81 mL/min (unitless)") # Table 3 beta Vpu/F,CrCL = -1.03 (RSE 40%)
    e_alb_kns <- 1.46; label("Power exponent of baseline albumin on theta_pb, reference 37 g/L (unitless)") # Table 3 beta theta_pb,HSA = 1.46 (RSE 15%)
    e_snp_abcb1_rs1045642_hom_clef <- -0.64; label("log-scale effect of ABCB1 3435TT on CLout/F (unitless; exp(-0.64) = 0.527)") # Table 3 beta CLout/F,ABCB1 = -0.64 (RSE 44%)

    # Inter-individual variability. Table 3 prints IIV as CV%; converted to
    # log-normal variances with omega^2 = log(CV^2 + 1). Vcu/F, Qu/F, CLin/F
    # and Vcell/F carry no IIV.
    etald1 ~ 0.17698 # Table 3 Tk0 IIV 44%: log(0.44^2 + 1)
    etalcl ~ 0.086178 # Table 3 CLu/F IIV 30%: log(0.30^2 + 1)
    etalvp ~ 0.39878 # Table 3 Vpu/F IIV 70%: log(0.70^2 + 1)
    etalkns ~ 0.025278 # Table 3 theta_pb IIV 16%: log(0.16^2 + 1)
    etalclef ~ 0.39878 # Table 3 CLout/F IIV 70%: log(0.70^2 + 1)

    # Inter-occasion variability (Methods Eq. 4, kappa_ik ~ N(0, Gamma)), one
    # eta per occasion; occasions 2-4 are fixed to the occasion-1 variance so
    # the four slots share one estimate. Table 3 prints IOV as CV%, converted
    # as for IIV.
    etaiov_d1_1 ~ 0.30748 # Table 3 Tk0 IOV 60%: log(0.60^2 + 1)
    etaiov_d1_2 ~ fixed(0.30748)
    etaiov_d1_3 ~ fixed(0.30748)
    etaiov_d1_4 ~ fixed(0.30748)
    etaiov_cl_1 ~ 0.051548 # Table 3 CLu/F IOV 23%: log(0.23^2 + 1)
    etaiov_cl_2 ~ fixed(0.051548)
    etaiov_cl_3 ~ fixed(0.051548)
    etaiov_cl_4 ~ fixed(0.051548)
    etaiov_q_1 ~ 1.20624 # Table 3 Qu/F IOV 153%: log(1.53^2 + 1)
    etaiov_q_2 ~ fixed(1.20624)
    etaiov_q_3 ~ fixed(1.20624)
    etaiov_q_4 ~ fixed(1.20624)
    etaiov_vp_1 ~ 0.58339 # Table 3 Vpu/F IOV 89%: log(0.89^2 + 1)
    etaiov_vp_2 ~ fixed(0.58339)
    etaiov_vp_3 ~ fixed(0.58339)
    etaiov_vp_4 ~ fixed(0.58339)
    etaiov_clef_1 ~ 0.60328 # Table 3 CLout/F IOV 91%: log(0.91^2 + 1)
    etaiov_clef_2 ~ fixed(0.60328)
    etaiov_clef_3 ~ fixed(0.60328)
    etaiov_clef_4 ~ fixed(0.60328)
    etaiov_vpbmc_1 ~ 0.10337 # Table 3 Vcell/F IOV 33%: log(0.33^2 + 1)
    etaiov_vpbmc_2 ~ fixed(0.10337)
    etaiov_vpbmc_3 ~ fixed(0.10337)
    etaiov_vpbmc_4 ~ fixed(0.10337)

    # Residual error: proportional on each concentration (Results section 3.2.1)
    propSd <- 0.29; label("Proportional residual error, total plasma MPA (fraction)") # Table 3 sigma_t = 29% (RSE 4%)
    propSd_Cu <- 0.28; label("Proportional residual error, unbound plasma MPA (fraction)") # Table 3 sigma_u = 28% (RSE 6%)
    propSd_Cpbmc <- 0.39; label("Proportional residual error, PBMC MPA (fraction)") # Table 3 sigma_cell = 39% (RSE 16%)
  })

  model({
    # Occasion indicators for the inter-occasion variability
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    iov_d1 <- occ1 * etaiov_d1_1 + occ2 * etaiov_d1_2 + occ3 * etaiov_d1_3 + occ4 * etaiov_d1_4
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 + occ4 * etaiov_cl_4
    iov_q <- occ1 * etaiov_q_1 + occ2 * etaiov_q_2 + occ3 * etaiov_q_3 + occ4 * etaiov_q_4
    iov_vp <- occ1 * etaiov_vp_1 + occ2 * etaiov_vp_2 + occ3 * etaiov_vp_3 + occ4 * etaiov_vp_4
    iov_clef <- occ1 * etaiov_clef_1 + occ2 * etaiov_clef_2 + occ3 * etaiov_clef_3 + occ4 * etaiov_clef_4
    iov_vpbmc <- occ1 * etaiov_vpbmc_1 + occ2 * etaiov_vpbmc_2 + occ3 * etaiov_vpbmc_3 + occ4 * etaiov_vpbmc_4

    # Individual parameters. Continuous covariates enter as powers of
    # COV / COV_REF (Eq. 5); the categorical ABCB1 effect is log-additive
    # (Eq. 6 read as exp(beta * indicator), see the vignette).
    d1 <- exp(ld1 + etald1 + iov_d1)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl + iov_cl) * (CRCL / 54.81)^e_crcl_cl
    q <- exp(lq + iov_q)
    vp <- exp(lvp + etalvp + iov_vp) * (CRCL / 54.81)^e_crcl_vp
    kns <- exp(lkns + etalkns) * (ALB / 37)^e_alb_kns
    clin <- exp(lclin)
    clef <- exp(lclef + etalclef + iov_clef + e_snp_abcb1_rs1045642_hom_clef * SNP_ABCB1_RS1045642_HOM)
    vpbmc <- exp(lvpbmc + iov_vpbmc)

    # The central compartment holds unbound MPA (Figure 2): its volume Vcu/F
    # relates the state to the UNBOUND concentration, and elimination,
    # distribution and PBMC uptake are all driven by that concentration.
    Cu <- central / vc
    Cpbmc <- pbmc / vpbmc

    d/dt(central) <- -cl * Cu - q * Cu + q * peripheral1 / vp - clin * Cu + clef * Cpbmc
    d/dt(peripheral1) <- q * Cu - q * peripheral1 / vp
    d/dt(pbmc) <- clin * Cu - clef * Cpbmc

    # Zero-order absorption straight into the central compartment; dose
    # records must carry rate = -2.
    dur(central) <- d1

    # Total plasma MPA = unbound + bound = (1 + theta_pb) * unbound (Eq. 1)
    Cc <- (1 + kns) * Cu

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
    Cpbmc ~ prop(propSd_Cpbmc)
  })
}
