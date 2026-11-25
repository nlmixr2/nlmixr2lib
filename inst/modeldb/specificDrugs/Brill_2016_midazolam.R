Brill_2016_midazolam <- function() {
  description <- paste(
    "Semi-physiologically based PK model for midazolam and its CYP3A",
    "metabolite 1-OH-midazolam (blood concentrations) after oral and",
    "intravenous midazolam in 20 morbidly obese patients before and 18 of them",
    "one year after bariatric surgery. Gut wall, portal vein and liver are at",
    "quasi-steady state: gut-wall extraction by the Qgut model and hepatic",
    "extraction by the well-stirred model, both from unbound intrinsic",
    "clearances, so that 1-OH-midazolam is formed pre-systemically as well as",
    "systemically. Midazolam has three disposition compartments and a",
    "depot-plus-five-transit oral absorption chain (Ka = Ktr); 1-OH-midazolam",
    "has two. Hepatic blood flow scales with total body weight (allometric",
    "scenario 1 of the paper). After surgery the hepatic intrinsic clearance",
    "and Ka are estimated separately and Q is multiplied by 3.22; before",
    "surgery body weight scales the midazolam central and peripheral volumes.",
    "CAUTION: the second midazolam inter-compartmental clearance Q2 is not",
    "reported in the paper; it is fixed at 0.5 L/min, the INITIAL estimate in",
    "the authors' published NONMEM code, not a fitted value."
  )
  reference <- paste(
    "Brill MJE, Valitalo PAJ, Darwich AS, van Ramshorst B, van Dongen HPA,",
    "Rostami-Hodjegan A, Danhof M, Knibbe CAJ (2016). Semiphysiologically",
    "based pharmacokinetic model for midazolam and CYP3A mediated metabolite",
    "1-OH-midazolam in morbidly obese and weight loss surgery patients.",
    "CPT Pharmacometrics Syst Pharmacol 5(1):20. doi:10.1002/psp4.12048."
  )
  vignette <- "Brill_2016_midazolam"
  units <- list(time = "min", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Brill 2016 S3 code TBW. Sets hepatic blood flow at every occasion,",
        "QH = 3.75 * WT^0.75 / 60 L/min (Table 1 scenario 1), and from it the",
        "portal-vein, hepatic-artery and villous flows. Before surgery",
        "(PRIOR_BARIATRIC = 0) it also scales the midazolam central volume",
        "linearly, Vc = 66.9 * (1 + 0.0435 * (WT - 127)), and the first",
        "peripheral volume (and with it the second) as a power,",
        "Vp = 31.0 * (WT / 127)^3.93. The S3 code centres both at 127 kg; Table",
        "3 prints 144 kg (see the article's Errata). The linear Vc term turns",
        "negative below 104 kg, outside the observed pre-surgery range of",
        "112-186 kg. Record the weight measured at each occasion."
      ),
      source_name = "TBW"
    ),
    PRIOR_BARIATRIC = list(
      description = "Bariatric surgery performed before this record (1 = post-surgery, 0 = before surgery)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (morbidly obese patient before bariatric surgery)",
      notes = paste(
        "Brill 2016 S3 code OBES (1 = morbidly obese, 2 = 'Bariatric surgery",
        "patients') and OCC (1, 2); PRIOR_BARIATRIC = OBES - 1 = OCC - 1.",
        "0 = the day of laparoscopic Roux-en-Y gastric bypass or sleeve",
        "gastrectomy, 1 = one year later. Time-varying within subject; the",
        "random effects are shared by both occasions (no inter-occasion",
        "variability in the model)."
      ),
      source_name = "OBES"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "midazolam", units = "ug", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "ug", specimen = "whole blood", verified = TRUE),
    peripheral2 = list(analyte = "midazolam", units = "ug", specimen = "whole blood", verified = TRUE),
    central_1ohm = list(
      analyte = "1-OH-midazolam (midazolam-mass equivalents)",
      units = "ug",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1_1ohm = list(
      analyte = "1-OH-midazolam (midazolam-mass equivalents)",
      units = "ug",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "Before surgery 26-57 years (mean 43.6, SD 7.6); after surgery 27-58 years (mean 45.5, SD 7.4)",
    weight_range = "Before surgery 112-186 kg (mean 144.4, SD 21.7); after surgery 62-138 kg (mean 98.3, SD 18.0)",
    sex_female_pct = 60,
    disease_state = paste(
      "Morbid obesity (BMI before surgery 40-68 kg/m^2, mean 47.1) in adults",
      "undergoing laparoscopic Roux-en-Y gastric bypass or sleeve gastrectomy;",
      "restudied one year after surgery (BMI 24-50 kg/m^2, mean 31.9; mean",
      "weight loss 44.5 kg). No CYP3A inducers or inhibitors, no renal",
      "insufficiency."
    ),
    dose_range = paste(
      "At each occasion 7.5 mg oral midazolam followed by 5 mg IV",
      "160 +/- 48 min later; a mean of 22 blood samples per patient and",
      "occasion, assayed for midazolam and 1-OH-midazolam."
    ),
    regions = "Netherlands",
    notes = paste(
      "Brill 2016 Table 2: 20 patients (12 F / 8 M) before surgery, 18 of them",
      "(11 F / 7 M) one year after. NTC01519726 / EudraCT 2011-003293-93.",
      "The same cohort as Brill 2015 (doi:10.1007/s11095-015-1752-9, midazolam",
      "only); occasion 1 was also analysed in Brill 2014 (Brill_2014_midazolam)."
    )
  )

  ini({
    # Final semi-PBPK model, hepatic blood flow scenario 1 (Brill 2016
    # Table 3, 'Value (RSE %)'); structure and fixed constants from the
    # authors' NONMEM code (Supplementary Information S3).

    # Oral absorption, Ka = Ktr, estimated separately per surgery status
    lka_presurg <- log(0.126); label("Ka = Ktr before bariatric surgery (1/min)") # Table 3 'Ka morbidly obese = Ktr' 0.126 (10%); S3 THETA(4)
    lka_postsurg <- log(0.242); label("Ka = Ktr after bariatric surgery (1/min)") # Table 3 'Ka weight loss patients = Ktr' 0.242 (9%); S3 THETA(22)

    # Midazolam unbound intrinsic clearances (blood basis)
    lcl_int_h_presurg <- log(16.8); label("Midazolam hepatic unbound intrinsic clearance before surgery (L/min)") # Table 3 'CLint,H morbidly obese' 16.8 (14%); S3 THETA(7)
    lcl_int_h_postsurg <- log(25.5); label("Midazolam hepatic unbound intrinsic clearance after surgery (L/min)") # Table 3 'CLint,H weight loss patients' 25.5 (15%); S3 THETA(23)
    lcl_int_g <- log(0.0199); label("Midazolam gut-wall unbound intrinsic clearance (L/min)") # Table 3 'CLint,G' 0.0199 (35%); S3 THETA(9)

    # Midazolam distribution
    lvc <- log(66.9); label("Midazolam central volume at WT = 127 kg before surgery, and after surgery (L)") # Table 3 'Vcentral weight loss patients' 66.9 (13%); S3 THETA(1)
    e_wt_vc <- fixed(0.0435); label("Linear WT slope on midazolam Vc before surgery (1/kg)") # Table 3 X 0.0435 FIX; S3 THETA(21)
    lvp <- log(31.0); label("Midazolam first peripheral volume at WT = 127 kg before surgery, and after surgery (L)") # Table 3 'Vperi 1 weight loss patients' 31.0 (19%); S3 THETA(3)
    e_wt_vp <- fixed(3.93); label("Power exponent of WT/127 on midazolam peripheral volumes before surgery (unitless)") # Table 3 Y 3.93 FIX; S3 THETA(20)
    # Table 3 'Vperi 2 = Vperi 1 * Z', Z = 10.8 (13%); S3 V12 = THETA(17) * V6.
    # Carried as an absolute volume (10.8 x 31.0 L) that shares the Vp random
    # effect and weight scaling, which reproduces Vp2 = Z * Vp exactly.
    lvp2 <- log(10.8 * 31.0); label("Midazolam second peripheral volume at WT = 127 kg before surgery, and after surgery (L)") # Table 3 Z 10.8 x Vperi 1 31.0 L
    lq <- log(1.41); label("Midazolam inter-compartmental clearance to peripheral1 before surgery (L/min)") # Table 3 'Q1' 1.41 (15%); S3 THETA(2)
    e_prior_bariatric_q <- fixed(3.22); label("Multiplicative factor on midazolam Q after bariatric surgery (unitless)") # Table 3 A 3.22 FIX; S3 THETA(24) 'fQ BA' applied as THETA(2)*THETA(24) when OBES = 2
    # NOT A FITTED VALUE. S3 'Q12=THETA(18)' is estimated, but its final
    # value is printed nowhere in the paper or supplements; 0.5 is the S3
    # $THETA initial estimate '(0, 0.5) ;18, Q12'.
    lq2 <- fixed(log(0.5)); label("Midazolam inter-compartmental clearance to peripheral2 (L/min; S3 initial estimate, final value unreported)") # S3 $THETA 18 initial 0.5; unreported final

    # 1-OH-midazolam (no IIV: S3 $OMEGA 8-10 '0 FIX')
    lvc_1ohm <- log(41.7); label("1-OH-midazolam central volume (L)") # Table 3 'Vcentral, 1-OH' 41.7 (11%); S3 THETA(10)
    lvp_1ohm <- log(16.4); label("1-OH-midazolam peripheral volume (L)") # Table 3 'Vperi, 1-OH' 16.4 (25%); S3 THETA(15)
    lq_1ohm <- log(0.652); label("1-OH-midazolam inter-compartmental clearance (L/min)") # Table 3 'Q1-OH' 0.652 (23%); S3 THETA(16)
    lcl_int_h_1ohm <- log(27.4); label("1-OH-midazolam hepatic unbound intrinsic clearance (L/min)") # Table 3 'CLint,H,1-OH' 27.4 (9%); S3 THETA(11)
    lcl_int_g_1ohm <- log(11.9); label("1-OH-midazolam gut-wall unbound intrinsic clearance (L/min)") # Table 3 'CLint,G,1-OH' 11.9 (180%); S3 THETA(13)

    # IIV, Table 3 '%'; omega^2 = log(1 + CV^2). S3 $OMEGA: no IIV on Q,
    # F1, ALAG1 or any 1-OH-midazolam parameter.
    etalka ~ log(1 + 0.44^2) # Table 3 IIV Ka 44% (20%); S3 ETA(11), shared by both surgery strata
    etalvc ~ log(1 + 0.63^2) # Table 3 IIV Vcentral 63% (49%); S3 ETA(1)
    etalvp ~ log(1 + 1.13^2) # Table 3 IIV Vperipheral 1 113% (24%); S3 ETA(3), also carried by Vp2 = Z * V6
    etalcl_int_h ~ log(1 + 0.48^2) # Table 3 IIV CLint,H 48% (21%); S3 ETA(6), shared by both surgery strata
    etalcl_int_g ~ log(1 + 4.93^2) # Table 3 IIV CLint,G 493% (35%); S3 ETA(7)

    # Proportional residual error, one SD per surgery stratum, shared by
    # midazolam and 1-OH-midazolam (S3 $ERROR Y = IPRED * (1 + ERR(1 or 2)))
    propSd_presurg <- 0.326; label("Proportional residual SD before surgery, both analytes (fraction)") # Table 3 'Morbidly obese patients' 32.6% (18%)
    propSd_postsurg <- 0.237; label("Proportional residual SD after surgery, both analytes (fraction)") # Table 3 'Weight loss patients' 23.7% (8%)
  })
  model({
    # PRIOR_BARIATRIC = 1 after surgery (S3 OBES = 2), 0 before (OBES = 1)
    pre <- 1 - PRIOR_BARIATRIC

    # Fixed physiological and drug constants (S3 $PK and $THETA FIX)
    fub <- 0.033 # S3 THETA(8) 0.033 FIX, midazolam unbound fraction in blood (Table 1 prints 0.0303)
    fu_gut <- 1 # Table 1 fuG = 1; S3 FUG = 1
    cl_per <- 0.1766 # S3 THETA(19) 0.1766 FIX, permeability clearance (L/min); Table 1 CLperm 0.177
    fub_1ohm <- 0.106 # S3 THETA(12) 0.106 FIX; Table 1 fuB,1-OH 0.106
    fu_gut_1ohm <- 1 # S3 THETA(14) 1 FIX; Table 1 fuG,1-OH = 1

    qh <- 3.75 * WT^0.75 / 60 # S3 QH, L/min (Table 1 scenario 1, Brown 1997)
    qpv <- 0.75 * qh # S3 QPV, portal vein
    qha <- 0.25 * qh # S3 QHA, hepatic artery
    qin <- 0.4 * qh # S3 QIN, small intestine
    qmu <- 0.8 * qin # S3 QMU, mucosa
    qvi <- 0.6 * qmu # S3 QVI, villi
    qgut <- qvi * cl_per / (qvi + cl_per) # Eq. 4

    # Individual parameters
    ka <- exp(lka_presurg * pre + lka_postsurg * PRIOR_BARIATRIC + etalka)
    ktr <- ka # S3 KTR = KA
    cl_int_h <- exp(lcl_int_h_presurg * pre + lcl_int_h_postsurg * PRIOR_BARIATRIC + etalcl_int_h)
    cl_int_g <- exp(lcl_int_g + etalcl_int_g)
    vc <- exp(lvc + etalvc) * (PRIOR_BARIATRIC + pre * (1 + e_wt_vc * (WT - 127)))
    wt_vp <- PRIOR_BARIATRIC + pre * (WT / 127)^e_wt_vp
    vp <- exp(lvp + etalvp) * wt_vp
    vp2 <- exp(lvp2 + etalvp) * wt_vp
    q <- exp(lq) * e_prior_bariatric_q^PRIOR_BARIATRIC
    q2 <- exp(lq2)
    vc_1ohm <- exp(lvc_1ohm)
    vp_1ohm <- exp(lvp_1ohm)
    q_1ohm <- exp(lq_1ohm)
    cl_int_h_1ohm <- exp(lcl_int_h_1ohm)
    cl_int_g_1ohm <- exp(lcl_int_g_1ohm)

    # Extraction ratios (Eqs. 1-3, 5, 6)
    eh <- cl_int_h * fub / (qh + cl_int_h * fub)
    fhep <- 1 - eh
    egut <- cl_int_g * fu_gut / (qgut + cl_int_g * fu_gut)
    fgut <- 1 - egut
    eh_1ohm <- cl_int_h_1ohm * fub_1ohm / (qh + cl_int_h_1ohm * fub_1ohm)
    fhep_1ohm <- 1 - eh_1ohm
    egut_1ohm <- cl_int_g_1ohm * fu_gut_1ohm / (qvi + cl_int_g_1ohm * fu_gut_1ohm)
    fgut_1ohm <- 1 - egut_1ohm

    # Quasi-steady-state gut wall, portal vein and liver (S3 $DES, S2).
    # The organ volumes are 1 L in the code, so these are concentrations.
    cb <- central / vc
    cb_1ohm <- central_1ohm / vc_1ohm
    c_gutwall <- ktr * transit5 / qvi
    c_portal <- (qvi * c_gutwall * fgut + qpv * cb) / qpv
    c_liver <- (qha * cb + qpv * c_portal) / qh
    c_gutwall_1ohm <- egut * c_gutwall
    c_portal_1ohm <- (qvi * c_gutwall_1ohm * fgut_1ohm + qpv * cb_1ohm) / qpv
    c_liver_1ohm <- (qha * cb_1ohm + qpv * c_portal_1ohm + qh * c_liver * eh) / qh

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(central) <- fhep * qh * c_liver - qh * cb -
      q * cb + q / vp * peripheral1 -
      q2 * cb + q2 / vp2 * peripheral2
    d/dt(peripheral1) <- q * cb - q / vp * peripheral1
    d/dt(peripheral2) <- q2 * cb - q2 / vp2 * peripheral2
    d/dt(central_1ohm) <- fhep_1ohm * qh * c_liver_1ohm - qh * cb_1ohm -
      q_1ohm * cb_1ohm + q_1ohm / vp_1ohm * peripheral1_1ohm
    d/dt(peripheral1_1ohm) <- q_1ohm * cb_1ohm - q_1ohm / vp_1ohm * peripheral1_1ohm

    # Blood concentrations, ug/L = ng/mL. 1-OH-midazolam is in
    # midazolam-mass equivalents: the S3 code applies no molecular-weight
    # conversion between parent and metabolite.
    propSd <- propSd_presurg * pre + propSd_postsurg * PRIOR_BARIATRIC
    Cc <- cb
    Cc_1ohm <- cb_1ohm
    Cc ~ prop(propSd)
    Cc_1ohm ~ prop(propSd)
  })
}
