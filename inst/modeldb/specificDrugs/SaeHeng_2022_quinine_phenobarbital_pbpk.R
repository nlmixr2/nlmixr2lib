SaeHeng_2022_quinine_phenobarbital_pbpk <- function() {
  description <- "PBPK (whole-body, SimBiology; 23 perfusion-limited tissues per drug). Quinine-phenobarbital drug-drug interaction in adults with cerebral malaria and seizures. Phenobarbital induces hepatic CYP3A4, UGT1A1, CYP2C19 and CYP2C9 (Emax/EC50 turnover on liver concentration), roughly doubling quinine clearance (CYP3A4 + UGT1A1) and requiring a doubled quinine regimen; quinine competitively inhibits CYP3A4 (Ki). Clearances, tissue partitioning (Poulin-Theil) and organ sizes (Bosgra anthropometry) are algebraic functions of age and body weight. Deterministic typical-value model (the variability switch novarsw is 0); the published simulations add Bosgra inter-individual variability. Both drugs given by i.v. infusion. Structure, parameters and dosing reproduced from the deposited SimBiology project (supplement s005) and Tables S6/S7."
  reference <- "Sae-heng T, Rajoli RKR, Siccardi M, Karbwang J, Na-Bangchang K. Physiologically based pharmacokinetic modeling for dose optimization of quinine-phenobarbital coadministration in patients with cerebral malaria. CPT Pharmacometrics Syst Pharmacol. 2022;11:104-115. doi:10.1002/psp4.12737 (PMC8752110). Model input parameters: Table S6 (physicochemistry, clearances, induction Emax/EC50); physiology: Table S7 (Bosgra 2012 anthropometry). Structural ODEs, rate laws, tissue-composition fractions and dosing decoded from the deposited MATLAB SimBiology project, supplement PSP4-11-104-s005.zip."
  vignette <- "SaeHeng_2022_quinine_phenobarbital"

  # 23 perfusion-limited tissue states per drug plus the IV-infusion holding
  # depot, urine and cumulative-metabolism sinks. The phenobarbital subsystem
  # carries the '_pb' suffix; the bare states are quinine. heartc_l / heartc_r
  # are the two functional sides of the heart (venous -> heartc_r -> lung ->
  # heartc_l -> arterial); gonads, sc_site (subcutaneous) and im_site
  # (intramuscular) are decoded SimBiology compartments (SC / IM routes unused
  # on the i.v. regimen). None is a registered canonical, so all are declared here.
  paper_specific_compartments <- c(
    "a_metabolized_pb",
    "a_remainder_pb",
    "adipose_pb",
    "arterial_pb",
    "bone_pb",
    "brain_pb",
    "depot_iv_pb",
    "gonads",
    "gonads_pb",
    "gut_pb",
    "heartc_l",
    "heartc_l_pb",
    "heartc_r",
    "heartc_r_pb",
    "im_site",
    "im_site_pb",
    "kidney_pb",
    "liver_pb",
    "lung_pb",
    "muscle_pb",
    "pancreas_pb",
    "sc_site",
    "sc_site_pb",
    "skin_pb",
    "spleen_pb",
    "stomach_pb",
    "urine_pb",
    "venous_pb"
  )

  compartmentData <- list(
    a_metabolized = list(
      analyte = "quinine (cumulative metabolized amount)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_metabolized_pb = list(
      analyte = "phenobarbital (cumulative metabolized amount)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_remainder = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    a_remainder_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    adipose_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "quinine", units = "mg/L", specimen = "whole blood", verified = TRUE),
    arterial_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "whole blood", verified = TRUE),
    bone = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    bone_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    brain_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    depot_iv = list(analyte = "quinine", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv_pb = list(analyte = "phenobarbital", units = "mg", specimen = "administration site", verified = TRUE),
    gonads = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    gonads_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    gut_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    heartc_l = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    heartc_l_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    heartc_r = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    heartc_r_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    im_site = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    im_site_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    kidney_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    liver_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    lung_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    muscle_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    pancreas = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    pancreas_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    sc_site = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    sc_site_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    skin_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    spleen_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    stomach = list(analyte = "quinine", units = "mg/L", specimen = "tissue", verified = TRUE),
    stomach_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "quinine", units = "mg", specimen = "urine", verified = TRUE),
    urine_pb = list(analyte = "phenobarbital", units = "mg", specimen = "urine", verified = TRUE),
    venous = list(analyte = "quinine", units = "mg/L", specimen = "whole blood", verified = TRUE),
    venous_pb = list(analyte = "phenobarbital", units = "mg/L", specimen = "whole blood", verified = TRUE)
  )

  units <- list(
    time = "h",
    dosing = "mg (quinine base and phenobarbital, i.v.)",
    concentration = "mg/L (Cc quinine plasma, Cc_pb phenobarbital plasma)"
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the Bosgra 2012 organ-weight equations and the age-dependent",
        "microsomal-protein-per-gram-liver (MPPGL) function; simulated cohort",
        "18-60 y (Virtual population). Enters the model as AGE."
      ),
      source_name = "XAge"
    )
  )

  # Body weight is fixed at 60 kg in the published simulations (XWeight is set
  # to 60 by an initial-assignment rule inside model()); it is therefore not a
  # covariate column here. CYP2C19 genotype enters the published simulations
  # only through the phenobarbital intrinsic-clearance scaling of the intermediate
  # / poor-metabolizer virtual sub-cohorts (Table 1/2 strata); the deposited
  # typical-value model carries the extensive-metabolizer (wild-type) chain, so
  # no CYP2C19 covariate column is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Fixed at 60 kg in the published simulations (XWeight <- 60 in model()); not fitted or varied."
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      notes = "Phenobarbital IM/PM strata (Tables 1-2) scale phenobarbital CLint; the deposited typical-value chain is the extensive-metabolizer form, so no phenotype column is referenced in model()."
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      notes = "See CYP2C19_IM; poor-metabolizer stratum reduces phenobarbital CLint by 21-27% (Discussion), not encoded in the deposited extensive-metabolizer typical-value model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    disease_state = "cerebral malaria with concurrent seizures (+/- acute renal failure, lactic acidosis)",
    age_range = "18-60 years",
    weight_median = "60 kg (fixed in simulations)",
    sex_female_pct = 50,
    dose_range = paste(
      "Quinine 2000 mg i.v. loading (rate 250 mg/h) then 1200 mg i.v. (rate 150 mg/h)",
      "8-hourly; phenobarbital 90 mg i.v. daily (30-min infusion) for 17 days,",
      "quinine started on day 14 at phenobarbital steady state."
    ),
    regions = "Thailand (malaria-endemic)",
    notes = paste(
      "Virtual cohort of 100 (50 M / 50 F). CYP2C19 extensive, intermediate and",
      "poor metabolizer sub-cohorts simulated (Tables 1-2); the deposited",
      "typical-value model is the extensive-metabolizer (wild-type) chain."
    )
  )

  ini({
    # Quinine physicochemistry and clearance -- supplement Table S6 + decoded SimBiology 'Quinine' variant
    fu <- fixed(0.07)
    label("Quinine fraction unbound in plasma; Table S6 'Fraction of unbound drug' = 0.109 (whole-blood-referenced 0.07 used as fu, R=1.13)")
    pKa <- fixed(8.5)
    label("Quinine pKa; Table S6 = 8.5")
    Kpc <- fixed(3.44)
    label("Quinine logP (Kpc = logP); Table S6 'Log P' = 3.44")
    CL_liver_CYP3A4 <- fixed(11.6)
    label("Quinine hepatic CYP3A4 intrinsic-clearance base (CLinvitro*Abundance); decoded model Note, back-calc from in-vivo CL 4.86 L/h (Table S6), fm,CYP3A4=0.44")
    CL_liver_UGT1A1 <- fixed(14.22)
    label("Quinine hepatic UGT1A1 intrinsic-clearance base; decoded model Note, fm,UGT1A1=0.56")
    Ki3A4 <- fixed(17.84)
    label("Quinine CYP3A4 competitive-inhibition constant; Table S6 Ki3A4 = 12.65 mg/L (decoded as-run 17.84)")
    Vd_correct_fact <- fixed(0.5)
    label("Quinine Vss correction factor; decoded SimBiology 'Vd_correct_fact' = 0.5")

    # Phenobarbital physicochemistry and clearance -- Table S6 + decoded 'Phenobarbital' variant
    fu_pb <- fixed(0.49)
    label("Phenobarbital fraction unbound; Table S6 = 0.49")
    pKa_pb <- fixed(7.3)
    label("Phenobarbital pKa; Table S6 = 7.3")
    Kpc_pb <- fixed(1.47)
    label("Phenobarbital logP (Kpc_pb); Table S6 'Log P' = 1.47")
    bpr_pb <- fixed(0.83)
    label("Phenobarbital blood:plasma ratio; Table S6 R_bp = 0.83")
    bpr <- fixed(1.13)
    label("Quinine blood:plasma ratio (bpr); Table S6 R_bp = 1.13")
    CLint_liver_pb <- fixed(0.124)
    label("Phenobarbital hepatic intrinsic-clearance base; decoded model, back-calc from in-vivo CL 0.31 L/h (Table S6)")
    Vd_correct_fact_pb <- fixed(0.8)
    label("Phenobarbital Vss correction factor; decoded SimBiology = 0.8")
    pKa_2 <- fixed(1.0)
    label("Phenobarbital second-ionization pKa term (decoded D_vo_w_pb formula); as-run = 1")

    # Phenobarbital enzyme-induction Emax / EC50 -- Table S6 (as-run values from the decoded model where noted)
    Emax3A4_pb <- fixed(9.2)
    label("Phenobarbital CYP3A4 induction Emax; Table S6 E_max = 20.1 (decoded as-run 9.2)")
    Ec503A4_pb <- fixed(75.0)
    label("Phenobarbital CYP3A4 induction EC50 (mg/L); Table S6 EC50 = 109.88 (decoded as-run 75)")
    ECmaxUGT1A1 <- fixed(2.8)
    label("Phenobarbital UGT1A1 induction Emax; Table S6 = 2.8")
    EC50UGT1A1 <- fixed(232.0)
    label("Phenobarbital UGT1A1 induction EC50 (mg/L); Table S6 = 232")
    ECmax2C19 <- fixed(1.6)
    label("Phenobarbital CYP2C19 induction Emax; Table S6 E_max = 8.0 (decoded as-run 1.6)")
    EC502C19 <- fixed(116.0)
    label("Phenobarbital CYP2C19 induction EC50 (mg/L); Table S6 = 11.6 (decoded as-run 116)")

    # Dosing-kinetics, subject and variability switches -- decoded SimBiology dose objects / variants
    kinf <- fixed(250.0)
    label("Quinine IV-infusion holding-compartment transfer rate (/h); decoded MassAction kf = 250")
    kinf_pb <- fixed(180.0)
    label("Phenobarbital IV-infusion holding-compartment transfer rate (/h); decoded kf = 180")
    kurine <- fixed(0.012)
    label("Quinine renal-excretion rate from kidney (fraction of kf); decoded 0.012")
    kurine_pb <- fixed(0.045)
    label("Phenobarbital renal-excretion rate; decoded 0.045")
    pH <- fixed(7.4)
    label("Physiological pH for ionization / Poulin-Theil partitioning; decoded model = 7.4")
    sub_area <- fixed(0.002)
    label("Subcutaneous depot fractional area of adipose (SC route unused here); decoded = 0.002")
    im_area <- fixed(0.004)
    label("Intramuscular depot fractional area of muscle (IM route unused here); decoded = 0.004")
    novarsw <- fixed(0.0)
    label("Variability switch (0 = typical-value deterministic solve; decoded stochastic runs set 1)")
    varsw <- fixed(1.0)
    label("Complementary variability switch (1 - novarsw); 1 selects the population-mean enzyme abundance")

    # Poulin-Theil tissue-composition fractions -- decoded 'TestSubject:Human' variant (Rodgers-Rowland)
    Vadn <- fixed(0.79)
    label("Poulin-Theil tissue-composition fraction 'Vadn'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vadp <- fixed(0.002)
    label("Poulin-Theil tissue-composition fraction 'Vadp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vadw <- fixed(0.18)
    label("Poulin-Theil tissue-composition fraction 'Vadw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vbon <- fixed(0.074)
    label("Poulin-Theil tissue-composition fraction 'Vbon'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vbop <- fixed(0.0011)
    label("Poulin-Theil tissue-composition fraction 'Vbop'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vbow <- fixed(0.439)
    label("Poulin-Theil tissue-composition fraction 'Vbow'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vbrn <- fixed(0.051)
    label("Poulin-Theil tissue-composition fraction 'Vbrn'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vbrp <- fixed(0.0565)
    label("Poulin-Theil tissue-composition fraction 'Vbrp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vbrw <- fixed(0.77)
    label("Poulin-Theil tissue-composition fraction 'Vbrw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vgun <- fixed(0.0487)
    label("Poulin-Theil tissue-composition fraction 'Vgun'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vgup <- fixed(0.0163)
    label("Poulin-Theil tissue-composition fraction 'Vgup'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vguw <- fixed(0.718)
    label("Poulin-Theil tissue-composition fraction 'Vguw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vhen <- fixed(0.0115)
    label("Poulin-Theil tissue-composition fraction 'Vhen'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vhep <- fixed(0.0166)
    label("Poulin-Theil tissue-composition fraction 'Vhep'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vhew <- fixed(0.758)
    label("Poulin-Theil tissue-composition fraction 'Vhew'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vkin <- fixed(0.0207)
    label("Poulin-Theil tissue-composition fraction 'Vkin'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vkip <- fixed(0.0162)
    label("Poulin-Theil tissue-composition fraction 'Vkip'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vkiw <- fixed(0.783)
    label("Poulin-Theil tissue-composition fraction 'Vkiw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vlin <- fixed(0.0348)
    label("Poulin-Theil tissue-composition fraction 'Vlin'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vlip <- fixed(0.0252)
    label("Poulin-Theil tissue-composition fraction 'Vlip'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vliw <- fixed(0.751)
    label("Poulin-Theil tissue-composition fraction 'Vliw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vlun <- fixed(0.003)
    label("Poulin-Theil tissue-composition fraction 'Vlun'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vlup <- fixed(0.009)
    label("Poulin-Theil tissue-composition fraction 'Vlup'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vluw <- fixed(0.811)
    label("Poulin-Theil tissue-composition fraction 'Vluw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vmun <- fixed(0.0238)
    label("Poulin-Theil tissue-composition fraction 'Vmun'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vmup <- fixed(0.0072)
    label("Poulin-Theil tissue-composition fraction 'Vmup'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vmuw <- fixed(0.76)
    label("Poulin-Theil tissue-composition fraction 'Vmuw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vnlp <- fixed(0.0035)
    label("Poulin-Theil tissue-composition fraction 'Vnlp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vphp <- fixed(0.00225)
    label("Poulin-Theil tissue-composition fraction 'Vphp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vskn <- fixed(0.0284)
    label("Poulin-Theil tissue-composition fraction 'Vskn'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vskp <- fixed(0.0111)
    label("Poulin-Theil tissue-composition fraction 'Vskp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vskw <- fixed(0.718)
    label("Poulin-Theil tissue-composition fraction 'Vskw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vspn <- fixed(0.0201)
    label("Poulin-Theil tissue-composition fraction 'Vspn'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vspp <- fixed(0.0198)
    label("Poulin-Theil tissue-composition fraction 'Vspp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vspw <- fixed(0.788)
    label("Poulin-Theil tissue-composition fraction 'Vspw'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")
    Vwp <- fixed(0.945)
    label("Poulin-Theil tissue-composition fraction 'Vwp'; decoded SimBiology TestSubject:Human variant (Rodgers-Rowland composition)")

    # Partition coefficients held at 1 (organs without a Poulin-Theil rule in the decoded model)
    kp_gonads_pb <- fixed(1.0)
    label("Gonads tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")
    kp_pancreas <- fixed(1.0)
    label("Pancreas tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")
    kp_pancreas_pb <- fixed(1.0)
    label("Pancreas tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")
    kp_remainder <- fixed(1.0)
    label("Remainder tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")
    kp_remainder_pb <- fixed(1.0)
    label("Remainder tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")
    kp_stomach <- fixed(1.0)
    label("Stomach tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")
    kp_stomach_pb <- fixed(1.0)
    label("Stomach tissue:plasma partition coefficient held at 1 (no Poulin-Theil rule for this organ in the decoded model)")

    # As-run defaults from the decoded model (unused enzyme pathways; do not affect either drug's clearance)
    Ec502B6_pb <- fixed(1.0)
    label("As-run default from the decoded model ('Ec502B6_pb' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    Ec502C8_pb <- fixed(1.0)
    label("As-run default from the decoded model ('Ec502C8_pb' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    Emax2B6_pb <- fixed(1.0)
    label("As-run default from the decoded model ('Emax2B6_pb' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    Emax2C8_pb <- fixed(1.0)
    label("As-run default from the decoded model ('Emax2C8_pb' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    HBD <- fixed(1.0)
    label("As-run default from the decoded model ('HBD' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    Ki2B6 <- fixed(1.0)
    label("As-run default from the decoded model ('Ki2B6' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    Ki2C8 <- fixed(1.0)
    label("As-run default from the decoded model ('Ki2C8' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    PSA <- fixed(1.0)
    label("As-run default from the decoded model ('PSA' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
    fInd2C19_1_total <- fixed(1.0)
    label("As-run default from the decoded model ('fInd2C19_1_total' = 1.0); does not affect quinine or phenobarbital clearance (unused enzyme pathway)")
  })

  model({
    # --- organ weights, volumes, blood flows, physchem, Kp, clearances (SimBiology rules) ---
    XBMI <- exp((3.294783))
    EP <- (bpr-(1-0.45))/0.45
    EP_pb <- (bpr_pb-(1-0.45))/0.45
    fut <- 1/(1+(((1-fu)/fu)*0.5))
    fut_pb <- 1/(1+(((1-fu_pb)/fu_pb)*0.5))
    Kt1 <- (2.2)
    Kt2 <- (2.15)
    Kt3 <- (2.1)
    PC <- 10^Kpc
    QCC <- (15)
    Abundance_CYP3A4gut <- (70.5)
    WBrain <- (0.405*exp(-AGE/629)*(3.68-2.68*exp(-AGE/0.89)))
    WGonads <- (3.3+53*(1-exp((AGE/17.5)^5.4*cos(5.4*pi))))/1000
    WThymus <- (14*((7.1-6.1*exp(-AGE/11.9))*((0.14+0.86*exp(-AGE/10.3)))))/1000
    kp_bone <- (PC*(Vbon+0.3*Vbop)+(1*Vbow+0.7*Vbop))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_brain <- (PC*(Vbrn+0.3*Vbrp)+(1*(Vbrw+0.7*Vbrp)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_gut <- (PC*(Vgun+0.3*Vgup)+(1*(Vguw+0.7*Vgup)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_heart <- (PC*(Vhen+0.3*Vhep)+(1*(Vhew+0.7*Vhep)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_kidney <- (PC*(Vkin+0.3*Vkip)+(1*(Vkiw+0.7*Vkip)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_liver <- (PC*(Vlin+0.3*Vlip)+(1*(Vliw+0.7*Vlip)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_lung <- (PC*(Vlun+0.3*Vlup)+(1*Vluw+0.7*Vlup))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_muscle <- (PC*(Vmun+0.3*Vmup)+(1*(Vmuw+0.7*Vmup)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_skin <- (PC*(Vskn+0.3*Vskp)+(1*(Vskw+0.7*Vskp)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    kp_spleen <- (PC*(Vspn+0.3*Vspp)+(1*(Vspw+0.7*Vspp)))/(PC*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu/fut)*Vd_correct_fact
    Dpc <- 10^((1.115*log10(PC))-1.35)
    Abundance_CYP3A4 <- (abs(141-80*0.5+259*0.5)*novarsw) + (141*varsw)
    Abundance_CYP2D6 <- abs(0.8+0.5*0.6)
    Abundance_CYP3A5 <- abs((16))
    Abundance_CYP2C19 <- abs((14))
    Peff <- 10^(-2.546-0.011*PSA-0.278*HBD)
    MMPGL <- abs((10^(1.407+0.0158*AGE-0.00038*AGE^2+0.0000024*AGE^3)))
    Abundance_CYP1A1 <- (0.58)
    Abundance_CYP1A2 <- (52)
    Ka <- 2*Peff*60*60/(1.5)
    Abundance_CYP2C8 <- abs((30.8))
    Abundance_CYP2B6 <- (17)
    PC_pb <- 10^Kpc_pb
    kp_bone_pb <- (PC_pb*(Vbon+0.3*Vbop)+(1*Vbow+0.7*Vbop))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_brain_pb <- (PC_pb*(Vbrn+0.3*Vbrp)+(1*(Vbrw+0.7*Vbrp)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_gut_pb <- (PC_pb*(Vgun+0.3*Vgup)+(1*(Vguw+0.7*Vgup)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_heart_pb <- (PC_pb*(Vhen+0.3*Vhep)+(1*(Vhew+0.7*Vhep)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_kidney_pb <- (PC_pb*(Vkin+0.3*Vkip)+(1*(Vkiw+0.7*Vkip)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_liver_pb <- (PC_pb*(Vlin+0.3*Vlip)+(1*(Vliw+0.7*Vlip)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_lung_pb <- (PC_pb*(Vlun+0.3*Vlup)+(1*Vluw+0.7*Vlup))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_muscle_pb <- (PC_pb*(Vmun+0.3*Vmup)+(1*(Vmuw+0.7*Vmup)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_skin_pb <- (PC_pb*(Vskn+0.3*Vskp)+(1*(Vskw+0.7*Vskp)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    kp_spleen_pb <- (PC_pb*(Vspn+0.3*Vspp)+(1*(Vspw+0.7*Vspp)))/(PC_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*(fu_pb/fut_pb)*Vd_correct_fact_pb
    Abundance_UGT1A1 <- abs(max(8.9,min(137.9,(34))))
    ClintUGT1A1 <- CL_liver_UGT1A1/Abundance_UGT1A1
    Abundance_CYP2C19gut <- min(3.9,max(0.6,abs((2.1))))
    Abundance_CYP2C9 <- abs(max(8.49,min(87.20,(37.53))))
    Clint2C9 <- (CLint_liver_pb/Abundance_CYP2C9)*0.1
    Abundance_2C19active <- Abundance_CYP2C19*fInd2C19_1_total
    XWeight <- 60
    v_brain <- WBrain/1.035
    v_im_site <- 0.004*1.041
    v_gonads <- 1
    v_sc_site <- 1
    kp_gonads <- 1
    D_pb <- 10^((1.115*log10(PC_pb))-1.35)
    D_vo_w <- 10^(log10(Dpc)-log10(1+10^(-pH+pKa)))
    qc_total <- QCC*(XWeight^0.75)
    q_adipose <- qc_total*0.052
    q_bone <- qc_total*0.042
    q_brain <- qc_total*0.11
    q_ha <- qc_total*0.12
    q_kidney <- qc_total*0.175
    Qlu <- qc_total*0.025
    q_muscle <- qc_total*0.19
    q_remainder <- qc_total*0.01
    q_skin <- qc_total*0.06
    q_pv <- 0.20*qc_total
    q_gonads <- qc_total*0.01
    q_hv <- q_ha+q_pv
    Qgu <- q_pv/4
    Qpa <- q_pv/4
    Qsp <- q_pv/4
    Qst <- q_pv/4
    WAdipose <- (((((1.20*XBMI)+(0.23*AGE)-16.2)*XWeight)/100))
    kp_adipose <- (D_vo_w*(Vadn+0.3*Vadp)+(1*(Vadw+0.7*Vadp)))*(fu)/(D_vo_w*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*Vd_correct_fact
    XHeight <- sqrt(XWeight/XBMI)
    D_vo_w_pb <- 10^(log10(D_pb)-log10(1+10^(pH+pKa_pb+pH-pKa_2)))
    kp_adipose_pb <- (D_vo_w_pb*(Vadn+0.3*Vadp)+(1*(Vadw+0.7*Vadp)))*(fu_pb)/(D_vo_w_pb*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))*Vd_correct_fact_pb
    v_adipose <- WAdipose/0.916
    WLiver <- exp((-0.6786+1.98*log(XHeight)))
    XBSA <- (71.84*(XWeight^0.425)*((XHeight*100)^0.725))/10000
    WBlood <- (3.33*XBSA-0.81)
    WHeart <- exp((-2.502+2.13*log(XHeight)))
    WIntestines <- exp((-1.351+2.47*log(XHeight)))
    WKidneys <- exp((-2.306+1.93*log(XHeight)))
    WLungs <- exp((-2.092+2.1*log(XHeight)))
    WPancreas <- exp((-3.431+2.43*log(XHeight)))
    WRemaining <- exp((-0.072+1.95*log(XHeight)))
    WSkin <- (exp(1.64*XBSA-1.93))
    WSpleen <- exp((-3.123+2.16*log(XHeight)))
    WStomach <- exp((-3.266+2.45*log(XHeight)))
    WBones <- exp((0.0689+2.67*log(XHeight)))
    WTotal_weight1 <- WLungs+WHeart+WBones+WKidneys+WStomach+WIntestines+WSpleen+WPancreas+WLiver+WRemaining+WBrain+WSkin+WBlood+WAdipose+WThymus+WGonads
    WMuscle <- 0.93*XWeight-WTotal_weight1
    ClintCYP2C9active_pb <- (Clint2C9*Abundance_CYP2C9*MMPGL*1000*60*WLiver/1000000)
    v_venous <- WBlood*0.67
    v_stomach <- WStomach/1.05
    v_spleen <- WSpleen/1.054
    v_skin <- WSkin/1.1
    v_pancreas <- WPancreas/1.045
    v_muscle <- WMuscle/1.041
    v_lung <- WLungs/1.05
    v_liver <- WLiver
    v_kidney <- WKidneys/1.05
    v_heart <- WHeart/1.03
    v_gut <- WIntestines/1.05
    v_bone <- WBones/1.2
    v_arterial <- WBlood*0.33
    v_remainder <- WRemaining
    fInh2C8_pb <- 1+(liver_pb/v_liver/Ki2C8)
    fInd2C8_pb <- 1+ (Emax2C8_pb*liver_pb/v_liver)/(Ec502C8_pb+(liver_pb/v_liver))
    Abundance_CYP2C8active <- Abundance_CYP2C8*fInd2C8_pb/fInh2C8_pb
    fInd2B6_pb <- 1+ (Emax2B6_pb*liver_pb/v_liver)/(Ec502B6_pb+(liver_pb/v_liver))
    fInh2B6_pb <- 1+(liver_pb/v_liver/Ki2B6)
    Abundance_CYP2B6active <- Abundance_CYP2B6*fInd2B6_pb/fInh2B6_pb
    fInd3A4_pb <- (1+ (Emax3A4_pb*liver_pb/v_liver)/(Ec503A4_pb+(liver_pb/v_liver)))
    fInh3A4 <- 1+(liver/v_liver/Ki3A4)
    fIndUGT1A1 <- (1+ (ECmaxUGT1A1*liver_pb/v_liver)/(EC50UGT1A1+(liver_pb/v_liver)))
    Abundance_CYP3A4active_pb <- Abundance_CYP3A4/fInh3A4
    fInd2C19 <- 1+ (ECmax2C19*liver_pb/v_liver)/(EC502C19+(liver_pb/v_liver))
    Abundance_CYP2C19active_pb <- Abundance_CYP2C19*fInd2C19
    Abundance_CYP3A4active <- (Abundance_CYP3A4*fInd3A4_pb)/(fInh3A4)
    Clint3A4 <- CL_liver_CYP3A4/Abundance_CYP3A4active_pb
    Abundance_UGT1A1active <- Abundance_UGT1A1*fIndUGT1A1
    ClintUGT1A1active <- ClintUGT1A1*Abundance_UGT1A1active*MMPGL*1000*60*WLiver/1000000
    Clint2C19_pb <- (CLint_liver_pb/Abundance_CYP2C19active_pb)*0.9
    Clint_gut_pb <- (Clint2C19_pb*Abundance_CYP2C19gut*1000*60/1000000)
    Clint_gut <- Clint3A4*Abundance_CYP3A4gut*1000*60/1000000
    ClintCYP3A4active <- (Clint3A4*Abundance_CYP3A4active*MMPGL*1000*60*WLiver/1000000)
    ClintCYP2C19active_pb <- (Clint2C19_pb*Abundance_CYP2C19active_pb*MMPGL*1000*60*WLiver/1000000)*0.65
    Clearance1 <- ClintCYP3A4active+ClintUGT1A1active
    cl_blood <- (q_hv*(fu/bpr)*Clearance1/(q_hv+Clearance1*(fu/bpr)))
    Clearance1_pb <- ClintCYP2C19active_pb+ClintCYP2C9active_pb
    cl_blood_pb <- (q_hv*(fu_pb/bpr_pb)*Clearance1_pb/(q_hv+Clearance1_pb*(fu_pb/bpr_pb)))

    # --- perfusion fluxes (SimBiology reactions) ---
    J1 <- q_hv*liver/v_liver*bpr/kp_liver-q_hv*venous/v_venous
    J2 <- q_pv/4*pancreas/v_pancreas*bpr/kp_pancreas-q_pv/4*liver/v_liver*bpr/kp_liver
    J3 <- q_pv/4*spleen/v_spleen*bpr/kp_spleen-q_pv/4*liver/v_liver*bpr/kp_liver
    J4 <- q_brain*brain/v_brain*bpr/kp_brain-q_brain*venous/v_venous
    J5 <- q_brain*arterial/v_arterial-q_brain*brain/v_brain*bpr/kp_brain
    J6 <- q_pv/4*arterial/v_arterial-q_pv/4*spleen/v_spleen*bpr/kp_spleen
    J7 <- q_pv/4*arterial/v_arterial-q_pv/4*pancreas/v_pancreas*bpr/kp_pancreas
    J8 <- q_muscle*(1-im_area/v_muscle)*arterial/v_arterial-q_muscle*(1-im_area/v_muscle)*muscle/v_muscle*bpr/kp_muscle
    J9 <- q_muscle*(1-im_area/v_muscle)*muscle/v_muscle*bpr/kp_muscle-q_muscle*(1-im_area/v_muscle)*venous/v_venous
    J10 <- q_adipose*(1-sub_area/v_adipose)*adipose/v_adipose*bpr/kp_adipose-q_adipose*(1-sub_area/v_adipose)*venous/v_venous
    J11 <- q_adipose*(1-sub_area/v_adipose)*arterial/v_arterial-q_adipose*(1-sub_area/v_adipose)*adipose/v_adipose*bpr/kp_adipose
    J12 <- q_skin*skin/v_skin*bpr/kp_skin-q_skin*venous/v_venous
    J13 <- q_bone*bone/v_bone*bpr/kp_bone-q_bone*venous/v_venous
    J14 <- q_skin*arterial/v_arterial-q_skin*skin/v_skin*bpr/kp_skin
    J15 <- q_bone*arterial/v_arterial-q_bone*bone/v_bone*bpr/kp_bone
    J16 <- q_kidney*arterial/v_arterial-q_kidney*kidney/v_kidney*bpr/kp_kidney
    J17 <- q_kidney*kidney/v_kidney*bpr/kp_kidney-q_kidney*venous/v_venous
    J18 <- q_pv/4*arterial/v_arterial-q_pv/4*gut/v_gut*bpr/kp_gut
    J19 <- q_pv/4*arterial/v_arterial-q_pv/4*stomach/v_stomach*bpr/kp_stomach
    J20 <- q_pv/4*gut/v_gut*bpr/kp_gut-q_pv/4*liver/v_liver*bpr/kp_liver
    J21 <- q_pv/4*stomach/v_stomach*bpr/kp_stomach-q_pv/4*liver/v_liver*bpr/kp_liver
    J22 <- q_ha*arterial/v_arterial-q_ha*liver/v_liver*bpr/kp_liver
    J23 <- cl_blood*venous/v_venous
    J24 <- q_remainder*a_remainder/v_remainder*bpr/kp_remainder-q_remainder*venous/v_venous
    J25 <- q_remainder*arterial/v_arterial-q_remainder*a_remainder/v_remainder*bpr/kp_remainder
    J26 <- cl_blood_pb*venous_pb/v_venous
    J27 <- q_hv*liver_pb/v_liver*bpr_pb/kp_liver_pb-q_hv*venous_pb/v_venous
    J28 <- q_ha*arterial_pb/v_arterial-q_ha*liver_pb/v_liver*bpr_pb/(kp_liver_pb)
    J29 <- q_pv/4*arterial_pb/v_arterial-(q_pv/4)*stomach_pb/v_stomach*bpr_pb/(kp_stomach_pb)
    J30 <- q_pv/4*stomach_pb/v_stomach*bpr_pb/(kp_stomach_pb)-(q_pv/4)*liver_pb/v_liver*bpr_pb/(kp_liver_pb)
    J31 <- q_pv/4*gut_pb/v_gut*bpr_pb/(kp_gut_pb)-(q_pv/4)*liver_pb/v_liver*bpr_pb/(kp_liver_pb)
    J32 <- q_pv/4*arterial_pb/v_arterial-(q_pv/4)*gut_pb/v_gut*bpr_pb/(kp_gut_pb)
    J33 <- q_pv/4*spleen_pb/v_spleen*bpr_pb/(kp_spleen_pb)-(q_pv/4)*liver_pb/v_liver*bpr_pb/(kp_liver_pb)
    J34 <- q_pv/4*arterial_pb/v_arterial-q_pv/4*spleen_pb/v_spleen*bpr_pb/(kp_spleen_pb)
    J35 <- q_pv/4*pancreas_pb/v_pancreas*bpr_pb/(kp_pancreas_pb)-(q_pv/4)*liver_pb/v_liver*bpr_pb/(kp_liver_pb)
    J36 <- q_pv/4*arterial_pb/v_arterial-(q_pv/4)*pancreas_pb/v_pancreas*bpr_pb/(kp_pancreas_pb)
    J37 <- q_brain*brain_pb/v_brain*bpr_pb/(kp_brain_pb)-q_brain*venous_pb/v_venous
    J38 <- q_brain*arterial_pb/v_arterial-q_brain*brain_pb/v_brain*bpr_pb/(kp_brain_pb)
    J39 <- q_muscle*(1-im_area/v_muscle)*muscle_pb/v_muscle*bpr_pb/kp_muscle_pb-q_muscle*(1-im_area/v_muscle)*venous_pb/v_venous
    J40 <- q_muscle*(1-im_area/v_muscle)*arterial_pb/v_arterial-q_muscle*(1-im_area/v_muscle)*muscle_pb/v_muscle*bpr_pb/kp_muscle_pb
    J41 <- q_adipose*(1-sub_area/v_adipose)*adipose_pb/v_adipose*bpr_pb/kp_adipose_pb-q_adipose*(1-sub_area/v_adipose)*venous_pb/v_venous
    J42 <- q_adipose*(1-sub_area/v_adipose)*arterial_pb/v_arterial-q_adipose*(1-sub_area/v_adipose)*adipose_pb/v_adipose*bpr_pb/kp_adipose_pb
    J43 <- q_skin*skin_pb/v_skin*bpr_pb/(kp_skin_pb)-q_skin*venous_pb/v_venous
    J44 <- q_skin*arterial_pb/v_arterial-q_skin*skin_pb/v_skin*bpr_pb/(kp_skin_pb)
    J45 <- q_bone*bone_pb/v_bone*bpr_pb/(kp_bone_pb)-q_bone*venous_pb/v_venous
    J46 <- q_bone*arterial_pb/v_arterial-q_bone*bone_pb/v_bone*bpr_pb/(kp_bone_pb)
    J47 <- q_remainder*a_remainder_pb/v_remainder*bpr_pb/(kp_remainder_pb)-q_remainder*venous_pb/v_venous
    J48 <- q_remainder*arterial_pb/v_arterial-q_remainder*a_remainder_pb/v_remainder*bpr_pb/(kp_remainder_pb)
    J49 <- q_kidney*kidney_pb/v_kidney*bpr_pb/(kp_kidney_pb)-q_kidney*venous_pb/v_venous
    J50 <- q_kidney*arterial_pb/v_arterial-q_kidney*kidney_pb/v_kidney*bpr_pb/(kp_kidney_pb)
    J51 <- qc_total*heartc_r/v_heart*bpr/kp_heart-qc_total*lung/v_lung*bpr/kp_lung
    J52 <- qc_total*heartc_l/v_heart*bpr/kp_heart-(qc_total)*arterial/v_arterial
    J53 <- qc_total*heartc_r_pb/v_heart*bpr_pb/kp_heart_pb-qc_total*lung_pb/v_lung*bpr_pb/kp_lung_pb
    J54 <- qc_total*heartc_l_pb/v_heart*bpr_pb/kp_heart_pb-(qc_total)*arterial_pb/v_arterial
    J55 <- qc_total*venous/v_venous-qc_total*heartc_r/v_heart*bpr/kp_heart
    J56 <- qc_total*venous_pb/v_venous-qc_total*heartc_r_pb/v_heart*bpr_pb/kp_heart_pb
    J57 <- qc_total*lung/v_lung*bpr/kp_lung-qc_total*heartc_l/v_heart*bpr/kp_heart
    J58 <- qc_total*lung_pb/v_lung*bpr_pb/kp_lung_pb-qc_total*heartc_l_pb/v_heart*bpr_pb/kp_heart_pb
    J59 <- q_adipose*(sub_area/v_adipose)*arterial/v_arterial-q_adipose*(sub_area/v_adipose)*sc_site/v_sc_site*bpr/kp_adipose
    J60 <- q_adipose*(sub_area/v_adipose)*sc_site/v_sc_site*bpr/kp_adipose-q_adipose*(sub_area/v_adipose)*venous/v_venous
    J61 <- q_gonads*arterial/v_arterial-q_gonads*gonads/v_gonads*bpr/kp_gonads
    J62 <- q_gonads*gonads/v_gonads*bpr/kp_gonads-q_gonads*venous/v_venous
    J63 <- q_muscle*(im_area/v_muscle)*arterial/v_arterial-q_muscle*(im_area/v_muscle)*im_site/v_im_site*bpr/kp_muscle
    J64 <- q_muscle*(im_area/v_muscle)*im_site/v_im_site*bpr/kp_muscle-q_muscle*(im_area/v_muscle)*venous/v_venous
    J65 <- q_gonads*gonads_pb/v_gonads*bpr_pb/kp_gonads_pb-q_gonads*venous_pb/v_venous
    J66 <- q_gonads*arterial_pb/v_arterial-q_gonads*gonads_pb/v_gonads*bpr_pb/kp_gonads_pb
    J67 <- q_muscle*(im_area/v_muscle)*arterial_pb/v_arterial-q_muscle*(im_area/v_muscle)*im_site_pb/v_im_site*bpr_pb/kp_muscle_pb
    J68 <- q_muscle*(im_area/v_muscle)*im_site_pb/v_im_site*bpr_pb/kp_muscle_pb-q_muscle*(im_area/v_muscle)*venous_pb/v_venous
    J69 <- q_adipose*(sub_area/v_adipose)*sc_site_pb/v_sc_site*bpr_pb/kp_adipose_pb-q_adipose*(sub_area/v_adipose)*venous_pb/v_venous
    J70 <- q_adipose*(sub_area/v_adipose)*arterial_pb/v_arterial-q_adipose*(sub_area/v_adipose)*sc_site_pb/v_sc_site*bpr_pb/kp_adipose_pb
    J71 <- kinf_pb*depot_iv_pb
    J72 <- kinf*depot_iv
    J73 <- 0.1*kurine*kidney
    J74 <- 0.1*kurine_pb*kidney_pb

    # --- ODEs ---
    d/dt(a_metabolized) <- +J23
    d/dt(a_metabolized_pb) <- +J26
    d/dt(a_remainder) <- -J24 +J25
    d/dt(a_remainder_pb) <- -J47 +J48
    d/dt(adipose) <- -J10 +J11
    d/dt(adipose_pb) <- -J41 +J42
    d/dt(arterial) <- -J5 -J6 -J7 -J8 -J11 -J14 -J15 -J16 -J18 -J19 -J22 -J25 +J52 -J59 -J61 -J63
    d/dt(arterial_pb) <- -J28 -J29 -J32 -J34 -J36 -J38 -J40 -J42 -J44 -J46 -J48 -J50 +J54 -J66 -J67 -J70
    d/dt(bone) <- -J13 +J15
    d/dt(bone_pb) <- -J45 +J46
    d/dt(brain) <- -J4 +J5
    d/dt(brain_pb) <- -J37 +J38
    d/dt(depot_iv) <- -J72
    d/dt(depot_iv_pb) <- -J71
    d/dt(gonads) <- +J61 -J62
    d/dt(gonads_pb) <- -J65 +J66
    d/dt(gut) <- +J18 -J20
    d/dt(gut_pb) <- -J31 +J32
    d/dt(heartc_l) <- -J52 +J57
    d/dt(heartc_l_pb) <- -J54 +J58
    d/dt(heartc_r) <- -J51 +J55
    d/dt(heartc_r_pb) <- -J53 +J56
    d/dt(im_site) <- +J63 -J64
    d/dt(im_site_pb) <- +J67 -J68
    d/dt(kidney) <- +J16 -J17 -J73
    d/dt(kidney_pb) <- -J49 +J50 -J74
    d/dt(liver) <- -J1 +J2 +J3 +J20 +J21 +J22
    d/dt(liver_pb) <- -J27 +J28 +J30 +J31 +J33 +J35
    d/dt(lung) <- +J51 -J57
    d/dt(lung_pb) <- +J53 -J58
    d/dt(muscle) <- +J8 -J9
    d/dt(muscle_pb) <- -J39 +J40
    d/dt(pancreas) <- -J2 +J7
    d/dt(pancreas_pb) <- -J35 +J36
    d/dt(sc_site) <- +J59 -J60
    d/dt(sc_site_pb) <- -J69 +J70
    d/dt(skin) <- -J12 +J14
    d/dt(skin_pb) <- -J43 +J44
    d/dt(spleen) <- -J3 +J6
    d/dt(spleen_pb) <- -J33 +J34
    d/dt(stomach) <- +J19 -J21
    d/dt(stomach_pb) <- +J29 -J30
    d/dt(urine) <- +J73
    d/dt(urine_pb) <- +J74
    d/dt(venous) <- +J1 +J4 +J9 +J10 +J12 +J13 +J17 -J23 +J24 -J55 +J60 +J62 +J64 +J72
    d/dt(venous_pb) <- -J26 +J27 +J37 +J39 +J41 +J43 +J45 +J47 +J49 -J56 +J65 +J68 +J69 +J71

    # Observations: plasma concentrations (venous blood / blood:plasma ratio)
    Cc <- venous / v_venous / bpr
    Cc_pb <- venous_pb / v_venous / bpr_pb
  })
}
