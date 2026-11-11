Lindauer_2025_risperidone_oral <- function() {
  description <- paste(
    "Population PK model of oral risperidone and 9-hydroxyrisperidone that",
    "Lindauer 2025 used as the oral comparator for its dopamine D2 receptor",
    "occupancy simulations. Two-compartment disposition for each analyte with",
    "a lagged, sequential zero- then first-order absorption from one depot;",
    "a fraction FP of the absorbed dose reaches the circulation as",
    "9-hydroxyrisperidone (first pass) and the rest as risperidone, which is",
    "converted systemically through a formation clearance. FP and the",
    "formation clearance depend on CYP2D6 metabolizer phenotype (a",
    "three-subpopulation mixture model in the source). Most parameters are",
    "fixed to Vermeulen 2007; the shared central volume was re-estimated on",
    "data on file. D2 receptor occupancy follows the Uchida 2011 Emax model",
    "of the active moiety (risperidone + 9-hydroxyrisperidone)."
  )
  reference <- paste(
    "Lindauer A, Snoeck E, Laveille C, Ayani I, Ochoa Diaz de Monasterioguren L,",
    "Almendros M, Martinez-Gonzalez J, Anta L, Gutierro I. Exposure-Efficacy",
    "Analysis and Dopamine D2 Receptor Occupancy in Adults with Schizophrenia",
    "after Treatment with the Monthly Intramuscular Injectable Risperidone ISM.",
    "J Clin Pharmacol. 2025;65(3):350-360. doi:10.1002/jcph.6152 (Supplemental",
    "Appendix 2). Absorption and disposition structure from Vermeulen A,",
    "Piotrovsky V, Ludwig EA. J Pharmacokinet Pharmacodyn. 2007;34(2):183-206",
    "(doi:10.1007/s10928-006-9040-2); D2 receptor occupancy model from Uchida H",
    "et al. J Clin Psychopharmacol. 2011;31(3):318-325",
    "(doi:10.1097/JCP.0b013e318218d339)."
  )
  vignette <- "Lindauer_2025_risperidone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CYP2D6_PM = list(
      description = "CYP2D6 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer when CYP2D6_IM is also 0)",
      notes = paste(
        "Selects the poor-metabolizer typical values of FP (0.0175) and of",
        "the risperidone-to-9-hydroxyrisperidone formation clearance CLPM",
        "(1.18 L/h). In the source the phenotype is a latent mixture-model",
        "subpopulation with prior probability P(PM) = 4.93%; supply the",
        "genotype-derived phenotype, or draw it with that probability."
      ),
      source_name = "mixture subpopulation PM (P(PM) = 4.93%)"
    ),
    CYP2D6_IM = list(
      description = "CYP2D6 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer when CYP2D6_PM is also 0)",
      notes = paste(
        "Selects the intermediate-metabolizer typical values of FP (0.109)",
        "and CLPM (4.37 L/h). Mixture prior probability P(IM) = 8.26%; the",
        "extensive-metabolizer reference (FP 0.413, CLPM 19.6 L/h) takes the",
        "remaining 86.81%."
      ),
      source_name = "mixture subpopulation IM (P(IM) = 8.26%)"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on F1, D1 and KA: a new draw of etaiov_fdepot, etaiov_d1 and etaiov_ka for every value.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Appendix S2 reports IOV on F1 (46.5% CV), D1 (158% CV) and KA",
        "(108% CV). The vignette starts a new occasion at each once-daily",
        "dose (OCC = dose number). Supply an integer per record; a",
        "single-occasion simulation can set OCC = 1 throughout."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "risperidone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "risperidone", units = "mg", specimen = "tissue", verified = TRUE),
    central_9oh = list(analyte = "9-hydroxyrisperidone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_9oh = list(analyte = "9-hydroxyrisperidone", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    disease_state = "Patients with schizophrenia or bipolar disorder receiving oral risperidone (comparator model; the fitted population is described in Vermeulen 2007 and in data on file)",
    dose_range = "Simulated at 3 and 4 mg once daily",
    cyp2d6_mixture = "P(PM) = 4.93%, P(IM) = 8.26%, P(EM) = 86.81% (Appendix S2 table)",
    notes = paste(
      "The source reports the parameter table and model diagram (Appendix S2)",
      "but no demographics for the population in which the shared central",
      "volume V2 = V4 was re-estimated (ROVI data on file). All other",
      "parameters, IIV and IOV are fixed to Vermeulen 2007."
    )
  )

  ini({
    # Absorption: dose enters the depot as a zero-order input of duration D1
    # after a lag ALAG1, then leaves it first-order with rate KA
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F1 (unitless)") # Appendix S2 'F1' = 1 fix
    ltlag <- fixed(log(0.165)); label("Absorption lag time ALAG1 (h)") # Appendix S2 'ALAG1 (h)' = 0.165 fix
    ld1 <- fixed(log(0.458)); label("Duration of the zero-order input into the depot D1 (h)") # Appendix S2 'D1 (h)' = 0.458 fix
    lka <- fixed(log(2.34)); label("First-order absorption rate constant KA out of the depot (1/h)") # Appendix S2 'KA (1/h)' = 2.34 fix

    # Fraction of the absorbed dose reaching the circulation as
    # 9-hydroxyrisperidone (FP), by CYP2D6 phenotype, on the logit scale
    logitfm_9oh_em <- fixed(log(0.413 / (1 - 0.413))); label("Logit of FP in CYP2D6 extensive metabolizers (unitless)") # Appendix S2 'FP (EM)' = 0.413 fix
    logitfm_9oh_im <- fixed(log(0.109 / (1 - 0.109))); label("Logit of FP in CYP2D6 intermediate metabolizers (unitless)") # Appendix S2 'FP (IM)' = 0.109 fix
    logitfm_9oh_pm <- fixed(log(0.0175 / (1 - 0.0175))); label("Logit of FP in CYP2D6 poor metabolizers (unitless)") # Appendix S2 'FP (PM)' = 0.0175 fix

    # Risperidone disposition
    lcl <- fixed(log(2.84)); label("Risperidone elimination clearance other than conversion to 9-hydroxyrisperidone CLP (L/h)") # Appendix S2 'CLP (L/h)' = 2.84 fix
    lvc <- log(118); label("Central volume of risperidone V2, shared with 9-hydroxyrisperidone V4 (L)") # Appendix S2 'V2=V4 (L)' = 118 (RSE 4%)
    lvp <- fixed(log(100)); label("Peripheral volume of risperidone V3 (L)") # Appendix S2 'V3 (L)' = 100 fix
    lq <- fixed(log(3.65)); label("Intercompartmental clearance of risperidone QP (L/h)") # Appendix S2 'QP (L/h)' = 3.65 fix

    # Systemic conversion of risperidone to 9-hydroxyrisperidone, by phenotype
    lcl_form_9oh_em <- fixed(log(19.6)); label("Risperidone-to-9-hydroxyrisperidone formation clearance CLPM, extensive metabolizers (L/h)") # Appendix S2 'CLPM (EM) (L/h)' = 19.6 fix
    lcl_form_9oh_im <- fixed(log(4.37)); label("Risperidone-to-9-hydroxyrisperidone formation clearance CLPM, intermediate metabolizers (L/h)") # Appendix S2 'CLPM (IM) (L/h)' = 4.37 fix
    lcl_form_9oh_pm <- fixed(log(1.18)); label("Risperidone-to-9-hydroxyrisperidone formation clearance CLPM, poor metabolizers (L/h)") # Appendix S2 'CLPM (PM) (L/h)' = 1.18 fix

    # 9-hydroxyrisperidone disposition (central volume shared with risperidone)
    lcl_9oh <- fixed(log(5.99)); label("9-hydroxyrisperidone elimination clearance CLM (L/h)") # Appendix S2 'CLM (L/h)' = 5.99 fix
    lvp_9oh <- fixed(log(91.8)); label("Peripheral volume of 9-hydroxyrisperidone V5 (L)") # Appendix S2 'V5 (L)' = 91.8 fix
    lq_9oh <- fixed(log(1.67)); label("Intercompartmental clearance of 9-hydroxyrisperidone QM (L/h)") # Appendix S2 'QM (L/h)' = 1.67 fix

    # Uchida 2011 D2 receptor occupancy Emax model of the active moiety
    emax <- fixed(100); label("Maximal dopamine D2 receptor occupancy (%)") # Methods 'Simulations of D2 Receptor Occupancy': Emax fixed to 100%
    lec50 <- fixed(log(8.2)); label("Active-moiety concentration giving 50% D2 receptor occupancy (ng/mL)") # Methods 'Simulations of D2 Receptor Occupancy': EC50 = 8.2 ng/mL

    # IIV, all fixed; variance = log(1 + CV^2) for log-normal parameters and
    # SD^2 for the logit-domain FP
    etaltlag ~ fixed(log(1 + 0.41^2)) # Appendix S2 'ALAG1' IIV 41% CV
    etald1 ~ fixed(log(1 + 1.13^2)) # Appendix S2 'D1' IIV 113% CV
    etalka ~ fixed(log(1 + 1.49^2)) # Appendix S2 'KA' IIV 149% CV
    etalogitfm_9oh ~ fixed(1.17^2) # Appendix S2 'FP' IIV 117 (footnote a: %SD in the logit domain), shared by the three phenotypes
    etalcl ~ fixed(log(1 + 1.84^2)) # Appendix S2 'CLP' IIV 184% CV
    etalvc ~ fixed(log(1 + 0.30^2)) # Appendix S2 'V2=V4' IIV 30% CV
    etalvp ~ fixed(log(1 + 0.54^2)) # Appendix S2 'V3' IIV 54% CV
    etalq ~ fixed(log(1 + 2.15^2)) # Appendix S2 'QP' IIV 215% CV
    etalcl_form_9oh ~ fixed(log(1 + 0.33^2)) # Appendix S2 'CLPM' IIV 33% CV, shared by the three phenotypes
    etalcl_9oh ~ fixed(log(1 + 0.20^2)) # Appendix S2 'CLM' IIV 20% CV
    etalvp_9oh ~ fixed(log(1 + 0.81^2)) # Appendix S2 'V5' IIV 81% CV

    # IOV, all fixed; a new draw for every value of the OCC column
    etaiov_fdepot ~ fixed(log(1 + 0.465^2)) | OCC # Appendix S2 'F1' IOV 46.5% CV
    etaiov_d1 ~ fixed(log(1 + 1.58^2)) | OCC # Appendix S2 'D1' IOV 158% CV
    etaiov_ka ~ fixed(log(1 + 1.08^2)) | OCC # Appendix S2 'KA' IOV 108% CV

    # Residual error is not reported (the source simulations exclude it)
    propSd <- fixed(0); label("Proportional residual error on risperidone (fraction; not reported)") # not reported in Appendix S2
    propSd_9oh <- fixed(0); label("Proportional residual error on 9-hydroxyrisperidone (fraction; not reported)") # not reported in Appendix S2
  })

  model({
    em <- 1 - CYP2D6_PM - CYP2D6_IM

    fdepot <- exp(lfdepot + etaiov_fdepot)
    tlag <- exp(ltlag + etaltlag)
    d1 <- exp(ld1 + etald1 + etaiov_d1)
    ka <- exp(lka + etalka + etaiov_ka)

    logitfm_9oh_ind <- em * logitfm_9oh_em + CYP2D6_IM * logitfm_9oh_im + CYP2D6_PM * logitfm_9oh_pm + etalogitfm_9oh
    fm_9oh <- expit(logitfm_9oh_ind)

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    cl_form_9oh <- exp(em * lcl_form_9oh_em + CYP2D6_IM * lcl_form_9oh_im + CYP2D6_PM * lcl_form_9oh_pm + etalcl_form_9oh)
    vc_9oh <- vc
    cl_9oh <- exp(lcl_9oh + etalcl_9oh)
    vp_9oh <- exp(lvp_9oh + etalvp_9oh)
    q_9oh <- exp(lq_9oh)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * (1 - fm_9oh) * depot - (cl + cl_form_9oh + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(central_9oh) <- ka * fm_9oh * depot + cl_form_9oh / vc * central - (cl_9oh + q_9oh) / vc_9oh * central_9oh + q_9oh / vp_9oh * peripheral1_9oh
    d/dt(peripheral1_9oh) <- q_9oh / vc_9oh * central_9oh - q_9oh / vp_9oh * peripheral1_9oh

    f(depot) <- fdepot
    alag(depot) <- tlag
    dur(depot) <- d1

    # Amounts in mg; mg/L * 1000 = ng/mL. The active moiety is the sum of the
    # two analytes, as in the source's D2 receptor occupancy simulations.
    Cc <- 1000 * central / vc
    Cc_9oh <- 1000 * central_9oh / vc_9oh
    AM <- Cc + Cc_9oh
    D2RO <- emax * AM / (exp(lec50) + AM)

    Cc ~ prop(propSd)
    Cc_9oh ~ prop(propSd_9oh)
  })
}
