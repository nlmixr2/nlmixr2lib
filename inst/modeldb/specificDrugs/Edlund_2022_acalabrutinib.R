Edlund_2022_acalabrutinib <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral acalabrutinib and",
    "its active metabolite ACP-5862 in adults with B-cell malignancies and",
    "healthy subjects (Edlund 2022). Acalabrutinib: two-compartment",
    "disposition with first-order elimination; absorption through a dosing",
    "depot, a chain of five transit compartments and an absorption depot",
    "(mean transit time MTT = (Ntr + 1) / ktr) followed by first-order",
    "absorption, with between-occasion variability on MTT and on the",
    "relative bioavailability F1. ACP-5862: two-compartment disposition",
    "with first-order elimination, formed at 0.4 x CL/F (fraction",
    "metabolised held at 0.4 from the human ADME study). Healthy-subject",
    "status acts on CL/F and Vp/F, ECOG performance status >= 2 on CL/F,",
    "and concomitant proton-pump-inhibitor use on F1. Separate log-scale",
    "residual errors for each analyte on rich-profile versus sparse",
    "sampling occasions."
  )
  reference <- paste(
    "Edlund H, Bellanti F, Liu H, Vishwanathan K, Tomkinson H, Ware J,",
    "Sharma S, Buil-Bruna N. Improved characterization of the",
    "pharmacokinetics of acalabrutinib and its pharmacologically active",
    "metabolite, ACP-5862, in patients with B-cell malignancies and in",
    "healthy subjects using a population pharmacokinetic approach.",
    "Br J Clin Pharmacol. 2022;88(2):846-852. doi:10.1111/bcp.14988"
  )
  vignette <- "Edlund_2022_acalabrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DIS_HEALTHY = list(
      description = "Healthy-subject indicator: 1 = healthy subject, 0 = patient with a B-cell malignancy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with a B-cell malignancy, the most common category)",
      notes = paste(
        "Health status in Edlund 2022. 138 healthy subjects (four phase 1",
        "studies) and 575 patients (eight phase 1-3 studies; Table S1-S2).",
        "Multiplies CL/F by (1 + 0.467) and Vp/F by (1 - 0.556) (Table 1,",
        "categorical form Eq. 2 of the Supplement)."
      ),
      source_name = "health status"
    ),
    ECOG_GE2 = list(
      description = "Baseline ECOG performance status >= 2: 1 = yes, 0 = ECOG 0 or 1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0-1)",
      notes = paste(
        "Table 1 row 'ECOG 2 - CL/F', Figure S10 'ECOG 2+'. 27 of 575",
        "patients had ECOG 2 or 3 (Table S2). ECOG was not recorded for",
        "healthy subjects; set ECOG_GE2 = 0 for them. Multiplies CL/F by",
        "(1 - 0.171)."
      ),
      source_name = "ECOG"
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI)",
      notes = paste(
        "Table 1 row 'PPI - F1'. 66 of 575 patients were PPI users",
        "(including one imputed; Table S2); omeprazole co-administration",
        "in the healthy-subject DDI studies ACE-HV-004 and ACE-HV-112 was",
        "part of the analysed dosing (Table S1). Multiplies F1 by",
        "(1 - 0.358)."
      ),
      source_name = "PPI"
    ),
    OCC = list(
      description = "Dosing-occasion index (1-4) for the between-occasion variability on MTT and F1",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Edlund 2022 Results 3.1 and Figure S2: occasions 1-3 are the",
        "rich-sampling (>= 3 samples) dosing occasions; sparse samples",
        "collected on different dosing days were lumped into one unique",
        "sparse occasion, here OCC = 4. Table 1 footnote b gives BOV",
        "shrinkage as 'the mean of 4 occasions'. Records with any other",
        "OCC value carry no BOV."
      ),
      source_name = "OCC"
    ),
    SAMPLE_INTENSIVE = list(
      description = "Rich-profile occasion indicator: 1 = observation from a rich-sampling occasion, 0 = sparse sample",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sparse sampling)",
      notes = paste(
        "Edlund 2022 Results 3.2: separate residual errors for samples",
        "from occasions with rich profiles vs sparse samples. A rich",
        "profile is >= 3 samples per occasion (Results 3.1). Selects",
        "expSdIntensive / expSdSparse (parent) and their ACP-5862",
        "counterparts."
      ),
      source_name = "rich vs sparse occasion"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the full covariate model (Figure S3) but not retained (95% CI within 0.8-1.25); median-normalised power form per Supplement Eq. 1.",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units = "mL/min/1.73m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the full covariate model (Figure S3) but not retained.",
      source_name = "eGFR"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other races)",
      notes = "Effect on Vc/F was outside the 0.8-1.25 interval but poorly estimated (RSE > 110%) and removed (Results 3.2).",
      source_name = "race"
    ),
    CONMED_H2RA = list(
      description = "Concomitant H2-receptor antagonist use",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no H2RA)",
      notes = "Tested in the full covariate model (Figure S3) but not retained.",
      source_name = "H2RA"
    ),
    HEPIMP = list(
      description = "Hepatic impairment (NCI-ODWG)",
      units = "(category)",
      type = "categorical",
      reference_category = "normal hepatic function",
      notes = "Tested in the full covariate model (Figure S3) but not retained.",
      source_name = "hepatic impairment"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(
      analyte = "acalabrutinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "The 'absorption depot' of Figure S10: receives the last ktr transfer and empties into central at Ka."
    ),
    central = list(analyte = "acalabrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acalabrutinib", units = "mg", specimen = "tissue", verified = TRUE),
    central_acp = list(analyte = "ACP-5862", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_acp = list(analyte = "ACP-5862", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 713L,
    n_studies = 12L,
    age_range = "18-90 years (healthy subjects mean 42.0, patients mean 66.9)",
    weight_range = "39.9-148.6 kg (overall mean 81.0)",
    sex_female_pct = 33.2,
    race_ethnicity = c(White = 89.6, Black = 6.8, Asian = 1.0, Other = 2.6),
    disease_state = paste(
      "138 healthy subjects and 575 adults with B-cell malignancies: chronic",
      "lymphocytic leukaemia (443), Waldenstrom macroglobulinaemia (50),",
      "mantle cell lymphoma (45), diffuse large B-cell lymphoma (15),",
      "multiple myeloma (13) and follicular lymphoma (9)"
    ),
    dose_range = "75-100 mg single doses (healthy subjects); 100-400 mg QD or 100-200 mg BID (patients)",
    n_observations = "8935 acalabrutinib samples from 712 subjects and 2394 ACP-5862 samples from 304 subjects",
    notes = "Edlund 2022 Results 3.1 and Supplement Tables S1-S3. NONMEM 7.3 SAEM with importance sampling; BLQ data handled with M3."
  )

  ini({
    # Acalabrutinib structural parameters (Edlund 2022 Table 1, final model;
    # reference: patient, ECOG 0-1, no PPI).
    lcl <- log(134); label("Apparent clearance of acalabrutinib (L/h)") # Table 1 'CL/F (L/h)' = 134 (RSE 1.97%)
    lvc <- log(31.0); label("Apparent central volume of acalabrutinib (L)") # Table 1 'Vc/F (L)' = 31.0 (RSE 15.4%)
    lq <- log(20.9); label("Apparent intercompartmental clearance of acalabrutinib (L/h)") # Table 1 'Q/F (L/h)' = 20.9 (RSE 3.94%)
    lvp <- log(110); label("Apparent peripheral volume of acalabrutinib (L)") # Table 1 'Vp/F (L)' = 110 (RSE 3.45%)
    lka <- log(1.48); label("First-order absorption rate constant from the absorption depot (1/h)") # Table 1 'Ka (1/h)' = 1.48 (RSE 2.63%)
    lmtt <- log(0.459); label("Mean transit time of the absorption chain (h)") # Table 1 'MTT (h)' = 0.459 (RSE 4.16%)
    ntr <- fixed(5); label("Number of transit compartments (count)") # Results 3.2 '5-transit compartment chain'; Figure S10 'Ntr = 5'
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, reference value (unitless)") # Table 1 lists no F1 typical value; F1 is relative with BOV only (Results 3.2)
    fm <- fixed(0.4); label("Fraction of acalabrutinib clearance forming ACP-5862 (fraction)") # Methods 2.2 'fixed to 0.4 based on ... human absorption, distribution, metabolism and excretion study'

    # ACP-5862 structural parameters (Table 1).
    lcl_acp <- log(21.8); label("Apparent clearance of ACP-5862 (L/h)") # Table 1 'CLM/F (L/h)' = 21.8 (RSE 1.98%)
    lvc_acp <- log(22.7); label("Apparent central volume of ACP-5862 (L)") # Table 1 'VcM/F (L)' = 22.7 (RSE 6.66%)
    lq_acp <- log(26.7); label("Apparent intercompartmental clearance of ACP-5862 (L/h)") # Table 1 'QM/F (L/h)' = 26.7 (RSE 8.65%)
    lvp_acp <- log(89.2); label("Apparent peripheral volume of ACP-5862 (L)") # Table 1 'VpM/F (L)' = 89.2 (RSE 4.58%)

    # Covariate effects (Table 1, footnote c 'Relative change (1 + estimate)';
    # categorical form theta_i = theta_pop * (1 + theta_cov)^cov, Supplement
    # Eq. 2). The PDF drops the minus signs; signs are restored from the
    # Supplement Discussion (healthy CL/F 196 L/h = 134 x 1.467; healthy Vp/F
    # 49 L = 110 x 0.444) and from the reported exposure changes (ECOG >= 2
    # +21% AUC = 1 / 0.829; PPI -36% AUC = 0.642).
    e_dis_healthy_cl <- 0.467; label("Relative change in CL/F for healthy subjects (fraction)") # Table 1 'Healthy subject - CL/F' = 0.467 (RSE 19.3%)
    e_dis_healthy_vp <- -0.556; label("Relative change in Vp/F for healthy subjects (fraction)") # Table 1 'Healthy subject - Vp/F' = -0.556 (RSE 5.59%)
    e_ecog_ge2_cl <- -0.171; label("Relative change in CL/F for ECOG performance status of 2 or more (fraction)") # Table 1 'ECOG 2 - CL/F' = -0.171 (RSE 36.2%)
    e_conmed_ppi_fdepot <- -0.358; label("Relative change in F1 with concomitant proton-pump inhibitor (fraction)") # Table 1 'PPI - F1' = -0.358 (RSE 7.01%)

    # Between-subject variability. Table 1 prints CV%; converted to the
    # log-normal variance omega^2 = log(1 + CV^2).
    etalcl ~ 0.0550978 # Table 1 'BSV CL/F (CV%)' = 23.8; log(1 + 0.238^2)
    etalvc ~ 2.1150500 # Table 1 'BSV Vc/F (CV%)' = 270; log(1 + 2.70^2)
    etalvp ~ 0.1075702 # Table 1 'BSV Vp/F (CV%)' = 33.7; log(1 + 0.337^2)
    etalcl_acp ~ 0.0138280 # Table 1 'BSV CLM/F (CV%)' = 11.8; log(1 + 0.118^2)
    etalvc_acp ~ 0.2011302 # Table 1 'BSV VcM/F (CV%)' = 47.2; log(1 + 0.472^2)
    etalq_acp ~ 0.1532780 # Table 1 'BSV QM/F (CV%)' = 40.7; log(1 + 0.407^2)
    etalvp_acp ~ 0.0362008 # Table 1 'BSV VpM/F (CV%)' = 19.2; log(1 + 0.192^2)

    # Between-occasion variability on MTT and F1, four occasions (three rich
    # occasions plus one lumped sparse occasion; Table 1 footnote b), one
    # reported variance each, encoded as an occasion-indicator expansion.
    # Occasions 2-4 repeat occasion 1's variance (NONMEM BLOCK(1) SAME).
    etaiov_mtt_1 ~ 0.8722970 # Table 1 'BOV MTT (CV%)' = 118; log(1 + 1.18^2)
    etaiov_mtt_2 ~ fixed(0.8722970) # same variance as occasion 1
    etaiov_mtt_3 ~ fixed(0.8722970) # same variance as occasion 1
    etaiov_mtt_4 ~ fixed(0.8722970) # same variance as occasion 1
    etaiov_fdepot_1 ~ 0.2736245 # Table 1 'BOV F1 (CV%)' = 56.1; log(1 + 0.561^2)
    etaiov_fdepot_2 ~ fixed(0.2736245) # same variance as occasion 1
    etaiov_fdepot_3 ~ fixed(0.2736245) # same variance as occasion 1
    etaiov_fdepot_4 ~ fixed(0.2736245) # same variance as occasion 1

    # Residual error: exponential (additive on log-transformed data), SD on
    # the log scale, separate for rich-profile and sparse occasions.
    expSdIntensive <- 0.586; label("Log-scale residual SD, acalabrutinib, rich-profile occasions (log units)") # Table 1 'Residual error (SD)' = 0.586 (RSE 0.558%)
    expSdSparse <- 0.856; label("Log-scale residual SD, acalabrutinib, sparse samples (log units)") # Table 1 'Residual error sparse (SD)' = 0.856 (RSE 1.16%)
    expSdIntensive_acp <- 0.334; label("Log-scale residual SD, ACP-5862, rich-profile occasions (log units)") # Table 1 'Residual error metabolite (SD)' = 0.334 (RSE 1.48%)
    expSdSparse_acp <- 0.234; label("Log-scale residual SD, ACP-5862, sparse samples (log units)") # Table 1 'Residual error metabolite sparse (SD)' = 0.234 (RSE 4.25%)
  })

  model({
    # Molecular weights (g/mol) from the molecular formulae: acalabrutinib
    # C26H23N7O2 and ACP-5862 C26H23N7O3. The model was fitted on the molar
    # scale (concentrations in nM; Supplement Figures S2-S8), so the
    # formation flux is converted from acalabrutinib mass to ACP-5862 mass.
    # Cross-check: Supplement LLOQs 1.0 ng/mL = 2.1 nM and 5 ng/mL = 10.4 nM.
    mw_parent <- 465.52
    mw_acp <- 481.52

    # Occasion indicators and between-occasion variability.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # Individual parameters (categorical covariates per Supplement Eq. 2).
    cl <- exp(lcl + etalcl) * (1 + e_dis_healthy_cl)^DIS_HEALTHY * (1 + e_ecog_ge2_cl)^ECOG_GE2
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (1 + e_dis_healthy_vp)^DIS_HEALTHY
    ka <- exp(lka)
    mtt <- exp(lmtt + iov_mtt)
    ktr <- (ntr + 1) / mtt
    fdepot <- exp(lfdepot + iov_fdepot) * (1 + e_conmed_ppi_fdepot)^CONMED_PPI

    cl_acp <- exp(lcl_acp + etalcl_acp)
    vc_acp <- exp(lvc_acp + etalvc_acp)
    q_acp <- exp(lq_acp + etalq_acp)
    vp_acp <- exp(lvp_acp + etalvp_acp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_acp <- cl_acp / vc_acp
    k12_acp <- q_acp / vc_acp
    k21_acp <- q_acp / vp_acp

    # Dosing depot -> five transit compartments -> absorption depot (transit6)
    # -> central at Ka (Figure S10). Six ktr transfers give MTT = (Ntr + 1) / ktr.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ka * transit6
    d/dt(central) <- ka * transit6 - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # ACP-5862 is formed at Fm x CL/F (molar), the remaining (1 - Fm) x CL/F
    # being other elimination of acalabrutinib.
    d/dt(central_acp) <- fm * kel * central * (mw_acp / mw_parent) -
      (kel_acp + k12_acp) * central_acp + k21_acp * peripheral1_acp
    d/dt(peripheral1_acp) <- k12_acp * central_acp - k21_acp * peripheral1_acp

    f(depot) <- fdepot

    # Concentrations: mg / L x 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc_acp <- 1000 * central_acp / vc_acp

    expSd <- expSdIntensive * SAMPLE_INTENSIVE + expSdSparse * (1 - SAMPLE_INTENSIVE)
    expSd_acp <- expSdIntensive_acp * SAMPLE_INTENSIVE + expSdSparse_acp * (1 - SAMPLE_INTENSIVE)
    Cc ~ lnorm(expSd)
    Cc_acp ~ lnorm(expSd_acp)
  })
}
