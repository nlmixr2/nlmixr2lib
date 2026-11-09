Yao_2025_acalabrutinib <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral acalabrutinib and",
    "its active metabolite ACP-5862 in adults with B-cell malignancies",
    "(including 105 Chinese and 6 Japanese patients) and healthy subjects",
    "(Yao 2025): the Edlund 2022 global model re-estimated with the East",
    "Asian data and an East Asian race effect. Acalabrutinib: two-compartment",
    "disposition with first-order elimination; absorption through a dosing",
    "depot, a chain of five transit compartments and an absorption depot",
    "(mean transit time MTT = (Ntr + 1) / ktr) followed by first-order",
    "absorption, with between-occasion variability on MTT and on the",
    "relative bioavailability F1. ACP-5862: two-compartment disposition",
    "with first-order elimination, formed at 0.4 x CL/F (fraction",
    "metabolised fixed at 0.4). Healthy-subject status acts on CL/F and",
    "Vp/F, ECOG performance status >= 2 on CL/F, concomitant",
    "proton-pump-inhibitor use on F1, and East Asian (Chinese or Japanese)",
    "race on F1 (+22.6%) and Vp/F (-38.1%). Separate log-scale residual",
    "errors for each analyte on rich-profile versus sparse sampling",
    "occasions."
  )
  reference <- paste(
    "Yao T, Ding J, Sharma S, Li Y, Xu Q, Ma P. Assessment of ethnic differences",
    "in pharmacokinetics and clinical responses of acalabrutinib between",
    "Chinese and White patients with B-cell malignancies.",
    "Br J Clin Pharmacol. 2025;91(7):2008-2019. doi:10.1002/bcp.70018"
  )
  vignette <- "Yao_2025_acalabrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    RACE_ASIAN_NORTHEAST = list(
      description = "East Asian race indicator: 1 = Chinese or Japanese patient, 0 = non-East Asian",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-East Asian: White, Black or African American, other)",
      notes = paste(
        "Yao 2025 Table 3 abbreviation 'EASTASIA, Chinese and Japanese",
        "patients' (105 Chinese from NCT03932331 and 6 Japanese from",
        "NCT03198650). Multiplies F1 by (1 + 0.226) and Vp/F by",
        "(1 - 0.381) (Table 3; Results 'Vp/F was 38% lower and F was 23%",
        "higher in East Asian individuals')."
      ),
      source_name = "EASTASIA"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-subject indicator: 1 = healthy subject, 0 = patient with a B-cell malignancy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with a B-cell malignancy)",
      notes = paste(
        "Table 3 'HEALTHGP, healthy individuals'; 138 healthy subjects from",
        "four phase 1 studies. Multiplies CL/F by (1 + 0.645) and Vp/F by",
        "(1 - 0.424)."
      ),
      source_name = "HEALTHGP"
    ),
    ECOG_GE2 = list(
      description = "Baseline ECOG performance status >= 2: 1 = yes, 0 = ECOG 0 or 1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0-1)",
      notes = paste(
        "Table 3 row 'CL/F_ECOG2+'. ECOG is not recorded for healthy",
        "subjects; set ECOG_GE2 = 0 for them. Multiplies CL/F by",
        "(1 - 0.128)."
      ),
      source_name = "ECOG2+"
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI)",
      notes = "Table 3 row 'F1_PPI'. Multiplies F1 by (1 - 0.307). One of the 105 Chinese patients took a PPI (omeprazole; Table 1).",
      source_name = "PPI"
    ),
    OCC = list(
      description = "Dosing-occasion index (1-4) for the between-occasion variability on MTT and F1",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Yao 2025 Table 3 footnote a: BOV shrinkage is 'the mean of the four",
        "occasions'. The model structure is carried over from Edlund 2022,",
        "where occasions 1-3 are the rich-sampling dosing occasions and",
        "sparse samples on different dosing days are lumped into one",
        "unique sparse occasion (here OCC = 4). Records with any other OCC",
        "value carry no BOV."
      ),
      source_name = "OCC"
    ),
    SAMPLE_INTENSIVE = list(
      description = "Rich-profile occasion indicator: 1 = observation from a rich-sampling occasion, 0 = sparse sample",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sparse sampling)",
      notes = paste(
        "Yao 2025 Methods: 'separate residual error parameters were included",
        "for samples obtained on occasions for which rich profiles vs.",
        "sparse samples had been collected'. Selects expSdIntensive /",
        "expSdSparse (parent) and their ACP-5862 counterparts."
      ),
      source_name = "rich vs sparse occasion"
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
      notes = "The absorption depot of Figure 1: receives the last ktr transfer and empties into central at Ka."
    ),
    central = list(analyte = "acalabrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acalabrutinib", units = "mg", specimen = "tissue", verified = TRUE),
    central_acp = list(analyte = "ACP-5862", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_acp = list(analyte = "ACP-5862", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 824L,
    n_studies = 14L,
    age_range = "Chinese patients mean 60.5 (SD 10.1) years; whole-dataset range not reported",
    weight_range = "Chinese patients mean 65.9 (SD 10.2) kg; whole-dataset range not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = c(White = 75.5, Black = 5.7, Chinese = 12.7, Japanese = 0.7, Other = 5.3),
    disease_state = paste(
      "138 healthy subjects and 686 adults with B-cell malignancies (chronic",
      "lymphocytic leukaemia, mantle cell lymphoma, follicular lymphoma,",
      "diffuse large B-cell lymphoma, small lymphocytic lymphoma, multiple",
      "myeloma, Waldenstrom macroglobulinaemia). Chinese cohort (n = 105):",
      "CLL 64, MCL 34, FL 4, SLL 2, DLBCL 1; 28.6% female; ECOG 0/1/2 =",
      "66/35/4; normal / mild / moderate hepatic function 83/19/3; mean",
      "creatinine clearance 74.7 mL/min."
    ),
    dose_range = "Global studies 100-400 mg daily (100 mg QD to 200 mg BID; single doses in healthy subjects); China and Japan studies 100 mg BID",
    regions = "Global (14 trials) plus China (NCT03932331) and Japan (NCT03198650)",
    n_observations = "11152 acalabrutinib and 3367 ACP-5862 plasma concentrations; East Asian subset 972 samples from 111 patients",
    notes = paste(
      "Yao 2025 Results 'Dataset', Tables 1-2 (Chinese cohort only) and the",
      "Figure 3 caption race counts (Black 47, Chinese 105, Japanese 6,",
      "other 44, White 622). NONMEM 7.3 SAEM with importance sampling;",
      "BLQ data (4.2% for acalabrutinib in Chinese patients) handled with",
      "M3. Run 525bCLL3."
    )
  )

  ini({
    # Acalabrutinib structural parameters (Yao 2025 Table 3, final model run
    # 525bCLL3; reference: patient, non-East Asian, ECOG 0-1, no PPI).
    lcl <- log(130.3); label("Apparent clearance of acalabrutinib (L/h)") # Table 3 'CL/F (L/h)' = 130.3 (RSE 0.4%)
    lvc <- log(28.8); label("Apparent central volume of acalabrutinib (L)") # Table 3 'Vc/F (L)' = 28.8 (RSE 4.5%)
    lq <- log(19.5); label("Apparent intercompartmental clearance of acalabrutinib (L/h)") # Table 3 'Q/F (L/h)' = 19.5 (RSE 1.3%)
    lvp <- log(105.6); label("Apparent peripheral volume of acalabrutinib (L)") # Table 3 'Vp/F (L)' = 105.6 (RSE 0.7%)
    lka <- log(1.4); label("First-order absorption rate constant from the absorption depot (1/h)") # Table 3 'Ka (h-1)' = 1.4 (RSE 6.6%)
    lmtt <- log(0.4); label("Mean transit time of the absorption chain (h)") # Table 3 'MTT (h)' = 0.4 (RSE 5%)
    ntr <- fixed(5); label("Number of transit compartments (count)") # Methods 'two-compartment model with five transit compartment'; Figure 1
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, reference value (unitless)") # Table 3 lists no F1 typical value; F1 is relative with BOV and covariates only
    fm <- fixed(0.4); label("Fraction of acalabrutinib clearance forming ACP-5862 (fraction)") # Methods 'the fraction of the drug metabolized having been fixed to 0.4'

    # ACP-5862 structural parameters (Table 3).
    lcl_acp <- log(21.5); label("Apparent clearance of ACP-5862 (L/h)") # Table 3 'CLM/F (L/h)' = 21.5 (RSE 0.6%)
    lvc_acp <- log(24.3); label("Apparent central volume of ACP-5862 (L)") # Table 3 'VcM/F (L)' = 24.3 (RSE 1.7%)
    lq_acp <- log(21.1); label("Apparent intercompartmental clearance of ACP-5862 (L/h)") # Table 3 'QM/F (L/h)' = 21.1 (RSE 2.9%)
    lvp_acp <- log(90.9); label("Apparent peripheral volume of ACP-5862 (L)") # Table 3 'VpM/F (L)' = 90.9 (RSE 1.3%)

    # Covariate effects: proportional shift P = TVP x (1 + Pcov) for a
    # subject in the covariate category (Methods equation; Table 3 footnote
    # 'Relative change (1 + estimate)').
    e_dis_healthy_cl <- 0.645; label("Relative change in CL/F for healthy subjects (fraction)") # Table 3 'CL/F_HEALTHGP' = 0.645 (RSE 17.1%)
    e_dis_healthy_vp <- -0.424; label("Relative change in Vp/F for healthy subjects (fraction)") # Table 3 'Vp/F_HEALTHGP' = -0.424 (RSE 5.9%)
    e_ecog_ge2_cl <- -0.128; label("Relative change in CL/F for ECOG performance status of 2 or more (fraction)") # Table 3 'CL/F_ECOG2+' = -0.128 (RSE 43.9%)
    e_conmed_ppi_fdepot <- -0.307; label("Relative change in F1 with concomitant proton-pump inhibitor (fraction)") # Table 3 'F1_PPI' = -0.307 (RSE 6.3%)
    e_race_asian_northeast_vp <- -0.381; label("Relative change in Vp/F for East Asian (Chinese or Japanese) patients (fraction)") # Table 3 'Vp/F_EASTASIA' = -0.381 (RSE 15%)
    e_race_asian_northeast_fdepot <- 0.226; label("Relative change in F1 for East Asian (Chinese or Japanese) patients (fraction)") # Table 3 'F1_EASTASIA' = 0.226 (RSE 44.3%)

    # Between-subject variability. Table 3 prints the BSV and BOV magnitudes
    # as CV% (mislabelled with the parameter units); converted to the
    # log-normal variance omega^2 = log(1 + CV^2).
    etalcl ~ 0.0524255 # Table 3 'BSV CL/F' = 23.2; log(1 + 0.232^2)
    etalvc ~ 2.3382997 # Table 3 'BSV Vc/F' = 306; log(1 + 3.06^2)
    etalvp ~ 0.1212269 # Table 3 'BSV Vp/F' = 35.9; log(1 + 0.359^2)
    etalcl_acp ~ 0.0172739 # Table 3 'BSV CLM/F' = 13.2; log(1 + 0.132^2)
    etalvc_acp ~ 0.1309195 # Table 3 'BSV VcM/F' = 37.4; log(1 + 0.374^2)
    etalq_acp ~ 0.3288547 # Table 3 'BSV QM/F' = 62.4; log(1 + 0.624^2)
    etalvp_acp ~ 0.0823622 # Table 3 'BSV VpM/F' = 29.3; log(1 + 0.293^2)

    # Between-occasion variability on MTT and F1, four occasions (Table 3
    # footnote a), one reported variance each, encoded as an
    # occasion-indicator expansion. Occasions 2-4 repeat occasion 1's
    # variance (NONMEM BLOCK(1) SAME).
    etaiov_mtt_1 ~ 0.9409833 # Table 3 'BOV MTT' = 125; log(1 + 1.25^2)
    etaiov_mtt_2 ~ fixed(0.9409833) # same variance as occasion 1
    etaiov_mtt_3 ~ fixed(0.9409833) # same variance as occasion 1
    etaiov_mtt_4 ~ fixed(0.9409833) # same variance as occasion 1
    etaiov_fdepot_1 ~ 0.2467359 # Table 3 'BOV F1' = 52.9; log(1 + 0.529^2)
    etaiov_fdepot_2 ~ fixed(0.2467359) # same variance as occasion 1
    etaiov_fdepot_3 ~ fixed(0.2467359) # same variance as occasion 1
    etaiov_fdepot_4 ~ fixed(0.2467359) # same variance as occasion 1

    # Residual error: additive on log-transformed concentrations (Methods),
    # SD on the log scale, separate for rich-profile and sparse occasions.
    expSdIntensive <- 0.641; label("Log-scale residual SD, acalabrutinib, rich-profile occasions (log units)") # Table 3 'RE parent (SD)' = 0.641 (RSE 0.5%)
    expSdSparse <- 0.925; label("Log-scale residual SD, acalabrutinib, sparse samples (log units)") # Table 3 'RE parent sparse samples (SD)' = 0.925 (RSE 0.9%)
    expSdIntensive_acp <- 0.36; label("Log-scale residual SD, ACP-5862, rich-profile occasions (log units)") # Table 3 'RE metabolite (SD)' = 0.36 (RSE 1.3%)
    expSdSparse_acp <- 0.401; label("Log-scale residual SD, ACP-5862, sparse samples (log units)") # Table 3 'RE metabolite sparse samples (SD)' = 0.401 (RSE 1.5%)
  })

  model({
    # Molecular weights (g/mol): acalabrutinib 465.5 (Yao 2025 Methods,
    # 'MW acalabrutinib = 465.5'; formula C26H23N7O2) and ACP-5862 from its
    # formula C26H23N7O3. Concentrations were modelled in molar units
    # (Methods, E-R analyses), so the formation flux is converted from
    # acalabrutinib mass to ACP-5862 mass.
    mw_parent <- 465.52
    mw_acp <- 481.52

    # Occasion indicators and between-occasion variability.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # Individual parameters (proportional categorical covariates, Methods).
    cl <- exp(lcl + etalcl) * (1 + e_dis_healthy_cl)^DIS_HEALTHY * (1 + e_ecog_ge2_cl)^ECOG_GE2
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (1 + e_dis_healthy_vp)^DIS_HEALTHY *
      (1 + e_race_asian_northeast_vp)^RACE_ASIAN_NORTHEAST
    ka <- exp(lka)
    mtt <- exp(lmtt + iov_mtt)
    ktr <- (ntr + 1) / mtt
    fdepot <- exp(lfdepot + iov_fdepot) * (1 + e_conmed_ppi_fdepot)^CONMED_PPI *
      (1 + e_race_asian_northeast_fdepot)^RACE_ASIAN_NORTHEAST

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
    # -> central at Ka (Figure 1). Six ktr transfers give MTT = (Ntr + 1) / ktr.
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
