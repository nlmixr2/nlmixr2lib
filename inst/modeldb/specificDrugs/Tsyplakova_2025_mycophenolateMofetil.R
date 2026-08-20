Tsyplakova_2025_mycophenolateMofetil <- function() {
  description <- "One-compartment population PK model with first-order absorption for mycophenolic acid after oral mycophenolate mofetil (MMF) in adult renal transplant recipients, with post-transplant time and total daily dose on apparent clearance and inter-occasion variability across monthly visits"
  reference <- paste(
    "Tsyplakova A, Catic-Djordjevic A, Stefanovic N, Karalis VD (2025)",
    "Optimizing Mycophenolate Therapy in Renal Transplant Patients Using",
    "Machine Learning and Population Pharmacokinetic Modeling.",
    "Med Sci (Basel) 13(4):235. doi:10.3390/medsci13040235.",
    sep = " "
  )
  vignette <- "Tsyplakova_2025_mycophenolate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    POD = list(
      description        = "Post-transplant time (days elapsed since renal transplantation)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tsyplakova 2025 reports this covariate in MONTHS ('Diff_in_months' in Eq. 11; 'PTP',",
        "post-transplant time). The canonical POD column is in DAYS, so model() divides by",
        "30.4375 days/month before forming the paper's ratio. The MMF normalisation constant is",
        "21 months (Eq. 11) -- materially lower than the 67 months used by the EC-MPS model",
        "(Eq. 8), consistent with the Eq. 1 definition of Cmean as the per-subgroup covariate",
        "mean and with the 13-patient MMF subgroup being less far post-transplant than the",
        "63-patient EC-MPS subgroup; the whole-cohort median in Table 1 is 70 months (IQR 84.3).",
        "All patients were at least 3 months post-transplant with stable graft function (Section",
        "2.1). No upper cap was applied by the authors. The effect is POSITIVE (exponent 0.33,",
        "roughly double the EC-MPS estimate): apparent clearance rises with time since",
        "transplantation.",
        sep = " "
      ),
      source_name        = "Diff_in_months / PTP"
    ),
    DOSE_MPA_MGD = list(
      description        = "Total daily dose of mycophenolic acid, on the MPA-equivalent mass scale",
      units              = "mg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tsyplakova 2025 Eq. 11 'TDD', normalised to 1500 mg/day -- the SAME constant the EC-MPS",
        "model uses (Eq. 8), which is why the constant is read as an MPA-equivalent cohort mean",
        "rather than a per-formulation one. MMF doses in this study were 500-2000 mg/day given",
        "twice daily and were converted onto the common MPA-equivalent scale with a factor of",
        "0.72 (Section 2.1: MMF 500-2000 mg/day maps onto EC-MPS 360-1440 mg/day). Populate this",
        "column on the same converted scale. The positive exponent (1.27) is",
        "indication-by-confounding in a therapeutically-monitored stable cohort -- the authors",
        "state higher-dose patients likely required those doses to compensate for greater",
        "clearance, and that TDD may partly proxy for body weight, which was not recorded.",
        sep = " "
      ),
      source_name        = "TDD"
    ),
    OCC = list(
      description        = "Occasion index; each monthly follow-up visit is a separate occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Tsyplakova 2025 Section 2.3: 'Inter-occasion variability (IOV) was also incorporated,",
        "treating each subsequent visit as a separate occasion.' Section 2.2: patients were",
        "monitored for six months with MPA measured monthly, so occasions run 1..6; the paper does",
        "not state the realised maximum number of occasions per patient. Decomposed inside model()",
        "into binary indicators oc1..oc6 that multiplex the per-occasion IOV etas on log-ka,",
        "log-Vc and log-CL. Monolix estimates a single IOV standard deviation per parameter shared",
        "across occasions (Table 3 gamma), so occasions 2-6 are fixed to the occasion-1 variance.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 13,
    n_studies      = 1,
    age_median     = "51 years (IQR 14; whole 76-patient cohort)",
    sex_female_pct = 34.2,
    disease_state  = "adult renal transplant recipients, at least 3 months post-transplant with stable graft function, on mycophenolate + tacrolimus + low-dose prednisone",
    dose_range     = "500-2000 mg/day mycophenolate mofetil, given twice daily (360-1440 mg/day MPA-equivalent after the 0.72 conversion)",
    regions        = "Serbia (University Clinical Centre of Nis)",
    notes          = paste(
      "Tsyplakova 2025 Table 1: 76 patients total contributed 209 MPA plasma samples and 65 saliva",
      "samples; 13 (17.1%) received MMF and 63 (82.9%) received EC-MPS, so this MMF model is based",
      "on the 13-patient subgroup and is the smaller and less precisely estimated of the two",
      "models the paper reports (the authors selected the EC-MPS model for their Monte Carlo",
      "dosing simulations 'due to the larger proportion of data derived from this formulation',",
      "Section 3.3). Median post-transplantation time 70 months (IQR 84.3) across the whole",
      "cohort; 54 (71%) living-donor and 18 (24%) deceased-donor grafts. Laboratory medians: urea",
      "7.8 mmol/L (IQR 5.4), creatinine 136 umol/L (IQR 60), WBC 7.9 x10^9/L, RBC 4.7 x10^12/L,",
      "haemoglobin 138 g/L, haematocrit 41.4%, platelets 225 x10^9/L. Body weight and serum",
      "albumin were NOT recorded, which is why the model carries no allometric or protein-binding",
      "term (Section 4). Only steady-state trough (C0) concentrations were available -- one",
      "measurement per occasion -- which is why the authors selected a one-compartment structure",
      "and estimated ka and V by maximum a posteriori estimation while CL was estimated by maximum",
      "likelihood (Sections 2.3 and 3.1.2).",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Tsyplakova 2025 Table 3 and Eqs. 9-11. Only oral data
    # were available, so V and CL are apparent (V/F and CL/F); the authors report
    # them as V and Cl and did not estimate bioavailability separately.
    lka <- log(0.23);    label("Absorption rate constant (1/h)")            # Table 3 kapop = 0.23 (RSE 22.6%); Eq. 9
    lvc <- log(196.43);  label("Apparent central volume of distribution V/F (L)")  # Table 3 Vpop = 196.43 L (RSE 28.9%); Eq. 10
    lcl <- log(9.3);     label("Apparent clearance CL/F at the reference covariate values (L/h)")  # Table 3 Clpop = 9.3 L/h (RSE 11.2%); Eq. 11

    # Covariate effects on apparent clearance. Eq. 1 gives the general continuous
    # form log(Pi) = log(Ppop) + betaC * log(Ci / Cmean) + eta + kappa, so each
    # coefficient is a power exponent on the covariate normalised to Cmean.
    e_pod_cl          <- 0.33; label("Power exponent on (post-transplant time / 21 months) for CL/F (unitless)")  # Table 3 beta(PTP) = 0.33 (RSE 21.8%); Eq. 11
    e_dose_mpa_mgd_cl <- 1.27; label("Power exponent on (total daily MPA dose / 1500 mg/day) for CL/F (unitless)")  # Table 3 beta(TDD) = 1.27 (RSE 23.1%); Eq. 11

    # Inter-individual variability. Table 3 reports the "Standard Deviation of the
    # Random Effects", so each omega below is SD^2 on the log scale.
    etalka ~ 0.27^2;  # Table 3 omega(ka) = 0.27 (RSE 26.5%)
    etalvc ~ 0.09^2;  # Table 3 omega(V)  = 0.09 (RSE 32.4%)
    etalcl ~ 0.32^2;  # Table 3 omega(Cl) = 0.32 (RSE 19.2%)

    # Inter-occasion variability across the six monthly visits. Table 3 reports a
    # single gamma per parameter shared across occasions (Monolix estimates one IOV
    # SD per parameter), so occasions 2-6 are fixed to the occasion-1 variance.
    etaiov_ka_1 ~ 0.48^2;         # Table 3 gamma(ka) = 0.48 (RSE 39.8%)
    etaiov_ka_2 ~ fix(0.48^2);    # shared IOV variance across occasions
    etaiov_ka_3 ~ fix(0.48^2);    # shared IOV variance across occasions
    etaiov_ka_4 ~ fix(0.48^2);    # shared IOV variance across occasions
    etaiov_ka_5 ~ fix(0.48^2);    # shared IOV variance across occasions
    etaiov_ka_6 ~ fix(0.48^2);    # shared IOV variance across occasions

    etaiov_vc_1 ~ 0.33^2;         # Table 3 gamma(V) = 0.33 (RSE 21.4%)
    etaiov_vc_2 ~ fix(0.33^2);    # shared IOV variance across occasions
    etaiov_vc_3 ~ fix(0.33^2);    # shared IOV variance across occasions
    etaiov_vc_4 ~ fix(0.33^2);    # shared IOV variance across occasions
    etaiov_vc_5 ~ fix(0.33^2);    # shared IOV variance across occasions
    etaiov_vc_6 ~ fix(0.33^2);    # shared IOV variance across occasions

    etaiov_cl_1 ~ 0.27^2;         # Table 3 gamma(Cl) = 0.27 (RSE 22.3%)
    etaiov_cl_2 ~ fix(0.27^2);    # shared IOV variance across occasions
    etaiov_cl_3 ~ fix(0.27^2);    # shared IOV variance across occasions
    etaiov_cl_4 ~ fix(0.27^2);    # shared IOV variance across occasions
    etaiov_cl_5 ~ fix(0.27^2);    # shared IOV variance across occasions
    etaiov_cl_6 ~ fix(0.27^2);    # shared IOV variance across occasions

    # Proportional-only residual error (Eq. 4), the structure selected for the MMF
    # model. Table 3 reports no additive component 'a' for this model.
    propSd <- 0.17; label("Proportional residual error (fraction)")  # Table 3 'b' = 0.17 (RSE 37.5%)
  })

  model({
    # Occasion indicators for IOV. OCC is the integer occasion column; each
    # monthly follow-up visit is a separate occasion (Section 2.3).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5 + oc6 * etaiov_vc_6
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # Post-transplant time on the paper's scale. POD is supplied in days (the
    # canonical unit); Eq. 11 normalises post-transplant time in MONTHS to 21.
    podMonth <- POD / 30.4375

    # Individual parameters. Eqs. 9-11.
    ka <- exp(lka + etalka + iov_ka)
    vc <- exp(lvc + etalvc + iov_vc)
    cl <- exp(lcl + etalcl + iov_cl) *
      (podMonth / 21)^e_pod_cl *
      (DOSE_MPA_MGD / 1500)^e_dose_mpa_mgd_cl

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
