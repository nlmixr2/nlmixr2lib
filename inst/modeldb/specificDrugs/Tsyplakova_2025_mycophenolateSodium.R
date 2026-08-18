Tsyplakova_2025_mycophenolateSodium <- function() {
  description <- "One-compartment population PK model with first-order absorption for mycophenolic acid after oral enteric-coated mycophenolate sodium (EC-MPS) in adult renal transplant recipients, with post-transplant time and total daily dose on apparent clearance and inter-occasion variability across monthly visits"
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
        "Tsyplakova 2025 reports this covariate in MONTHS ('Diff_in_months' in Eq. 8; 'PTP',",
        "post-transplant time). The canonical POD column is in DAYS, so model() divides by",
        "30.4375 days/month before forming the paper's ratio. The EC-MPS normalisation constant",
        "is 67 months (Eq. 8), the EC-MPS-subgroup mean per the Eq. 1 definition of Cmean; the",
        "whole-cohort median in Table 1 is 70 months (IQR 84.3). All patients were at least 3",
        "months post-transplant with stable graft function (Section 2.1 inclusion criteria), so",
        "the power term is only supported in the late post-transplant window. No upper cap was",
        "applied by the authors. The effect is POSITIVE (exponent 0.16): apparent clearance rises",
        "with time since transplantation, which Tsyplakova 2025 Discussion attributes to",
        "progressively improving renal function, normalisation of albumin binding and recovery of",
        "hepatic glucuronidation.",
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
        "Tsyplakova 2025 Eq. 8 'TDD', normalised to 1500 mg/day. The two formulations were placed",
        "on a common MPA-equivalent scale by converting MMF doses with a factor of 0.72 (Section",
        "2.1), and both the EC-MPS and the MMF model use the same 1500 mg/day constant, which is",
        "why the constant is read as an MPA-equivalent cohort mean rather than a per-formulation",
        "one. Protocol EC-MPS doses in this study were 360-1440 mg/day given twice daily; the",
        "Discussion notes some participants exceeded the 2160 mg/day maximum therapeutic dose, so",
        "the observed TDD range extends above the protocol range. The positive exponent (0.77) is",
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
        "not state the realised maximum number of occasions per patient (209 plasma samples across",
        "76 patients, so most patients contributed fewer than 6). Decomposed inside model() into",
        "binary indicators oc1..oc6 that multiplex the per-occasion IOV etas on log-ka, log-Vc and",
        "log-CL. Monolix estimates a single IOV standard deviation per parameter shared across",
        "occasions (Table 2 gamma), so occasions 2-6 are fixed to the occasion-1 variance.",
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
    n_subjects     = 63,
    n_studies      = 1,
    age_median     = "51 years (IQR 14; whole 76-patient cohort)",
    sex_female_pct = 34.2,
    disease_state  = "adult renal transplant recipients, at least 3 months post-transplant with stable graft function, on mycophenolate + tacrolimus + low-dose prednisone",
    dose_range     = "360-1440 mg/day enteric-coated mycophenolate sodium, given twice daily",
    regions        = "Serbia (University Clinical Centre of Nis)",
    notes          = paste(
      "Tsyplakova 2025 Table 1: 76 patients total contributed 209 MPA plasma samples and 65 saliva",
      "samples; 63 (82.9%) received EC-MPS and 13 (17.1%) received MMF, so this EC-MPS model is",
      "based on the 63-patient subgroup. Median post-transplantation time 70 months (IQR 84.3);",
      "54 (71%) living-donor and 18 (24%) deceased-donor grafts. Laboratory medians: urea 7.8",
      "mmol/L (IQR 5.4), creatinine 136 umol/L (IQR 60), WBC 7.9 x10^9/L, RBC 4.7 x10^12/L,",
      "haemoglobin 138 g/L, haematocrit 41.4%, platelets 225 x10^9/L. Body weight and serum albumin",
      "were NOT recorded, which is why the model carries no allometric or protein-binding term and",
      "why creatinine clearance could not be computed (Section 4). Only steady-state trough (C0)",
      "concentrations were available -- one measurement per occasion -- which is why the authors",
      "selected a one-compartment structure and estimated ka and V by maximum a posteriori",
      "estimation while CL was estimated by maximum likelihood (Sections 2.3 and 3.1.1).",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Tsyplakova 2025 Table 2 and Eqs. 6-8. Only oral data
    # were available, so V and CL are apparent (V/F and CL/F); the authors report
    # them as V and Cl and did not estimate bioavailability separately.
    lka <- log(0.18);    label("Absorption rate constant (1/h)")            # Table 2 kapop = 0.18 (RSE 15.7%); Eq. 6
    lvc <- log(192.42);  label("Apparent central volume of distribution V/F (L)")  # Table 2 Vpop = 192.42 L (RSE 18.8%); Eq. 7
    lcl <- log(9.3);     label("Apparent clearance CL/F at the reference covariate values (L/h)")  # Table 2 Clpop = 9.3 L/h (RSE 8.34%); Eq. 8

    # Covariate effects on apparent clearance. Eq. 1 gives the general continuous
    # form log(Pi) = log(Ppop) + betaC * log(Ci / Cmean) + eta + kappa, so each
    # coefficient is a power exponent on the covariate normalised to Cmean.
    e_pod_cl          <- 0.16; label("Power exponent on (post-transplant time / 67 months) for CL/F (unitless)")  # Table 2 beta(PTP) = 0.16 (RSE 24.7%); Eq. 8
    e_dose_mpa_mgd_cl <- 0.77; label("Power exponent on (total daily MPA dose / 1500 mg/day) for CL/F (unitless)")  # Table 2 beta(TDD) = 0.77 (RSE 19.2%); Eq. 8

    # Inter-individual variability. Table 2 reports the "Standard Deviation of the
    # Random Effects", so each omega below is SD^2 on the log scale.
    etalka ~ 0.36^2;  # Table 2 omega(ka) = 0.36 (RSE 21.2%)
    etalvc ~ 0.52^2;  # Table 2 omega(V)  = 0.52 (RSE 38.8%)
    etalcl ~ 0.27^2;  # Table 2 omega(Cl) = 0.27 (RSE 20.9%)

    # Inter-occasion variability across the six monthly visits. Table 2 reports a
    # single gamma per parameter shared across occasions (Monolix estimates one IOV
    # SD per parameter), so occasions 2-6 are fixed to the occasion-1 variance.
    etaiov_ka_1 ~ 0.28^2;         # Table 2 gamma(ka) = 0.28 (RSE 48.3%)
    etaiov_ka_2 ~ fix(0.28^2);    # shared IOV variance across occasions
    etaiov_ka_3 ~ fix(0.28^2);    # shared IOV variance across occasions
    etaiov_ka_4 ~ fix(0.28^2);    # shared IOV variance across occasions
    etaiov_ka_5 ~ fix(0.28^2);    # shared IOV variance across occasions
    etaiov_ka_6 ~ fix(0.28^2);    # shared IOV variance across occasions

    etaiov_vc_1 ~ 0.52^2;         # Table 2 gamma(V) = 0.52 (RSE 29.2%)
    etaiov_vc_2 ~ fix(0.52^2);    # shared IOV variance across occasions
    etaiov_vc_3 ~ fix(0.52^2);    # shared IOV variance across occasions
    etaiov_vc_4 ~ fix(0.52^2);    # shared IOV variance across occasions
    etaiov_vc_5 ~ fix(0.52^2);    # shared IOV variance across occasions
    etaiov_vc_6 ~ fix(0.52^2);    # shared IOV variance across occasions

    etaiov_cl_1 ~ 0.31^2;         # Table 2 gamma(Cl) = 0.31 (RSE 11.2%)
    etaiov_cl_2 ~ fix(0.31^2);    # shared IOV variance across occasions
    etaiov_cl_3 ~ fix(0.31^2);    # shared IOV variance across occasions
    etaiov_cl_4 ~ fix(0.31^2);    # shared IOV variance across occasions
    etaiov_cl_5 ~ fix(0.31^2);    # shared IOV variance across occasions
    etaiov_cl_6 ~ fix(0.31^2);    # shared IOV variance across occasions

    # Combined additive + proportional residual error (Eq. 5), the structure
    # selected for the EC-MPS model.
    addSd  <- 0.04; label("Additive residual error (mg/L)")        # Table 2 'a' = 0.04 (RSE 26%)
    propSd <- 0.06; label("Proportional residual error (fraction)")  # Table 2 'b' = 0.06 (RSE 24.6%)
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
    # canonical unit); Eq. 8 normalises post-transplant time in MONTHS to 67.
    podMonth <- POD / 30.4375

    # Individual parameters. Eqs. 6-8.
    ka <- exp(lka + etalka + iov_ka)
    vc <- exp(lvc + etalvc + iov_vc)
    cl <- exp(lcl + etalcl + iov_cl) *
      (podMonth / 67)^e_pod_cl *
      (DOSE_MPA_MGD / 1500)^e_dose_mpa_mgd_cl

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
