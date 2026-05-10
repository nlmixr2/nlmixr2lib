Ekhart_2008_carboplatin <- function() {
  description <- "Two-compartment population PK model for free (ultrafilterable) carboplatin in adult cancer patients (Ekhart 2008)"
  reference <- "Ekhart C, Rodenhuis S, Schellens JHM, Beijnen JH, Huitema ADR. Carboplatin dosing in overweight and obese patients with normal renal function, does weight matter? Cancer Chemother Pharmacol. 2009;64(1):115-122. doi:10.1007/s00280-008-0856-x"
  vignette <- "Ekhart_2008_carboplatin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  # Final structural model has no covariate effects: allometric scaling on CL
  # with ABW / IBW / AIBW / Benezet / FFM / LBM each gave delta-OFV < 6.63
  # versus the basic model (Table 6), so weight was not retained.
  covariateData <- list()

  population <- list(
    n_subjects        = 240L,
    n_courses         = 380L,
    n_observations    = 4478L,
    age_range         = "16-75 years",
    age_median        = "47 years",
    weight_range      = "46-170 kg",
    weight_median     = "70 kg",
    height_range      = "153-210 cm",
    height_median     = "171 cm",
    bmi_range         = "16-46 kg/m^2",
    bmi_median        = "24 kg/m^2",
    bsa_range         = "1.49-2.94 m^2",
    bsa_median        = "1.81 m^2",
    sex_female_pct    = 67.1,
    bmi_categories    = "underweight (BMI < 18.5) 7 (3%), normal (18.5-25) 146 (61%), overweight (25-30) 72 (30%), obese (>=30) 15 (6%)",
    serum_creatinine_range_uM           = "18-124 (median 57)",
    creatinine_clearance_range_mL_min   = "55-451 (median 126; Cockcroft-Gault)",
    albumin_range_g_L                   = "18-52 (median 42)",
    disease_state     = "Adult cancer patients with normal renal function receiving carboplatin in combination chemotherapy (NSCLC, ovarian cancer, high-risk and metastatic breast cancer, refractory germ cell cancer, epithelial breast cancer)",
    dose_range        = "Carboplatin 267-600 mg/m^2/day or AUC 6-20 mg.min/mL (Calvert formula); 30 min to 1 h IV infusions; conventional and high-dose CTC / tCTC / miniCTC / paclitaxel-carboplatin regimens",
    regions           = "Netherlands (Antoni van Leeuwenhoek Hospital / Slotervaart Hospital, Amsterdam)",
    notes             = "Pooled across previously published studies (Ekhart 2008 references 10-14). Free (ultrafilterable) platinum measured by flameless atomic absorption spectrometry."
  )

  ini({
    # Structural parameters - basic 2-compartment model fitted with FOCE-INTERACTION
    # after logarithmic data transformation. Values are typical-individual
    # estimates of the population PK parameters (Table 5).
    lcl  <- log(8.38);   label("Carboplatin clearance CL (L/h)")                       # Ekhart 2008 Table 5, CL row
    lvc  <- log(15.4);   label("Central volume of distribution V (L)")                 # Ekhart 2008 Table 5, V row
    lk12 <- log(0.135);  label("Distribution rate constant k12, central->peripheral (1/h)")  # Ekhart 2008 Table 5, k12 row
    lk21 <- log(0.215);  label("Distribution rate constant k21, peripheral->central (1/h)")  # Ekhart 2008 Table 5, k21 row

    # Inter-individual variability (omega^2 = log(1 + CV^2)).
    # CV% values are taken from Table 5 column "% IIV (RSE %)".
    etalcl  ~ 0.0369  # log(1 + 0.194^2); IIV CL = 19.4% CV
    etalvc  ~ 0.0208  # log(1 + 0.145^2); IIV V  = 14.5% CV
    etalk12 ~ 0.2089  # log(1 + 0.482^2); IIV k12 = 48.2% CV
    etalk21 ~ 0.1740  # log(1 + 0.436^2); IIV k21 = 43.6% CV
    # Inter-occasion variability is reported in Table 5 (CL: 9.14% CV, V: 10.8% CV)
    # but is not encoded here - capturing IOV would require an OCC covariate
    # column and per-occasion eta multiplexing. See vignette Errata.

    # Residual error (proportional, on linear scale; the source fitted a
    # log-additive residual after logarithmic data transformation, which is
    # equivalent to a proportional error in linear space for nlmixr2).
    propSd <- 0.197;     label("Proportional residual error (fraction)")              # Ekhart 2008 Table 5, residual error row
  })

  model({
    # Individual PK parameters preserve the paper's micro-constant
    # parameterization (CL, V, k12, k21) so the IIV magnitudes from Table 5
    # apply on their native scale. Q and Vp are derivable as
    #   q  = k12 * vc,  vp = q / k21
    # but were not separately estimated.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, central volume in L -> mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
