Raigani_2017_alteplase_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous-bolus model for the",
    "plasma fibrinolytic activity of Alteplase (wild-type recombinant tissue",
    "plasminogen activator) in female Wistar rats after an i.v. bolus dose of",
    "4000 IU, fit as the reference comparator alongside the novel mt-PA",
    "variant (Raigani 2017; 6 rats). Plasma activity was fit to the",
    "bi-exponential Cp(t) = D*(A*exp(-alpha*t) + B*exp(-beta*t)) with a",
    "nonlinear least-squares program; Alteplase (Dose = 4000) is the",
    "reference level of the treatment covariate (Table 3), converted to",
    "canonical CL/Vc/Q/Vp by the standard bi-exponential relationships. The",
    "companion mt-PA fit from the same paper is Raigani_2017_mtpa_rat.",
    "Volumes are per-animal (mL) and clearances per-animal (mL/min) for a",
    "230-250 g rat; the dose and concentration are expressed in",
    "international units (IU) of fibrinolytic activity, so no molar amount or",
    "body-weight normalization is defined. No inter-individual variance or",
    "residual error model was reported, so the model is typical-value only."
  )
  reference <- paste(
    "Raigani M, Rouini MR, Golabchifar AA, Mirabzadeh E, Vaziri B,",
    "Barkhordari F, Davami F, Mahboudi F (2017). Scale up and",
    "pharmacokinetic study of a novel mutated chimeric tissue plasminogen",
    "activator (mt-PA) in rats. Scientific Reports 7:43028.",
    "doi:10.1038/srep43028.",
    sep = " "
  )
  vignette <- "Raigani_2017_mtpa"
  units <- list(time = "min", dosing = "IU", concentration = "IU/mL")

  # Issue #482: the measured signal is plasma fibrinolytic activity (IU/mL),
  # so the state amount is in IU of activity rather than a mass.
  compartmentData <- list(
    central = list(analyte = "alteplase", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alteplase", units = "IU", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar, female)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "12-14 weeks",
    weight_range = "230-250 g",
    sex_female_pct = 100,
    disease_state = "Healthy female Wistar rats (Pasteur Institute of Iran).",
    dose_range = paste(
      "Single i.v. bolus of Alteplase 4000 IU in 0.5 mL into the tail vein",
      "(positive control). Blood sampled before and at 2, 5, 10, 15, 20, 30",
      "and 45 min after injection."
    ),
    notes = paste(
      "Plasma fibrinolytic activity measured by an amidolytic (Chromolize",
      "t-PA) assay and fit to a two-compartment bi-exponential model (Table",
      "1 micro-constants and volumes; Table 3 macro-constant regression",
      "estimates). Alteplase is the reference treatment level; the novel",
      "mt-PA variant is extracted as Raigani_2017_mtpa_rat. Individual-rat",
      "parameter spread is shown graphically in Figure 5 but was not reported",
      "as variance components; no residual error model was reported."
    )
  )

  ini({
    # Two-compartment i.v.-bolus disposition for the Alteplase reference fit.
    #
    # The paper estimated the bi-exponential macro-constants (Table 3,
    # Dose = 4000, the reference level) of
    #   Cp(t) = D*(A*exp(-alpha*t) + B*exp(-beta*t))                (Eq. 1)
    # with A = 0.00176 1/mL, B = 0.000191 1/mL, alpha = 0.383 1/min,
    # beta = 0.0811 1/min. Because Cp carries the dose factor D, A and B have
    # units 1/mL and
    #   Vc = 1/(A+B)                                                (Eq. 2)
    #   k21 = (A*beta + B*alpha)/(A+B);  k10 = alpha*beta/k21;  k12 = alpha+beta-k21-k10
    #   CL = k10*Vc;  Q = k12*Vc;  Vp = Q/k21
    # giving the canonical values below (Vc = 512.6 mL, Vp = 336.9 mL,
    # Vss = 849.5 mL, CL = 143.9 mL/min, t1/2,beta = 8.5 min), reproducing the
    # Alteplase entries of Table 1 (V1 509.3, V2 328.7, Vss 838.0, CLp 142.6,
    # k10 0.3, k12 0.07, k21 0.12, t1/2,beta 8.3 min).
    lcl <- log(143.88); label("Clearance CL (mL/min)") # Derived from Table 3 Alteplase macro-constants via Eq. 2-4; cf. Table 1 CLp 142.6 mL/min
    lvc <- log(512.6); label("Central volume Vc (mL)") # Derived from Table 3 Alteplase: Vc = 1/(A+B); cf. Table 1 V1 509.3 mL
    lq <- log(37.28); label("Intercompartmental clearance Q (mL/min)") # Derived from Table 3 Alteplase: Q = k12*Vc
    lvp <- log(336.9); label("Peripheral volume Vp (mL)") # Derived from Table 3 Alteplase: Vp = Q/k21; cf. Table 1 V2 328.7 mL

    # No inter-individual variance or residual error model was reported (the
    # fit is a nonlinear least-squares regression; individual spread appears
    # only graphically in Figure 5). Encoded as a typical-value model with a
    # fixed-zero proportional residual, following Betts_2018_mAb_tg32mouse.R.
    propSd <- fixed(0)
    label("Proportional residual error (fraction)")
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
