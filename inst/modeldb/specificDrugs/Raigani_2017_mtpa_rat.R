Raigani_2017_mtpa_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous-bolus model for the",
    "plasma fibrinolytic activity of mt-PA, a novel mutated chimeric tissue",
    "plasminogen activator (deletion of the finger, EGF and kringle-1",
    "domains of t-PA, insertion of a GHRP fibrin-binding tetrapeptide, and",
    "a KHRR 296-299 AAAA substitution conferring PAI-1 resistance), in",
    "female Wistar rats after i.v. bolus doses of 100, 200 or 400 IU",
    "(Raigani 2017; 13 rats). The authors fit plasma activity to the",
    "bi-exponential Cp(t) = D*(A*exp(-alpha*t) + B*exp(-beta*t)) with a",
    "nonlinear least-squares program; treatment (mt-PA vs the Alteplase",
    "comparator) and mt-PA dose were tested as log-additive covariates on A,",
    "B, alpha and beta (Table 3). mt-PA dose was not significant across the",
    "three doses, so the values here are the pooled mt-PA estimate",
    "(Alteplase reference plus the mt-PA treatment coefficient), converted to",
    "canonical CL/Vc/Q/Vp by the standard bi-exponential relationships. The",
    "companion Alteplase fit from the same paper is",
    "Raigani_2017_alteplase_rat. Volumes are per-animal (mL) and clearances",
    "per-animal (mL/min) for a 230-250 g rat; the dose and the concentration",
    "are expressed in international units (IU) of fibrinolytic activity, so",
    "no molar amount or body-weight normalization is defined. No",
    "inter-individual variance or residual error model was reported, so the",
    "model is typical-value only."
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

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The measured signal is plasma fibrinolytic activity
  # (IU/mL), so the state amount is in IU of activity rather than a mass.
  compartmentData <- list(
    central = list(analyte = "mtPA", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mtPA", units = "IU", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar, female)",
    n_subjects = 13L,
    n_studies = 1L,
    age_range = "12-14 weeks",
    weight_range = "230-250 g",
    sex_female_pct = 100,
    disease_state = "Healthy female Wistar rats (Pasteur Institute of Iran).",
    dose_range = paste(
      "Single i.v. bolus of mt-PA 100 IU (n = 3), 200 IU (n = 4) or 400 IU",
      "(n = 6) in 0.5 mL into the tail vein. Blood sampled before and at 2,",
      "5, 10, 15, 20, 30, 45, 60, 75, 90 and 120 min after injection."
    ),
    notes = paste(
      "Plasma fibrinolytic activity measured by an amidolytic (Chromolize",
      "t-PA) assay and fit to a two-compartment bi-exponential model (Table",
      "1 micro-constants and volumes; Table 3 macro-constant regression",
      "estimates). Covariate testing (Wald) found no significant difference",
      "among the three mt-PA doses, so the three dose groups are pooled. The",
      "paper's companion drug Alteplase (4000 IU, n = 6) is extracted as",
      "Raigani_2017_alteplase_rat. Individual-rat parameter spread is shown",
      "graphically in Figure 5 but was not reported as variance components;",
      "no residual error model was reported."
    )
  )

  ini({
    # Two-compartment i.v.-bolus disposition for the pooled mt-PA fit.
    #
    # The paper estimated the bi-exponential macro-constants (Table 3) of
    #   Cp(t) = D*(A*exp(-alpha*t) + B*exp(-beta*t))                (Eq. 1)
    # for an Alteplase reference group (Dose = 4000) and applied log-additive
    # covariate effects for the mt-PA treatment (Eq. 7:
    # log P_test = log P_reference + beta_cov). mt-PA dose was not significant
    # (Table 3, p > 0.05 for the dose coefficients on alpha; the three mt-PA
    # dose groups do not differ), so the pooled mt-PA macro-constants are the
    # reference values shifted by the mt-PA coefficient:
    #   alpha = exp(log(0.383)  + (-0.0362)) = 0.3694 1/min   # Table 3 alpha(4000)=0.383; beta(alpha,Dose=mt-PA)=-0.0362
    #   beta  = exp(log(0.0811) + (-0.981))  = 0.03041 1/min  # Table 3 beta(4000)=0.0811; beta(beta,Dose=mt-PA)=-0.981
    #   A     = exp(log(0.00176)+  1.96)     = 0.01249 1/mL   # Table 3 A(4000)=0.00176; beta(A,Dose=mt-PA)=1.96
    #   B     = exp(log(0.000191)+ 3.31)     = 0.005231 1/mL  # Table 3 B(4000)=0.000191; beta(B,Dose=mt-PA)=3.31
    # Because Cp carries the dose factor D, A and B have units 1/mL and
    #   Vc = 1/(A+B)                                           (Eq. 2, D/(A+B) with Cp = D*(...))
    #   k21 = (A*beta + B*alpha)/(A+B);  k10 = alpha*beta/k21;  k12 = alpha+beta-k21-k10
    #   CL = k10*Vc;  Q = k12*Vc;  Vp = Q/k21
    # giving the canonical values below (Vc = 56.4 mL, Vp = 79.3 mL, Vss =
    # 135.7 mL, CL = 4.86 mL/min, t1/2,beta = 22.8 min), which reproduce the
    # pooled mt-PA entries of Table 1 (Vc 48-59, Vss 91-164, CLp 3.9-5.9,
    # t1/2,beta 19.1-26.1 across the three doses).
    lcl <- log(4.86); label("Clearance CL (mL/min)") # Derived from Table 3 pooled mt-PA macro-constants via Eq. 2-4; cf. Table 1 CLp 3.9-5.9 mL/min
    lvc <- log(56.4); label("Central volume Vc (mL)") # Derived from Table 3 pooled mt-PA: Vc = 1/(A+B); cf. Table 1 V1 48.3-58.7 mL
    lq <- log(10.34); label("Intercompartmental clearance Q (mL/min)") # Derived from Table 3 pooled mt-PA: Q = k12*Vc
    lvp <- log(79.3); label("Peripheral volume Vp (mL)") # Derived from Table 3 pooled mt-PA: Vp = Q/k21; cf. Table 1 V2 42.9-105.2 mL

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
