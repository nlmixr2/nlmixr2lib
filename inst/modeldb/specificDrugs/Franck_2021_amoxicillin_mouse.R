Franck_2021_amoxicillin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Two-compartment population PK submodel for oral",
    "amoxicillin in mice infected intranasally with Streptococcus pneumoniae",
    "serotype 1 (Franck 2021), with first-order absorption after a lag time,",
    "first-order elimination, and a pharmacokinetic interaction in which",
    "intraperitoneal monophosphoryl lipid A (MPLA, 2 mg/kg) coadministration",
    "lowers amoxicillin clearance linearly with the amoxicillin dose",
    "(CL = 124 - 0.145 * dose in ug; -40.9% at 14 mg/kg, -1.17% at 0.4 mg/kg).",
    "Interindividual variability on CL and Q and a proportional residual error.",
    "This is the stand-alone PK submodel of Supplementary Table S1; its typical",
    "values were then fixed in the sequential PK/PD-survival model",
    "Franck_2021_amoxicillin_mpla_mouse.",
    sep = " "
  )
  reference <- paste(
    "Franck S, Michelet R, Casilag F, Sirard JC, Wicha SG, Kloft C.",
    "A Model-Based Pharmacokinetic/Pharmacodynamic Analysis of the Combination",
    "of Amoxicillin and Monophosphoryl Lipid A Against S. pneumoniae in Mice.",
    "Pharmaceutics. 2021;13(4):469. doi:10.3390/pharmaceutics13040469.",
    "PK submodel parameters from Supplementary Table S1 (pharmaceutics-13-00469-s001.pdf).",
    sep = " "
  )
  vignette <- "Franck_2021_amoxicillin_mpla_pneumonia"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  covariateData <- list(
    DOSE_AMOXICILLIN_UG = list(
      description = "Administered single oral amoxicillin dose in ug (absolute amount per mouse)",
      units = "ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through the MPLA pharmacokinetic interaction on clearance",
        "(Supplementary Eq. 1, P = theta1 + theta2 * DOSE), so it has no effect",
        "when CONMED_MPLA = 0. The paper's arithmetic (Results: 73.3 mL/h at",
        "14 mg/kg and 123 mL/h at 0.4 mg/kg with MPLA) reproduces exactly with",
        "the mg/kg dose multiplied by a 25 g body weight (14 mg/kg = 350 ug,",
        "0.4 mg/kg = 10 ug); Supplementary Section S1 gives body weight ~25 g.",
        sep = " "
      ),
      source_name = "DOSE"
    ),
    CONMED_MPLA = list(
      description = "Coadministration of monophosphoryl lipid A (1 = MPLA 2.0 mg/kg IP given with the amoxicillin dose, 0 = amoxicillin alone)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Binary covariate (Supplementary Section S2, 'Pharmacokinetic Submodel').",
        "Only a single MPLA dose level (2.0 mg/kg) was studied and MPLA PK was",
        "not measured, so the interaction is a yes/no effect scaled by the",
        "amoxicillin dose.",
        sep = " "
      ),
      source_name = "MPLA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "amoxicillin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amoxicillin", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "amoxicillin", units = "ug", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "mouse (RjOrl:Swiss / CD-1, female, S. pneumoniae serotype 1 pneumonia model)",
    n_subjects = 106,
    n_studies = 1,
    age_range = "6-8 weeks",
    weight_range = "~25 g",
    sex_female_pct = 100,
    disease_state = paste(
      "Pneumonia after intranasal infection with 1-4 x 10^6 CFU Streptococcus",
      "pneumoniae serotype 1 (clinical isolate E1586, amoxicillin MIC 0.016 mg/L);",
      "treatment given 12 h after infection.",
      sep = " "
    ),
    dose_range = paste(
      "Amoxicillin 0.4 or 14 mg/kg single oral gavage, with or without",
      "monophosphoryl lipid A 2.0 mg/kg intraperitoneally.",
      sep = " "
    ),
    regions = "France (Institut Pasteur de Lille)",
    notes = paste(
      "Supplementary Section S1: 106 total serum amoxicillin concentrations",
      "(15.1% below the LLOQ of 0.01 ug/mL) from RjOrl:Swiss mice, 1-2",
      "retro-orbital samples per mouse at 0.167, 0.5, 1, 2, 3, 6 and 12 h after",
      "dosing, 3-4 samples per time point and group. Total (not unbound)",
      "concentrations were modelled; protein binding is ~17%. Mouse type",
      "(RjOrl:Swiss vs Balb/cJRj) was tested and not retained.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Supplementary Table S1 ('Structural Submodel').
    # All are apparent (/F) values; F was fixed to 1.
    lka <- fixed(log(5.04))
    label("First-order absorption rate constant ka (1/h)") # Table S1: ka = 5.04 1/h, fixed (footnote *: fixed during model development for stability)
    ltlag <- log(0.125)
    label("Absorption lag time tlag (h)") # Table S1: tlag = 0.125 h (RSE 10.0%)
    lvc <- log(15.4)
    label("Apparent central volume of distribution Vc/F (mL)") # Table S1: Vc/F = 15.4 mL (RSE 25.5%)
    lvp <- log(50.7)
    label("Apparent peripheral volume of distribution Vp/F (mL)") # Table S1: Vp/F = 50.7 mL (RSE 10.1%)
    lq <- log(71.9)
    label("Apparent intercompartmental clearance Q/F (mL/h)") # Table S1: Q/F = 71.9 mL/h (RSE 17.4%)
    lcl <- log(124)
    label("Apparent amoxicillin clearance without MPLA, CL/F (mL/h)") # Table S1: CL_AMX/F = 124 mL/h (RSE 5.90%)
    lfdepot <- fixed(log(1))
    label("Oral bioavailability F (fraction)") # Table S1 abbreviations: 'F: Bioavailability of AMX fixed to 1'

    # MPLA pharmacokinetic interaction (Supplementary Eq. 1, P = theta1 + theta2 * DOSE)
    e_dose_mpla_cl <- -0.145
    label("Additive change in CL/F per ug amoxicillin dose when MPLA is coadministered (mL/h/ug)") # Table S1: FC_AMX+MPLA = -0.145 mL/h/ug (RSE 18.1%)

    # Interindividual variability -- Table S1 reports %CV; omega^2 = log(CV^2 + 1)
    etalcl ~ 0.05331 # Table S1: IIV CL 23.4 %CV -> log(0.234^2 + 1) = 0.05331
    etalq ~ 0.06392 # Table S1: IIV Q 25.7 %CV -> log(0.257^2 + 1) = 0.06392

    # Residual error
    propSd <- 0.282
    label("Proportional residual error on serum amoxicillin concentration (fraction)") # Table S1: proportional RUV 28.2 %CV; Supplementary Section S3
  })
  model({
    # Individual parameters. Clearance carries the additive MPLA x dose
    # interaction of Supplementary Eq. 1 inside the exponential IIV.
    ka <- exp(lka)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq + etalq)
    cl <- (exp(lcl) + e_dose_mpla_cl * DOSE_AMOXICILLIN_UG * CONMED_MPLA) * exp(etalcl)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)
    alag(depot) <- exp(ltlag)

    # Total amoxicillin serum concentration (ug in mL = ug/mL = mg/L)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
