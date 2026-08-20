Betts_2018_mAb_tg32mouse <- function() {
  description <- "Preclinical (hFcRn Tg32 transgenic mouse). Class-level typical two-compartment population PK model for monoclonal antibodies with linear clearance in mice expressing human neonatal Fc receptor (Betts 2018, n=11 mAbs)."
  reference <- "Betts A, Keunecke A, van Steeg TJ, van der Graaf PH, Avery LB, Jones H, Berkhout J. Linear pharmacokinetic parameters for monoclonal antibodies are similar within a species and across different pharmacological targets: A comparison between human, cynomolgus monkey and hFcRn Tg32 transgenic mouse using a population-modeling approach. MAbs. 2018;10(5):751-64. doi:10.1080/19420862.2018.1462429"
  vignette <- "Betts_2018_mAb_linear_PK"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "mAb tg32mouse", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb tg32mouse", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "hFcRn Tg32 transgenic mouse",
    n_subjects     = NA,
    n_studies      = NA,
    disease_state  = "Homozygous hFcRn Tg32 transgenic mice (human FcRn knock-in) receiving single-dose IV monoclonal antibody at 3.5 or 5 mg/kg (n=5-6 per mAb) as an early PK-screening tool for human PK prediction.",
    dose_range     = "IV: 3.5 mg/kg (1 mAb) or 5 mg/kg (10 mAbs) per Betts 2018 Materials and methods.",
    notes          = "Class-level typical values pooled from Pfizer internal data for 11 monoclonal antibodies with linear PK in hFcRn Tg32 mice. Reference body weight assumed 0.02 kg (20 g) per Betts 2018 Materials and methods."
  )

  ini({
    # Structural parameters at 0.02 kg reference. Betts 2018 Table 2 Tg32 hFcRn
    # Transgenic Mouse column reports CL, V1, Q, V2 in weight-normalized units
    # (mL/h/kg or mL/kg). Values converted here to absolute per 0.02 kg mouse,
    # day and L units:
    #   CL = 0.35 mL/h/kg * 0.02 kg * 24 h/day / 1000 = 0.000168 L/day
    #   V1 = 59.28 mL/kg * 0.02 kg / 1000              = 0.0011856 L
    #   Q  = 4.40 mL/h/kg * 0.02 kg * 24 / 1000        = 0.002112 L/day
    #   V2 = 60.54 mL/kg * 0.02 kg / 1000              = 0.0012108 L
    lcl <- log(0.000168);  label("Clearance (L/day)")                              # Betts 2018 Table 2, Tg32 Mouse column: CL 0.35 mL/h/kg
    lvc <- log(0.0011856); label("Central volume of distribution (L)")             # Betts 2018 Table 2, Tg32 Mouse column: V1 59.28 mL/kg
    lq  <- log(0.002112);  label("Inter-compartmental clearance (L/day)")          # Betts 2018 Table 2, Tg32 Mouse column: Q 4.40 mL/h/kg
    lvp <- log(0.0012108); label("Peripheral volume of distribution (L)")          # Betts 2018 Table 2, Tg32 Mouse column: V2 60.54 mL/kg

    # IIV block on CL and V1 (variances and covariance on log scale).
    # Betts 2018 Table 2 Tg32 Mouse column: IIV CL = 0.41, COV CL-V1 = 0.11, IIV V1 = 0.12.
    etalcl + etalvc ~ c(0.41,
                        0.11, 0.12)

    # Per-compound residual error is in Supplementary Table 1 (not on disk); fixed at 0
    # per the skill's standing policy for unreported RUV (documented in vignette Errata).
    propSd <- fixed(0); label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
