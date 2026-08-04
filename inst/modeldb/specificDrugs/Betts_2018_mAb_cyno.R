Betts_2018_mAb_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Class-level typical two-compartment population PK model for monoclonal antibodies with linear clearance in cynomolgus monkeys (Betts 2018, n=23 mAbs)."
  reference <- "Betts A, Keunecke A, van Steeg TJ, van der Graaf PH, Avery LB, Jones H, Berkhout J. Linear pharmacokinetic parameters for monoclonal antibodies are similar within a species and across different pharmacological targets: A comparison between human, cynomolgus monkey and hFcRn Tg32 transgenic mouse using a population-modeling approach. MAbs. 2018;10(5):751-64. doi:10.1080/19420862.2018.1462429"
  vignette <- "Betts_2018_mAb_linear_PK"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "mAb", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "cynomolgus monkey",
    n_subjects     = NA,
    n_studies      = NA,
    disease_state  = "Cynomolgus monkeys receiving single-dose IV monoclonal antibody at 1-3 dose levels (n=2 per dose) as part of preclinical PK studies.",
    dose_range     = "IV: 0.01-200 mg/kg across the 23 mAbs (individual mAb linear dose ranges in Betts 2018 Table 1).",
    notes          = "Class-level typical values pooled from Pfizer internal data for 23 monoclonal antibodies with linear PK in cynomolgus monkeys. Reference body weight assumed 3 kg per Betts 2018 Materials and methods. Non-linear (target-mediated) doses were excluded from the source fit."
  )

  ini({
    # Structural parameters at 3 kg reference. Betts 2018 Table 2 Cynomolgus Monkey
    # column reports CL, V1, Q, V2 in weight-normalized units (mL/h/kg or mL/kg).
    # Values converted here to absolute per 3 kg cynomolgus monkey, day and L units:
    #   CL = 0.27 mL/h/kg * 3 kg * 24 h/day / 1000 mL/L = 0.01944 L/day
    #   V1 = 39.29 mL/kg * 3 kg / 1000                   = 0.11787 L
    #   Q  = 1.00 mL/h/kg * 3 kg * 24 / 1000             = 0.072 L/day
    #   V2 = 27.56 mL/kg * 3 kg / 1000                   = 0.08268 L
    lcl <- log(0.01944); label("Clearance (L/day)")                                # Betts 2018 Table 2, Cynomolgus Monkey column: CL 0.27 mL/h/kg
    lvc <- log(0.11787); label("Central volume of distribution (L)")               # Betts 2018 Table 2, Cynomolgus Monkey column: V1 39.29 mL/kg
    lq  <- log(0.072);   label("Inter-compartmental clearance (L/day)")            # Betts 2018 Table 2, Cynomolgus Monkey column: Q 1.00 mL/h/kg
    lvp <- log(0.08268); label("Peripheral volume of distribution (L)")            # Betts 2018 Table 2, Cynomolgus Monkey column: V2 27.56 mL/kg

    # IIV block on CL and V1 (variances and covariance on log scale).
    # Betts 2018 Table 2 Cynomolgus Monkey column: IIV CL = 0.38, COV CL-V1 = 0.09, IIV V1 = 0.10.
    etalcl + etalvc ~ c(0.38,
                        0.09, 0.10)

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
