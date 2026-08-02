Betts_2018_mAb_human <- function() {
  description <- "Class-level typical two-compartment population PK model for monoclonal antibodies with linear clearance in adult humans (Betts 2018, n=18 mAbs)."
  reference <- "Betts A, Keunecke A, van Steeg TJ, van der Graaf PH, Avery LB, Jones H, Berkhout J. Linear pharmacokinetic parameters for monoclonal antibodies are similar within a species and across different pharmacological targets: A comparison between human, cynomolgus monkey and hFcRn Tg32 transgenic mouse using a population-modeling approach. MAbs. 2018;10(5):751-64. doi:10.1080/19420862.2018.1462429"
  vignette <- "Betts_2018_mAb_linear_PK"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA,
    n_studies      = NA,
    disease_state  = "Adult healthy volunteers or patients receiving therapeutic monoclonal antibodies with linear (non-target-mediated) clearance.",
    dose_range     = "IV: 0.1-800 mg/kg or 22-800 mg total across the 18 mAbs (individual mAb linear dose ranges in Betts 2018 Table 1).",
    notes          = "Class-level typical values pooled from Pfizer internal data for 18 monoclonal antibodies with linear PK in humans (single-dose IV, n=3-24 individuals per dose level per mAb). Reference body weight assumed 70 kg. Non-linear (target-mediated) doses were excluded from the source fit; see Betts 2018 Materials and methods."
  )

  ini({
    # Structural parameters at 70 kg reference. Betts 2018 Table 2 Human column
    # reports CL, V1, Q, V2 in weight-normalized units (mL/h/kg or mL/kg).
    # Values converted here to absolute per 70 kg human, day and L units:
    #   CL = 0.15 mL/h/kg * 70 kg * 24 h/day / 1000 mL/L = 0.252 L/day
    #   V1 = 46.31 mL/kg * 70 kg / 1000 mL/L            = 3.2417 L
    #   Q  = 0.27 mL/h/kg * 70 kg * 24 h/day / 1000     = 0.4536 L/day
    #   V2 = 31.47 mL/kg * 70 kg / 1000                 = 2.2029 L
    lcl <- log(0.252); label("Clearance (L/day)")                                  # Betts 2018 Table 2, Human column: CL 0.15 mL/h/kg
    lvc <- log(3.2417); label("Central volume of distribution (L)")                # Betts 2018 Table 2, Human column: V1 46.31 mL/kg
    lq  <- log(0.4536); label("Inter-compartmental clearance (L/day)")             # Betts 2018 Table 2, Human column: Q 0.27 mL/h/kg
    lvp <- log(2.2029); label("Peripheral volume of distribution (L)")             # Betts 2018 Table 2, Human column: V2 31.47 mL/kg

    # IIV block on CL and V1 (variances and covariance on log scale).
    # Betts 2018 Table 2 Human column: IIV CL = 0.48, COV CL-V1 = 0.09, IIV V1 = 0.09.
    etalcl + etalvc ~ c(0.48,
                        0.09, 0.09)

    # Proportional residual error was estimated per compound (Betts 2018 Table 2 footnote:
    # "Residual errors per compound are shown in Supplementary Table 1"); the supplement
    # is not on disk, so no class-level typical value is available. Fixed at 0 per the
    # skill's standing policy for unreported RUV (documented in vignette Errata).
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
