Betts_2018_mAb_combined <- function() {
  description <- "Preclinical + human. Class-level typical two-compartment population PK model with body-weight allometric scaling for monoclonal antibodies with linear clearance, jointly fitted across human, cynomolgus monkey and hFcRn Tg32 mouse (Betts 2018 combined all-species dataset, n=27 mAbs)."
  reference <- "Betts A, Keunecke A, van Steeg TJ, van der Graaf PH, Avery LB, Jones H, Berkhout J. Linear pharmacokinetic parameters for monoclonal antibodies are similar within a species and across different pharmacological targets: A comparison between human, cynomolgus monkey and hFcRn Tg32 transgenic mouse using a population-modeling approach. MAbs. 2018;10(5):751-64. doi:10.1080/19420862.2018.1462429"
  vignette <- "Betts_2018_mAb_linear_PK"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "mAb combined", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb combined", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor for CL, V1, Q and V2 with reference weight 70 kg. Species-representative body weights used in the Betts 2018 fit were 70 kg (human, subject-specific), 3 kg (cynomolgus monkey, assumed), and 0.02 kg (hFcRn Tg32 transgenic mouse, assumed).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human + cynomolgus monkey + hFcRn Tg32 transgenic mouse",
    n_subjects     = NA,
    n_studies      = NA,
    disease_state  = "Combined-species dataset pooling 18 mAbs with clinical data, 23 mAbs with cynomolgus monkey data, and 11 mAbs with hFcRn Tg32 mouse data (total 27 mAbs).",
    dose_range     = "IV: mAb-specific linear dose ranges (see Betts 2018 Table 1). Non-linear doses were excluded.",
    notes          = "Combined-species jointly-fitted 2-compartment model with allometric scaling on all four disposition parameters (CL, V1, Q, V2). Reference body weight 70 kg (human). Species-representative body weights in the source fit: 70 kg human, 3 kg cynomolgus monkey, 0.02 kg hFcRn Tg32 mouse."
  )

  ini({
    # Structural parameters at 70 kg reference. Betts 2018 Table 4 combined
    # 'Tg32 Mouse, Cyno and Human (n=27 mAbs)' column reports CL, V1, Q, V2 in
    # weight-normalized units (mL/h/kg or mL/kg) at the 70 kg human reference:
    #   CL = 0.16 mL/h/kg * 70 kg * 24 h/day / 1000 mL/L = 0.2688 L/day
    #   V1 = 45.19 mL/kg * 70 kg / 1000                   = 3.1633 L
    #   Q  = 0.28 mL/h/kg * 70 kg * 24 / 1000             = 0.4704 L/day
    #   V2 = 30.81 mL/kg * 70 kg / 1000                   = 2.1567 L
    lcl <- log(0.2688); label("Clearance at 70 kg reference (L/day)")              # Betts 2018 Table 4, all-species column: CL 0.16 mL/h/kg
    lvc <- log(3.1633); label("Central volume of distribution at 70 kg (L)")       # Betts 2018 Table 4, all-species column: V1 45.19 mL/kg
    lq  <- log(0.4704); label("Inter-compartmental clearance at 70 kg (L/day)")    # Betts 2018 Table 4, all-species column: Q 0.28 mL/h/kg
    lvp <- log(2.1567); label("Peripheral volume of distribution at 70 kg (L)")    # Betts 2018 Table 4, all-species column: V2 30.81 mL/kg

    # Allometric exponents (unitless). Betts 2018 Table 4 note gives the scaling
    # equation Y_species = Y_typical * (BW / BW_ref)^exponent. All exponents are
    # estimated (not held constant at 0.75 / 1); reported with narrow CIs.
    e_wt_cl <- 0.89; label("Allometric exponent on CL (unitless)")                  # Betts 2018 Table 4, all-species column: alpha
    allo_v1 <- 0.98; label("Allometric exponent on V1 (unitless)")                  # Betts 2018 Table 4, all-species column: beta
    e_wt_q  <- 0.67; label("Allometric exponent on Q (unitless)")                   # Betts 2018 Table 4, all-species column: gamma
    allo_v2 <- 0.95; label("Allometric exponent on V2 (unitless)")                  # Betts 2018 Table 4, all-species column: delta

    # IIV block on CL and V1 (variances and covariance on log scale).
    # Betts 2018 Table 4 all-species column: IIV CL = 0.47, COV CL-V1 = 0.08, IIV V1 = 0.11.
    etalcl + etalvc ~ c(0.47,
                        0.08, 0.11)

    # Per-compound residual error is in Supplementary Table 2 (not on disk); fixed at 0
    # per the skill's standing policy for unreported RUV (documented in vignette Errata).
    propSd <- fixed(0); label("Proportional residual error (fraction)")
  })
  model({
    # WT is body weight in kg. Allometric scaling with reference 70 kg per Betts 2018 Table 4.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_v1
    q  <- exp(lq)           * (WT / 70)^e_wt_q
    vp <- exp(lvp)          * (WT / 70)^allo_v2

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
