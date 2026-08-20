Yamakawa_2011_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral imatinib in Japanese adults with ",
    "chronic myeloid leukemia (Yamakawa 2011). No covariate is retained: ",
    "the source paper found single-nucleotide polymorphisms in the influx ",
    "transporter SLCO1B3 and the efflux transporter ABCB1 to be associated ",
    "with imatinib clearance, but did not build a covariate model on the ",
    "population PK model. Inter-individual variability is estimated ",
    "independently on CL/F, Vc/F and ka; residual error is additive only. ",
    "TRANSCRIBED FROM A SECONDARY SOURCE: the parameter values come from ",
    "Table 1 of Yang 2025, an external evaluation of 15 published imatinib ",
    "population PK models, not from the primary publication. Re-extract ",
    "from Yamakawa 2011 when that paper is obtained."
  )
  reference <- paste0(
    "Yamakawa Y, Hamada A, Nakashima R, Yuki M, Hirayama C, Kawaguchi T, ",
    "Saito H. Association of genetic polymorphisms in the influx ",
    "transporter SLCO1B3 and the efflux transporter ABCB1 with imatinib ",
    "pharmacokinetics in patients with chronic myeloid leukemia. Ther Drug ",
    "Monit. 2011;33(2):244-250. doi:10.1097/FTD.0b013e31820beb02. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Yamakawa et al. ",
    "(2011)' and Table 1 footnote d. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SNP_SLCO1B3 = list(
      description = "SLCO1B3 (OATP1B3) genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "Yang 2025 Table 1 footnote d: 'The paper showed that single ",
        "nucleotide polymorphisms in SLCO1B3 and ABCB1 were associated ",
        "with imatinib clearance but did not build a covariate model on ",
        "the population PK model.' The association is therefore reported ",
        "in the primary but carries no coefficient that could be encoded, ",
        "and Yang 2025 Table 1 lists this row's covariates as 'None'. The ",
        "specific variant identifiers are not given by the secondary ",
        "source, so no SNP_<GENE>_RS<rsid> canonical is asserted; read the ",
        "primary to recover them. SLCO1B3 encodes the hepatic uptake ",
        "transporter OATP1B3, for which imatinib is a substrate."
      )
    ),
    SNP_ABCB1 = list(
      description = "ABCB1 (P-glycoprotein / MDR1) genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "See the SNP_SLCO1B3 note above: associated with imatinib ",
        "clearance in the primary but not carried into the population PK ",
        "covariate model, per Yang 2025 Table 1 footnote d. ABCB1 encodes ",
        "the efflux transporter P-glycoprotein, for which imatinib is both ",
        "a substrate and an inhibitor."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    n_observations = "622 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "21-80 years",
    disease_state  = "Japanese adults with chronic myeloid leukemia (CML)",
    dose_range     = "Oral imatinib 100-600 mg total daily dose",
    regions        = "Japan",
    bioanalytical  = "HPLC; limit of quantification not reported in Yang 2025 Table 1",
    notes          = paste0(
      "Demographic detail beyond the row above (weight range, sex split) ",
      "is not reported by the secondary source and must be read from the ",
      "primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Yamakawa row) -----
    lka <- log(2.06); label("First-order absorption rate constant ka (1/h)")  # Yang 2025 Table 1: Ka = 2.06
    lcl <- log(8.7); label("Apparent oral clearance CL/F (L/h)")  # Yang 2025 Table 1: CL/F = 8.7
    lvc <- log(430); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 430

    # ----- Inter-individual variability -----
    # This is the only row of Yang 2025 Table 1 that reports its BSV using
    # the symbol omega rather than CV% or Var(eta): 'omega_CL: 0.363,
    # omega_Vc: 0.457, omega_Ka: 0.747'. Omega is the conventional symbol
    # for the STANDARD DEVIATION of the log-scale random effect, so the
    # variances below are the squares of the tabulated values. The
    # secondary source's deliberate use of three different labels across
    # the table -- Var(eta) for Judson, CV% for most rows, omega here --
    # indicates each row is labelled with the scale it actually reports.
    # No covariances are tabulated, so the three etas are independent.
    etalcl ~ 0.131769  # Yang 2025 Table 1: omega_CL = 0.363 -> omega^2
    etalvc ~ 0.208849  # Yang 2025 Table 1: omega_Vc = 0.457 -> omega^2
    etalka ~ 0.558009  # Yang 2025 Table 1: omega_Ka = 0.747 -> omega^2

    # ----- Residual unexplained variability -----
    # This is the only model among the 15 with a PURELY ADDITIVE residual
    # error and no proportional component. The additive term is tabulated
    # as 0.4 mg/L; Cc is reported in ng/mL, so 0.4 mg/L x 1000 =
    # 400 ng/mL. That is a large absolute error relative to the 1000-3000
    # ng/mL imatinib trough target range, and it means this model's
    # residual is essentially constant across the concentration range
    # rather than scaling with it.
    addSd <- 400; label("Additive residual error (ng/mL)")  # Yang 2025 Table 1: Add 0.4 mg/L = 400 ng/mL
  })

  model({
    # ----- 1. Individual parameters -----
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL,
    # the unit in which the additive residual error is expressed.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
