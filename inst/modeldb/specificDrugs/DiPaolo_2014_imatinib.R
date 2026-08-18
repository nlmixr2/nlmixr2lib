DiPaolo_2014_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral imatinib in Italian adults with ",
    "chronic myeloid leukemia (Di Paolo 2014). Both CL/F and Vc/F are ",
    "reduced by a THRESHOLD effect of alpha-1 acid glycoprotein above 0.94 ",
    "g/L (94 mg/dL), and CL/F is additionally reduced by 11.8% in carriers ",
    "of the SLC22A1 (hOCT1) c.480C>G rs683369 variant allele, the ",
    "pharmacogenetic association after which the source paper is titled. ",
    "Inter-individual variability is estimated on CL/F only; residual ",
    "error is proportional. TRANSCRIBED FROM A SECONDARY SOURCE: the ",
    "parameter values come from Table 1 of Yang 2025, an external ",
    "evaluation of 15 published imatinib population PK models, not from ",
    "the primary publication. Re-extract from Di Paolo 2014 when that ",
    "paper is obtained."
  )
  reference <- paste0(
    "Di Paolo A, Polillo M, Capecchi M, Cervetti G, Barate C, Angelini S, ",
    "Guerrini F, Fontanelli G, Arici R, Ciabatti E, Grassi S, Bocci G, ",
    "Hrelia P, Danesi R, Petrini M, Galimberti S. The c.480C>G ",
    "polymorphism of hOCT1 influences imatinib clearance in patients ",
    "affected by chronic myeloid leukemia. Pharmacogenomics J. ",
    "2014;14(4):328-335. doi:10.1038/tpj.2014.7. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Di Paolo et al. ",
    "(2014)' and Table 1 footnote f. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    AAG = list(
      description        = "Alpha-1 acid glycoprotein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "UNIT CONVERSION IS REQUIRED RELATIVE TO THE PRINTED THRESHOLD. ",
        "Yang 2025 Table 1 prints the effect as ",
        "CL/F = 13.093 x (0.775 if AGP > 94 mg/dL) and ",
        "Vc/F = 359.238 x (0.793 if AGP > 94 mg/dL) -- note the threshold ",
        "is given in mg/dL for this row, whereas the same table's ",
        "abbreviation list gives the general AGP unit as g/L (used by ",
        "Petain_2008_imatinib.R). 94 mg/dL = 0.94 g/L, so the model tests ",
        "AAG > 0.94 against the canonical g/L column. That threshold sits ",
        "just above the upper end of the healthy adult reference range ",
        "(roughly 0.5-1.0 g/L), which is a sensible cut point for ",
        "identifying the acute-phase-elevated subgroup.",
        " ",
        "This is a THRESHOLD (step) effect, not the power-form continuous ",
        "effect used by Petain_2008_imatinib.R: above the cut point CL/F ",
        "is multiplied by 0.775 and Vc/F by 0.793, and below it both ",
        "factors are exactly 1. The direction matches the mechanism -- ",
        "imatinib is extensively bound to AAG, so a higher AAG lowers the ",
        "unbound fraction and therefore both the apparent clearance and ",
        "the apparent volume measured on total plasma concentrations. ",
        "Because the two factors are nearly equal, the effect is close to ",
        "a pure shift in both parameters, leaving the elimination ",
        "half-life almost unchanged while raising total-drug exposure by ",
        "about 29%."
      ),
      source_name        = "AGP"
    ),
    SNP_SLC22A1_RS683369 = list(
      description        = "SLC22A1 (OCT1) rs683369 c.480C>G L160F variant carrier indicator; 1 = CG or GG, 0 = CC",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CC wild-type homozygote, whose CL/F factor is 1)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). Yang 2025 Table 1 ",
        "footnote f gives the coefficient explicitly: 'When SLC22A1 CG or ",
        "GG theta_SLC22A1 = 0.882, else theta_SLC22A1 = 1 when SLC22A1 ",
        "CC'. This is a DOMINANT genetic model -- heterozygotes and ",
        "variant homozygotes share one coefficient -- so a single carrier ",
        "indicator is sufficient and paired _HET / _HOM indicators (as ",
        "used for ABCG2 in Jiang_2023_imatinib.R) would be unidentifiable ",
        "here. Carriers have an 11.8% lower apparent oral clearance. ",
        "SLC22A1 encodes OCT1, the principal hepatic uptake transporter ",
        "for imatinib, so a reduced-function variant lowering hepatocellular ",
        "uptake and therefore apparent clearance is mechanistically ",
        "coherent. Yang 2025 Table 1 lists this row's covariate under the ",
        "protein alias 'hOCT1'; the primary is titled after this ",
        "association."
      ),
      source_name        = "hOCT1 c.480C>G (theta_SLC22A1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    n_observations = "117 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "27-79 years",
    disease_state  = "Adults with chronic myeloid leukemia (CML)",
    dose_range     = "Oral imatinib 200, 400 or 600 mg total daily dose",
    regions        = "Italy",
    bioanalytical  = "HPLC, limit of quantification 25 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "The sparsest sampling of the 15 models evaluated by Yang 2025 -- ",
      "117 samples from 60 patients, under two per subject -- which is ",
      "consistent with the absence of any random effect on Vc/F. ",
      "Demographic detail beyond the row above (weight range, sex split) ",
      "is not reported by the secondary source and must be read from the ",
      "primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Di Paolo row) -----
    # Typical values are for the REFERENCE subject: AAG at or below 0.94
    # g/L and SLC22A1 rs683369 CC wild-type, at which both covariate
    # factors equal 1.
    lka <- log(1.29); label("First-order absorption rate constant ka (1/h)")  # Yang 2025 Table 1: Ka = 1.29
    lcl <- log(13.093); label("Apparent oral clearance CL/F at the reference subject (L/h)")  # Yang 2025 Table 1: CL/F = 13.093 x (0.775 if AGP > 94 mg/dL) x theta_SLC22A1
    lvc <- log(359.238); label("Apparent central volume Vc/F at the reference subject (L)")  # Yang 2025 Table 1: Vc/F = 359.238 x (0.793 if AGP > 94 mg/dL)

    # ----- Threshold effect of alpha-1 acid glycoprotein -----
    # Multiplicative factors applied only when AAG exceeds 0.94 g/L
    # (94 mg/dL as printed in the source).
    e_aag_cl <- 0.775; label("CL/F multiplicative factor when AAG exceeds 0.94 g/L (unitless)")  # Yang 2025 Table 1: (0.775 if AGP > 94 mg/dL)
    e_aag_vc <- 0.793; label("Vc/F multiplicative factor when AAG exceeds 0.94 g/L (unitless)")  # Yang 2025 Table 1: (0.793 if AGP > 94 mg/dL)

    # ----- SLC22A1 rs683369 genotype effect on CL/F (dominant model) -----
    e_slc22a1_cl <- 0.882; label("CL/F multiplicative factor for SLC22A1 rs683369 CG or GG carriers (unitless)")  # Yang 2025 Table 1 footnote f: theta_SLC22A1 = 0.882 for CG or GG, 1 for CC

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports only 'CV%(CL): 19.62%' for this row: no
    # random effect on Vc/F or ka is tabulated, consistent with the very
    # sparse sampling (117 samples from 60 patients). The tabulated CV% is
    # taken as omega (the log-scale standard deviation), so the variance is
    # (CV/100)^2 -- the convention used throughout this transcription and
    # in the shipped Jiang_2023_imatinib.R.
    etalcl ~ 0.038494  # Yang 2025 Table 1: CV%(CL) 19.62% -> omega^2 = 0.1962^2

    # ----- Residual unexplained variability -----
    propSd <- 0.2898; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 28.98%
  })

  model({
    # ----- 1. Covariate factors -----
    # Threshold indicator: 1 when AAG exceeds 0.94 g/L (the published 94
    # mg/dL cut point converted to the canonical g/L column), 0 otherwise.
    # Raising the factor to the power of the indicator collapses it to
    # exactly 1 below the threshold.
    aag_high <- (AAG > 0.94)
    aag_factor_cl <- e_aag_cl^aag_high
    aag_factor_vc <- e_aag_vc^aag_high

    # Dominant genetic model: the factor applies to CG and GG carriers
    # alike and collapses to 1 for CC wild-type homozygotes.
    slc22a1_factor <- e_slc22a1_cl^SNP_SLC22A1_RS683369

    # ----- 2. Individual parameters -----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * aag_factor_cl * slc22a1_factor
    vc <- exp(lvc) * aag_factor_vc

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 5. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
