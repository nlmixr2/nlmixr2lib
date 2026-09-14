Tegenge_2020_immunoglobulin <- function() {
  description <- "Two-compartment population PK model for intravenous polyclonal immunoglobulin G in very low birth-weight neonates (Tegenge 2020)"
  reference   <- "Tegenge MA, Mahmood I. Population pharmacokinetics of immunoglobulin intravenous preparation in very low birth weight neonates. Int Immunopharmacol. 2020;80:106192. doi:10.1016/j.intimp.2019.106192 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 35)"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (birth weight in this cohort)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling with an ESTIMATED exponent of 3.5e-07 -- numerically indistinguishable from zero, so the final model is effectively weight-independent. van der Zeeuw 2026 section 3.2.1.5 reports the exponent and the Discussion (section 4) attributes the near-zero value to the narrow birth-weight range (0.78-1.38 kg), which left the exponent unidentifiable. The reference weight is not reported by the review; 1.08 kg (the midpoint of the birth-weight range in Table 1) is used here, and because the exponent is ~0 the choice of reference is numerically immaterial.",
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 20L,
    n_studies        = 1L,
    age_range        = "1-6 days postnatal",
    gestational_age  = "26-35 weeks",
    weight_range     = "0.78-1.38 kg (birth weight)",
    weight_median    = "Not reported",
    sex_female_pct   = round(100 * 8 / 20, 1),
    race_ethnicity   = "Not reported",
    disease_state    = "Very low birth-weight preterm neonates receiving intravenous immunoglobulin",
    dose_range       = "IVIg 500 or 750 mg/kg",
    regions          = "United States",
    notes            = "Single study (van der Zeeuw 2026 Table 1; the review cites doi:10.1016/s0022-3476(88)80070-2 as the source of the PK data). Baseline IgG not reported in Table 1. This is the only model in the review developed exclusively in neonates, and the review EXCLUDED it from the cross-model Figure 2 simulation for that reason (section 3.2.3). Clearance scaled to 70 kg is roughly twice that of older children and adults (0.0027 vs 0.0014 L/day per kg), prompting the authors to suggest higher dosing in this population."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Tegenge et al.
    # (2020) [35]'. Absolute values for a typical very low birth-weight
    # neonate, in L and L/day. IVIg only, so no depot, Ka or bioavailability.
    lcl     <- log(0.0027); label("Clearance for a typical VLBW neonate (L/day)")       # van der Zeeuw 2026 Table 4: CL = 0.0027
    lvc     <- log(0.008);  label("Central volume of distribution (L)")                 # van der Zeeuw 2026 Table 4: Vc = 0.008
    lq      <- log(0.045);  label("Intercompartmental clearance (L/day)")               # van der Zeeuw 2026 Table 4: Q = 0.045
    lvp     <- log(0.055);  label("Peripheral volume of distribution (L)")              # van der Zeeuw 2026 Table 4: Vp = 0.055

    # Allometric exponent -- estimated, essentially zero. Retained explicitly
    # so the reported structure is preserved even though it contributes a
    # factor of ~1 at every weight in the observed range.
    e_wt    <- 3.5e-07; label("Allometric exponent on CL, Q, Vc and Vp (unitless)")     # van der Zeeuw 2026 section 3.2.1.5: estimated allometric scaling component 3.5e-07

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'Tegenge et al. set
    # the endogenous IgG concentration to 5 g/L based on observed pre-infusion
    # IgG concentrations in neonates [35]'. Held constant, so wrapped in fixed().
    bl_igg  <- fixed(5); label("Endogenous (pre-infusion) IgG concentration (g/L)")     # van der Zeeuw 2026 section 3.2.1.6, held constant at 5 g/L

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2). Values confirmed twice in the source --
    # Table 4 and the prose of section 3.2.1.5.
    etalcl ~ 0.239332  # 52% CV; van der Zeeuw 2026 Table 4 and section 3.2.1.5 IIV 'CL = 52'
    etalvc ~ 0.047265  # 22% CV; van der Zeeuw 2026 Table 4 and section 3.2.1.5 IIV 'Vc = 22'
    etalvp ~ 0.080750  # 29% CV; van der Zeeuw 2026 Table 4 and section 3.2.1.5 IIV 'Vp = 29'

    # Residual error. van der Zeeuw 2026 Table 4 prints '9%' in the
    # proportional column and '-' in the additive column.
    propSd <- 0.09; label("Proportional residual error (fraction)")                     # van der Zeeuw 2026 Table 4: Prop = 9%
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 1.08)^e_wt
    vc <- exp(lvc + etalvc) * (WT / 1.08)^e_wt
    q  <- exp(lq)           * (WT / 1.08)^e_wt
    vp <- exp(lvp + etalvp) * (WT / 1.08)^e_wt

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Intravenous administration only: doses go directly to `central`.
    # States hold EXOGENOUS (therapeutic) IgG only.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
