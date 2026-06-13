Garot_2011_ceftriaxone <- function() {
  description <- "Two-compartment IV-infusion population PK model for ceftriaxone in critically ill adult ICU patients with sepsis, severe sepsis, or septic shock (Garot 2011)"
  reference <- "Garot D, Respaud R, Lanotte P, Simon N, Mercier E, Ehrmann S, Perrotin D, Dequin PF, Le Guellec C. Population pharmacokinetics of ceftriaxone in critically ill septic patients: a reappraisal. Br J Clin Pharmacol. 2011;72(5):758-767. doi:10.1111/j.1365-2125.2011.04005.x"
  vignette <- "Garot_2011_ceftriaxone"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Measured creatinine clearance (24-hour urine collection, raw mL/min; NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLcr. Measured by 24-hour urine collection (not BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form; precedent: Delattre_2010_amikacin.R). The paper's Table 3 final-model equation centres the covariate on the population median CLcr expressed in L/h: CL = theta1 + theta2 * (CLcr / 4.26), with CLcr in L/h. 4.26 L/h is equivalent to 71 mL/min (4.26 * 1000 / 60). The model below carries CRCL in mL/min and divides by 71 mL/min to match the paper's normalization. The text reports a population median CLcr of 68.5 mL/min (range 5.5-214 mL/min), slightly different from the 71 mL/min normalization constant; the discrepancy is preserved as published rather than re-centred. Haemofiltration (n=12) was tested and had no effect on PK parameters, so HF status is not carried as a covariate.",
      source_name        = "CLcr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 1L,
    age_range      = "35-86 years",
    age_median     = "68 years (mean)",
    weight_range   = "Not reported in the paper",
    weight_median  = "Not reported in the paper",
    sex_female_pct = 28,
    race_ethnicity = "Not reported (single French university-hospital ICU population)",
    disease_state  = "Sepsis (n=19), severe sepsis (n=9), or septic shock (n=26); 40 mechanically ventilated; 12 on haemofiltration",
    dose_range     = "1 g or 2 g ceftriaxone IV infusion over 20 minutes, typically once daily (median daily dose 2 g)",
    regions        = "France (single-centre: CHRU de Tours)",
    saps_ii        = "mean 50 (range 9-87)",
    renal_function = "Measured creatinine clearance (24-h urine collection) median 68.5 mL/min (range 5.5-214); 12 patients on haemofiltration",
    sepsis_origin  = "Lung 33, urinary 10, intra-abdominal 6, skin/soft-tissue 2, CNS 1, ENT 1, undetermined 1",
    notes          = "Baseline demographics per Garot 2011 Table 1. Single-centre prospective study (July 2006 - March 2008) at CHRU Tours, France (NCT00449800). PK sampling on two occasions (PK1 day 2; PK2 day 5 or 48 h after catecholamine withdrawal for septic shock). 20 patients had full sampling (10 samples per occasion); 34 patients had sparse sampling (6 samples per occasion). 709 total ceftriaxone concentrations contributed to model building. The model did not support interoccasion variability."
  )

  ini({
    # Structural parameters at the typical septic patient (median CLcr 71 mL/min,
    # i.e. 4.26 L/h); Garot 2011 Table 3 final-estimate column.
    lcl <- log(0.56); label("Non-renal component of CL (intercept) (L/h)") # Garot 2011 Table 3: theta1 = 0.56 L/h (intercept of additive linear CL ~ CLcr term)
    lvc <- log(10.3); label("Central volume of distribution V1 (L)")      # Garot 2011 Table 3: V1 = 10.3 L
    lvp <- log(7.35); label("Peripheral volume of distribution V2 (L)")   # Garot 2011 Table 3: V2 = 7.35 L
    lq  <- log(5.28); label("Intercompartmental clearance Q (L/h)")       # Garot 2011 Table 3: Q = 5.28 L/h

    # Covariate effect: additive linear CL ~ CLcr with divisive normalization to
    # the population median CLcr expressed as 4.26 L/h (= 71 mL/min). Mirrors the
    # form CL = theta1 + theta2 * (CLcr/4.26) reported in Garot 2011 Table 3.
    # The CL slope is interpreted as L/h gained per unit of (CRCL / 71 mL/min);
    # at the median patient (CRCL = 71 mL/min) CL = 0.56 + 0.32 = 0.88 L/h, which
    # matches the population-mean CL of 0.88 L/h reported in the Results.
    e_crcl_cl <- 0.32; label("Renal CL slope per (CRCL / 71 mL/min) (L/h)") # Garot 2011 Table 3: theta2 = 0.32 L/h

    # Inter-individual variability (Garot 2011 Table 3 final-model BSV omega^2
    # column). The paper reports variance values 0.24, 0.23, 0.42 with the
    # small-omega approximation CV% = sqrt(omega^2)*100 in parentheses (49%, 48%,
    # 65%). The variance values are used directly; the corresponding CV under
    # the exact log-normal CV = sqrt(exp(omega^2)-1) is slightly higher (52%,
    # 51%, 75%) than the approximation the authors reported. No IIV on Q.
    etalcl ~ 0.24 # Garot 2011 Table 3: omega^2(CL) = 0.24 (49% CV approx)
    etalvc ~ 0.23 # Garot 2011 Table 3: omega^2(V1) = 0.23 (48% CV approx)
    etalvp ~ 0.42 # Garot 2011 Table 3: omega^2(V2) = 0.42 (65% CV approx)

    # Proportional residual error. Garot 2011 Table 3 reports sigma^2 proportional
    # (%) = 24, SE 0.009, 95% CI 21-28%. The "(%)" notation indicates the
    # estimate is reported as a percent CV (24% CV residual); converting back to
    # variance gives sigma^2 ~= 0.0576 (matched by the SE on sigma^2 = 0.009
    # propagating to a CV 95% CI of approximately 21-28%). In nlmixr2 the
    # proportional residual is encoded as the SD of relative error, i.e. propSd
    # = 0.24 maps to 24% CV.
    propSd <- 0.24; label("Proportional residual error (fraction)") # Garot 2011 Table 3: sigma^2 proportional (%) = 24 (CV = 24%)
  })
  model({
    # Individual PK parameters. CL has an additive linear covariate term on CRCL
    # with divisive normalization (TVCL = intercept + slope * CRCL / 71 mL/min),
    # wrapped in exp(etalcl) to give log-normal IIV. V1, V2 follow the standard
    # exp(lparam + etalparam) form. No IIV on Q.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL / 71)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
