Beal_2001_iv1cmt_bql <- function() {
  description <- "Methodology reference. One-compartment IV-bolus generative toy model from Beal 2001 (the M1-M7 BQL methods paper). Typical CL is fixed at 0.693 and Vd at 1 by the author so the time unit equals one half-life; a unit dose is administered at time 0. Encoded as the SI1 population simulation model with geometric SD 0.2 on both CL and Vd and homogeneous residual variance 0.02."
  reference <- "Beal SL. Ways to fit a PK model with some data below the quantification limit. J Pharmacokinet Pharmacodyn. 2001 Oct;28(5):481-504. doi:10.1023/a:1012299115260. PMID 11768292."
  vignette <- "Beal_2001_iv1cmt_bql"
  units <- list(time = "half_life", dosing = "dose_unit", concentration = "dose_unit/vol_unit")

  covariateData <- list()

  population <- list(
    species        = "None (methodology paper; simulation-only toy model with no drug, no patients, no fitted estimates).",
    n_subjects     = 20L,
    n_studies      = 1L,
    disease_state  = "N/A (Monte Carlo simulation study; not a fit of any real molecule).",
    dose_range     = "Single unit-valued IV bolus dose at t = 0 (Beal 2001 section 2.2, page 486).",
    regions        = "N/A",
    scope_note     = paste(
      "Filed under inst/modeldb/pharmacokinetics/ (not specificDrugs/) at the operator's direction",
      "(sidecar zotero-077-beal_2001_unknown request-001 q1=B, response 2026-06-21).",
      "The paper is a methodology reference introducing the M1-M7 estimation methods for handling",
      "concentration measurements below the assay's quantification limit; it is not a fit of any",
      "specific drug. The typical values CL = 0.693 and Vd = 1 are scale-fixed teaching constants",
      "chosen by Beal so that time is measured in half-lives and the dose is dimensionless."
    ),
    notes          = paste(
      "Section 2.2 (page 486): 'The clearance CL is taken to be .693, so that the units of time",
      "may be regarded as half-lives, and the volume of distribution Vd is taken to be 1.",
      "A single unit-valued dose is given at time 0. For each individual, the observation times",
      "are taken to be .5, 1, 1.5, 2, 2.5, and 3.' Population data section 2.2 (page 488):",
      "'each data set is comprized of data from 20 individuals, where individual-specific values",
      "of CL and Vd are drawn independently from lognormal distributions with geometric means",
      ".693 and 1, respectively.' The paper studies two omega_CL sets (0.2 and 0.4); the primary",
      "set (omega_CL = omega_Vd = 0.2, page 488 first set) is encoded here. Residual variance follows",
      "SI1 (page 487): g(t) = 0.02 for all t (homogeneous variance, matched by the DA1 analytic",
      "model in section 2.3, page 489)."
    )
  )

  ini({
    lcl <- fixed(log(0.693)); label("Clearance (dose_unit / vol_unit / half_life)")  # Beal 2001 p. 486 section 2.2 (clearance fixed at .693 so time is in half-lives)
    lvc <- fixed(log(1));     label("Volume of distribution (vol_unit)")             # Beal 2001 p. 486 section 2.2 (Vd fixed at 1)

    etalcl ~ fixed(0.04)  # variance = 0.2^2; Beal 2001 p. 488 section 2.2 first omega set (geometric SD on CL = 0.2, no correlation)
    etalvc ~ fixed(0.04)  # variance = 0.2^2; Beal 2001 p. 488 section 2.2 first omega set (geometric SD on Vd = 0.2, no correlation)

    addSd <- fixed(sqrt(0.02)); label("Additive residual error (SI1 homogeneous variance g(t) = 0.02)")  # Beal 2001 p. 487 section 2.2 SI1 residual variance
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
