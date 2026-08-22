Cendros_2025_enflicoxib <- function() {
  description <- paste(
    "Veterinary (dog). Joint parent + metabolite population PK model for",
    "oral enflicoxib (a COX-2 selective NSAID dosed once weekly) and its",
    "active pyrazol metabolite in dogs. Enflicoxib is described by a",
    "two-compartment model with first-order absorption and an absorption",
    "lag time; relative bioavailability F is not identifiable from oral-only",
    "data and is FIXED at 1 while still carrying inter-individual",
    "variability. Parent elimination is split into two parallel apparent",
    "clearance arms: CL2/F, the irreversible biotransformation clearance",
    "that forms the pyrazol metabolite, and CL1/F, clearance through other",
    "elimination pathways. The pyrazol metabolite is described by a",
    "three-compartment model (central plus a shallow and a deep peripheral",
    "compartment) with first-order elimination, fed by the CL2/F formation",
    "flux out of the parent central compartment. Total body weight is the",
    "only retained covariate and enters every clearance and volume as a",
    "fixed-exponent allometric power of (WGT / 9.9 kg), with exponent 0.75",
    "on clearances and 1.0 on volumes. The metabolite central volume VM/F",
    "shares the same THETA as the parent central volume V/F. The structural",
    "parameters, allometry, inter-individual variability and residual error",
    "were estimated in healthy Beagle dogs by Cendros 2022; Cendros 2025",
    "reproduces that parameter set (its Table 1) and externally validates it",
    "against sparse plasma samples from 83 client-owned dogs of any breed",
    "with naturally occurring osteoarthritis treated weekly for 6 months.",
    "No covariate other than body weight influenced the PK, and no",
    "time-dependent PK or over-accumulation was observed.")

  reference <- paste(
    "Cendros JM, Salichs M, Encina G, Vela JM, Homedes J. Enflicoxib for the",
    "long-term management of canine osteoarthritis - External validation of a",
    "population pharmacokinetic model in dogs with osteoarthritis.",
    "Front Vet Sci. 2025;12:1645857. doi:10.3389/fvets.2025.1645857.",
    "Parameter estimates reproduced in Cendros 2025 Table 1 originate from:",
    "Cendros JM, Salichs M, Encina G, Vela JM, Homedes JM. Pharmacology of",
    "enflicoxib, a new coxib drug: efficacy and dose determination by clinical",
    "and PK-guided approach for the treatment of osteoarthritis in dogs based",
    "on an acute arthritis induction model. Vet Med Sci. 2022;8:31-45.",
    "doi:10.1002/vms3.670.",
    sep = " "
  )

  vignette <- "Cendros_2025_enflicoxib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot               = list(analyte = "enflicoxib", units = "mg", specimen = "administration site", verified = FALSE),
    central             = list(analyte = "enflicoxib", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1         = list(analyte = "enflicoxib", units = "mg", specimen = "plasma", verified = FALSE),
    central_pyrazol     = list(analyte = "enflicoxib pyrazol metabolite", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1_pyrazol = list(analyte = "enflicoxib pyrazol metabolite", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2_pyrazol = list(analyte = "enflicoxib pyrazol metabolite", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight of the dog.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the model. Enters every apparent",
        "clearance and every apparent volume as a fixed-exponent allometric",
        "power of (WGT / median(WGT)) per the Cendros 2025 allometric equation",
        "(P_POP = theta1 * (WGT / median(WGT))^theta2). median(WGT) = 9.9 kg is",
        "the median body weight of the healthy Beagle population in which the",
        "model was estimated; it is read off the explicit '(WGT/9.9)' terms in",
        "every Table 1 parameter equation. Exponents are 0.75 on clearances and",
        "1.0 on volumes and were fixed a priori, not estimated (Cendros 2025",
        "Discussion: 'already incorporates allometric exponents in clearances and",
        "volumes of distribution (0.75 and 1.00, respectively)'). Body weight was",
        "the only covariate showing a clear trend with plasma levels in the",
        "external-validation cohort (Cendros 2025 Results, 'Relationship between",
        "plasma levels and covariates')."),
      source_name        = "WGT"
    )
  )

  # Screened in the Cendros 2025 exploratory covariate analysis but NOT
  # retained in the model: "No relationship or very weak correlations were
  # observed for all covariates, except for body weight" (Results,
  # 'Relationship between plasma levels and covariates'). No point estimates
  # are reported for any of these, so they are documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Dog age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; no relationship with enflicoxib or pyrazol metabolite plasma levels. Cendros 2025 Discussion: 'In agreement with other NSAIDS of the same class, no effect has been seen for age or sex'."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical covariate (35/83 = 42.2% female); not retained."
    )
  )

  population <- list(
    species        = "dog (Beagle for parameter estimation; client-owned dogs of any breed for external validation)",
    n_subjects     = 83L,
    n_studies      = 1L,
    age_range      = "2-16 years (mean +/- SD 8.7 +/- 3)",
    weight_range   = "4.9-64.9 kg (mean +/- SD 27.0 +/- 15)",
    weight_median  = "9.9 kg in the healthy Beagle population used for parameter estimation (the allometric reference)",
    sex_female_pct = 42.2,
    disease_state  = "naturally occurring osteoarthritis with clinical signs (pain and lameness) for at least 3 weeks plus radiographic evidence in at least one pelvic or thoracic limb joint; baseline clinical sum score >= 4",
    dose_range     = "oral Daxocox tablets, 8 mg/kg loading dose on day 0 then 4 mg/kg once weekly for 26 weeks (27 administrations); actual mean administered doses were 10.4 mg/kg loading and 5.2 mg/kg maintenance",
    regions        = "Portugal and Hungary",
    n_observations = "142 plasma samples (75 on day 44, 2 on day 93, 65 on day 189); 2 samples per dog by design",
    breeds         = "42 purebred (50.6%) and 41 mixed-bred (49.4%), more than 25 breeds represented; Labrador Retriever and German Shepherd most frequent",
    notes          = paste(
      "Two-tier population. The structural model, parameter estimates, IIV and",
      "residual error in ini() were estimated by Cendros 2022 in young healthy",
      "Beagle dogs (hence the 9.9 kg allometric reference). Cendros 2025",
      "externally validated that parameter set, unchanged, against the sparse",
      "field-study cohort described by the demographics above, using VPC / pcVPC",
      "/ NPDE and a maximum a posteriori Bayesian (POSTHOC, MAXEVAL=0) fit; no",
      "parameter was re-estimated. Concentrations below the limit of",
      "quantification (5.0 ng/mL enflicoxib, 2.5 ng/mL pyrazol metabolite) were",
      "handled by the M3 likelihood method. Dosing was with food, which",
      "increases absorption.")
  )

  ini({
    # ------------------------------------------------------------------
    # All fixed effects are the final population estimates in Cendros 2025
    # Table 1 ("Pharmacokinetic parameter estimates and bootstrap results
    # for the popPK model of enflicoxib and pyrazol metabolite (3)"), which
    # reproduces the Cendros 2022 final model. Each Table 1 row is written
    # as "<PARAM>/F = theta_n * (WGT/9.9)^exponent" with the theta value on
    # the following row; the theta value is what is carried here and the
    # allometric term is applied in model().
    #
    # Parameter roles are taken verbatim from the Table 1 footnotes:
    #   dagger  (parent enflicoxib): relative bioavailability (F); absorption
    #     rate (Ka); lag time (Tlag); apparent volume of distribution in the
    #     central compartment (V/F); apparent volume of distribution in the
    #     peripheral compartment (V5/F); apparent clearance to other
    #     elimination pathways (CL1/F); apparent clearance of irreversible
    #     biotransformation to pyrazol metabolite (CL2/F); and apparent
    #     distribution clearance between central and peripheral compartment
    #     (CL5/F).
    #   double-dagger (pyrazol metabolite): apparent volume of distribution
    #     in the central compartment (VM/F); in the shallow peripheral
    #     compartment (VMP/F); in the deep peripheral compartment (VM2/F);
    #     apparent clearance (CLM/F); apparent distribution clearance between
    #     central and shallow peripheral compartments (CLMP/F); and between
    #     central and deep peripheral compartments (CLM2/F).
    # ------------------------------------------------------------------

    # --- Enflicoxib (parent) absorption ---
    lfdepot      <- fixed(log(1))      ; label("Log relative oral bioavailability F of enflicoxib (unitless); not identifiable from oral-only data, anchored at 1, but still carries IIV") # Cendros 2025 Table 1 (enflicoxib): "Relative bioavailability (F) ... 1 FIX", IIV 39%
    ltlag        <- log(0.249)         ; label("Log absorption lag time Tlag of enflicoxib (h)")                                     # Cendros 2025 Table 1 (enflicoxib): Lag Time (Tlag) = 0.249 hr
    lka          <- log(0.220)         ; label("Log first-order absorption rate constant ka of enflicoxib (1/h)")                    # Cendros 2025 Table 1 (enflicoxib): Absorption rate (ka) = 0.220 hr^-1

    # --- Enflicoxib (parent) disposition ---
    lvc          <- log(6.59)          ; label("Log apparent enflicoxib central volume V/F at the 9.9 kg allometric reference (L)")  # Cendros 2025 Table 1 (enflicoxib): V/F = theta4 * (WGT/9.9)^1.0, theta4 = 6.59 L
    lvp          <- log(27.70)         ; label("Log apparent enflicoxib peripheral volume V5/F at the 9.9 kg allometric reference (L)") # Cendros 2025 Table 1 (enflicoxib): V5/F = theta11 * (WGT/9.9)^1.0, theta11 = 27.70 L
    lq           <- log(12.1)          ; label("Log apparent enflicoxib inter-compartmental clearance CL5/F at the 9.9 kg allometric reference (L/h)") # Cendros 2025 Table 1 (enflicoxib): CL5/F = theta10 * (WGT/9.9)^0.75, theta10 = 12.1 L/hr

    # --- Enflicoxib (parent) elimination: two parallel apparent arms ---
    # Table 1 splits the parent's total apparent clearance into a
    # metabolite-formation arm (CL2/F) and a non-formation arm (CL1/F);
    # the parent's total loss is (cl_met + cl_nonmet) / vc * central and
    # the metabolite's input is cl_met / vc * central. This is the
    # lcl_met / lcl_nonmet canonical decomposition (cf. Lehr_2010_tesofensine.R).
    lcl_met      <- log(0.520)         ; label("Log apparent enflicoxib biotransformation clearance CL2/F to the pyrazol metabolite at the 9.9 kg allometric reference (L/h)") # Cendros 2025 Table 1 (enflicoxib), Elimination: CL2/F = theta6 * (WGT/9.9)^0.75, theta6 = 0.520 L/hr; footnote dagger "apparent clearance of irreversible biotransformation to pyrazol metabolite (CL2/F)"
    lcl_nonmet   <- log(0.193)         ; label("Log apparent enflicoxib clearance CL1/F through elimination pathways other than pyrazol-metabolite formation at the 9.9 kg allometric reference (L/h)") # Cendros 2025 Table 1 (enflicoxib), Elimination: CL1/F = theta5 * (WGT/9.9)^0.75, theta5 = 0.193 L/hr; footnote dagger "apparent clearance to other elimination pathways (CL1/F)"

    # --- Pyrazol metabolite disposition ---
    # VM/F carries the SAME theta (theta4 = 6.59 L) as the parent V/F.
    # This is not a transcription slip: Table 1 lists thetas 1-13 with no
    # gaps once F / Tlag / ka take theta1-theta3, so a separate metabolite
    # central-volume theta does not exist. It is confirmed numerically -
    # the metabolite terminal half-life implied by VM/F = 6.59 L is 13.81 d,
    # matching the 13.8 d the paper reports (Results, "Enflicoxib and
    # pyrazol metabolite plasma levels"), and the metabolite half-life is
    # strongly sensitive to VM/F.
    lvc_pyrazol  <- log(6.59)          ; label("Log apparent pyrazol metabolite central volume VM/F at the 9.9 kg allometric reference (L); shares theta4 with the parent V/F") # Cendros 2025 Table 1 (pyrazol metabolite): VM/F = theta4 * (WGT/9.9)^1.0, theta4 = 6.59 L
    lvp_pyrazol  <- log(29.3)          ; label("Log apparent pyrazol metabolite shallow peripheral volume VMP/F at the 9.9 kg allometric reference (L)") # Cendros 2025 Table 1 (pyrazol metabolite): VMP/F = theta9 * (WGT/9.9)^1.0, theta9 = 29.3 (units column prints "-"; VMP is a volume so L, cf. the sibling theta13 row which prints L)
    lvp2_pyrazol <- log(13.4)          ; label("Log apparent pyrazol metabolite deep peripheral volume VM2/F at the 9.9 kg allometric reference (L)") # Cendros 2025 Table 1 (pyrazol metabolite): VM2/F = theta13 * (WGT/9.9)^1.0, theta13 = 13.4 L
    lq_pyrazol   <- log(13.2)          ; label("Log apparent pyrazol metabolite inter-compartmental clearance CLMP/F to the shallow peripheral compartment at the 9.9 kg allometric reference (L/h)") # Cendros 2025 Table 1 (pyrazol metabolite): CLMP/F = theta8 * (WGT/9.9)^0.75, theta8 = 13.2 L/hr
    lq2_pyrazol  <- log(0.131)         ; label("Log apparent pyrazol metabolite inter-compartmental clearance CLM2/F to the deep peripheral compartment at the 9.9 kg allometric reference (L/h)") # Cendros 2025 Table 1 (pyrazol metabolite): CLM2/F = theta12 * (WGT/9.9)^0.75, theta12 = 0.131 L/hr
    lcl_pyrazol  <- log(0.111)         ; label("Log apparent pyrazol metabolite elimination clearance CLM/F at the 9.9 kg allometric reference (L/h)") # Cendros 2025 Table 1 (pyrazol metabolite), Elimination: CLM/F = theta7 * (WGT/9.9)^0.75, theta7 = 0.111 L/hr

    # --- Allometric exponents (fixed a priori, not estimated) ---
    e_wt_cl      <- fixed(0.75)        ; label("Allometric exponent of (WT / 9.9 kg) on every apparent clearance (unitless)") # Cendros 2025 Table 1: every clearance row is written "(WGT/9.9)^0.75"; Discussion confirms the exponents were fixed a priori ("already incorporates allometric exponents in clearances and volumes of distribution (0.75 and 1.00, respectively)")
    e_wt_vc      <- fixed(1.0)         ; label("Allometric exponent of (WT / 9.9 kg) on every apparent volume (unitless)")    # Cendros 2025 Table 1: every volume row is written "(WGT/9.9)^1.0"; Discussion confirms the exponents were fixed a priori

    # ------------------------------------------------------------------
    # Inter-individual variability. Cendros 2025 reports IIV as a percent
    # CV in the Table 1 "IIV" columns and states the random effects are
    # log-normal (Methods, "Established pop PK model in healthy Beagle
    # dogs": P_i = P_POP * exp(eta_i), eta ~ N(0, omega^2), "The magnitude
    # of IIV is expressed as the coefficient of variation (%CV)"), so
    # omega^2 = log(1 + CV^2).
    #
    # No covariance / correlation between random effects is reported, so
    # the omega matrix is encoded as diagonal.
    # ------------------------------------------------------------------
    etalfdepot     ~ 0.1417   # Table 1 enflicoxib F      IIV  39% -> log(1 + 0.39^2)  = 0.1417
    etaltlag       ~ 0.3800   # Table 1 enflicoxib Tlag   IIV  68% -> log(1 + 0.68^2)  = 0.3800
    etalka         ~ 0.7034   # Table 1 enflicoxib ka     IIV 101% -> log(1 + 1.01^2)  = 0.7034
    etalvc         ~ 0.6331   # Table 1 enflicoxib V/F    IIV  94% -> log(1 + 0.94^2)  = 0.6331
    etalcl_met     ~ 0.2153   # Table 1 enflicoxib CL2/F  IIV  49% -> log(1 + 0.49^2)  = 0.2153
    etalcl_pyrazol ~ 0.0917   # Table 1 metabolite CLM/F  IIV  31% -> log(1 + 0.31^2)  = 0.0917

    # ------------------------------------------------------------------
    # Residual error. Cendros 2025 Methods gives the residual model as
    # ln(C_ij) = ln(C_pred,ij) + eps_ij (additive on the natural-log scale),
    # which is proportional error on the linear concentration scale in
    # nlmixr2, and reports its magnitude as a percent CV in the Table 1
    # "Residual variability" rows. The CV is carried directly as the
    # proportional SD (for these magnitudes the log-scale SD differs by
    # <3%: sqrt(log(1 + 0.34^2)) = 0.331 and sqrt(log(1 + 0.20^2)) = 0.198).
    # ------------------------------------------------------------------
    propSd         <- 0.34    ; label("Proportional residual error on enflicoxib plasma concentration (fraction)")           # Cendros 2025 Table 1 (enflicoxib): Residual variability = 34%
    propSd_pyrazol <- 0.20    ; label("Proportional residual error on pyrazol metabolite plasma concentration (fraction)")   # Cendros 2025 Table 1 (pyrazol metabolite): Residual variability = 20%
  })

  model({
    # Allometric scaling of every apparent clearance and volume on total
    # body weight, normalised to the 9.9 kg median weight of the healthy
    # Beagle population in which the parameters were estimated
    # (Cendros 2025: P_POP = theta1 * (WGT / median(WGT))^theta2, and every
    # Table 1 row spells the term out as (WGT/9.9)).
    allo_cl_factor <- (WT / 9.9)^e_wt_cl
    allo_v  <- (WT / 9.9)^e_wt_vc

    # --- Individual enflicoxib (parent) parameters ---
    fdepot     <- exp(lfdepot + etalfdepot)
    tlag       <- exp(ltlag + etaltlag)
    ka         <- exp(lka + etalka)
    vc         <- exp(lvc + etalvc) * allo_v
    vp         <- exp(lvp) * allo_v
    q          <- exp(lq) * allo_cl_factor
    cl_met     <- exp(lcl_met + etalcl_met) * allo_cl_factor
    cl_nonmet  <- exp(lcl_nonmet) * allo_cl_factor

    # --- Individual pyrazol metabolite parameters ---
    vc_pyrazol  <- exp(lvc_pyrazol) * allo_v
    vp_pyrazol  <- exp(lvp_pyrazol) * allo_v
    vp2_pyrazol <- exp(lvp2_pyrazol) * allo_v
    q_pyrazol   <- exp(lq_pyrazol) * allo_cl_factor
    q2_pyrazol  <- exp(lq2_pyrazol) * allo_cl_factor
    cl_pyrazol  <- exp(lcl_pyrazol + etalcl_pyrazol) * allo_cl_factor

    # --- Absorption ---
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # --- Enflicoxib: two-compartment with first-order absorption ---
    # Total parent loss is the sum of the two apparent elimination arms;
    # only the cl_met arm feeds the metabolite.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl_met + cl_nonmet) / vc * central -
                          q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # --- Pyrazol metabolite: three-compartment with first-order elimination ---
    # Input is the parent's irreversible biotransformation flux
    # cl_met / vc * central (mass flux; the apparent-clearance form absorbs
    # the molar-mass ratio and the metabolite's own availability, as every
    # metabolite parameter in Table 1 is likewise reported as an apparent
    # "/F" quantity).
    d/dt(central_pyrazol)     <-  cl_met / vc * central -
                                  cl_pyrazol / vc_pyrazol * central_pyrazol -
                                  q_pyrazol / vc_pyrazol * central_pyrazol +
                                  q_pyrazol / vp_pyrazol * peripheral1_pyrazol -
                                  q2_pyrazol / vc_pyrazol * central_pyrazol +
                                  q2_pyrazol / vp2_pyrazol * peripheral2_pyrazol
    d/dt(peripheral1_pyrazol) <-  q_pyrazol / vc_pyrazol * central_pyrazol -
                                  q_pyrazol / vp_pyrazol * peripheral1_pyrazol
    d/dt(peripheral2_pyrazol) <-  q2_pyrazol / vc_pyrazol * central_pyrazol -
                                  q2_pyrazol / vp2_pyrazol * peripheral2_pyrazol

    # Plasma concentrations. Doses are in mg and volumes in L, so
    # central / vc has units mg/L = ug/mL; multiply by 1000 to report in
    # ng/mL, the scale used throughout Cendros 2025.
    Cc         <- central         / vc         * 1000
    Cc_pyrazol <- central_pyrazol / vc_pyrazol * 1000

    Cc         ~ prop(propSd)
    Cc_pyrazol ~ prop(propSd_pyrazol)
  })
}
