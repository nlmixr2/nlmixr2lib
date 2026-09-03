Geldof_2008_fluvoxamine_sert_plasma_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Pharmacokinetic-pharmacodynamic model",
    "relating ex vivo 5-HT transporter (SERT) occupancy in rat frontal cortex",
    "to fluvoxamine PLASMA concentration, fit by Geldof et al. (2008, Br J",
    "Pharmacol 154:1369-1378). This is the first of the three PK-PD models",
    "reported side by side in Table 3 of that paper; the ECF- and brain-tissue-",
    "driven variants are Geldof_2008_fluvoxamine_sert_ecf_rat and",
    "Geldof_2008_fluvoxamine_sert_brain_rat. SERT occupancy was measured once",
    "per animal by ex vivo [3H]citalopram autoradiography after a single 30 min",
    "IV infusion of 1 or 7.3 mg/kg fluvoxamine. Occupancy is related to the",
    "driving concentration by a hyperbolic Bmax (Emax) function,",
    "B = Bmax * C / (EC50 + C), with no hysteresis (the paper found none, so no",
    "effect compartment is used). The plasma PK is a three-compartment linear",
    "disposition whose parameters are FIXED at the mean post-hoc estimates",
    "reported in Table 1 of the same paper (inherited from the upstream Geldof",
    "2007 rat population PK model) -- they were not re-estimated here.",
    "Inter-individual variability is on EC50 only; residual variability on",
    "occupancy is additive.",
    sep = " "
  )
  reference <- paste(
    "Geldof M, Freijer JI, van Beijsterveldt L, Langlois X, Danhof M.",
    "Pharmacokinetic-pharmacodynamic modelling of fluvoxamine 5-HT transporter",
    "occupancy in rat frontal cortex.",
    "Br J Pharmacol. 2008;154(6):1369-1378.",
    "doi:10.1038/bjp.2008.179. PMCID: PMC2483389.",
    "Plasma PK parameters fixed at the mean post-hoc estimates of the upstream",
    "Geldof 2007 rat population three-compartment PK model",
    "(Eur J Pharm Sci 2007;30(1):45-55; doi:10.1016/j.ejps.2006.10.001),",
    "as tabulated in Table 1 of the 2008 paper.",
    sep = " "
  )
  vignette <- "Geldof_2008_fluvoxamine_sert_occupancy"

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: analyte and specimen confirmed against
  # the source paper (Methods, "Data analysis"; Table 1 caption).
  compartmentData <- list(
    central     = list(analyte = "fluvoxamine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluvoxamine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "fluvoxamine", units = "ng", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Demographics and dose covariates were screened by the authors ("On the
  # basis of a covariate analysis, no differences in the PK of the different
  # dose groups or between the rats in the brain sampling studies and the
  # microdialysis studies could be detected", Results p1373; "No significant
  # correlation was observed between fluvoxamine dose and the PD parameter
  # estimates", Results p1374) but none was retained in the final model.
  covariatesDataExcluded <- list(
    DOSE_FLUVOXAMINE_MGKG = list(
      description = "Administered fluvoxamine dose level (1, 3.7 or 7.3 mg/kg, 30 min IV infusion)",
      units       = "mg/kg",
      type        = "continuous",
      notes       = paste(
        "Screened against the PD parameter estimates and against the plasma PK",
        "parameters; no significant correlation with either was found, so dose",
        "was not retained as a covariate (paper Results p1373-p1374)."
      )
    ),
    STUDY = list(
      description = "Source protocol (brain-sampling study vs microdialysis study)",
      units       = NA,
      type        = "categorical",
      notes       = paste(
        "Screened as a covariate on the plasma PK parameters; no difference",
        "between the two protocols could be detected (paper Results p1373).",
        "Table 1 nevertheless reports the per-protocol mean post-hoc estimates",
        "separately; the pooled 'Brain sampling + microdialysis' row is used here."
      )
    )
  )

  population <- list(
    species        = "rat (male Wistar, Charles River Wiga GmbH, Sulzfeld, Germany)",
    n_subjects     = 47L,
    n_studies      = 2L,
    age_range      = "adult (specific age not reported; group-housed 1 week on arrival, then 2 days individually after cannulation surgery)",
    weight_range   = "226-250 g body weight at the start of the experiments",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy rats prepared with a permanent right-jugular-vein cannula for",
      "fluvoxamine administration and a permanent left-femoral-artery cannula",
      "for serial arterial blood sampling. SERT occupancy in the frontal cortex",
      "was determined ex vivo by [3H]citalopram autoradiography on 20 um",
      "cryosections, expressed as a percentage of the labelling in the",
      "corresponding brain area of untreated control animals. Because only",
      "unoccupied transporters bind the radioligand, labelling is inversely",
      "related to occupancy."
    ),
    dose_range     = paste(
      "Single 30 min IV infusion of fluvoxamine free base into the right",
      "jugular vein at 1 mg/kg (24 rats) or 7.3 mg/kg (23 rats) in the",
      "brain-sampling study that provided the SERT occupancy observations",
      "(flow rate 20 uL/min, BAS BeeHive pump). The companion microdialysis",
      "study (26 rats) additionally used 3.7 mg/kg and contributed the brain",
      "ECF concentration data."
    ),
    regions        = "preclinical (in-vivo rat); Leiden University, The Netherlands",
    notes          = paste(
      "n_subjects = 47 counts the brain-sampling animals that contributed the",
      "ex vivo SERT occupancy observations used to fit the PD model (paper",
      "Figure 4 caption: 'the total population of 47 rats'). A further 26",
      "animals in the companion microdialysis study contributed brain ECF and",
      "plasma concentration data; n_studies = 2 counts both protocols.",
      "Only ONE SERT occupancy observation could be obtained per animal",
      "(destructive sampling), so inter- and intra-individual variability could",
      "not be separated -- the paper reports both an IIV term on EC50 and an",
      "additive residual error, and explicitly notes that no distinction",
      "between the two random effects is identifiable (Discussion p1377).",
      "Doses are reported per kg but the structural plasma PK parameters are",
      "absolute (mL, mL/min) and reflect the typical mid-range Wistar rat; the",
      "body-weight range is narrow (226-250 g) and weight was not a covariate.",
      "The LOQ for fluvoxamine was 1 ng/mL in plasma, brain ECF and brain",
      "tissue; the estimated plasma EC50 (0.48 ng/mL) lies BELOW that LOQ, so",
      "the lower limb of the plasma concentration-effect curve is informed by",
      "model-predicted rather than directly measured concentrations",
      "(paper Results p1374 and Figure 3a)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma PK -- three-compartment linear disposition. All six parameters
    # are FIXED at the mean post-hoc estimates reported in Geldof 2008 BJP
    # Table 1, row "Brain sampling + microdialysis" (the pooled cohort).
    # They originate in the upstream Geldof 2007 rat popPK model and were
    # NOT re-estimated in the present paper: the individual plasma
    # concentration at each occupancy sampling time was supplied to the PD
    # fit as a post-hoc empirical-Bayes prediction (paper Data analysis,
    # p1373). Population THETAs, IIV and residual error for the plasma PK
    # are therefore not available in this source and are not encoded.
    # ------------------------------------------------------------------
    lcl  <- fixed(log(29.6));  label("Plasma systemic clearance CL (mL/min)")                              # Geldof 2008 BJP Table 1, 'Brain sampling + microdialysis' row: CL = 29.6 mL/min
    lvc  <- fixed(log(294.4)); label("Plasma central volume of distribution V1 (mL)")                      # Geldof 2008 BJP Table 1, row 1: V1 = 294.4 mL
    lvp  <- fixed(log(858.1)); label("Plasma peripheral 1 volume of distribution V2 (mL)")                 # Geldof 2008 BJP Table 1, row 1: V2 = 858.1 mL
    lq   <- fixed(log(31.8));  label("Inter-compartmental clearance central <-> peripheral1 Q2 (mL/min)")  # Geldof 2008 BJP Table 1, row 1: Q2 = 31.8 mL/min
    lvp2 <- fixed(log(136.3)); label("Plasma peripheral 2 volume of distribution V3 (mL)")                 # Geldof 2008 BJP Table 1, row 1: V3 = 136.3 mL (IIV not identifiable upstream; the n = 187 population estimate is used, identical across all three Table 1 rows)
    lq2  <- fixed(log(1.0));   label("Inter-compartmental clearance central <-> peripheral2 Q3 (mL/min)")  # Geldof 2008 BJP Table 1, row 1: Q3 = 1.0 mL/min (IIV not identifiable upstream; identical across all three Table 1 rows)

    # ------------------------------------------------------------------
    # PD -- hyperbolic Bmax (Emax) model for SERT occupancy, estimated in
    # THIS paper (Geldof 2008 BJP Table 3, 'Plasma' column):
    #   B = Bmax * C / (EC50 + C)
    # with C = the plasma fluvoxamine concentration. The paper found no
    # hysteresis between plasma concentration and occupancy, so occupancy is
    # related directly to Cc with no effect compartment (Discussion p1376).
    # Bmax is expressed in percent occupancy; EC50 in ng/mL.
    # ------------------------------------------------------------------
    lemax <- log(94.5); label("Maximum SERT occupancy Bmax (% of control [3H]citalopram labelling)")  # Geldof 2008 BJP Table 3, Plasma column: Bmax = 94.5 % (CV 1.1%)
    lec50 <- log(0.48); label("Plasma fluvoxamine concentration at half-maximal SERT occupancy EC50 (ng/mL)")  # Geldof 2008 BJP Table 3, Plasma column: EC50 = 0.48 ng/mL (CV 11.6%)

    # ------------------------------------------------------------------
    # Inter-individual variability. The paper models IIV exponentially
    # (P_i = theta * exp(eta_i), Methods p1373) and reports the NONMEM
    # variance omega^2. IIV could be identified on EC50 only; it was not
    # estimable on Bmax and was fixed to zero there.
    # ------------------------------------------------------------------
    etalec50 ~ 0.34  # Geldof 2008 BJP Table 3, Plasma column: omega^2(eta_EC50) = 0.34 (CV 36.2%)

    # ------------------------------------------------------------------
    # Residual variability -- ADDITIVE on the percent-occupancy scale
    # (paper Methods p1373: Bo_ij = B_ij + eps1_ij, eps ~ N(0, sigma^2)).
    # Table 3 reports the VARIANCE sigma^2 = 33.5; nlmixr2 takes the SD:
    #   sqrt(33.5) = 5.78792 -> 5.7879 percentage points.
    # ------------------------------------------------------------------
    addSd <- 5.7879; label("Additive residual SD on SERT occupancy (percentage points)")  # Geldof 2008 BJP Table 3, Plasma column: sigma^2(eps1ij) = 33.5 (CV 27.2%); SD = sqrt(33.5) = 5.7879
  })

  model({
    # 1. Individual parameters
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    vp2  <- exp(lvp2)
    q2   <- exp(lq2)

    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system -- three-compartment plasma disposition. The dose is a
    #    30 min IV infusion into the central compartment (AMT + DUR/RATE in
    #    the event table).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 4. Observations
    Cc <- central / vc  # plasma fluvoxamine concentration (ng/mL); the PD driver. No residual error is reported for plasma PK in this paper, so Cc is a diagnostic output only.

    # Hyperbolic Bmax model (paper Methods p1373, unnumbered equation
    # immediately above the inter-individual variability equation).
    sertOccupancy <- emax * Cc / (ec50 + Cc)

    sertOccupancy ~ add(addSd)
  })
}
