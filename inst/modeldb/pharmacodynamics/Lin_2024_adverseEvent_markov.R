Lin_2024_adverseEvent_markov <- function() {
  description <- paste(
    "Joint population PK and adverse-event (AE) grade model for a",
    "de-identified oral compound, used by Lin 2024 as the data-generating",
    "model for a simulation study of time-averaged exposure metrics (CavTE)",
    "in exposure-response analyses. PK is one-compartment with first-order",
    "absorption, an absorption lag and first-order elimination. AE grade",
    "(0 / 1 / 2) evolves as a first-order Markov chain whose per-visit",
    "transition probabilities follow a proportional-odds (cumulative-logit)",
    "model conditional on the previous grade. The cumulative logits are",
    "shifted by an Emax function of the time-averaged central-compartment",
    "concentration, with a larger maximum effect when the previous grade was",
    "0 (Emax0) than when it was >= 1 (Emax1). Grades 1 and 2 are pooled into",
    "an 'Any Grade' event in the source analysis. The compound is",
    "de-identified in the source, so the file stem names the endpoint rather",
    "than a drug (cf. Cardilin_2018_radiation_radiosensitizer_mouse).",
    sep = " "
  )
  reference <- paste(
    "Lin Y-W, Largajolli A, Edwards AY, Cheung SYA, Patel K, Hennig S.",
    "Impact of using time-averaged exposure metrics on binary endpoints in",
    "exposure-response analyses.",
    "Front Pharmacol. 2025;15:1487062.",
    "doi:10.3389/fphar.2024.1487062.",
    "All parameter values from Table 1 ('Parameter values used for the",
    "Example Scenario').",
    "Mixed-effects Markov + proportional-odds form adapted from",
    "Zingmark PH, Kagedal M, Karlsson MO.",
    "J Pharmacokinet Pharmacodyn. 2005;32(2):261-281,",
    "doi:10.1007/s10928-005-0021-7, as cited by the source Methods.",
    "Same cumulative-logit encoding as",
    "modellib('Hansson_2013_sunitinib_hfs').",
    sep = " "
  )
  vignette <- "Lin_2024_adverseEvent_markov"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # `auc` is a running exposure accumulator (not a physical compartment); it
  # integrates Cc so the model can form the time-averaged concentration that
  # drives the AE logits. Declared here to keep the compartment lint clean.
  paper_specific_compartments <- c("auc")

  compartmentData <- list(
    depot = list(
      analyte = "de-identified study drug", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "de-identified study drug", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    auc = list(
      analyte = "de-identified study drug", units = "ng/mL*h",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 200L,
    n_studies     = 1L,
    age_range     = "not reported (source cohort is de-identified and modified)",
    weight_range  = "not reported (source cohort is de-identified and modified)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state = paste(
      "Not disclosed. The source states the model is 'based on a real data",
      "example' from 'a clinical trial dataset, de-identified and modified",
      "for ethical and confidentiality reasons'; neither the drug, the",
      "indication nor the baseline demographics are reported.",
      sep = " "
    ),
    dose_range    = "60 mg orally once daily for four 28-day cycles (112 days total).",
    regions       = "not reported",
    biomarkers    = paste(
      "Adverse-event grade (0 = no event, 1, 2). Grades 1 and 2 are pooled",
      "into an 'Any Grade' category and only the first event per subject is",
      "retained for the exposure-response analysis; Grade 0 subjects are",
      "censored.",
      sep = " "
    ),
    notes = paste(
      "Virtual populations of n = 50, 100 or 200 were simulated (Methods 2.1);",
      "n_subjects records the largest. The observed 'Any Grade' incidence",
      "rate for each simulated scenario is reported in Supplementary Table S1",
      "and is used as the validation answer key in the vignette.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # PK -- Table 1, 'Parameters for the PopPK and proportional odds
    # models' block. One-compartment, first-order absorption with lag,
    # first-order elimination (Methods 2.1 + Figure 1A).
    # ------------------------------------------------------------------
    lka   <- log(4.23)  ; label("Absorption rate constant (1/h)")              # Table 1 Ka = 4.23 /h
    lcl   <- log(17.7)  ; label("Apparent clearance (L/h)")                    # Table 1 CL = 17.7 L/h
    lvc   <- log(229)   ; label("Apparent central volume of distribution (L)") # Table 1 V = 229 L
    ltlag <- log(0.154) ; label("Absorption lag time (h)")                     # Table 1 Alag = 0.154 h

    # IIV reported in Table 1 as CV%; omega^2 = log(CV^2 + 1).
    etalcl ~ 0.2558818  # Table 1 CL IIV 54.0% CV -> log(1 + 0.540^2)
    etalvc ~ 0.1100026  # Table 1 V  IIV 34.1% CV -> log(1 + 0.341^2)
    etalka ~ 0.6481629  # Table 1 Ka IIV 95.5% CV -> log(1 + 0.955^2)
    # Table 1 gives no IIV for Alag ('--'), so ltlag carries no eta.

    # ------------------------------------------------------------------
    # Proportional-odds Markov AE model -- Table 1 transition parameters.
    #
    # Table 1 reports, per previous grade, a pair (B<prev>1, B<prev>2).
    # These are the Zingmark 2005 / Hansson 2013 cumulative-logit
    # parameterisation: the FIRST member is the baseline cumulative logit
    # for grade >= 1, and the SECOND is a NEGATIVE INCREMENT added to it
    # to give the cumulative logit for grade >= 2. That guarantees
    # P(grade >= 1) >= P(grade >= 2) by construction.
    #
    # Reading the pair as two independent cumulative logits is impossible
    # for previous grade 0: B01 = -6.59 < B02 = -1.80 would make
    # P(>= 2) > P(>= 1). Three independent checks fix the reading above:
    #   (1) all three second members (-1.80, -6.70, -0.684) are negative,
    #       exactly as the increment parameterisation requires;
    #   (2) Hansson 2013 (cited by the same Zingmark lineage) reports its
    #       Table 3 in the identical B1 / B2 / B3 increment form -- see
    #       modellib('Hansson_2013_sunitinib_hfs');
    #   (3) only B01 = -6.59 as the grade >= 1 baseline reproduces the
    #       published incidence rates of Supplementary Table S1 (the
    #       alternative, -1.80, gives ~100% incidence at every drug-effect
    #       level against an observed 8.5%-99.5%). See vignette.
    # ------------------------------------------------------------------
    b01 <- -6.59  ; label("Baseline cumulative logit for AE grade >= 1 given previous grade 0")                  # Table 1 B01 = -6.59
    b02 <- -1.80  ; label("Increment to the cumulative logit for AE grade >= 2 given previous grade 0")          # Table 1 B02 = -1.80
    b11 <-  0.311 ; label("Baseline cumulative logit for AE grade >= 1 given previous grade 1")                  # Table 1 B11 = 0.311
    b12 <- -6.70  ; label("Increment to the cumulative logit for AE grade >= 2 given previous grade 1")          # Table 1 B12 = -6.70
    b21 <- -0.563 ; label("Baseline cumulative logit for AE grade >= 1 given previous grade 2")                  # Table 1 B21 = -0.563
    b22 <- -0.684 ; label("Increment to the cumulative logit for AE grade >= 2 given previous grade 2")          # Table 1 B22 = -0.684

    # Emax drug effect on the logit scale (Methods 2.1: "maximal drug effect
    # implemented on a logit scale, as described previously (Zingmark 2005)").
    # The Emax applied depends on the previous event grade.
    lemax_px0 <- log(4.73) ; label("Maximum drug effect on the AE cumulative logits when previous grade = 0 (log-odds units)")   # Table 1 Emax0 = 4.73
    lemax_px1 <- log(1.09) ; label("Maximum drug effect on the AE cumulative logits when previous grade >= 1 (log-odds units)")  # Table 1 Emax1 = 1.09
    lec50     <- log(6.05) ; label("Concentration at half-maximal drug effect on the AE logits (ng/mL)")        # Table 1 EC50 = 6.05 ng/mL

    # Table 1 reports 10.0% IIV (CV%) on every proportional-odds parameter.
    # No distribution or equation is given; encoded as exponential (log-normal)
    # IIV with omega^2 = log(1 + 0.10^2), applied multiplicatively so each
    # parameter keeps its sign (the b-parameters are logits and may be
    # negative). See vignette Assumptions and deviations.
    etab01       ~ 0.0099503  # Table 1 B01   IIV 10.0% CV -> log(1 + 0.10^2)
    etab02       ~ 0.0099503  # Table 1 B02   IIV 10.0% CV
    etab11       ~ 0.0099503  # Table 1 B11   IIV 10.0% CV
    etab12       ~ 0.0099503  # Table 1 B12   IIV 10.0% CV
    etab21       ~ 0.0099503  # Table 1 B21   IIV 10.0% CV
    etab22       ~ 0.0099503  # Table 1 B22   IIV 10.0% CV
    etalemax_px0 ~ 0.0099503  # Table 1 Emax0 IIV 10.0% CV
    etalemax_px1 ~ 0.0099503  # Table 1 Emax1 IIV 10.0% CV
    etalec50     ~ 0.0099503  # Table 1 EC50  IIV 10.0% CV

    # The source simulates without residual error (mrgsolve simulation study,
    # Methods 2.2) and reports no residual-error estimate for the PK.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # not reported: the source performs simulation only
  })

  model({
    # --- 1. Individual PK parameters -------------------------------------
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    tlag <- exp(ltlag)

    # --- 2. One-compartment PK with first-order absorption and lag -------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central
    alag(depot)   <- tlag

    # Table 1 gives EC50 in ng/mL while dose is in mg and vc in L, so the
    # amount/volume ratio (mg/L) is scaled by 1000 to reach ng/mL.
    Cc <- 1000 * central / vc

    # --- 3. Time-averaged central concentration driving the AE logits ----
    # Methods 2.1 states the logits are "positively dependent on
    # concentration in the central compartment" but does not state the
    # transition grid or whether the driver is instantaneous or averaged.
    # A weekly grid reproduces Supplementary Table S1 far better than any
    # other grid tried, whichever driver is used. Given a weekly grid the
    # driver is only weakly identified: the cumulative time-average used
    # here (the paper's own CavTE definition, Introduction) gives a mean
    # absolute error of ~1.4 percentage points over the six published
    # drug-effect levels, the mid-interval and interval-average
    # instantaneous concentrations give ~2.7 and ~1.5, and only the
    # dose-time trough is excluded (~9.6, under-predicting by up to 19
    # points). The cumulative average is preferred because it alone does
    # not depend on a within-day sampling phase the paper never states.
    # See vignette Assumptions and deviations.
    d/dt(auc) <- Cc
    cavg      <- auc / (t + 1e-10)

    # --- 4. Emax drug effect on the cumulative logits --------------------
    emax_px0 <- exp(lemax_px0 + etalemax_px0)
    emax_px1 <- exp(lemax_px1 + etalemax_px1)
    ec50     <- exp(lec50 + etalec50)

    eff_px0 <- emax_px0 * cavg / (ec50 + cavg)
    eff_px1 <- emax_px1 * cavg / (ec50 + cavg)

    # --- 5. Individual transition parameters (multiplicative 10% CV IIV) --
    i_b01 <- b01 * exp(etab01)
    i_b02 <- b02 * exp(etab02)
    i_b11 <- b11 * exp(etab11)
    i_b12 <- b12 * exp(etab12)
    i_b21 <- b21 * exp(etab21)
    i_b22 <- b22 * exp(etab22)

    # --- 6. Cumulative logits per previous grade -------------------------
    # Second Table 1 member is an increment, so the grade >= 2 logit is the
    # sum. Emax1 applies to both previous grade 1 and previous grade 2.
    clge1_px0 <- i_b01
    clge2_px0 <- i_b01 + i_b02
    clge1_px1 <- i_b11
    clge2_px1 <- i_b11 + i_b12
    clge1_px2 <- i_b21
    clge2_px2 <- i_b21 + i_b22

    pge1_px0 <- expit(clge1_px0 + eff_px0)
    pge2_px0 <- expit(clge2_px0 + eff_px0)
    pge1_px1 <- expit(clge1_px1 + eff_px1)
    pge2_px1 <- expit(clge2_px1 + eff_px1)
    pge1_px2 <- expit(clge1_px2 + eff_px1)
    pge2_px2 <- expit(clge2_px2 + eff_px1)

    # --- 7. Transition probabilities (each previous-grade row sums to 1) --
    p00 <- 1 - pge1_px0
    p01 <- pge1_px0 - pge2_px0
    p02 <- pge2_px0

    p10 <- 1 - pge1_px1
    p11 <- pge1_px1 - pge2_px1
    p12 <- pge2_px1

    p20 <- 1 - pge1_px2
    p21 <- pge1_px2 - pge2_px2
    p22 <- pge2_px2

    # --- 8. Observation --------------------------------------------------
    # The published likelihood is the joint Markov + proportional-odds form,
    # which rxode2 / nlmixr2 cannot express as a "previous DV"-conditioned
    # discrete likelihood (same deviation as
    # modellib('Hansson_2013_sunitinib_hfs')). The transition probabilities
    # above are emitted as deterministic outputs and the vignette advances
    # the Markov chain from them; the fitted observation is the PK.
    Cc ~ prop(propSd)
  })
}
