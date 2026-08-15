Zhang_2025_cefiderocol <- function() {
  description <- "Two-compartment population PK model for cefiderocol in healthy Chinese adults after 2 g intravenous infusion over 3 h, parameterised by the hybrid disposition rate constants alpha and beta"
  reference <- paste(
    "Zhang C, Yu S, Li S, Wu X, Wei Q, He J, Cao G, Yang H, Wang J,",
    "Fujitani K, Katsube T, Zhang J, Dou H.",
    "Pharmacokinetic, Pharmacokinetic/Pharmacodynamic, and Safety",
    "Investigations of Cefiderocol in Chinese Healthy Subjects.",
    "Adv Ther. 2025;42(5):2285-2297. doi:10.1007/s12325-025-03147-1",
    sep = " "
  )
  vignette <- "Zhang_2025_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # The model carries no covariate effects: all 12 participants were healthy
  # Chinese adults with normal renal function studied at a single dose level,
  # and the source reports no covariate screen.
  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 12,
    n_studies = 1,
    age_range = "18-45 years",
    age_median = "24.8 years (mean)",
    weight_range = "47.6-76.2 kg",
    weight_median = "59.94 kg (mean)",
    sex_female_pct = 50,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "2 g intravenous infusion over 3 h; single dose on day 1 then q8h on days 2-5",
    regions = "China (single centre, Huashan Hospital, Shanghai)",
    renal_function = "normal; participants with clinically relevant renal disease were excluded",
    notes = paste(
      "Demographics from Zhang 2025 Table 1. Body mass index 19.1-25.4 kg/m2",
      "(mean 22.67). Eleven of 12 participants were ethnic Han.",
      "Open-label single-centre phase 1 study, ChiCTR2300076607."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural disposition parameters.
    #
    # The source fitted a two-compartment model with first-order elimination
    # in WinNonlin 8.1 and reports it in the MACRO-CONSTANT parameterisation:
    # the central volume V1 plus the transfer constant K21 and the two hybrid
    # disposition rate constants Alpha and Beta (Supplementary Table S2).
    # Those four numbers are transcribed verbatim below, and model() derives
    # kel, k12 and the clearance-parameterisation equivalents (cl, q, vp) from
    # them using the standard identities
    #
    #     alpha + beta  = kel + k12 + k21
    #     alpha * beta  = kel * k21
    #
    # so that every value in this ini() block appears literally in the source.
    # The derived typical values are cl = 5.197 L/h, q = 1.437 L/h,
    # vp = 3.422 L, kel = 0.5452 /h, k12 = 0.1508 /h; these reproduce the
    # published noncompartmental summaries closely (see the vignette source
    # trace) but are NOT themselves printed in the paper.
    # ---------------------------------------------------------------------
    lvc <- log(9.532); label("Central volume of distribution V1 (L)") # Zhang 2025 Supplementary Table S2, V1 = 9.532 (SD 1.455)
    lk21 <- log(0.420); label("Peripheral to central transfer rate constant k21 (1/h)") # Zhang 2025 Supplementary Table S2, K21 = 0.420 (SD 0.109)
    lalpha <- log(0.845); label("Hybrid distribution rate constant alpha of the two-compartment disposition (1/h)") # Zhang 2025 Supplementary Table S2, Alpha = 0.845 (SD 0.187)
    lbeta <- log(0.271); label("Hybrid terminal elimination rate constant beta of the two-compartment disposition (1/h)") # Zhang 2025 Supplementary Table S2, Beta = 0.271 (SD 0.039)

    # ---------------------------------------------------------------------
    # Interindividual variability.
    #
    # Supplementary Table S2 reports an arithmetic SD alongside each of the
    # four parameters above (n = 12): V1 1.455 L (CV 15.3%), K21 0.109 /h
    # (CV 26.0%), Alpha 0.187 /h (CV 22.1%), Beta 0.039 /h (CV 14.4%). Those
    # are MARGINAL spreads of the 12 individual WinNonlin estimates, not a
    # population omega, and the Monte Carlo simulation that produced the
    # published PTA/CFR results sampled them jointly while, in the authors'
    # words, "considering correlation among PK parameters" -- a correlation
    # matrix that is reported nowhere in the paper or its supplement.
    #
    # The marginals cannot be used on their own: a valid two-compartment
    # parameter set requires alpha > k21 > beta (otherwise
    # k12 = (alpha - k21)(k21 - beta)/k21 turns negative), and drawing the
    # three rate constants independently from the reported marginals violates
    # that ordering in roughly 9% of draws (8.6% reading the reported SDs as
    # lognormal, 12.4% as normal; see the vignette, which recomputes this). Encoding the marginals as
    # independent etas would therefore not reproduce the authors' simulation
    # and would generate structurally invalid subjects.
    #
    # These are declared as fixed(0) to record that the source DID report
    # between-subject spread on exactly these four parameters, while making
    # clear that no usable variance-covariance structure was published.
    # fixed(0) here means "joint distribution unreported", NOT "the authors
    # observed zero variability". Do not treat simulations from this model as
    # reproducing the published PTA or CFR spread without supplying an omega.
    # See the vignette Errata.
    # ---------------------------------------------------------------------
    etalvc ~ fixed(0) # Zhang 2025 Supplementary Table S2 (SD 1.455 L reported; correlation with the other parameters not reported)
    etalk21 ~ fixed(0) # Zhang 2025 Supplementary Table S2 (SD 0.109 /h reported; correlation with the other parameters not reported)
    etalalpha ~ fixed(0) # Zhang 2025 Supplementary Table S2 (SD 0.187 /h reported; correlation with the other parameters not reported)
    etalbeta ~ fixed(0) # Zhang 2025 Supplementary Table S2 (SD 0.039 /h reported; correlation with the other parameters not reported)

    # No residual error model was reported: the disposition parameters come
    # from an individual WinNonlin curve fit summarised as mean (SD), not from
    # a likelihood-based population fit, so no sigma exists to transcribe.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported by the source)") # Zhang 2025 - no residual error model reported
  })

  model({
    # 1. Individual parameters, in the source's macro-constant parameterisation.
    #    The IIV terms are declared by the source but carry zero variance here
    #    (see ini()).
    vc <- exp(lvc + etalvc)
    k21 <- exp(lk21 + etalk21)
    alpha <- exp(lalpha + etalalpha)
    beta <- exp(lbeta + etalbeta)

    # 2. Micro-constants, from the standard two-compartment identities
    #    alpha * beta = kel * k21 and alpha + beta = kel + k12 + k21.
    kel <- alpha * beta / k21
    k12 <- alpha + beta - k21 - kel

    # 3. Clearance-parameterisation equivalents, carried as derived outputs so
    #    that the familiar quantities are available without re-deriving them.
    cl <- kel * vc
    q <- k12 * vc
    vp <- q / k21

    # 4. ODE system. Intravenous administration into central; the source gave
    #    2 g as a 3-h infusion.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
