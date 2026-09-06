Nagy_2017_obiltoxaximab_survival <- function() {
  description <- paste(
    "Preclinical (NZW rabbit + cynomolgus macaque).",
    "Weibull cure-rate (mixture) survival model for obiltoxaximab treatment of inhalational",
    "anthrax, fit simultaneously to infected rabbit and macaque survival data and used to",
    "justify the 16 mg/kg human dose under the US FDA Animal Rule.",
    "The cure fraction psurv (probability of surviving to the Day-28 end of study) follows an",
    "Emax dose-response on the logit scale minus an exponential log10 prior-to-treatment (PTT)",
    "bacteremia penalty; the Weibull death-rate lambda for the non-cured fraction is log-linear in",
    "log10 PTT bacteremia. Species was investigated as a covariate but is not retained in the",
    "final model. The model is algebraic and deterministic (no ODE state, no IIV, no residual",
    "error): it returns the survivor function directly."
  )
  reference <- paste(
    "Nagy CF, Mondick J, Serbina N, Casey LS, Carpenter SE, French J, Guttendorf R.",
    "Animal-to-human dose translation of obiltoxaximab for treatment of inhalational anthrax",
    "under the US FDA animal rule.",
    "Clin Transl Sci. 2017;10(1):12-19. doi:10.1111/cts.12433.",
    "Equations are from Results, 'Animal survival modeling'; parameter estimates are from",
    "Supplementary Table S2 (file CTS-10-12-s004).",
    "The same fitted survival model is reported a second time, with its logit(psurv) equation",
    "printed in full, by Yamamoto BJ, Shadiack AM, Carpenter S, Sanford D, Henning LN,",
    "O'Connor E, Gonzales N, Mondick J, French J, Stark GV, Fisher AC, Casey LS, Serbina NV.",
    "Efficacy projection of obiltoxaximab for treatment of inhalational anthrax across a range",
    "of disease severity. Antimicrob Agents Chemother. 2016;60(10):5787-5795.",
    "doi:10.1128/AAC.00972-16, Materials and Methods, 'Survival modeling'.",
    "That printing is the source for the exp() wrapping the bacteremia term (see the model block).",
    "The companion population PK model from the same paper is available as",
    "modellib('Nagy_2017_obiltoxaximab').",
    sep = " "
  )
  vignette <- "Nagy_2017_obiltoxaximab"
  units <- list(
    time = "day",
    dosing = "mg/kg",
    concentration = "probability (the model output sur is the probability of surviving beyond time t, not a drug concentration)"
  )

  covariateData <- list(
    DOSE_OBILTOXAXIMAB_MGKG = list(
      description        = "Administered obiltoxaximab dose level, in mg/kg, as a per-subject constant",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the Emax term Emax * dose / (ED50 + dose) on logit(psurv).",
        "Set to 0 for placebo animals, which reduces the survival model to its",
        "untreated baseline. Dose levels in the survival data set were 0 (placebo),",
        "1, 4, 8, 16 and 32 mg/kg i.v. across rabbit Study 2 and macaque Studies 2-5",
        "(Nagy 2017 Table 1). Obiltoxaximab was given as a single i.v. dose triggered by",
        "a significant body-temperature rise (rabbits) or a positive serum",
        "protective-antigen signal (macaques), so a single dose level applies per animal.",
        "This is a covariate carrying the assigned DOSE LEVEL, not the rxode2 event",
        "column `amt` -- the model has no PK compartment."
      ),
      source_name        = "DOSE"
    ),
    BACT_PTT_LOG10CFU = list(
      description        = "Prior-to-treatment (PTT) quantitative bacteremia, expressed as log10 colony-forming units per mL of blood",
      units              = "log10 CFU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The paper's disease-severity covariate. Enters logit(psurv) through the",
        "exponential penalty exp((theta1 * BACT_PTT_LOG10CFU)^theta2) and enters the",
        "Weibull death rate through log(lambda) = lambda0 + lambda1 * BACT_PTT_LOG10CFU.",
        "Supplementary Figure S2 stratifies the observed and predicted dose-response by",
        "quartiles of this covariate: [BLQ, 3.02], [3.03, 3.95], [3.96, 4.87] and",
        "(4.87, 8.56] log10 CFU. A value of 0 represents no detectable PTT bacteremia,",
        "for which the paper reports a predicted 0.98 probability of cure at 16 mg/kg.",
        "Higher values mean more advanced disease at the time of treatment; the paper",
        "describes ~5 log10 CFU/mL as a 'point of no return' beyond which survival",
        "probability is near zero regardless of dose."
      ),
      source_name        = "PTT"
    )
  )

  population <- list(
    species        = "rabbit (New Zealand White) + cynomolgus macaque",
    n_subjects     = 261L,
    n_studies      = 5L,
    disease_state  = paste(
      "Inhalational anthrax following aerosol challenge with a target 200 LD50 of",
      "Bacillus anthracis (Ames strain) spores. Randomized, blinded, parallel-group,",
      "placebo-controlled trigger-to-treat studies conducted under GLP.",
      "Treatment was triggered by a significant increase in body temperature (rabbits)",
      "or a positive serum protective-antigen ECL signal (macaques); if neither occurred,",
      "obiltoxaximab was given at a predetermined time (rabbits 72 h, macaques 54 h",
      "post-challenge)."
    ),
    dose_range     = "Placebo, 1, 4, 8, 16, or 32 mg/kg obiltoxaximab as a single i.v. dose",
    notes          = paste(
      "Survival analyses used infected rabbit Study 2 (N = 70) and infected cynomolgus",
      "macaque Studies 2-5 (N = 44 + 48 + 48 + 51), giving n_subjects = 261 from the",
      "Table 1 enrolment counts; the paper does not report a pooled analysis-set size.",
      "The study end point is survival to Day 28 post-challenge.",
      "Species (rabbit vs macaque) was investigated as a covariate on the survival",
      "function (Methods) but does not appear in the final parameter table, so the",
      "final model pools the two species.",
      "The paper also developed an exposure (AUC)-response survival model as part of the",
      "initial analysis and reports that it led to identical inferences; only the",
      "dose-response model is parameterised in the supplement and therefore extracted here."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Weibull cure-rate model, Nagy 2017 Supplementary Table S2.
    #
    # Survivor function (Results, 'Animal survival modeling'):
    #   P(T > t) = psurv + (1 - psurv) * exp(-(lambda * t)^alpha)
    #
    #   logit(psurv) = theta0 - exp((theta1 * log10 PTT)^theta2)
    #                  + Emax * dose / (ED50 + dose)
    #   log(lambda)  = lambda0 + lambda1 * log10 PTT
    #
    # These three equations render as "formula-not-decoded" in the docling
    # conversion of the PDF; they were recovered verbatim with
    # `pdftotext -layout` from the on-disk PDF (page 4, right column).
    #
    # NOTE: Nagy's DISPLAYED logit(psurv) equation omits the exp() shown above,
    # but its own prose one paragraph earlier says the model used "an
    # exponential effect of log10 (PTT bacteremia) on logit(psurv)", and the
    # survivor function in that same display block is itself mis-typeset
    # (printed as "exp[-(lambda t)]alpha", with the Weibull shape exponent
    # stranded outside the bracket). The exp() form is taken from the second
    # published report of this same fit, Yamamoto 2016
    # doi:10.1128/AAC.00972-16; see the model block for the numerical evidence.
    #
    # No IIV and no residual error are reported (Table S2 lists point
    # estimates and Wald 95% CIs only), which is expected for a
    # population-level parametric time-to-event model.
    # ------------------------------------------------------------------

    theta0 <- 0.105; label("Baseline logit for P(cure) (unitless)")                                   # Table S2 theta0 = 0.105 (95% CI -0.961, 1.17)
    emax   <- 4.060; label("Maximum obiltoxaximab effect on P(cure) on the logit scale (unitless)")   # Table S2 Emax = 4.060 (95% CI 2.56, 5.56)
    ed50   <- 1.640; label("Obiltoxaximab dose giving half-maximal effect on P(cure) (mg/kg)")        # Table S2 ED50 = 1.640 (95% CI 0.515, 5.22)
    theta1 <- 0.296; label("Rate for log10(PTT bacteremia) on P(cure) (per log10 CFU)")               # Table S2 theta1 = 0.296 (95% CI 0.206, 0.386)
    theta2 <- 1.320; label("Exponent for the quantitative bacteremia effect (unitless)")              # Table S2 theta2 = 1.320 (95% CI 0.534, 2.11)
    lambda0 <- -2.240; label("Intercept of log Weibull death rate lambda (log 1/day)")                # Table S2 lambda0 = -2.240 (95% CI -2.5, -1.98)
    lambda1 <- 0.171;  label("Slope of log10(PTT bacteremia) on log Weibull death rate (per log10 CFU)") # Table S2 lambda1 = 0.171 (95% CI 0.166, 0.226)
    alpha   <- 2.830;  label("Weibull shape parameter (unitless)")                                    # Table S2 alpha = 2.830 (95% CI 2.52, 3.18)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Cure fraction: probability of surviving to the Day-28 end of
    #    study. Emax dose-response minus an exponential log10 PTT
    #    bacteremia penalty, both on the logit scale.
    # ------------------------------------------------------------------
    # The exp() around the bacteremia term is printed explicitly by Yamamoto
    # 2016 (doi:10.1128/AAC.00972-16), Materials and Methods, 'Survival
    # modeling', which reports this same fitted model:
    #   logit(psurv) = theta0 - exp[(theta1 x log10(PTT))^theta2]
    #                  + Emax x dose/(ED50 + dose)
    # and describes it twice in the same paragraph as "an exponential effect of
    # log10(PTT bacteremia) on logit(psurv)". Recovered with `pdftotext -layout`
    # (the docling conversion mangles the sub/superscripts).
    #
    # Nagy 2017's own displayed equation omits this exp(), so the two printings
    # of this single fit disagree. Yamamoto is followed because (a) BOTH papers'
    # prose says "exponential effect", (b) Nagy's display block mis-typesets the
    # adjacent survivor function too, and (c) with all eight Table S2 parameters
    # held fixed, the exp() form is closest on 12 of 13 published anchors
    # (probability-scale RMS 0.037 vs 0.334). See the vignette Errata.
    #
    # Without the exp() the model predicts 85-98% survival at every dose and
    # every bacteremia level -- including 85% survival for UNTREATED animals at
    # 10^6 CFU/mL, against 0-14.3% observed (Yamamoto 2016 Table 2) -- and so
    # contradicts the bacteremia "point of no return" that is the central claim
    # of both papers.
    bactTerm <- exp((theta1 * BACT_PTT_LOG10CFU)^theta2)
    doseTerm <- emax * DOSE_OBILTOXAXIMAB_MGKG / (ed50 + DOSE_OBILTOXAXIMAB_MGKG)

    logitPsurv <- theta0 - bactTerm + doseTerm
    psurv <- 1 / (1 + exp(-logitPsurv))

    # ------------------------------------------------------------------
    # 2. Weibull death rate for the non-cured fraction.
    # ------------------------------------------------------------------
    lambda <- exp(lambda0 + lambda1 * BACT_PTT_LOG10CFU)

    # ------------------------------------------------------------------
    # 3. Cure-rate (mixture) survivor function. As t grows the Weibull
    #    term decays to zero and sur approaches the cure fraction psurv,
    #    which is the quantity the paper's dose-response figures and ED50
    #    / ED90 statements refer to.
    # ------------------------------------------------------------------
    sur <- psurv + (1 - psurv) * exp(-(lambda * t)^alpha)

    # Cumulative hazard of the mixture, for users who want to drive a
    # time-to-event simulation from this model.
    cumhazard <- -log(sur)
  })
}
