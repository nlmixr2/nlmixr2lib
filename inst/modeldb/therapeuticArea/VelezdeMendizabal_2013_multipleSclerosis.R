VelezdeMendizabal_2013_multipleSclerosis <- function() {
  description <- "Discrete count disease-dynamics model for the number of contrast-enhancing lesions (CELs) per monthly T1-weighted post-contrast MRI in relapsing-remitting multiple sclerosis. Negative binomial with first- and second-order Markov dependence on the previous monthly CEL counts ('NB nested MAK2 with steroid effects' from Velez de Mendizabal 2013): per-month expected CEL count lambda = lambda_0 + PDV * theta_pdv_eff + PPDV * theta_ppdv, where PDV and PPDV are the observed CEL counts at the previous and second-previous monthly MRIs and theta_pdv_eff is the typical theta_PDV (= 0.447) in months with no corticosteroid administration and theta_PDV_S (= 0.145, about 67 percent lower) in months in which the patient received corticosteroids for a clinical relapse. The source publication uses a negative-binomial likelihood with mean lambda and overdispersion OVDP (variance = lambda * (1 + OVDP * lambda)); following the ddmore/Plan_2012_pain.R and ddmore/Schoemaker_2018_levetiracetam.R precedents, the nlmixr2 observation is declared as Poisson(lambda) for fitting-API compatibility, and OVDP is exposed as a derived model variable so a downstream user can post-process Poisson samples through a negative-binomial correction if they need the full source dispersion. No drug PK is dosed in this model; the corticosteroid effect enters as a binary per-record covariate."
  reference <- paste(
    "Velez de Mendizabal N, Hutmacher MM, Troconiz IF, Goni J, Villoslada P,",
    "Bagnato F, Bies RR. (2013).",
    "Predicting Relapsing-Remitting Dynamics in Multiple Sclerosis Using",
    "Discrete Distribution Models: A Population Approach.",
    "PLoS ONE 8(9):e73361.",
    "doi:10.1371/journal.pone.0073361.",
    sep = " "
  )
  vignette <- "VelezdeMendizabal_2013_multipleSclerosis"
  units <- list(
    time          = "month",
    dosing        = "(none; corticosteroid effect enters as a binary per-record covariate)",
    concentration = "(count of contrast-enhancing lesions per monthly MRI; unitless integer)"
  )

  covariateData <- list(
    PDV = list(
      description        = "Observed CEL count at the previous monthly MRI; first-order Markov-state covariate for the negative-binomial count likelihood. Convention: PDV = 0 at the first monthly observation (no prior month) so the first-order Markov contribution is zero; for every subsequent observation, set PDV to the observed CEL count from the immediately preceding month.",
      units              = "(count, non-negative integer)",
      type               = "count",
      reference_category = NULL,
      notes              = "Source column PDV (canonical name matches the paper's PDV). Per the Velez de Mendizabal 2013 Methods section vi (Negative Binomial Markov elements) and equation 5, PDV 'takes the value of the previous dependent variable'. Enters the model as a linear coefficient on theta_pdv (or, during corticosteroid-treated months, on theta_pdv_s); see model() block. Existing canonical entry from inst/references/covariate-columns.md; PDV was originally registered for the Schoemaker 2018 levetiracetam seizure-count model and is reused here for monthly MS CEL counts (the canonical concept -- a per-record Markov-feedback observed count -- generalises across count-likelihood Markov-feedback PD models).",
      source_name        = "PDV"
    ),
    PPDV = list(
      description        = "Observed CEL count two monthly MRIs prior; second-order Markov-state covariate for the negative-binomial count likelihood. Convention: PPDV = 0 at the first two monthly observations (no second-prior month) so the second-order Markov contribution is zero; for every subsequent observation, set PPDV to the observed CEL count two months prior.",
      units              = "(count, non-negative integer)",
      type               = "count",
      reference_category = NULL,
      notes              = "Source column PPDV (canonical name matches the paper's PPDV). New canonical entry registered alongside this model in inst/references/covariate-columns.md, parallel to the existing PDV entry; PPDV captures the second-order Markov state of the count likelihood. Enters the model as a linear coefficient on theta_ppdv (equation 5). Velez de Mendizabal 2013 observed a decreasing magnitude pattern theta_PDV > theta_PPDV > theta_PPPDV across the first, second, and third Markov orders (Table 1), with the third-order fit improvement no longer statistically significant; the final selected model includes the first and second orders only.",
      source_name        = "PPDV"
    ),
    CONMED_STEROID = list(
      description        = "Indicator for systemic corticosteroid administration in the current monthly record. In the Velez de Mendizabal 2013 cohort, six of the nine patients received corticosteroids (intravenous methylprednisolone 1 g/day for 3-5 days or oral prednisone taper) for the treatment of clinical relapses; the column is 1 in months in which such a course was administered and 0 otherwise.",
      units              = "(binary, 0 / 1)",
      type               = "binary",
      reference_category = "0 (no corticosteroid administration this month)",
      notes              = "Source column STEROID. Canonical name CONMED_STEROID from inst/references/covariate-columns.md. The existing CONMED_STEROID register entry describes a baseline / time-fixed corticosteroid-use indicator (Narwal 2013 / Zheng 2016 sifalimumab SLE / asthma cohorts where systemic steroids are standard of care at study entry); in Velez de Mendizabal 2013 the same canonical concept is used time-varying (per monthly MRI record, on / off depending on whether the patient received a relapse-treatment course that month). The register entry's description was generalised to cover both temporal grains alongside this extraction. Effect: switches the first-order Markov coefficient from theta_pdv (no corticosteroids, 0.447) to theta_pdv_s (with corticosteroids, 0.145) in equation 5 -- about a 67 percent reduction; the source authors interpret this as steroids contributing to the inflammatory resolution of persistent (older) CELs without affecting the appearance of newly active lesions in the same month. Source paper Table 3, RSE 32.06 percent for theta_PDV_S. The authors evaluated whether the effect carried over to the immediately following month and found no significant lag-1 carry-over; that extension is not encoded here.",
      source_name        = "STEROID"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9L,
    n_studies      = 1L,
    age_range      = "(adult MS cohort; specific age range not tabulated in the paper -- refer to Bagnato et al. 2003 reference [25] for the full demographic detail)",
    age_median     = "(not reported in the Velez de Mendizabal 2013 publication)",
    weight_range   = "(not reported)",
    weight_median  = "(not reported)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Relapsing-remitting multiple sclerosis. Patients were immunomodulator- and immunosuppressant-naive at enrollment except for intravenous methylprednisolone (1 g/day for 3-5 days) or oral prednisone taper for clinical relapses; required to have been steroid-free for at least one month at study entry. Six of the nine model-building patients received corticosteroids during the 48-month observation window for relapse treatment.",
    dose_range     = "(no planned drug regimen; corticosteroids administered on a per-relapse basis only, encoded here as the binary CONMED_STEROID covariate)",
    regions        = "United States (NIH Bethesda, MD)",
    notes          = "Model-building cohort: n = 9 MS patients sequentially enrolled at the NIH Intramural Program (study approved by the Intramural Research Board of NINDS); patients underwent monthly 1.5 T T1-weighted post-contrast MRI for 48 months. The total CEL count per month was 'the sum of all the CELs that were enhancing at that month for the last time' (each CEL counted only once across the study; Methods, Patients and MRI Scans). External-validation cohort (Figure S1): n = 14 relapsing-remitting MS patients with monthly MRIs during a 6-month pre-therapy phase, none of whom received immunosuppressive therapy before the first scan; mean CEL count per patient per month 4.08 in the validation cohort vs 3.26 in the model-building cohort (Discussion). EDSS time-courses are shown alongside the CEL counts in Figure 1 but EDSS was not retained as a covariate in the final NB nested MAK2 model; the only retained covariate effect is the binary corticosteroid indicator on the first-order Markov term."
  )

  ini({
    # Parameter estimates from Velez de Mendizabal 2013 Table 3 ('NB nested
    # MAK2 with steroid effects model parameters'). All values are the
    # publication's final point estimates; the paper reports them on the
    # linear scale (lambda_0 = 0.923 etc.) and we apply log-transforms in
    # ini() for the log-normal IIV parameterization (positive-constrained
    # parameters; standard pharmacometric convention).

    llambda0 <- log(0.923)
    label("Log of typical baseline monthly CEL count rate lambda_0 (CELs/month) -- intercept of the equation 5 linearized Markov mean")
    # Table 3 lambda_0 = 0.923, RSE 26.54 percent.

    lovdp <- log(0.132)
    label("Log of typical negative-binomial overdispersion OVDP (unitless) -- exposed for downstream NB post-processing but not used by the nlmixr2 Poisson observation likelihood; see model() and the validation vignette's Assumptions and deviations section")
    # Table 3 OVDP = 0.132, RSE 25.37 percent. No IIV reported in Table 3.

    ltheta_pdv <- log(0.447)
    label("Log of typical first-order Markov coefficient theta_PDV (CELs added per CEL observed at the immediately preceding monthly MRI) in months with no corticosteroid administration")
    # Table 3 theta_PDV = 0.447, RSE 21.40 percent.

    ltheta_pdv_s <- log(0.145)
    label("Log of typical first-order Markov coefficient theta_PDV_S that replaces theta_PDV in months with corticosteroid administration (about 67 percent reduction relative to the no-steroid value)")
    # Table 3 theta_PDV_S = 0.145, RSE 32.06 percent. No IIV reported in
    # Table 3. The Discussion text quotes 'theta_PDV is diminished 66.44
    # percent when the patient was treated with immunosuppressive drugs';
    # (0.447 - 0.145) / 0.447 = 0.676 = 67.6 percent, consistent with the
    # quoted reduction within rounding.

    ltheta_ppdv <- log(0.150)
    label("Log of typical second-order Markov coefficient theta_PPDV (CELs added per CEL observed two monthly MRIs ago)")
    # Table 3 theta_PPDV = 0.150, RSE 48.00 percent. No IIV reported in
    # Table 3.

    # Inter-individual variability. Source Table 3 reports ISV(CV percent)
    # only for lambda_0 (66.18 percent) and theta_PDV (35.63 percent); the
    # other three parameters have no random effect in the published model.
    # Table 1 (same row, 'NB_nested MAK2 steroids') tabulates the omega
    # values 0.438 and 0.127 directly; the published CV percent maps to the
    # tabulated omega via CV percent = sqrt(omega^2) * 100 (the small-omega
    # approximation: sqrt(0.438) * 100 = 66.18 percent and sqrt(0.127) * 100
    # = 35.64 percent, both consistent with Table 3). Treating the tabulated
    # omega as a variance is the standard NONMEM $OMEGA convention and is
    # carried over here. Under the chosen log-normal IIV parameterization
    # the exact CV percent will be moderately larger than the published
    # sqrt approximation (sqrt(exp(0.438) - 1) * 100 = 74 percent for
    # lambda_0 and sqrt(exp(0.127) - 1) * 100 = 37 percent for theta_PDV);
    # the tabulated omega values are preserved verbatim because they
    # correspond directly to the NONMEM $OMEGA diagonal entries.
    etallambda0   ~ 0.438
    # Table 1 ('NB_nested MAK2 steroids' row) omega lambda = 0.438;
    # corresponds to Table 3 ISV(CV percent) lambda_0 = 66.18.

    etaltheta_pdv ~ 0.127
    # Table 1 ('NB_nested MAK2 steroids' row) omega PDV = 0.127;
    # corresponds to Table 3 ISV(CV percent) theta_PDV = 35.63.
  })

  model({
    # Subject-typical parameter back-transforms. Log-normal IIV on lambda_0
    # and theta_PDV (the only two parameters that carry a random effect in
    # the source); OVDP, theta_PDV_S, and theta_PPDV are population-typical
    # with no eta per Table 3.
    lambda0    <- exp(llambda0    + etallambda0)
    theta_pdv  <- exp(ltheta_pdv  + etaltheta_pdv)
    theta_pdv_s <- exp(ltheta_pdv_s)
    theta_ppdv <- exp(ltheta_ppdv)
    ovdp       <- exp(lovdp)

    # Corticosteroid switch on the first-order Markov coefficient. The
    # source paper replaces theta_PDV by the (smaller) theta_PDV_S in months
    # in which the patient received corticosteroids; Table 3 reports no IIV
    # on theta_PDV_S, so the steroid branch uses the population-typical
    # reduced value with no subject random effect. The literal binary
    # switch (1 - CONMED_STEROID) * theta_pdv + CONMED_STEROID * theta_pdv_s
    # is exact for CONMED_STEROID in {0, 1}.
    theta_pdv_eff <- (1 - CONMED_STEROID) * theta_pdv + CONMED_STEROID * theta_pdv_s

    # Equation 5 of the source paper, with the steroid covariate effect on
    # the first-order Markov term. lambda is the expected CEL count for the
    # current monthly MRI given the observed CEL counts from the previous
    # month (PDV), two months prior (PPDV), and the binary CONMED_STEROID
    # indicator for the current month. Every term on the right-hand side
    # is non-negative, so lambda > 0 and the Poisson likelihood is well
    # defined.
    lambda <- lambda0 + PDV * theta_pdv_eff + PPDV * theta_ppdv

    # Observation. The source publication uses a negative-binomial
    # likelihood with mean lambda and overdispersion OVDP (equation 10:
    # variance = lambda * (1 + OVDP * lambda)). Per the
    # ddmore/Plan_2012_pain.R and ddmore/Schoemaker_2018_levetiracetam.R
    # precedents for count-likelihood deviations, the nlmixr2 observation
    # is declared as Poisson(lambda) so the model fits with the standard
    # ini() / model() API. The deterministic typical-value trajectory of
    # lambda is unaffected by the choice of observation distribution; the
    # difference is only in the dispersion of stochastic VPC samples. OVDP
    # is emitted as a model variable above so a downstream user can
    # post-process Poisson samples to negative-binomial samples by drawing
    # N ~ NegBinom(lambda, OVDP) instead of N ~ Poisson(lambda). See the
    # validation vignette's Assumptions and deviations section for the
    # deviation rationale and a worked NB-correction example.
    cel_count <- lambda
    cel_count ~ pois(lambda)
  })
}
