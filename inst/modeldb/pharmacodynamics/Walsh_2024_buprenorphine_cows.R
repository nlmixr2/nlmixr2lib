Walsh_2024_buprenorphine_cows <- function() {
  description <- paste(
    "Bounded-integer maximum-inhibition (Imax) exposure-response model of",
    "the Clinical Opiate Withdrawal Scale (COWS) total score under",
    "subcutaneous long-acting buprenorphine (CAM2038 Q1W) in",
    "non-treatment-seeking adults with moderate-to-severe opioid use",
    "disorder. The driver is the time-matched buprenorphine plasma",
    "concentration CP_BPN_NGML. The COWS total score (0-44, 45 ordered",
    "categories) is modelled with the bounded-integer approach of Ueckert",
    "(2021) in the numerically-stable implementation the authors cite: the",
    "prediction lives on the scale of the standard-normal quantile",
    "function and is mapped onto the score range by the normal CDF.",
    "Structure (Supplementary Methods, 'Model development: COWS score',",
    "third and final variant):",
    "IPRED(Cp) = BASEb.int - (BASEb.int - LBASEb.int) * Cp^gamma /",
    "(IC50^gamma + Cp^gamma), with gamma fixed to 1 and IC90 estimated in",
    "place of IC50 (IC50 = IC90 / 9 for gamma = 1). The drug effect",
    "therefore drags the latent score from a pre-treatment baseline down",
    "to a lower asymptote rather than to zero. Back-transforming the two",
    "asymptotes reproduces clinically sensible scores: 45 * pnorm(-0.887)",
    "= 8.4 (mild withdrawal, pre-treatment) and 45 * pnorm(-2.11) = 0.78",
    "(essentially no withdrawal, fully blocked).",
    "Encoded here as probitNorm(addSd, 0, 45), the continuous relaxation",
    "of the bounded-integer likelihood -- the integer flooring step is not",
    "applied (see vignette Errata).",
    "The model has NO PK ODE and no compartments: it is a static algebraic",
    "exposure-response relationship, and CP_BPN_NGML is a required",
    "time-varying input covariate supplied per observation record. The",
    "CAM2038 popPK model is not reported in this paper and is not packaged",
    "in nlmixr2lib at extraction time.",
    sep = " "
  )
  reference <- paste(
    "Walsh SL, Comer SD, Aguiar Zdovc J, Sarr C, Bjornsson M,",
    "Strandgarden K, Hjelmstrom P, Tiberg F.",
    "Pharmacokinetic-pharmacodynamic analysis of drug liking blockade by",
    "buprenorphine subcutaneous depot (CAM2038) in participants with opioid",
    "use disorder.",
    "Neuropsychopharmacology. 2024;49(7):1050-1057.",
    "doi:10.1038/s41386-023-01793-z.",
    "Parameter estimates from Supplementary Table S4; structural equations",
    "from the Supplementary Appendix (41386_2023_1793_MOESM1_ESM.pdf).",
    "Bounded-integer model form:",
    "Ueckert S. Modeling composite assessment data using item response",
    "theory. CPT Pharmacometrics Syst Pharmacol. 2018;7(4):205-218;",
    "and the improved-stability implementation cited by Walsh 2024.",
    "Underlying phase 2 study NCT02611752, reported in",
    "Walsh SL, Comer SD, Lofwall MR, et al. JAMA Psychiatry.",
    "2017;74(9):894-902. doi:10.1001/jamapsychiatry.2017.1874.",
    "Sister endpoint models from the same paper:",
    "modellib('Walsh_2024_buprenorphine_drugLiking'),",
    "modellib('Walsh_2024_buprenorphine_desireToUse').",
    sep = " "
  )
  vignette <- "Walsh_2024_buprenorphine_opioidBlockade"

  depends <- c("CP_BPN_NGML")

  units <- list(
    time          = "day",
    dosing        = "N/A (PD-only; buprenorphine plasma concentration is a required input covariate)",
    concentration = "COWS total score (observation, 0-44 ordinal); ng/mL (CP_BPN_NGML input covariate)"
  )

  covariateData <- list(
    CP_BPN_NGML = list(
      description        = "Time-matched buprenorphine plasma concentration driving the Imax inhibition of the COWS total score. Time-varying; supplied per observation record in the event table rather than computed from a coupled PK model.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Supplementary Methods, 'Model development: COWS score': Cp is the time-matched BPN concentration. COWS was evaluated on the same day before the hydromorphone challenge, so all challenge levels contributed. Observed BPN plasma concentration range across the phase 2 study was 0.636-12.3 ng/mL; the paper's application simulation swept 0-10 ng/mL (Supplementary Methods, 'Model application'). Member of the canonical CP_<DRUG>_<UNITS> plasma-PD-driver family (siblings CP_MORPH_NGML, CP_OXY_NGML, CP_FBX_NGML). The upstream CAM2038 popPK model is not reported in this paper and is not packaged in nlmixr2lib; concentrations must be supplied externally.",
      source_name        = "Cp (Walsh 2024 Supplementary Methods, COWS Imax equation)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 1L,
    n_observations = 598L,
    age_range      = "18-54 years (median 36.0; mean 35.8, SD 9.1)",
    age_median     = "36.0 years",
    weight_range   = "53.0-110.0 kg (median 73.4; mean 75.9, SD 14.0)",
    weight_median  = "73.4 kg",
    sex_female_pct = 25.5,
    race_ethnicity = c(Black = 51, White = 47, Other = 2),
    disease_state  = "Non-treatment-seeking adults with moderate to severe opioid use disorder (DSM-5), physically dependent on opioids.",
    dose_range     = "CAM2038 (subcutaneous long-acting buprenorphine) 24 mg (n = 22) or 32 mg (n = 25) once weekly on days 0 and 7. Intramuscular hydromorphone challenges of 0, 6 or 18 mg; COWS was assessed before the challenge so all challenge levels contributed.",
    regions        = "USA (multisite: University of Kentucky; Columbia University / New York State Psychiatric Institute)",
    biomarkers     = "Clinical Opiate Withdrawal Scale (COWS) total score: 11 opioid-withdrawal signs each rated 0-4 by a trained observer, total range 0-44, where 5-12 represents mild symptoms.",
    notes          = "Baseline demographics from Supplementary Table S1; observation counts from Supplementary Table S2 (COWS: 598 PD records, 284 + 314 across the 24 mg and 32 mg arms). Study NCT02611752. No covariate analysis is reported for this endpoint. Estimation used Monte Carlo importance sampling (IMP) in NONMEM 7.5 with MATRIX=R standard errors."
  )

  ini({
    # =====================================================================
    # Structural parameters. Walsh 2024 Supplementary Table S4 ("Parameter
    # estimates of the final COWS score PK/PD model").
    #
    # The COWS total score is an integer 0-44 (45 ordered categories) and
    # is modelled with the BOUNDED-INTEGER approach: the prediction lives
    # on the "scale of quantile function of standard distribution" (the
    # standard-normal quantile scale, i.e. the probit scale) and is mapped
    # onto the observable score range by the normal CDF. Hence the `probit`
    # transform prefix on the two baseline parameters below.
    #
    # Three Imax variants were tested (Supplementary Methods, "Model
    # development: COWS score"). The FINAL model is the third:
    #
    #                                    (BASEb.int - LBASEb.int) * Cp^gamma
    #   IPRED(Cp) = BASEb.int - ( ------------------------------------------- )
    #                                       IC50^gamma + Cp^gamma
    #
    # so the drug drags the latent score from BASEb.int down towards
    # LBASEb.int (a lower asymptote), not towards zero. gamma was fixed to
    # 1 and "the IC90 was estimated instead of the IC50".
    #
    # Back-checks that validate this encoding:
    #   * IC50 = IC90 / 9^(1/gamma) = 0.109 / 9 = 0.0121 ng/mL, matching the
    #     IC50 = 0.012 ng/mL that Table S4 reports as a derived quantity
    #     (RSE "N/A"). Confirms gamma = 1.
    #   * 45 * pnorm(-0.887) = 8.44 -> pre-treatment COWS ~8, inside the
    #     5-12 "mild symptoms" band expected for this OUD cohort.
    #   * 45 * pnorm(-2.11)  = 0.78 -> fully-blocked COWS <1, consistent
    #     with the paper's simulation result that at 0.100 ng/mL "80.4% of
    #     the simulated patients had a COWS score <= 1" (Results).
    # =====================================================================

    # Baseline COWS on the probit (standard-normal quantile) scale, before
    # any buprenorphine exposure (Cp = 0).
    probitbase <- -0.887
    label("BASEb.int: baseline COWS total score on the probit (bounded-integer) scale (unitless)")
    # Walsh 2024 Table S4: BASEb.int = -0.887 (RSE 8.67%).

    # Lower asymptotic COWS on the same probit scale, reached after the
    # maximum buprenorphine-driven decrease from baseline.
    probitbase_low <- -2.11
    label("LBASEb.int: lower asymptotic COWS total score on the probit (bounded-integer) scale (unitless)")
    # Walsh 2024 Table S4: LBASEb.int = -2.11 (RSE 4.58%).

    # Buprenorphine plasma concentration producing 90% of the maximum
    # decrease. Estimated in place of IC50 by the authors.
    lic90 <- log(0.109)
    label("IC90: BPN plasma concentration giving 90% of the maximal COWS decrease (ng/mL)")
    # Walsh 2024 Table S4: IC90 = 0.109 ng/mL (RSE 43.1%).

    # Sigmoidicity (Hill) coefficient, fixed to 1 by the authors.
    hill <- fixed(1)
    label("gamma: sigmoidicity (Hill) coefficient of the Imax relationship (unitless)")
    # Walsh 2024 Supplementary Methods: "the sigmoidicity, gamma, was fixed to 1".

    # =====================================================================
    # Bounded-integer scaling parameter. In the bounded-integer model the
    # latent variable is normal with this SD on the probit scale; together
    # with the 45-category range it sets how sharply the latent prediction
    # maps onto discrete COWS scores. Encoded as the SD of the
    # probit-normal residual, i.e. probitNorm(addSd, 0, 45).
    # =====================================================================
    addSd <- 0.265
    label("SDb.int: bounded-integer scaling SD on the probit scale (unitless)")
    # Walsh 2024 Table S4: SDb.int = 0.265 (RSE 5.66%).

    # =====================================================================
    # Inter-individual variability. Supplementary Methods, "Model
    # development: COWS score": "The additive IIV on BASEb.int., exponential
    # IIV on LBASEb.int., and the correlation between these two parameters
    # was also estimated." The two IIVs therefore enter differently in
    # model(): additively on probitbase, multiplicatively (exp) on
    # probitbase_low.
    #
    # Table S4 reports the two IIV magnitudes and their CORRELATION (0.0203)
    # rather than a covariance, so the off-diagonal is reconstructed as
    #   cov = corr * omega_1 * omega_2 = 0.0203 * 0.454 * 0.0948 = 0.000874.
    # The correlation's RSE is 478%, i.e. it is not distinguishable from
    # zero; it is retained as reported rather than dropped. See vignette
    # "Assumptions and deviations".
    # =====================================================================
    etaprobitbase + etaprobitbase_low ~ c(0.454^2,
                                          0.0203 * 0.454 * 0.0948, 0.0948^2)
    # Walsh 2024 Table S4: IIV BASEb.int = 0.454 (RSE 15.2%); IIV LBASEb.int = 0.0948 (RSE 14.3%); correlation = 0.0203 (RSE 478%).
  })

  model({
    # -------------------------------------------------------------------
    # Individual parameters. Note the two DIFFERENT IIV forms the paper
    # specifies: additive on the baseline, exponential on the lower
    # asymptote. Exponential IIV on a negative parameter scales its
    # magnitude while preserving its sign.
    # -------------------------------------------------------------------
    pbase     <- probitbase + etaprobitbase
    pbase_low <- probitbase_low * exp(etaprobitbase_low)

    # -------------------------------------------------------------------
    # The authors estimated IC90; recover IC50 via the sigmoid-Imax
    # identity IC90 = IC50 * 9^(1/gamma). With gamma fixed to 1 this is
    # simply IC50 = IC90 / 9 = 0.0121 ng/mL (Table S4's derived value).
    # -------------------------------------------------------------------
    ic90 <- exp(lic90)
    ic50 <- ic90 / 9^(1 / hill)

    # -------------------------------------------------------------------
    # Imax exposure-response between the two probit-scale asymptotes. This
    # is the paper's IPRED, which lives "on the scale of quantile function
    # of standard distribution" -- NOT on the observable 0-44 COWS scale.
    # -------------------------------------------------------------------
    latentcows <- pbase -
      (pbase - pbase_low) * CP_BPN_NGML^hill / (ic50^hill + CP_BPN_NGML^hill)

    # -------------------------------------------------------------------
    # Map the latent probit-scale prediction onto the observable score
    # range: probitInv(x, 0, 45) == 45 * pnorm(x). This back-transform is
    # REQUIRED because rxode2's probitNorm(sd, low, hi) expects the
    # prediction on the untransformed (observed) scale and applies
    # probit((y - low) / (hi - low)) internally -- feeding it the latent
    # value directly puts the argument outside (0, 45) and silently yields
    # NA. Because probitInv is the exact inverse of that internal probit,
    # the round trip is lossless and the residual error lands on the
    # paper's quantile scale, which is the bounded-integer construction.
    #
    # Back-check: probitInv(-0.887, 0, 45) = 8.44 (pre-treatment COWS,
    # inside the paper's 5-12 "mild withdrawal" band) and
    # probitInv(-2.11, 0, 45) = 0.78 (fully blocked).
    # -------------------------------------------------------------------
    cows <- probitInv(latentcows, 0, 45)

    # -------------------------------------------------------------------
    # Bounded-integer observation. The published likelihood additionally
    # floors the mapped value to an integer; that step is not applied here
    # (see vignette Errata) -- downstream users reproducing the paper's
    # category proportions should apply floor() to the simulated score.
    # -------------------------------------------------------------------
    cows ~ probitNorm(addSd, 0, 45)
  })
}
