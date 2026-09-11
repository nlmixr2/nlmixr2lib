Mahomed_2025_CAP256V2LS <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for the broadly neutralizing anti-HIV-1 monoclonal antibody CAP256V2LS in young South African women, supporting 1200 mg fixed versus 5-20 mg/kg weight-based dosing; recombinant human hyaluronidase (ENHANZE drug product) coadministration lowers clearance, central volume and bioavailability, and relative bioavailability is fixed to 1 on the second dosing occasion (Mahomed 2025)."
  reference <- "Mahomed S, Beliveau M, Heredia-Ortiz R, Osman F, Letsoalo M, Garrett N, Gengiah TN, Archary D, Wang J, Narpala S, Castro M, Serebryannyy L, Carlton K, Koup RA, Moore PL, Morris L, Abdool Karim Q, Abdool Karim SS. Population pharmacokinetics of weight-based compared with fixed dosing of CAP256V2LS, a broadly neutralizing antibody for HIV prevention in women. J Antimicrob Chemother. 2025; 80: 2135-2144. doi:10.1093/jac/dkaf181"
  vignette <- "Mahomed_2025_CAP256V2LS"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "CAP256V2LS", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "CAP256V2LS", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "CAP256V2LS", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_HYALURONIDASE = list(
      description        = "Coadministration of recombinant human hyaluronidase (rHuPH20; ENHANZE drug product, EDP) with the subcutaneous CAP256V2LS injection",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no EDP)",
      notes              = paste(
        "Time-fixed per dosing occasion. EDP was used to permit the full weight-based dose",
        "to be delivered as a single subcutaneous injection (Mahomed 2025 Methods, Trial design;",
        "Table S1 marks the EDP arms with an asterisk). 28/52 (53.8%) of the PK population",
        "received EDP (Table 3). Supplementary Table S6 carries three separate EDP effects:",
        "clearance x 0.378, central volume x 0.41 and an effective relative-bioavailability",
        "reduction factor of 0.452. Although hyaluronidase is a local-delivery adjuvant, the",
        "published model applies the clearance and volume effects to every record irrespective",
        "of route, and the paper's own IV +EDP simulation scenarios (Tables S3-S5) confirm that",
        "reading; the bioavailability effect applies only to subcutaneous (depot) doses.",
        "The trial administered EDP only with subcutaneous doses, so IV +EDP is a simulated",
        "hypothetical scenario rather than an observed one."
      ),
      source_name        = "EDP"
    ),
    OCC = list(
      description        = "Integer-valued dosing-occasion index (1 = first dose, 2 = second dose)",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (first dose)",
      notes              = paste(
        "Time-varying; constant within a dosing occasion. Only groups 2d and 2f received a",
        "second CAP256V2LS dose (16 or 24 weeks after the first; Table S1), so OCC = 2 exists",
        "only for those participants. Inter-occasion variability on F1 was tested and then",
        "replaced by a categorical occasion effect that improved the objective function",
        "substantially (Supplementary Methods, Structural model). Supply OCC >= 1 on every",
        "record; records with OCC >= 2 take relative bioavailability of 1."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened in the stepwise covariate analysis and NOT retained in the final model.",
        "Mahomed 2025 Discussion: 'body weight did not significantly impact the PK parameters",
        "of CAP256V2LS, despite the expected allometric relationship ... likely due to the",
        "relatively narrow range of body weight in the study population (an ~2-fold",
        "difference)'. Body weight still enters the paper's simulations through the",
        "weight-based mg/kg dose amount, not through any PK parameter."
      )
    ),
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis and not retained; Supplementary Results states 'No significant individual residual effects of covariates were observed on PK parameters, indicating that no additional covariate was required.'"
    ),
    CREAT = list(
      description = "Baseline serum creatinine concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Reported as a baseline PK-population characteristic (Table 3, median 58.0 umol/L) and screened but not retained in the final model."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase activity",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Reported as a baseline PK-population characteristic (Table 3, median 16.0 IU/L) and screened but not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 52L,
    n_studies      = 1L,
    age_range      = "18-43 years",
    age_median     = "24.5 years",
    weight_range   = "45.3-93.6 kg",
    weight_median  = "64.9 kg",
    sex_female_pct = 100.0,
    race_ethnicity = "Not reported by the source; the CAPRISA 012B trial enrolled young women recruited in Durban and the surrounding KwaZulu-Natal Province, South Africa.",
    disease_state  = "Healthy HIV-negative young women enrolled in the CAPRISA 012B first-in-human Phase 1 dose-escalation trial of the broadly neutralizing antibody CAP256V2LS for HIV prevention (PACTR202003767867253)",
    dose_range     = "5, 10 and 20 mg/kg IV or SC (weight-based, n = 44; two of the SC groups received a second dose at 16 or 24 weeks) and 1200 mg SC fixed dose (n = 8), given alone or sequentially with VRC07-523LS",
    regions        = "South Africa (CAPRISA eThekwini Clinical Research Site, Durban, KwaZulu-Natal)",
    route_split    = "IV 8 (15.4%); SC 44 (84.6%)",
    edp_split      = "EDP coadministered in 28 (53.8%); no EDP in 24 (46.2%)",
    sampling       = "767 CAP256V2LS PK sampling records, of which 119 were below the limit of quantitation and excluded, leaving 648 measurable observations. Plasma concentrations were measured by an electrochemiluminescence sandwich immunoassay on the Meso Scale Discovery platform.",
    notes          = "Demographics from Mahomed 2025 Table 3 (Participant characteristics, Overall N = 52 column). Groups 1c and 1d (HIV-positive participants dosed IV at 20 mg/kg) were still in follow-up and are not represented in this analysis (Table S1 footnote). No PK model was developed for VRC07-523LS in this paper."
  )

  ini({
    # Structural parameters -- Supplementary Table S6 (final model). Footnote a:
    # 'These model estimates were computed in log space'.
    lka <- log(0.033)    ; label("First-order subcutaneous absorption rate constant KA (1/h)")   # Table S6: KA = 0.033 1/h (RSE 13.8%)
    lcl <- log(0.00919)  ; label("Systemic clearance CL (L/h)")                                  # Table S6: CL = 0.00919 L/h (RSE 1.4%); main text Results reports 9.2 mL/h
    lvc <- log(4.76)     ; label("Central volume of distribution V2 (L)")                        # Table S6: V2 = 4.76 L (RSE 8.7%); main text Results reports 4.76 L
    lq  <- log(0.00694)  ; label("Inter-compartmental clearance Q (L/h)")                        # Table S6: Q = 0.00694 L/h (RSE 13%)
    lvp <- log(2.94)     ; label("Peripheral volume of distribution V3 (L)")                     # Table S6: V3 = 2.94 L (RSE 19.1%)

    # Relative bioavailability of the subcutaneous depot, estimated on the logit
    # scale so that it is bounded by 1 (Table S6 footnote b: 'Converted to
    # exp(x)/[1+exp(x)] since the actual value estimated x was fitted to be
    # bounded to 1').
    logitfdepot <- log(0.712 / (1 - 0.712))  ; label("Logit of subcutaneous relative bioavailability F1 on the first dosing occasion, no EDP (fraction)")  # Table S6: F1 = 0.712 (RSE 54.1%); logit(0.712) = 0.9050

    # Covariate effects.
    #
    # EDP (recombinant human hyaluronidase) effects on CL and V2 are reported as
    # multiplicative factors ('x if EDP') on the linear scale.
    e_edp_cl <- 0.378 ; label("Multiplicative effect of EDP coadministration on CL (unitless)")  # Table S6: 'x if EDP' 0.378 (RSE 12.3%); Supplementary Results, 'reduction in central clearance of CAP256V2LS (x0.38)'
    e_edp_vc <- 0.41  ; label("Multiplicative effect of EDP coadministration on V2 (unitless)")  # Table S6: 'x if EDP' 0.41 (RSE 31.4%); Supplementary Results, 'volume of distribution (0.41)'

    # The EDP effect on F1 is reported in Table S6 as an EFFECTIVE reduction
    # factor of 0.452 on the back-transformed fraction (Supplementary Results:
    # 'an effective reduction factor of 0.45 is obtained for F1'), i.e.
    # F1(EDP) = 0.712 * 0.452 = 0.3218. Because F1 itself is carried on the
    # logit scale, the coefficient stored here is the logit-scale difference
    # logit(0.712 * 0.452) - logit(0.712) = -1.6504, which reproduces exactly
    # that back-transformed fraction. Confirmed against the paper's own
    # simulations: Table S3 gives a single-dose SC AUC of 152.2 (no EDP) versus
    # 166.3 (EDP) ug/mL*week, which the multiplicative reading reproduces to
    # within 1% while a logit-additive reading of 0.452 would predict ~380.
    e_edp_fdepot <- log(0.712 * 0.452 / (1 - 0.712 * 0.452)) - log(0.712 / (1 - 0.712))  ; label("Effect of EDP coadministration on logit F1 (logit units)")  # Table S6: F1 'x if EDP' 0.452 (RSE 12.3%)

    # Occasion effect on F1: Table S6 reports F1 'if OCC=2' as 1, FIXED. The
    # underlying fitted quantity is a fixed additive offset on the logit scale;
    # the base model (Table S2) prints it as +20, which saturates
    # expit(logit(F1) + 20) to 1 to well beyond the reported precision. Encoded
    # with that offset so the back-transformed value equals the 1 that Table S6
    # reports. Confirmed against Table S3: the steady-state SC AUC equals
    # Dose/CL with F1 = 1 (217.7 vs 213.7 no-EDP; 567.8 vs 565.4 with EDP),
    # i.e. both the first-occasion F1 and the EDP factor are washed out on the
    # second and later occasions.
    e_occ2_fdepot <- fixed(20) ; label("Additive offset on logit F1 for dosing occasions 2 and later (logit units)")  # Table S6: F1 'if OCC=2' = 1 FIXED; Table S2 base model prints the logit-scale offset as '+ if OCC=2  20  FIXED'

    # Inter-individual variability. The final model carries IIV on CL and F1
    # only; the IIV on V2 quoted in the main text belongs to the base model
    # (Supplementary Results: 'The final population PK model included
    # inter-individual variability (IIV) in central parameters (CL and F1)').
    etalcl          ~ 0.019  # Abstract: 'Inter-individual variability in bioavailability and clearance was 0.212 and 0.019'; Table S6 IIV(CL) = 13.8% (shrinkage 38.2%), 0.138^2 = 0.0190
    etalogitfdepot  ~ 0.212  # Abstract, as above; Table S6 IIV(F1) = 46% (shrinkage 37.9%), 0.46^2 = 0.2116; carried on the logit scale because F1 is

    # Residual error: additive on log-transformed concentrations (Supplementary
    # Methods, Structural model: 'The structural population PK model included an
    # additive error model on log-transformed concentrations'), which is the
    # nlmixr2 lnorm() endpoint.
    expSd <- 0.196 ; label("Log-scale additive residual error SD (log ug/mL)")  # Table S6: Log-Additive Error = 0.196 (RSE 1.8%)
  })

  model({
    # 1. Derived covariate terms.
    occ2 <- (OCC >= 2)

    # 2. Individual parameters. EDP effects enter as the multiplicative factors
    #    printed in Table S6 ('x if EDP'), raised to the 0/1 indicator.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * e_edp_cl^CONMED_HYALURONIDASE
    vc <- exp(lvc)          * e_edp_vc^CONMED_HYALURONIDASE
    q  <- exp(lq)
    vp <- exp(lvp)

    fdepot <- expit(logitfdepot +
                      e_edp_fdepot  * CONMED_HYALURONIDASE +
                      e_occ2_fdepot * occ2 +
                      etalogitfdepot)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system. Subcutaneous doses go to 'depot'; intravenous doses go
    #    directly to 'central' and bypass both KA and F1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # 5. Bioavailability.
    f(depot) <- fdepot

    # 6. Observation and error. Dose in mg and vc in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
