# Joint population pharmacokinetic model of chloroquine and its active
# metabolite desethylchloroquine in adults with Plasmodium vivax
# mono-infection treated with the standard 3-day chloroquine regimen on
# the Thai-Myanmar border (Hoglund 2016, Malaria Journal 15:129;
# doi:10.1186/s12936-016-1181-1).

Hoglund_2016_chloroquine <- function() {
  description <- paste(
    "Joint parent + metabolite population PK model for oral chloroquine",
    "and its active metabolite desethylchloroquine in adults with",
    "Plasmodium vivax mono-infection treated with the standard 3-day",
    "25 mg base/kg chloroquine regimen (Hoglund 2016). A",
    "one-transit-compartment absorption chain with ktr = 2 / MTT feeds a",
    "two-compartment chloroquine disposition model; a fixed fraction",
    "fm = 0.18 of systemic chloroquine clearance is routed, with a molar",
    "correction, into a two-compartment desethylchloroquine disposition",
    "model, and the remaining 82% leaves as other elimination. No",
    "covariates were retained in the final model: body weight (allometric",
    "exponents fixed at 0.75 / 1.0), age, sex, parasite clearance time and",
    "fever clearance time were screened and dropped. Relative",
    "bioavailability is fixed at 1 and carries between-subject",
    "variability. NONMEM additive residual error on log-transformed",
    "observations is encoded here as a proportional residual in the linear",
    "concentration space. Table 1 of the source is internally inconsistent",
    "for two desethylchloroquine parameters; the values used here are",
    "reconstructed from the paper's own bootstrap CIs and reported",
    "terminal half-lives (see the ini() comments and the vignette Errata).",
    sep = " "
  )
  reference <- paste(
    "Hoglund R, Moussavi Y, Ruengweerayut R, Cheomung A, Abelo A,",
    "Na-Bangchang K (2016). Population pharmacokinetics of a three-day",
    "chloroquine treatment in patients with Plasmodium vivax infection on",
    "the Thai-Myanmar border. Malaria Journal 15:129.",
    "doi:10.1186/s12936-016-1181-1.",
    sep = " "
  )
  vignette <- "Hoglund_2016_chloroquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are in mg of the named analyte base: the
  # chloroquine states hold mg of chloroquine base, and the mass flux
  # entering central_dcq carries the mwDCQ / mwCQ molar correction so the
  # desethylchloroquine states hold mg of desethylchloroquine base.
  # Whole blood is the assayed matrix (Methods, Blood sampling and drug
  # analysis: "Whole blood concentrations of chloroquine and
  # desethylchloroquine were measured").
  compartmentData <- list(
    depot            = list(analyte = "chloroquine",          units = "mg", specimen = "administration site", verified = TRUE),
    transit1         = list(analyte = "chloroquine",          units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "chloroquine",          units = "mg", specimen = "whole blood",         verified = TRUE),
    peripheral1      = list(analyte = "chloroquine",          units = "mg", specimen = "whole blood",         verified = TRUE),
    central_dcq      = list(analyte = "desethylchloroquine",  units = "mg", specimen = "whole blood",         verified = TRUE),
    peripheral1_dcq  = list(analyte = "desethylchloroquine",  units = "mg", specimen = "whole blood",         verified = TRUE)
  )

  # No covariates are referenced in model(): the final model is
  # covariate-free (Results, Population pharmacokinetic models: "No
  # covariates were added in the final model").
  covariateData <- list()

  # Screened in the stepwise covariate search but NOT retained in the
  # final model (Methods, Population pharmacokinetics: "Relationships
  # between all parameters estimated in the base model and covariates,
  # i.e., body weight (BW), age, sex, parasite clearance time (PCT) and
  # fever clearance time (FCT) were evaluated"). Documentation only;
  # checkModelConventions() does not validate names that appear only
  # here. PCT and FCT are deliberately NOT registered as canonical
  # covariate columns -- nothing in the library carries them as a model
  # input, and registering an unused canonical fails the
  # structural-cleanliness check in checkNamingRegisters().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight. Applied as an allometric power model on all CL and V terms with exponents fixed at 0.75 (clearances) and 1.0 (volumes), normalised to the median body weight (Methods, Eq for Pt). Not retained in the final model.",
      units = "kg",
      type = "continuous",
      notes = "Screened as an allometric covariate on every clearance and volume; dropped in the backward elimination step."
    ),
    AGE = list(
      description = "Age at enrolment. Screened in the stepwise covariate search; not retained.",
      units = "years",
      type = "continuous",
      notes = "Cohort age range 17-52 years (Methods, Patients and study design)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male). Screened in the stepwise covariate search; not retained.",
      units = "unitless",
      type = "categorical",
      notes = "Cohort was 39 female / 36 male (Methods, Patients and study design)."
    ),
    PCT = list(
      description = "Parasite clearance time, the time taken for the parasite count to fall below the level of microscopic detection. Screened as a continuous covariate via the centred linear model Pt = theta1 * (1 + theta2 * (PCT - median PCT)); not retained.",
      units = "h",
      type = "continuous",
      notes = "Documentation-only key (not a registered canonical covariate column). Cohort median (95% CI) 30 (18-36) h (Results)."
    ),
    FCT = list(
      description = "Fever clearance time, the time taken for the temperature to return to normal (< 37.3 C). Screened as a continuous covariate via the centred linear model VPt = theta1 * (1 + theta2 * (FCT - median FCT)); significant on VP CQ/F in the forward addition step but NOT retained after backward elimination.",
      units = "h",
      type = "continuous",
      notes = "Documentation-only key (not a registered canonical covariate column). Cohort median (95% CI) 24 (12-42) h (Results). The Table 1 footnote still defines a leftover 'KFCT' symbol for this effect, but no KFCT row appears in Table 1 and the Results text states the effect was dropped in the backward step."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 75L,
    n_studies      = 1L,
    age_range      = "17-52 years (Methods, Patients and study design)",
    sex_female_pct = 52.0,
    race_ethnicity = "8 Thai and 67 Burmese migrant workers (Methods, Patients and study design)",
    disease_state  = paste(
      "Acute Plasmodium vivax mono-infection. Median (95% CI) admission",
      "parasitaemia 4898 (1206-29,480) parasites/uL. All 75 patients",
      "completed the 42-day follow-up with a 100% cure rate; neither",
      "recurrence of P. vivax parasitaemia nor appearance of P. falciparum",
      "occurred. Median (95% CI) parasite clearance time 30 (18-36) h and",
      "fever clearance time 24 (12-42) h (Results)."
    ),
    dose_range     = paste(
      "Standard 3-day chloroquine regimen, 25 mg base/kg body weight in",
      "total, given as 250 mg chloroquine phosphate tablets (Government",
      "Pharmaceutical Organization of Thailand): 10 mg base/kg at 0 h and",
      "5 mg base/kg at 6-12 h on day 0, then 5 mg base/kg on each of day 1",
      "and day 2. All doses were supervised and taken with 250 mL of water;",
      "patients were observed for at least 30 min after ingestion.",
      "Primaquine 15 mg base daily for 14 days was co-administered from",
      "day 1 onwards (Methods, Patients and study design)."
    ),
    regions        = "Mae Tao Clinic for migrant workers, Tak Province, Thailand (Thai-Myanmar border); samples collected during the 2010-2011 clinical efficacy study",
    notes          = paste(
      "Whole-blood chloroquine and desethylchloroquine quantified by HPLC",
      "with UV detection; LOQ 2 ng/mL for both analytes, assay accuracy",
      "0.25-5.7% relative error and precision < 5% CV. Samples below the",
      "LOQ (< 5% of samples) were excluded from the analysis. Sampling was",
      "pre-dose and at 1, 6, 12, 24, 25, 36, 48 and 49 h after the first",
      "dose, then on days 7, 14, 21, 28, 35 and 42. Estimation used FOCE in",
      "NONMEM 7.12 on natural-log-transformed concentrations. The number of",
      "observations is reported inconsistently by the source: the Abstract",
      "says 1045 observations from 75 participants while the Results say",
      "1405 -- the two differ by a digit transposition and the paper does",
      "not resolve which is correct."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters: Hoglund 2016 Table 1, "Final covariate
    # model (RSE)" column. All clearance and volume terms are apparent
    # (divided by bioavailability F) because "oral dosing was not
    # accompanied by an intravenous dose" (Methods, end of Population
    # pharmacokinetics). Reported on the linear scale; log() is applied
    # here for the nlmixr2 internal log scale.
    #
    # VERIFICATION OF THE WHOLE TABLE. Hoglund 2016 Table 1 also reports
    # the terminal half-lives computed from these same estimates
    # ("The pharmacokinetic population parameters estimated from the
    # final covariate model were used to calculate terminal half-life").
    # Those half-lives are an independent, non-circular check on the
    # four-parameter disposition set for each analyte, via the standard
    # two-compartment terminal slope
    #   beta = (a - sqrt(a^2 - 4*k10*k21)) / 2,  a = k10 + k12 + k21.
    # Chloroquine, using the Table 1 values verbatim, gives
    # t1/2 = 10.72 days against the reported 10.7 days -- so the
    # chloroquine row set is confirmed as printed, and CL CQ/F is
    # confirmed to be TOTAL apparent elimination clearance (it enters as
    # k10 = CL/Vc), not the formation clearance alone. See the
    # desethylchloroquine block below for the same check applied there.
    # ---------------------------------------------------------------

    lmtt <- log(0.773)
    label("Mean transit time of the 1-transit-compartment absorption chain MTT (h)")
    # Hoglund 2016 Table 1: MTT = 0.773 h (RSE 43.1%).
    # NOTE: the bootstrap 95% CI printed alongside is 0.809-2.38, which
    # does NOT contain the point estimate. The final-model estimate is
    # used here because ini() encodes the final model, and the bootstrap
    # converged for only 603 of 1000 runs (Table 1 footnote a) on a
    # parameter the authors separately found poorly determined (adding
    # BSV on MTT gave a 149% RSE and was dropped; Results). Unlike the
    # disposition parameters, MTT has no independent falsifier in the
    # paper -- the terminal half-lives do not constrain absorption -- so
    # the discrepancy is recorded rather than resolved. See vignette
    # Errata.
    # Number of transit compartments = 1 (Results: "a one transit
    # compartment model for the absorption of chloroquine"), so the
    # chain depot -> transit1 -> central has NN + 1 = 2 equal-rate
    # transfers and ktr = (NN + 1) / MTT = 2 / MTT (Savic & Karlsson
    # 2007 convention; same idiom as the sibling Hoglund_2018_mefloquine.R
    # with NN = 1 and Ali_2018_amodiaquine.R with NN = 2).

    lcl <- log(6.13)
    label("Apparent total chloroquine elimination clearance CL CQ/F (L/h)")
    # Hoglund 2016 Table 1: CL CQ/F = 6.13 L/h (RSE 3.40%; bootstrap
    # 95% CI 5.74-6.55). This is TOTAL apparent clearance out of the
    # chloroquine central compartment -- confirmed by the terminal
    # half-life check above (k10 = 6.13 / 468 reproduces the reported
    # 10.7 days) and by Fig. 1, which draws k10 (elimination) and k23
    # (transformation to DCQ) as two SEPARATE arrows leaving Central CQ.
    # The Table 1 footnote glosses this symbol as "apparent clearance of
    # CQ for transformation into desethylchloroquine", which conflicts
    # with both; the footnote is a drafting leftover, and fm below
    # carries the transformation fraction.

    lvc <- log(468)
    label("Apparent chloroquine central volume of distribution Vc CQ/F (L)")
    # Hoglund 2016 Table 1: VC CQ/F = 468 L (RSE 16.0%; bootstrap
    # 95% CI 137-529)

    lq <- log(37.7)
    label("Apparent chloroquine inter-compartmental clearance Q CQ/F (L/h)")
    # Hoglund 2016 Table 1: QCQ/F = 37.7 L/h (RSE 18.9%; bootstrap
    # 95% CI 31.1-69.0)

    lvp <- log(1600)
    label("Apparent chloroquine peripheral volume of distribution Vp CQ/F (L)")
    # Hoglund 2016 Table 1: VP CQ/F = 1600 L (RSE 5.21%; bootstrap
    # 95% CI 1470-1800)

    # ---------------------------------------------------------------
    # Desethylchloroquine disposition. Two of the four values printed in
    # Hoglund 2016 Table 1 are internally impossible, and both are
    # repaired here from the paper's OWN data. The repairs are not
    # substitutions from outside the source.
    #
    #   printed VP DCQ/F = 566,257 L with bootstrap 95% CI 198-341
    #   printed QDCQ/F   =   31.46 L/h with bootstrap 95% CI 1.11-1.83
    #
    # Each printed point estimate falls outside its own bootstrap CI --
    # VP DCQ by more than three orders of magnitude. Substituting the
    # four printed values into the two-compartment terminal-slope
    # formula gives a desethylchloroquine half-life of about 8500 days
    # against the 8.74 days Table 1 itself reports.
    #
    # Taking VP DCQ/F = 257 L and QDCQ/F = 1.46 L/h instead reproduces
    # the reported half-life to within rounding: beta = 0.0033061 /h,
    # t1/2 = 8.736 days against the reported 8.74 days. Both repaired
    # values lie inside their own bootstrap CIs (257 in 198-341; 1.46 in
    # 1.11-1.83, essentially at its centre), and the digit strings "257"
    # and "1.46" are literally present in the mangled cells -- the
    # failure is a typesetting mash ("566" + "257"; a stray leading "3"),
    # not a different set of numbers. Three independent constraints
    # (bootstrap CI, reported half-life, printed digits) agree on both.
    # The same half-life formula reproduces the chloroquine half-life
    # exactly from the unmodified chloroquine row set, which is what
    # establishes that this IS the authors' own calculation.
    # Recorded prominently in the vignette Errata.
    # ---------------------------------------------------------------

    lcl_dcq <- log(2.04)
    label("Apparent desethylchloroquine elimination clearance CL DCQ/F (L/h)")
    # Hoglund 2016 Table 1: CL DCQ/F = 2.04 L/h (RSE 3.50%; bootstrap
    # 95% CI 1.90-2.18). Printed value used as-is.

    lvc_dcq <- log(2.27)
    label("Apparent desethylchloroquine central volume of distribution Vc DCQ/F (L)")
    # Hoglund 2016 Table 1: VC DCQ/F = 2.27 L (RSE 14.1%; bootstrap
    # 95% CI 1.62-2.90). Printed value used as-is -- it is consistent
    # with its own bootstrap CI and is required by the half-life check
    # above. The value is small for a central volume; the authors
    # address this directly, attributing the low desethylchloroquine
    # volume estimates to the fixation of the formed fraction at 18%
    # (Discussion).

    lq_dcq <- log(1.46)
    label("Apparent desethylchloroquine inter-compartmental clearance Q DCQ/F (L/h)")
    # Hoglund 2016 Table 1 prints "31.46 (12.3)" with bootstrap 95% CI
    # 1.11-1.83. Repaired to 1.46 L/h: see the block comment above.

    lvp_dcq <- log(257)
    label("Apparent desethylchloroquine peripheral volume of distribution Vp DCQ/F (L)")
    # Hoglund 2016 Table 1 prints "566,257 (14.4)" with bootstrap 95% CI
    # 198-341. Repaired to 257 L: see the block comment above.

    fm <- fixed(0.18)
    label("Fraction of chloroquine clearance forming desethylchloroquine (unitless)")
    # Hoglund 2016 Results: "The parameter describing the transformation
    # of chloroquine into desethylchloroquine (CLm) was fixed to 18% of
    # the transformation clearance from parent drug to metabolite [26]".
    # The Discussion states what the 18% is: "the fixation of CLm to
    # 18%, an estimation based on the fraction desethylchloroquine of
    # the total chloroquine dose recovered in urine" -- i.e. a fraction
    # metabolised taken from external data, not estimated here. This is
    # the FIXED-fm resolution of the fm / metabolite-volume
    # identifiability problem described in the `fm` register entry, so
    # central_dcq holds a true amount and the desethylchloroquine
    # volumes above are true (not fm-scaled) apparent volumes.

    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability F (unitless, reference value 1)")
    # Hoglund 2016 Methods, Population pharmacokinetics: "Relative
    # bioavailability was added with a typical value of 100% with an
    # estimate of between-subject variability." The typical value is
    # therefore fixed at 1 and only the BSV (etalfdepot below) is
    # estimated; Table 1 lists a BSV F row and no F point-estimate row.

    # ---------------------------------------------------------------
    # Between-subject variability. Hoglund 2016 Table 1 "Interindividual
    # variability" block, footnote b: "Listed as coefficient of
    # variation (CV; %) and their RSE (%) in parenthesis". BSV is
    # exponential (Methods, Eq: Pi = Pp * exp(eta_i)), so the parameters
    # are log-normal and the internal log-scale variance is recovered as
    # omega^2 = log(CV^2 + 1) -- the same conversion used by the sibling
    # Hoglund_2018_mefloquine.R.
    #
    #   BSV VC DCQ  48.7% -> log(0.487^2 + 1) = 0.212826
    #   BSV VP CQ   20.0% -> log(0.200^2 + 1) = 0.039221
    #   BSV VP DCQ  86.8% -> log(0.868^2 + 1) = 0.561570
    #   BSV F       19.4% -> log(0.194^2 + 1) = 0.036945
    #
    # BSV was retained on exactly these four parameters (Results: "BSV
    # was kept on the apparent volume of distribution of
    # desethylchloroquine, relative bioavailability and the apparent
    # volume of distribution in the peripheral compartment of
    # chloroquine and desethylchloroquine"). BSV on MTT was tested and
    # rejected (149% RSE), and no BSV is reported on any clearance.
    # No correlation structure is reported, so the etas are independent.
    # ---------------------------------------------------------------

    etalvc_dcq ~ 0.212826
    # Hoglund 2016 Table 1: BSV VC DCQ = 48.7% CV (RSE 47.5%; bootstrap 95% CI 17.9-71.1)

    etalvp ~ 0.039221
    # Hoglund 2016 Table 1: BSV VP CQ = 20.0% CV (RSE 61.6%; bootstrap 95% CI 8.25-31.9)

    etalvp_dcq ~ 0.561570
    # Hoglund 2016 Table 1: BSV VP DCQ = 86.8% CV (RSE 30.5%; bootstrap 95% CI 49.1-116)

    etalfdepot ~ 0.036945
    # Hoglund 2016 Table 1: BSV F = 19.4% CV (RSE 31.7%; bootstrap 95% CI 13.1-25.4)

    # ---------------------------------------------------------------
    # Residual error. Hoglund 2016 Methods: concentrations were
    # "transformed into their natural logarithms" and an additive
    # residual model Cobs = Cp + eps_ad was applied, with the paper
    # noting "An additive model on log-transformed data is equivalent to
    # an exponential model". An additive-on-log-scale residual maps to a
    # proportional residual in nlmixr2's linear concentration space (the
    # same mapping used by Hoglund_2012_piperaquine.R,
    # Hoglund_2017_piperaquine.R and Hoglund_2018_mefloquine.R).
    #
    # Table 1 reports these as SDs on the log scale, not variances: the
    # printed RSE of 5.34% on 0.401 implies a 95% interval of
    # 0.401 * (1 +/- 1.96 * 0.0534) = 0.359-0.443, which matches the
    # printed bootstrap CI of 0.360-0.444 exactly. So 0.401 is the
    # parameter itself and is used directly as the SD.
    # ---------------------------------------------------------------

    propSd <- 0.401
    label("Proportional residual SD for chloroquine whole-blood concentration (SD on log scale)")
    # Hoglund 2016 Table 1: Proportional error CQ = 0.401 (RSE 5.34%;
    # bootstrap 95% CI 0.360-0.444)

    propSd_dcq <- 0.431
    label("Proportional residual SD for desethylchloroquine whole-blood concentration (SD on log scale)")
    # Hoglund 2016 Table 1: 'Proporional error DCQ' = 0.431 (RSE 4.97%;
    # bootstrap 95% CI 0.393-0.479). Row label misspelled in the source.
  })

  model({
    # Molecular weights of the free base of each analyte (g/mol).
    # Chloroquine is C18H26ClN3 and desethylchloroquine is C16H22ClN3
    # (chloroquine less one ethyl group, C2H4, 28.05 g/mol). The
    # transformation is 1:1 in moles, and the authors fitted on the
    # molar scale (the observed-concentration axes of Figs. 2, 3 and 4
    # are all umole/L), so the mass flux leaving the chloroquine central
    # compartment is multiplied by mwDCQ / mwCQ = 0.91231 to become the
    # mass flux of desethylchloroquine entering central_dcq. This is the
    # same molar-correction idiom used by the sibling antimalarial
    # parent + desethyl-metabolite models Ali_2018_amodiaquine.R and
    # Ding_2024_amodiaquine.R, and is what lets both analytes be dosed
    # and read out in mass units here. The molecular weights are
    # standard chemical constants, not fitted values.
    mwCQ        <- 319.87
    mwDCQ       <- 291.82
    molarFactor <- mwDCQ / mwCQ

    # Absorption chain rate constant. NN = 1 transit compartment, so the
    # chain depot -> transit1 -> central has NN + 1 = 2 equal-rate
    # transfers and ktr = (NN + 1) / MTT = 2 / MTT.
    mtt <- exp(lmtt)
    ktr <- 2 / mtt

    # Individual PK parameters. The final model is covariate-free. BSV
    # is carried only on Vp CQ, Vc DCQ, Vp DCQ and F (Results); every
    # other parameter is a typical value.
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    cl_dcq <- exp(lcl_dcq)
    vc_dcq <- exp(lvc_dcq + etalvc_dcq)
    q_dcq  <- exp(lq_dcq)
    vp_dcq <- exp(lvp_dcq + etalvp_dcq)

    # Two-compartment micro-rate constants for each analyte (1/h).
    kel     <- cl     / vc
    k12     <- q      / vc
    k21     <- q      / vp
    kel_dcq <- cl_dcq / vc_dcq
    k12_dcq <- q_dcq  / vc_dcq
    k21_dcq <- q_dcq  / vp_dcq

    # ODE system, replicating Hoglund 2016 Fig. 1. Compartment amounts
    # are in mg of the analyte base and volumes are in L, so amount /
    # volume is mg/L and is scaled by 1000 below to give ng/mL.

    # Absorption: depot -> transit1 -> chloroquine central.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1

    # Chloroquine central + one peripheral. Total clearance out of
    # central is kel * central; Fig. 1 splits this into k10 (other
    # elimination) and k23 (transformation to DCQ), and fm is the
    # fraction taken by k23. The fm split does NOT change the total
    # rate of loss from central, which is what the reported chloroquine
    # terminal half-life of 10.7 days constrains.
    d/dt(central)     <-  ktr * transit1 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Desethylchloroquine central + one peripheral, fed by the fm
    # fraction of chloroquine elimination after the molar correction.
    d/dt(central_dcq)     <-  molarFactor * fm * kel * central -
                              kel_dcq * central_dcq -
                              k12_dcq * central_dcq + k21_dcq * peripheral1_dcq
    d/dt(peripheral1_dcq) <-  k12_dcq * central_dcq - k21_dcq * peripheral1_dcq

    # Relative oral bioavailability on the depot compartment: typical
    # value fixed at 1 with between-subject variability only.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Whole-blood concentrations in ng/mL: amount (mg) / volume (L) is
    # mg/L, which is 1000 ng/mL.
    Cc     <- 1000 * central     / vc
    Cc_dcq <- 1000 * central_dcq / vc_dcq

    # Additive residual on the log scale in the source maps to a
    # proportional residual here (see the ini() comment on propSd).
    Cc     ~ prop(propSd)
    Cc_dcq ~ prop(propSd_dcq)
  })
}
