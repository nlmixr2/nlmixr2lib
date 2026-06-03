Syvanen_2011_verapamil_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Sprague-Dawley). Population mixed-effects",
    "popPK model for (R)-[11C]verapamil in plasma and whole-brain PET",
    "tissue, fit by Syvanen et al. (2011, BMC Med Imaging) as part of a",
    "PET study comparing P-glycoprotein (P-gp) functionality at the",
    "blood-brain barrier between kainate-induced post-status-epilepticus",
    "rats (n = 22) and saline-treated controls (n = 20), with paired",
    "tariquidar (15 mg/kg IV) vs vehicle co-administration arms. The",
    "structural model is a three-compartment plasma disposition",
    "(central + 2 peripherals) coupled to a two-compartment brain model",
    "(brain_csf = fast-exchange brain compartment connected to plasma",
    "via Qin in / Qout out; brain_deep = deep-brain compartment exchanging",
    "with brain_csf via Qbr). Plasma curves are complete-metabolite-",
    "corrected before fitting, so the model describes intact",
    "(R)-[11C]verapamil kinetics only. Body weight is the only continuous",
    "covariate (allometric on plasma CL, reference weight 0.3084 kg).",
    "Tariquidar co-administration multiplies Vp1 by 1.20, Vbr1 by 2.41,",
    "and Qin by 12.0; the kainate-induced post-SE state multiplies Vbr1",
    "by 1.32 (no significant effect on Qin or Qout); both categorical",
    "effects use the paper's theta^COV multiplicative form (Equation 5)."
  )
  reference <- paste(
    "Syvanen S, Luurtsema G, Molthoff CFM, Windhorst AD, Huisman MC,",
    "Lammertsma AA, Voskuyl RA, de Lange ECM.",
    "(R)-[11C]verapamil PET studies to assess changes in P-glycoprotein",
    "expression and functionality in rat blood-brain barrier after",
    "exposure to kainate-induced status epilepticus.",
    "BMC Med Imaging. 2011;11:1. doi:10.1186/1471-2342-11-1.",
    sep = " "
  )
  vignette <- "Syvanen_2011_verapamil_rat"

  paper_specific_etas <- c("etalvbr")
  paper_specific_residual_sds <- c("propSd_Cbrain")

  units <- list(
    time          = "min",
    dosing        = "MBq",
    concentration = "MBq/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at the time of PET scanning.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per rat in the source study (one body-weight value per",
        "animal recorded on the experimental day, Table 1). Used for",
        "allometric scaling of plasma clearance with reference weight",
        "0.3084 kg (the mean body weight across all 42 successfully scanned",
        "animals; paper Population mixed effects modelling section and",
        "Table 3 footer). The four rat-group mean weights are 327, 292,",
        "326, and 293 g; the overall mean is 308 g, rounded to 0.3084 kg",
        "as the model's reference. Allometric form: CL * (WT / 0.3084)^1.98",
        "(paper Equation 4)."
      ),
      source_name        = "WT"
    ),
    CONMED_TARIQUIDAR = list(
      description        = paste(
        "Indicator for tariquidar (P-glycoprotein inhibitor) pre-",
        "co-administration: 1 = animal received a single 15 mg/kg IV bolus",
        "of tariquidar 20-30 minutes prior to the (R)-[11C]verapamil",
        "injection; 0 = animal received only vehicle (3 mL/kg of 5%",
        "glucose in saline)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (vehicle co-administration; no tariquidar)",
      notes              = paste(
        "Time-fixed per animal in the source study (each rat is allocated",
        "to either the tariquidar or vehicle arm; 21 of the 42 successfully",
        "scanned rats received tariquidar). Tariquidar plasma concentrations",
        "2 h post-dose were 1.8 +/- 0.3 ug/mL (saline) and 1.6 +/- 0.3",
        "ug/mL (kainate), considered equivalent across arms so the indicator",
        "is treated as binary rather than continuous (paper Results section",
        "'P-gp blockage by Tariquidar co-administration'). Multiplicative",
        "fractional-change form 'theta^COVARIATE' applied to Vp1 (theta =",
        "1.20), Vbr1 (theta = 2.41), and Qin (theta = 12.0) -- the paper",
        "screened it as a covariate on Qout but did not retain it; only",
        "the three retained effects are encoded here (paper Table 3 and",
        "Equation 5)."
      ),
      source_name        = "tariquidar (categorical 0/1)"
    ),
    DIS_POSTSE_KAINATE = list(
      description        = paste(
        "Indicator for the post-status-epilepticus state induced by kainic",
        "acid 7 days prior to PET scanning: 1 = animal received repeated",
        "IP injections of kainic acid (10 mg/kg followed by 5 mg/kg every",
        "30-60 min until stage IV-V seizures by Racine scale, or 30 mg/kg",
        "total) 7 days before scanning; 0 = animal received an equivalent",
        "volume of saline (control)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (saline control; no kainate-induced status epilepticus)",
      notes              = paste(
        "Time-fixed per animal in the source study (each rat is allocated",
        "to either the kainate or saline arm; 22 of the 42 successfully",
        "scanned rats were kainate-treated). The 7-day post-induction time",
        "point was chosen because previous studies indicated P-gp",
        "expression peaks around this time; at scanning the rats had not",
        "yet developed spontaneous epilepsy (typically begins 3-4 weeks",
        "after SE) but were progressing toward the chronic epileptic",
        "state. Multiplicative fractional-change form 'theta^COVARIATE'",
        "applied only to Vbr1 (theta = 1.32); the paper screened it as a",
        "covariate on Qin and Qout but did not retain it (paper Table 3",
        "and Equation 5, plus Discussion p1789 explaining that the Vbr1",
        "increase reflects increased intracellular brain distribution",
        "rather than altered BBB transport)."
      ),
      source_name        = "rat group (kainate vs saline; categorical 0/1)"
    )
  )

  population <- list(
    species        = "rat (male Sprague-Dawley, Harlan, Horst, The Netherlands)",
    n_subjects     = 42L,
    n_studies      = 1L,
    age_range      = "Adult; specific age not reported (animals housed for ~1 week of habituation after arrival, then 7 days of post-SE / post-saline observation before scanning)",
    weight_range   = "200-224 g at arrival; 308 +/- 31 g (mean +/- SD) on the experimental day across all 42 successfully scanned rats (Table 1)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Four treatment arms: saline + vehicle (n = 10), kainate + vehicle",
      "(n = 11), saline + tariquidar (n = 10), kainate + tariquidar",
      "(n = 11). Kainate-treated rats had received repeated IP kainic acid",
      "injections 7 days prior to scanning to induce status epilepticus",
      "(mimicking the chronic-epilepsy onset phase); the kainate arm",
      "represents an early post-SE pre-epileptic state. P-gp",
      "immunohistochemistry showed only a small non-significant trend",
      "toward higher P-gp expression in kainate-treated rats vs controls."
    ),
    dose_range     = paste(
      "Single intravenous bolus injection of no-carrier-added",
      "(R)-[11C]verapamil, 20.8 +/- 4.5 MBq (0.15-0.3 mL) at PET scan",
      "start; specific activity 56 +/- 26 GBq/umol (Table 1).",
      "Tariquidar pre-dose (when applicable): 15 mg/kg IV bolus in 3",
      "mL/kg vehicle (5% glucose in saline), administered 20-30 min",
      "before the (R)-[11C]verapamil injection."
    ),
    regions        = "preclinical (in-vivo rat); Leiden University, The Netherlands",
    notes          = paste(
      "Of an original 50 animals, 8 were excluded (4 died from kainate",
      "treatment, 1 did not respond to kainate, 3 had technical issues",
      "such as scanner movement or extravascular tracer injection),",
      "leaving 42 successfully scanned rats. Anaesthesia: isoflurane",
      "(4% induction / 2% maintenance) in 1 L/min oxygen via nose mask.",
      "Femoral artery and vein cannulated 1-2 h pre-scan for tracer",
      "administration and arterial blood sampling. PET acquisition on a",
      "HRRT (CTI / Siemens) scanner for 60 min in pairs of rats; frame",
      "sequence 6 x 10, 2 x 30, 3 x 60, 2 x 150, 2 x 300, 4 x 600",
      "seconds. Arterial blood sampling at 0.5, 1, 3, 5, 10, 15, 20,",
      "30, 45, 60 min for total plasma activity; subset metabolite",
      "analysis at 10, 30, 60 min. Plasma curves used for fitting are",
      "complete-metabolite-corrected (intact-fraction multiplied",
      "throughout), so the model describes intact",
      "(R)-[11C]verapamil only. Brain time-activity curves derived from",
      "a 127 uL whole-brain VOI drawn on the co-acquired [18F]FDG image,",
      "corrected for an intravascular contribution assuming 3% cerebral",
      "blood volume. Model fit using NONMEM VI with ADVAN9 + FOCE-I."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Final NONMEM estimates from Table 3 of
    # Syvanen et al. 2011. Reference units retained from the paper:
    # volumes in mL, clearances in mL/min, time in min, dose in MBq.
    # Plasma compartments use canonical names (central / peripheral1 /
    # peripheral2); brain compartments use brain_csf (Vbr1, fast-exchange
    # brain compartment connected to plasma) and brain_deep (Vbr2, the
    # deep brain compartment), following the Xie 2000 m3g_rat precedent.
    # ------------------------------------------------------------------

    # Plasma central + 2 peripherals (paper Table 3, "Plasma" block)
    lvc  <- log(22.8);  label("Plasma central volume Vc (mL)")                              # Syvanen 2011 Table 3: Vc = 22.8 mL (SE 2.41)
    lvp  <- log(504);   label("Plasma peripheral 1 volume Vp1 (mL)")                        # Syvanen 2011 Table 3: Vp1 = 504 mL (SE 40.7)
    lvp2 <- log(70.2);  label("Plasma peripheral 2 volume Vp2 (mL)")                        # Syvanen 2011 Table 3: Vp2 = 70.2 mL (SE 9.91)
    lcl  <- log(14.7);  label("Plasma systemic clearance CL (mL/min)")                      # Syvanen 2011 Table 3: CL = 14.7 mL/min (SE 0.681)
    lq   <- log(16.1);  label("Inter-compartmental clearance central <-> peripheral1 Q1 (mL/min)") # Syvanen 2011 Table 3: Q1 = 16.1 mL/min (SE 1.75)
    lq2  <- log(22.7);  label("Inter-compartmental clearance central <-> peripheral2 Q2 (mL/min)") # Syvanen 2011 Table 3: Q2 = 22.7 mL/min (SE 3.34)

    # Brain two-compartment model (paper Table 3, "Brain" block).
    # Paper-specific PK-parameter names (lvbr1, lvbr2, lqin, lqout, lqbr)
    # because the brain-side parameterisation does not collapse to the
    # canonical PK V / Q register.
    lvbr1 <- log(7.23);   label("Brain compartment 1 (brain_csf) volume Vbr1 (mL)")         # Syvanen 2011 Table 3: Vbr1 = 7.23 mL (SE 2.00)
    lvbr2 <- log(10.7);   label("Brain compartment 2 (brain_deep) volume Vbr2 (mL)")        # Syvanen 2011 Table 3: Vbr2 = 10.7 mL (SE 2.22)
    lqin  <- log(1.75);   label("BBB influx clearance central -> brain_csf Qin (mL/min)")   # Syvanen 2011 Table 3: Qin = 1.75 mL/min (SE 0.26)
    lqout <- log(1.81);   label("BBB efflux clearance brain_csf -> central Qout (mL/min)")  # Syvanen 2011 Table 3: Qout = 1.81 mL/min (SE 0.283)
    lqbr  <- log(0.692);  label("Inter-compartmental brain clearance brain_csf <-> brain_deep Qbr (mL/min)") # Syvanen 2011 Table 3: Qbr = 0.692 mL/min (SE 0.143)

    # ------------------------------------------------------------------
    # Covariate effects. The paper uses two functional forms (Equation 5
    # and Table 3 footer):
    #   - Categorical (tariquidar, kainate post-SE): theta^COVARIATE,
    #     applied as a multiplicative ratio when COVARIATE = 1. Stored
    #     here as log(theta) so the model() block can add the effect on
    #     the log scale (cleaner than power form for numerical stability
    #     and consistent with the Xie 2000 precedent).
    #   - Continuous (body weight): (WT/0.3084)^theta, applied as a
    #     power-law on the body-weight ratio. Stored here as theta
    #     (the allometric exponent) directly.
    # All covariate-effect THETAs were estimated, not fixed, so no
    # fixed() wrappers (paper Discussion p1789 -- stepwise forward
    # addition / backward deletion identified these as the retained
    # significant covariates).
    # ------------------------------------------------------------------
    e_conmed_tariquidar_vp   <- log(1.20)
    label("Log fractional change in Vp1 with tariquidar (paper theta = 1.20, ratio form)")  # Syvanen 2011 Table 3: tariquidar covariate on Vp1 = 1.20 (SE 0.119)
    e_conmed_tariquidar_vbr1 <- log(2.41)
    label("Log fractional change in Vbr1 with tariquidar (paper theta = 2.41, ratio form)") # Syvanen 2011 Table 3: tariquidar covariate on Vbr1 = 2.41 (SE 0.505)
    e_conmed_tariquidar_qin  <- log(12.0)
    label("Log fractional change in Qin with tariquidar (paper theta = 12.0, ratio form)")  # Syvanen 2011 Table 3: tariquidar covariate on Qin = 12.0 (SE 0.554)
    e_dis_postse_kainate_vbr1 <- log(1.32)
    label("Log fractional change in Vbr1 with kainate post-SE state (paper theta = 1.32, ratio form)") # Syvanen 2011 Table 3: kainate covariate on Vbr1 = 1.32 (SE 0.416)
    e_wt_cl <- 1.98
    label("Allometric body-weight exponent on plasma CL (unitless, paper Equation 4)")       # Syvanen 2011 Table 3: weight covariate on CL = 1.98 (SE 0.362), reference weight 0.3084 kg

    # ------------------------------------------------------------------
    # Inter-animal variability (paper "Inter-individual variability"
    # column, Table 3). Paper Equation 3:
    #   theta_i = theta_pop * exp(eta_i),  eta_i ~ N(0, omega)
    # The values reported in Table 3 are NONMEM-convention variances
    # (omega^2 on the log scale); SDs would be sqrt(value). Paper text:
    # "Inter-individual variation was investigated for all parameters,
    # but incorporated only for those parameters for which it
    # significantly (Delta-OFV > 3.83) improved the model."
    #
    # Shared etas: Q1 and Q2 are reported with identical variance estimates
    # (0.220) and identical SEs (0.073); Vbr1 and Vbr2 are reported with
    # identical variance estimates (0.115) and identical SEs (0.054). The
    # identical SE-to-the-same-decimal-place is essentially impossible if
    # two independent etas were estimated separately, so the encoding here
    # interprets these as shared etas in the NONMEM idiom
    # (Q1 = TVQ1 * exp(ETA_Q); Q2 = TVQ2 * exp(ETA_Q); same for Vbr).
    # This shared-eta interpretation is documented in the vignette's
    # Assumptions and deviations section.
    # ------------------------------------------------------------------
    etalvc  ~ 0.175  # Syvanen 2011 Table 3: IIV(Vc)  variance = 0.175 (SE 0.085)
    etalvp2 ~ 0.036  # Syvanen 2011 Table 3: IIV(Vp2) variance = 0.036 (SE 0.034)
    etalcl  ~ 0.065  # Syvanen 2011 Table 3: IIV(CL)  variance = 0.065 (SE 0.016)
    etalq   ~ 0.220  # Syvanen 2011 Table 3: IIV(Q1) == IIV(Q2) = 0.220 (SE 0.073); encoded as a shared eta applied to both q1 and q2 in model()
    etalvbr ~ 0.115  # Syvanen 2011 Table 3: IIV(Vbr1) == IIV(Vbr2) = 0.115 (SE 0.054); encoded as a shared eta applied to both vbr1 and vbr2 in model()

    # ------------------------------------------------------------------
    # Residual error. Paper "Residual errors" block of Table 3, plus the
    # Methods text: "Proportional error models were included for the
    # residual variability." The values 0.118 (blood) and 0.226 (brain)
    # are reported with their SEs (0.0152 and 0.0354) and are interpreted
    # as proportional SDs in the model Y_obs = Y_pred * (1 + eps),
    # eps ~ N(0, sigma^2). Mapped directly to propSd (plasma Cc) and
    # propSd_Cbrain (whole-brain VOI observation).
    # ------------------------------------------------------------------
    propSd        <- 0.118
    label("Proportional residual SD on plasma (R)-[11C]verapamil concentration Cc (fraction)") # Syvanen 2011 Table 3: residual error 'blood' = 0.118 (SE 0.0152)
    propSd_Cbrain <- 0.226
    label("Proportional residual SD on whole-brain (R)-[11C]verapamil concentration Cbrain (fraction)") # Syvanen 2011 Table 3: residual error 'brain' = 0.226 (SE 0.0354)
  })

  model({
    # ------------------------------------------------------------------
    # Derived individual parameters. The categorical covariate effects
    # use the paper's theta^COVARIATE multiplicative form (Equation 5),
    # which on the log scale is log(theta) * COVARIATE -- equivalent to
    # the Xie 2000 e_<cov>_<param> exponential-effect idiom. The weight
    # covariate uses (WT / 0.3084)^e_wt_cl directly per Equation 4.
    # ------------------------------------------------------------------
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + e_conmed_tariquidar_vp * CONMED_TARIQUIDAR)
    vp2 <- exp(lvp2 + etalvp2)
    cl  <- exp(lcl  + etalcl) * (WT / 0.3084)^e_wt_cl
    q   <- exp(lq   + etalq)
    q2  <- exp(lq2  + etalq)  # shared eta with q (Q1 and Q2 share IIV per Table 3)

    vbr1 <- exp(lvbr1 + e_conmed_tariquidar_vbr1 * CONMED_TARIQUIDAR +
                e_dis_postse_kainate_vbr1 * DIS_POSTSE_KAINATE + etalvbr)
    vbr2 <- exp(lvbr2 + etalvbr)  # shared eta with vbr1 (Vbr1 and Vbr2 share IIV per Table 3)
    qin  <- exp(lqin  + e_conmed_tariquidar_qin * CONMED_TARIQUIDAR)
    qout <- exp(lqout)
    qbr  <- exp(lqbr)

    # ------------------------------------------------------------------
    # Concentrations driving the ODE flux terms. Paper-faithful form:
    # transfer between compartments is rate = Q * (concentration). The
    # paper's K1 macroconstant (PET notation) maps to Qin / (Vbr1 + Vbr2)
    # via the table footer relation 'K1 = Qin/Vc * (Vc/(Vbr1+Vbr2))'.
    # ------------------------------------------------------------------
    cpl    <- central     / vc
    cbr1   <- brain_csf   / vbr1
    cbr2   <- brain_deep  / vbr2

    # ODE system. Plasma central exchanges with two peripherals (Q1, Q2)
    # and clears systemically (CL); brain_csf exchanges with plasma via
    # asymmetric Qin / Qout and with brain_deep via symmetric Qbr.
    # Direction conventions: central -> brain_csf at Qin * cpl;
    # brain_csf -> central at Qout * cbr1; brain_csf <-> brain_deep at
    # Qbr * (cbr1 - cbr2).
    d/dt(central)     <- -cl  * cpl -
                          q   * (cpl - peripheral1 / vp) -
                          q2  * (cpl - peripheral2 / vp2) -
                          qin * cpl +
                          qout * cbr1
    d/dt(peripheral1) <-  q   * (cpl - peripheral1 / vp)
    d/dt(peripheral2) <-  q2  * (cpl - peripheral2 / vp2)
    d/dt(brain_csf)   <-  qin  * cpl -
                          qout * cbr1 -
                          qbr  * (cbr1 - cbr2)
    d/dt(brain_deep)  <-  qbr  * (cbr1 - cbr2)

    # ------------------------------------------------------------------
    # Observations. Cc is the plasma activity concentration (paper-named
    # 'blood' in Table 3; the plasma curve is complete-metabolite-corrected
    # before model fitting). Cbrain is the whole-brain VOI activity
    # concentration as observed by PET, computed as total activity in
    # both brain compartments per total brain volume (Vbr1 + Vbr2).
    # ------------------------------------------------------------------
    Cc     <- cpl
    Cbrain <- (brain_csf + brain_deep) / (vbr1 + vbr2)

    Cc     ~ prop(propSd)
    Cbrain ~ prop(propSd_Cbrain)
  })
}
