Hashimoto_1994_zonisamide <- function() {
  description <- "Steady-state Michaelis-Menten population PK model for zonisamide in 68 Japanese epileptic patients (pediatric + adult) on chronic oral zonisamide. A power-of-weight body-size factor scales both volume of distribution and Vmax; concomitant carbamazepine multiplicatively increases Vmax (Hashimoto 1994 Eqs. 1-4)."
  reference   <- "Hashimoto Y, Odani A, Tanigawara Y, Yasuhara M, Okuno T, Hori R. Population analysis of the dose-dependent pharmacokinetics of zonisamide in epileptic patients. Biol Pharm Bull. 1994;17(2):323-326. doi:10.1248/bpb.17.323"
  vignette    <- "Hashimoto_1994_zonisamide"
  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 33 kg per Hashimoto 1994 Eq. 2 (page 324): SIZE_i = 33 * (WT_i/33)^theta1 with theta1 = 0.741 (Table II). The same SIZE factor scales both V (Eq. 3) and Vmax (Eq. 4 rearrangement). Time-fixed at baseline in the TDM dataset (chronic oral maintenance therapy).",
      source_name        = "WT"
    ),
    CONMED_CBZ = list(
      description        = "Indicator for concomitant carbamazepine coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = "Hashimoto 1994 Eq. 4 (page 324) defines CBZ = 1 in patients receiving carbamazepine concomitantly, otherwise CBZ = 0. When CONMED_CBZ = 1, Vmax is multiplied by theta2 = 1.13 (Table II), i.e., Vmax is 13 percent higher in patients on carbamazepine. Valproate and phenytoin coadministration were tested as covariates but did not significantly affect zonisamide PK (page 325) and are not included. Phenobarbital was not tested because only 5 of 68 patients received it.",
      source_name        = "CBZ"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in Hashimoto 1994 preliminary analysis (page 325, Figs. 1-2 stratified younger <=10 yr vs older >10 yr) but not retained as a covariate in the final model. The age effect on zonisamide PK is fully captured by the power-of-weight SIZE factor (Eq. 2); the age-stratified fits of Eq. 6 were identical across the two age groups."
    ),
    CONMED_VPA = list(
      description = "Indicator for concomitant valproate coadministration",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in Hashimoto 1994 (page 325, 'data not shown') but not retained as a covariate; valproate did not significantly affect the size-corrected dose vs serum concentration relationship."
    ),
    CONMED_PHT = list(
      description = "Indicator for concomitant phenytoin coadministration",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in Hashimoto 1994 (page 325, 'data not shown') but not retained as a covariate; phenytoin did not significantly affect the size-corrected dose vs serum concentration relationship."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 68L,
    n_observations   = 266L,
    n_studies        = 1L,
    age_range        = "30 patients <=10 yr, 5 patients >20 yr (remainder 11-20 yr)",
    age_median       = "mean 11.2 (SD 6.4) years",
    weight_range     = "approx 10-70 kg (mean 33.4 SD 17.7)",
    weight_median    = "mean 33.4 (SD 17.7) kg",
    sex_female_pct   = 45.6,
    race_ethnicity   = "Japanese (single-centre cohort at Kyoto University Hospital)",
    disease_state    = "Epileptic patients on chronic oral zonisamide (Excegran tablet or powder, Dainippon Pharmaceutical Co., Osaka). 2 patients on zonisamide alone; 37 also on carbamazepine, 32 on valproate, 27 on phenytoin, 5 on phenobarbital (overlap permitted).",
    dose_range       = "Mean daily dose 135 (SD 104) mg/d, administered orally at 12-hour intervals to 62 of 68 patients (tau = 0.5 day). 26 patients sampled at a single dose level; 17 at three or more dose levels.",
    regions          = "Japan (Kyoto University Hospital, single centre, November 1989-July 1992).",
    notes            = "Hashimoto 1994 Table I baseline demographics. 60 outpatients, 31 females, 37 males. 78 percent (208/266) of samples taken 2-6 h post-dose (approximate peak); 12 percent (33/266) at trough 12 h post-dose. All samples obtained at steady-state after more than one month of stable therapy. Bioavailability assumed F = 1 (page 326)."
  )

  ini({
    # Structural Michaelis-Menten parameters - Hashimoto 1994 Table II (page 325).
    # Reference covariates: 33-kg patient without concomitant carbamazepine.
    # V is reported as L/kg (Table II): typical Vc at 33 kg = 1.27 L/kg * 33 kg
    # = 41.91 L (matches the page-326 statement "V was estimated to be 1.27 l/kg
    # in a typical 33-kg patient"). Vmax is reported as mg/d/kg: typical Vmax
    # at 33 kg = 27.6 mg/d/kg * 33 kg = 910.8 mg/d.
    lvc   <- log(41.91); label("Volume of distribution Vc at reference 33 kg (L)")  # Table II: V = 1.27 L/kg; 1.27 * 33 = 41.91 L
    lvmax <- log(910.8); label("Maximum elimination rate Vmax at reference 33 kg (mg/day)")  # Table II: Vmax = 27.6 mg/d/kg; 27.6 * 33 = 910.8 mg/d
    lkm   <- log(45.9);  label("Michaelis-Menten constant Km (mg/L)")  # Table II: Km = 45.9 ug/mL = 45.9 mg/L

    # Covariate effects - Hashimoto 1994 Table II. Both fixed: the paper provides
    # point estimates and 95 percent CIs but does not separately estimate IIV
    # on either coefficient. Hypothesis testing in Table III rejected theta1 = 1
    # (linear weight, p << 0.005) and theta2 = 1 (no CBZ effect, p < 0.01).
    e_wt_vc_vmax      <- fixed(0.741); label("Shared power exponent of (WT/33) on Vc and Vmax (unitless)")  # Table II: theta1 = 0.741 (95% CI 0.665-0.817)
    e_conmed_cbz_vmax <- fixed(1.13);  label("Multiplicative factor on Vmax when CONMED_CBZ = 1 (unitless)")  # Table II: theta2 = 1.13 (95% CI 1.02-1.24)

    # IIV - Hashimoto 1994 Eq. 4 and Table II: omega_CL = 29.7 percent CV
    # (95% CI 25.2-33.6) applied additively to the steady-state effective
    # clearance CL = (Vmax - D/tau)/Km via CL_i = CL_typ * (1 + eta_CL).
    # Translated to log-normal IIV on Vmax for time-resolved ODE simulation:
    # since Km has no IIV in the source, IIV on Vmax preserves the implied
    # CL variability structure because at any concentration C the effective
    # CL = Vmax/(Km + C) is proportional to Vmax for fixed Km. The numerical
    # variance follows the (CV/100)^2 convention used by related Japanese
    # popPK models in this registry (Yukawa 1990 phenytoin). The vignette
    # Assumptions section discusses this re-parameterization.
    etalvmax ~ 0.0882  # Table II: omega_CL = 29.7 percent CV; 0.297^2 = 0.0882

    # Residual error - Hashimoto 1994 Table II: sigma = 17.8 percent CV
    # intra-individual (Eq. 5: SZC_ij = SZC_ij* * (1 + epsilon_ij) with
    # epsilon ~ N(0, sigma^2)).
    propSd <- 0.178; label("Proportional residual SD on Cc (fraction)")  # Table II: sigma = 17.8 percent (95% CI 14.5-20.7)
  })

  model({
    # Individual M-M parameters. Hashimoto 1994 Eqs. 2-4 (page 324):
    #   SIZE_i = 33 * (WT/33)^theta1
    #   V_i    = V * SIZE_i              -> equivalently V_typ_33kg * (WT/33)^theta1
    #   Vmax_i = Vmax * SIZE_i * theta2^CBZ -> Vmax_typ_33kg * (WT/33)^theta1 * theta2^CBZ
    # exp(lvc) and exp(lvmax) already encode the 33-kg typical value, so the
    # SIZE multiplier appears here only as (WT/33)^theta1.
    vc   <- exp(lvc)              * (WT / 33)^e_wt_vc_vmax
    vmax <- exp(lvmax + etalvmax) * (WT / 33)^e_wt_vc_vmax * e_conmed_cbz_vmax^CONMED_CBZ
    km   <- exp(lkm)

    # One-compartment PK with Michaelis-Menten elimination. The paper assumes
    # no significant absorption phase (page 324, samples drawn 2-12 h post-dose
    # at steady state) and treats each oral dose as an instantaneous bolus
    # into the central compartment; bioavailability F = 1 (page 326). The time
    # unit is days to match the published Vmax units (mg/d). Dose events should
    # specify cmt = central and amt in mg.
    Cc <- central / vc

    d/dt(central) <- -vmax * Cc / (km + Cc)

    Cc ~ prop(propSd)
  })
}
