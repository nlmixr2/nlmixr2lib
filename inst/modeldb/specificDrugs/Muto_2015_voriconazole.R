Muto_2015_voriconazole <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption (lag time, oral bioavailability) and parallel linear plus time-dependent Michaelis-Menten elimination for voriconazole in 21 immunocompromised Japanese pediatric subjects (Muto 2015). Vmax declines with time after the first dose toward Vmax * (1 - Vmax_inh) with half-time T50; the maximum inhibition fraction Vmax_inh is fixed to 1 (full inhibition) for CYP2C19 heterozygous-extensive-metabolizer or poor-metabolizer subjects and modeled on the logit scale otherwise. Allometric scaling on all clearances (exponent 0.75) and all volumes (exponent 1) to a 70 kg reference; oral bioavailability F1 is modeled on the logit scale with a Manly-transformed log-normal random effect."
  reference   <- "Muto C, Shoji S, Tomono Y, Liu P. Population pharmacokinetic analysis of voriconazole from a pharmacokinetic study with immunocompromised Japanese pediatric subjects. Antimicrob Agents Chemother. 2015;59(6):3216-3223. doi:10.1128/AAC.04993-14"
  vignette    <- "Muto_2015_voriconazole"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling reference 70 kg per Muto 2015 Methods (page 3217) and Appendix equations: clearances (Vmax1, CL, Q) scale with (WT/70)^0.75; volumes (V2, V3) scale with WT/70. Cohort body-weight range 11.5-55.2 kg (Table 2, n = 21); typical-value tables for 20 kg and 50 kg pediatric subjects are reported in Table 4.",
      source_name        = "wt"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator (heterozygous extensive metabolizer, HEM, in Muto 2015 nomenclature)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive / ultra-rapid metabolizer; CYP2C19_PM = 0 must also be 0 for the EM/UM reference)",
      notes              = "Muto 2015 calls the heterozygous *1/*2 (or similar) genotype a heterozygous extensive metabolizer (HEM); the canonical CYP2C19_IM phenotype indicator is the same construct. In Muto 2015 Final Model, the CYP2C19 effect on Vmax_inh dichotomizes UM/EM versus HEM/PM, with HEM/PM subjects fixed to Vmax_inh = 1 (100% inhibition of the saturable elimination arm at large t). The model () block therefore uses the OR of CYP2C19_IM and CYP2C19_PM to switch between the logit-estimated typical value and the fixed value of 1. Cohort distribution: 9 EM, 10 HEM, 2 PM (Table 2).",
      source_name        = "CYP2C19 genotype (EM / HEM / PM)"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator (homozygous loss-of-function alleles, e.g., *2/*2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive / ultra-rapid / intermediate metabolizer; pairs with CYP2C19_IM)",
      notes              = "Muto 2015 reports 2 of 21 subjects (9.5%) as CYP2C19 poor metabolizers (PM). The covariate effect on Vmax_inh treats HEM and PM identically in the final model (paper Discussion page 3221: 'The CYP2C19 PM subjects were combined with the CYP2C19 HEM subjects during the covariate evaluation process'); both groups share Vmax_inh = 1.",
      source_name        = "CYP2C19 genotype (EM / HEM / PM)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 21L,
    n_studies       = 1L,
    n_observations  = 276L,
    age_range       = "3-14 years",
    age_median      = "10 years",
    weight_range    = "11.5-55.2 kg",
    weight_median   = "31.5 kg",
    bmi_range       = "11.9-22.9 kg/m^2",
    sex_female_pct  = 57.1,
    race_ethnicity  = "Japanese (single ethnicity)",
    disease_state   = "Immunocompromised pediatric subjects at high risk for systemic fungal infection (largely hematological-malignancy and post-transplant indications). Concomitant medications and complicated background therapy are common and contribute to the high IIV in oral bioavailability and exposure (Muto 2015 Discussion).",
    dose_range      = "Intravenous loading 9 mg/kg q12h (children and adolescents <50 kg) or 6 mg/kg (adolescents >=50 kg) on day 1; intravenous maintenance 8 mg/kg q12h (children and adolescents <50 kg) or 4 mg/kg q12h (adolescents >=50 kg) on days 2-7; oral suspension 9 mg/kg q12h capped at 350 mg per dose (children and adolescents <50 kg) or 200 mg q12h (adolescents >=50 kg) on days 8-14. IV infusions at 3 mg/kg/h. Oral doses at least 1 h before or after a meal.",
    regions         = "Six centers in Japan",
    cyp2c19_status  = c(EM_pct = 42.9, HEM_pct = 47.6, PM_pct = 9.5),
    notes           = "Population from Table 2 (demographics) and Table 1 (dosing regimens). The current popPK analysis applies normal-inverse-Wishart Bayesian priors derived from the earlier non-Japanese pediatric + adult model of Friberg 2012 (Antimicrob Agents Chemother 56:3032-3042) to stabilize estimation in the small Japanese cohort (n = 21). Of the 21 subjects, 18 contributed oral-period concentrations (3 subjects discontinued at the i.v.-to-oral switch); 152 i.v. and 124 oral concentration records contributed to the final model fit."
  )

  ini({
    # Structural parameters - Muto 2015 Table 3 'Estimate for theta' column.
    # Reference body weight is 70 kg (Appendix equations); typical-value tables
    # in Table 4 are reported at 50 kg and 20 kg after allometric scaling.

    # Saturable elimination (time-dependent Vmax auto-inhibition).
    lkm           <- log(0.922); label("Michaelis-Menten constant Km (ug/mL)")                                                 # Table 3 Km estimate 0.922 (RSE 30%)
    lvmax         <- log(118);   label("Maximum elimination rate Vmax at 1 h after first dose, Vmax,1 (mg/h per 70 kg)")       # Table 3 Vmax,1 estimate 118 (RSE 14%)
    logitvmaxinh  <- 2.61;       label("Logit of maximum fractional Vmax inhibition Vmax_inh for CYP2C19 UM/EM (unitless)")    # Table 3 Vmax_inh logit estimate 2.61 (RSE 19%); expit(2.61) = 0.932
    lt50          <- log(2.45);  label("Half-time of the time-dependent Vmax decay T50 (h)")                                   # Table 3 T50 estimate 2.45 (RSE 6.3%)
    lvmaxscale    <- log(1.25);  label("Scaling factor for the shared Km/Vmax,1 random effect on Vmax,1 (theta_Vmax,scale; unitless)") # Table 3 theta_Vmax,scale estimate 1.25 (RSE 12%)

    # Linear disposition.
    lcl <- log(6.02); label("Linear clearance CL (L/h per 70 kg)")                # Table 3 CL estimate 6.02 (RSE 11%)
    lvc <- log(75.0); label("Central volume of distribution V2 (L per 70 kg)")    # Table 3 V2 estimate 75.0 (RSE 3.2%)
    lvp <- log(101);  label("Peripheral volume of distribution V3 (L per 70 kg)") # Table 3 V3 estimate 101 (RSE 6.1%)
    lq  <- log(24.6); label("Inter-compartmental clearance Q (L/h per 70 kg)")    # Table 3 Q estimate 24.6 (RSE 4.4%)

    # Absorption.
    lka         <- log(1.38);  label("First-order absorption rate constant ka (1/h)")              # Table 3 ka estimate 1.38 (RSE 14%)
    ltlag       <- log(0.121); label("Absorption lag time Alag (h)")                               # Table 3 Alag estimate 0.121 (RSE 2.8%)
    logitfdepot <- 0.597;      label("Logit of oral bioavailability F1 for the depot (unitless)")  # Table 3 F1 logit estimate 0.597 (RSE 13%); expit(0.597) = 0.645

    # Allometric exponents (fixed at canonical values per Muto 2015 Methods
    # page 3217: 'All clearance terms ... were scaled allometrically, using
    # weight to a power of 0.75, and all volume terms ... were scaled
    # allometrically, using weight to a power of 1.').
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on clearances Vmax,1/CL/Q (unitless)") # Methods page 3217 (fixed)
    e_wt_vc <- fixed(1);    label("Allometric exponent of WT on volumes V2/V3 (unitless)")         # Methods page 3217 (fixed)

    # Manly transformation parameter for the F1 random effect.
    # ETATR = ((exp(etalogitfdepot))^lambda - 1) / lambda;
    # logit(F1,i) = logit(F1) + ETATR. See Table 3 footnote f and Table 4
    # footnote d for the transform definition.
    lbcflambda <- log(0.330); label("Manly-transformation power for the F1 random effect, theta_BC-F (unitless)") # Table 3 theta_BC-F estimate 0.330 (RSE 23%)

    # IIV - Muto 2015 Table 3 'Estimate for Omega SD' column.
    # Reported as omega SDs (log-scale); variances are (omega SD)^2.
    # Block correlations on Km/Vmax, V3, CL, Q come from Table 3 estimates:
    #   Corr(Km,V3)=-0.52, Corr(Km,CL)=0.26, Corr(V3,CL)=0.15,
    #   Corr(Km,Q) =-0.61, Corr(V3,Q) =0.88, Corr(CL,Q) =0.097.
    # The shared random effect etalkm_vmax drives both Km and Vmax,1
    # (with Vmax,1 scaled by exp(lvmaxscale)); declaring the eta with the
    # 'etal<p1>_<p2>' shared-eta naming convention so checkModelConventions
    # accepts it without flagging the missing 1-to-1 lvmax pairing.
    etalkm_vmax + etalvp + etalcl + etalq ~ c(
      1.8496,
      -0.5544, 0.6147,
      0.2462, 0.0819, 0.4844,
      -0.3602, 0.2995, 0.0293, 0.1884
    )                                                                          # Table 3: SDs 1.36, 0.784, 0.696, 0.434 and listed correlations
    etalvc         ~ 0.020164  # Table 3: SD 0.142, var = 0.142^2 = 0.020164
    etalka         ~ 0.799236  # Table 3: SD 0.894, var = 0.894^2 = 0.799236
    etalogitfdepot ~ 2.8561    # Table 3: SD 1.69, var = 1.69^2 = 2.8561 (Manly-transformed; see lbcflambda)

    # Residual error - Muto 2015 Table 3 reports a single residual error
    # estimate of 0.239 (RSE 5.8%) modeled as additive on log-transformed
    # concentrations (Methods page 3218: 'Within-subject variability ...
    # was modeled as additive errors on the log-transformed concentrations
    # (analogous to the proportional-error model on the untransformed
    # concentrations).'). expSd is the canonical SD name for the
    # ~ lnorm(...) log-normal residual.
    expSd <- 0.239; label("Log-scale residual SD (additive on log-transformed concentrations)") # Table 3 residual error 0.239 (RSE 5.8%)
  })

  model({
    # Derived covariate terms.
    # WT-based allometric scaling factors (reference 70 kg).
    wt_cl <- (WT / 70)^e_wt_cl
    wt_vc <- (WT / 70)^e_wt_vc

    # CYP2C19 HEM/PM combined indicator: 1 if either CYP2C19_IM or
    # CYP2C19_PM is 1 (mutually exclusive in the source paper), 0 for
    # the EM/UM reference. Boolean OR encoded as 1 - (1 - a)*(1 - b)
    # so rxode2 evaluates the expression without an ifelse.
    cyp_hempm <- 1 - (1 - CYP2C19_IM) * (1 - CYP2C19_PM)

    # Individual parameters.
    # Shared random effect etalkm_vmax drives Km directly and Vmax,1 via
    # the multiplicative scale exp(lvmaxscale). theta_Vmax,scale = 1.25
    # means Vmax,1's effective IIV SD on the log scale is 1.36 * 1.25 = 1.70.
    vmaxscale <- exp(lvmaxscale)
    km    <- exp(lkm   + etalkm_vmax)
    vmax1 <- exp(lvmax + etalkm_vmax * vmaxscale) * wt_cl
    cl    <- exp(lcl   + etalcl) * wt_cl
    q     <- exp(lq    + etalq)  * wt_cl
    vc    <- exp(lvc   + etalvc) * wt_vc
    vp    <- exp(lvp   + etalvp) * wt_vc
    ka    <- exp(lka   + etalka)
    tlag  <- exp(ltlag)
    t50   <- exp(lt50)

    # Time-dependent Vmax (auto-inhibition decay).
    # Vmax(t) = Vmax,1 * (1 - Vmax_inh * (t-1) / ((t-1) + (T50-1)));
    # t is time after the first dose. At t = 1 h, Vmax = Vmax,1
    # (reference); as t grows, Vmax decays toward Vmax,1 * (1 - Vmax_inh).
    # Vmax_inh = 1 (full inhibition) for CYP2C19 HEM/PM subjects; otherwise
    # logit-estimated typical value expit(logitvmaxinh) = 0.932 for the
    # CYP2C19 UM/EM reference. The formula is the model as written in the
    # Muto 2015 Appendix and is intended for t >= 0 (PK simulation starts
    # at the first dose); the (t-1) factor in the numerator yields the
    # paper's stated identity Vmax(1 h) = Vmax,1.
    vmax_inh <- (1 - cyp_hempm) * expit(logitvmaxinh) + cyp_hempm
    vmax     <- vmax1 * (1 - vmax_inh * (t - 1) / ((t - 1) + (t50 - 1)))

    # Bioavailability - Manly-transformed random effect on logit(F1):
    # ETATR = ((exp(eta))^lambda - 1) / lambda. lambda = exp(lbcflambda)
    # = 0.330; at lambda -> 0 the transform reduces to ETATR = eta (the
    # Box-Cox limit), at lambda = 1 ETATR = exp(eta) - 1. The transform is
    # applied as an additive shift on the logit scale; see Muto 2015 Table
    # 3 footnote f and Table 4 footnote d for the definition.
    bcflambda  <- exp(lbcflambda)
    manly_eta  <- ((exp(etalogitfdepot))^bcflambda - 1) / bcflambda
    fdepot     <- expit(logitfdepot + manly_eta)

    # Micro-constants for the two-compartment disposition.
    k12 <- q / vc
    k21 <- q / vp

    # ODE system - two-compartment with first-order oral absorption,
    # parallel linear (CL) and Michaelis-Menten (Vmax(t), Km) elimination
    # from the central compartment. Concentration Cc = central / vc
    # (mg in central, L in vc, so Cc in mg/L = ug/mL).
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central - vmax * Cc / (km + Cc) -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability applies to the depot compartment only (i.v. doses go
    # directly into central with F = 1). Absorption lag time is set on the
    # depot via alag().
    f(depot)    <- fdepot
    alag(depot) <- tlag

    Cc ~ lnorm(expSd)
  })
}
