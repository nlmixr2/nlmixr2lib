Philippe_2015_cyclosporine <- function() {
  description <- "Pediatric PK-PD-time-to-event model for oral cyclosporine in children with severe aplastic anemia (Philippe 2015). PK is a two-compartment model with first-order absorption, lag time, and linear elimination; absorption parameters (F, Tlag, ka) are fixed from the literature, and V1, V2, Cl, Q are allometrically scaled to body weight (reference 34 kg; fixed exponents 0.75 on clearance and 1 on volume). The pharmacodynamic interface model (Eq. 5) describes an effective concentration Ce driven by the predicted trough concentration Ctrough, with production active only when Ctrough lies inside an effective range (lower bound gamma1 = 87 ng/mL, upper bound gamma2 = 120 ng/mL) and first-order elimination at rate alpha. The instantaneous hazard of neutrophil response (Eq. 6) is lambda(t) = lambda0 * (1 + slope * Ce); cumhaz and sur are exposed as derived outputs. In this implementation the predicted Cc (multiplied by 1000 to convert mg/L to ng/mL) is used as the Ctrough input to the interface model; see vignette Assumptions and deviations for the full justification."
  reference <- "Philippe M, Henin E, Bertrand Y, Plantaz D, Goutelle S, Bleyzac N. Model-based determination of effective blood concentrations of cyclosporine for neutrophil response in the treatment of severe aplastic anemia in children. AAPS J. 2015;17(5):1157-1166. doi:10.1208/s12248-015-9779-8"
  vignette <- "Philippe_2015_cyclosporine"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L (= ug/mL); the pharmacodynamic bounds gamma1 and gamma2 are in ng/mL and Cc is rescaled inside model() by a factor of 1000 to compare with them"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg). Time-fixed in the source paper (baseline weight only); WT is the canonical column for body weight (baseline or time-varying).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 34 kg, equal to the pooled-cohort mean of the 23 included treatment courses (Table I). Allometric scaling per Methods Eq. 4 with fixed exponents 0.75 on Cl and Q and 1 on V1 and V2 (no estimated uncertainty was reported on these exponents; they were held at standard physiological values). The pooled-cohort weight range is 9.8-79.3 kg (Table I).",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "1-15 years (mean 8.5)",
    weight_range   = "9.8-79.3 kg (mean 34.0)",
    sex_female_pct = 43.5,
    race_ethnicity = "Not reported (French pediatric centres)",
    disease_state  = "Pediatric severe aplastic anemia (SAA) defined by two of: absolute neutrophil count < 0.5e9/L, platelets < 20e9/L, hemoglobin < 80 g/L with reticulocytes < 20e9/L. Patients received immunosuppressive therapy with anti-thymocyte globulin (160 mg/kg horse ATG or 18.75 mg/kg rabbit ATG) plus oral cyclosporine.",
    dose_range     = "Oral cyclosporine 5 mg/kg every 12 h (initial regimen); trough concentrations subsequently adjusted to an initial target of 150 ng/mL.",
    administration = "Oral, twice daily (per os)",
    regions        = "France (single-centre Lyon n = 22 courses + Grenoble n = 1 course; pediatric hematology/oncology units, retrospective 1998-2013)",
    notes          = "21 patients corresponding to 23 treatment courses (three patients received a second IST course after first-course failure and were modelled as independent individuals). 341 trough blood concentrations with median 14 measurements per patient (range 7-27). 15/23 courses (65.2%) achieved neutrophil response (ANC > 0.5e9/L on two consecutive occasions) with mean time-to-response 69 days (range 19-182). Three additional prospectively followed patients were used as the external validation cohort (Table III)."
  )

  ini({
    # Absorption parameters fixed from literature (Methods page 1158:
    # "Since no pharmacokinetic information in the absorption phase was
    # available, parameters related to the absorption (bioavailability,
    # lag time, and first-order absorption rate constant) were fixed,
    # according to the literature (24)"). Wrapped in fixed() to preserve
    # the held-constant status.
    lka     <- fixed(log(0.829)); label("First-order oral absorption rate constant ka (1/h)")  # Table II: ka = 0.829 1/h (fixed)
    ltlag   <- fixed(log(0.648)); label("Oral absorption lag time Tlag (h)")                   # Table II: Tlag = 0.648 h (fixed)
    lfdepot <- fixed(log(0.386)); label("Oral bioavailability F (dimensionless)")              # Table II: F = 0.386 (fixed)

    # Two-compartment disposition; typical values at the reference weight
    # of 34 kg (Methods Eq. 4, where BWpop is the cohort mean from
    # Table I). All four disposition parameters are allometrically scaled
    # to (WT/34) with fixed exponents 0.75 on Cl and Q and 1 on V1 and V2.
    lcl <- log(7.2);  label("Clearance Cl/F at WT = 34 kg (L/h)")                              # Table II: Cl = 7.2 L/h at 34 kg (rse 13%)
    lvc <- log(37.5); label("Central volume V1/F at WT = 34 kg (L)")                           # Table II: V1 = 37.5 L at 34 kg (rse 25%)
    lq  <- log(2.7);  label("Inter-compartment clearance Q/F at WT = 34 kg (L/h)")             # Table II: Q = 2.7 L/h at 34 kg (rse 26%)
    lvp <- log(1690); label("Peripheral volume V2/F at WT = 34 kg (L)")                        # Table II: V2 = 1690 L at 34 kg (rse 25%)

    # Allometric exponents fixed at the standard physiological values per
    # Methods Eq. 4 (no uncertainty reported). Cl and Q share exponent
    # 0.75; V1 and V2 share exponent 1. Recorded as four parameters
    # (one per canonical disposition parameter) so the lint accepts each
    # covariate-effect name.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/34) for Cl (unitless)")          # Methods Eq. 4: exponent fixed to 0.75 for clearance parameters
    e_wt_q  <- fixed(0.75); label("Allometric exponent on (WT/34) for Q (unitless)")           # Methods Eq. 4: exponent fixed to 0.75 for clearance parameters (same value as e_wt_cl)
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/34) for V1 (Vc) (unitless)")     # Methods Eq. 4: exponent fixed to 1 for volume parameters
    e_wt_vp <- fixed(1);    label("Allometric exponent on (WT/34) for V2 (Vp) (unitless)")     # Methods Eq. 4: exponent fixed to 1 for volume parameters (same value as e_wt_vc)

    # Inter-individual variability (Monolix exponential model, Eq. 1).
    # The reported IIV% CV is converted to log-normal variance via
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.09405  # Table II: Cl IIV = 31.4 %CV; log(1 + 0.314^2) = 0.09405
    etalvc ~ 0.08618  # Table II: V1 IIV = 30   %CV; log(1 + 0.30^2)  = 0.08618
    etalq  ~ 0.36363  # Table II: Q  IIV = 66.2 %CV; log(1 + 0.662^2) = 0.36363
    etalvp ~ 0.08618  # Table II: V2 IIV = 30   %CV; log(1 + 0.30^2)  = 0.08618

    # Residual error: combined additive + proportional model per Methods
    # page 1158 ("combined model wij = sqrt(a^2 + b^2 * Cpred^2)"), with
    # a and b given in Table II.
    addSd  <- 0.03; label("Additive residual SD on Cc (mg/L = ug/mL)")                          # Table II: a = 0.03 mg/L (rse 9%)
    propSd <- 0.25; label("Proportional residual SD on Cc (fraction)")                          # Table II: b = 25% (rse 9%)

    # Time-to-response and interface-model parameters (Table II). No IIV
    # could be estimated on these parameters: "No IIV could be estimated
    # on these parameters (IIV fixed to 0 for all the parameters)"
    # (Results page 1162). Typical values were estimated, not fixed -
    # they carry %rse in Table II - so they are NOT wrapped in fixed();
    # only the IIV is omitted.
    lambda0 <- 0.0027; label("Baseline hazard of response without treatment (1/day)")           # Table II: lambda0 = 0.0027 day-1 (rse 58%)
    alphaIE <- 0.028;  label("Elimination rate constant of effective concentration (1/day)")    # Table II: alpha = 0.028 day-1 (rse 107%); implied half-life ln(2)/alpha = 25 days
    slopeIE <- 17.2;   label("Hazard slope on Ce (mL/ng)")                                      # Table II: slope = 17.2 mL/ng (rse 102%)
    gamma1  <- 87;     label("Lower bound of optimal CsA trough concentration range (ng/mL)")   # Table II: gamma1 = 87 ng/mL (rse 3%)
    gamma2  <- 120;    label("Upper bound of optimal CsA trough concentration range (ng/mL)")   # Table II: gamma2 = 120 ng/mL (rse 3%)
  })
  model({
    # Individual PK parameters (allometric scaling per Methods Eq. 4 with
    # fixed exponents 0.75 / 1 as recorded in ini()).
    cl <- exp(lcl + etalcl) * (WT / 34)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 34)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 34)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 34)^e_wt_vp
    ka         <- exp(lka)
    tlag_depot <- exp(ltlag)
    fdepot     <- exp(lfdepot)

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Oral absorption with lag and first-order ka.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag_depot

    # Cyclosporine whole-blood concentration. Dose in mg, vc in L gives
    # Cc in mg/L (= ug/mL), matching Fig. 3 of the paper.
    Cc <- central / vc

    # Convert Cc to ng/mL for comparison with the interface bounds
    # gamma1 and gamma2 reported in ng/mL throughout the paper.
    Cc_ngmL <- 1000 * Cc

    # Indicator variables H1 and H2 per Eq. 5: H1 = 1 if Ctrough > gamma1,
    # H2 = 1 if Ctrough < gamma2, both 0 otherwise. We use the predicted
    # Cc as the Ctrough input (see vignette Assumptions and deviations).
    H1 <- (Cc_ngmL > gamma1)
    H2 <- (Cc_ngmL < gamma2)

    # Interface model (Eq. 5):
    #   dCe/dt = -alpha * Ce + sqrt((Ctrough - gamma1) * (gamma2 - Ctrough) * H1 * H2)
    # The paper expresses dCe/dt in (ng/mL)/day (alpha is in 1/day).
    # This model's TIME axis is in hours (consistent with the PK ka in
    # 1/h and Tlag in h); we therefore divide alpha and the sqrt
    # production term by 24 to express the equation per hour.
    alphaIE_h    <- alphaIE / 24
    production_h <- sqrt((Cc_ngmL - gamma1) * (gamma2 - Cc_ngmL) * H1 * H2) / 24
    d/dt(effect) <- -alphaIE_h * effect + production_h
    effect(0)    <- 0

    # Time-to-event hazard (Eq. 6):
    #   lambda(t) = lambda0 * (1 + slope * Ce(t))
    # lambda0 is converted from 1/day to 1/h to match the model TIME
    # axis (hours). slope multiplied by Ce (ng/mL * mL/ng) is
    # dimensionless and needs no conversion.
    lambda0_h <- lambda0 / 24
    hazard    <- lambda0_h * (1 + slopeIE * effect)

    # Cumulative hazard and survival probability (Eq. 7). For t in
    # hours, cumhaz integrates the per-hour hazard so sur = exp(-cumhaz)
    # is the probability of no neutrophil response up to time t.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)

    # Concentration observation model (combined additive + proportional;
    # Methods page 1158).
    Cc ~ add(addSd) + prop(propSd)
  })
}
