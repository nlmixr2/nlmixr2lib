Mo_2014_azd1208_mouse <- function() {
  description <- "Preclinical (mouse, MOLM-16 acute myeloid leukemia xenograft). Lifespan-based tumor growth inhibition (LS TGI) model of the oral pan-Pim kinase inhibitor AZD1208 with a non-cycle-specific Emax drug effect (Mo 2014 Case Study 2, Tables 5 and 6). Tumor cells have a common lifespan td_div at the end of which they divide. The tumors grew only exponentially, so Mo 2014 set the size threshold wth to a value far above the observed volumes; in that limit the division efficiency p(w) of Eq. 13 is the constant p0, and the recursive lifespan product of Eq. 22 reduces to p0^n(t), where n(t) is the number of elapsed lifespans. Following the authors' published MATLAB implementation, n(t) = round(t / td_div) once t >= td_div (and 0 before), not the floor INT(t / td_div) printed in Equations 7 and 22 (see the vignette). AZD1208 removes proliferating cells at the saturable rate Emax * C / (EC50 + C), the Emax replacement of k2 * C in Eq. 18, into a non-proliferating apoptotic pool that dies after td_apo (Eq. 24, encoded with rxode2 delay()); the surviving fraction of the proliferating cells, exp(-integral of the kill rate), is carried as the cumulative kill hazard cumhaz. One-compartment oral AZD1208 PK with the elimination rate constant set equal to the absorption rate constant (Table 5). Time is in days, doses in mg/kg into depot, plasma concentration in ng/mL and tumor volume in mm^3. Solve with method = 'dop853' and maxsteps = 1e7."
  reference <- paste(
    "Mo G, Gibbons F, Schroeder P, Krzyzanski W. Lifespan based",
    "pharmacokinetic-pharmacodynamic model of tumor growth inhibition by",
    "anticancer therapeutics. PLoS One. 2014;9(10):e109747.",
    "doi:10.1371/journal.pone.0109747"
  )
  vignette <- "Mo_2014_lifespan_tgi"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "AZD1208", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "AZD1208", units = "mg/kg", specimen = "plasma", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells M", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "apoptotic (non-proliferating) tumor cells A",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative drug-induced kill hazard of proliferating cells",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (MOLM-16 human acute myeloid leukemia xenograft)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    disease_state = "subcutaneous MOLM-16 human acute myeloid leukemia xenograft",
    dose_range = "AZD1208 0.3, 1, 3, 10 or 30 mg/kg orally (tumor study, Figure 9B); single 3, 10 or 30 mg/kg oral doses (PK study, 3 mice per dose)",
    notes = "Data provided by AstraZeneca (Mo 2014 Case Study 2; Keeton et al. Blood 2014;123:905-913). Day 1 is the first dose. The control data were too sparse to estimate the doubling time, which was fixed to the 3-day in vitro doubling time of MOLM-16 cells. Mean data fitted with MATLAB nlinfit. Initial tumor volume 170 mm^3, the first measured value (Figure 9 caption). The dosing regimen of the efficacy study is not restated in Mo 2014; see the vignette."
  )

  ini({
    # Lifespan TGI parameters: Mo 2014 Table 6.
    ltd_div <- fixed(log(3)); label("Tumor cell lifespan (doubling time) T (day)") # Table 6: T = 3 day, set to the MOLM-16 in vitro doubling time
    ltd_apo <- fixed(log(1.92)); label("Apoptosis duration TA (day)") # Table 6: TA = 1.92 day, grid search
    lp0 <- log(1.57); label("Cell division efficiency p0 (unitless)") # Table 6: p0 = 1.57
    lkin0 <- log(27.8); label("Tumor production rate before the experiment kin0 (mm^3/day)") # Table 6: kin0 = 27.8 mm^3/day
    lemax <- log(0.159); label("Maximum drug-induced kill rate Emax (1/day; printed as unitless)") # Table 6: Emax = 0.159
    lec50 <- log(182); label("AZD1208 plasma concentration at half-maximal kill rate EC50 (ng/mL)") # Table 6: EC50 = 182 ng/mL
    lrbase_tumor <- fixed(log(170)); label("Tumor volume at the start of the experiment w0 (mm^3)") # Figure 9 caption: initial tumor volume fixed to 170 mm^3

    # AZD1208 PK: Mo 2014 Table 5 (one-compartment, first-order absorption;
    # kel set equal to ka, Case Study 2 text).
    lka <- log(5.52); label("AZD1208 absorption rate constant ka, also the elimination rate constant (1/day)") # Table 5: ka = 5.52 1/day
    lvc <- log(4.86); label("AZD1208 apparent volume of distribution V/F (L/kg)") # Table 5: V/F = 4.86 L/kg

    # Mo 2014 fitted the mean profiles by unweighted nonlinear least
    # squares (MATLAB nlinfit) and reports no residual-error estimate.
    addSd_tumor_vol <- fixed(0); label("Additive residual error on tumor volume (mm^3; not reported)") # not reported (nlinfit least squares)
  })

  model({
    td_div <- exp(ltd_div)
    td_apo <- exp(ltd_apo)
    p0 <- exp(lp0)
    kin0 <- exp(lkin0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    rbase <- exp(lrbase_tumor)
    ka <- exp(lka)
    kel <- ka
    vc <- exp(lvc)

    # AZD1208 PK (dose in mg/kg into depot; F is absorbed into V/F;
    # mg/L x 1000 = ng/mL)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc * 1000

    # Tumor volume w = M + A (Eq. 17)
    tumor_vol <- cycling_cells + damaged_cells1
    # Emax drug effect replacing k2 * C in Eq. 18 (Case Study 2 text)
    eff <- emax * Cc / (ec50 + Cc)

    # Number of elapsed lifespans in the lifespan product, as in the
    # supplement S5 code: round(t / T) from t = T on, 0 before.
    n_div <- 0
    if (t >= td_div) {n_div <- round(t / td_div)}

    # Apoptotic outflow E(C(t - TA)) M(t - TA) (Eq. 24)
    c_apo <- 0
    m_apo <- 0
    if (t >= td_apo) {
      c_apo <- delay(central, td_apo) / vc * 1000
      m_apo <- delay(cycling_cells, td_apo)
    }
    eff_apo <- emax * c_apo / (ec50 + c_apo)

    # Proliferating pool (Eq. 22 with p(w) = p0), apoptotic pool (Eq. 24)
    # and the surviving fraction exp(-cumhaz) of Eq. 22
    cycling_cells(0) <- rbase
    d/dt(cycling_cells) <- kin0 * (p0 - 1) * exp(-cumhaz) * p0^n_div - eff * cycling_cells
    d/dt(damaged_cells1) <- eff * cycling_cells - eff_apo * m_apo
    d/dt(cumhaz) <- eff

    tumor_vol ~ add(addSd_tumor_vol)
  })
}
