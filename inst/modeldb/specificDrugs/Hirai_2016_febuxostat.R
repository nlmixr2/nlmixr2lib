Hirai_2016_febuxostat <- function() {
  description <- "Semi-mechanistic PK-PD model for the hypouricemic effect of febuxostat (non-purine xanthine oxidase inhibitor) in hyperuricemic patients with or without chronic kidney disease. Combines a 2-compartment first-order absorption PK model for febuxostat (literature-digitized, WinNonlin-fit) with (i) a reversible mixed-type inhibition of xanthine oxidase (parameters transferred from Takano 2005 in vitro bovine XO study), (ii) a 2-compartment endogenous uric acid disposition model (adapted from Scott 1969), and (iii) a 1-compartment endogenous xanthine model. Renal uric acid clearance is a power function of creatinine clearance (Tykarski 1991). Non-renal uric acid clearance is assumed constant. All parameters are literature-derived and fixed (deterministic mechanistic model; no population-PK fit). Reproduces Figure 2 of Hirai 2016 (10 mg febuxostat QD in a normouricemic 60 kg patient with CrCl 100 mL/min)."
  reference <- "Hirai T, Kimura T, Echizen H. Modeling and simulation for estimating the influence of renal dysfunction on the hypouricemic effect of febuxostat in hyperuricemic patients due to overproduction or underexcretion of uric acid. Biol Pharm Bull. 2016;39(6):1013-1021. doi:10.1248/bpb.b15-01031"
  vignette <- "Hirai_2016_febuxostat"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  paper_specific_compartments <- c("uapool1", "uapool2", "xanthine")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "febuxostat", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "febuxostat", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "febuxostat", units = "mg", specimen = "plasma", verified = FALSE),
    uapool1     = list(analyte = "uric acid", units = "mg", specimen = "urine", verified = FALSE),
    uapool2     = list(analyte = "uric acid", units = "mg", specimen = "plasma", verified = FALSE),
    xanthine    = list(analyte = "xanthine", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Enters as a multiplicative factor on all body-weight-normalized volumes (Vc, Vp for febuxostat; Vd1, Vd2 for uric acid; Vd for xanthine) and thereby on all body-weight-normalized clearances (CL, Q, kel_R_UA). Reference values used by Hirai 2016 across data sets: 60 kg for Japanese and 70 kg for Caucasian patients (Table 1 and Methods p. 1015).",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault, actual body weight; the paper treats estimated GFR mL/min per 1.73 m^2 and CLcr mL/min as exchangeable per Results p. 1017)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Sole driver of renal uric acid clearance via Hirai 2016 Eq. 9 (Tykarski 1991): CL_R_UA (mL/min) = 1.23 * CRCL^0.433. Renal impairment classes reported by the paper: normal (>80 mL/min), mild CKD (50-80), moderate CKD (30-49), severe CKD (<29). Febuxostat AUC was reported to be elevated by 48 percent in moderate CKD and 74 percent in severe CKD compared with normal renal function (Results p. 1017), but the paper does not provide per-CKD-class point estimates for the febuxostat PK parameters -- see vignette Errata.",
      source_name        = "CLcr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 735,
    n_studies      = 5,
    age_range      = "adults (specific ages not reported per data set)",
    weight_range   = "59-88 kg across data sets (Table 1); reference weights 60 kg (Japanese) and 70 kg (Caucasian)",
    sex_female_pct = NA_real_,
    disease_state  = "healthy volunteers (n = 156) and hyperuricemic patients (n = 579) across normal renal function and mild / moderate / severe chronic kidney disease",
    dose_range     = "10-240 mg febuxostat orally once daily for at least 7 days (Table 1 across data sets 1-39)",
    regions        = "Japan and United States",
    notes          = "Data set assembled from Medline searches, the Japanese Pharmaceuticals and Medical Devices Agency Common Technical Documents summary, and the febuxostat interview form and prescribing information issued by Teijin Pharma (Methods p. 1013-1014, references 9-12, 18-20). 39 data sets in total, listed in Table 1. Mean plasma concentrations of febuxostat were digitized from published figures using UnGraph 5 (BIOSOFT) and fit with an open two-compartment first-order-absorption model in Phoenix WinNonlin (Certara). The febuxostat PK parameters below are the range of point-estimate values reported by the paper for the four renal function classes (Vd1: 0.23-0.30 L/kg, k10: 0.25-0.46/h, k12: 0.13-0.32/h, k21: 0.16-0.20/h; Results p. 1018). This model encodes the NORMAL-renal-function endpoints of those ranges (matching the representative 60 kg patient in Figure 2). Fu of febuxostat was reported as 0.9 percent (normal / mild), 0.8 percent (moderate CKD), and 1.2 percent (severe CKD); the normal / mild value is used here."
  )

  ini({
    # ============================================================
    # Febuxostat 2-compartment PK with first-order absorption
    # Source: Hirai 2016 Results p. 1018 (WinNonlin refit of digitized
    # literature concentration-time data from Mayer 2005 [ref 9] and
    # PMDA CTD [ref 12]). Foral and ka are assumed universally
    # (Methods p. 1014); Vc, CL, Q, Vp are BW-normalized and given as
    # ranges across four renal-function classes. This file encodes the
    # NORMAL-renal-function endpoints (Vd1 = 0.30 L/kg, k10 = 0.46/h,
    # k12 = 0.22/h [midpoint of 0.13-0.32; no monotone CKD trend],
    # k21 = 0.20/h) and reparameterizes as CL, Q, Vc, Vp:
    #   CL = k10 * Vd1 = 0.46 * 0.30 = 0.138 L/h/kg
    #   Q  = k12 * Vd1 = 0.22 * 0.30 = 0.066 L/h/kg
    #   Vp = Q / k21   = 0.066 / 0.20 = 0.33 L/kg
    # The paper multiplies Foral separately in the absorption input
    # (Eq. 1), so Vd1/Vd2 are TRUE (not V/F) volumes; encoded here via
    # f(depot) = Foral so central-compartment concentrations divide by
    # the true Vc.
    # ============================================================
    lka       <- fixed(log(2.18));    label("Febuxostat absorption rate constant ka (1/h)")                       # Hirai 2016 Methods p. 1014: assumed ka = 2.18 h^-1 (Mayer 2005 ref 9, PMDA CTD ref 12); no BSV
    lfdepot   <- fixed(log(0.65));    label("Febuxostat oral bioavailability Foral (fraction)")                    # Hirai 2016 Methods p. 1014: Foral = 65 percent (Mayer 2005 ref 9, PMDA CTD ref 12 human mass-balance); no BSV
    lvc       <- fixed(log(0.30));    label("Febuxostat central volume of distribution Vd1 (L/kg)")                # Hirai 2016 Results p. 1018: range 0.23-0.30 L/kg; upper endpoint (normal renal function) used
    lcl       <- fixed(log(0.138));   label("Febuxostat clearance CL (L/h/kg) = k10 * Vd1")                        # Hirai 2016 Results p. 1018: k10 range 0.25-0.46/h; upper endpoint (normal renal function) used, times Vd1 = 0.30 L/kg
    lq        <- fixed(log(0.066));   label("Febuxostat intercompartmental clearance Q (L/h/kg) = k12 * Vd1")      # Hirai 2016 Results p. 1018: k12 range 0.13-0.32/h; midpoint 0.22/h used (no monotone CKD trend reported), times Vd1 = 0.30 L/kg
    lvp       <- fixed(log(0.33));    label("Febuxostat peripheral volume of distribution Vp (L/kg) = Q / k21")    # Hirai 2016 Results p. 1018: k21 range 0.16-0.20/h; upper endpoint (normal renal function) used, so Vp = 0.066 / 0.20 = 0.33 L/kg
    fu_febx   <- fixed(0.009);        label("Febuxostat plasma unbound fraction (dimensionless)")                  # Hirai 2016 Results p. 1018 (Mayer 2005 ref 9): fu = 0.9 percent for normal / mild renal function, 0.8 percent for moderate CKD, 1.2 percent for severe CKD

    # ============================================================
    # Xanthine oxidase reversible mixed-type inhibition kinetic
    # constants for febuxostat. The paper adopts the bovine XO values
    # of Takano 2005 (ref 13) and assumes they translate to human XO.
    # Enzyme kinetics: v = [S]*Vmax/(Km + [S]) at baseline;
    # v' = [S]*Vmax/(Km*(1 + [I]/Ki) + [S]*(1 + [I]/Ki_prime)) with
    # febuxostat present. Only the ratio v'/v is used as a
    # dimensionless multiplier on UA formation (Hirai 2016 Eq. 4).
    # ============================================================
    km        <- fixed(2.7);          label("Michaelis-Menten constant Km of bovine XO for xanthine (uM)")         # Hirai 2016 Methods p. 1014 (Takano 2005 ref 13): Km = 2.7 uM
    ki        <- fixed(0.6);          label("Competitive inhibition constant Ki of febuxostat on XO (nM)")         # Hirai 2016 Methods p. 1014 (Takano 2005 ref 13): Ki = 0.6 nM
    ki_prime  <- fixed(3.1);          label("Uncompetitive inhibition constant Ki' of febuxostat on XO (nM)")      # Hirai 2016 Methods p. 1014 (Takano 2005 ref 13): Ki' = 3.1 nM

    # ============================================================
    # Endogenous uric acid two-compartment disposition (Hirai 2016
    # Eqs. 5-8). The kinetic parameters were transferred from Scott
    # 1969 (ref 14), a human UA pool study on healthy subjects and
    # patients (approximate mean BW ~80 kg). Because Hirai 2016
    # applies these values to non-Caucasian patients as well, the
    # volumes are body-weight-normalized (L/kg).
    # ============================================================
    lvd_ua_1     <- fixed(log(0.25));  label("Uric acid central compartment volume Vd1(UA) (L/kg)")               # Hirai 2016 Results p. 1018 (Scott 1969 ref 14): Vd1(UA) = 0.25 L/kg
    lvd_ua_2     <- fixed(log(0.06));  label("Uric acid peripheral compartment volume Vd2(UA) (L/kg)")            # Hirai 2016 Results p. 1018 (Scott 1969 ref 14): Vd2(UA) = 0.06 L/kg
    lk12_ua      <- fixed(log(0.0068));label("Uric acid central-to-peripheral transfer rate k12(UA) (1/h)")        # Hirai 2016 Results p. 1018 (Scott 1969 ref 14): k12(UA) = 0.0068/h
    lk21_ua      <- fixed(log(0.031)); label("Uric acid peripheral-to-central transfer rate k21(UA) (1/h)")        # Hirai 2016 Results p. 1018 (Scott 1969 ref 14): k21(UA) = 0.031/h
    lkel_nr_ua   <- fixed(log(0.010)); label("Uric acid non-renal elimination rate kel_NR(UA) (1/h)")              # Hirai 2016 Results p. 1018 (Scott 1969 ref 14): kel_NR(UA) = 0.010/h (assumed constant irrespective of renal function)

    # ============================================================
    # Renal uric acid clearance as a function of creatinine clearance
    # Hirai 2016 Eq. 9 (adopted from Tykarski 1991 ref 17):
    #   CL_R(UA) [mL/min] = a * CRCL^b   with a = 1.23, b = 0.433
    # kel_R(UA) [1/h] = CL_R(UA) / (BW * Vd1(UA)) via Eq. 10 (with
    # a mL/min -> L/h conversion of 60/1000).
    # ============================================================
    a_cl_r_ua       <- fixed(1.23);   label("Renal UA clearance regression intercept a (mL/min)")                  # Hirai 2016 Eq. 9 (Tykarski 1991 ref 17): a = 1.23 mL/min
    e_crcl_cl_r_ua  <- fixed(0.433);  label("Renal UA clearance regression exponent b (dimensionless)")            # Hirai 2016 Eq. 9 (Tykarski 1991 ref 17): b = 0.433

    # ============================================================
    # Xanthine one-compartment endogenous kinetic model (Hirai 2016
    # Eqs. 11-12). Vd_xanthine is ASSUMED equal to the sum of Vd1(UA)
    # + Vd2(UA) = 0.31 L/kg (Methods p. 1015). Baseline xanthine
    # concentration 0.29 mg/L is the average reported by Mayer 2005
    # (ref 9). The elimination rate constant kel_xanthine is derived
    # inside model() from the steady-state balance
    # Xanthine_syn = Pool_xanthine(0) * kel_xanthine.
    # ============================================================
    lvd_xanthine   <- fixed(log(0.31));  label("Xanthine volume of distribution Vd_xanthine (L/kg)")               # Hirai 2016 Methods p. 1015: Vd_xanthine assumed = Vd(UA) = 0.31 L/kg (sum of Vd1 + Vd2)
    bl_xanthine    <- fixed(0.29);       label("Baseline serum xanthine concentration (mg/L)")                     # Hirai 2016 Methods p. 1015 (Mayer 2005 ref 9): baseline serum xanthine 0.29 mg/L

    # ============================================================
    # Baseline serum uric acid (mg/dL). Enters as the initial
    # condition for UApool1 and UApool2 via Hirai 2016 Eqs. 6-7 and
    # sets UAsyn via Eq. 8. Default 6.0 mg/dL is the normouricemic
    # reference used throughout the paper (Methods p. 1015 virtual
    # subject; Fig. 4 hyperuricemic-CKD scenarios use 9.0 mg/dL).
    # ============================================================
    bl_ua          <- fixed(6.0);        label("Baseline serum uric acid (mg/dL)")                                 # Hirai 2016 Methods p. 1015: normouricemic reference 6.0 mg/dL (target UA level); overridable per simulation via nlmixr2::ini() or rxode2::rxSolve inits argument
  })

  model({
    # ------------------------------------------------------------
    # Physical constants (spelled out in model() per convention).
    # ------------------------------------------------------------
    mw_febx <- 316.375   # febuxostat molecular weight (g/mol; C16H16N2O3S)
    mw_xan  <- 152.11    # xanthine molecular weight  (g/mol; C5H4N4O2)

    # ------------------------------------------------------------
    # Febuxostat PK parameters (per-subject; body-weight scaled).
    # ------------------------------------------------------------
    ka       <- exp(lka)                    # /h
    fdepot   <- exp(lfdepot)                # bioavailability (dimensionless)
    vc_L     <- exp(lvc) * WT               # L
    cl_Lh    <- exp(lcl) * WT               # L/h
    q_Lh     <- exp(lq)  * WT               # L/h
    vp_L     <- exp(lvp) * WT               # L
    kel      <- cl_Lh / vc_L                # /h (= k10 of the paper)
    k12      <- q_Lh  / vc_L                # /h
    k21      <- q_Lh  / vp_L                # /h

    # ------------------------------------------------------------
    # UA subsystem parameters (per-subject).
    # ------------------------------------------------------------
    vd_ua_1_L    <- exp(lvd_ua_1) * WT      # L (central UA compartment)
    vd_ua_2_L    <- exp(lvd_ua_2) * WT      # L (peripheral UA compartment)
    k12_ua       <- exp(lk12_ua)            # /h
    k21_ua       <- exp(lk21_ua)            # /h
    kel_nr_ua    <- exp(lkel_nr_ua)         # /h

    # Renal UA clearance from Cockcroft-Gault CLcr via Tykarski 1991.
    cl_r_ua_mlmin <- a_cl_r_ua * CRCL^e_crcl_cl_r_ua   # mL/min
    cl_r_ua_Lh    <- cl_r_ua_mlmin * 60 / 1000         # L/h
    kel_r_ua      <- cl_r_ua_Lh / vd_ua_1_L            # /h  (Hirai 2016 Eq. 10)

    # Initial UA amounts (mg) and endogenous UA synthesis rate (mg/h).
    # Hirai 2016 Eqs. 6-7: UApool_i(0) = 10 * bl_ua [mg/dL] * BW * Vd_i(UA) [L/kg]
    #                    = bl_ua [mg/dL] * 10 [dL/L] * Vd_i(UA) [L]
    # (The '10' converts mg/dL to mg/L before multiplying by volume in L.)
    uapool1_ic <- 10 * bl_ua * vd_ua_1_L    # mg
    uapool2_ic <- 10 * bl_ua * vd_ua_2_L    # mg
    uasyn      <- uapool1_ic * (kel_r_ua + kel_nr_ua)  # mg/h  (Hirai 2016 Eq. 8)

    # ------------------------------------------------------------
    # Xanthine subsystem parameters (per-subject).
    # ------------------------------------------------------------
    vd_xanthine_L <- exp(lvd_xanthine) * WT             # L
    xanthine_ic   <- bl_xanthine * vd_xanthine_L        # mg  (Hirai 2016 Eq. 12)
    xanthine_syn  <- uasyn                              # mg/h; Hirai 2016 Methods p. 1015: Xanthine_syn = UAsyn at steady state
    kel_xanthine  <- xanthine_syn / xanthine_ic         # /h  (steady-state constraint on Eq. 11)

    # ------------------------------------------------------------
    # Instantaneous concentrations feeding XO inhibition.
    # [I] = unbound plasma febuxostat (nM)
    # [S] = serum xanthine              (uM)
    # ------------------------------------------------------------
    fx_conc_mgL      <- central / vc_L                                     # mg/L total plasma febuxostat
    fx_unbound_nM    <- fu_febx * fx_conc_mgL * 1e6 / mw_febx              # (mg/L) / (g/mol) * 1e6 = nM
    xan_conc_mgL     <- xanthine / vd_xanthine_L                           # mg/L
    xan_conc_uM      <- xan_conc_mgL * 1000 / mw_xan                       # uM

    # XO inhibition ratio v'/v (Hirai 2016 Eqs. 2-4). Dimensionless.
    # v'/v = (Km + [S]) / (Km * (1 + [I]/Ki) + [S] * (1 + [I]/Ki_prime))
    inhib_ratio <- (km + xan_conc_uM) /
                   (km * (1 + fx_unbound_nM / ki) +
                    xan_conc_uM * (1 + fx_unbound_nM / ki_prime))

    # ------------------------------------------------------------
    # ODE system.
    # Febuxostat PK (Hirai 2016 Eq. 1):
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    f(depot)          <-  fdepot

    # Uric acid two-compartment turnover (Hirai 2016 Eq. 5):
    d/dt(uapool1) <-  uasyn * inhib_ratio +
                       uapool2 * k21_ua -
                       uapool1 * (kel_r_ua + kel_nr_ua + k12_ua)
    d/dt(uapool2) <-  uapool1 * k12_ua - uapool2 * k21_ua

    # Xanthine one-compartment turnover (Hirai 2016 Eq. 11).
    # Synthesis Xanthine_syn is constant; elimination is by XO (inhibited
    # by febuxostat via inhib_ratio, since xanthine is the XO substrate
    # upstream of uric acid).
    d/dt(xanthine) <- xanthine_syn - xanthine * kel_xanthine * inhib_ratio

    # Initial conditions (drug-free steady state; Hirai 2016 Eqs. 6-7, 12).
    uapool1(0)  <- uapool1_ic
    uapool2(0)  <- uapool2_ic
    xanthine(0) <- xanthine_ic

    # ------------------------------------------------------------
    # Observation-time algebraic outputs.
    # Cc      : total plasma febuxostat  (mg/L; = ng/mL / 1000)
    # sUA     : serum uric acid          (mg/dL) -- inverse of Eq. 6
    # sXan    : serum xanthine           (mg/L)
    # UAsynRate: instantaneous UA formation rate through XO (mg/h)
    # ------------------------------------------------------------
    Cc         <- fx_conc_mgL
    sUA        <- uapool1 / (10 * vd_ua_1_L)
    sXan       <- xan_conc_mgL
    UAsynRate  <- uasyn * inhib_ratio
  })
}
