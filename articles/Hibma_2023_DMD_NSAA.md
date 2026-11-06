# DMD North Star Ambulatory Assessment and loss of ambulation (Hibma 2023)

## Model and source

- Citation: Hibma JE, Jayachandran P, Neelakantan S, Harnisch LO.
  Disease progression modeling of the North Star Ambulatory Assessment
  for Duchenne Muscular Dystrophy. CPT Pharmacometrics Syst Pharmacol.
  2023;12(3):375-386. <doi:10.1002/psp4.12921>
- Description: Joint disease-progression model for the North Star
  Ambulatory Assessment (NSAA, total score 0-34) and loss of ambulation
  (LOA) in boys with Duchenne muscular dystrophy (DMD), fit to the
  placebo and domagrozumab arms of the phase II trial B5161002 pooled
  with the CINRG natural-history study B5161011 (Hibma 2023). The NSAA
  is a modified indirect-response state integrated over age from birth
  (NSAA = 0 at age 0): a zero-order production rate kin and a
  dissipation rate kout \* age that grows linearly with age, so the
  score rises to a peak near age 5.5 years and then declines. LOA is a
  parametric time-to-event model on age with a log-normal baseline
  hazard; the individual NSAA dissipation rate constant enters the
  hazard as a power covariate (kout / median kout)^theta. Correlated
  exponential IIV on kin and kout; additive residual error on the NSAA
  total score. No covariates were retained and the model carries no drug
  effect.
- Article: <https://doi.org/10.1002/psp4.12921> (open access,
  PMC10014057)

The article and its Supporting Information (Appendix S1, Figures S1-S7,
Tables S1-S6) were reviewed in full. No erratum or corrigendum was found
on the journal page or in PubMed as of 2026-10-26.

## Population

The model was fit to 2001 NSAA total scores from 276 boys with Duchenne
muscular dystrophy: 1404 observations from all 120 participants of the
Pfizer phase II domagrozumab trial B5161002 (placebo and active
sequences, ages 6 to \<16 years at randomization) and 597 observations
from 156 participants of the CINRG natural-history study B5161011
(NCT00468832). Median baseline age was 9.1 years (range 6.0-16.0) in
B5161002 and 7.4 years (4.0-20.3) in B5161011; median baseline NSAA
total score was 21.5 and 23.0 (Table S3). Twenty-nine participants lost
ambulation during follow-up (Tables S2 and S4). All participants were
male. B5161002 participants were assumed to be on stable corticosteroids
throughout (Appendix S1).

Data handling that shaped the fit: NSAA scores above the upper limit of
33 were treated with the M3 censoring method; 29 scores of zero (below
the lower limit of 1) were excluded from the NSAA likelihood and
recorded as loss-of-ambulation (LOA) events; a score of zero was taken
as irreversible, so later non-zero scores were set to zero (Methods,
Appendix S1).

The same information is available programmatically via
`rxode2::rxode(readModelDb("Hibma_2023_DMD_NSAA"))$population`.

## Source trace

| Element | Value in model | Source |
|----|----|----|
| NSAA ODE `d/dt(nsaa) = kin - kout * age * nsaa` | `d/dt(nsaa) <- kin - kout * t * nsaa` | Equation 1; Table 1 footnote a |
| Initial condition NSAA = 0 at age 0 | `nsaa(0) = 0` (default) | Discussion (“setting the initial condition to zero (age)”); Figure 1 |
| `kin` | 8.757 NSAA points/year | Table 1 |
| `kout` | 0.05767 1/year (multiplied by age) | Table 1 |
| IIV `kin` | 39.75% -\> omega^2 = 0.15801 | Table 1 |
| IIV `kout` | 92.82% -\> omega^2 = 0.86156 | Table 1 |
| Correlation `kin`-`kout` | 96.97% -\> covariance 0.35778 | Table 1 |
| Additive residual error | 2.782 points | Table 1; Results (“an additive error model”) |
| Log-normal baseline hazard on age | `h0 = phi-pdf(z) / (t * sigma * (1 - Phi(z)))`, `z = (log t - mu) / sigma` | Equation 8 |
| `mu` (location) | `ltmed_loa` = 3.216 | Table 1 |
| `sigma` (scale) | `sdlogt_loa` = 0.3551 | Table 1 |
| Hazard covariate `(kout / median kout)^theta` | `cov_kout_haz`, theta = `e_kout_haz` = 2.476 | Equations 9-10; Table 1 |
| Survival `S = exp(-cumhaz)` | `sur` | Equations 5-6 |

## Typical-value trajectory

``` r

mod <- rxode2::rxode(readModelDb("Hibma_2023_DMD_NSAA"))
#> ℹ parameter labels from comments will be replaced by 'label()'
mod_typ <- rxode2::zeroRe(mod)

ages <- seq(0, 20, by = 0.05)
ev_typ <- rxode2::et(ages, cmt = "nsaa")
typ <- as.data.frame(rxode2::rxSolve(mod_typ, ev_typ))
#> ℹ omega/sigma items treated as zero: 'etalkin', 'etalkout'

peak <- typ[which.max(typ$nsaa), ]
knitr::kable(
  data.frame(
    Quantity = c("Age at peak NSAA (years)", "Peak NSAA (points)"),
    Model = signif(c(peak$time, peak$nsaa), 3)
  )
)
```

| Quantity                 | Model |
|:-------------------------|------:|
| Age at peak NSAA (years) |  5.45 |
| Peak NSAA (points)       | 27.90 |

The typical subject peaks at 5.45 years. Hibma 2023 reports that the
peak NSAA total score “occurred on average at ~5.9 years (B5161011) and
5.5 years (B5161002)”; those are means over the individual fits, so the
typical-value peak is expected to sit close to, not exactly on, them.

### Dimensional / closed-form check

Equation 1 is linear in NSAA with a time-varying rate, so it has the
closed form
`NSAA(t) = kin * exp(-kout t^2 / 2) * integral_0^t exp(kout s^2 / 2) ds`.
Both sides use the same parameters, so the difference is pure ODE-solver
error and a tight bound is correct.

``` r

kin_typ <- exp(mod$theta[["lkin"]])
kout_typ <- exp(mod$theta[["lkout"]])
nsaa_closed <- function(a) {
  kin_typ * exp(-kout_typ * a^2 / 2) *
    integrate(function(s) exp(kout_typ * s^2 / 2), 0, a)$value
}
chk_ages <- c(2, 4, 6, 8, 10, 12, 15, 18)
chk_nsaa <- typ |>
  filter(round(time, 2) %in% chk_ages) |>
  mutate(closed_form = vapply(time, nsaa_closed, numeric(1)),
         rel_diff = (nsaa - closed_form) / closed_form)
knitr::kable(chk_nsaa |> select(time, nsaa, closed_form, rel_diff), digits = 5)
```

| time |     nsaa | closed_form | rel_diff |
|-----:|---------:|------------:|---------:|
|    2 | 16.22743 |    16.22743 |        0 |
|    4 | 26.00522 |    26.00522 |        0 |
|    6 | 27.66626 |    27.66626 |        0 |
|    8 | 24.13095 |    24.13095 |        0 |
|   10 | 19.23913 |    19.23913 |        0 |
|   12 | 15.14874 |    15.14874 |        0 |
|   15 | 11.23118 |    11.23118 |        0 |
|   18 |  8.99524 |     8.99524 |        0 |

``` r

stopifnot(nrow(chk_nsaa) == length(chk_ages), all(abs(chk_nsaa$rel_diff) < 1e-4))
```

## Virtual cohort

The model has no covariates, so a virtual cohort is defined only by its
random effects. Each subject is followed from birth to age 20 years.

``` r

rxode2::rxSetSeed(20231)
n_sub <- 200
ev_pop <- rxode2::et(seq(0, 20, by = 0.25), cmt = "nsaa")
sim <- as.data.frame(rxode2::rxSolve(mod, ev_pop, nSub = n_sub))
```

### Closed-form survival check

For a subject with fixed `kout`, Equations 5, 8, 9 and 10 integrate to
`S(t) = (1 - Phi(z))^COV`. The ODE-integrated `sur` must match it
subject by subject (same drawn `kout`, so the tolerance is solver error
only).

``` r

mu <- mod$theta[["ltmed_loa"]]
sigma <- mod$theta[["sdlogt_loa"]]
theta_haz <- mod$theta[["e_kout_haz"]]
chk_sur <- sim |>
  filter(time > 0) |>
  mutate(cov_kout = (kout / kout_typ)^theta_haz,
         sur_closed = (1 - pnorm((log(time) - mu) / sigma))^cov_kout,
         abs_diff = abs(sur - sur_closed))
stopifnot(max(chk_sur$abs_diff) < 1e-5)
max(chk_sur$abs_diff)
#> [1] 4.339213e-07
```

## Replicate Figure 4: NSAA versus age

Figure 4 of Hibma 2023 is the NSAA visual predictive check by study.
Below, the simulated median and 2.5th/97.5th percentiles of NSAA (with
additive residual error) are shown versus age for the virtual cohort
over the observed age range. The source VPC is conditioned on each
study’s observed design (ages, dropout and the exclusion of zero
scores), so the bands here are an unconditional analogue, not a
reproduction of its shaded areas.

``` r

vpc <- sim |>
  filter(time >= 4, time <= 18) |>
  group_by(time) |>
  summarise(
    p025 = quantile(sim, 0.025),
    p50 = median(sim),
    p975 = quantile(sim, 0.975),
    .groups = "drop"
  )
ggplot(vpc, aes(time)) +
  geom_ribbon(aes(ymin = p025, ymax = p975), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = p50), linewidth = 1) +
  geom_line(data = filter(typ, time >= 4, time <= 18), aes(y = nsaa),
            linetype = "dashed", colour = "red") +
  labs(x = "Age (years)", y = "NSAA total score",
       caption = "Solid: simulated median; dashed red: typical value; band: 2.5th-97.5th percentile") +
  theme_bw()
```

![Replicates Figure 4 of Hibma 2023 (unconditional analogue): simulated
NSAA median and 95% prediction interval versus
age.](Hibma_2023_DMD_NSAA_files/figure-html/fig4-1.png)

Replicates Figure 4 of Hibma 2023 (unconditional analogue): simulated
NSAA median and 95% prediction interval versus age.

Read from Figure 4, the observed and simulated medians fall from about
25 points near age 5 to about 19-20 at age 10 and about 12 at age 15
(B5161011 panel). The simulated median at those ages is:

``` r

fig4_chk <- vpc |>
  filter(time %in% c(5, 10, 15)) |>
  mutate(figure4_approx = c(25, 19.5, 12))
knitr::kable(fig4_chk |> select(time, p50, figure4_approx), digits = 1)
```

| time |  p50 | figure4_approx |
|-----:|-----:|---------------:|
|    5 | 24.1 |           25.0 |
|   10 | 19.8 |           19.5 |
|   15 | 11.2 |           12.0 |

``` r

stopifnot(all(abs(fig4_chk$p50 - fig4_chk$figure4_approx) < 4))
```

## Replicate Figure 3: probability of ambulation

The population probability of remaining ambulatory is the mean of the
individual survival curves. It is computed two ways: as the mean of
`sur` over the simulated cohort, and exactly by integrating the
closed-form individual survival over the log-normal `kout` random effect
(`etalkout ~ N(0, 0.86156)`), which needs no simulation.

``` r

omega_kout <- sqrt(mod$omega["etalkout", "etalkout"])
sur_pop_exact <- function(a) {
  integrate(function(e) {
    (1 - pnorm((log(a) - mu) / sigma))^exp(theta_haz * e) * dnorm(e, 0, omega_kout)
  }, -Inf, Inf)$value
}
km <- sim |>
  group_by(time) |>
  summarise(sur_sim = mean(sur), .groups = "drop") |>
  filter(time > 0) |>
  mutate(sur_exact = vapply(time, sur_pop_exact, numeric(1)))
ggplot(km, aes(time)) +
  geom_line(aes(y = sur_exact), linewidth = 1) +
  geom_line(aes(y = sur_sim), colour = "steelblue", linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Age (years)", y = "Probability of ambulation",
       caption = "Solid: exact population mean; dashed: mean over the simulated cohort") +
  theme_bw()
```

![Replicates Figure 3 of Hibma 2023: population probability of
ambulation versus
age.](Hibma_2023_DMD_NSAA_files/figure-html/fig3-1.png)

Replicates Figure 3 of Hibma 2023: population probability of ambulation
versus age.

Figure 3’s simulated band reads at roughly 0.95 at age 12, 0.75-0.80 at
age 15 and 0.55-0.62 at age 19.

``` r

fig3_chk <- km |>
  filter(time %in% c(12, 15, 19)) |>
  mutate(figure3_approx = c(0.95, 0.78, 0.58))
knitr::kable(fig3_chk, digits = 3)
```

| time | sur_sim | sur_exact | figure3_approx |
|-----:|--------:|----------:|---------------:|
|   12 |   0.883 |     0.897 |           0.95 |
|   15 |   0.763 |     0.776 |           0.78 |
|   19 |   0.621 |     0.628 |           0.58 |

``` r

stopifnot(
  all(abs(fig3_chk$sur_exact - fig3_chk$figure3_approx) < 0.08),
  # The cohort mean estimates the exact curve; 200 subjects give a
  # standard error of at most ~0.035, so this envelope tolerates any draw.
  all(abs(fig3_chk$sur_sim - fig3_chk$sur_exact) < 0.12)
)
```

## Replicate Figure 5: one-year NSAA trajectories by baseline age

Figure 5 overlays one year of observed NSAA on simulated trajectories
for boys entering at ages 4, 5, 6 and 7 with a baseline NSAA between 16
and 30. The typical-value change over one year from each entry age is:

``` r

fig5 <- typ |>
  filter(round(time, 2) %in% c(4, 5, 6, 7, 8)) |>
  select(age = time, nsaa)
fig5 |>
  mutate(change_next_year = lead(nsaa) - nsaa) |>
  filter(age < 8) |>
  knitr::kable(digits = 2)
```

| age |  nsaa | change_next_year |
|----:|------:|-----------------:|
|   4 | 26.01 |             1.73 |
|   5 | 27.74 |            -0.07 |
|   6 | 27.67 |            -1.38 |
|   7 | 26.29 |            -2.15 |

The typical trajectory is near its peak between ages 4 and 8 and changes
by about two points or less per year, consistent with the nearly flat
median lines of Figure 5 and with the paper’s conclusion that a one-year
trial in this age window needs a treatment that restores function to
reach 80% power (Figure 6).

## Assumptions and deviations

- **IIV scale.** Table 1 prints IIV as a percentage. Its 95% CIs are
  exactly the square roots of a Wald interval on omega^2 (for `kin`:
  sqrt(0.15801 -/+ 1.96 x 0.02034 x 0.15801) = 0.3895, 0.4053, printed
  38.95-40.53; for `kout`: 0.9192-0.9371, printed 91.92-93.72), so the
  percentage is 100 x sqrt(omega^2) and omega^2 = (IIV/100)^2. The
  log(1 + CV^2) reading gives 0.147 and 0.623 instead and does not
  reproduce the printed CIs. The printed correlation CI (97.51, 96.51)
  has its limits in reverse order; the point estimate 96.97% is used.
- **Median kout in the hazard.** Equation 10 normalizes `kout` by
  “median(kout)” without printing the value. With log-normal IIV the
  population median of `kout` equals its typical value, so the model
  uses `exp(lkout)` and the ratio reduces to `exp(etalkout)`. If the
  authors used the median of the empirical Bayes estimates, the
  difference is absorbed by `mu` and small (kout eta shrinkage 8.8%).
- **Hazard time scale.** The hazard is a function of age (the model’s
  independent variable), and LOA is evaluated from birth, consistent
  with the age axis of Figures 3 and S5. The log-normal hazard is set to
  its limit of zero at age 0 to keep `log(t)` finite.
- **Residual error.** Table 1 reports the residual error as a theta in
  NSAA units, i.e. the standard deviation of an additive error (the
  Results name an additive error model).
- **Score bounds.** The model does not truncate NSAA at 0 or 34. The
  authors let predictions exceed 34 and capped them at 34 only in
  post-processing (Discussion); the NSAA likelihood excluded zero
  scores, and LOA is described by the hazard rather than by NSAA
  reaching zero.
- **M3 censoring and the LOA likelihood** are estimation details (data
  flags), not model structure; the packaged model declares only the NSAA
  endpoint and exposes `hazard`, `cumhaz` and `sur` for simulation.
- **No drug effect.** The B5161002 treatment sequences were tested as
  covariates and not retained; the model describes natural progression.
- **Table S1 columns.** Some Table S1 rows appear to have the two study
  columns transposed (for example, 120 participants with missing
  mutation status are listed under B5161011, but mutation status was not
  recorded for B5161002, which has 120 participants); the population
  metadata uses the per-study counts from the Results and Table S3.
- **Base models.** The separately fit base NSAA model (Table S5) and
  base LOA model (Table S6) are intermediate steps of the joint model
  and are not packaged.
