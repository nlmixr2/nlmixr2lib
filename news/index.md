# Changelog

## Version 0.3.2

- Add Kim 2006 model for IgG metabolism
- Add Xie 2019 agomelatine PK model
- Drop `qs` since it will be archived
- Update tumor growth inhibition models

## Version 0.3.1

CRAN release: 2025-10-22

- Bug fix for replacement of multiplicative expressions in `nlmixr2lib`
- phenylalanine_charbonneau_2021 had its net protein breakdown parameter
  corrected
- Kyhl_2016_nalmefene model was added

## Version 0.3.0

CRAN release: 2024-10-07

- Added ability to choose style type when modifying models. Currently
  supported styles are: “camel” for `variablesLikeThis`, “snake” for
  `variables_like_this`, “dot” for `variables.like.this` and “blank” for
  `variableslikethis`. This can be selected with
  [`setCombineType()`](https://nlmixr2.github.io/nlmixr2lib/reference/setCombineType.md).

- With the new combination style, you can change how `eta` variables are
  constructed with the `option(nlmixr2lib.etaCombineType="camel")` or
  whatever you wish it to the variable style to be.

- Added new model building framework for building models

  - **PK model building functions**

    - [`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md)/[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md)
      which were present before, but now modified and made a bit more
      robust, more closely matching literature method of transit
      compartments.

    - [`addDepot()`](https://nlmixr2.github.io/nlmixr2lib/reference/addDepot.md)/[`removeDepot()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeDepot.md)
      which were present before, but modified to be a bit more robust.

    - [`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md)
      which adds a Weibull absorption to a PK model

    - [`convertMM()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertMM.md)
      converts linear elimination to Michaelis-Menten elimination

    - `transPK()` converts the `cl` style parameter transformations to
      various other PK transformations like `k`, `aob`, `alpha`, `k12`

  - **PD model building functions**

  - [`addIndirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addIndirectLin.md)
    – this adds an indirect effect model to a PK model that has a
    concentration `Cc` in the model. This purposely uses a simple linear
    effect of `Cc*Ek` or `Cc*Ik` so it will be easy to parse and turn
    into other functional forms (like `Emax` or `Hill`). If the PK model
    is not present it will use `Cc` as a covariate in a purely PD
    models.

  - [`addIndirect()`](https://nlmixr2.github.io/nlmixr2lib/reference/addIndirect.md)
    – this builds on
    [`addIndirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addIndirectLin.md)
    and adds `Emax` or `Hill` models to a PK model. You can also set
    `imax=1` or `emax=1` to drop these parameters from being estimated
    in the model. Additionally `hill=TRUE` will add a Hill coefficient
    to the sigmoid model.

  - [`addEffectCmtLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addEffectCmtLin.md)
    – this adds an effect compartment based on the `Cc` in the model.
    The linear effect can be modified into other function forms.

  - [`addDirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addDirectLin.md)
    – this adds a direct effect model based on the `Cc` in the model.

  - **Changing functional forms of Effect models**

    - [`convertEmax()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertEmax.md)
      changes linear effect models to Emax models

    - [`convertEmaxHill()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertEmaxHill.md)
      changes linear effect models to Hill models

    - [`convertQuad()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertQuad.md)
      changes linear effect models to quadratic models

    - [`convertLogLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertLogLin.md)
      changes linear effect models to log-linear models

  - **Changing functional forms of Baselines in non-indirect response
    models**

    - [`addBaselineConst()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineConst.md)
      changes the zero baseline to a estimated constant

    - [`addBaselineLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineLin.md)
      changes the zero baseline to a estimated constant and a linear
      constant with respect to `time`.

    - [`addBaselineExp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineExp.md)
      changes the zero baseline to a exponential decay with respect to
      time

    - [`addBaseline1exp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaseline1exp.md)
      – the baseline effect is changed from zero to to an exponential
      approaching to a constant (with respect to time).

  - **Changing model properties** (all use
    [`addCmtProp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addCmtProp.md))

    - [`addBioavailability()`](https://nlmixr2.github.io/nlmixr2lib/reference/addCmtProp.md)
      adds bioavailability property to a compartment

    - [`addRate()`](https://nlmixr2.github.io/nlmixr2lib/reference/addCmtProp.md)
      adds a modeled rate to a compartment

    - [`addDur()`](https://nlmixr2.github.io/nlmixr2lib/reference/addCmtProp.md)
      adds modeled duration to a compartment

    - [`addIni()`](https://nlmixr2.github.io/nlmixr2lib/reference/addCmtProp.md)
      adds an initial value to a compartment

    - [`addLag()`](https://nlmixr2.github.io/nlmixr2lib/reference/addCmtProp.md)
      adds a lag time to the a compartment

- Add Carlsson Petri (2021) liraglutide PK model

- Add Cirincione (2017) exenatide immediate-release PK model

- Add a variety of indirect response models

- Add a variety of tumor growth inhibition models and move all oncology
  models into a new model database directory

- Add a variety of double-absorption PK models

- `cp` and related `cpddSd` and `cppropSd` were renamed to `Cc`,
  `CcAddSd` and `CcPropSd` (fix
  [\#70](https://github.com/nlmixr2/nlmixr2lib/issues/70)).

- Multiple-endpoint models will have the `DV` column in the modeldb
  separated by commas.

## Version 0.2.0

CRAN release: 2023-03-29

- Work with the new `rxode2` version 2.0.12
  [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html) and
  [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html)
  assignment methods.
- Therapeutic-area specific models have begun being added.
- Models can now give the user some additional information load via the
  `message` meta-data.
- Models can now be in different directories. The change is for ease of
  maintaining the library, it is not a change that affects users.
- A regression where
  [`addEta()`](https://nlmixr2.github.io/nlmixr2lib/reference/addEta.md)
  did not change the parameter, related to a change in `rxode2`, was
  fixed.
- [`addEta()`](https://nlmixr2.github.io/nlmixr2lib/reference/addEta.md)
  detects where to add etas more robustly when covariates are on the
  parameter.

### Models added

- Add Davda (2014) mAb consensus model
- Add Liu (2017) time-dependent clearance model based on nivolumab
- Add Kovalenko (2020) dupilumab PK model
- Add Soehoel (2022) tralokinumab PK model
- Add Zhu (2017) lebrikizumab PK model

## Version 0.1.0

CRAN release: 2022-10-31

- Initial version
