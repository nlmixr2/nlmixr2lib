# Convert double absorption to sequential delays

Reparameterizes a double-absorption model so the second path's delay is
expressed relative to the first path's zero-order duration (\`Tlag2 =
Tk01\`), matching Monolix's "sequential" double-absorption models. Only
valid when the first path absorbs zero-order; raises an error otherwise.

## Usage

``` r
convertAbsSequential(
  ui,
  central = "central",
  depot = "depot",
  depot2 = "depot2"
)
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- central:

  central compartment name

- depot:

  depot compartment name

- depot2:

  name of the second depot compartment

## Value

a model with the second delay tied to the first duration

## See also

Other absorption:
[`addSecondAbsorption()`](https://nlmixr2.github.io/nlmixr2lib/reference/addSecondAbsorption.md),
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md),
[`addZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addZeroOrderAbs.md),
[`convertAbsForceLongerDelay()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsForceLongerDelay.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md),
[`removeZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeZeroOrderAbs.md)

## Author

Matthew L. Fidler
