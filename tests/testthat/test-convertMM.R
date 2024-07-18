test_that(".replaceMultC works", {

  e1 <- str2lang("d/dt(central) <- ka * depot - kel * central")

  e2 <- str2lang("d/dt(central) <- ka * depot - (vm * central/vc)/(km + central/vc)")

  expect_equal(.replaceMultC(e1, str2lang("kel"),
                             str2lang("central"),
                             str2lang("(vm*central/vc)/(km+central/vc)")),
               e2)

  expect_equal(.replaceMultC(e1,
                             str2lang("central"),
                             str2lang("kel"),
                             str2lang("(vm*central/vc)/(km+central/vc)")),
               e2)

  expect_equal(.replaceMultC(e1,
                             str2lang("funny"),
                             str2lang("kel"),
                             str2lang("(vm*central/vc)/(km+central/vc)")),
               e1)

})
