test_that("invalid sample size",{
  expect_error(rBivarCop(n = 0, marg1 = 1, m1 = 1, c1 = 0.5, marg2 = 2,
                         m2 = 1.2, c2 = 0.65, cop = 2, tau = 0.5),
               "Please specify a valid input for sample size n.")
})

test_that("invalid sample size hier",{
  expect_error(dataGen(n = 19, gs = c(6, 4, 4, 4), marg1 = 1, marg2 = 1,
                       cop = 1, hyper_mu = c(1, 0.5, 0.25, 0.2, 0.1, 0.25),
                       hyper_prec = c(400, 400, 20, 800, 800, 20), seed = 1),
               "The total sample size for all groups in gs must sum to n.")
})
