test_that("check_continuous works", {
  expect_true(check_continuous(runif(5)))
  expect_true(check_continuous(c(NA,runif(5))))
  expect_false(check_continuous(round(c(NA,runif(5)))))
})
