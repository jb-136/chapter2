test_that("check_continuous works", {
  funky <- data.frame(
    runif(3),
    round(runif(3)),
    as.character(runif(3)),
    as.character(round(runif(3))),
    c("a", "b", "3"),
    c(NA, 1.4, 2.5),
    c(NA, "a", "b"),
    c(NA,4,5),
    stringsAsFactors=FALSE)
  funky[,1] <- as.numeric(funky[,1])
  funky[,2] <- as.numeric(funky[,2])
  funky[,6] <- as.numeric(funky[,6])
  colnames(funky) <- NULL
  expect_true(all(check_continuous(funky)[c(1,3,6)]))

})

test_that("chapter2_fitGeiger works", {
  data(geospiza, package="geiger")
  geospiza$dat[4,3] <- NA
  chapter2 <- match_data(geospiza$phy, geospiza$dat)
  fc_result <- chapter2_fitGeiger(chapter2, keep="continuous")
  expect_true(inherits(fc_result[[1]][[1]], "gfit"))
})
