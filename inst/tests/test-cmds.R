context("Check cmds")

test_that("verbose works", {
  a <- alply(replicate(5,c(0,1,1,1,0,1,1,1,0)), 2, matrix, nrow = 3, ncol = 3)

  expect_that(cmds(a,v=TRUE), prints_text("Iteration"))
  expect_that(cmds(QuadCurves,k=1,l=3,v=T), prints_text("Init.+1.87234"))
})
  
test_that("quadratic example works", {  
  expect_that(cmds(QuadCurves,k=1,l=3), prints_text("The algorithm converged.+0"))
  res <- cmds(QuadCurves,k=1,l=0)
  DL.X <- Dist.List(res$XL,res$params)
  expect_that(DL.X[-1], is_equivalent_to(QuadCurves[-1]))
})

test_that("expanding triangle example works", {
  res <- cmds(ExpandingTriangle,k=2,l=0)
  DL.X <- Dist.List(res$XL,res$params)
  expect_that(DL.X, is_equivalent_to(ExpandingTriangle))
})

test_that("expanding circle example works", {
  res <- cmds(ExpandingCircle,k=2,l=0)
  DL.X <- Dist.List(res$XL,res$params)
  expect_that(DL.X, is_equivalent_to(ExpandingCircle))
})

test_that("T=1 works", {
  DL <- list(a.m(dist(seq(0,1,l=5))))
  expect_that(cmds(DL), prints_text("Total cost C:.+0"))
})
              
