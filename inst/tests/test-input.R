context("Input checking for cmds")

test_that("incorrect D input raises errors", {
  a <- matrix(0,2,2)
  expect_that(cmds(a),throws_error("'D' is not a list."))
  
  a <- llply(rep(0, 3), matrix, nrow = 3, ncol = 2) 
  expect_that(cmds(a),throws_error("All elements of 'D' must be square matrices."))
 
  a <- alply(replicate(3, c(0,1,2,0)), 2, matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a),throws_error("All elements of 'D' must be symmetric matrices."))
  
  a <- alply(replicate(3, c(1,1,1,0)), 2, matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a),throws_error("All elements of 'D' must have zero diagonal."))

  a <- llply(rep(0, 3), matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a),throws_error("All rows of 'D' are equal to zero or NA."))
  
  a <- alply(replicate(3, c(0,NA,NA,0)), 2, matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a),gives_warning("Some rows of 'D' are equal to zero or NA."))

  a <- alply(replicate(3, c(0,Inf,Inf,0)), 2, matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a), throws_error("'D' contains infinite values."))

  a <- alply(replicate(3, c(0,-1,-1,0)), 2, matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a), throws_error("'D' contains negative values."))
  
  a <- list(matrix(0,2,2),matrix(0,3,3),matrix(0,2,2))
  expect_that(cmds(a), throws_error("All elements of 'D' must have the same dimensions."))
  
  a <- alply(replicate(3, c(0,1,NA,1,0,0,NA,0,0)), 2, matrix, nrow = 3, ncol = 3)
  #expect_that(cmds(a), gives_warning("There are some missing values in 'D'."))
  #expect_that(cmds(a)$W[[1]][1,3], equals(0))
  
})

test_that("Incorrect k or l raises errors", {
  a <- alply(replicate(3, c(0,1,1,1,0,1,1,1,0)), 2, matrix, nrow = 3, ncol = 3)
  expect_that(cmds(a,k=0.2), throws_error("'k' must be in 1,2, ..., N-1."))
  expect_that(cmds(a,k=3), throws_error("'k' must be in 1,2, ..., N-1."))

  expect_that(cmds(a,k=1,l=-0.5), throws_error("l must be positive."))
})

test_that("Incorrect W raises errors", {
  a <- alply(replicate(3, c(0,1,1,1,0,1,1,1,0)), 2, matrix, nrow = 3, ncol = 3)
  W <- llply(rep(1, 3), matrix, nrow = 2, ncol = 2)
  expect_that(cmds(a,W=W), throws_error("Your custom weights should be of the same dimensions as D."))
  W <- matrix(1,2,2)
  expect_that(cmds(a,W=W), throws_error("W must be either 'classic', 'kamada-kawai' or a custom list of weights."))
})

test_that("Incorrect optional parameters raise errors", {
  a <- alply(replicate(3, c(0,1,1,1,0,1,1,1,0)), 2, matrix, nrow = 3, ncol = 3)
  expect_that(cmds(a,M=matrix(1,2,2)), throws_error("The custom penalty matrix must be of size TxT."))

  expect_that(cmds(a,init="custom"), throws_error("'init' must be 'average' or 'random' or a custom list."))
  expect_that(cmds(a,init=matrix(0,3,1)), throws_error("'init' must be 'average' or 'random' or a custom list."))
  expect_that(cmds(a,init=replicate(3,list(matrix(0,4,1)))), throws_error("Your custom initialization should be a list with matrices of size Nxk."))
})

test_that("Init method 'average' yields correct format of XL", {
  a <- alply(replicate(5,c(0,1,1,1,0,1,1,1,0)), 2, matrix, nrow = 3, ncol = 3)

  expect_that(length(cmds(a,k=1)$XL), equals(5))
  expect_that(dim(cmds(a,k=1)$XL[[1]]), equals(c(1,3)))

  expect_that(length(cmds(a,k=2)$XL), equals(5))
  expect_that(dim(cmds(a,k=2)$XL[[1]]), equals(c(2,3)))
})

test_that("Init method 'random' yields correct format of XL", {
  a <- alply(replicate(5,c(0,1,1,1,0,1,1,1,0)), 2, matrix, nrow = 3, ncol = 3)

  expect_that(length(cmds(a,k=1,l=0,init="random")$XL), equals(5))
  expect_that(dim(cmds(a,k=1,l=0,init="random")$XL[[1]]), equals(c(1,3)))

  expect_that(length(cmds(a,k=2,l=0,init="random")$XL), equals(5))
  expect_that(dim(cmds(a,k=2,l=0,init="random")$XL[[1]]), equals(c(2,3)))
})
