context("ESCO simulations")

test.params <- newescoParams(nGenes = 50, nCells = 20)

test_that("escoSimulate takes input correctly", {
  sim<-escoSimulate(test.params, type = "single", nCells = 100)
  params = metadata(sim)$Params
  newnCells = getParam(params, "nCells")
  expect_equal(newnCells,100)
})

test_that("escoSimulate output is valid", {
  expect_true(validObject(escoSimulate(test.params)))
  expect_true(validObject(escoSimulate(test.params, type = "group", group.prob = c(0.5, 0.5))))
})

test_that("escoSimulate warns correctly", {
  expect_warning(escoSimulate(test.params, type = "group",
                               group.prob = c(1)),
                 "nGroups is 1, switching to single mode")
  expect_warning(escoSimulate(test.params, trials=2),
                 "Detect calls for multiple trials, but no directory to save files...")
  expect_warning(escoSimulate(test.params, dropout.type = c('nope')),
                 "Detect corrupted parameter dropout.type: should be either 'zeroinfalte' or 'downsample', or '' for no dropout...")
  expect_warning(escoSimulate(test.params, dropout.mid = c(1,2)),
                 "Detect calls for simulating multiple configurations, but no directory to save files...")
  expect_warning(escoSimulate(test.params, alpha_mean = c(0.1, 0.2)),
                 "Detect calls for simulating multiple configurations, but no directory to save files...")
})


