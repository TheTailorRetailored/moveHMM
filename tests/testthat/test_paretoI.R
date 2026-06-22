context("Pareto Type I step lengths")

test_that("density and distribution use the documented parameter order", {
    xmin <- 2
    mu <- 3

    expect_equal(dparetoI(xmin, xmin, mu), mu / xmin)
    expect_equal(dparetoI(4, xmin, mu), mu / xmin * (4 / xmin)^(-mu - 1))
    expect_equal(dparetoI(1, xmin, mu), 0)
    expect_equal(pparetoI(xmin, xmin, mu), 0)
    expect_equal(pparetoI(4, xmin, mu), 1 - (xmin / 4)^mu)
    expect_equal(dparetoI_rcpp(c(xmin, 4), xmin, mu), dparetoI(c(xmin, 4), xmin, mu))
})

test_that("random Pareto steps respect xmin", {
    set.seed(42)
    steps <- rparetoI(1000, xmin = 2, mu = 3)

    expect_length(steps, 1000)
    expect_true(all(steps >= 2))
})

test_that("simData can generate Pareto step lengths", {
    set.seed(42)
    data <- simData(
        nbAnimals = 1,
        nbStates = 1,
        stepDist = "paretoI",
        angleDist = "none",
        stepPar = c(2, 3),
        obsPerAnimal = 25
    )

    expect_s3_class(data, "moveData")
    expect_true(all(stats::na.omit(data$step) >= 2))
})
