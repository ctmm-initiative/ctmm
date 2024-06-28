test_that("check dexp2", {
  input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -10, -100)
  expected <- c(
    1,
    0.9999999979388464,
    0.8646647167633873,
    0.18126924692201812,
    0.0198013266932447,
    0,
    -0.02020134002675581,
    -0.22140275816016985,
    -485165194.40979034,
    -7.22597376812575e+86
  )
  for(idx in 1:length(input)){
    expect_identical(dexp2(input[idx]), expected[idx], label = paste0('dexp2 of ', input[idx]))
  }
})
test_that("check dexp1", {
  input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -10, -100)
  expected <- c(
    1,
    0.9999546000702375,
    0.6321205588285577,
    0.09516258196404043,
    0.009950166250831947,
    0,
    -0.010050167084168057,
    -0.10517091807564763,
    -22025.465794806714,
    -26881171418161356094253400435962903554686976
  )
  for(idx in 1:length(input)){
    expect_identical(dexp1(input[idx]), expected[idx], label = paste0('dexp1 of ', input[idx]))
  }
})