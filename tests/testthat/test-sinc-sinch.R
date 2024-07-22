test_that("check sinc", {
  input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -1, -10, -100)
  expected <- c(
    -0.005063656411097588,
    -0.05440211108893698,
    0.8414709848078965,
    0.9983341664682815,
    0.9999833334166665,
    1,
    0.9999833334166665,
    0.9983341664682815,
    0.8414709848078965,
    -0.05440211108893698,
    -0.005063656411097588
  )
  for(idx in 1:length(input)){
    expect_identical(sinc(input[idx]), expected[idx], label = paste0('sinc of ', input[idx]))
  }
  for(idx in 1:length(input)){
    expect_identical(sinc(input[idx], sin(input[idx])), expected[idx], label = paste0('sinc of ', input[idx]))
  }

  expect_identical(sinc(input), expected, label = 'sinc works for vectors')
  expect_identical(sinc(input, sin(input)), expected, label = 'sinc works for vectors (precalculated)')

  
})

test_that("check sinch", {
  input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -1, -10, -100)
  expected <- c(
    134405857090806786660967198606715892269056,
    1101.3232874703394,
    1.1752011936438014,
    1.0016675001984403,
    1.0000166667500003,
    1,
    1.0000166667500003,
    1.0016675001984403,
    1.1752011936438014,
    1101.3232874703394,
    134405857090806786660967198606715892269056
  )
  for(idx in 1:length(input)){
    expect_identical(sinch(input[idx]), expected[idx], label = paste0('sinch of ', input[idx]))
  }
  for(idx in 1:length(input)){
    expect_identical(sinch(input[idx], sinh(input[idx])), expected[idx], label = paste0('sinch of ', input[idx]))
  }

  expect_identical(sinch(input), expected, label = 'sinch works for vectors')
  expect_identical(sinch(input, sinh(input)), expected, label = 'sinch works for vectors (precalculated)')
})