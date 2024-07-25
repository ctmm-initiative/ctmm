testthat::test_that("Langevin Cilla DIM=1 works as expected",{
    data('buffalo')
    load(testthat::test_path('testdata', 'test_Langevin.rda'))

    t <- buffalo$Cilla$t/3600
    dt <- c(Inf, diff(t))
    dti <- sort(dt,index.return=TRUE)
    dtl <- unique(dti$x) # dt levels
    dti <- dti$ix # sort indices


    for(idx in 1:length(CTMMs)){
        testthat::expect_identical(
            Langevin(
                t = t,
                CTMM = c(
                    CTMMs[[idx]],
                    list(
                        dt = dt,
                        dti = dti,
                        dtl = dtl
                    )
                ),
                DIM = 1
             ),
            expected[[idx]]
        )
    }
})