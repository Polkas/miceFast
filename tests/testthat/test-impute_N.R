
context("miceFast-impute_N")


test_that("impute_N",
          {
            set.seed(1234)

            data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality))
            weights = rgamma(nrow(data),3,3) # a numeric vector - positive values
            #groups = as.numeric(airquality[,5]) # a numeric vector not integers - positive values
            groups = as.numeric(sample(1:3,nrow(data),replace=T)) # a numeric vector not integers - positive values

            model = new(miceFast)
            model$set_data(data) # providing by a reference
            model$set_w(weights) # providing by a reference
            model$set_g(groups)  # providing by a reference
            #impute adapt to provided parmaters like w or g
            #Warning - if data is not sorted increasingly by the g then it would be done automatically
            #during a first imputation
            #Simple mean - permanent imputation at the object and data

            imp_miceFast = model$impute_N("lm_bayes",2,c(3,4,5,6),10000)$imputations
            imp_miceFast_pred = model$impute("lm_pred",2,c(3,4,5,6))$imputations

            rm(model)

            expect_true(mean((imp_miceFast-imp_miceFast_pred)**2)<1)

          })
