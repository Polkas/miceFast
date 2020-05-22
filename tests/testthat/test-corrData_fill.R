

context("corrData-fill")


test_that("fill-data",
          {
            set.seed(1234)

            #parameters

            power = 4 # power of 10 - number of observations - should be adjusted to a computer capabilities

            nr_var = 7 #CHANGE - only if you generate a bigger corr matrix:  number of variables - independent and one dependent

            grs = max(c(10**(power-3),10)) # grouping variable - number of groups

            iters = 10 # number of iterations for benchmarking

            ## generete example - data

            ##positive-defined correlation matrix

            cors = matrix(c(1,0.6,0.7,0.4,0.4,0.5,0.35,
                            NA,1,0.2,0.05,0.1,0.12,0.15,
                            NA,NA,1,0.15,0.15,0.1,0.08,
                            NA,NA,NA,1,0.12,0.15,0.1,
                            NA,NA,NA,NA,1,0.15,0.2,
                            NA,NA,NA,NA,NA,1,0.15,
                            NA,NA,NA,NA,NA,NA,1),7,7,byrow = T)

            cors[lower.tri(cors)] = t(cors)[lower.tri(cors)]

            # automatic corr matrix - close to diagonal

            #cors = stats::rWishart(100,nr_var+10,diag(nr_var))

            #cors = apply(cors,1:2,mean)/(nr_var+10)

            #cors

            ##

            model = new(corrData,10,10^power,rep(0,nr_var),cors)

            data_bin = model$fill("binom")
            data_disc = model$fill("discrete")
            data_con = model$fill("contin")

            rm(model)

            cor_b=cor(data_bin)
            cor_d=cor(data_disc)
            cor_c=cor(data_con)

            cor_cor_a = cor(cbind(as.vector(cors),as.vector(cor_c)))[1,2] > 0.99
            cor_cor_b = cor(cbind(as.vector(cors),as.vector(cor_c)))[1,2] > 0.99
            cor_cor_c = cor(cbind(as.vector(cors),as.vector(cor_c)))[1,2] > 0.99

            expect_true(all(cor_cor_a,cor_cor_b,cor_cor_c))

          })

