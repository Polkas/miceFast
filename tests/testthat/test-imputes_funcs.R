
context("miceFast-imputes_funcs")


test_that("imputes",
          {

            set.seed(1234)

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

            n_vars = ncol(cors)

            posit_y = 1
            posit_x = 2:(n_vars-2)
            posit_w = n_vars-1
            posit_grs = n_vars
            posit_NA = n_vars+1

            ## NA index

            index_NA = 1:nrow(data_con) %in% sample(1:nrow(data_con),10^(power-1))

            fill_by_NA = function(v,index_NA){

              v[index_NA] = NA

              v
            }

            ######################

            group_d = floor(pnorm(data_disc[,posit_grs])*grs)
            group_c = floor(pnorm(data_con[,posit_grs])*grs)
            group_b = floor(pnorm(data_bin[,posit_grs])*grs)

            w_d = abs(data_disc[,posit_w])
            w_c = abs(data_con[,posit_w])
            w_b = abs(data_bin[,posit_w])

            data_disc_NA = cbind(fill_by_NA(data_disc[,posit_y],index_NA),data_disc[,posit_x],w_d,group_d,index_NA)
            data_con_NA = cbind(fill_by_NA(data_con[,posit_y],index_NA),data_con[,posit_x],w_c,group_c,index_NA)
            data_bin_NA = cbind(fill_by_NA(data_bin[,posit_y],index_NA),data_bin[,posit_x],w_b,group_b,index_NA)

            colnames(data_bin_NA) = c("y",paste0("x",posit_x),"weights","group","index_NA")
            colnames(data_disc_NA) = c("y",paste0("x",posit_x),"weights","group","index_NA")
            colnames(data_con_NA) = c("y",paste0("x",posit_x),"weights","group","index_NA")

            pred_lm_bayes = fill_NA(x=as.matrix(data_con_NA),
                                   model="lm_bayes",
                                   posit_y=1,
                                   posit_x=c(2,3,4))

            test0 = cor(cbind(pred_lm_bayes[index_NA],data_con[index_NA,1]))[1,2] > 0.5

            pred_lm_noise_N = fill_NA_N(x=as.matrix(data_con_NA),
                                       model="lm_noise",
                                       posit_y=1,
                                       posit_x=c(2,3,4),
                                       times=10)

            test1 = cor(cbind(pred_lm_noise_N[index_NA],data_con[index_NA,1]))[1,2] > 0.5

            pred_lda = fill_NA(x=as.matrix(data_bin_NA),
                              model="lda",
                              posit_y=1,
                              posit_x=c(2,3,4),w=data_bin_NA[,6])

            test2 = cor(cbind(pred_lda[index_NA],data_bin[index_NA,1]))[1,2] > 0.5

            vifs = VIF(x=as.matrix(data_con_NA)*5,
                        posit_y=1,
                        posit_x=c(2,3,4))

            collin_x = function(y,x){

              full_rows = !is.na(y)&complete.cases(x)

              XXinv=solve(crossprod(scale(x[full_rows,],scale=F)))

              vifs = sapply(1:ncol(XXinv),function(i) det(XXinv[i,i,drop=FALSE])*det(XXinv[-i,-i,drop=FALSE])/det(XXinv))

              vifs2 = sapply(1:ncol(XXinv),function(i) Reduce("*",svd(XXinv[i,i,drop=FALSE])$d)*Reduce("*",svd(XXinv[-i,-i,drop=FALSE])$d)/Reduce("*",svd(XXinv)$d))

              list(vifs,vifs2)

            }

          test3 = cor(cbind(vifs,collin_x(data_con_NA[,1]*5,data_con_NA[,c(2,3,4)]*5)[[1]]))[1,1] > 0.9


            expect_true(all(c(test0,test1,test2,test3)))

          })
