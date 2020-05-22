
context("miceFast-impute")


test_that("impute",
          {
            #library(dplyr)
            #library(mice)
            #library(broom)

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


            ######################Discrete

            #mice.impute.lda = mice.impute.lda(data_disc[,posit_y],!index_NA,data_disc[,posit_x])

            model = new(miceFast)
            data = data_disc_NA[,c(posit_y,posit_x)]
            model$set_data(data)
            pred_miceFast =  model$impute("lda",posit_y,posit_x)

            rm(model)

            a = table(pred_miceFast$imputations[index_NA] ,data_disc[index_NA,posit_y])
            #b = table(as.numeric(mice.impute.lda),data_disc[index_NA,posit_y])

            acc_a = sum(diag(a))/sum(a)

            #acc_b = sum(diag(b))/sum(b)
            #test_0 = (acc_a-acc_b)/acc_b >=  0

            test_0 = acc_a  >=  0.3

            ### grouping variable

            index_sort = sort(data_disc_NA[,posit_grs],index.return=TRUE)$ix

            index_rev = sort(sort(data_disc_NA[,posit_grs],index.return=TRUE)$ix,index.return=TRUE)$ix

            data_disc_NA_sort = data_disc_NA[index_sort,]

            #pred_Rbase = NULL
            #for(i in unique(data_disc_NA_sort[,posit_grs])){
            #  sub = data_disc_NA_sort[,posit_grs]==i
            #  temp = data_disc_NA_sort[sub,]
            #  pred = mice.impute.lda(temp[,posit_y],!temp[,posit_NA],temp[,posit_x])
            #  pred_Rbase = c(pred_Rbase,as.numeric(pred))
            #}

            #pred_dplyr = data_disc_NA_sort %>%
            #  as.data.frame() %>%
            #  group_by(group) %>%
            #  do(im = mice.impute.lda(as.matrix(.[,posit_y]),!.$index_NA,as.matrix(.[,posit_x]))) %>%
            #  tidy(im) %>%  arrange(group) %>% ungroup()%>% select(x) %>% unlist() %>% as.numeric()

            data = data_disc_NA_sort[,c(posit_y,posit_x)]
            g = data_disc_NA_sort[,posit_grs]

            model = new(miceFast)
            model$set_data(data)
            model$set_g(g)
            pred_miceFast =  model$impute("lda",posit_y,posit_x)

            rm(model)

            true_y = data_disc[index_sort,][index_NA[index_sort],posit_y]

            a = table(pred_miceFast$imputations[as.logical(pred_miceFast$index_imp)],true_y)
            #b = table(pred_dplyr,true_y)
            #c = table(pred_Rbase,true_y)

            acc_a = sum(diag(a))/sum(a)

            #acc_b = sum(diag(b))/sum(b)

            #test_1 = (acc_a-acc_b)/acc_b >= 0

            test_1 = acc_a >= 0.3

            #######################Binom

            #mice.impute.lda = mice.impute.lda(data_bin[,posit_y],!index_NA,data_bin[,posit_x])

            data = data_bin_NA[,c(posit_y,posit_x)]

            model = new(miceFast)
            model$set_data(data)
            pred_miceFast =  model$impute("lda",posit_y,posit_x)

            rm(model)

            a = table(pred_miceFast$imputations[index_NA] ,data_bin[index_NA,posit_y])
            #b = table(mice.impute.lda,data_bin[index_NA,posit_y])

            acc_a = sum(diag(a))/sum(a)

            #acc_b = sum(diag(b))/sum(b)

            #test_2 = (acc_a-acc_b)/acc_b >=  0
            test_2 = acc_a >=  0.5

            #####################Continous - LM Noise


            #mice.impute.norm.nob = mice.impute.norm.nob(data_con[,posit_y],!index_NA,data_con[,posit_x])

            data = data_con_NA[,c(posit_y,posit_x)]

            model = new(miceFast)
            model$set_data(data)
            pred_miceFast =  model$impute("lm_noise",posit_y,posit_x)

            rm(model)

            a = sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
            #b = sum((mice.impute.norm.nob - data_con[index_NA,posit_y])^2)

            #test_3 = abs((round(a)-round(b)))/round(a) < 0.1
            test_3 = round(a)/sum(index_NA) < 0.5

            #####################Continous - LM Bayes


            #mice.impute.norm.bayes = mice.impute.norm(data_con[,posit_y],!index_NA,data_con[,posit_x])

            data = data_con_NA[,c(posit_y,posit_x)]

            model = new(miceFast)
            model$set_data(data)
            pred_miceFast =  model$impute("lm_bayes",posit_y,posit_x)

            rm(model)

            a=sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
            #b=sum((mice.impute.norm.bayes - data_con[index_NA,posit_y])^2)

            #test_4 = abs((round(a)-round(b)))/round(a) < 0.1
            test_4 = round(a)/sum(index_NA) < 0.5

            #####################Continous - LM Predict


            #mice.impute.norm.pred = mice.impute.norm.predict(data_con[,posit_y],!index_NA,data_con[,posit_x])

            data = data_con_NA[,c(posit_y,posit_x)]

            model = new(miceFast)
            model$set_data(data)
            pred_miceFast =  model$impute("lm_pred",posit_y,posit_x)

            rm(model)

            a=sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
            #b=sum((mice.impute.norm.pred - data_con[index_NA,posit_y])^2)

            #test_5 = abs((round(a)-round(b)))/round(a) < 0.01
            test_5 = round(a)/sum(index_NA) < 0.5

            ## grouping variable

            index_sort = sort(data_con_NA[,posit_grs],index.return=TRUE)$ix

            index_rev = sort(sort(data_con_NA[,posit_grs],index.return=TRUE)$ix,index.return=TRUE)$ix

            data_con_NA_sort = data_con_NA[index_sort,]

            #pred_Rbase = NULL
            #for(i in unique(data_con_NA_sort[,posit_grs])){
            #  sub = data_con_NA_sort[,posit_grs]==i
             # temp = data_con_NA_sort[sub,]
             # pred = mice.impute.norm.predict(as.matrix(temp[,posit_y]),!temp[,posit_NA],as.matrix(temp[,posit_x]))
             # pred_Rbase = c(pred_Rbase,pred)
            #}

            #pred_dplyr = data_con_NA_sort %>%
            #  as.data.frame() %>%
            #  group_by(group) %>%
            #  do(im = mice.impute.norm.predict(as.matrix(.[,posit_y]),!.$index_NA,as.matrix(.[,posit_x]))) %>%
            #  tidy(im)  %>% ungroup()%>% select(x) %>% unlist() %>% as.numeric()

            data = cbind(data_con_NA_sort[,c(posit_y,posit_x)],1)
            g = data_con_NA_sort[,posit_grs]

            model = new(miceFast)
            model$set_data(data)
            model$set_g(g)
            pred_miceFast =  model$impute("lm_pred",posit_y,c(posit_x,max(posit_x)+1))

            rm(model)

            true_y = data_con[index_sort,][index_NA[index_sort],posit_y]

            a = sum((pred_miceFast$imputations[as.logical(pred_miceFast$index_imputed)]-true_y)^2)
            # = sum((pred_dplyr-true_y)^2)
            #c = sum((pred_Rbase-true_y)^2)

            #test_6 = abs((round(a)-round(b)))/round(a) < 0.01
            test_6 = round(a)/sum(index_NA) < 0.5

            expect_true(all(c(test_0,test_1,test_2,test_3,test_4,test_5,test_6)))

          })
