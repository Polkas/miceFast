#library(devtools)
#install_github("Rcpp","RcppCore")
#install.packages("pacman")
library(pacman)

p_load(Rcpp,
       mice,
       tidyverse,
       ggthemes,
       broom,
       miceFast,
       data.table,
       microbenchmark,
       car)

set.seed(1234)
#parameters
power = 5 # power of 10 - number of observations - should be adjusted to a computer capabilities
nr_var = 20 #CHANGE - only if you generate a bigger corr matrix:  number of variables - independent and one dependent
grs = max(c(10**(power-3),10)) # grouping variable - number of groups
iters = 30 # number of iterations for benchmarking
## generete example - data

##positive-defined correlation matrix
#cors = matrix(c(1,0.6,0.7,0.4,0.4,0.5,0.35,
#                NA,1,0.2,0.05,0.1,0.12,0.15,
#                NA,NA,1,0.15,0.15,0.1,0.08,
#                NA,NA,NA,1,0.12,0.15,0.1,
#                NA,NA,NA,NA,1,0.15,0.2,
#                NA,NA,NA,NA,NA,1,0.15,
#                NA,NA,NA,NA,NA,NA,1),7,7,byrow = T)

#cors[lower.tri(cors)] = t(cors)[lower.tri(cors)]
# automatic corr matrix

vars_mat = diag(nr_var)
vars_mat[vars_mat==0] = 0.3
covs = stats::rWishart(3,nr_var,vars_mat)
covs = apply(covs,1:2,mean)/(nr_var+10)
cors = diag(1/sqrt(diag(covs))) %*% covs%*%diag(1/sqrt(diag(covs)))

#diag(cors) = 1
#cors
##

model = new(corrData,10,10^power,rep(0,nr_var),cors)
data_bin = model$fill("binom")
data_disc = model$fill("discrete")
data_con = model$fill("contin")

n_vars = ncol(cors)
posit_y = 1
posit_x = 2:(n_vars-2)
posit_w = n_vars-1
posit_grs = n_vars
posit_NA = n_vars+1
posit_index = n_vars+2

## NA index

index_NA = 1:nrow(data_con) %in% sample(1:nrow(data_con),2*10^(power-1))
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

index = 1:(10**power)

data_disc_NA = cbind(fill_by_NA(data_disc[,posit_y],index_NA),data_disc[,posit_x],w_d,group_d,index_NA,index)
data_con_NA = cbind(fill_by_NA(data_con[,posit_y],index_NA),data_con[,posit_x],w_c,group_c,index_NA,index)
data_bin_NA = cbind(fill_by_NA(data_bin[,posit_y],index_NA),data_bin[,posit_x],w_b,group_b,index_NA,index)

colnames(data_bin_NA) = c("y",paste0("x",posit_x),"weights","group","index_NA","index")
colnames(data_disc_NA) = c("y",paste0("x",posit_x),"weights","group","index_NA","index")
colnames(data_con_NA) = c("y",paste0("x",posit_x),"weights","group","index_NA","index")


######################Discrete

mice.impute.lda = mice.impute.lda(data_disc[,posit_y],!index_NA,data_disc[,posit_x])

model = new(miceFast)
data = data_disc_NA[,c(posit_y,posit_x)]
model$set_data(data)
pred_miceFast =  model$impute("lda",posit_y,posit_x)
rm(model)

table(pred_miceFast$imputations[index_NA] ,data_disc[index_NA,posit_y])
table(as.numeric(mice.impute.lda),data_disc[index_NA,posit_y])

m1 = microbenchmark::microbenchmark(R_mice=mice.impute.lda(data_disc[,posit_y],!index_NA,data_disc[,posit_x]),
                               miceFast={
                                 model = new(miceFast)
                                 model$set_data(data_disc_NA[,c(posit_y,posit_x)])
                                 pred_miceFast =  model$impute("lda",posit_y,posit_x)
                                 rm(model)
                               },
                               times=iters)
m1

g1 = autoplot(m1,log=FALSE)+theme_economist()+ ggtitle("LDA discrete - without grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g1.png",g1)

### grouping variable

index_sort = sort(data_disc_NA[,posit_grs],index.return=TRUE)$ix

index_rev = sort(sort(data_disc_NA[,posit_grs],index.return=TRUE)$ix,index.return=TRUE)$ix

data_disc_NA_sort = data_disc_NA[index_sort,]

data_disc_NA_sort_DF = as.data.frame(data_disc_NA_sort)

true_y = data_disc[index_sort,][index_NA[index_sort],posit_y]

pred_Rbase = NULL
for(i in unique(data_disc_NA_sort[,posit_grs])){
    sub = data_disc_NA_sort[,posit_grs]==i
    temp = data_disc_NA_sort[sub,]
    pred = mice.impute.lda(temp[,posit_y],!temp[,posit_NA],temp[,posit_x])
    pred_Rbase = c(pred_Rbase,as.numeric(pred))
}

table(pred_Rbase,true_y)

pred_dplyr = data_disc_NA_sort_DF %>%
  group_by(group) %>%
  do(im = mice.impute.lda(as.matrix(.[,posit_y]),!.$index_NA,as.matrix(.[,posit_x]))) %>%
  tidy(im) %>%  arrange(group) %>% ungroup()%>% select(x) %>% unlist()

table(pred_dplyr,true_y)

data_disc_NA_sort_DT = data.table(data_disc_NA_sort)
pred_datatable = data_disc_NA_sort_DT[,{im=mice.impute.lda(as.matrix(.SD[['y']]),!index_NA,as.matrix(.SD[,posit_x,with=F]))},by=.(group)]

table(pred_datatable[['V1']],true_y)

pred_datatable_miceFast = data_disc_NA_sort_DT[,{im=fill_NA(as.matrix(.SD),"lda",posit_y,posit_x)},by=.(group)]

table(pred_datatable_miceFast[['V1']][index_NA[index_sort]],true_y)

data = data_disc_NA_sort
g = data_disc_NA_sort[,posit_grs]

model = new(miceFast)
model$set_data(data)
model$set_g(g)
pred_miceFast =  model$impute("lda",posit_y,posit_x)
rm(model)

table(pred_miceFast$imputations[as.logical(pred_miceFast$index_imp)],true_y)

data = data_disc_NA
g = data_disc_NA[,posit_grs]
model = new(miceFast)
model$set_data(data_disc_NA)
model$set_g(g)
pred_miceFast_rotate =  model$impute("lda",posit_y,posit_x)

table(pred_miceFast_rotate$imputations[order(model$get_index())][as.logical(pred_miceFast$index_imp)[order(model$get_index())]],data_disc[index_NA,posit_y])

##Performance

m2 = microbenchmark::microbenchmark(
  dplyr_mice={
  pred_dplyr = data_disc_NA_sort_DF %>%
    group_by(group) %>%
    do(im = mice.impute.lda(as.matrix(.[,posit_y]),!.$index_NA,as.matrix(.[,posit_x]))) %>%
    tidy(im)  %>%
    ungroup()%>%
    select(x) %>%
    unlist() },
  DT_mice={pred_datatable= data_disc_NA_sort_DT[,{im=mice.impute.lda(as.matrix(.SD[['y']]),!index_NA,as.matrix(.SD[,posit_x,with=F]))},by=.(group)]},
  DT_miceFast = {pred_datatable_miceFast = data_disc_NA_sort_DT[,{im=fill_NA(as.matrix(.SD),"lda",posit_y,posit_x)},by=.(group)]},
  R_mice={
      pred_all = NULL
      for(i in unique(data_disc_NA_sort[,nr_var])){
        sub = data_disc_NA_sort[,posit_grs]==i
        temp = data_disc_NA_sort[sub,]
        pred = mice.impute.lda(temp[,posit_y],!temp[,posit_NA],temp[,posit_x])
        pred_Rbase = c(pred_Rbase,as.numeric(pred))}},
  miceFast={
    model = new(miceFast)
    model$set_data(data)
    model$set_g(g)
    pred_miceFast =  model$impute("lda",posit_y,posit_x)
    rm(model)
    },
  times=iters)

m2

g2 = autoplot(m2,log=FALSE)+theme_economist()+ ggtitle("LDA discrete - with grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g2.png",g2)


#######################Binom

mice.impute.lda = mice.impute.lda(data_bin[,posit_y],!index_NA,data_bin[,posit_x])

data = data_bin_NA[,c(posit_y,posit_x)]

model = new(miceFast)
model$set_data(data)
pred_miceFast =  model$impute("lda",posit_y,posit_x)
rm(model)

table(pred_miceFast$imputations[index_NA] ,data_bin[index_NA,posit_y])
table(mice.impute.lda,data_bin[index_NA,posit_y])

m3 = microbenchmark::microbenchmark(R_mice={mice.impute.lda = mice.impute.lda(data_bin[,posit_y],!index_NA,data_bin[,posit_x])},
                               miceFast={
                                 model = new(miceFast)
                                 model$set_data(data)
                                 pred_miceFast =  model$impute("lda",posit_y,posit_x)
                                 rm(model)

                               },
                               times=iters)

m3

g3 = autoplot(m3,log=FALSE)+theme_economist()+ ggtitle("LDA binom - without grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g3.png",g3)


#####################Continous - LM Noise


mice.impute.norm.nob = mice.impute.norm.nob(data_con[,posit_y],!index_NA,data_con[,posit_x])

data = data_con_NA[,c(posit_y,posit_x)]

model = new(miceFast)
model$set_data(data)
pred_miceFast =  model$impute("lm_noise",posit_y,posit_x)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
sum((mice.impute.norm.nob - data_con[index_NA,posit_y])^2)

m4 = microbenchmark::microbenchmark(R_mice ={mice.impute.norm.nob = mice.impute.norm.nob(data_con[,posit_y],!index_NA,data_con[,posit_x])},
                               miceFast={
                                 model = new(miceFast)
                                 model$set_data(data)
                                 pred_miceFast =  model$impute("lm_noise",posit_y,posit_x)
                                 rm(model)

                               },
                               times=iters)
m4

g4 = autoplot(m4,log=FALSE)+theme_economist()+ ggtitle("linear regression noise - without grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g4.png",g4)

#####################Continous - LM Bayes

mice.impute.norm.bayes = mice.impute.norm(data_con[,posit_y],!index_NA,data_con[,posit_x])

data = data_con_NA[,c(posit_y,posit_x)]

model = new(miceFast)
model$set_data(data)
pred_miceFast =  model$impute("lm_bayes",posit_y,posit_x)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
sum((mice.impute.norm.bayes - data_con[index_NA,posit_y])^2)

m5 = microbenchmark::microbenchmark(R_mice = mice.impute.norm(data_con[,posit_y],!index_NA,data_con[,posit_x]),
                               miceFast={
                                 model = new(miceFast)
                                 model$set_data(data)
                                 pred_miceFast =  model$impute("lm_bayes",posit_y,posit_x)
                                 rm(model)

                               },
                               times=iters)
m5

g5 = autoplot(m5,log=FALSE)+theme_economist()+ ggtitle("linear regression bayes - without grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g5.png",g5)

#####################Continous - LM Predict


mice.impute.norm.pred = mice.impute.norm.predict(data_con[,posit_y],!index_NA,data_con[,posit_x])

data = cbind(data_con_NA[,c(posit_y,posit_x)],1)

model = new(miceFast)
model$set_data(data)
pred_miceFast =  model$impute("lm_pred",posit_y,c(posit_x,max(posit_x)+1))
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
sum((mice.impute.norm.pred - data_con[index_NA,posit_y])^2)

m6 = microbenchmark::microbenchmark(R_mice = {
  mice.impute.norm.pred = mice.impute.norm.predict(data_con[,posit_y],!index_NA,data_con[,posit_x])
},
miceFast={
  model = new(miceFast)
  model$set_data(data)
  pred_miceFast =  model$impute("lm_pred",posit_y,posit_x)
  rm(model)

},times=iters)

m6

g6 = autoplot(m6,log=FALSE)+theme_economist()+ ggtitle("linear regression predict - without grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g6.png",g6)


## grouping variable

index_sort = sort(data_con_NA[,posit_grs],index.return=TRUE)$ix

index_rev = sort(sort(data_con_NA[,posit_grs],index.return=TRUE)$ix,index.return=TRUE)$ix

data_con_NA_sort = data_con_NA[index_sort,]

pred_Rbase = NULL
for(i in unique(data_con_NA_sort[,posit_grs])){
  sub = data_con_NA_sort[,posit_grs]==i
  temp = data_con_NA_sort[sub,]
  pred = mice.impute.norm.predict(as.matrix(temp[,posit_y]),!temp[,posit_NA],as.matrix(temp[,posit_x]))
  pred_Rbase = c(pred_Rbase,pred)
}

pred_dplyr = data_con_NA_sort %>%
    as.data.frame() %>%
    group_by(group) %>%
    do(im = mice.impute.norm.predict(as.matrix(.[,posit_y]),!.$index_NA,as.matrix(.[,posit_x]))) %>%
    tidy(im)  %>% ungroup()%>% select(x) %>% unlist() %>% as.numeric()

data_con_NA_sort_DT = data.table(data_con_NA_sort)

pred_datatable = data_con_NA_sort_DT[,{im=mice.impute.norm.predict(as.matrix(.SD[,1]),!index_NA,as.matrix(.SD[,posit_x,with=F]))},by=.(group)]

pred_datatable_miceFast = data_con_NA_sort_DT[,{im=fill_NA(cbind(as.matrix(.SD),1),"lm_pred",1,posit_x)},by=.(group)]

data = cbind(data_con_NA_sort[,c(posit_y,posit_x)],1)
g = data_con_NA_sort[,posit_grs]

model = new(miceFast)
model$set_data(data)
model$set_g(g)
pred_miceFast =  model$impute("lm_pred",posit_y,posit_x)
rm(model)

true_y = data_con[index_sort,][index_NA[index_sort],posit_y]

sum((pred_miceFast$imputations[as.logical(pred_miceFast$index_imputed)]-true_y)^2) # no intercept
sum((pred_dplyr-true_y)^2)
sum((pred_Rbase-true_y)^2)
sum((pred_datatable$V1-true_y)^2)
sum((pred_datatable_miceFast$V1[index_NA[index_sort]]-true_y)^2) #no intercept

##Performance

m7 = microbenchmark::microbenchmark(
  dplyr_mice={
    pred_dplyr = data_con_NA_sort %>%
      as.data.frame() %>%
      group_by(group) %>%
      do(im = mice.impute.norm.predict(as.matrix(.[,posit_y]),!.$index_NA,as.matrix(.[,posit_x]))) %>%
      tidy(im)  %>%
      ungroup()%>%
      select(x) %>%
      unlist() %>%
      as.numeric()
    },
  R_mice={
    pred_Rbase = NULL
    for(i in unique(data_con_NA_sort[,posit_grs])){
      sub = data_con_NA_sort[,posit_grs]==i
      temp = data_con_NA_sort[sub,]
      pred = mice.impute.norm.predict(as.matrix(temp[,posit_y]),!temp[,posit_NA],as.matrix(temp[,posit_x]))
      pred_Rbase = c(pred_Rbase,pred)
    }
    },
  miceFast={
    model = new(miceFast)
    model$set_data(data)
    model$set_g(g)
    pred_miceFast =  model$impute("lm_pred",posit_y,posit_x)
    rm(model)
  },
  times=iters,
  DT_mice={pred_datatable= data_disc_NA_sort_DT[,{im=mice.impute.norm.predict(as.matrix(.SD[,1]),!index_NA,as.matrix(.SD[,posit_x,with=F]))},by=.(group)]},
  DT_miceFast = {pred_datatable_DT=data_disc_NA_sort_DT[,{im=fill_NA(as.matrix(.SD),"lm_pred",posit_y,posit_x)},by=.(group)]}
)

m7

g7 = autoplot(m7,log=FALSE)+
  theme_economist()+
  ggtitle("linear regression predict - with grouping")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g7.png",g7)

####
####Multiple Imputations
####

mice.impute.norm.nob = rowMeans(sapply(1:10,function(x) mice.impute.norm.nob(data_con[,posit_y],!index_NA,data_con[,posit_x])))

data = data_con_NA[,c(posit_y,posit_x)]

model = new(miceFast)
model$set_data(data)
pred_miceFast =  model$impute_N("lm_noise",posit_y,posit_x,10)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
sum((mice.impute.norm.nob - data_con[index_NA,posit_y])^2)

m8 = microbenchmark::microbenchmark(R_mice ={mice.impute.norm.nob = rowMeans(sapply(1:10,function(x) mice.impute.norm.nob(data_con[,posit_y],!index_NA,data_con[,posit_x])))},
                                    miceFast={
                                      model = new(miceFast)
                                      model$set_data(data)
                                      pred_miceFast =  model$impute_N("lm_noise",posit_y,posit_x,10)
                                      rm(model)

                                    },
                                    times=iters)
m8

g8 = autoplot(m8,log=FALSE)+theme_economist()+ ggtitle("linear regression noise - without grouping - multiple 10")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g8.png",g8)

#####################Continous - LM Bayes - multiple

mice.impute.norm.bayes = rowMeans(sapply(1:10,function(x) mice.impute.norm(data_con[,posit_y],!index_NA,data_con[,posit_x])))

data = data_con_NA[,c(posit_y,posit_x)]

model = new(miceFast)
model$set_data(data)
pred_miceFast =  model$impute_N("lm_bayes",posit_y,posit_x,10)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,posit_y])^2)
sum((mice.impute.norm.bayes - data_con[index_NA,posit_y])^2)

m9 = microbenchmark::microbenchmark(R_mice = {mice.impute.norm.bayes = rowMeans(sapply(1:10,function(x) mice.impute.norm(data_con[,posit_y],!index_NA,data_con[,posit_x])))},
                                    miceFast={
                                      model = new(miceFast)
                                      model$set_data(data)
                                      pred_miceFast =  model$impute_N("lm_bayes",posit_y,posit_x,10)
                                      rm(model)

                                    },
                                    times=iters)
m9

g9 = autoplot(m9,log=FALSE)+theme_economist()+ ggtitle("linear regression bayes - without grouping - multiple 10")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g9.png",g9)

#####################Continous - VIFS

data.vifs = data.frame(data_con[,c(posit_y,posit_x)])
colnames(data.vifs) = letters[1:(nr_var-2)]
vifs_car = car::vif(lm(a~.,data=data.vifs))

model = new(miceFast)
model$set_data(data_con)
vifs_miceFast =  model$vifs(posit_y,posit_x)
rm(model)

vifs_car
vifs_miceFast

m10 = microbenchmark::microbenchmark(car = { car::vif(lm(a~.,data=data.vifs))},
                                    miceFast={
                                      model = new(miceFast)
                                      model$set_data(data_con)
                                      vifs_miceFast =  model$vifs(posit_y,posit_x)
                                      rm(model)

                                    },
                                    times=iters)
m10

g10 = autoplot(m10,log=FALSE)+theme_economist()+ ggtitle("vifs")

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g10.png",g10)

#plot for README/Intro

dats = bind_rows(list(
data.frame(m1) %>% mutate(model = "LDA discrete - without grouping"),
data.frame(m2) %>% mutate(model = "LDA discrete - with grouping"),
data.frame(m3) %>% mutate(model = "LDA binom - without grouping"),
data.frame(m4) %>% mutate(model = "linear regression noise - without grouping"),
data.frame(m5) %>% mutate(model = "linear regression bayes - without grouping"),
data.frame(m6) %>% mutate(model = "linear regression predict - without grouping"),
data.frame(m7) %>% mutate(model = "linear regression predict - with grouping"),
data.frame(m8) %>% mutate(model = "linear regression noise - without grouping - multiple 10"),
data.frame(m9) %>% mutate(model = "linear regression bayes - without grouping - multiple 10")
))

dats_plot = dats %>% group_by(model,expr) %>% summarise(mean_time_sec=mean(time/10**9)) %>%
  group_by(model) %>%
  mutate(relative_time=mean_time_sec/min(mean_time_sec)) %>%
  rename('package'='expr')

g_summary = ggplot(dats_plot,aes(model,relative_time,fill=package)) +
  geom_bar(stat="identity",position="dodge") +
  theme(axis.text.x= element_text(angle=90)) +
  ggtitle(paste0("Benchmarks - 10^" ,power, "obs (20% NA) ",nr_var,"vars - 10^3 groups"))

ggsave("C:/Users/user/Desktop/burk/own_R_packages/miceFast/inst/extdata/images/g_summary.png",g_summary)

