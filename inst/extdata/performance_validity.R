#library(devtools)
#install_github("Rcpp","RcppCore")
#install.packages("pacman)
library(pacman)

p_load(Rcpp,
       mice,
       tidyverse,
       ggthemes,
       broom,
       miceFast)

set.seed(12345)

#parameters

power = 5 # power of 10 - number of observations - should be adjusted to a computer capabilities

nr_var = 7 #CHANGE - only if you generate a bigger corr matrix:  number of variables - independent and one dependent

grs = 10**(power-3) # grouping variable - number of groups

iters = 100 # number of iterations for benchmarking

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

colnames(data_bin) = c("y",paste0("x",1:(nr_var-2)),"group")
colnames(data_disc) = c("y",paste0("x",1:(nr_var-2)),"group")
colnames(data_con) = c("y",paste0("x",1:(nr_var-2)),"group")

## NA index

index_NA = 1:nrow(data_con) %in% sample(1:nrow(data_con),10^(power-1))

fill_NA = function(v,index_NA){

  v[index_NA] = NA

  v
}

######################Discrete

mice.impute.lda = mice.impute.lda(data_disc[,1],!index_NA,data_disc[,c(2:nr_var)])

model = new(miceFast,cbind(fill_NA(data_disc[,1],index_NA),data_disc[,2:nr_var]))
pred_miceFast =  model$impute("lda",1,2:nr_var)
#index_NA = pred_miceFast$index_NA

table(pred_miceFast$imputations[index_NA] ,data_disc[index_NA,1])
table(as.numeric(mice.impute.lda),data_disc[index_NA,1])

m1 = microbenchmark::microbenchmark(R=mice.impute.lda(data_disc[,1],!index_NA,data_disc[,c(2:nr_var)]),
                               miceFast={
                                 model = new(miceFast,cbind(fill_NA(data_disc[,1],index_NA),data_disc[,2:nr_var]))
                                 pred_miceFast =  model$impute("lda",1,2:nr_var)
                               },
                               times=iters)
m1

g1 = autoplot(m1,log=FALSE)+theme_economist()+ ggtitle("LDA discrete - without grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g1.png",g1)

### grouping variable


g = data_disc[,nr_var]

data_disc[,nr_var] = floor(pnorm(g)*grs)

data_disc = cbind(data_disc,index_NA)

data_disc = data_disc[order(data_disc[,nr_var]),]

gr= data_disc[,nr_var]

index_NA = as.logical(data_disc[,(nr_var+1)])

pred_Rbase = NULL
for(i in unique(data_disc[,nr_var])){
    temp = data_disc[data_disc[,nr_var]==i,]
    pred = mice.impute.lda(as.vector(temp[,1]),!temp[,(nr_var+1)],as.matrix(temp[,c(2:(nr_var-1))]))
    pred_Rbase = c(pred_Rbase,as.numeric(pred))
}

pred_dplyr = data_disc %>%
  as.data.frame() %>%
  group_by(group) %>%
  do(im = mice.impute.lda(as.matrix(.[,1]),!.$index_NA,as.matrix(.[,c(2:(nr_var-1))]))) %>%
  tidy(im)  %>% ungroup()%>% select(x) %>% unlist() %>% as.numeric()

model = new(miceFast,cbind(fill_NA(data_disc[,1],index_NA),data_disc[,2:(nr_var-1)]),gr,TRUE)
pred_miceFast =  model$impute("lda",1,2:(nr_var-1))
#index_NA = pred_miceFast$index_NA


table(pred_miceFast$imputations[index_NA],data_disc[index_NA,1])
table(pred_dplyr,data_disc[index_NA,1])
table(pred_Rbase,data_disc[index_NA,1])

##Performance

m2 = microbenchmark::microbenchmark(
  dplyr={
  pred_dplyr = data_disc %>%
    as.data.frame() %>%
    group_by(group) %>%
    do(im = mice.impute.lda(.$y,!.$index_NA,as.matrix(.[,c(2:(nr_var-1))]))) %>%
    tidy(im)  %>%
    ungroup()%>%
    select(x) %>%
    unlist() %>%
    as.numeric()},
  R_base={
      pred_all = NULL
      for(i in unique(data_disc[,nr_var])){
        temp = data_disc[data_disc[,nr_var]==i,]
        pred = mice.impute.lda(as.matrix(temp[,1]),!temp[,(nr_var+1)],as.matrix(temp[,c(2:(nr_var-1))]))
        pred_all = c(pred_all,as.numeric(pred))}},
  miceFast={
    model = new(miceFast,cbind(fill_NA(data_disc[,1],index_NA),data_disc[,2:(nr_var-1)]),gr,TRUE)
    pred_miceFast =  model$impute("lda",1,2:(nr_var-1))},
  times=iters)

m2

g2 = autoplot(m2,log=FALSE)+theme_economist()+ ggtitle("LDA discrete - with grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g2.png",g2)


#######################Binom

mice.impute.lda = mice.impute.lda(data_bin[,1],!index_NA,data_bin[,c(2:(nr_var-1))])

model = new(miceFast,cbind(fill_NA(data_bin[,1],index_NA),data_bin[,2:(nr_var-1)]))
pred_miceFast =  model$impute("lda",1,2:(nr_var-1))
#index_NA = pred_miceFast$index_NA

table(pred_miceFast$imputations[index_NA] ,data_bin[index_NA,1])
table(mice.impute.lda,data_bin[index_NA,1])

m3 = microbenchmark::microbenchmark(R=mice.impute.lda(data_bin[,1],!index_NA,data_bin[,c(2:(nr_var-1))]),
                               miceFast={
                                 model = new(miceFast,cbind(fill_NA(data_bin[,1],index_NA),data_bin[,2:(nr_var-1)]))
                                 pred_miceFast =  model$impute("lda",1,2:(nr_var-1))
                               },
                               times=iters)

m3

g3 = autoplot(m3,log=FALSE)+theme_economist()+ ggtitle("LDA binom - without grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g3.png",g3)


#####################Continous - LM Noise


mice.impute.norm.nob = mice.impute.norm.nob(data_con[,1],!index_NA,data_con[,c(2:(nr_var-1))])

model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),data_con[,2:(nr_var-1)]))
pred_miceFast =  model$impute("lm_noise",1,2:(nr_var-1))
#index_NA = pred_miceFast$index_NA

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,1])^2)
sum((mice.impute.norm.nob - data_con[index_NA,1])^2)

m4 = microbenchmark::microbenchmark(R = mice.impute.norm.nob(data_con[,1],!index_NA,data_con[,c(2:(nr_var-1))]),
                               miceFast={
                                 model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),data_con[,2:(nr_var-1)]))
                                 pred_miceFast =  model$impute("lm_noise",1,2:(nr_var-1))
                               },
                               times=iters)
m4

g4 = autoplot(m4,log=FALSE)+theme_economist()+ ggtitle("linear regression noise - without grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g4.png",g4)



#####################Continous - LM Bayes


mice.impute.norm.bayes = mice.impute.norm(data_con[,1],!index_NA,data_con[,c(2:(nr_var-1))])

model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),data_con[,2:(nr_var-1)]))
pred_miceFast =  model$impute("lm_bayes",1,2:(nr_var-1))
#index_NA = pred_miceFast$index_NA

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,1])^2)
sum((mice.impute.norm.bayes - data_con[index_NA,1])^2)

m5 = microbenchmark::microbenchmark(R = mice.impute.norm(data_con[,1],!index_NA,data_con[,c(2:(nr_var-1))]),
                               miceFast={
                                 model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),data_con[,2:(nr_var-1)]))
                                 pred_miceFast =  model$impute("lm_bayes",1,2:(nr_var-1))
                               },
                               times=iters)
m5

g5 = autoplot(m5,log=FALSE)+theme_economist()+ ggtitle("linear regression bayes - without grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g5.png",g5)



#####################Continous - LM Predict


mice.impute.norm.pred = mice.impute.norm.predict(data_con[,1],!index_NA,data_con[,c(2:(nr_var-1))])

model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),1,data_con[,c(2:(nr_var-1))]))
pred_miceFast =  model$impute("lm_pred",1,2:(nr_var))
#index_NA = pred_miceFast$index_NA

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA,1])^2)
sum((mice.impute.norm.pred - data_con[index_NA,1])^2)

m6 = microbenchmark::microbenchmark(R = {
  mice.impute.norm.predict(data_con[,1],!index_NA,data_con[,c(2:(nr_var-1))])
  }
,
miceFast={
  model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),1,data_con[,c(2:(nr_var-1))]))
  pred_miceFast =  model$impute("lm_pred",1,2:(nr_var-1))
},times=iters)

m6

g6 = autoplot(m6,log=FALSE)+theme_economist()+ ggtitle("linear regression predict - without grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g6.png",g6)


## grouping variable

g = data_con[,nr_var]

data_con[,nr_var] = floor(pnorm(g)*grs)

data_con = cbind(data_con,index_NA)

data_con = data_con[order(data_con[,nr_var]),]

gr= data_con[,nr_var]

index_NA = as.logical(data_con[,(nr_var+1)])

pred_Rbase = NULL
for(i in unique(data_con[,nr_var])){
  temp = data_con[data_con[,nr_var]==i,]
  pred = mice.impute.norm.predict(as.matrix(temp[,1]),!temp[,(nr_var+1)],as.matrix(temp[,c(2:(nr_var-1))]))
  pred_Rbase = c(pred_Rbase,pred)
}

pred_dplyr = data_con %>%
    as.data.frame() %>%
    group_by(group) %>%
    do(im = mice.impute.norm.predict(as.matrix(.[,1]),!.$index_NA,as.matrix(.[,c(2:(nr_var-1))]))) %>%
    tidy(im)  %>% ungroup()%>% select(x) %>% unlist() %>% as.numeric()

model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),1,data_con[,c(2:(nr_var-1))]),gr,TRUE)
pred_miceFast =  model$impute("lm_pred",1,2:(nr_var))
#index_NA = pred_miceFast$index_NA

sum((pred_miceFast$imputations[index_NA]-data_con[index_NA,1])^2)
sum((pred_dplyr-data_con[index_NA,1])^2)
sum((pred_Rbase-data_con[index_NA,1])^2)

##Performance

m7 = microbenchmark::microbenchmark(
  dplyr={
    pred_dplyr = data_con %>%
      as.data.frame() %>%
      group_by(group) %>%
      do(im = mice.impute.norm.predict(as.matrix(.[,1]),!.$index_NA,as.matrix(.[,c(2:(nr_var-1))]))) %>%
      tidy(im)  %>%
      ungroup()%>%
      select(x) %>%
      unlist() %>%
      as.numeric()
    },
  R_base={
    pred_Rbase = NULL
    for(i in unique(data_con[,nr_var])){
      temp = data_con[data_con[,nr_var]==i,]
      pred = mice.impute.norm.predict(as.matrix(temp[,1]),!temp[,(nr_var+1)],as.matrix(temp[,c(2:(nr_var-1))]))
      pred_Rbase = c(pred_Rbase,pred)
    }},
  miceFast={
    model = new(miceFast,cbind(fill_NA(data_con[,1],index_NA),1,data_con[,c(2:(nr_var-1))]),gr,TRUE)
    pred_miceFast =  model$impute("lm_pred",1,(nr_var))},
  times=iters)

m7

g7 = autoplot(m7,log=FALSE)+
  theme_economist()+
  ggtitle("linear regression predict - with grouping")

ggsave("C:/Users/user/Desktop/own_R_packages/miceFast/inst/extdata/images/g7.png",g7)

