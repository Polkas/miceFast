
context("miceFast-vif")


test_that("VIF",
          {

  collin_x = function(y,x){

  full_rows = !is.na(y)&complete.cases(x)

  XXinv=solve(crossprod(scale(x[full_rows,],scale=F)))

  vifs = sapply(1:ncol(XXinv),function(i) det(XXinv[i,i,drop=FALSE])*det(XXinv[-i,-i,drop=FALSE])/det(XXinv))

  vifs2 = sapply(1:ncol(XXinv),function(i) Reduce("*",svd(XXinv[i,i,drop=FALSE])$d)*Reduce("*",svd(XXinv[-i,-i,drop=FALSE])$d)/Reduce("*",svd(XXinv)$d))

  list(vifs,vifs2)

}

airquality2 = airquality

airquality2$Temp2 = airquality2$Temp**2

#vif_car = car::vif(lm(Ozone ~ ., data=airquality2))

vif_R = collin_x(airquality2[,1],as.matrix(airquality2[,-1]))

model = new(miceFast)
model$set_data(as.matrix(airquality2))

vif_det_miceFast = model$vifs(1,c(2,3,4,5,6,7))

all_vifs_sd = apply(do.call(rbind,list(vif_R[[1]],vif_R[[2]],as.vector(vif_det_miceFast))),2,sd)

expect_true(all(all_vifs_sd<0.01))

          })
