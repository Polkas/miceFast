#include <RcppArmadillo.h>
#include "miceFast.h"

//functions for getting indexes

arma::uvec get_index_full_R(arma::mat &x,int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  arma::uvec index_Full = arma::find((full_y + full_x) == 2);

  return index_Full;
}

arma::uvec get_index_NA_R(arma::mat &x,int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  arma::uvec index(x.n_rows,arma::fill::zeros);

  for(unsigned int i=0;i<x.n_rows;i++){

    if( (full_y(i)==0) && (full_x(i)==1) ){index(i) = 1;}

  }

  arma::uvec index_NA = arma::find(index==1);

  return index_NA;
}

arma::colvec impute_raw_R(arma::mat &x,std::string s, int posit_y,arma::uvec posit_x,int times){

  typedef arma::colvec (*pfunc)(arma::colvec&,arma::mat&,arma::mat&,int);
  std::map<std::string, pfunc> funMap = {
    {"lda",fastLda},
    {"lm_pred",fastLm_pred},
    {"lm_noise",fastLm_noise},
    {"lm_bayes",fastLm_bayes}};

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  arma::uvec index_full = get_index_full_R(x,posit_y, posit_x);
  arma::uvec index_NA = get_index_NA_R(x,posit_y, posit_x);

  if(!(index_NA.n_elem==0)){

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.rows(index_full);

  pfunc f = funMap[s];
  arma::colvec pred = (*f)(Y_full,X_full,X_NA,times);

  Y.rows(index_NA) = pred;

  }

  return Y;
}


arma::colvec imputeW_R(arma::mat &x,std::string s,int posit_y,arma::uvec posit_x,arma::colvec w,int times){

  typedef arma::colvec (*pfuncw)(arma::colvec&,arma::mat&,arma::colvec&,arma::mat&,int);
  std::map<std::string, pfuncw> funMapw = {{"lm_pred",fastLm_weighted},
  {"lm_noise",fastLm_weighted_noise},
  {"lm_bayes",fastLm_weighted_bayes}};

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  arma::uvec index_full = get_index_full_R(x,posit_y, posit_x);
  arma::uvec index_NA = get_index_NA_R(x,posit_y, posit_x);

  if(w.has_nan()){Rcpp::stop("There is NA values for weights");}
  if(arma::any(w<0)){Rcpp::stop("There are negative values for the weights variable");}

  if(!(index_NA.n_elem==0)){

  //dividing data to NA and full

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.elem(index_full);
  arma::colvec w_full = w.elem(index_full);

  pfuncw f = funMapw[s];
  arma::colvec pred = (*f)(Y_full,X_full,w_full,X_NA,times);

  Y.rows(index_NA) = pred;

  }

  return Y;
}

//eval vifs

//' \code{VIF} function for assessing VIF which helps to identify a collinearity problem
//'
//' @param x a numeric matrix - a numeric matrix with variables
//' @param posit_y an integer - a position of dependent variable
//' @param posit_x an integer vector - positions of independent variables
//'
//' @return load a numeric vector with VIF for all variables provided by posit_x
//'
//' @seealso \code{\link{fill_NA}} \code{\link{fill_NA_N}}
//'
//' @examples
//' \dontrun{
//' library(miceFast)
//' library(data.table)
//'
//' airquality2 = airquality
//' airquality2$Temp2 = airquality2$Temp**2
//' #install.packages("car")
//' #car::vif(lm(Ozone ~ ., data=airquality2))
//'
//'
//' data_DT = data.table(airquality2)
//' # VIF for variables at 1,3,4 positions - you include a y position to consider its NA values
//' data_DT[,.(vifs=VIF(x=as.matrix(.SD),
//'                     posit_y=1,
//'                     posit_x=c(2,3,4,5,6,7)))]
//'
//' ######################
//' #OR using OOP miceFast
//' ######################
//'
//' airquality2_mat = as.matrix(airquality2)
//' model = new(miceFast)
//' model$set_data(airquality2_mat)
//'
//' as.vector(model$vifs(1,c(2,3,4,5,6,7)))
//'
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec VIF(arma::mat &x,int posit_y,arma::uvec posit_x){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  arma::mat x_cols = x.cols(posit_x - 1);

  arma::uvec full_rows = get_index_full_R(x,posit_y - 1,posit_x - 1);

  arma::mat full_x = x_cols.rows(full_rows);

  arma::mat col_means = arma::mean(full_x,0);

  full_x.each_row() -= col_means;

  arma::mat XXinv = arma::inv(full_x.t()*full_x);

  int Ncol = full_x.n_cols;

  arma::vec vifs(Ncol);

  double det_XXinv = arma::det(XXinv);

  for(int i=0;i<Ncol;i++){

    arma::mat XXinv_small = XXinv;

    XXinv_small.shed_row(i);

    XXinv_small.shed_col(i);

    vifs(i)= XXinv(i,i) * arma::det(XXinv_small)/det_XXinv;

  }

  return vifs;
};


//' \code{fill_NA_N} function for an multiple imputations purpose. Multiple imputations to fill the missing data.
//' Non missing independent variables are used to approximate a missing observations for a dependent variable.
//' Quantitative models were built under Rcpp packages and the C++ library Armadillo.
//'
//' @param x a numeric matrix - a numeric matrix with variables
//' @param model a character - posibble options ("lm_bayes","lm_noise")
//' @param posit_y an integer - a position of dependent variable
//' @param posit_x an integer vector - positions of independent variables
//' @param w  a numeric vector - a weighting variable - only positive values
//' @param times an integer - a number of multiple imputations - default 10
//'
//' @return load average of N imputations in a numeric vector format
//'
//' @seealso \code{\link{fill_NA}} \code{\link{VIF}}
//'
//' @examples
//' \dontrun{
//' library(miceFast)
//' library(data.table)
//' library(magrittr)
//'
//' data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality),
//'              # a numeric vector - positive values
//'              weights = round(rgamma(nrow(airquality),3,3),1),
//'              # as.numeric is needed only for OOP miceFast - see on next pages
//'              groups = airquality[,5])
//' data_DT = data.table(data)
//'
//' # simple mean imputation - intercept at position 6
//' data_DT[,Ozone_imp:=fill_NA(x=as.matrix(.SD),
//'                            model="lm_pred",
//'                            posit_y=1,
//'                            posit_x=c(6),w=.SD[['weights']]),by=.(groups)] %>%
//' # avg of 10 multiple imputations - last posit_x equal to 9 not 10
//' # because the groups variable is not included in .SD
//' .[,Solar_R_imp:=fill_NA_N(as.matrix(.SD),
//'                          model="lm_bayes",
//'                          posit_y=2,
//'                          posit_x=c(3,4,5,6,9),w=.SD[['weights']],times=10),by=.(groups)]
//'
//' head(data_DT,10)
//'
//' ######################
//' #OR using OOP miceFast
//' ######################
//'
//' data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality))
//' weights = rgamma(nrow(data),3,3) # a numeric vector - positive values
//' #a numeric vector not integers - positive values - sorted increasingly
//' groups = as.numeric(airquality[,5])
//' #a numeric vector not integers - positive values - not sorted
//' #groups = as.numeric(sample(1:8,nrow(data),replace=T))
//'
//' model = new(miceFast)
//' model$set_data(data) # providing data by a reference
//' model$set_w(weights) # providing by a reference
//' model$set_g(groups)  # providing by a reference
//'
//' #impute adapt to provided parmaters like w or g
//' #Simple mean - permanent imputation at the object and data
//' #variable will be replaced by imputations
//' model$update_var(1,model$impute("lm_pred",1,c(6))$imputations)
//'
//' model$update_var(2,model$impute_N("lm_bayes",2,c(1,3,4,5,6),10)$imputations)
//'
//' #Printing data and retrieving an old order if data was sorted by the grouping variable
//' head(cbind(model$get_data(),model$get_g(),model$get_w())[order(model$get_index()),],3)
//' #the same
//' head(cbind(data,groups,weights)[order(model$get_index()),],3)
//'
//'
//' }
//'
//' @export
//'
// [[Rcpp::export]]
arma::colvec fill_NA_N(arma::mat &x, std::string model, int posit_y,arma::uvec posit_x,arma::colvec w=0,int times=10){

  if( !(model.compare("lm_bayes") == 0) && !(model.compare("lm_noise") == 0)){Rcpp::stop("Works only for `lm_bayes` and `lm_noise` models");}
  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg;

  if(w(0)==0 || (model.compare("lda") == 0)){
    pred_avg = impute_raw_R(x,model,posit_y,posit_x,times);
  } else{
    pred_avg = imputeW_R(x,model,posit_y,posit_x,w,times);
  }

  //index
  return pred_avg;
}

//' \code{fill_NA} function for an imputations purpose. Regular imputations to fill the missing data.
//' Non missing independent variables are used to approximate a missing observations for a dependent variable.
//' Quantitative models were built under Rcpp packages and the C++ library Armadillo.
//'
//' @param x a numeric matrix - a numeric matrix with variables
//' @param model a character - posibble options ("lda","lm_pred","lm_bayes","lm_noise")
//' @param posit_y an integer - a position of dependent variable
//' @param posit_x an integer vector - positions of independent variables
//' @param w  a numeric vector - a weighting variable - only positive values
//'
//' @return load imputations in a numeric vector format
//'
//' @seealso \code{\link{fill_NA_N}}  \code{\link{VIF}}
//'
//' @examples
//' \dontrun{
//' library(miceFast)
//' library(data.table)
//' library(magrittr)
//'
//' data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality),
//'              # a numeric vector - positive values
//'              weights = round(rgamma(nrow(airquality),3,3),1),
//'              # as.numeric is needed only for OOP miceFast - see on next pages
//'              groups = airquality[,5])
//' data_DT = data.table(data)
//'
//' # simple mean imputation - intercept at position 6
//' data_DT[,Ozone_imp:=fill_NA(x=as.matrix(.SD),
//'                            model="lm_pred",
//'                            posit_y=1,
//'                            posit_x=c(6),w=.SD[['weights']]),by=.(groups)] %>%
//' # avg of 10 multiple imputations - last posit_x equal to 9 not 10
//' # because the groups variable is not included in .SD
//' .[,Solar_R_imp:=fill_NA_N(as.matrix(.SD),
//'                          model="lm_bayes",
//'                          posit_y=2,
//'                          posit_x=c(3,4,5,6,9),w=.SD[['weights']],times=10),by=.(groups)]
//'
//' head(data_DT,10)
//'
//' ######################
//' #OR using OOP miceFast
//' ######################
//'
//' data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality))
//' weights = rgamma(nrow(data),3,3) # a numeric vector - positive values
//' #a numeric vector not integers - positive values - sorted increasingly
//' groups = as.numeric(airquality[,5])
//' #a numeric vector not integers - positive values - not sorted
//' #groups = as.numeric(sample(1:8,nrow(data),replace=T))
//'
//' model = new(miceFast)
//' model$set_data(data) # providing data by a reference
//' model$set_w(weights) # providing by a reference
//' model$set_g(groups)  # providing by a reference
//'
//' #impute adapt to provided parmaters like w or g
//' #Simple mean - permanent imputation at the object and data
//' #variable will be replaced by imputations
//' model$update_var(1,model$impute("lm_pred",1,c(6))$imputations)
//'
//' model$update_var(2,model$impute_N("lm_bayes",2,c(1,3,4,5,6),10)$imputations)
//'
//' #Printing data and retrieving an old order if data was sorted by the grouping variable
//' head(cbind(model$get_data(),model$get_g(),model$get_w())[order(model$get_index()),],3)
//' #the same
//' head(cbind(data,groups,weights)[order(model$get_index()),],3)
//'
//' }
//'
//' @export
// [[Rcpp::export]]
arma::colvec fill_NA(arma::mat &x,std::string model, int posit_y,arma::uvec posit_x,arma::colvec w=0){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg;

  if(w(0)==0 || (model.compare("lda") == 0)){
    pred_avg = impute_raw_R(x,model,posit_y,posit_x,1);
  } else {
    pred_avg = imputeW_R(x,model,posit_y,posit_x,w,1);
  }

  //index
  return pred_avg;
}
