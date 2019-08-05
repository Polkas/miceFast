// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

class corrData{

  int nr_cat = 2;
  int n_row;
  int n_col;
  arma::vec mu;
  arma::mat cors;

public:

  corrData(int _n_row,arma::vec _mu,arma::mat _cors):
      n_row(_n_row),
      n_col(_mu.n_elem),
      mu(_mu),
      cors(_cors){};

  corrData(int _nr_cat,int _n_row,arma::vec _mu,arma::mat _cors):
      nr_cat(_nr_cat),
      n_row(_n_row),
      n_col(_mu.n_elem),
      mu(_mu),
      cors(_cors){};

  arma::mat fill(std::string type);

};


std::map<std::string, int> types = {
    {"binom", 1},
    {"discrete", 2},
    {"contin", 3}
    };


arma::mat corrData::fill(std::string type){

  arma::mat X(n_row, n_col);

  std::mt19937 engine;  // Mersenne twister random number engine

  std::normal_distribution<double> distr(0.0, 1.0);

  X.imbue( [&]() { return distr(engine); } );

  arma::mat X_new(n_row, n_col);

  switch(types[type]){

  case 1:
      {

    X_new = X * arma::inv(arma::chol(arma::cor(X))) * arma::chol(cors);

    arma::colvec col0 = X_new.col(0);
    for(int i=0;i<n_row;i++){
            col0(i) = (Rcpp::stats::pnorm_0(col0(i),1,0))*2;
    }
    X_new.col(0) = arma::ceil(col0);

    break;
    }

  case 2:
      {
    X_new = X * arma::inv(arma::chol(arma::cor(X))) * arma::chol(cors);

    arma::colvec col0 = X_new.col(0);
    for(int i=0;i<n_row;i++){
            col0(i) = (Rcpp::stats::pnorm_0(col0(i),1,0))*nr_cat;
    }
    X_new.col(0) = arma::ceil(col0);

    break;
    }

  case 3:
      {
        X_new = X * arma::inv(arma::chol(arma::cor(X))) * arma::chol(cors);
        break;
    }

  default: {
      Rcpp::stop("binom,discrete or contin");
      }
  }

  for(int i=1; i<n_col;i++){X_new.col(i)=X_new.col(i)+mu(i);}

  return X_new;
}


RCPP_MODULE(corrData){
  using namespace Rcpp ;

     class_<corrData>("corrData")
    .constructor<int,arma::vec,arma::mat>()
    .constructor<int,int,arma::vec,arma::mat>()
    .method("fill", &corrData::fill)

  ;}


