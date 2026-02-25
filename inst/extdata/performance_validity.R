# Performance & Validity Benchmarks for miceFast
# Run from the package root: Rscript inst/extdata/performance_validity.R
#
# Compares miceFast vs mice on every imputation model (LDA, lm_pred,
# lm_noise, lm_bayes, PMM, weighted variants, grouped, multiple) and VIF.
# Saves benchmark plots to inst/extdata/images/.

suppressPackageStartupMessages({
  library(Rcpp)
  library(mice)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(miceFast)
  library(data.table)
  library(microbenchmark)
  library(car)
  library(rprojroot)
})

set.seed(1234)

# Portable base_path: works from package root or from this script's dir
base_path <- tryCatch(
  rprojroot::find_package_root_file(),
  error = function(e) {
    normalizePath(
      file.path(dirname(sys.frame(1)$ofile), "..", ".."),
      mustWork = FALSE
    )
  }
)
dir.create(
  file.path(base_path, "inst", "extdata", "images"),
  recursive = TRUE,
  showWarnings = FALSE
)
cat("Saving plots to:", file.path(base_path, "inst", "extdata", "images"), "\n")

# parameters
power <- 5 # power of 10 - number of observations - should be adjusted to a computer capabilities
nr_var <- 10 # CHANGE - only if you generate a bigger corr matrix:  number of variables - independent and one dependent
grs <- max(c(10**(power - 3), 10)) # grouping variable - number of groups
iters <- 30 # number of iterations for benchmarking
## generete example - data

## positive-defined correlation matrix
# cors = matrix(c(1,0.6,0.7,0.4,0.4,0.5,0.35,
#                NA,1,0.2,0.05,0.1,0.12,0.15,
#                NA,NA,1,0.15,0.15,0.1,0.08,
#                NA,NA,NA,1,0.12,0.15,0.1,
#                NA,NA,NA,NA,1,0.15,0.2,
#                NA,NA,NA,NA,NA,1,0.15,
#                NA,NA,NA,NA,NA,NA,1),7,7,byrow = T)

# cors[lower.tri(cors)] = t(cors)[lower.tri(cors)]
# automatic corr matrix

vars_mat <- diag(nr_var)
vars_mat[vars_mat == 0] <- 0.3
covs <- stats::rWishart(3, nr_var, vars_mat)
covs <- apply(covs, 1:2, mean) / (nr_var + 10)
cors <- diag(1 / sqrt(diag(covs))) %*% covs %*% diag(1 / sqrt(diag(covs)))

# diag(cors) = 1
# cors
##

model <- new(corrData, 10, 10^power, rep(0, nr_var), cors)
data_bin <- model$fill("binom")
data_disc <- model$fill("discrete")
data_con <- model$fill("contin")

n_vars <- ncol(cors)
posit_y <- 1
posit_x <- 2:(n_vars - 2)
posit_w <- n_vars - 1
posit_grs <- n_vars
posit_NA <- n_vars + 1
posit_index <- n_vars + 2

## NA index

index_NA <- 1:nrow(data_con) %in% sample(1:nrow(data_con), 2 * 10^(power - 1))
fill_by_NA <- function(v, index_NA) {
  v[index_NA] <- NA
  v
}

######################

group_d <- floor(pnorm(data_disc[, posit_grs]) * grs)
group_c <- floor(pnorm(data_con[, posit_grs]) * grs)
group_b <- floor(pnorm(data_bin[, posit_grs]) * grs)

w_d <- abs(data_disc[, posit_w])
w_c <- abs(data_con[, posit_w])
w_b <- abs(data_bin[, posit_w])

index <- 1:(10**power)

data_disc_NA <- cbind(
  fill_by_NA(data_disc[, posit_y], index_NA),
  data_disc[, posit_x],
  w_d,
  group_d,
  index_NA,
  index
)
data_con_NA <- cbind(
  fill_by_NA(data_con[, posit_y], index_NA),
  data_con[, posit_x],
  w_c,
  group_c,
  index_NA,
  index
)
data_bin_NA <- cbind(
  fill_by_NA(data_bin[, posit_y], index_NA),
  data_bin[, posit_x],
  w_b,
  group_b,
  index_NA,
  index
)

col_names <- c(
  "y",
  paste0("x", posit_x),
  "weights",
  "group",
  "index_NA",
  "index"
)
colnames(data_bin_NA) <- col_names
colnames(data_disc_NA) <- col_names
colnames(data_con_NA) <- col_names

###################### Discrete

pred_lda_disc_mice <- mice::mice.impute.lda(
  data_disc[, posit_y],
  !index_NA,
  data_disc[, posit_x]
)

model <- new(miceFast)
data <- data_disc_NA[, c(posit_y, posit_x)]
model$set_data(data)
pred_miceFast <- model$impute("lda", posit_y, posit_x)
rm(model)

table(pred_miceFast$imputations[index_NA], data_disc[index_NA, posit_y])
table(as.numeric(pred_lda_disc_mice), data_disc[index_NA, posit_y])

m1 <- microbenchmark::microbenchmark(
  R_mice = mice::mice.impute.lda(
    data_disc[, posit_y],
    !index_NA,
    data_disc[, posit_x]
  ),
  miceFast = {
    model <- new(miceFast)
    model$set_data(data_disc_NA[, c(posit_y, posit_x)])
    pred_miceFast <- model$impute("lda", posit_y, posit_x)
    rm(model)
  },
  times = iters
)
m1

g1 <- autoplot(m1, log = FALSE) +
  theme_economist() +
  ggtitle("LDA discrete - without grouping")

ggsave(sprintf("%s/inst/extdata/images/g1.png", base_path), g1)

### grouping variable

index_sort <- sort(data_disc_NA[, posit_grs], index.return = TRUE)$ix

data_disc_NA_sort <- data_disc_NA[index_sort, ]

data_disc_NA_sort_DF <- as.data.frame(data_disc_NA_sort)

true_y <- data_disc[index_sort, ][index_NA[index_sort], posit_y]

pred_Rbase <- NULL
for (i in unique(data_disc_NA_sort[, posit_grs])) {
  sub <- data_disc_NA_sort[, posit_grs] == i
  temp <- data_disc_NA_sort[sub, ]
  pred <- mice::mice.impute.lda(
    temp[, posit_y],
    !temp[, posit_NA],
    temp[, posit_x]
  )
  pred_Rbase <- c(pred_Rbase, as.numeric(pred))
}

table(pred_Rbase, true_y)

pred_dplyr <- tibble(data_disc_NA_sort_DF) %>%
  group_by(group) %>%
  do(
    im = mice::mice.impute.lda(
      as.matrix(.[, posit_y]),
      !.$index_NA,
      as.matrix(.[, posit_x])
    )
  ) %>%
  unnest(im) %>%
  arrange(group) %>%
  ungroup() %>%
  select(im) %>%
  unlist()

table(pred_dplyr, true_y)

data_disc_NA_sort_DT <- data.table(data_disc_NA_sort)
pred_datatable <- data_disc_NA_sort_DT[,
  {
    im <- mice::mice.impute.lda(
      as.matrix(.SD[["y"]]),
      !index_NA,
      as.matrix(.SD[, posit_x, with = F])
    )
  },
  by = .(group)
]

table(pred_datatable[["V1"]], true_y)

pred_datatable_miceFast <- data_disc_NA_sort_DT[,
  {
    im <- fill_NA(as.matrix(.SD), "lda", posit_y, posit_x)
  },
  by = .(group)
]

table(pred_datatable_miceFast[["V1"]][index_NA[index_sort]], true_y)

data <- data_disc_NA_sort
g <- data_disc_NA_sort[, posit_grs]

model <- new(miceFast)
model$set_data(data)
model$set_g(g)
pred_miceFast <- model$impute("lda", posit_y, posit_x)
rm(model)

table(pred_miceFast$imputations[as.logical(pred_miceFast$index_imp)], true_y)

data <- data_disc_NA
g <- data_disc_NA[, posit_grs]
model <- new(miceFast)
model$set_data(data_disc_NA)
model$set_g(g)
pred_miceFast_rotate <- model$impute("lda", posit_y, posit_x)

table(
  pred_miceFast_rotate$imputations[order(model$get_index())][as.logical(
    pred_miceFast$index_imp
  )[order(model$get_index())]],
  data_disc[index_NA, posit_y]
)

## Performance

m2 <- microbenchmark::microbenchmark(
  dplyr_mice = {
    pred_dplyr <- data_disc_NA_sort_DF %>%
      group_by(group) %>%
      do(
        im = mice::mice.impute.lda(
          as.matrix(.[, posit_y]),
          !.$index_NA,
          as.matrix(.[, posit_x])
        )
      ) %>%
      unnest(im) %>%
      ungroup() %>%
      select(im) %>%
      unlist()
  },
  DT_mice = {
    pred_datatable <- data_disc_NA_sort_DT[,
      {
        im <- mice::mice.impute.lda(
          as.matrix(.SD[["y"]]),
          !index_NA,
          as.matrix(.SD[, posit_x, with = F])
        )
      },
      by = .(group)
    ]
  },
  DT_miceFast = {
    pred_datatable_miceFast <- data_disc_NA_sort_DT[,
      {
        im <- fill_NA(as.matrix(.SD), "lda", posit_y, posit_x)
      },
      by = .(group)
    ]
  },
  R_mice = {
    pred_all <- NULL
    for (i in unique(data_disc_NA_sort[, nr_var])) {
      sub <- data_disc_NA_sort[, posit_grs] == i
      temp <- data_disc_NA_sort[sub, ]
      pred <- mice::mice.impute.lda(
        temp[, posit_y],
        !temp[, posit_NA],
        temp[, posit_x]
      )
      pred_all <- c(pred_all, as.numeric(pred))
    }
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    model$set_g(g)
    pred_miceFast <- model$impute("lda", posit_y, posit_x)
    rm(model)
  },
  times = iters
)

m2

g2 <- autoplot(m2, log = FALSE) +
  theme_economist() +
  ggtitle("LDA discrete - with grouping")

ggsave(sprintf("%s/inst/extdata/images/g2.png", base_path), g2)


####################### Binom

pred_lda_bin_mice <- mice::mice.impute.lda(
  data_bin[, posit_y],
  !index_NA,
  data_bin[, posit_x]
)

data <- data_bin_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data)
pred_miceFast <- model$impute("lda", posit_y, posit_x)
rm(model)

table(pred_miceFast$imputations[index_NA], data_bin[index_NA, posit_y])
table(pred_lda_bin_mice, data_bin[index_NA, posit_y])

m3 <- microbenchmark::microbenchmark(
  R_mice = {
    mice::mice.impute.lda(
      data_bin[, posit_y],
      !index_NA,
      data_bin[, posit_x]
    )
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pred_miceFast <- model$impute("lda", posit_y, posit_x)
    rm(model)
  },
  times = iters
)

m3

g3 <- autoplot(m3, log = FALSE) +
  theme_economist() +
  ggtitle("LDA binom - without grouping")

ggsave(sprintf("%s/inst/extdata/images/g3.png", base_path), g3)


##################### Continous - LM Noise

pred_noise_mice <- mice::mice.impute.norm.nob(
  data_con[, posit_y],
  !index_NA,
  data_con[, posit_x]
)

data <- data_con_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data)
pred_miceFast <- model$impute("lm_noise", posit_y, posit_x)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2)
sum((pred_noise_mice - data_con[index_NA, posit_y])^2)

m4 <- microbenchmark::microbenchmark(
  R_mice = {
    mice::mice.impute.norm.nob(
      data_con[, posit_y],
      !index_NA,
      data_con[, posit_x]
    )
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pred_miceFast <- model$impute("lm_noise", posit_y, posit_x)
    rm(model)
  },
  times = iters
)
m4

g4 <- autoplot(m4, log = FALSE) +
  theme_economist() +
  ggtitle("linear regression noise - without grouping")

ggsave(sprintf("%s/inst/extdata/images/g4.png", base_path), g4)

##################### Continous - LM Bayes

pred_bayes_mice <- mice::mice.impute.norm(
  data_con[, posit_y],
  !index_NA,
  data_con[, posit_x]
)

data <- data_con_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data)
pred_miceFast <- model$impute("lm_bayes", posit_y, posit_x)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2)
sum((pred_bayes_mice - data_con[index_NA, posit_y])^2)

m5 <- microbenchmark::microbenchmark(
  R_mice = mice::mice.impute.norm(
    data_con[, posit_y],
    !index_NA,
    data_con[, posit_x]
  ),
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pred_miceFast <- model$impute("lm_bayes", posit_y, posit_x)
    rm(model)
  },
  times = iters
)
m5

g5 <- autoplot(m5, log = FALSE) +
  theme_economist() +
  ggtitle("linear regression bayes - without grouping")

ggsave(sprintf("%s/inst/extdata/images/g5.png", base_path), g5)

##################### Continous - LM Predict

pred_lmpred_mice <- mice::mice.impute.norm.predict(
  data_con[, posit_y],
  !index_NA,
  data_con[, posit_x]
)

data <- cbind(data_con_NA[, c(posit_y, posit_x)], 1)

model <- new(miceFast)
model$set_data(data)
pred_miceFast <- model$impute("lm_pred", posit_y, c(posit_x, max(posit_x) + 1))
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2)
sum((pred_lmpred_mice - data_con[index_NA, posit_y])^2)

m6 <- microbenchmark::microbenchmark(
  R_mice = {
    mice::mice.impute.norm.predict(
      data_con[, posit_y],
      !index_NA,
      data_con[, posit_x]
    )
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pred_miceFast <- model$impute("lm_pred", posit_y, posit_x)
    rm(model)
  },
  times = iters
)

m6

g6 <- autoplot(m6, log = FALSE) +
  theme_economist() +
  ggtitle("linear regression predict - without grouping")

ggsave(sprintf("%s/inst/extdata/images/g6.png", base_path), g6)


## grouping variable

index_sort <- sort(data_con_NA[, posit_grs], index.return = TRUE)$ix

data_con_NA_sort <- data_con_NA[index_sort, ]

pred_Rbase <- NULL
for (i in unique(data_con_NA_sort[, posit_grs])) {
  sub <- data_con_NA_sort[, posit_grs] == i
  temp <- data_con_NA_sort[sub, ]
  pred <- mice::mice.impute.norm.predict(
    as.matrix(temp[, posit_y]),
    !temp[, posit_NA],
    as.matrix(temp[, posit_x])
  )
  pred_Rbase <- c(pred_Rbase, pred)
}

pred_dplyr <- data_con_NA_sort %>%
  as.data.frame() %>%
  group_by(group) %>%
  do(
    im = mice::mice.impute.norm.predict(
      as.matrix(.[, posit_y]),
      !.$index_NA,
      as.matrix(.[, posit_x])
    )
  ) %>%
  unnest(im) %>%
  ungroup() %>%
  select(im) %>%
  unlist() %>%
  as.numeric()

data_con_NA_sort_DT <- data.table(data_con_NA_sort)

pred_datatable <- data_con_NA_sort_DT[,
  {
    im <- mice::mice.impute.norm.predict(
      as.matrix(.SD[, 1]),
      !index_NA,
      as.matrix(.SD[, posit_x, with = F])
    )
  },
  by = .(group)
]

pred_datatable_miceFast <- data_con_NA_sort_DT[,
  {
    im <- fill_NA(cbind(as.matrix(.SD), 1), "lm_pred", 1, posit_x)
  },
  by = .(group)
]

data <- cbind(data_con_NA_sort[, c(posit_y, posit_x)], 1)
g <- data_con_NA_sort[, posit_grs]

model <- new(miceFast)
model$set_data(data)
model$set_g(g)
pred_miceFast <- model$impute("lm_pred", posit_y, posit_x)
rm(model)

true_y <- data_con[index_sort, ][index_NA[index_sort], posit_y]

sum(
  (pred_miceFast$imputations[as.logical(pred_miceFast$index_imputed)] -
    true_y)^2
) # no intercept
sum((pred_dplyr - true_y)^2)
sum((pred_Rbase - true_y)^2)
sum((pred_datatable$V1 - true_y)^2)
sum((pred_datatable_miceFast$V1[index_NA[index_sort]] - true_y)^2) # no intercept

## Performance

m7 <- microbenchmark::microbenchmark(
  dplyr_mice = {
    pred_dplyr <- data_con_NA_sort %>%
      as.data.frame() %>%
      group_by(group) %>%
      do(
        im = mice::mice.impute.norm.predict(
          as.matrix(.[, posit_y]),
          !.$index_NA,
          as.matrix(.[, posit_x])
        )
      ) %>%
      unnest(im) %>%
      ungroup() %>%
      select(im) %>%
      unlist() %>%
      as.numeric()
  },
  R_mice = {
    pred_Rbase <- NULL
    for (i in unique(data_con_NA_sort[, posit_grs])) {
      sub <- data_con_NA_sort[, posit_grs] == i
      temp <- data_con_NA_sort[sub, ]
      pred <- mice::mice.impute.norm.predict(
        as.matrix(temp[, posit_y]),
        !temp[, posit_NA],
        as.matrix(temp[, posit_x])
      )
      pred_Rbase <- c(pred_Rbase, pred)
    }
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    model$set_g(g)
    pred_miceFast <- model$impute("lm_pred", posit_y, posit_x)
    rm(model)
  },
  times = iters,
  DT_mice = {
    pred_datatable <- data_con_NA_sort_DT[,
      {
        im <- mice::mice.impute.norm.predict(
          as.matrix(.SD[, 1]),
          !index_NA,
          as.matrix(.SD[, posit_x, with = F])
        )
      },
      by = .(group)
    ]
  },
  DT_miceFast = {
    pred_datatable_DT <- data_con_NA_sort_DT[,
      {
        im <- fill_NA(as.matrix(.SD), "lm_pred", posit_y, posit_x)
      },
      by = .(group)
    ]
  }
)

m7

g7 <- autoplot(m7, log = FALSE) +
  theme_economist() +
  ggtitle("linear regression predict - with grouping")

ggsave(sprintf("%s/inst/extdata/images/g7.png", base_path), g7)

####
#### Multiple Imputations
####

mice_multi_noise <- rowMeans(sapply(
  1:10,
  function(x) {
    mice::mice.impute.norm.nob(
      data_con[, posit_y],
      !index_NA,
      data_con[, posit_x]
    )
  }
))

data <- data_con_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data)
pred_miceFast <- model$impute_N("lm_noise", posit_y, posit_x, 10)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2)
sum((mice_multi_noise - data_con[index_NA, posit_y])^2)

m8 <- microbenchmark::microbenchmark(
  R_mice = {
    rowMeans(sapply(
      1:10,
      function(x) {
        mice::mice.impute.norm.nob(
          data_con[, posit_y],
          !index_NA,
          data_con[, posit_x]
        )
      }
    ))
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pred_miceFast <- model$impute_N("lm_noise", posit_y, posit_x, 10)
    rm(model)
  },
  times = iters
)
m8

g8 <- autoplot(m8, log = FALSE) +
  theme_economist() +
  ggtitle("linear regression noise - without grouping - multiple 10")

ggsave(sprintf("%s/inst/extdata/images/g8.png", base_path), g8)

##################### Continous - LM Bayes - multiple

mice_multi_bayes <- rowMeans(sapply(
  1:10,
  function(x) {
    mice::mice.impute.norm(data_con[, posit_y], !index_NA, data_con[, posit_x])
  }
))

data <- data_con_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data)
pred_miceFast <- model$impute_N("lm_bayes", posit_y, posit_x, 10)
rm(model)

sum((pred_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2)
sum((mice_multi_bayes - data_con[index_NA, posit_y])^2)

m9 <- microbenchmark::microbenchmark(
  R_mice = {
    rowMeans(sapply(
      1:10,
      function(x) {
        mice::mice.impute.norm(
          data_con[, posit_y],
          !index_NA,
          data_con[, posit_x]
        )
      }
    ))
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pred_miceFast <- model$impute_N("lm_bayes", posit_y, posit_x, 10)
    rm(model)
  },
  times = iters
)
m9

g9 <- autoplot(m9, log = FALSE) +
  theme_economist() +
  ggtitle("linear regression bayes - without grouping - multiple 10")

ggsave(sprintf("%s/inst/extdata/images/g9.png", base_path), g9)

##################### Continous - VIFS

data.vifs <- data.frame(data_con[, c(posit_y, posit_x)])
colnames(data.vifs) <- letters[1:(nr_var - 2)]
vifs_car <- car::vif(lm(a ~ ., data = data.vifs))

model <- new(miceFast)
model$set_data(data_con)
vifs_miceFast <- model$vifs(posit_y, posit_x)
rm(model)

vifs_car
vifs_miceFast

m10 <- microbenchmark::microbenchmark(
  car = {
    car::vif(lm(a ~ ., data = data.vifs))
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data_con)
    vifs_miceFast <- model$vifs(posit_y, posit_x)
    rm(model)
  },
  times = iters
)
m10

g10 <- autoplot(m10, log = FALSE) + theme_economist() + ggtitle("vifs")

ggsave(sprintf("%s/inst/extdata/images/g10.png", base_path), g10)

##################### PMM

data <- data_con_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data)
pmm_miceFast <- model$impute_N("pmm", posit_y, posit_x, k = 5)
rm(model)

pred_pmm_mice <- mice::mice.impute.pmm(
  data_con[, posit_y],
  !index_NA,
  data_con[, posit_x]
)

sum((pmm_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2)
sum((pred_pmm_mice - data_con[index_NA, posit_y])^2)

m11 <- microbenchmark::microbenchmark(
  R_mice = {
    mice::mice.impute.pmm(data_con[, posit_y], !index_NA, data_con[, posit_x])
  },
  miceFast = {
    model <- new(miceFast)
    model$set_data(data)
    pmm_miceFast <- model$impute_N("pmm", posit_y, posit_x, k = 5)
    rm(model)
  },
  times = iters
)
m11

g11 <- autoplot(m11, log = FALSE) + theme_economist() + ggtitle("pmm")

ggsave(sprintf("%s/inst/extdata/images/g11.png", base_path), g11)

##################### Weighted LM Predict

data_w <- data_con_NA[, c(posit_y, posit_x)]

model <- new(miceFast)
model$set_data(data_w)
model$set_w(w_c)
pred_w_miceFast <- model$impute("lm_pred", posit_y, posit_x)
rm(model)

cat(
  "Weighted lm_pred SSE:",
  sum((pred_w_miceFast$imputations[index_NA] - data_con[index_NA, posit_y])^2),
  "\n"
)

m12 <- microbenchmark::microbenchmark(
  miceFast_weighted = {
    model <- new(miceFast)
    model$set_data(data_w)
    model$set_w(w_c)
    model$impute("lm_pred", posit_y, posit_x)
    rm(model)
  },
  miceFast_unweighted = {
    model <- new(miceFast)
    model$set_data(data_w)
    model$impute("lm_pred", posit_y, posit_x)
    rm(model)
  },
  times = iters
)
m12

g12 <- autoplot(m12, log = FALSE) +
  theme_economist() +
  ggtitle("weighted vs unweighted lm_pred")

ggsave(sprintf("%s/inst/extdata/images/g12.png", base_path), g12)

##################### fill_NA / fill_NA_N functional interface

data_fill <- data_con_NA[, c(posit_y, posit_x)]

res_fill <- fill_NA(data_fill, "lm_pred", posit_y, posit_x)
res_fill_N <- fill_NA_N(data_fill, "lm_noise", posit_y, posit_x, k = 10)

cat(
  "fill_NA  lm_pred SSE:",
  sum((res_fill[index_NA] - data_con[index_NA, posit_y])^2),
  "\n"
)
cat(
  "fill_NA_N lm_noise (k=10) SSE:",
  sum((res_fill_N[index_NA] - data_con[index_NA, posit_y])^2),
  "\n"
)

m13 <- microbenchmark::microbenchmark(
  fill_NA_lm_pred = fill_NA(data_fill, "lm_pred", posit_y, posit_x),
  fill_NA_N_lm_noise = fill_NA_N(
    data_fill,
    "lm_noise",
    posit_y,
    posit_x,
    k = 10
  ),
  fill_NA_N_pmm = fill_NA_N(data_fill, "pmm", posit_y, posit_x, k = 5),
  times = iters
)
m13

g13 <- autoplot(m13, log = FALSE) +
  theme_economist() +
  ggtitle("fill_NA / fill_NA_N functional interface")

ggsave(sprintf("%s/inst/extdata/images/g13.png", base_path), g13)

# plot for README/Intro

dats <- bind_rows(list(
  data.frame(m1) %>% mutate(model = "LDA discrete - without grouping"),
  data.frame(m2) %>% mutate(model = "LDA discrete - with grouping"),
  data.frame(m3) %>% mutate(model = "LDA binom - without grouping"),
  data.frame(m4) %>%
    mutate(model = "linear regression noise - without grouping"),
  data.frame(m5) %>%
    mutate(model = "linear regression bayes - without grouping"),
  data.frame(m6) %>%
    mutate(model = "linear regression predict - without grouping"),
  data.frame(m7) %>%
    mutate(model = "linear regression predict - with grouping"),
  data.frame(m8) %>%
    mutate(model = "linear regression noise - without grouping - multiple 10"),
  data.frame(m9) %>%
    mutate(model = "linear regression bayes - without grouping - multiple 10"),
  data.frame(m10) %>% mutate(model = "VIF"),
  data.frame(m11) %>% mutate(model = "pmm"),
  data.frame(m12) %>% mutate(model = "weighted vs unweighted lm_pred"),
  data.frame(m13) %>% mutate(model = "fill_NA / fill_NA_N")
))

dats_plot <- dats %>%
  group_by(model, expr) %>%
  summarise(mean_time_sec = mean(time / 10**9), .groups = "drop") %>%
  group_by(model) %>%
  mutate(relative_time = mean_time_sec / min(mean_time_sec)) %>%
  rename("package" = "expr")

g_summary <- ggplot(dats_plot, aes(model, relative_time, fill = package)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle(paste0(
    "Benchmarks - 10^",
    power,
    "obs (20% NA) ",
    nr_var,
    "vars and ",
    grs,
    " groups"
  ))

ggsave(sprintf("%s/inst/extdata/images/g_summary.png", base_path), g_summary)
