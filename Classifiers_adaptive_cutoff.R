set.seed(1)

library(keras)
library(sn)
library(abind)
library(nortest)
library(tseries)
library(dplyr)
library(moments)

N <- 60 # Two dimensional unit square regular grid of size 60×60
nu <- 0.5 # Smoothness parameter: 0.5 or 1.0 

load(paste("Train-data-beta30-p4-nu", nu, "-gridsize", N, ".Rdata", sep=""))
load(paste("Test-data-beta50-p10-nu", nu, "-gridsize", N, ".Rdata", sep=""))

train.data <- abind(train.normal, train.nonnormal)
test.data <- abind(test.normal, test.nonnormal)

y.train <- c(rep(0, dim(train.normal)[3]), rep(1, dim(train.nonnormal)[3]))
y.test <- c(rep(0, dim(test.normal)[3]), rep(1, dim(test.nonnormal)[3]))

### Individual normality tests
descriptor <- function(x){
  return(c(shapiro.test(x)$stat, lillie.test(x)$stat, jarque.bera.test(x)$stat, ad.test(x)$stat,
           skewness(x), kurtosis(x)))
}

desc_f <- function(x){
  x <- as.vector(x)
  return(c(descriptor(x)))
}

train.desc <- t(apply(train.data, 3, desc_f))
test.normal.desc <- t(apply(test.normal, 3, desc_f))
test.nonnormal.desc <- t(apply(test.nonnormal, 3, desc_f))


### Learning neural networks and linear classifiers
batch_size <- 128
epochs <- 500

# Neural Networks
model <- keras_model_sequential() %>%
  layer_dense(units = 256, activation = 'relu', input_shape=c(6)) %>%
  layer_dropout(rate = 0.3) %>% 
  layer_batch_normalization() %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>% 
  layer_batch_normalization() %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = .0001),
  metrics = c('accuracy')
)

model %>% fit(
  train.desc, y.train,
  batch_size = batch_size,
  epochs = epochs#,
)

# Linear
model.linear <- keras_model_sequential() %>%
  layer_dense(units = 1,activation = 'sigmoid',  input_shape=c(6)) 

model.linear %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(learning_rate = .0001),
  metrics = c('accuracy')
)

model.linear %>% fit(
  train.desc, y.train,
  batch_size = batch_size,
  epochs = epochs#,
)

### Cutoffs from Neural Networks
train.preds <- model %>% predict(train.desc)
train.preds <- train.preds[1:dim(train.normal)[3]]
cutoff <- quantile(train.preds, .95)

# Adaptive cutoffs
nn <- length(train.preds) / length(betas)
cutoffs <- c()
for(k in 1:length(betas)){
  cutoffs[k] <- quantile(train.preds[((k-1)*nn+1):(k*nn)], 0.95)  
}

test.stats.normal <- model %>% predict(test.normal.desc)
test.stats.nonnormal <- model %>% predict(test.nonnormal.desc)


### Cutoffs from Linear classifier
train.preds.linear <- model.linear %>% predict(train.desc)
train.preds.linear <- train.preds.linear[1:dim(train.normal)[3]]
cutoff.linear <- quantile(train.preds.linear, .95)

# Adaptive cutoffs
cutoffs.linear <- c()
for(k in 1:length(betas)){
  cutoffs.linear[k] <- quantile(train.preds.linear[((k-1)*nn+1):(k*nn)], 0.95)  
}


test.stats.normal.linear <- model.linear %>% predict(test.normal.desc)
test.stats.nonnormal.linear <- model.linear %>% predict(test.nonnormal.desc)


save(train.desc, test.normal.desc, test.nonnormal.desc,
     test.stats.normal, test.stats.nonnormal,
     test.stats.normal.linear, test.stats.nonnormal.linear,
     cutoff, cutoff.linear,
     cutoffs, cutoffs.linear,
     file = paste("NN-Linear-train.beta", length(betas), "-test.beta", length(betas.test),"-p",length(powers.test),
                  "-nu", nu, "-gridsize", N, ".Rdata", sep=""))
