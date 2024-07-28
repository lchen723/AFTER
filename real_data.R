library(penalized)
library(tidyverse)
library(ggpubr)
W <- function(r_star,g_x){
  a <- data.frame(r_star); b <- data.frame(g_x)
  data <- cbind(a,b)
  colnames(data) <- c('Y',"X")
  w <- lm(Y~X+0,data = data)
  return(w$coefficients)
}
interval<- function(lower_bound,upper_bound){         
  interval_size = (upper_bound-lower_bound)/7
  interval_points <- c(lower_bound)
  for(i in (1:7)){
    interval_points <- c(interval_points,lower_bound+interval_size*i)
  }
  return(interval_points)
}

# import data
data(nki70)
real_data <- nki70
data <- real_data %>% select(c(1:2),c(8:77))
data$time <- log(data$time)


# high data classifications 0.6, high variance 0.75
# generate v, which follows which V follows N(0,0.75), and constant = 1
set.seed(1)
v <- data.frame(rnorm(dim(data)[1],0,0.75))
constant <- data.frame(rep(1,dim(data)[1]))


# ME correction
# correction for survival time(y)
# assume y_with_measurement = y + 1 + 1*TSPYLS + V, which V follows N(0,0.75)
# correct y: y_hat = y_with_measurement_error-E(1 + 1*1*TSPYLS + V) = y_with_measurement_error-E(w)
# use sample mean to E(x) and V(x)
estimated_w <- 1+ 1*mean(data[,3])+mean(v[,1])
correction_last_part <- data[,3]-mean(data[,3])
y_hat <- data.frame(data[,1] - estimated_w-correction_last_part)
colnames(y_hat) <- c("y_hat")


# correction for censoring indicator
#correct indicator:p(censoring indicator hat = 1) = (p(indicator_with_error = 1) - 0.6) / 1 - 0.6 - 0.6
#indicator_hat_probability > 0, indicator_hat = 1, otherwise = 0
indicator_hat_probability <- (data[,2]-0.6)/(1-0.6-0.6)
indicator_hat <- NULL
for (i in c(1:length(indicator_hat_probability))){
  if (indicator_hat_probability[i]<0){
    indicator_hat[i] <- 0
  }else{
    indicator_hat[i] <- 1
  }
}


#boosting for corrected data
#corrected data
corrected_data<- cbind(y_hat,indicator_hat,data[,3:72])
corrected_data<- data
colnames(corrected_data)[1:2] <- c('Y','censoring_indicator')


# number of sample
n = dim(corrected_data)[1]
# dimension
p = 70


censoring_indicator <- corrected_data$censoring_indicator
variable_catch <- c()
y <- corrected_data$Y


#step 0.
sum_of_every_fitted_value <- rep(0,times=n)
y_star <- corrected_data$Y
r <- y - sum_of_every_fitted_value
r_star <- y_star - sum_of_every_fitted_value
df3 <- data.frame(corrected_data$censoring_indicator,r_star,r)
colnames(df3)[1] <- "censoring_indicator"
survival_probability <- c()

for (i in c(1:dim(df3)[1])){
  yy <- df3[df3$r<=df3$r[i],]
  u <- c(yy$r)
  data_producted <- 1
  for (j in c(1:length(u))){
    total_sum_of_denominator <- sum((df3$r>=u[j])*1)
    censoring_indicator_1<- df3[df3$censoring_indicator==1,]
    total_sum_of_numerator <- sum((censoring_indicator_1$r==u[j])*1)
    data_producted = data_producted * (1-total_sum_of_numerator/total_sum_of_denominator)
  }
  survival_probability[i] <- data_producted
}


survival_probability[which(survival_probability==0)] <- 0.000001


sum_riemann_for_every_n <- c()
upper = max(r)
for (i in c(1:length(df3$r))){
  point <- interval(df3$r[i],upper)
  riemann_x <- c()
  data_producted <- 1
  for (j in c(1:length(point))){
    xx <- df3[df3$r<=point[j],]
    u <- c(xx$r)
    for (k in c(1:length(u))){
      total_sum_of_denominator <- sum((df3$r>=u[k])*1)
      censoring_indicator_1<- df3[df3$censoring_indicator==1,]
      total_sum_of_numerator <- sum((censoring_indicator_1$r==u[k])*1)
      data_producted = data_producted *  (1-total_sum_of_numerator/total_sum_of_denominator)
    }
    riemann_x[j]<- data_producted
  }

  riemann_minus <- c()
  for (l in c(1:length(riemann_x)-1)){
    riemann_minus[l] <- riemann_x[l+1]-riemann_x[l]
  }


  sum_integrate <-c(0)
  for (m in c(1:length(riemann_minus))){
    sum_integrate <- sum_integrate + point[m+1] * riemann_minus[m]
  }
  sum_riemann_for_every_n[i] <- sum_integrate
}


y_star_update <- c()
for (i in c(1:length(y))){
  y_star_update[i] <- censoring_indicator[i] * y[i]+
    (1-censoring_indicator[i])*(sum_of_every_fitted_value[i] - sum_riemann_for_every_n[i]/survival_probability[i])
}


y_star <- y_star_update


add_g_every_time <- data.frame()
df <- as.data.frame(matrix(numeric(0),ncol = p, nrow = n))
df[is.na(df)] <- 0
for (i in c(1:p)){
  colnames(df)[i] <- i
}


iterations = 0
stop_times = 50


while (TRUE) {
  r_star <- y_star - sum_of_every_fitted_value
  residual<- c(0)
  add_f <- c()
  times <- c(0)
  variable <- c()


  for (i in c(1:p)){
    times = times + 1
    f_variables <- smooth.spline(x=corrected_data[,i+2],y=r_star,cv=FALSE,all.knots=c(0,0.2,0.4,0.6,0.8,1))
    if (residual==0){
      add_f <- f_variables
      residual <- f_variables$pen.crit
      variable <- times
    } else if(f_variables$pen.crit < residual){
      add_f <- f_variables
      variable <- times
      residual <- f_variables$pen.crit
    }
  }


  variable_catch <- c(variable_catch,variable)
  fit = fitted(add_f)
  fit[which(is.na(fit))]=0

  w <- W(r_star,fit)
  add_g <- as.data.frame(w * fit, ncol = 1)

  stop_value <- 1e-3
  number_of_added_value <- sum((abs(add_g) < stop_value)*1)

  if (iterations == stop_times){
    break
  }

  if (number_of_added_value != n && iterations!= stop_times){
    df[as.character(variable)] <- df[,variable] + add_g

    sum_of_every_fitted_value <- sum_of_every_fitted_value + w * fit
    r_star = y_star - sum_of_every_fitted_value
    r <- y - sum_of_every_fitted_value
    df1 <- data.frame(corrected_data$censoring_indicator,r_star,r)
    colnames(df1)[1] <- "censoring_indicator"
    survival_probability <- c()

    for (i in c(1:dim(df1)[1])){
      yy <- df1[df1$r<=df1$r[i],]
      u <- c(yy$r)
      data_producted <- 1
      for (j in c(1:length(u))){
        total_sum_of_denominator <- sum((df1$r>=u[j])*1)
        censoring_indicator_1<- df1[df1$censoring_indicator==1,]
        total_sum_of_numerator <- sum((censoring_indicator_1$r==u[j])*1)
        data_producted = data_producted * (1-total_sum_of_numerator/total_sum_of_denominator)
      }
      survival_probability[i] <- data_producted
    }

    survival_probability[which(survival_probability==0)] <- 0.000001

    sum_riemann_for_every_n <- c()
    upper = max(r)
    for (i in c(1:length(df1$r))){
      point <- interval(df1$r[i],upper)
      riemann_x <- c()
      data_producted <- 1
      for (j in c(1:length(point))){
        xx <- df1[df1$r<=point[j],]
        u <- c(xx$r)
        for (k in c(1:length(u))){
          total_sum_of_denominator <- sum((df1$r>=u[k])*1)
          total_sum_of_denominator
          censoring_indicator_1<- df1[df1$censoring_indicator==1,]
          total_sum_of_numerator <- sum((censoring_indicator_1$r==u[k])*1)
          data_producted = data_producted *  (1-total_sum_of_numerator/total_sum_of_denominator)
        }
        riemann_x[j]<- data_producted
      }

      riemann_minus <- c()
      for (l in c(1:length(riemann_x)-1)){
        riemann_minus[l] <- riemann_x[l+1]-riemann_x[l]
      }


      sum_integrate <-c(0)
      for (m in c(1:length(riemann_minus))){
        sum_integrate <- sum_integrate + point[m+1] * riemann_minus[m]
      }
      sum_riemann_for_every_n[i] <- sum_integrate
    }


    y_star_update <- c()

    for (i in c(1:length(y))){
      y_star_update[i] <- censoring_indicator[i] * y[i]+
        (1-censoring_indicator[i])*(sum_of_every_fitted_value[i] - sum_riemann_for_every_n[i]/survival_probability[i])
    }

    y_star <- y_star_update

    iterations = iterations + 1

  }else{
    break
  }
}


a <- table(variable_catch)
a[a>=2]  
aa <- as.numeric(names(a[a>=2]))




for (i in 1:length(aa)){
  name <- paste0(names(a[a>=2])[i],'.jpeg')
  jpeg(name)
  x_1 <- as.data.frame(cbind(corrected_data[,aa[i]+2],df[,aa[i]]+0.0))
  colnames(x_1) <- c('x','y')
  x_1 <- x_1[order(x_1$x),]
  plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[aa[i]+2],type='l', lwd=2,lty=1)
  dev.off()
  
  
}


# 2
x_1 <- as.data.frame(cbind(corrected_data[,4],df[,2]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[4],type='l', lwd=2,lty=1)
dev.off()


jpeg("h2.jpeg")
# 10
x_1 <- as.data.frame(cbind(corrected_data[,12],df[,10]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[12],type='l', lwd=2,lty=1)
dev.off()


jpeg("h3.jpeg")
# 13
x_1 <- as.data.frame(cbind(corrected_data[,15],df[,13]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y,col='black', lwd=2,lty=1,type='l',xlab="x",ylab='f', main= colnames(corrected_data)[15])
dev.off()


jpeg("h4.jpeg")
# 19
x_1 <- as.data.frame(cbind(corrected_data[,21],df[,19]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y,col='black', lwd=2,lty=1,type='l',xlab="x",ylab='f', main= colnames(corrected_data)[21])
dev.off()


jpeg("h5.jpeg")
# 33
x_1 <- as.data.frame(cbind(corrected_data[,35],df[,33]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y,col='black', lwd=2,lty=1,type='l',xlab="x",ylab='f', main= colnames(corrected_data)[35])
dev.off()


jpeg("h6.jpeg")
# 35
x_1 <- as.data.frame(cbind(corrected_data[,37],df[,35]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y,col='black', lwd=2,lty=1,type='l', xlab="x",ylab='f', main= colnames(corrected_data)[37])
dev.off()


jpeg("h7.jpeg")
# 47
x_1 <- as.data.frame(cbind(corrected_data[,41],df[,39]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[41],type='l', lwd=2,lty=1)
dev.off()


jpeg("h8.jpeg")
# 52
x_1 <- as.data.frame(cbind(corrected_data[,54],df[,52]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[54],type='l', lwd=2,lty=1)
dev.off()


jpeg("h9.jpeg")
# 53
x_1 <- as.data.frame(cbind(corrected_data[,55],df[,53]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[55],type='l', lwd=2,lty=1)
dev.off()


jpeg("h10.jpeg")
# 62
x_1 <- as.data.frame(cbind(corrected_data[,64],df[,62]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[64],type='l', lwd=2,lty=1)
dev.off()


jpeg("h11.jpeg")
# 64
x_1 <- as.data.frame(cbind(corrected_data[,66],df[,64]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[66],type='l', lwd=2,lty=1)
dev.off()


jpeg("h12.jpeg")
# 67
x_1 <- as.data.frame(cbind(corrected_data[,69],df[,67]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[69],type='l', lwd=2,lty=1)
dev.off()


survival_data_a <- df
sur_time <-apply(survival_data_a,1,sum)
sur_time <-exp(sur_time)
sur_time <-sort(sur_time)
sur_time <-data.frame(sur_time)
sur_time
pro <- seq(1,144)/144
pro <- sort(pro, decreasing= TRUE)
pro <- data.frame(pro)


# survival probability
ddf1 <- cbind(sur_time, pro)
plot(ddf1$sur_time,ddf1$pro, xlab="time",ylab='survival probability', main='survival curve',type='l', lwd=2,lty=1)


# import data
data(nki70)
real_data <- nki70
data <- real_data %>% select(c(1:2),c(8:77))
dim(data)[1]
data$time <- log(data$time)
# low data classifications 0.1, high variance 0.15
# generate v, which follows which V follows N(0,0.15), and constant = 1
set.seed(1)
v <- data.frame(rnorm(dim(data)[1],0,0.15))
constant <- data.frame(rep(1,dim(data)[1]))

# ME correction
# correction for survival time(y)
# assume y_with_measurement = y + 1 + 1*TSPYLS + V, which V follows N(0,0.15)
# correct y: y_hat = y_with_measurement_error-E(1 + 1*1*TSPYLS + V) = y_with_measurement_error-E(w)
estimated_w <- 1+ 1*mean(data[,3])+mean(v[,1])
correction_last_part <- data[,3]-mean(data[,3])
y_hat <- data.frame(data[,1] - estimated_w-correction_last_part)
colnames(y_hat) <- c("y_hat")



# correction for censoring indicator
#correct indicator:p(censoring indicator hat = 1) = (p(indicator_with_error = 1) - 0.1) / 1 - 0.1- 0.1
#indicator_hat_probability > 0, indicator_hat = 1, otherwise = 0
indicator_hat_probability <- (data[,2]-0.1)/(1-0.1-0.1)
indicator_hat <- NULL
for (i in c(1:length(indicator_hat_probability))){
  if (indicator_hat_probability[i]<0){
    indicator_hat[i] <- 0
  }else{
    indicator_hat[i] <- 1
  }
}



#boosting for corrected data
#corrected data
corrected_data<- cbind(y_hat,indicator_hat,data[,3:72])
colnames(corrected_data)[1:2] <- c('Y','censoring_indicator')



# number of sample
n = dim(corrected_data)[1]
# dimension
p = 70


censoring_indicator <- corrected_data$censoring_indicator
variable_catch <- c()
y <- corrected_data$Y

#step 0.
sum_of_every_fitted_value <- rep(0,times=n)
y_star <- corrected_data$Y
r <- y - sum_of_every_fitted_value
r_star <- y_star - sum_of_every_fitted_value
df3 <- data.frame(corrected_data$censoring_indicator,r_star,r)
colnames(df3)[1] <- "censoring_indicator"
survival_probability <- c()

for (i in c(1:dim(df3)[1])){
  yy <- df3[df3$r<=df3$r[i],]
  u <- c(yy$r)
  data_producted <- 1
  for (j in c(1:length(u))){
    total_sum_of_denominator <- sum((df3$r>=u[j])*1)
    censoring_indicator_1<- df3[df3$censoring_indicator==1,]
    total_sum_of_numerator <- sum((censoring_indicator_1$r==u[j])*1)
    data_producted = data_producted * (1-total_sum_of_numerator/total_sum_of_denominator)
  }
  survival_probability[i] <- data_producted
}


survival_probability[which(survival_probability==0)] <- 0.000001


sum_riemann_for_every_n <- c()
upper = max(r)
for (i in c(1:length(df3$r))){
  point <- interval(df3$r[i],upper)
  riemann_x <- c()
  data_producted <- 1
  for (j in c(1:length(point))){
    xx <- df3[df3$r<=point[j],]
    u <- c(xx$r)
    for (k in c(1:length(u))){
      total_sum_of_denominator <- sum((df3$r>=u[k])*1)
      censoring_indicator_1<- df3[df3$censoring_indicator==1,]
      total_sum_of_numerator <- sum((censoring_indicator_1$r==u[k])*1)
      data_producted = data_producted *  (1-total_sum_of_numerator/total_sum_of_denominator)
    }
    riemann_x[j]<- data_producted
  }

  riemann_minus <- c()
  for (l in c(1:length(riemann_x)-1)){
    riemann_minus[l] <- riemann_x[l+1]-riemann_x[l]
  }


  sum_integrate <-c(0)
  for (m in c(1:length(riemann_minus))){
    sum_integrate <- sum_integrate + point[m+1] * riemann_minus[m]
  }
  sum_riemann_for_every_n[i] <- sum_integrate
}


y_star_update <- c()
for (i in c(1:length(y))){
  y_star_update[i] <- censoring_indicator[i] * y[i]+
    (1-censoring_indicator[i])*(sum_of_every_fitted_value[i] - sum_riemann_for_every_n[i]/survival_probability[i])
}

y_star <- y_star_update

add_g_every_time <- data.frame()
df <- as.data.frame(matrix(numeric(0),ncol = p, nrow = n))
df[is.na(df)] <- 0
for (i in c(1:p)){
  colnames(df)[i] <- i
}


iterations = 0
stop_times = 50


while (TRUE) {
  r_star <- y_star - sum_of_every_fitted_value
  residual<- c(0)
  add_f <- c()
  times <- c(0)
  variable <- c()


  for (i in c(1:p)){
    times = times + 1
    f_variables <- smooth.spline(x=corrected_data[,i+2],y=r_star,cv=FALSE,all.knots=c(0,0.2,0.4,0.6,0.8,1))
    if (residual==0){
      add_f <- f_variables
      residual <- f_variables$pen.crit
      variable <- times
    } else if(f_variables$pen.crit < residual){
      add_f <- f_variables
      variable <- times
      residual <- f_variables$pen.crit
    }
  }


  variable_catch <- c(variable_catch,variable)
  fit = fitted(add_f)
  fit[which(is.na(fit))]=0

  w <- W(r_star,fit)
  add_g <- as.data.frame(w * fit, ncol = 1)


  stop_value <- 1e-3
  number_of_added_value <- sum((abs(add_g) < stop_value)*1)

  if (iterations == stop_times){
    break
  }

  if (number_of_added_value != n && iterations!= stop_times){
    df[as.character(variable)] <- df[,variable] + add_g

    sum_of_every_fitted_value <- sum_of_every_fitted_value + w * fit
    r_star = y_star - sum_of_every_fitted_value
    r <- y - sum_of_every_fitted_value
    df1 <- data.frame(corrected_data$censoring_indicator,r_star,r)
    colnames(df1)[1] <- "censoring_indicator"
    survival_probability <- c()

    for (i in c(1:dim(df1)[1])){
      yy <- df1[df1$r<=df1$r[i],]
      u <- c(yy$r)
      data_producted <- 1
      for (j in c(1:length(u))){
        total_sum_of_denominator <- sum((df1$r>=u[j])*1)
        censoring_indicator_1<- df1[df1$censoring_indicator==1,]
        total_sum_of_numerator <- sum((censoring_indicator_1$r==u[j])*1)
        data_producted = data_producted * (1-total_sum_of_numerator/total_sum_of_denominator)
      }
      survival_probability[i] <- data_producted
    }

    survival_probability[which(survival_probability==0)] <- 0.000001

    sum_riemann_for_every_n <- c()
    upper = max(r)

    for (i in c(1:length(df1$r))){
      point <- interval(df1$r[i],upper)
      riemann_x <- c()
      data_producted <- 1
      for (j in c(1:length(point))){
        xx <- df1[df1$r<=point[j],]
        u <- c(xx$r)
        for (k in c(1:length(u))){
          total_sum_of_denominator <- sum((df1$r>=u[k])*1)
          total_sum_of_denominator
          censoring_indicator_1<- df1[df1$censoring_indicator==1,]
          total_sum_of_numerator <- sum((censoring_indicator_1$r==u[k])*1)
          data_producted = data_producted *  (1-total_sum_of_numerator/total_sum_of_denominator)
        }
        riemann_x[j]<- data_producted
      }

      riemann_minus <- c()
      for (l in c(1:length(riemann_x)-1)){
        riemann_minus[l] <- riemann_x[l+1]-riemann_x[l]
      }


      sum_integrate <-c(0)
      for (m in c(1:length(riemann_minus))){
        sum_integrate <- sum_integrate + point[m+1] * riemann_minus[m]
      }
      sum_riemann_for_every_n[i] <- sum_integrate
    }


    y_star_update <- c()

    for (i in c(1:length(y))){
      y_star_update[i] <- censoring_indicator[i] * y[i]+
        (1-censoring_indicator[i])*(sum_of_every_fitted_value[i] - sum_riemann_for_every_n[i]/survival_probability[i])
    }

    y_star <- y_star_update

    iterations = iterations + 1

  }else{
    break
  }
}




a <- table(variable_catch)
a[a>=2]  # 4, 5, 7, 25, 30, 35, 40, 51, 64, 66, 68
aa <- as.numeric(names(a[a>=2]))



for (i in 1:length(aa)){
  name <- paste0(names(a[a>=2])[i],'.jpeg')
  jpeg(name)
  x_1 <- as.data.frame(cbind(corrected_data[,aa[i]+2],df[,aa[i]]+0.0))
  colnames(x_1) <- c('x','y')
  x_1 <- x_1[order(x_1$x),]
  plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[aa[i]+2],type='l', lwd=2,lty=1)
  dev.off()
  
  
}



jpeg("l1.jpeg")
# 4
x_1 <- as.data.frame(cbind(corrected_data[,6],df[,4]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
c <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[6],type='l', lwd=2,lty=1)
dev.off()


jpeg("l2.jpeg")
# 5
x_1 <- as.data.frame(cbind(corrected_data[,7],df[,5]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
d <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[7],type='l', lwd=2,lty=1)
dev.off()


jpeg("l3.jpeg")
# 7
x_1 <- as.data.frame(cbind(corrected_data[,9],df[,7]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
e <- plot(x_1$x,x_1$y,col='black', lwd=2,lty=1,type='l',xlab="x",ylab='f', main= colnames(corrected_data)[9])
dev.off()


jpeg("l4.jpeg")
# 25
x_1 <- as.data.frame(cbind(corrected_data[,27],df[,25]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
f <- plot(x_1$x,x_1$y,col='black', lwd=2,lty=1,type='l', xlab="x",ylab='f', main= colnames(corrected_data)[27])
dev.off()



jpeg("l5.jpeg")
# 30
x_1 <- as.data.frame(cbind(corrected_data[,32],df[,30]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
g <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[32],type='l', lwd=2,lty=1)
dev.off()



jpeg("l6.jpeg")
# 35
x_1 <- as.data.frame(cbind(corrected_data[,37],df[,35]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
h <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[37],type='l', lwd=2,lty=1)
dev.off()



jpeg("l7.jpeg")
# 51
x_1 <- as.data.frame(cbind(corrected_data[,53],df[,51]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
i <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[53],type='l', lwd=2,lty=1)
dev.off()



jpeg("l8.jpeg")
# 54
x_1 <- as.data.frame(cbind(corrected_data[,56],df[,54]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
j <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[56],type='l', lwd=2,lty=1)
dev.off()



jpeg("l9.jpeg")
# 64
x_1 <- as.data.frame(cbind(corrected_data[,66],df[,64]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
k <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[66],type='l', lwd=2,lty=1)
dev.off()



jpeg("l10.jpeg")
# 66
x_1 <- as.data.frame(cbind(corrected_data[,68],df[,66]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
l <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[68],type='l', lwd=2,lty=1)
dev.off()



jpeg("l11.jpeg")
# 68
x_1 <- as.data.frame(cbind(corrected_data[,70],df[,68]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
m <- plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[70],type='l', lwd=2,lty=1)
dev.off()


survival_data_b <- df
sur_time <- apply(survival_data_b ,1,sum)
sur_time <-exp(sur_time)
sur_time <-sort(sur_time)
sur_time  <- data.frame(sur_time)
sur_time
pro <- seq(1,144)/144
pro <- sort(pro, decreasing= TRUE)
pro <- data.frame(pro)

# survival probability
ddf2 <- cbind(sur_time, pro)
plot(ddf2$sur_time,ddf2$pro, xlab="time",ylab='survival probability', main='survival curve',type='l', lwd=2,lty=1)







#####  naive data
#boosting for naive data
colnames(data)[1:2] <- c('Y','censoring_indicator')
data$time <- log(data$time)
# number of sample
n = dim(data)[1]
# dimension
p = 70


censoring_indicator <- data$censoring_indicator
variable_catch <- c()
y <- data$Y


#step 0.
sum_of_every_fitted_value <- rep(0,times=n)
y_star <- data$Y
r <- y - sum_of_every_fitted_value
r_star <- y_star - sum_of_every_fitted_value
df3 <- data.frame(data$censoring_indicator,r_star,r)
colnames(df3)[1] <- "censoring_indicator"
survival_probability <- c()

for (i in c(1:dim(df3)[1])){
  yy <- df3[df3$r<=df3$r[i],]
  u <- c(yy$r)
  data_producted <- 1
  for (j in c(1:length(u))){
    total_sum_of_denominator <- sum((df3$r>=u[j])*1)
    censoring_indicator_1<- df3[df3$censoring_indicator==1,]
    total_sum_of_numerator <- sum((censoring_indicator_1$r==u[j])*1)
    data_producted = data_producted * (1-total_sum_of_numerator/total_sum_of_denominator)
  }
  survival_probability[i] <- data_producted
}



survival_probability[which(survival_probability==0)] <- 0.000001




sum_riemann_for_every_n <- c()
upper = max(r)
for (i in c(1:length(df3$r))){
  point <- interval(df3$r[i],upper)
  riemann_x <- c()
  data_producted <- 1
  for (j in c(1:length(point))){
    xx <- df3[df3$r<=point[j],]
    u <- c(xx$r)
    for (k in c(1:length(u))){
      total_sum_of_denominator <- sum((df3$r>=u[k])*1)
      censoring_indicator_1<- df3[df3$censoring_indicator==1,]
      total_sum_of_numerator <- sum((censoring_indicator_1$r==u[k])*1)
      data_producted = data_producted *  (1-total_sum_of_numerator/total_sum_of_denominator)
    }
    riemann_x[j]<- data_producted
  }

  riemann_minus <- c()
  for (l in c(1:length(riemann_x)-1)){
    riemann_minus[l] <- riemann_x[l+1]-riemann_x[l]
  }


  sum_integrate <-c(0)
  for (m in c(1:length(riemann_minus))){
    sum_integrate <- sum_integrate + point[m+1] * riemann_minus[m]
  }
  sum_riemann_for_every_n[i] <- sum_integrate
}


y_star_update <- c()
for (i in c(1:length(y))){
  y_star_update[i] <- censoring_indicator[i] * y[i]+
    (1-censoring_indicator[i])*(sum_of_every_fitted_value[i] - sum_riemann_for_every_n[i]/survival_probability[i])
}


y_star <- y_star_update


add_g_every_time <- data.frame()
df <- as.data.frame(matrix(numeric(0),ncol = p, nrow = n))
df[is.na(df)] <- 0
for (i in c(1:p)){
  colnames(df)[i] <- i
}


iterations = 0
stop_times = 50



while (TRUE) {
  r_star <- y_star - sum_of_every_fitted_value
  residual<- c(0)
  add_f <- c()
  times <- c(0)
  variable <- c()


  for (i in c(1:p)){
    times = times + 1
    f_variables <- smooth.spline(x=data[,i+2],y=r_star,cv=FALSE,all.knots=c(0,0.2,0.4,0.6,0.8,1))
    if (residual==0){
      add_f <- f_variables
      residual <- f_variables$pen.crit
      variable <- times
    } else if(f_variables$pen.crit < residual){
      add_f <- f_variables
      variable <- times
      residual <- f_variables$pen.crit
    }
  }


  variable_catch <- c(variable_catch,variable)
  fit = fitted(add_f)
  fit[which(is.na(fit))]=0

  w <- W(r_star,fit)

  add_g <- as.data.frame(w * fit, ncol = 1)
  stop_value <- 1e-3

  number_of_added_value <- sum((abs(add_g) < stop_value)*1)

  if (iterations == stop_times){
    break
  }

  if (number_of_added_value != n && iterations!= stop_times){
    df[as.character(variable)] <- df[,variable] + add_g

    sum_of_every_fitted_value <- sum_of_every_fitted_value + w * fit
    r_star = y_star - sum_of_every_fitted_value
    r <- y - sum_of_every_fitted_value
    df1 <- data.frame(corrected_data$censoring_indicator,r_star,r)
    colnames(df1)[1] <- "censoring_indicator"
    survival_probability <- c()

    for (i in c(1:dim(df1)[1])){
      yy <- df1[df1$r<=df1$r[i],]
      u <- c(yy$r)
      data_producted <- 1
      for (j in c(1:length(u))){
        total_sum_of_denominator <- sum((df1$r>=u[j])*1)
        censoring_indicator_1<- df1[df1$censoring_indicator==1,]
        total_sum_of_numerator <- sum((censoring_indicator_1$r==u[j])*1)
        data_producted = data_producted * (1-total_sum_of_numerator/total_sum_of_denominator)
      }
      survival_probability[i] <- data_producted
    }

    survival_probability[which(survival_probability==0)] <- 0.000001

    sum_riemann_for_every_n <- c()
    upper = max(r)

    for (i in c(1:length(df1$r))){
      point <- interval(df1$r[i],upper)
      riemann_x <- c()
      data_producted <- 1
      for (j in c(1:length(point))){
        xx <- df1[df1$r<=point[j],]
        u <- c(xx$r)
        for (k in c(1:length(u))){
          total_sum_of_denominator <- sum((df1$r>=u[k])*1)
          total_sum_of_denominator
          censoring_indicator_1<- df1[df1$censoring_indicator==1,]
          total_sum_of_numerator <- sum((censoring_indicator_1$r==u[k])*1)
          data_producted = data_producted *  (1-total_sum_of_numerator/total_sum_of_denominator)
        }
        riemann_x[j]<- data_producted
      }

      riemann_minus <- c()
      for (l in c(1:length(riemann_x)-1)){
        riemann_minus[l] <- riemann_x[l+1]-riemann_x[l]
      }


      sum_integrate <-c(0)
      for (m in c(1:length(riemann_minus))){
        sum_integrate <- sum_integrate + point[m+1] * riemann_minus[m]
      }
      sum_riemann_for_every_n[i] <- sum_integrate
    }


    y_star_update <- c()

    for (i in c(1:length(y))){
      y_star_update[i] <- censoring_indicator[i] * y[i]+
        (1-censoring_indicator[i])*(sum_of_every_fitted_value[i] - sum_riemann_for_every_n[i]/survival_probability[i])
    }

    y_star <- y_star_update

    iterations = iterations + 1

  }else{
    break
  }
}



variable_catch
a <- table(variable_catch)
a[a>=2]  # 4, 5, 30, 33, 35, 64, 66, 68
aa <- as.numeric(names(a[a>=2]))


for (i in 1:length(aa)){
  name <- paste0(names(a[a>=2])[i],'.jpeg')
  jpeg(name)
  x_1 <- as.data.frame(cbind(data[,aa[i]+2],df[,aa[i]]+0.0))
  colnames(x_1) <- c('x','y')
  x_1 <- x_1[order(x_1$x),]
  plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(data)[aa[i]+2],type='l', lwd=2,lty=1)
  dev.off()
  
  
}




jpeg("n1.jpeg")
# 4
x_1 <- as.data.frame(cbind(corrected_data[,6],df[,4]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[6],type='l', lwd=2,lty=1)
dev.off()


jpeg("n2.jpeg")
# 5
x_1 <- as.data.frame(cbind(corrected_data[,7],df[,5]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[7],type='l', lwd=2,lty=1)
dev.off()


jpeg("n3.jpeg")
# 30
x_1 <- as.data.frame(cbind(corrected_data[,32],df[,30]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[32],type='l', lwd=2,lty=1)
dev.off()



jpeg("n4.jpeg")
# 33
x_1 <- as.data.frame(cbind(corrected_data[,35],df[,33]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[35],type='l', lwd=2,lty=1)
dev.off()


jpeg("n5.jpeg")
# 35
x_1 <- as.data.frame(cbind(corrected_data[,37],df[,35]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[37],type='l', lwd=2,lty=1)
dev.off()


jpeg("n6.jpeg")
# 64
x_1 <- as.data.frame(cbind(corrected_data[,66],df[,64]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[66],type='l', lwd=2,lty=1)
dev.off()


jpeg("n7.jpeg")
# 66
x_1 <- as.data.frame(cbind(corrected_data[,68],df[,66]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[68],type='l', lwd=2,lty=1)
dev.off()


jpeg("n8.jpeg")
# 68
x_1 <- as.data.frame(cbind(corrected_data[,70],df[,68]+0.0))
colnames(x_1) <- c('x','y')
x_1 <- x_1[order(x_1$x),]
plot(x_1$x,x_1$y, xlab="x",ylab='f', main= colnames(corrected_data)[70],type='l', lwd=2,lty=1)
dev.off()


survival_data_c <- df
sur_time <- apply(survival_data_c,1,sum)
sur_time <-exp(sur_time)
sur_time <-sort(sur_time)
sur_time  <- data.frame(sur_time)
pro <- seq(1,144)/144
pro <- sort(pro, decreasing= TRUE)
pro <- data.frame(pro)


# survival probability
ddf3 <- cbind(sur_time, pro)
plot(ddf3$sur_time,ddf3$pro, xlab="time",ylab='survival probability', main='survival curve',type='l', lwd=2,lty=1)
points(ddf1$sur_time,ddf1$pro,type="l", col="red", lwd=2)
points(ddf2$sur_time,ddf1$pro,type="l", col="blue", lwd=2)

jpeg("survival curve.jpeg")
# 68
ddf3 <- cbind(sur_time, pro)
plot(ddf3$sur_time,ddf3$pro, xlab="time",ylab='survival probability', main='survival curve',type='l', lwd=2,lty=1)
lines(ddf1$sur_time,ddf1$pro,lty=2, col="red", lwd=2)
lines(ddf2$sur_time,ddf1$pro,lty=3, col="blue", lwd=2)
dev.off()


