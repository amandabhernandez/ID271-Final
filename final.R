### AUTHORS: AHz, RH, AZ
### Written in: R version 4.2.2
### Purpose: ID 271 Midterm


################################################################################
# LOAD DEPENDENCIES  ####
################################################################################

pacman::p_load(tidyverse, lubridate, mgcv, splines, broom, ggeffects, gtsummary, 
               flextable, cowplot)

#set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat_columns <- c("WARMCOLD", "WKDAY", "RACE", "TWODRINK",
                 "SMK", "CHD", "FBG", "STATIN")

nas_bp <- readRDS("nas_bp.rds") %>% 
  mutate(across(any_of(cat_columns), factor))

nas_bp_long <- nas_bp %>% 
  pivot_longer(names_to = "variable", values_to = "value", -ID)



################################################################################
# What's in the data  ####
################################################################################

skimr::skim(nas_bp)

nas_bp_long %>% 
  subset(variable %in% c("SBP", "TEMPC24H", "BC24H")) %>% 
           ggplot(aes(x = value)) + 
           geom_histogram(aes(y = after_stat(density)),
                          color = "black",
                          fill = "white",
                          bins = 30) +
           geom_density(fill = "red", alpha = 0.25)+ 
           facet_wrap(~variable, scales = "free")
         
ggplot(nas_bp, aes(x = log(SBP))) +
  geom_histogram(aes(y = after_stat(density)),
                 color = "black",
                 fill = "white",
                 bins = 30) +
  geom_density(fill = "red", alpha = 0.25)




################################################################################
# SOME MODELING  ####
################################################################################

# bc we have multiple visits per person, we are going to want to use a mixed effect
# model. SBP appears slightly non-normally distributed, so I've log transformed it, 
# and we can use a gaussian distribution 



# model 1 -- exp + confounders only ---------------------------------------

mod1 <- gamm(log(SBP) ~ s(TEMPC24H, bs = 'cr', fx = FALSE) + 
               s(BC24H, bs = 'cr', fx = FALSE) + 
               s(RHUM24H, bs = 'cr', fx = FALSE) + 
               # control for seasonality 
               SINETIME + COSTIME + WKDAY, 
             random = list(ID =  ~ 1),
             family = gaussian(), data = nas_bp)

summary(mod1$gam)
# all edfs = 1, so they are linear and we can remove the splines
as_flextable(mod1$gam)

# model 2 -- exp + confounders + predictors --------------------------------

# let's examine non-linearity with continuous variable predictors and SBP 
#remember to make a spline for TEMPC24H and BC24H 

mod2 <- gamm(log(SBP) ~ TEMPC24H + BC24H + RHUM24H + 
               s(AGE, bs = 'cr', fx = FALSE) + 
               RACE + TWODRINK + CHD + DIABETE +
               s(BMI, bs = 'cr', fx = FALSE) + SMK + FBG + 
               s(NEDUC, bs = 'cr', fx = FALSE) + STATIN + 
               # control for seasonality 
               SINETIME + COSTIME + WKDAY, 
             random = list(ID =  ~ 1),
             family = gaussian(), data = nas_bp)

#should we be including all the predictors? should we justify our decision based 
#on p-values? 


#get edfs
summary(mod2$gam)
as_flextable(mod2$gam)

#show plots for splines 
plot(mod1$gam,scale=0, pages = 1)

# model 3 -- exp + confounders + predictors (wo splines) -------------------

# we don't need a spline for BMI or NEDUC-- its linear

mod3 <- gamm(log(SBP) ~ TEMPC24H + BC24H + RHUM24H + 
               s(AGE, bs = 'cr', fx = FALSE) + 
               RACE + TWODRINK + CHD + DIABETE +
               BMI + SMK + FBG + 
               NEDUC + STATIN + 
               # control for seasonality 
               SINETIME + COSTIME + WKDAY,
             random = list(ID =  ~ 1),
             family = gaussian(), data = nas_bp)

summary(mod3$gam)
as_flextable(mod3$gam)

# dose response curves
ggpredict(mod3$gam, terms = "TEMPC24H") %>% 
  plot()

ggpredict(mod3$gam, terms = "BC24H") %>% 
  plot()

# hist(mod3$lme$coef$random$ID[,1] + mod3$lme$coef$fixed[2], 
#      main="Histogram of random slopes by ID",xlab="Random slopes" )



# alternate way... 
# pred_wtemp <- exp(predict(mod2$gam,type='terms')[,1])
# 
# plot(nas_bp$TEMPC24H, pred_wtemp)
# 
# pred_wbc <- exp(predict(mod2$gam,type='terms',se.fit=TRUE)$fit[,2])
# 
# plot(nas_bp$BC24H, pred_wbc)
# 
# plot(mod2$gam)

# 
# mod1 <- gamm(log(SBP) ~ TEMPC24H + BC24H + RHUM24H + SINETIME + COSTIME + as.factor(WKDAY), 
#              random = list(ID=~1),
#              #treat ID as random intercept to control for
#              #highly-correlated nature of repeat measurements within subjects 
#              family=gaussian(),
#              #we have normally distributed continuous data 
#              data=nas_bp_test)
# 
# mod1_pred_TEMPC24H <- ggpredict(mod1$gam, terms =  "TEMPC24H", type = "fe")
# mod1_pred_BC24H <- ggpredict(mod1$gam, terms= "BC24H")




# get effect estimate and 95% CIs
data.frame(intervals(mod3$lme)$fixed)[2:3,] %>% 
  rownames_to_column(var = "term") 



################################################################################
# CHECK FOR EFFECT MODIFICATION  ####
################################################################################


# check for effect modification by diabetes status

mod3_diab <- gamm(log(SBP) ~ TEMPC24H + BC24H + RHUM24H + 
                    s(AGE, bs = 'cr', fx = FALSE) + 
                    RACE + TWODRINK + CHD + DIABETE +
                    BMI + SMK + FBG + 
                    NEDUC + STATIN + 
                    DIABETE*BC24H + 
                    # control for seasonality 
                    SINETIME + COSTIME + WKDAY,
                  random = list(ID =  ~ 1),
                  family = gaussian(), data = nas_bp)

tidy(mod3_diab$gam, parametric = TRUE)

# 
# # get estimates 
# # for people wo diabetes (DIABETE = 0) get the effect of same day BC on SBP 
# data.frame((mod3_diab$lme$coef$fixed)) %>% 
#   rownames_to_column(var = "term")  %>% 
#   subset(term == "XBC24H")
# 
# # for people w diabetes (DIABETE = 1) get the effect of same day BC on SBP
# 
# data.frame((mod3_diab$lme$coef$fixed)) %>% 
#   rownames_to_column(var = "term")  %>% 
#   filter(term %in% c("XBC24H", "XBC24H:DIABETE")) %>% 
#   group_by(.) %>% 
#   summarize(sum(X.mod3_diab.lme.coef.fixed.))
# 

# get estimates and std errors 
mod3_diab_interaction_stderrs <- data.frame(mod3_diab$lme$varFix) %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "term2", values_to = "stderr", -term) %>% 
  filter(str_detect(term, "BC24H") & str_detect(term2, "BC24H"))

# get the estimates and std errors 
broom.mixed::tidy(mod3_diab$gam, parametric = TRUE, conf.int = TRUE) %>% 
  filter(term %in% c("BC24H", "BC24H:DIABETE")) %>% 
  select(term, estimate) %>% 
  group_by(.) %>% 
  mutate(estimate = case_when(term == "BC24H:DIABETE" ~ sum(estimate),
                              TRUE ~ estimate)) %>% 
  left_join(data.frame(term = c("BC24H", "BC24H:DIABETE"),
                       std.err = c(sqrt(mod3_diab_interaction_stderrs$stderr[which(mod3_diab_interaction_stderrs$term == "XBC24H" & mod3_diab_interaction_stderrs$term2 == "XBC24H")]),
                                   sqrt(sum(mod3_diab_interaction_stderrs$stderr))))) 



# do this for all categorical subgroups? 

model_fun <- list()
mod3_emm <- list()
int_term <- list()
for(i in cat_columns){
  
  model_fun[[i]] <- paste0("log(SBP) ~ TEMPC24H + BC24H + RHUM24H + s(AGE, bs = 'cr', fx = FALSE) + RACE + TWODRINK + CHD + DIABETE + BMI + SMK + FBG + NEDUC + STATIN + ",
                      i, "*BC24H + SINETIME + COSTIME + WKDAY")
  
  mod3_emm[[i]] <- gamm(as.formula(model_fun[[i]]),
                    random = list(ID =  ~ 1),
                    family = gaussian(), data = nas_bp)
  
  int_term[[i]] <- tidy(mod3_emm[[i]]$gam, parametric = TRUE) %>% 
    filter(str_detect(term, ":"))
  
}

inspect_int_terms <- bind_rows(int_term, .id = "model") %>% 
  mutate(p = case_when(p.value < 0.001 ~ "***",
                       p.value < 0.01 ~ "**",
                       p.value < 0.05 ~ "*",
                              TRUE ~ paste0(round(p.value, 4))))



# stratified model by WARMCOLD
# 
# mod2_warmcold0 <- gamm(log(SBP) ~ TEMPC24H + BC24H + RHUM24H + 
#                     # control for seasonality 
#                     SINETIME + COSTIME + WKDAY,
#                   random = list(ID =  ~ 1),
#                   family = gaussian(), data = nas_bp[nas_bp$WARMCOLD == "0",])
# 
# mod2_warmcold1 <- gamm(log(SBP) ~ TEMPC24H + BC24H + RHUM24H + 
#                     # control for seasonality 
#                     SINETIME + COSTIME + WKDAY,
#                   random = list(ID =  ~ 1),
#                   family = gaussian(), data = nas_bp[nas_bp$WARMCOLD == "1",])
# 
# 
# as_flextable(mod2_warmcold0$gam)
# as_flextable(mod2_warmcold1$gam)
# 
# warmcold_strat_0 <- ggpredict(mod2_warmcold0$gam, terms = "TEMPC24H")
# 
# warmcold_strat_1 <- ggpredict(mod2_warmcold1$gam, terms = "TEMPC24H")
# 
# ggplot(warmcold_strat_0) +
#   geom_point(nas_bp, mapping = aes(x = TEMPC24H, y = SBP, color = WARMCOLD), alpha = 0.1, show.legend = FALSE) +
#   geom_line(aes(x, y = predicted, color = "Cold season")) +
#   geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = .1) + 
#   geom_line(warmcold_strat_1, mapping = aes(x, y = predicted, color = "Warm season")) + 
#   geom_ribbon(warmcold_strat_1, mapping = aes(x = x, ymin = conf.low, ymax = conf.high), fill = "red", alpha = .1) + 
#   scale_color_manual(name = "Season", values = c("Cold season" = "blue", "Warm season" = "red", "0" = "blue", "1" = "red"))
# check: 
# mod2_emm$lme$coef$fixed[3] + mod2_emm$lme$coef$fixed[15]
# sqrt(mod2_emm$lme$varFix[9,9])
# sqrt(sum(mod2_emm$lme$varFix[c(9,15), c(9,15)]))

# yep, we got the same thing 
# 
# 
# check for effect modification by WARMCOLD status
# check for significance for p value 

## make some histograms for the effect heterogeneity 

## lagged effects 






# subgroups -- should we do interaction by race 
# 
# 
# 
# 
# GARBAGE
# # mod2_emm <- gamm(log(SBP) ~ TEMPC24H  + 
#                    DIABETE + factor(SMK) + 
#                    s(AGE, bs = 'cr', fx = FALSE) + factor(RACE) + 
#                    s(BMI, bs = 'cr', fx = FALSE) + 
#                    DIABETE*BC24H + 
#                    # control for seasonality 
#                    SINETIME + COSTIME + as.factor(WKDAY), 
#                  random = list(ID =  ~ 1),
#                  family = gaussian(), data = nas_bp)



temp_dl <- as.matrix(mod3_lag_tidy$estimate[str_detect(mod3_lag_tidy$term, "MPC")])
rownames(temp_dl) <- mod3_lag_tidy$term[str_detect(mod3_lag_tidy$term, "MPC")]

temp_beta_matrix <- ll%*%as.matrix(temp_dl)
temp_cov <-  ll%*%as.matrix(vcov(mod3_lag$gam)[2:4,2:4])%*%t(ll)
temp_se <- sqrt(diag(temp_cov))

temp_est_se <- data.frame(var=c('temp_lag3','temp_lag1','temp_lag0'),
                          coef=temp_beta_matrix,
                          se=temp_se)



bc_dl <- as.matrix(mod3_lag_tidy$estimate[str_detect(mod3_lag_tidy$term, "BC")])
rownames(bc_dl) <- mod3_lag_tidy$term[str_detect(mod3_lag_tidy$term, "BC")]


temp_est_se <- data.frame(var=c('temp_lag2','temp_lag1','temp_lag0'),
                          coef=temp_beta_matrix,
                          se=temp_se)


mod3_lag_tidy <- tidy(mod3_lag$gam, parametric = TRUE)

temp_dl <- mod3_lag_tidy[str_detect(mod3_lag_tidy$term, "MPC"),]

lags <- c(0:2)


gplots::plotCI(y = temp_dl$estimate, x = lags, uiq = 1.96*(temp_dl$std.error), err = "y")


ll <- matrix(nrow=4,ncol=3)
ll[1,] <- rep(1,3) #intercept
ll[2,] <- 0:2 #lag0 to lag2, i term
ll[3,] <- ll[2,]^2 #i^2 term
ll[4,] <- ll[2,]^3 #i^3 term


ll <- t(ll) #transpose





lagm <- function(x, n) {
  nn <- length(x)
  xn <- c(rep(NA, n), x[1:(nn - n)])
  return(xn)
}

ll <- matrix(nrow=4,ncol=6)
ll[1,] <- rep(1,6) #intercept
ll[2,] <- 0:5 #lag0 to lag5, i term
ll[3,] <- ll[2,]^2 #i^2 term
ll[4,] <- ll[2,]^3 #i^3 term

ll <- t(ll) #transpose

lent <- length(nas_bp$ID)
xx <- matrix(nrow = lent, ncol = 6)
yy <- matrix(nrow = lent, ncol = 6)
for (i in 1:6) {
  xx[, i] <- lagm(nas_bp$TEMPC24H, (i - 1))
  yy[,i] <-  lagm(nas_bp$BC24H, (i - 1))
}

## Then the Z's are

zz<-xx%*%ll
zz2 <- xx%*%ll

#these are the same -- the lag part worked right 
summary(xx[,2])
summary(nas_bp$TMPCD1)


mod2_lag <- gamm(log(SBP) ~ zz[,1] + zz[,2] + zz[,3] + 
                   zz2[,1] + zz2[,2] + zz2[,3] +
                   RHUM24H +
                   # control for seasonality 
                   SINETIME + COSTIME + WKDAY,
                 random = list(ID =  ~ 1),
                 family = gaussian(), data = nas_bp)

mod2_lag_tidy <- tidy(mod2_lag$gam, parametric = TRUE)


temp_dl <- as.matrix(mod2_lag_tidy$estimate[str_detect(mod2_lag_tidy$term, "MPC")])
rownames(temp_dl) <- mod2_lag_tidy$term[str_detect(mod2_lag_tidy$term, "MPC")]

temp_beta_matrix <- ll%*%as.matrix(temp_dl)
temp_cov <-  ll%*%as.matrix(vcov(mod2_lag_tidy$gam)[2:4,2:4])%*%t(ll)
temp_se <- sqrt(diag(temp_cov))

temp_est_se <- data.frame(var=c('temp_lag3','temp_lag1','temp_lag0'),
                          coef=temp_beta_matrix,
                          se=temp_se)


bc_dl <- as.matrix(mod2_lag_tidy$estimate[str_detect(mod2_lag_tidy$term, "BC")])
rownames(bc_dl) <- mod2_lag_tidy$term[str_detect(mod2_lag_tidy$term, "BC")]




# create empty matrix
x_weights <- matrix(nrow=4,ncol=3)

# fill in matrix
x_weights[1,] <- rep(1,3) # for the intercept
x_weights[2,] <- 0:2 # for lag0 to lag5 of ozone, i term
x_weights[3,] <- x_weights[2,]^2 # i^2 term
x_weights[4,] <- x_weights[2,]^3 # i^3 term

# We transpose the matrix so that we can conduct matrix multiplication
x_weights <- t(x_weights)


# create matrix with the lagged ozone
lag_mat <- matrix(nrow=length(nas_bp$ID), ncol=3)
lag_mat[,1] <- nas_bp$TEMPC24H
lag_mat[,2] <- nas_bp$TMPCD1
lag_mat[,3] <- nas_bp$TMPCD2

for (i in 1:3) {
  lag_mat[,i]<-lag_ozone(nas_bp$TEMPC24H,(i-1))
}

# Let's calculate the 'z' coefficients
z_coeff<-lag_mat%*%x_weights

mod2_lag <- gamm(log(SBP) ~ z_coeff[,1] + z_coeff[,2] + z_coeff[,3] + 
                   BC24H + BCD1 + BCD2 +
                   RHUM24H +
                   # control for seasonality 
                   SINETIME + COSTIME + WKDAY,
                 random = list(ID =  ~ 1),
                 family = gaussian(), data = nas_bp)



# Extract 'c' coefficients
c_coeff <- as.matrix(mod2_lag$gam$coefficients[2:4])

# Now multiple them by the x weights to get the 'b's, which are what we are interested in.
b_coeff <- x_weights%*%as.matrix(c_coeff)

# Typically, we are also interested in inference (e.g., 95% CI). So let's first extract the covariance matrix from the model. Then, we need to multiply them by the weights squared
ozone_covar <- x_weights%*%as.matrix(vcov(mod)[8:11, 8:11])%*%t(x_weights)

# We square root the variance, which are in the diagonal of the covariance matrix, to get the standard error
ozone_se <- sqrt(diag(ozone_covar))

# Last step! Let's put the beta estimates and the standard errors together in a single dataframe for easy plotting latter.
ozone_results <- data.frame(lag=c(5:0),
                            coef=b_coeff, se=ozone_se)


mod_tidy <- broom::tidy(mod2_lag$gam, parametric = TRUE)
o3_dl <- as.matrix(mod_tidy$estimate[str_detect(mod_tidy$term, "z_")])
rownames(o3_dl) <- mod_tidy$term[str_detect(mod_tidy$term, "z_")]


#now multiple them by the weights to become b's
o3_est <- lag_mat%*%as.matrix(o3_dl) #beta matrix

# o3_cov <- ll%*%as.matrix(vcov(mod)[8:11,8:11])%*%t(ll)
# replace the code above with this
o3_cov <- lag_mat%*%as.matrix(vcov(mod2_lag$lme)[str_detect(rownames(vcov(mod2_lag$lme)), "z_"), str_detect(colnames(vcov(mod2_lag$lme)), "z_")])%*%t(ll)

o3_se <- sqrt(diag(o3_cov))



# mod2_bc_warmcold <- gamm(SBP ~ TEMPC24H  + 
#                        BC24H*WARMCOLD + 
#                        RHUM24H + 
#                        BMI + 
#                        s(AGE, bs = 'cr', fx = FALSE) + 
#                        NEDUC + 
#                        SMK + TWODRINK + STATIN + FBG + 
#                        # control for seasonality 
#                        SINETIME + COSTIME + WKDAY, 
#                      random = list(ID =  ~ 1),
#                      family = gaussian(), data = nas_bp)


# as_flextable(mod2_bc_warmcold$gam)


# warmcold0_bc <- ggpredict(mod2_bc_warmcold$gam, terms = "BC24H", condition = c(WARMCOLD = "0"))
# warmcold1_bc <- ggpredict(mod2_bc_warmcold$gam, terms = "BC24H", condition = c(WARMCOLD = "1"))

# ggplot(warmcold0_bc) +
#   #geom_point(nas_bp, mapping = aes(x = TEMPC24H, y = SBP), alpha = 0.1) +
#   geom_line(aes(x, y = predicted, color = "Cold season")) +
#   geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = "blue", alpha = .1) + 
#   geom_line(warmcold1_bc, mapping = aes(x, y = predicted, color = "Warm season")) + 
#   geom_ribbon(warmcold1_bc, mapping = aes(x = x, ymin = conf.low, ymax = conf.high), fill = "red", alpha = .1) + 
#   scale_color_manual(name = "Season", values = c("Cold season" = "blue", "Warm season" = "red")) + 
#   labs(x = "TEMPC24H", y = "Predicted SBP")

