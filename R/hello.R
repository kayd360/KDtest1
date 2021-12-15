
#' @title overview_tab
#'
#' @description Provides an overview table for the time and scope conditions of
#'     a data set
#'
#' @param dat A data set object
#' @param id Scope (e.g., country codes or individual IDs)
#' @param time Time (e.g., time periods are given by years, months, ...)
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#'
#'  # ---------------------------------------------------------
#'  # Libraries and Data preparation
#'  # ---------------------------------------------------------
#'
#'  ## --------------------------------
#'  ## Basic Libraries
#'  ## --------------------------------
#'  rm(list = ls())
#'  library(lattice)
#'  library(rgl)
#'  library(car)
#'  library(rsm)
#'  library(effects)
#'  library(splines)
#'  library(relaimpo)
#'  library(lmtest)
#'
#'  ## --------------------------------
#'  ## Data Load
#'  ## --------------------------------
#'  # loading from txt
#'  Sal = read.table("P130.txt", header = TRUE, sep = "\t", dec = ".", col.names = c("salary", "exper", "educ", "manag"))
#'
#'  # loading from .RData file
#'  load("USTemp.RData")
#'
#'  # is dataset complete
#'  all(complete.cases(USTemp))
#'
#'  ## --------------------------------
#'  ## Factoring the Data and Revelling
#'  ## --------------------------------
#'  Sal$educ = factor(Sal$educ, levels = c(1, 2, 3), labels = c("HighS", "Bachelor", "Master"))
#'  Sal$manag = factor(Sal$manag, levels = c(0, 1), labels = c("NonMan", "Man"))
#'  Sal$educ = relevel(Sal$educ, ref = "Master")
#'
#'
#'  # ---------------------------------------------------------
#'  # Modelling
#'  # ---------------------------------------------------------
#'
#'  ## --------------------------------
#'  ## Creating Models
#'  ## --------------------------------
#'
#'  fit = lm(salary ~ educ + manag + exper, data = Sal) # without interaction
#'  fitInter = lm(salary ~ educ + manag + exper + educ:manag, data = Sal) # with interaction terms
#'
#'  fit = lm(temp ~ poly(long, degree = 3, raw = TRUE) + poly(lat, degree = 3, raw = TRUE), data = USTemp) # polynomial model
#'  fitI = lm(temp ~ long + I(long^2) + I(long^3) + lat + I(lat^2) + I(lat^3), data = USTemp) # polynomial using Identitiy equation
#'  fitR = lm(temp ~ poly(long, degree = 3, raw = TRUE) + lat, data = USTemp) # poly longitude + linear lattitude
#'
#'  # with knots / splines
#'  fit = lm(temp ~ bs(long, degree = 1, knots = c(-80, -100)) + bs(lat, degree = 1, knots = c(30, 40)), data = USTemp)
#'  # function bs() will calculate the B-splines
#'  # argument "degree" sets the degree of a polynomial function; degree = 1 for a straight line regression function
#'  # argument "knots" sets values of knots into particular values of the explanatory variables (long and lat)
#'
#'  # alternative splines with df
#'  fit2 = lm(temp ~ bs(long, degree = 1, df = 1 + 2) + bs(lat, degree = 1, df = 1 + 2), data = USTemp)
#'
#'  # slicing the data for splines
#'  new.data = expand.grid(lat = 35, long = seq(-71, -120, by = -1))
#'  plot(seq(-71, -120, by = -1), predict(fit, new.data), type = "l") # plotting slice data
#'
#'  coef(fit)
#'  head(model.matrix(fit))
#'
#'  ## --------------------------------
#'  ## Plotting
#'  ## --------------------------------
#'
#'  # basic plots
#'  plot(resid(fit) ~ Sal$exper)
#'  plot(resid(fit) ~ interaction(Sal$educ, Sal$manag))
#'
#'  # x-y plots
#'  xyplot(salary ~ exper | manag, group = educ, auto.key = TRUE, data = Sal)
#'  scatter3d(temp ~ long + lat, data = USTemp, surface = TRUE) # scatter3d plot has Y axis
#'
#'  # contour plots
#'  contour(fit, lat ~ long, image = TRUE) # contour plot has Y-axis and then X-axis
#'  text(USTemp$long, USTemp$lat, labels=USTemp$city, cex= 0.6) # but everywhere else it is X-axis to Y-axis
#'
#'  # perspective or surface plots
#'  persp(fit, lat ~ long, col = rev(heat.colors(5)), contours = list(z = "top", col="blue"), theta = 60, phi = 45)
#'
#'  # effect plots
#'  eff1 = Effect(c("lat", "long"), fit, xlevels = list(long = seq(-120, -70, by = 1), lat = seq(25, 45, by = 5)))
#'  plot(eff1, x.var = "long", lines = list(multiline = TRUE)) #Plotting longitude holding latitude constant
#'
#'  eff2 = Effect(c("lat", "long"), fit, xlevels = list(long = seq(-120, -70, by = 10), lat = seq(25, 45, by = 1)))
#'  plot(eff2, x.var = "lat", lines = list(multiline = TRUE))  #Plotting lattitude holding longitude constant
#'
#'  plot(distance/speed ~ speed, data = brake)
#'  points(brake$speed, fitted(fit_full), pch = 4, col = "blue")
#'
#'
#'
#'  ## --------------------------------
#'  ## TESTS
#'  ## --------------------------------
#'
#'  # get : RSS, p value, degrees of freedom,
#'  anova(fit, fitInter)
#'
#'  # residual sum of squares from model
#'  (Qe = sum(resid(fit)^2))
#'  (QeR = sum(resid(fitR)^2))
#'
#'  # value of the F test statistic
#'  Ftest = anova(fit, fitInter)
#'  Ftest$F[2]
#'
#'  # squared value of the t test statistic
#'  ttest = as.data.frame(coef(summary(fitInter)))
#'  ttest$`t value`[8]^2
#'
#'  # Assess the relative importance of the explanatory variables
#'  (relImp = calc.relimp(fit_CC, type = "lmg"))
#'
#'  # Externally studentized Residuals
#'  plot(rstudent(fit_temp) ~ fitted(fit_temp))
#'  abline(h = 0, lty = 3)
#'  residualPlots(fit_temp, tests = FALSE, quadratic = FALSE, id = list(n=3), type = "rstudent")
#'
#'  # Cooks Distance
#'  head(sort(cooks.distance(fit_temp), decreasing = TRUE))
#'
#'  ## --------------------------------
#'  ## INTERPRETATIONS AND CONCLUSIONS
#'  ## --------------------------------
#'
#'  # EXTERNAL STUDENTIZED RESIDUALS
#'  #---------------------------
#'  # : using these charts, we may subjectivelly judge the 1st two assumptions
#'  # 1: correct specification of the regression model;
#'  # 2: constant variance of the error term)
#'  # from the strong set of assumptions
#'  # the externally studentized residuals seem to be quite randomly distributed around zero in the case of plots with longitude and fitted values on the x-axis,
#'  # while there may be seen slightly "quadratic" tendency in the case of the chart with latitude on the x-axis;
#'  # CONCLUSION: we found no obvious support against the assumptions, but a model with latitude in the form of quadratic polynomial could be also studied
#'
#'  # COOKS DISTANCE
#'  #---------------------------
#'  # we sort the observations by Cook's distance (in the descending order)
#'  # using the "rule of thumb", as no observations have the Cook's distance higher than 1, there are no influential observations
#'
#'  # BONEFERRI'S CORRECTION REGRESSION OUTLIERS
#'  #---------------------------
#'  n = dim(model.matrix(fit_temp))[1] # number of observations (US cities); n = 56
#'
#'  p = dim(model.matrix(fit_temp))[2]
#'  # number of parameters of the regression model; p = 5 (Intercept, 3 parameters of the cubic polynomial for longitude, 1 parameter for the latitude)
#'
#'  max(abs(rstudent(fit_temp))) > qt(1 - 0.05/(2 * n), n - p - 1)
#'  # "at least one" => it is enough to judge the observation with the highest residual (if this one is not a regression outlier that there is no regression outlier at all)
#'  # we perform the test at significance level 0.05; Bonferroni correction ("* n") is applied since we are testing multiple hypotheses (because we are testing the residual which is the largest in the absolute value; imagine we "do not know", which one it is). If the observation (residual) to be tested was chosen in advance, we would be testing just one hypothesis and no correction would be applied.
#'  # we do not reject the null hypothesis that there are no regression outliers (ie it seems there are none regression outliers)
#'
#'  # ADDITIVE RELATIONSHIP
#'  #---------------------------
#'  # a regression model providing a correct description of reality: a model with satisfied weak set of assumptions (about the error term, ie epsilon)
#'  # in practice, you will never know (for sure); you may only find some support either in favour or against the assumptions
#'  # for solving this task we plot residuals (y-axis) against the explanatory variables (x-axis)
#'  # a support that regression model provides a correct description of reality: residuals randomly distributed around zero
#'
#'  # firstly, we plot the residuals against the exper variable (exper is a quantitative variable: function plot() creates a scatterplot)
#'  # there are no obvious orderlinesses
#'  # secondly, we may plot the residuals against the manag and educ variables (these are qualitative variables: function plot() creates boxplots)
#'  # for a precision, we construct 6 boxplots (one for each combination of educ and manag levels)
#'  # the assumption that the regression model provides a correct description of reality is not correct sinc
#'
#'  # HYPOTHESIS TESTING
#'  #---------------------------
#'  # anova(fit,fit_Inter)
#'  # Statistical hypothesis test about beta7=0 (ie the interaction term is not important in the model; ie the relationship is additive)
#'  # as the p-value is higher than 0.01, at 1% significance level the model assuming additivity (ie the model without the interaction term) seems acceptable
#'
#'  # HYPOTHESIS F and T TESTING
#'  #---------------------------
#'  # anova(fit, fitInter) --- T Test
#'  # coef(summary(fitInter)) -- F test
#'  # Task 3) Illustrate the equivalence between the general linear hypothesis F test and the standard t-test making use of appropriate outputs from the R software
#'  # the p-value of the general linear hypothesis F test "anova(fit, fitInter)" equals 0.08469
#'  # the p-value of the standard t-test for the beta7 (the interaction term) equals 8.468718e-02, ie the p-values are indeed the same
#'
#'
#'  # LACK OF FIT
#'  #---------------------------
#'  # anova(fitSL, fit_full)
#'  # to assess the straight line model, we run (at a 5% significance level) lack of fit statistical hypothesis test; we compare the straight line model with the "full" model
#'  # (with many dummy variables)
#'  # note that the lack of fit test is a special case (example) of the general linear hypothesis test
#'  # Regression line (Model 1 in R output) residual sum of squares = 12.3012
#'  # Sum of squares due to Pure Error: Q_PE(y) = 7.6932; this is the residual sum of squares of the "full" model (with many dummy variables;
#'  # Model 2 in R output; this is the lowest possible residual sum of squares, which can be achieved for the modelled relationship
#'  # Sum of squares due to Lack of Fit:  Q_LOF(y) = 4.608; mathematically this is squared norm of vector of differences between fitted values of the two compared models;
#'  # sum of squared distances of the blue crosses from the straight line in above chart
#'  # null (tested) hypothesis: straight line model provides a correct description of the modelled relationship
#'  # Conclusion: as the p-value=0.7728 is higher than our significance level 0.05, we do not reject H0; this result of the test suggests that relationship between speed and distance/speed truly has a nature of straight line
#'
#'  ## Analysis of Variance Table
#'  ##
#'  ## Model 1: distance/speed ~ speed
#'  ## Model 2: distance/speed ~ factor(speed)
#'  ##   Res.Df     RSS Df Sum of Sq      F Pr(>F)
#'  ## 1     61 12.3012
#'  ## 2     34  7.6932 27     4.608 0.7543 0.7728
#' @export
#' @importFrom dplyr "%>%"

hello<-function()
{

}

