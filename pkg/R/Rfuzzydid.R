#' @title wald_did
#' @export
#' @description  Computes the Wald-DiD estimator, equivalent to an IV regression of y on d with
#' g and t as included instruments and g:t as the excluded instrument for d
#' @param i Rows of the data frame to operate on
#' @param df Data frame to operate on
#' @param y_name Outcome variable
#' @param g_name Group variable
#' @param t_name Time period variable
#' @param d_name Treatment variable, can be any ordered variable
#' @param X_name Covariates to include
wald_did = function(i, df, y_name, g_name, t_name, d_name, X_name = NULL){
	df = df[i,]
#	select_df = function(g, t, var) {
#		return(df[df[[g_name]] == g & df[[t_name]] == t, var])
#	}
#	if (is.null(X_name)) {
#		numerator = mean(select_df(1, 1, y_name), na.rm=TRUE) - mean(select_df(1, 0, y_name), na.rm=TRUE) - mean(select_df(0, 1, y_name), na.rm=TRUE) + mean(select_df(0, 0, y_name), na.rm=TRUE)
#		denominator = mean(select_df(1, 1, d_name), na.rm=TRUE) - mean(select_df(1, 0, d_name), na.rm=TRUE) - mean(select_df(0, 1, d_name), na.rm=TRUE) + mean(select_df(0, 0, d_name), na.rm=TRUE)
#		return(numerator/denominator)
#	}
	if(is.null(X_name)) {
		X_name = "wdidconstwdid"
		df[[X_name]] = 1
	}

	X_formula = paste(X_name, collapse = "+")
	yX_formula = paste0(y_name, " ~ 0 + ", X_formula)
	yX_formula = as.formula(yX_formula)
	dX_formula = paste0(d_name, " ~ 0 + ", X_formula)
	dX_formula = as.formula(dX_formula)

	reg_y_10 = lm(yX_formula, data = df[df[[g_name]] == 1 & df[[t_name]] == 0, ])
	reg_y_01 = lm(yX_formula, data = df[df[[g_name]] == 0 & df[[t_name]] == 1, ])
	reg_y_00 = lm(yX_formula, data = df[df[[g_name]] == 0 & df[[t_name]] == 0, ])

	reg_d_10 = lm(dX_formula, data = df[df[[g_name]] == 1 & df[[t_name]] == 0, ])
	reg_d_01 = lm(dX_formula, data = df[df[[g_name]] == 0 & df[[t_name]] == 1, ])
	reg_d_00 = lm(dX_formula, data = df[df[[g_name]] == 0 & df[[t_name]] == 0, ])

	y11 = df[df[[g_name]] == 1 & df[[t_name]] == 1, y_name]
	X11 = as.data.frame(df[df[[g_name]] == 1 & df[[t_name]] == 1, X_name])
	colnames(X11) = X_name
	d11 = df[df[[g_name]] == 1 & df[[t_name]] == 1, d_name]
	numerator = sum(y11 - predict(reg_y_10, X11) - predict(reg_y_01, X11) + predict(reg_y_00, X11))
	denominator = sum(d11 - predict(reg_d_10, X11) - predict(reg_d_01, X11) + predict(reg_d_00, X11))

	wdid = numerator / denominator
	return(wdid)
}

#' @title wald_tc
#' @description  Computes the time-corrected Wald-DiD estimator as in CH2018 <doi:10.1093/restud/rdx049>
#' @param i Rows of the data frame to operate on
#' @param df Data frame to operate on
#' @param y_name Outcome variable
#' @param g_name Group variable
#' @param t_name Time period variable
#' @param d_name Treatment variable, can be any ordered variable
#' @param X_name Covariates to include
wald_tc = function(i, df, y_name, g_name, t_name, d_name, X_name = NULL){
	df = df[i,]
	if (is.null(X_name)) {
		X_name = "wtcconstwtc"
		df[[X_name]] = 1
	}

	X_formula = paste(X_name, collapse = "+")
	yX_formula = paste0(y_name, " ~ 0 + ", X_formula)
	yX_formula = as.formula(yX_formula)
	dX_formula = paste0(d_name, " ~ 0 + ", X_formula)
	dX_formula = as.formula(dX_formula)
	
	reg_y_10 = lm(yX_formula, data = df[df[[g_name]] == 1 & df[[t_name]] == 0,])
	reg_y_d01 = list()
	reg_y_d00 = list()
	delta_d = list()
	for (d in unique(df[[d_name]])) {
		n_d01 = nrow(df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 1,])
		n_d00 = nrow(df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 0,])
		D = as.character(d)
		reg_y_d01[[D]] = lm(yX_formula, data = df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 1,])
		reg_y_d00[[D]] = lm(yX_formula, data = df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 0,])
		tc_1 = sum(df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 1, y_name]) / n_d01
		tc_2 = sum(df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 0, y_name]) / n_d00
		delta_d[[D]] = tc_1 - tc_2
	}


	y11 = df[df[[g_name]] == 1 & df[[t_name]] == 1, y_name]
	d11 = df[df[[g_name]] == 1 & df[[t_name]] == 1, d_name]
	n11 = length(y11)
	y10 = df[df[[g_name]] == 1 & df[[t_name]] == 0, y_name]
	d10 = df[df[[g_name]] == 1 & df[[t_name]] == 0, d_name]
	n10 = length(y10)
	sum_corrections = 0
	for (i in 1:n10) {
		elem = y10[i] + delta_d[[as.character(d10[i])]]
		sum_corrections = sum_corrections + elem
	}
	numerator = mean(y11) - (sum_corrections / n10)
	denominator = mean(d11) - mean(d10)

	wtc = numerator / denominator
	return(wtc)
}

#' @title fuzzydid
#' @description Computes the corrected LATE estimators as in de Chaisemartin and D'Haultfoeuille (2018)
#' <doi:10.1093/restud/rdx049>
#' @param df Data frame to operate on
#' @param y_name Outcome variable
#' @param g_name Group variable
#' @param t_name Time period variable
#' @param d_name Treatment variable, can be any ordered variable
#' @param X_name Covariates to include
#' @param est Estimators to compute. Should contain at least one of "did", "tc"
#' @param nboot Number of bootstrap samples for standard errors
fuzzydid = function(df, y_name, g_name, t_name, d_name, X_name = NULL, est = c("did", "tc"), nboot = 50) {
	b = list()
	boot_se = list()
	ci95 = list()

	if ("did" %in% est) {
		print("Computing WDID...")
		b[["did"]] = wald_did(1:nrow(df), df, y_name, g_name, t_name, d_name)
		boot = bootstrap(1:nrow(df), 50, wald_did, df, y_name, g_name, t_name, d_name)
		boot_se[["did"]] = sd(boot$thetastar)
		ci95[["did"]] = quantile(boot$thetastar, c(0.05, 0.95))
	}

	if ("tc" %in% est) {
		print("Computing WTC...")
		b[["tc"]] = wald_tc(1:nrow(df), df, y_name, g_name, t_name, d_name)
		boot = bootstrap(1:nrow(df), 50, wald_tc, df, y_name, g_name, t_name, d_name)
		boot_se[["tc"]] = sd(boot$thetastar)
		ci95[["tc"]] = quantile(boot$thetastar, c(0.05, 0.95))
	}
	result = list("b" = b, "boot_se" = boot_se, "ci95" = ci95)
	class(result) = "fuzzydid"
	return(result)
}

#' @title summary.fuzzydid
#' @export
#' @description Displays the estimations in a table
#' @param object A fuzzydid object
summary.fuzzydid <- function(object) {
  b <- object$b
  boot_se <- object$boot_se
  ci95 <- object$ci95

  summary_df <- data.frame(
    Estimator = paste0("W_", toupper(names(b))),
    Point_Estimate = unlist(b),
    Bootstrap_SE = unlist(boot_se),
    CI_Lower = sapply(ci95, `[[`, 1),
    CI_Upper = sapply(ci95, `[[`, 2)
  )
	rownames(summary_df) = NULL

  knitr::kable(summary_df, caption = "Summary of Estimators")
}

