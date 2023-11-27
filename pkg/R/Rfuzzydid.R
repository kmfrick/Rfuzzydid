#' @title wald_did
#' @description  Computes the Wald-DiD estimator, equivalent to an IV regression of y on d with
#' g and t as included instruments and g:t as the excluded instrument for d
#' @param i Rows of the data frame to operate on
#' @param df Data frame to operate on
#' @param y_name Outcome variable
#' @param g_name Group variable
#' @param t_name Time period variable
#' @param d_name Treatment variable, can be any ordered variable
#' @param X_name Covariates to include
#' @importFrom stats "as.formula" "lm" "predict"
wald_did = function(i, df, y_name, g_name, t_name, d_name, X_name = NULL){
	df = df[i,]
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

	#print(paste0("E[y11] = ", mean(y11), "; E[y10+delta_d] = ", sum_corrections / n10))
	#print(paste0("E[d10] = ", mean(d11), "; E[d10] = ", mean(d10)))
	numerator = mean(y11) - (sum_corrections / n10)
	denominator = mean(d11) - mean(d10)

	wtc = numerator / denominator
	return(wtc)
}



#' @title wald_cic
#' @description  Computes the Wald-CiC estimator
#' @param i Rows of the data frame to operate on
#' @param df Data frame to operate on
#' @param y_name Outcome variable
#' @param g_name Group variable
#' @param t_name Time period variable
#' @param d_name Treatment variable, can be any ordered variable
#' @param X_name Covariates to include
#' @importFrom stats "ecdf" "quantile"
wald_cic = function(i, df, y_name, g_name, t_name, d_name, X_name = NULL){
	df = df[i,]
	if (is.null(X_name)) {
		X_name = "wcicconstwcic"
		df[[X_name]] = 1
	}


	y11 = df[df[[g_name]] == 1 & df[[t_name]] == 1, y_name]
	d11 = df[df[[g_name]] == 1 & df[[t_name]] == 1, d_name]
	n11 = length(y11)
	y10 = df[df[[g_name]] == 1 & df[[t_name]] == 0, y_name]
	d10 = df[df[[g_name]] == 1 & df[[t_name]] == 0, d_name]
	n10 = length(y10)

	Fy_d00 = list()
	Fy_d01 = list()
	Qy_d01 = list()
	for (d in unique(df[[d_name]])) {
		D = as.character(d)
		y_d00 = df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 0, y_name]
		y_d01 = df[df[[d_name]] == d & df[[g_name]] == 0 & df[[t_name]] == 1, y_name]
		Fy_d00[[D]] = ecdf(y_d00)
		Fy_d01[[D]] = ecdf(y_d01)
		Qy_d01[[D]] = function(y) {
			cdf_values = Fy_d00[[D]](y)
			quantiles_d01 = quantile(y_d01, cdf_values, names = FALSE)
			quantiles_d01[cdf_values < min(Fy_d01[[D]](y_d01))] = min(y_d01)
			return(quantiles_d01)
		}
	}

	sum_corrections = 0
	for (i in 1:n10) {
		D = as.character(d10[i])
		Qd = Qy_d01[[D]](y10[i])
		sum_corrections = sum_corrections + Qd
	}
	#print(paste0("E[y11] = ", mean(y11), "; E[Qd] = ", sum_corrections / n10))
	#print(paste0("E[d10] = ", mean(d11), "; E[d10] = ", mean(d10)))
	numerator = mean(y11) - (sum_corrections / n10)
	denominator = mean(d11) - mean(d10)

	wcic = numerator / denominator
	return(wcic)
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
#' @param est Estimators to compute. Should contain at least one of "did", "tc", "cic"
#' @param nboot Number of bootstrap samples for standard errors
#' @importFrom bootstrap "bootstrap"
#' @importFrom stats "quantile" "sd"
#' @export
fuzzydid = function(df, y_name, g_name, t_name, d_name, X_name = NULL, est = c("did", "tc", "cic"), nboot = 50) {
	b = list()
	boot_se = list()
	ci95 = list()
	d_outp = list()

	if ("did" %in% est) {
		print("Computing W_DID...")
		d_outp[["did"]] = d_name
		b[["did"]] = wald_did(1:nrow(df), df, y_name, g_name, t_name, d_name)
		boot = bootstrap(1:nrow(df), 50, wald_did, df, y_name, g_name, t_name, d_name)
		boot_se[["did"]] = sd(boot$thetastar)
		ci95[["did"]] = quantile(boot$thetastar, c(0.05, 0.95))
	}

	if ("tc" %in% est) {
		print("Computing W_TC...")
		d_outp[["tc"]] = d_name
		b[["tc"]] = wald_tc(1:nrow(df), df, y_name, g_name, t_name, d_name)
		boot = bootstrap(1:nrow(df), 50, wald_tc, df, y_name, g_name, t_name, d_name)
		boot_se[["tc"]] = sd(boot$thetastar)
		ci95[["tc"]] = quantile(boot$thetastar, c(0.05, 0.95))
	}

	if ("cic" %in% est) {
		print("Computing W_CIC...")
		d_outp[["cic"]] = d_name
		b[["cic"]] = wald_cic(1:nrow(df), df, y_name, g_name, t_name, d_name)
		boot = bootstrap(1:nrow(df), 50, wald_cic, df, y_name, g_name, t_name, d_name)
		boot_se[["cic"]] = sd(boot$thetastar)
		ci95[["cic"]] = quantile(boot$thetastar, c(0.05, 0.95))
	}
	result = list("d_outp" = d_outp, "b" = b, "boot_se" = boot_se, "ci95" = ci95, n = nrow(df),
	n11 = sum(df[[g_name]] == 1 & df[[t_name]] == 1),
	n10 = sum(df[[g_name]] == 1 & df[[t_name]] == 0),
	n01 = sum(df[[g_name]] == 0 & df[[t_name]] == 1),
	n00 = sum(df[[g_name]] == 0 & df[[t_name]] == 0)
	)
	class(result) = "fuzzydid"
	result
}

#' @title summary.fuzzydid
#' @description Displays the estimations in a table
#' @param object A fuzzydid object
#' @param ... Extra arguments are ignored
#' @importFrom knitr "kable"
#' @exportS3Method Rfuzzydid::summary
summary.fuzzydid = function(object, ...) {
  b = object$b
  boot_se = object$boot_se
  ci95 = object$ci95

  summary_df = data.frame(
    Estimator = paste0("W_", toupper(names(b))),
    Point_Estimate = unlist(b),
    Bootstrap_SE = unlist(boot_se),
    CI_Lower = sapply(ci95, `[[`, 1),
    CI_Upper = sapply(ci95, `[[`, 2)
  )
	rownames(summary_df) = NULL

	kable(summary_df, caption = "Summary of Estimators")
}

#' @title tidy.fuzzydid
#' @description tidy for modelsummary
#' @param x A fuzzydid object
#' @param ... Extra arguments are ignored
#' @importFrom broom "tidy"
#' @export
tidy.fuzzydid = function(x, ...) {
    ret = data.frame(
				model = paste0("W_", toupper(names(x$b))),
				term = paste0("fit_", unlist(x$d_outp)),
        estimate = unlist(x$b),
        std.error = unlist(x$boot_se),
        conf.low = sapply(x$ci95, `[[`, 1),
        conf.high = sapply(x$ci95, `[[`, 2)
    )
    ret
}

#' @title glance.fuzzydid
#' @description glance for modelsummary
#' @param x A fuzzydid object
#' @param ... Extra arguments are ignored
#' @importFrom broom "glance"
#' @export
glance.fuzzydid = function(x, ...) {
	ret = data.frame("Num.Obs." = x$n,
									 "N.11" = x$n11,
	"N.10" = x$n10,
	"N.01" = x$n01,
	"N.00" = x$n00
	)
	ret
}


