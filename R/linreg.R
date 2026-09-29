#' Linear regression using the QR decomposition
#'
#' @description
#' `linreg` is a reference class (RC) for multiple linear regression. Creating
#' an object with `linreg$new(formula, data)` fits the model: the coefficients
#' are estimated with the QR decomposition of the design matrix, and the
#' results are stored in the object's fields.
#'
#' @details
#' `linreg$new()` takes two arguments:
#' * `formula`: a formula such as `y ~ x1 + x2`. The left side is the
#'   dependent variable; the right side lists the independent variables.
#' * `data`: a data frame containing every variable in `formula`.
#'
#' The design matrix \eqn{X} is split as \eqn{X = QR}. The coefficients solve
#' \eqn{R\hat{\beta} = Q^T y}{R beta = Q'y}, and their variance is
#' \eqn{\hat{\sigma}^2 (R^T R)^{-1}}{sigma^2 (R'R)^(-1)}.
#'
#' @section Methods:
#' * `print()`: prints the call and the coefficients.
#' * `plot()`: draws two diagnostic plots with ggplot2: residuals against
#'   fitted values, and the square root of the absolute standardized
#'   residuals against fitted values (Scale-Location).
#' * `resid()`: returns the vector of residuals.
#' * `pred()`: returns the vector of fitted values.
#' * `coef()`: returns the named vector of coefficients.
#' * `summary()`: prints each coefficient with its standard error, t-value,
#'   p-value and significance stars, followed by the residual standard error
#'   and the degrees of freedom.
#'
#' @field formula The model formula.
#' @field data_name Name of the data frame passed as `data`, used when
#'   printing.
#' @field X Design matrix built with [stats::model.matrix()].
#' @field y Values of the dependent variable.
#' @field beta Named vector of estimated coefficients.
#' @field y_hat Fitted values.
#' @field t_value t-values of the coefficients.
#' @field p_value p-values of the coefficients.
#' @field e Residuals.
#' @field V_e Residual variance.
#' @field dof Degrees of freedom: number of observations minus number of
#'   coefficients.
#' @field se Standard errors of the coefficients.
#'
#' @seealso [stats::lm()], which fits the same model.
#'
#' @examples
#' m <- linreg$new(Petal.Length ~ Sepal.Width + Sepal.Length, data = iris)
#' m$print()
#' m$coef()
#' head(m$pred())
#' head(m$resid())
#' m$summary()
#'
#' m2 <- linreg$new(Petal.Length ~ Species, data = iris)
#' m2$plot()
#'
#' @export linreg
#' @exportClass linreg
linreg <- setRefClass(
  "linreg", # "linreg" is the class name
  fields  = list(
    formula = "formula", data_name = "character", X = "matrix", y = "numeric",
    beta = "numeric", y_hat = "numeric", t_value = "numeric",
    p_value = "numeric", e = "numeric", V_e = "numeric", dof = "numeric",
    se = "numeric"
  ), # stored data, and its type
  methods = list(
    initialize = function(formula, data) { # runs automatically on $new()
      # Validates inputs
      stopifnot(
        "`formula`` is not a formula!" = inherits(formula, "formula"),
        "`data` is not a data.frame!" = is.data.frame(data)
      )

      formula <<- formula
      data_name <<- deparse(substitute(data))

      # Gets data
      X <<- model.matrix(formula, data = data)
      y <<- data[[all.vars(formula)[1]]]

      # Performs QR decomposition
      QR <- qr(X)
      Q <- qr.Q(QR)
      R <- qr.R(QR)

      # Solves for the beta coefficients
      beta <<- drop(backsolve(r = R, x = crossprod(Q, y)))
      names(beta) <<- colnames(X)

      # Calculates fitted values
      y_hat <<- drop(X %*% beta)

      # Calculates residuals
      e <<- y - y_hat

      # Calculates the degrees of freedom
      dof <<- nrow(X) - ncol(X)

      # Calculates the variance of the residuals
      V_e <<-sum(e^2) / dof

      # Calculates the variance of the beta coefficients
      V_beta <- V_e * chol2inv(R)

      # Calculates the standard errors
      se <<- sqrt(diag(V_beta))

      # Calculates the t-values of the beta coefficients
      t_value <<- beta / se

      # Calculates the p-values of the beta coefficients
      p_value <<- 2 * pt(abs(t_value), dof, lower.tail = FALSE)
    },
    print = function() {
      cat("Call:\n")
      cat(
        "linreg(formula = ", deparse1(formula), ", data = ", data_name, ")\n\n",
        sep = ""
      )
      cat("Coefficients:\n")
      base::print(beta)
    },
    plot = function() {
      # Figure 1
      df <- data.frame(fitted = y_hat, resid = e) # columns to plot
      p1 <- ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = resid)) +
        ggplot2::geom_point(shape = 1) + # hollow circles
        ggplot2::stat_summary(
          fun = stats::median, geom = "line", colour = "red"
        ) +
        ggplot2::labs(
          title = "Residuals vs Fitted",
          x = paste0("Fitted values\nlinreg(", deparse1(formula), ")"),
          y = "Residuals"
        )
      base::print(p1)

      # Figure 2
      df <- data.frame(fitted = y_hat, std_resid = sqrt(abs(e / sqrt(V_e))))
      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = std_resid)) +
        ggplot2::geom_point(shape = 1) + # hollow circles
        ggplot2::stat_summary(
          fun = stats::median, geom = "line", colour = "red"
        ) +
        ggplot2::labs(
          title = "Scale-Location",
          x = paste0("Fitted values\nlinreg(", deparse1(formula), ")"),
          y = expression(sqrt(abs("Standardized residuals")))
        )
      base::print(p2)
    },
    resid = function() {
      e
    },
    pred = function() {
      y_hat
    },
    coef = function() {
      beta
    },
    summary = function() {
      # Calculates stars according to p-value
      stars <- ifelse(
        p_value < 0.001, "***", ifelse(
          p_value < 0.01,  "**", ifelse(
            p_value < 0.05,  "*", ifelse(
              p_value < 0.1,   ".", ""
            )
          )
        )
      )

      # Calculates beta coefficients table and prints it
      tab <- data.frame(beta, se, t_value, p_value, stars)
      names(tab) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
      base::print(tab)

      # Prints last line
      cat(
        "\nResidual standard error:", sqrt(V_e), "on", dof, 
        "degrees of freedom\n"
      )

    },
    show = function() {
      print()
    }
  )
)