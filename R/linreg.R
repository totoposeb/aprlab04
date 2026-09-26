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