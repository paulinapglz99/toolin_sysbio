# Función de error para Michaelis-Menten
error_f <- function(params, x, y) {
  c <- params[1]
  a <- params[2]
  b <- params[3]
  y_pred <- f(x, c, a, b)
  sum((y - y_pred)^2)
}

# Función de error para la función tangente
error_g <- function(params, x, y) {
  a <- params[1]
  b <- params[2]
  c <- params[3]
  y_pred <- g(x, a, b, c)
  sum((y - y_pred)^2)
}


# Datos observados
x_vals <- ref1$log10dil
y_vals <- ref1$signal

# Ajuste de Michaelis-Menten
initial_params_f <- c(c = 1, a = 1, b = 1)
fit_f <- optim(par = initial_params_f, fn = error_f, x = x_vals, y = y_vals)
params_f <- fit_f$par  # Parámetros ajustados

# Ajuste de la función tangente
initial_params_g <- c(a = 1, b = 1, c = 1)
fit_g <- optim(par = initial_params_g, fn = error_g, x = x_vals, y = y_vals)
params_g <- fit_g$par  # Parámetros ajustados


# Graficar los puntos originales
plot(x = x_vals, y = y_vals,
     xlab = "Log10 dilution",
     ylab = "Signal",
     pch = 16, cex = 1, type = "b", lwd = 2)

# Agregar ajuste de Michaelis-Menten
y_pred_f <- f(x_vals, params_f[1], params_f[2], params_f[3])
lines(x_vals, y_pred_f, col = "blue", lwd = 2, lty = 2)

# Agregar ajuste de la función tangente
y_pred_g <- g(x_vals, params_g[1], params_g[2], params_g[3])
lines(x_vals, y_pred_g, col = "red", lwd = 2, lty = 3)
