library(tidyverse)

# ________________________________________________________
# ==================== 1. EXPONENCIAL ====================
# ________________________________________________________

# Lê dados
load("01_data//res_exp_n1000_B1000_ncens100.RData")

# Prepara dados para visualização
dados <- rbind(
  data.frame(c = res_exp$c, bias = res_exp$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_exp$c, bias = res_exp$bias_lambda, tipo = "lambda")
)

# Percentual de Censura x Viés relativo de S(t)
plot_exp <- ggplot(res_exp) +
  aes(x = c, y = bias_S_t) +
  geom_line(colour = "#112446") +
  theme_minimal()
plotly::ggplotly(plot_exp)

# Percentual de Censura x Viés relativo de S(t) e parâmetros
plot_bias <- ggplot(dados) +
  aes(x = c, y = bias, colour = tipo) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()
plotly::ggplotly(plot_bias)



# ________________________________________________________
# ====================== 2. WEIBULL ======================
# ________________________________________________________

# Lê dados
load("01_data//res_wb_n1000_B1000_ncens100.RData")

# Prepara dados para visualização
dados <- rbind(
  data.frame(c = res_wb$c, bias = res_wb$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_wb$c, bias = res_wb$bias_shape, tipo = "shape"),
  data.frame(c = res_wb$c, bias = res_wb$bias_scale, tipo = "scale"))

# Percentual de Censura x Viés relativo de S(t)
plot_wb <- ggplot(res_wb) +
  aes(x = c, y = bias_S_t) +
  geom_line(colour = "#112446") +
  theme_minimal()
plotly::ggplotly(plot_wb)

# Percentual de Censura x Viés relativo de S(t) e parâmetros
plot_bias <- ggplot(dados) +
  aes(x = c, y = bias, colour = tipo) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()
plotly::ggplotly(plot_bias)



# ________________________________________________________
# ==================== 3. LOG-NORMAL =====================
# ________________________________________________________

# Lê dados
load("01_data//res_ln_n1000_B1000_ncens100.RData")

# Prepara dados para visualização
dados <- rbind(
  data.frame(c = res_ln$c, bias = res_ln$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_ln$c, bias = res_ln$bias_meanlog, tipo = "meanlog"),
  data.frame(c = res_ln$c, bias = res_ln$bias_sdlog, tipo = "sdlog"))

# Percentual de Censura x Viés relativo de S(t)
plot_ln <- ggplot(res_ln) +
  aes(x = c, y = bias_S_t) +
  geom_line(colour = "#112446") +
  theme_minimal()
plotly::ggplotly(plot_ln)

# Percentual de Censura x Viés relativo de S(t) e parâmetros
plot_bias <- ggplot(dados) +
  aes(x = c, y = bias, colour = tipo) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()
plotly::ggplotly(plot_bias)


# ________________________________________________________
# ================= 4. GAMA GENERALIZADA =================
# ________________________________________________________

# Lê dados
load("01_data//res_gg_n1000_B1000_ncens100.RData")

# Prepara dados para visualização
dados <- rbind(
  data.frame(c = res_gg$c, bias = res_gg$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_gg$c, bias = res_gg$bias_mu, tipo = "mu"),
  data.frame(c = res_gg$c, bias = res_gg$bias_sigma, tipo = "sigma"),
  data.frame(c = res_gg$c, bias = res_gg$bias_Q, tipo = "Q"))

# Percentual de Censura x Viés relativo de S(t)
plot_gg <- ggplot(res_gg) +
  aes(x = c, y = bias_S_t) +
  geom_line(colour = "#112446") +
  theme_minimal()
plotly::ggplotly(plot_gg)

# Percentual de Censura x Viés relativo de S(t) e parâmetros
plot_bias <- ggplot(dados) +
  aes(x = c, y = bias, colour = tipo) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()
plotly::ggplotly(plot_bias)