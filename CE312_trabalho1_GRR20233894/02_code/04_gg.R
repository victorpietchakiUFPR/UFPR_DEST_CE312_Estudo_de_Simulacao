# ESTUDO APLICADO DO MÉTODO DE SIMULAÇÃO DE MONTE CARLO
# ________________________________________________________
# AVALIANDO O VIÉS DE ESTIMADORES DA FUNÇÃO DE SOBREVIVÊNCIA
# EM AMOSTRAS DE DIFERENTES PERCENTUAIS DE CENSURA:
# UMA APLICAÇÃO DO MÉTODO DE MONTE CARLO

# ________________________________________________________
# ================ 0. PACOTES NECESSÁRIOS ================
# ________________________________________________________

# Carregando pacotes...
library(tidyverse)
library(survival)
library(flexsurv)
# install.packages("future.apply")
library(future.apply)
library(future)
plan(multisession, workers = 2)




# ______________________________________________________________
# =========== 1. DEFINIÇÃO DE PARÂMETROS DA SIMULAÇÃO ==========
# ______________________________________________________________

# Tamanho amostral
n <- 1e3

# Número de simulações
B <- 1e3

# Tempos até o evento ~ Distribuição Gamma Generalizada
mu <- log(0.3)
sigma <- 1/3
Q <- 3

# Número de percentuais de censura
n_cens_perc <- 1e2

# Tempos de censura ~ Distribuição Uniforme Contínua
# theta <- seq(1e3, 1e6, length.out = n_cens_perc) # Controle % de censura
calc_cens <- function(theta, n = 20000){
  T_ <- rgengamma(n, mu, sigma, Q)
  C_ <- runif(n, 0, theta)
  mean(C_ < T_)
}

# Gera grid para theta a fim o obter percentuais de censura
# distribuídos de forma quase equilibrada entre 0 e 1.
theta_grid <- exp(seq(log(0.01), log(1e2), length.out = n_cens_perc))
cens_grid <- sapply(theta_grid, calc_cens)
df_calib <- data.frame(
  cens = cens_grid,
  theta = theta_grid
) |>
  arrange(cens) |>
  group_by(cens) |>
  summarise(theta = mean(theta), .groups = "drop")
theta_from_cens <- function(target_cens){
  approx(
    x = df_calib$cens,
    y = df_calib$theta,
    xout = target_cens,
    rule = 2
  )$y
}
cens_targets <- seq(0.05, 0.95, length.out = n_cens_perc)
theta <- sapply(cens_targets, theta_from_cens)



# ______________________________________________________________
# =============== 2. GERANDO AMOSTRAS ALEATÓRIAS ===============
# ______________________________________________________________

# Função para gerar amostras aleatórias de dados de sobrevivência com censura
# para a distribuição Gamma Generalizada
aa.gengamma <- function(n, mu, sigma, Q, theta) {
  
  # Tempos até o evento
  T_ <- rgengamma(n = n, mu = mu, sigma = sigma, Q = Q)
  
  # Tempos de censura
  C_ <- runif(n, min = 0, max = theta)
  
  # Tempos Observados
  tempos <- pmin(T_, C_)
  # Indicador de censura
  cens <- as.integer(T_ <= C_)
  # Tempos observados com indicador de censura
  X <- cbind(tempos, cens) |> as.data.frame()
  colnames(X) <- c("time", "event")
  
  return(X)
}


# ________________________________________________________
# ================= 3. ESTIMAÇÃO DE S(t) =================
# ________________________________________________________

fit.gengamma <- function(n = n, mu = mu, sigma = sigma, Q = Q, theta = theta){
  
  # Gera amostra aleatória
  X <- aa.gengamma(n = n, mu = mu, sigma = sigma, Q = Q, theta = theta)
  
  # Modelo Paraamétrico: Distribuição Gamma Generalizada
  # Ajusta parâmetros usando MLE
  fit <- flexsurvreg(
    Surv(time, event) ~ 1
    ,data = X
    ,dist = "gengamma"
  )
  
  # Extraindo estimativas pontuais
  params_fit <- fit$res
  mu_hat    <- params_fit["mu"   , "est"]
  sigma_hat <- params_fit["sigma", "est"]
  Q_hat     <- params_fit["Q"    , "est"]
  params <- c(mu_hat, sigma_hat, Q_hat)
  names(params) <- c("mu", "sigma", "Q")
  
  # Percentual de censura
  cens_perc <- mean(1 - X$event)
  
  # Retorna as estimativas pontuais
  return(c(c = cens_perc, params))
}


# ________________________________________________________
# ================= 3. CALCULANDO O VIÉS =================
# ________________________________________________________
# Definindo um grid de valores de t
t_gg <- seq(0.05, quantile(rgengamma(1e5, mu, sigma, Q), 0.90), length.out = 1e4)

# Função de sobrevivência teórica
S_t_teorica_gg <- pgengamma(q = t_gg, mu = mu, sigma = sigma, Q = Q, lower.tail = FALSE)

# Função de sobrevivência estimada
S_t_estimada_gg <- function(t, params){
  pgengamma(q = t, mu = params["mu"], sigma = params["sigma"], Q = params["Q"], lower.tail = FALSE)
}

# Calcula viés para cada percentual de censura
sim_gg <- function(B, n, mu, sigma, Q, theta){
  
  params_ams_gg <- replicate(
    B,
    tryCatch(
      fit.gengamma(n, mu, sigma, Q, theta),
      error = function(e) rep(NA, 4)
    )
  ) |> t()
  
  params_ams_gg <- na.omit(params_ams_gg)
  
  if(nrow(params_ams_gg) == 0) return(c(c = NA, bias = NA))
  
  S_t_estimado_ams_gg <- apply(
    params_ams_gg[, 2:ncol(params_ams_gg)],
    1,
    function(params) S_t_estimada_gg(t_gg, params)
  )
  
  # Calculando o Erro Relativo Integrado (IMRE) para cada coluna (simulação)
  if(ncol(S_t_estimado_ams_gg) == length(t_gg) && nrow(S_t_estimado_ams_gg) != length(t_gg)){
    S_t_estimado_ams_gg <- t(S_t_estimado_ams_gg)
  }
  
  # Métrica: Média do Erro Relativo ao longo da curva (Independente da escala do tempo)
  imre_por_simulacao <- apply(S_t_estimado_ams_gg, 2, function(S_estimada) {
    mean((S_estimada - S_t_teorica_gg) / S_t_teorica_gg, na.rm = TRUE)
  })
  
  bias_S_t <- mean(imre_por_simulacao, na.rm = TRUE)
  bias_mu <- mean(params_ams_gg[, "mu"]/mu - 1, na.rm = TRUE)
  bias_sigma <- mean(params_ams_gg[, "sigma"]/sigma - 1, na.rm = TRUE)
  bias_Q <- mean(params_ams_gg[, "Q"]/Q - 1, na.rm = TRUE)
  cens_perc <- mean(params_ams_gg[, "c"], na.rm = TRUE)
  
  return(c(c = cens_perc,
           bias_S_t = bias_S_t,
           bias_mu = bias_mu,
           bias_sigma = bias_sigma,
           bias_Q = bias_Q))
}

# Executa a simulação para cada percentual de censura
res_gg <- future_sapply(theta, function(theta) {
  sim_gg(B, n, mu, sigma, Q, theta)
})

# Salva dados em .RData
res_gg <- as.data.frame(t(res_gg))
path <- paste0("01_data//res_gg_", "n", n, "_B", B, "_ncens", n_cens_perc, ".RData")
save(res_gg, file = path)
colnames(res_gg) <- c("c", "bias_S_t", "bias_mu", "bias_sigma", "bias_Q")

# Organiza dados para visualização
dados <- rbind(
  data.frame(c = res_gg$c, bias = res_gg$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_gg$c, bias = res_gg$bias_mu, tipo = "lambda"),
  data.frame(c = res_gg$c, bias = res_gg$bias_sigma, tipo = "sigma"),
  data.frame(c = res_gg$c, bias = res_gg$bias_Q, tipo = "Q")
)