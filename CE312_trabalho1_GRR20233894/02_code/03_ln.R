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

# Tempos até o evento ~ Distribuição Log Normal
meanlog <- 3
sdlog <- 0.3

# Número de percentuais de censura
n_cens_perc <- 1e2

# Tempos de censura ~ Distribuição Uniforme Contínua
# theta <- seq(1e3, 1e6, length.out = n_cens_perc) # Controle % de censura
calc_cens <- function(theta, n = 20000){
  T_ <- rlnorm(n, meanlog, sdlog)
  C_ <- runif(n, 0, theta)
  mean(C_ < T_)
}


# Gera grid para theta a fim o obter percentuais de censura
# distribuídos de forma quase equilibrada entre 0 e 1.
theta_grid <- exp(seq(log(0.001), log(1e3), length.out = n_cens_perc))

# Calcula percentuais de censura para cada valor de theta
cens_grid <- sapply(theta_grid, calc_cens)

# Organiza dados para obter a função de calibração entre theta e censura
df_calib <- data.frame(
  cens = cens_grid,
  theta = theta_grid
) |>
  arrange(cens) |>
  group_by(cens) |>
  summarise(theta = mean(theta), .groups = "drop")

# Obtém o theta correspondente a cada percentual de censura
# Obs.: Usa-se interpolação linear
theta_from_cens <- function(target_cens){
  approx(
    x = df_calib$cens,
    y = df_calib$theta,
    xout = target_cens,
    rule = 2
  )$y
}

# Obtém o theta linearmante interpolado correspondente
# Obs.: para cada percentual de censura desejado
cens_targets <- seq(0.05, 0.95, length.out = n_cens_perc)
theta <- sapply(cens_targets, theta_from_cens)




# ______________________________________________________________
# =============== 2. GERANDO AMOSTRAS ALEATÓRIAS ===============
# ______________________________________________________________

# Função para gerar amostras aleatórias de dados de sobrevivência com censura
# para a distribuição Gamma Generalizada
aa.ln <- function(n, meanlog, sdlog, theta) {
  
  # Tempos até o evento
  T_ <- rlnorm(n = n, meanlog = meanlog, sdlog = sdlog)
  
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

fit.ln <- function(n = n, meanlog = meanlog, sdlog = sdlog, theta = theta){
  
  # Gera amostra aleatória
  X <- aa.ln(n = n, meanlog = meanlog, sdlog = sdlog, theta = theta)
  
  # Modelo Paramétrico: Distribuição Log-Normal
  # Ajusta parâmetros usando MLE
  fit <- survreg(
    Surv(time, event) ~ 1
    ,data = X
    ,dist = "lognorm"
  )
  
  # Extraindo estimativas pontuais
  meanlog_hat <- fit$coef
  sdlog_hat <- fit$scale
  params <- c(meanlog_hat, sdlog_hat)
  names(params) <- c("meanlog", "sdlog")
  
  # Percentual de censura
  cens_perc <- mean(1 - X$event)
  
  # Retorna as estimativas pontuais
  return(c(c = cens_perc, params))
}


# ________________________________________________________
# ================= 3. CALCULANDO O VIÉS =================
# ________________________________________________________
# Definindo um grid de valores de t
t_ln <- seq(0.05, qlnorm(0.99, meanlog, sdlog), length.out = 1e3)

# Função de sobrevivência teórica
S_t_teorica_ln <- plnorm(q = t_ln, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE)

# Função de sobrevivência estimada
S_t_estimada_ln <- function(t, params){
  plnorm(q = t, meanlog = params["meanlog"], sdlog = params["sdlog"], lower.tail = FALSE)
}

# Calcula viés para cada percentual de censura
sim_ln <- function(B, n, meanlog, sdlog, theta){
  
  params_ams_ln <- replicate(
    B,
    tryCatch(
      fit.ln(n, meanlog, sdlog, theta),
      error = function(e) rep(NA, 4)
    )
  ) |> t()
  
  params_ams_ln <- na.omit(params_ams_ln)
  
  if(nrow(params_ams_ln) == 0) return(c(c = NA, bias = NA))
  
  S_t_estimado_ams_ln <- apply(
    params_ams_ln[, 2:ncol(params_ams_ln)],
    1,
    function(params) S_t_estimada_ln(t_ln, params)
  )
  
  # Calculando o Erro Relativo Integrado (IMRE) para cada coluna (simulação)
  if(ncol(S_t_estimado_ams_ln) == length(t_ln) && nrow(S_t_estimado_ams_ln) != length(t_ln)){
    S_t_estimado_ams_ln <- t(S_t_estimado_ams_ln)
  }
  
  # Métrica: Média do Erro Relativo ao longo da curva (Independente da escala do tempo)
  imre_por_simulacao <- apply(S_t_estimado_ams_ln, 2, function(S_estimada) {
    mean((S_estimada - S_t_teorica_ln) / S_t_teorica_ln, na.rm = TRUE)
  })
  
  bias_S_t <- mean(imre_por_simulacao, na.rm = TRUE)
  bias_meanlog <- mean(params_ams_ln[, "meanlog"]/meanlog - 1, na.rm = TRUE)
  bias_sdlog <- mean(params_ams_ln[, "sdlog"]/sdlog - 1, na.rm = TRUE)
  cens_perc <- mean(params_ams_ln[, "c"], na.rm = TRUE)
  
  return(c(c = cens_perc,
           bias_S_t = bias_S_t,
           bias_meanlog = bias_meanlog,
           bias_sdlog = bias_sdlog))
}

# Executa a simulação para cada percentual de censura
res_ln <- future_sapply(theta, function(theta) {
  sim_ln(B, n, meanlog, sdlog, theta)},
  future.seed = TRUE)

# Salva dados em .RData
res_ln <- as.data.frame(t(res_ln))
path <- paste0("01_data//res_ln_", "n", n, "_B", B, "_ncens", n_cens_perc, ".RData")
save(res_ln, file = path)
colnames(res_ln) <- c("c", "bias_S_t", "bias_meanlog", "bias_sdlog")

# Organiza dados para visualização
dados <- rbind(
  data.frame(c = res_ln$c, bias = res_ln$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_ln$c, bias = res_ln$bias_meanlog, tipo = "meanlog"),
  data.frame(c = res_ln$c, bias = res_ln$bias_sdlog, tipo = "sdlog"))