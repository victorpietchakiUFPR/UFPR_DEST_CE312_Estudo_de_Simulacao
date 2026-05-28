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
library(future)
library(future.apply)

plan(multisession, workers = 4)


# ______________________________________________________________
# =========== 1. DEFINIÇÃO DE PARÂMETROS DA SIMULAÇÃO ==========
# ______________________________________________________________

# Tamanho amostral
n <- 1e3

# Número de simulações
B <- 1e3

# Número de percentuais de censura
n_cens_perc <- 1e2

# Tempos até o evento ~ Distribuição Exponencial
lambda <- 0.3

# # Tempos de censura ~ Distribuição Uniforme Contínua
calc_cens <- function(theta, n = 20000){
  T_ <- rexp(n, rate = 1/lambda)
  C_ <- runif(n, 0, theta)
  mean(C_ < T_)
}

# Gera grid para theta a fim o obter percentuais de censura
# distribuídos de forma quase equilibrada entre 0 e 1.
theta_grid <- exp(seq(log(0.001), log(1e3), length.out = n_cens_perc))
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
aa.exp <- function(n, lambda, theta) {
  
  # Tempos até o evento
  T_ <- rexp(n = n, rate = 1/lambda)
  
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

fit.exp <- function(n = n, lambda = lambda, theta = theta){
  
  # Gera amostra aleatória
  X <- aa.exp(n = n, lambda = lambda, theta = theta)
  
  # Modelo Paraamétrico: Distribuição Gamma Generalizada
  # Ajusta parâmetros usando MLE
  fit <- flexsurvreg(
    Surv(time, event) ~ 1
    ,data = X
    ,dist = "exponential"
  )
  
  # Extraindo estimativas pontuais
  lambda_hat <- 1/exp(fit$coefficients[1])
  
  # lambda_hat <- sum(X$event) / sum(X$time)
  
  # Percentual de censura
  cens_perc <- mean(1 - X$event)
  
  # Retorna as estimativas pontuais
  return(c(c = cens_perc, lambda = lambda_hat))
}


# ________________________________________________________
# ================= 3. CALCULANDO O VIÉS =================
# ________________________________________________________
# Definindo um grid de valores de t
t_exp <- seq(0, quantile(rexp(1e5, rate = 1/lambda), 0.99), length.out = 1e4)

# Função de sobrevivência teórica
S_t_teorica_exp <- pexp(q = t_exp, rate = 1/lambda, lower.tail = FALSE)

# Função de sobrevivência estimada
S_t_estimada_exp <- function(t, lambda){
  pexp(q = t, rate = 1/lambda, lower.tail = FALSE)
}

# Calcula viés para cada percentual de censura
sim_exp <- function(B, n, lambda, theta){
  
  # Gerando amostras e estimando
  params_ams_exp <- replicate(
    B,
    fit.exp(
      n = n,
      lambda = lambda,
      theta = theta
    )
  ) |> t()
  
  # Estimando S(t) para cada simulação
  S_t_estimado_ams_exp <- sapply(
    params_ams_exp[, "lambda"]
    ,function(lambda) S_t_estimada_exp(t_exp, lambda)
  )
  S_t_estimado_ams_mean_exp <- rowMeans(S_t_estimado_ams_exp)
  
  # Estimando o viés
  cens_perc <- mean(params_ams_exp[, "c"])
  bias_S_t <- mean(S_t_estimado_ams_mean_exp/S_t_teorica_exp - 1)
  bias_lambda <- mean(params_ams_exp[, "lambda"]/lambda - 1)
  
  return(c(c = cens_perc, bias_S_t = bias_S_t, bias_lambda = bias_lambda))
  
}

# Executa a simulação para cada percentual de censura
res_exp <- future_sapply(
  theta,
  function(x) {
    sim_exp(
      B = B,
      n = n,
      lambda = lambda,
      theta = x
    )
  },
  future.seed = TRUE
)

# Salva dados em .RData
res_exp <- as.data.frame(t(res_exp))
path <- paste0("01_data//res_exp_", "n", n, "_B", B, "_ncens", n_cens_perc, ".RData")
save(res_exp, file = path)

# Organiza dados para visualização
dados <- rbind(
  data.frame(c = res_exp$c, bias = res_exp$bias_S_t, tipo = "S(t)"),
  data.frame(c = res_exp$c, bias = res_exp$bias_lambda, tipo = "lambda")
)