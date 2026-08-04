res <- postsim_bh_ui[[1]]$sim_out
alloc_mat   <- res$raw_allocation_age   # (A x T) 각 연령 x 주차에 배포된 도즈
eff_mat     <- res$vacc_delayed         # S에게 간 도즈 (delay 걸리기 전)
wasted_mat  <- res$wasted_dose   

weeks <- 1:ncol(alloc_mat)

df_dose <- tibble(
  week      = weeks,
  effective = colSums(eff_mat),       # 모든 연령 합계
  ineffective = colSums(wasted_mat)   # 모든 연령 합계
)

df_long <- df_dose |>
  pivot_longer(cols = c(effective, ineffective),
               names_to = "status",
               values_to = "doses")

library(dplyr)
library(tidyr)
library(ggplot2)

make_vacc_by_state_df <- function(S, E, I, R, alloc, target_idx = NULL) {
  # target_idx가 주어지면 그 연령만 선택
  if (!is.null(target_idx)) {
    S     <- S[target_idx, ]
    E     <- E[target_idx, ]
    I     <- I[target_idx, ]
    R     <- R[target_idx, ]
    alloc <- alloc[target_idx, ]
  }
  
  T_len <- ncol(alloc)
  out <- vector("list", T_len - 1)
  
  for (t in 2:T_len) {
    alloc_t <- alloc[, t]
    
    S_prev <- S[, t - 1]
    E_prev <- E[, t - 1]
    I_prev <- I[, t - 1]
    R_prev <- R[, t - 1]
    
    N_prev <- S_prev + E_prev + I_prev + R_prev
    
    prop_S <- ifelse(N_prev > 0, S_prev / N_prev, 0)
    prop_E <- ifelse(N_prev > 0, E_prev / N_prev, 0)
    prop_I <- ifelse(N_prev > 0, I_prev / N_prev, 0)
    prop_R <- ifelse(N_prev > 0, R_prev / N_prev, 0)
    
    out[[t - 1]] <- tibble(
      week = t,
      S = sum(alloc_t * prop_S),
      E = sum(alloc_t * prop_E),
      I = sum(alloc_t * prop_I),
      R = sum(alloc_t * prop_R)
    )
  }
  
  bind_rows(out) |>
    pivot_longer(
      cols      = c(S, E, I, R),
      names_to  = "status",
      values_to = "doses"
    )
}
res <- postsim_bh_ui[[1]]$sim_out

df_state <- make_vacc_by_state_df(
  S     = res$S,
  E     = res$E,
  I     = res$I,
  R     = res$R,                 # ⬅️ 여기서도 R 이라고 맞춰줌
  alloc = res$raw_allocation_age # target_idx는 일단 전체 연령 사용
)

ggplot(df_state, aes(x = week, y = doses, fill = status)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Week",
    y = "Number of vaccinated individuals",
    fill = "Infection status at vaccination",
    title = "Vaccinated individuals by infection status at time of vaccination"
  )


df_state_prop <- df_state |>
  group_by(week) |>
  mutate(
    total = sum(doses),
    prop  = if_else(total > 0, doses / total, NA_real_)  # 0/0 방지
  ) |>
  ungroup() |>
  filter(total > 0)   

# Panel A
p1 <- ggplot(df_state_prop, aes(x = week, y = prop, fill = status)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      S = "#0072B2",   # Blue
      E = "#E69F00",   # Orange
      I = "#D55E00",   # Red
      R = "#009E73"    # Dark Green
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "Week",
    y = "Proportion of vaccinated individuals",
    fill = "Infection status at vaccination",
    title = "A. Proportion of vaccinated individuals by infection status"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


# Panel B
p2 <- ggplot(df_state_prop %>% filter(status != "S"),
             aes(week, prop, fill = status)) +
  geom_col(position = "stack") +
  coord_cartesian(ylim = c(0, 0.12)) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c(
    E = "#E69F00",
    I = "#D55E00",
    R = "#009E73"
  )) +
  theme_bw(base_size = 11) +
  labs(
    x = "Week",
    y = "Proportion of vaccinated individuals",
    fill = "Infection status at vaccination",
    title = "B. Proportions for E, I, and R"
  )

p1 + p2 + plot_layout(ncol = 1)
