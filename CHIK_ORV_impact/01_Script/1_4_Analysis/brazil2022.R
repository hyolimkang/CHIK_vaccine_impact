totcases_brazil <- read_excel("00_Data/0_1_Raw/epi_report_plisa.xlsx", 
                               sheet = "totcase_time_series")

totcases_brazil$cases <- as.numeric(totcases_brazil$cases)
totcases_brazil$cases[is.na(totcases_brazil$cases)] <- 0
totcases_brazil$week <- 1:nrow(totcases_brazil)

totcases_brazil_21 <- totcases_brazil %>% filter(year == 2021)
totcases_brazil_22 <- totcases_brazil %>% filter(year == 2022)
totcases_brazil_2122 <- totcases_brazil[417:460, ]

ggplot(totcases_brazil)+
  geom_line(aes(x = week, y = cases, color = as.factor(year))) +
  labs(color = "Year")+
  theme_minimal()

ggplot(totcases_brazil_2122)+
  geom_line(aes(x = week, y = cases))+
  theme_minimal()