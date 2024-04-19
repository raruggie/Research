# Run NWIS.R first

# make plot of boxplots of C for obsevaed CQ pairs and predicted C for each site

x <- df.predC %>% select(Name, Date, X_00060_00003, pred.C)

y <- df.RP %>% select(Name, sample_dt, Q, C)  %>% rename(Date = 2)

z <- left_join(y, x, by = c('Name', 'Date')) %>% select(-5) %>% 
  pivot_longer(cols = c(C, pred.C))

ggplot(z, aes(x = Q, y = value, color = name))+
  geom_point()+
  facet_wrap(~Name, scales = 'free')+
  scale_x_log10()+
  scale_y_log10()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
























