
# load data:



#### Hypothesis testing ####

# set up list of dataframes for response variables and predictors

# OLS slope
# OLS slope stormflow samples
# AANY.hydrosep
# medC

resp.vars <- c("OLS.Slope.1s", "OLS.Slope.2s_hydrosep_post", "OLS.AANY.2s.hydrosep.method2", "medC")

# set up dataframes for the 4 response variables::

l.resp.allpred <- lapply(resp.vars, \(i) df.PCA %>% filter(Name %in% df.g$Name) %>% select(i, v.pred.cols) %>% rename(term = 1)) %>% purrr::set_names(resp.vars)

# lets look at the distributions of the response variables:

l.resp.allpred[[1]] %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(., aes(value))+
  geom_dotplot()+
  facet_wrap(~name, scales = 'free')

# lets look at histogram of term colored by percent forest:

l.resp.allpred[[2]] %>% 
  ggplot(., aes(term))+
  geom_histogram(aes(fill = R_FORESTNLCD06, group = R_FORESTNLCD06))

# using the NYS data, test the following hypotheses:

#### ~ H1 ####

# H1: Unless nearly 100% forested, most catchments display nutrient mobilization as Q increases

# what we see with the histograms:

bind_rows(l.resp.allpred[[1]] %>% pivot_longer(-term) %>% mutate(resp = 'Overall Slope', .before = 1), l.resp.allpred[[2]] %>% pivot_longer(-term) %>% mutate(resp = 'Stormflow Slope', .before = 1)) %>% 
  filter(name == 'R_FORESTNLCD06') %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  ggplot(., aes(term))+
  geom_histogram(aes(fill = R_FORESTNLCD06, group = R_FORESTNLCD06))+
  facet_wrap(~resp)+
  ylab('Number of Sites')+
  xlab('Slope Value')

# is that lower forest does not equate with more slope...

#### ~ H2 ####

# H2: Aggregate measures of the fraction developed land (ie. total ag, total urban) provide moderate explanation of watershed variations in nutrient load

# pick just these to try to predict yield and medC:

l.H2 <- lapply(l.resp.allpred[3:4], \(i) i %>% select(term, R_PLANTNLCD06, R_DEVNLCD06) %>% lm(term ~ ., data = .))

summary(l.H2[[1]]) # none are significant in predicting the yield
summary(l.H2[[2]]) # agriculture is significant at 5% and developed is signficant at 10% for predicting medC

# lets try PLSR:

l.H2 <- lapply(l.resp.allpred[3:4], \(i) i %>% select(term, R_PLANTNLCD06, R_DEVNLCD06) %>% pls::plsr(term ~ ., data = ., scale = T, validation = 'CV',segments = 10,jackknife = TRUE))

varImp(l.H2[[1]])
varImp(l.H2[[2]])

RMSEP(l.H2[[1]])
RMSEP(l.H2[[2]])

R2(l.H2[[1]])
R2(l.H2[[2]])

#### ~ H3 ####

# H3: Incorporating measures of the fraction of developed land in close proximity to streams (i.e. land use within stream buffers) improves prediction of watershed variations in nutrient load

# pick just these to try to predict yield and medC:

l.H3 <- lapply(l.resp.allpred[3:4], \(i) i %>% select(term, RIP.CSA.100, R_RIP100_DEV, R_RIP100_PLANT) %>% lm(term ~ ., data = .))

summary(l.H3[[1]]) # none are signficant in predicting the yield
summary(l.H3[[2]]) # agriculture is significant at 10% for predicting medC

# lets try plsr:

l.H3 <- lapply(l.resp.allpred[3:4], \(i) i %>% select(term, RIP.CSA.100, R_RIP100_DEV, R_RIP100_PLANT) %>% pls::plsr(term ~ ., data = ., scale = T, validation = 'CV',segments = 10,jackknife = TRUE))

varImp(l.H3[[1]])
varImp(l.H3[[2]])

RMSEP(l.H3[[1]])
RMSEP(l.H3[[2]])

R2(l.H3[[1]])
R2(l.H3[[2]])

# lets compare l.H2 and l.H3:

RMSEP(l.H2[[2]])
RMSEP(l.H3[[2]])

# RMSEP did not go down