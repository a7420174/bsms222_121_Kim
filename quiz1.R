m <- readRDS("data_BrainSpan.rds")
e <- m$e1
s <- m$s %>% mutate(col_id = paste(structure_acronym, donor_id, gender, age, sep='_'))
g <- m$g
colnames(e) <- s$col_id
e$gene_symbol <- g$gene_symbol
devstage <- function(age){
  case_when(grepl('pcw', age) ~ as.numeric(gsub(' pcw', '', age)),
            grepl('mos', age) ~ as.numeric(gsub(' mos', '', age)) * 30 + 365,
            grepl('yrs', age) ~ as.numeric(gsub(' yrs', '', age)) * 365 + 365)
}

# Q1
e %>% gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, structure_acronym)) %>% 
  ggplot(aes(structure_acronym, expr)) + geom_boxplot()

# Q2
merge(g %>% select(gene_symbol,gene_type), e) %>% 
  gather(col_id, expr, -gene_symbol, -gene_type) %>% 
  ggplot(aes(gene_type, expr)) + geom_boxplot()

# Q3
e1 <- e %>% filter(gene_symbol == "SCN2A") %>% 
  gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, age, region = structure_acronym))
e1$age2 <- devstage(e1$age)
e1 %>% ggplot(aes(age2, expr)) + geom_point() + scale_x_log10()

# Q4
e1 %>% ggplot(aes(age2, expr, col = region)) + 
  geom_point() + scale_x_log10()

# Q5
e2 <- e %>% filter(grepl('^SCN', gene_symbol)) %>% 
  gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, age))
e2$age2 <- devstage(e2$age)
e2 %>% ggplot(aes(age2, expr)) + geom_point() + scale_x_log10()

# Q6
e3 <- e %>% filter(grepl('^KCN', gene_symbol)) %>% 
  gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, age))
e3$age2 <- devstage(e3$age)
e3 %>% ggplot(aes(age2, expr)) + geom_point() + scale_x_log10()

# Q7
e4 <- e %>% filter(gene_symbol == "XIST") %>% 
  gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, age, gender))
e4$age2 <- devstage(e4$age)
e4 %>% ggplot(aes(age2, expr, col = gender)) + 
  geom_point() + scale_x_log10()

# Q8
top3 <- e %>% gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, gender)) %>% 
  group_by(gene_symbol, gender) %>% summarise(avg = mean(expr)) %>% 
  spread(gender, avg) %>% mutate(diff = abs(`M`-`F`)) %>% 
  arrange(desc(diff)) %>% head(3)
e5 <- e %>% filter(gene_symbol %in% top3$gene_symbol) %>% 
  gather(col_id, expr, -gene_symbol) %>% 
  merge(s %>% select(col_id, age, gender))
e5$age2 <- devstage(e5$age)
e5 %>% ggplot(aes(age2, expr, col = gender)) + 
  geom_point() + scale_x_log10() + 
  facet_wrap(~gene_symbol,nrow = 1)
