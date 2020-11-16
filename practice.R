rbind.data.frame(d_luad %>% 
                   count(tumor_stage) %>%
                   mutate(Type = 'LUAD'),
                 d_lusc %>% 
                   count(tumor_stage) %>%
                   mutate(Type = 'LUSC')) %>%
  ggplot(aes(tumor_stage, n, fill=tumor_stage)) + 
  geom_bar(stat="identity", position = position_dodge(), show.legend = F) + 
  labs(x = 'Short Letter Code', 
       y = 'Sample Count', 
       title = 'Number of Short Letter Code in TCGA-LUAD dataset') +
  facet_wrap(~ Type, nrow = 2) + coord_flip() + scale_fill_manual(values = rev(c(colfunc(10),'grey')))
rev
colfunc<-colorRampPalette(c("red","blue"))

d_luad %>% ggplot(aes(gender, fill = gender)) + geom_bar()

rbind.data.frame(d_luad %>% count(gender) %>% mutate(type = 'luad'),
                 d_lusc %>% count(gender) %>% mutate(type = 'lusc')) %>% 
  ggplot(aes(gender, n, fill = gender)) + 
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(~ type)

rbind.data.frame(d_luad %>% select(gender, cigarettes_per_day) %>% mutate(type = 'laud'),
                 d_lusc %>% select(gender, cigarettes_per_day) %>% mutate(type = 'lusc')) %>% 
  ggplot(aes(gender, cigarettes_per_day, color = type)) + geom_point() + 
  labs(x = 'sex', y = 'smoking days')

quantile(d_lusc$age_at_diagnosis/365, na.rm = T)
d_lusc %>% ggplot(aes(age_at_diagnosis, fill = gender)) + geom_histogram(bins = 100)
a <- d_lusc %>% ggplot(aes(age_at_diagnosis, fill = gender)) + geom_density(alpha = 0.2)
b <- d_luad %>% ggplot(aes(age_at_diagnosis, fill = gender)) + geom_density(alpha = 0.2)
plot_grid(a, b)

p <- ggplot(d_luad, aes(gender, age_at_diagnosis/365, fill=gender)) + 
  geom_violin()

m <- readRDS("data_BrainSpan.rds")
n <- read.delim("non_alt_loci_set.txt")
n1 <- merge(m$g, n %>% select(name, ensembl_gene_id),
            by = "ensembl_gene_id", all.x = T)
m$s %>% ggplot(aes(structure_acronym)) + geom_boxplot()

cbind.data.frame(m$g %>% select(gene_type), m$e1) %>% 
  gather(sample, expression, -gene_type) %>% 
  ggplot(aes(gene_type, expression)) + geom_boxplot() + scale_y_log10()

MFC <- m$s %>% filter(structure_acronym == "MFC") %>% .$column_num %>% paste0("V", .)
MD <- m$s %>% filter(structure_acronym == "MD") %>% .$column_num %>% paste0("V", .)
m$e %>% gather(sample, expression) %>% 
  mutate(region = case_when(sample %in% MFC ~ "MFC", sample %in% MD ~ "MD")) %>% 
  ggplot(aes(region, expression)) + geom_boxplot() + scale_y_log10()

m1 <- cbind.data.frame(m$g, m$e)
data.frame(SCN2A = as.numeric(m1 %>% filter(gene_symbol == "SCN2A") %>% select(V9:V513)), stage = m$s$age) %>% 
  ggplot(aes(fct_inorder(stage), SCN2A)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

data.frame(SCN2A = ,stage = m$s$age) %>% 
  ggplot(aes(fct_inorder(stage), SCN2A)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

S4 <- read_excel("aav8130_Data_S4.xls")
View(S4)
S4 %>% ggplot(aes(`Fold change`, -log10(`q value`))) + geom_point()
S4_4 %>% mutate(sum = rowSums(.[-(1:6)])) %>% filter(`Cell type` == 'OPC')
df %>% ggplot(aes(DEGs_num)) + geom_histogram(binwidth = 10)
df %>% ggplot(aes(A, DEGs_num)) + geom_point()
library(pscl)
summary(p1 <- glm(DEGs_num ~ A, family = poisson, data = df %>% filter(`Cell type` == 'L2/3')))

5144 5978
S1_2 <- read_excel("aav8130_Data_S1.xlsx", sheet = 2)
S1_2 %>% mutate(Individual = as.character(Individual)) %>% filter(Individual == c("5144", "5978")) %>%
  ggplot(aes(Bulk_reads, `Sample name`)) + geom_point()
S4_4 %>% select(`ASD patients with DE signal`) %>% separate_rows(`ASD patients with DE signal`, sep = ",") %>% View()
separate_rows()
family <- read_delim("family.txt", delim = "\t")
test <- left_join(S4_4, family %>% select(`Gene group name`, "Gene name" = `Approved symbol`))
S4_4 %>% count(`Gene ID`) %>% filter(n>1)
S4_4 %>% nrow()
test %>% nrow()
test %>% count(`Gene group name`)

S4_1 %>% filter(`Gene name` %in% c("PRRC2B","ZMYND8","C1orf115","KPNB1","CLMP","VWA8","AC012501.2","NDUFAF6","SNX32","TXNDC16"))
S4_4 %>% filter(`Gene name` %in% c("PRRC2B","ZMYND8","C1orf115","KPNB1","CLMP","VWA8","AC012501.2","NDUFAF6","SNX32","TXNDC16"))

S2_1 <- read_excel("aav8130_Data_S2.xlsx")
S2_1 %>% count(individual)

test <- read.delim("phs000126.pha002865.txt.gz", sep = "\t", comment.char = "#")
