startd <- function(d){
  str_extract(d, "\\d+")
}

startp <- function(p){
  str_extract(p, "\\d+")
}

d1 <- d[d$Inheritence %in% c("DeNovo", "Inherited") & d$Effect == "Missense",]
dna <- as.numeric(startd(d1$c.DNA))
pro <- as.numeric(startp(d1$p.Protein))

plot(dna, pro)

table(d$Effect)
max(as.numeric(d$SeizureOnsetDays[d$Effect == "Frameshift"]), na.rm = TRUE)
max(as.numeric(d$SeizureOnsetDays[d$Effect == "Missense"]), na.rm = TRUE)
max(as.numeric(d$SeizureOnsetDays[d$Effect == "Nonsense"]), na.rm = TRUE)

d %>% group_by(Effect) %>%
summarise(max = max(as.numeric(gsub('[^0-9]','',SeizureOnsetDays)), na.rm = TRUE))

d1 <- d[d$Effect == "Missense",]
cDNA <- as.numeric(gsub('[^0-9]','',d1$c.DNA))
protein <- as.numeric(gsub('[^0-9]','',d1$p.Protein))
plot(cDNA,protein)

## Plot the relationship between seizure onset days and age at assessment.
S_Day <- as.numeric(str_extract(d$SeizureOnsetDays, '[0-9]+'))
Age_at_Ass <- ifelse(str_detect(d$PatientAgeAtAssessment, 'm'),
                     as.numeric(str_extract(d$PatientAgeAtAssessment, "[0-9]+"))/12,
                     as.numeric(str_extract(d$PatientAgeAtAssessment, "[0-9]+")))
plot(S_Day,Age_at_Ass, log = "xy")

