violindata <- maindata %>%
  select(SimpleSortPheno, IndexedPheno, AICDA, MKI67, GC_FL_Comp.CXCR4.APC) %>%
  filter(SimpleSortPheno %in% c("GC", "FL"))
#filter(MKI67 > 0)
#filter(AICDA > 0)

MKI67plus <- violindata %>%
  filter(MKI67 > 0) %>%
  drop_na()

MKI67negative <- violindata %>%
  filter(MKI67 == 0) %>%
  drop_na()

ggplot(MKI67plus, aes(x = IndexedPheno, y = GC_FL_Comp.CXCR4.APC)) +
  geom_violin() +
  geom_jitter(aes(color = IndexedPheno), width = 0.15, alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(
    "GC LZ" = "#213891",
    "GC DZ" = "#F50021",
    "GC other" = "#C41ED4",
    "FL" = "#25853B")) +
  labs(y = "CXCR4 (APC)", x = "Phenotype")

ggplot(MKI67negative, aes(x = IndexedPheno, y = GC_FL_Comp.CXCR4.APC)) +
  geom_violin() +
  geom_jitter(aes(color = IndexedPheno), width = 0.15, alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(
    "GC LZ" = "#213891",
    "GC DZ" = "#F50021",
    "GC other" = "#C41ED4",
    "FL" = "#25853B")) +
  labs(y = "CXCR4 (APC)", x = "Phenotype")

scatterdataA <- maindata %>%
  filter(SimpleSortPheno %in% c("FL")) %>%
  filter(AICDA > 0)

ggplot(scatterdataA, aes(x = AICDA, y = GC_FL_Comp.CXCR4.APC)) +
  geom_point()


scatterdataM <- maindata %>%
  filter(SimpleSortPheno %in% c("FL")) %>%
  filter(MKI67 > 0)

ggplot(scatterdataM, aes(x = MKI67, y = GC_FL_Comp.CXCR4.APC)) +
  geom_point()