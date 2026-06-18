# Carregando os pacotes necessários 

library("terra")
library("sf")
library("landscapemetrics")
library("tmap")
library("ggplot2")
library("vegan")
library("tidyverse")
library("rgbif")
library("geodata")
library("geobr")
library("tidyterra")
library("dplyr")

cat("✓ Todos os pacotes carregados com sucesso!\n")

#Estabelecendo pasta de input e output

setwd("C:/faculdade/ProvaEcopais")

# Carregando o Raster de Bom Jesus do Tocantins (PA)
mun_2024 <- rast("toca_2024.tif")

# Buscando dados do IBGE
geobr::read_municipality(1501576)

# Configurando coordenadas (Lat-long --> metros)
mun_2024_utm <- project(mun_2024, "EPSG:31982", method = "near")

# Verficando configuração
print(mun_2024_utm)

# Identificando 
unique(mun_2024_utm)

# criando legendas e definindo cores
legenda <- data.frame(
  value = c(3,4,6,9,11,12,15,24,25,29,31,33,39,41),
  label = c(
    "Formação Florestal",
    "Formação Savânica",
    "Floresta Alagável",
    "Floresta Plantada",
    "Campo Alagado e Área Pantanosa",
    "Formação Campestre",
    "Pastagem",
    "Área Urbanizada",
    "Outras Áreas Não Vegetadas",
    "Afloramento Rochoso",
    "Aquicultura",
    "Rio, Lago e Oceano",
    "Soja",
    "Outras Lavouras Temporárias"
  )
)
cores <- c(
  "#1f8d49", # 3  Formação Florestal
  "#7dc975", # 4  Formação Savânica
  "#04381d", # 6  Floresta Alagável
  "#7a5900", # 9  Floresta Plantada
  "#519799", # 11 Campo Alagado e Área Pantanosa
  "#d6bc74", # 12 Formação Campestre
  "#edde8e", # 15 Pastagem
  "#d4271e", # 24 Área Urbanizada
  "#db4d4f", # 25 Outras Áreas Não Vegetadas
  "#ffaa5f", # 29 Afloramento Rochoso
  "#9c0027", # 31 Aquicultura
  "#2532e4", # 33 Rio, Lago e Oceano
  "#f5b3c8", # 39 Soja
  "#ffefc3"  # 41 Outras Lavouras Temporárias
)
legenda$Cor <- cores
levels(mun_2024_utm) <- legenda[,1:2]

# Vendo o resultado 
ggplot() +
  geom_spatraster(data = mun_2024_utm) +
  scale_fill_manual(values = cores, na.value = "transparent") +
  labs(title = "Bom Jesus do Tocantins 2024") +
  theme_minimal()

# Calculando métricas paisagem
check_landscape(mun_2024_utm)

# Métricas para Landscape (paisagem)
# diversidade 
shannon <- lsm_l_shdi(mun_2024_utm)
head(shannon)
# Área total
area_total <- lsm_l_ta(mun_2024_utm)
print(area_total)

# Métricas para class (classe) 
# Proporção de cada classe
proporcao <- lsm_c_pland(mun_2024_utm)
print(proporcao)
### Criando tabela 
legendatab <- c(
  "3"  = "Formação Florestal",
  "4"  = "Formação Savânica",
  "6"  = "Floresta Alagável",
  "9"  = "Floresta Plantada",
  "11" = "Campo Alagado e Área Pantanosa",
  "12" = "Formação Campestre",
  "15" = "Pastagem",
  "24" = "Área Urbanizada",
  "25" = "Outras Áreas Não Vegetadas",
  "29" = "Afloramento Rochoso",
  "31" = "Aquicultura",
  "33" = "Rio, Lago e Oceano",
  "39" = "Soja",
  "41" = "Outras Lavouras Temporárias"
)

proporcao_tab <- proporcao %>%
  mutate(classe_nome = legendatab[as.character(class)]) %>%
  select(classe_nome, value) %>%
  rename(
    Classe = classe_nome,
    PLAND = value
  ) %>%
  arrange(desc(PLAND))

proporcao_tab

# Borda por CLASSE
borda <- lsm_c_te(mun_2024_utm)
print(borda)

#criando tabela
borda_tab <- borda %>%
  mutate(
    classe_nome = legendatab[as.character(class)]
  ) %>%
  select(
    classe_nome,
    value
  ) %>%
  rename(
    Classe = classe_nome,
    Borda_m = value
  ) %>%
  arrange(desc(Borda_m))

borda_tab

# Métricas para Patch (Nivel de Fragmentos)
# Área dos fragmentos
areas <- lsm_p_area(mun_2024_utm)
head(areas)
### Criando tabela 

areas %>%
  mutate(
    classe_nome = "Formação Florestal"
  ) %>%
  select(
    classe_nome,
    id,
    value
  ) %>%
  rename(
    Classe = classe_nome,
    Fragmento = id,
    Area = value
  )

# Métricas para classe florestal
areas %>%
  filter(class == 3) %>%
  summarise(
    area_total = sum(value),
    n_fragmentos = n(),
    maior_fragmento = max(value)
  )

# Distância para o vizinho mais próximo
isolamento <- lsm_p_enn(mun_2024_utm)
head(isolamento)

# Realizando o sorteio de pontos
# 15 pontos em 2024
# Definindo área de sorteio 
vetor_total <- as.polygons(mun_2024_utm, dissolve = TRUE)
municipio_limite <- aggregate(vetor_total)
limite_seguro <- buffer(municipio_limite, width = -2000)

floresta <- vetor_total[vetor_total$label == "Formação Florestal", ]
area_sorteio <- crop(floresta, limite_seguro)

if (is.list(area_sorteio)) {
  area_sorteio <- do.call(rbind, area_sorteio)
}

if (nrow(area_sorteio) == 0) {
  stop("Erro: Nenhuma floresta encontrada a mais de 1km da borda.")
}

# Definindo distância minima de 2000m e sorteando pontos

set.seed(1)
pontos_finais <- NULL
tentativas <- 0

while(is.null(pontos_finais) || nrow(pontos_finais) < 15) {
  tentativas <- tentativas + 1
  cand <- spatSample(area_sorteio, size = 1, method = "random")
  
  if (nrow(cand) > 0) {
    if (is.null(pontos_finais)) {
      pontos_finais <- cand
    } else {
      dists <- distance(cand, pontos_finais)
      if (min(dists) >= 2000) { 
        pontos_finais <- rbind(pontos_finais, cand)
      }
    }
  }
  
  if (tentativas > 5000) break
}

# Aplicando buffers
# Verificar CRS
if (!identical(crs(pontos_finais), crs(mun_2024_utm))) {
  pontos_finais <- project(pontos_finais, crs(mun_2024_utm))
}

buffers <- buffer(pontos_finais, width = 1500)
buffers$id_paisagem <- 1:nrow(buffers)

# Extração e processamento
extracao <- terra::extract(mun_2024_utm, buffers)
nome_col_raster <- names(extracao)[2]

tabela_final <- extracao %>%
  rename(id_buffer = ID, categoria_bruta = !!sym(nome_col_raster)) %>%
  group_by(id_buffer) %>%
  mutate(total_px = n()) %>%
  group_by(id_buffer, categoria_bruta) %>%
  summarise(
    pixels = n(),
    percentagem = (pixels / first(total_px)) * 100,
    .groups = "drop"
  ) %>%
  mutate(Categoria_Limpa = as.character(categoria_bruta)) %>%
  left_join(legenda %>% select(label, Cor), by = c("Categoria_Limpa" = "label"))

write.csv(tabela_final, "uso_solo_murici_final.csv", row.names = FALSE)
writeVector(pontos_finais, "pontos_final.shp", overwrite=TRUE)

# Vizualizando resultado 
ggplot() +
  geom_spatraster(data = mun_2024_utm) +
  scale_fill_manual(values = legenda$Cor, na.value = "white", name = "Uso do Solo") +
  geom_spatvector(data = buffers, fill = NA, color = "black", linewidth = 1) +
  geom_spatvector(data = pontos_finais, color = "red", size = 2, shape = 19) +
  geom_spatvector_text(data = pontos_finais, aes(label = 1:15), 
                       color = "white", size = 3, vjust = -0.8) +
  labs(title = "Pontos amostrais e buffers de 1500m") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Painel 5X3
lista_recortes <- list()
for(i in 1:nrow(buffers)) {
  crop_i <- crop(mun_2024_utm, buffers[i, ])
  mask_i <- mask(crop_i, buffers[i, ])
  df_i <- as.data.frame(mask_i, xy = TRUE, cells = TRUE)
  df_i$id_buffer <- paste("Paisagem", i)
  lista_recortes[[i]] <- df_i
}
df_painel <- bind_rows(lista_recortes)

ggplot(df_painel) +
  geom_tile(aes(x = x, y = y, fill = label)) +
  scale_fill_manual(values = setNames(legenda$Cor, legenda$label)) +
  facet_wrap(~id_buffer, nrow = 3, ncol = 5, scales = "free") +
  theme_minimal() +
  labs(title = "Painel de Amostragem: 15 Paisagens de Murici-AL",
       subtitle = "Recortes circulares de 1500m de raio (Interior de Floresta)",
       fill = "Uso do Solo") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

agregacao <- lsm_c_ai(mun_2024_utm)
print(agregacao)

agregacao_tab <- agregacao %>%
  mutate(
    classe_nome = legendatab[as.character(class)]
  ) %>%
  select(
    classe_nome,
    value
  ) %>%
  rename(
    Classe = classe_nome,
    AI = value
  ) %>%
  arrange(desc(AI))

agregacao_tab
