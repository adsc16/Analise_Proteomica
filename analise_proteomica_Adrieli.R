# 00. Instalação e carregamento de pacotes
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DEP", force = TRUE)
BiocManager::install(c("MSnbase", "MsCoreUtils", "DEP"), ask = FALSE, update = TRUE)
library(MSnbase)
library(MsCoreUtils)
library("DEP")
library("dplyr")


# Pacotes usados
library(vsn) 
library(MSnbase) 
library(SummarizedExperiment) 
library(dplyr) 
library(stringr) 
library(ggplot2) 
library(RColorBrewer) 
library(ggrepel) 
library(tidyverse) 
library(clusterProfiler)
library(org.Mm.eg.db) #para mus musculus (org.Hs.eg.db para humano)

# ------------------------------------------------------------------------------
# 01. Importação e preparo dos dados
# ------------------------------------------------------------------------------
# caminho do arquivo
caminho_arquivo <- "C:/Users/User/Desktop/Mestrado - FMRP/Projeto/Resultados/MS/BaF3 Tratamentos/Tratamentos/tratamentos.csv"

# arquivo CSV
input <- read.csv(caminho_arquivo, header = TRUE, sep = ";", dec = ",", stringsAsFactors = FALSE)

input$Gene.names %>% duplicated() %>% any()

# Garantir que todas as colunas de intensidade sejam numéricas
intensity_cols <- 2:ncol(input)
for(col in intensity_cols) {
  input[, col] <- as.numeric(input[, col])
}

# A função make_unique() cria nomes únicos para cada proteína.
# Gera as colunas 'name' e 'ID' 
data_unique <- make_unique(input, "Protein", "Protein", delim = ";")

# ------------------------------------------------------------------------------
# 02. Desenho experimental
# ------------------------------------------------------------------------------
# O pacote DEP exige um data.frame com colunas obrigatórias: label, condition, replicate
sample_names <- colnames(data_unique)[intensity_cols]

# Extraindo informações dos nomes das amostras 
tipo_celular <- sapply(strsplit(sample_names, "_"), "[[", 1) # WT ou T315I
tratamento   <- sapply(strsplit(sample_names, "_"), "[[", 2) # CTR, CZ, IMT, etc.
replicata    <- as.numeric(sapply(strsplit(sample_names, "_"), "[[", 3)) # 01, 02, etc.

experimental_design <- data.frame(
  label = sample_names,
  condition = paste0(tipo_celular, "_", tratamento), # Ex: WT_CTR
  replicate = replicata,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 03. Objeto Summarized experiment
# ------------------------------------------------------------------------------
data_se <- make_se(
  data_unique,
  columns = intensity_cols,
  expdesign = experimental_design
)

print(data_se)
# ------------------------------------------------------------------------------
# 04. Filtragem 
# ------------------------------------------------------------------------------
my_colors <- c(
  brewer.pal(8, "Set1"),
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(9, "Set3"),
  brewer.pal(8, "Accent")
)[1:42]

# Frequência de quantificação das proteínas por amostra
p <- plot_frequency(data_se)
p + scale_fill_manual(values = my_colors)

# Filtragem de valores faltantes:
# thr = 0 significa que a proteína deve ser detectada em TODAS as replicatas de pelo menos UMA condição.
data_filt <- filter_missval(data_se, thr = 0)

# Visualizar número de proteínas detectadas após filtro
plot_numbers(data_filt) +
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sobreposição de proteínas entre condições
plot_coverage(data_filt) +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

# ------------------------------------------------------------------------------
# 05. Normalização e imputação
# ------------------------------------------------------------------------------
# Normalização VSN (Variance Stabilizing Normalization)
data_norm <- normalize_vsn(data_filt)

# Comparar antes e depois da normalização
plot_normalization(data_filt, data_norm)

# Visualizar distribuição de valores faltantes (Missing values)
plot_missval(data_filt)
plot_detect(data_filt)

# Imputação de dados faltantes (MinProb é ideal para valores faltantes não aleatórios - MNAR)
# Simula valores baixos baseados em uma distribuição gaussiana
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Visualizar efeito da imputação
plot_imputation(data_norm, data_imp)

# ------------------------------------------------------------------------------
# 06. Análise de expressão diferencial
# ------------------------------------------------------------------------------
# Extrair GENOTYPE (WT ou T315I)
genotype <- str_extract(colnames(data_imp), "WT|T315I")
# Extrair CONDITION (CTR, CZ, IMT, NT, AMS, CH)
condition <- str_extract(colnames(data_imp), "CTR|CZ|IMT|NT|AMS|CH")

# Extrair REPLICATE 
replicate <- str_extract(colnames(data_imp), "[0-9]+$")

# Criar condição combinada
combined_condition <- paste(genotype, condition, sep = "_")

# Verificar se combinou corretamente
combined_condition

# Atualizando colData com as variáveis separadas para testes, se necessário
colData(data_imp)$genotype <- genotype
colData(data_imp)$condition <- combined_condition
colData(data_imp)$replicate <- replicate

# Teste de expressão diferencial (todos os contrastes possíveis)
data_diff <- test_diff(data_imp, type = "all")

# Adicionando rejeições
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

# Obter tabela completa de resultados
df_results <- get_results(dep)

# ------------------------------------------------------------------------------
# 07. Visualizações (PCA, Heatmap, Volcano)
# ------------------------------------------------------------------------------
cores_12 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", 
              "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#000000", "#FFD700")

# PCA (exemplo com 500 proteínas)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4) +
  scale_color_manual(values = cores_12) +
  theme(plot.title = element_blank())

# Correlação de Pearson
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "RdYlBu")

# Heatmap de proteínas significativas
plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 5, 
             show_row_names = FALSE, indicate = c("condition", "replicate"))

# Volcano Plot de um contraste específico (Exemplo: T315I CTR vs WT CTR)
plot_volcano(dep, contrast = "T315I_CTR_vs_WT_CTR", label_size = 3, add_names = TRUE)

# Plot de uma proteína específica, exemplo sf3b1
plot_single(dep, proteins = "Sf3b1", type = "centered")

# ------------------------------------------------------------------------------
# 08. Gráficos (Dispersão)
# ------------------------------------------------------------------------------
# 8.1 Dispersão - exemplo: Comparação do efeito do NT157 entre WT e T315I
ggplot(df_results, aes(x = WT_CTR_vs_WT_NT_ratio, y = T315I_CTR_vs_T315I_NT_ratio)) +
  geom_point(alpha = 0.6, color = "gray50") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_point(data = subset(df_results, abs(WT_CTR_vs_WT_NT_ratio) > 2 & abs(T315I_CTR_vs_T315I_NT_ratio) < 0.5),
             aes(x = WT_CTR_vs_WT_NT_ratio, y = T315I_CTR_vs_T315I_NT_ratio), color = "red", size = 2) +
  geom_text_repel(data = subset(df_results, abs(WT_CTR_vs_WT_NT_ratio) > 2 & abs(T315I_CTR_vs_T315I_NT_ratio) < 0.5),
                  aes(label = name), size = 3, max.overlaps = 15) +
  labs(title = "Comparação do Efeito do NT157: WT vs T315I",
       x = "Log2 FC na WT (Controle / NT157)",
       y = "Log2 FC na T315I (Controle / NT157)") +
  theme_bw()

# ------------------------------------------------------------------------------
# 9. Cálculo do z-score (Exportação para excel)
# ------------------------------------------------------------------------------
# Utilizando a matriz original (removendo a primeira coluna de nomes)
nomes_genes <- input[, 1] 
matriz_bruta <- as.matrix(input[, -1])
class(matriz_bruta) <- "numeric"

# Aplicar Log2 (+1 para evitar log(0))
matriz_log <- log2(matriz_bruta + 1)

# Calcular Z-Score por linha (Proteína). scale() age nas colunas, por isso transpomos (t) duas vezes.
matriz_zscore <- t(scale(t(matriz_log)))

# Salvar o resultado
df_zscore <- data.frame(Genes = nomes_genes, matriz_zscore, check.names = FALSE)
write.csv2(df_zscore, "Todas_Proteinas_ZScore.csv", row.names = FALSE)

# ==========================================================================
# 10.Análise de Enriquecimento Funcional (Processos Biológicos)
# Comparar clusters hierárquicos e gerar visualizações de vias
# ==========================================================================

# --- 0: Carregamento das Bibliotecas ---
library(clusterProfiler)
library(org.Hs.eg.db) #humano
library(enrichplot)
library(ggplot2)
library(dplyr)

library(ComplexHeatmap)

# 1. Heatmap 

meu_heatmap <- plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, 
                            col_limit = 5, show_row_names = FALSE, 
                            indicate = c("condition", "replicate"))

# 2. Extrair a ordem das linhas (proteínas) separadas pelas 6 fatias (clusters) do gráfico
lista_clusters <- row_order(meu_heatmap)

# 4. Pegar os nomes das proteínas significativas 
# (Eles estão na mesmíssima ordem que o DEP envia para o gráfico)
nomes_sig <- df_results$name[df_results$significant]

# 5. Montar a 'tabela_final' amarrando o nome da proteína ao cluster do gráfico
tabela_final <- data.frame(name = character(), Cluster = character(), stringsAsFactors = FALSE)

for (i in seq_along(lista_clusters)) {
  # Extrai os índices matemáticos daquele pedaço do heatmap
  indices_do_cluster <- lista_clusters[[i]]
  
  # Puxa os nomes reais das proteínas usando esses índices
  proteinas_do_cluster <- nomes_sig[indices_do_cluster]
  
  # Cola tudo na nossa tabela final
  tabela_final <- rbind(tabela_final, data.frame(
    name = proteinas_do_cluster,
    Cluster = paste0("Cluster_", i),
    stringsAsFactors = FALSE
  ))
}

# 6. Verificar quantas proteínas caíram em cada cluster do Heatmap
print("Quantidade de proteínas por cluster extraídas do Heatmap:")
table(tabela_final$Cluster)

# --- PASSO 6: Preparação dos Dados ---
# 6.1 Converter Símbolos de Genes (Symbols) para ENTREZ ID (usando o banco de camundongo)
ids <- bitr(tabela_final$name, 
            fromType = "SYMBOL", 
            toType = "ENTREZID", 
            OrgDb = "org.Mm.eg.db")

# 6.2 Unir os IDs numéricos de volta na tabela original de clusters
dados_com_ids <- merge(tabela_final, ids, by.x = "name", by.y = "SYMBOL")

# --- PASSO 7: Rodar a Comparação de Clusters ---
comparacao_go <- compareCluster(
  ENTREZID ~ Cluster, 
  data = dados_com_ids, 
  fun = "enrichGO",       
  OrgDb = "org.Mm.eg.db", 
  ont = "BP",             
  pAdjustMethod = "BH",   
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE         
)

# --- PASSO 8: Visualizações ---
plot_geral <- dotplot(comparacao_go, showCategory = 5) + 
  ggtitle("Comparação Funcional Geral dos Clusters (Processos Biológicos)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot_geral)

# ==========================================================================
# 11. Gráfico de Barras de Proteínas Diferencialmente Expressas (DEP)
# Contar e plotar o número de proteínas UP e DOWN por tratamento
# ==========================================================================

# 1. Extrair resultados do objeto DEP
df_results <- get_results(dep)

# 2. Configurar nomes e códigos dos tratamentos
treatments <- c("IMT", "NT", "AMS", "CZ", "CH")
names(treatments) <- c("Imatinibe", "NT157", "Amsacrina", "Clotrimazol", "Clorexidina")

# 3. Função de extração e contagem
# Lógica: Ratio < 0 = Proteína aumentada no Tratamento (UP)
extract_counts <- function(trt_code, trt_name, lineage) {
  
  # Monta o nome das colunas geradas pelo DEP (ex: WT_CTR_vs_WT_IMT)
  col_prefix <- paste0(lineage, "_CTR_vs_", lineage, "_", trt_code)
  col_sig    <- paste0(col_prefix, "_significant")
  col_ratio  <- paste0(col_prefix, "_ratio")
  
  # Se a coluna não existir, retorna nulo para evitar erros
  if(!col_sig %in% colnames(df_results)) return(NULL)
  
  # Contagem dos significativos
  # Ratio < 0 --> Maior no Tratamento (UP/Vermelho)
  # Ratio > 0 --> Menor no Tratamento (DOWN/Azul)
  n_up   <- sum(df_results[[col_sig]] & df_results[[col_ratio]] < 0, na.rm = TRUE)
  n_down <- sum(df_results[[col_sig]] & df_results[[col_ratio]] > 0, na.rm = TRUE)
  
  # Retorna um dataframe com o resumo
  return(data.frame(
    Treatment  = trt_name,
    Lineage    = lineage,
    Type       = c("Regulada positivamente", "Regulada negativamente"),
    Count      = c(n_up, n_down),
    Plot_Value = c(n_up, -n_down) # O valor negativo joga a barra do DOWN para baixo do eixo X
  ))
}

# 4. Iterar sobre os tratamentos e gerar os dados
plot_data_list <- list()

for (trt in names(treatments)) {
  plot_data_list[[paste0("WT_", trt)]]    <- extract_counts(treatments[trt], trt, "WT")
  plot_data_list[[paste0("T315I_", trt)]] <- extract_counts(treatments[trt], trt, "T315I")
}

# Combina a lista em um único dataframe
plot_data <- bind_rows(plot_data_list)

# --- AQUI ESTÁ O TRUQUE DA ORDEM ---
# 4.1 Definimos a ordem exata desejada no eixo X
ordem_desejada <- c("Imatinibe", "NT157", "Amsacrina", "Clotrimazol", "Clorexidina")

# 4.2 Transformamos a coluna Treatment em FATOR com os LEVELS fixos
plot_data$Treatment <- factor(plot_data$Treatment, levels = ordem_desejada)

# 4.3 Ajuste estético dos nomes das linhagens (WT e T315I)
plot_data$Lineage <- factor(plot_data$Lineage, 
                            levels = c("WT", "T315I"),
                            labels = c("Ba/F3 WT", "Ba/F3 T315I"))

# 5. Plotar o Gráfico de Barras Espelhado
grafico_barras <- ggplot(plot_data, aes(x = Treatment, y = Plot_Value, fill = Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  
  # Cores personalizadas
  scale_fill_manual(values = c("Regulada negativamente" = "#377EB8", 
                               "Regulada positivamente" = "#E41A1C")) +
  
  # Adiciona os números absolutos em cima/embaixo das barras
  geom_text(aes(label = abs(Count), 
                vjust = ifelse(Plot_Value > 0, -0.5, 1.5)), 
            size = 3.5, fontface = "bold") +
  
  # Linha preta marcando o zero
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # Separa o gráfico por linhagem
  facet_grid(~Lineage, scales = "free_x", space = "free_x", switch = "x") +
  
  # Textos e Títulos
  labs(
    title    = "Proteínas Diferencialmente Expressas (DEP)",
    subtitle = "Comparação: Controle vs Tratado (Log2FC cutoff: 2 | FC = 4)",
    y        = "Número de Proteínas", 
    x        = NULL, 
    fill     = "Efeito do tratamento:"
  ) +
  
  # Estética e Limpeza do tema
  theme_minimal() +
  theme(
    legend.position    = "top",
    strip.placement    = "outside",
    strip.text         = element_text(size = 12, face = "bold"),
    strip.background   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    panel.grid.major.x = element_blank()
  )

# Mostra o gráfico final
print(grafico_barras)
