# Analise_Proteomica
Script base em R para análise de expressão diferencial em proteômica (LFQ) utilizando pacotes do Bioconductor (ex. DEP, clusterProfiler).Utilizado para análise de dados de projeto de Mestrado.

## Contexto e Uso do Script
Este código funciona como uma base para análise estatística e visualização dos dados obtidos por espectrometria de massas. 

**Nota sobre reprodutibilidade:** O script aqui presente representa a estrutura metodológica base. Durante o desenvolvimento do projeto, filtragens pontuais, alterações de *cutoffs* (ex: Log2FC e p-valor), testes de diferentes contrastes e outras modificações visuais foram executados a partir deste script, quando necessário.

## Principais Etapas do Pipeline
1. **Importação dos dados:** Importação, filtragem de *missing values* e definição do desenho experimental (Condições, Linhagens, Tratamentos e Replicatas).
2. **Normalização e Imputação:** Estabilização de variância (VSN) e imputação de dados ausentes (método *MinProb* para perdas não-aleatórias).
3. **Análise de Expressão Diferencial:** Identificação de proteínas significativamente reguladas (UP/DOWN) utilizando testes estatísticos padronizados pelo pacote DEP.
4. **Visualização de Dados:** Geração de *PCA*, *Volcano Plots*, *Heatmaps* (via ComplexHeatmap) e Gráficos de Dispersão personalizados.
5. **Análise de Enriquecimento Funcional:** Extração de clusters hierárquicos e análise de Gene Ontology (Processos Biológicos) das proteínas diferenciais utilizando clusterProfiler e visualização via enrichplot.

## Principais Dependências
Certifique-se de ter as seguintes bibliotecas instaladas no R antes de executar:
* DEP
* clusterProfiler
* ComplexHeatmap
* ggplot2 e dplyr

## Contato 
Desenvolvido por adsc16 (https://github.com/adsc16) para o análise de dados do projeto de Mestrado em USP-FMRP.
