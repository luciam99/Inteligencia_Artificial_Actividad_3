############################################################
### Análisis de un conjunto de datos de origen biológico mediante técnicas de machine learning supervisadas y no supervisadas
### 1) Reducción dimensionalidad: PCA, t-SNE, Isomap,LLE,LE,MVU,UMAP
### 2) Clusterización: K-means, Dendogramas, Heatmap
### 3) Supervisado: SVM, KNN, SVM-Kernel,Árboles de decisión,Random Forest + métricas (CM, Precision, Sensitivity, Specificity, F1)
############################################################

library(tidyverse)
library(caret)
library(factoextra)
library(randomForest)
library(ggplot2)
library(stats)  
library(Rtsne)
library(kernlab)

# ###############################################################################
# Generación de un único dataframe
# ###############################################################################

f_c <- list.files(pattern="classes.csv$", ignore.case=TRUE)
df_c<-read.csv(f_c, header = FALSE,sep = ";",col.names = c("ID","Clase"))


colnames_g <- readLines("column_names.txt")
f_g <- list.files(pattern="gene_expression.csv$", ignore.case=TRUE)

df <- df_c %>% # Se unifica el dataframe que contiene las clases con el que contiene la expresión génica
  mutate(read.csv(
    f_g,
    header = FALSE,
    # Se asigna el título (nombre de gen) de cada columna del dataframe que contiene los datos de expresión génica
    col.names = colnames_g,
    sep = ";"
  ))

# ###############################################################################
# Imputación de NAs
# ###############################################################################

df_NAs <- df %>%
  select(where( ~ any(is.na(.)))) # La imputación no procede puesto que el dataset no contiene valores faltantes

# ###############################################################################
# Exploración de los datos
# ###############################################################################

table (df$Clase) # Evaluar el balance las clases

genes <- df %>% select(DUOXA1:TTC31) # dataframe con los datos de expresión génica de cada gen
range(genes, na.rm = TRUE) # Consultar el intervalo de valores de expresión génica
res_gen<-t(summary(genes)) # Resumen de estadistícos descriptivos para cada gen


# Gráfico de Densidad para evaluar la distribución global de los valores de expresión génica

g_lar <- genes %>% #Convertir a formato largo para ggplot
  pivot_longer(everything(), names_to = "Gen", values_to = "Expresion")

ggplot(g_lar, aes(x = Expresion)) +
  geom_density(fill = "orange", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribución de la Expresión Génica", 
       x = "Nivel de Expresión", y = "Densidad") 

# Histograma de los ceros

porcentaje_ceros <- colMeans(genes == 0) * 100 # Cálculolo del % de ceros para cada gen

# Histograma para identificar genes con altos porcentaje de ceros
hist(porcentaje_ceros, 
     main = "Histograma de Sparsity (Ceros)",
     xlab = "% de ceros por gen", 
     col = "lightgreen")


# Se genera un nuevo dataframe donde se excluyen los genes con % de ceros superior a 80 (por considerarse menos informativos)

df_filtrado <- bind_cols(df %>% select(ID, Clase), genes[, porcentaje_ceros <=80])

# Escalado de datos 

X_scaled <- df_filtrado %>%
  select(where(is.numeric)) %>% scale() #Aplicación de Z-score (Media=0, SD=1) a la expresión génica

# ###############################################################################
# Métodos de aprendizaje no supervisado
# ###############################################################################

############################Reducción dimensionalidad###########################

# PCA--------------------------------------------------------------------------

# Se aplica PCA a la matriz de datos previamente escalada
pca <- prcomp(X_scaled, center=FALSE)

# Se extraen coordenadas de las muestras en los dos primeros componentes
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Clase = df_filtrado$Clase #Agrupar según clase
)

# Visualización de la estructura global de los datos con PCA
ggplot(pca_df, aes(PC1, PC2, color = Clase)) +
  geom_point(size = 2, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA (PC1 vs PC2)", x = "PC1", y = "PC2")

# t-SNE-------------------------------------------------------------------------

# t-SNE no permite duplicados 
X_unique <- unique(X_scaled)

# Establecer semilla para reproducibilidad 
set.seed(42)

# Ejecución de t-SNE
tsne_res <- Rtsne(X_scaled, 
                  dims = 2, 
                  perplexity = 30, 
                  verbose = FALSE, 
                  max_iter = 500,
                  check_duplicates = FALSE) 

# Dataframe para visualización
tsne_df <- data.frame(
  TSNE1 = tsne_res$Y[, 1],
  TSNE2 = tsne_res$Y[, 2],
  Clase = df_filtrado$Clase
)

# Visualización
ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = Clase)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "t-SNE: Visualización de clústeres biológicos",
       subtitle = "Basado en expresión génica filtrada y escalada",
       x = "t-SNE dimensión 1", 
       y = "t-SNE dimensión 2")

# ################################Clusterización################################

set.seed(123)

# K-means ----------------------------------------------------------------------

# Elbow method para evaluar el número de centroides óptimo
fviz_nbclust(X_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 5, linetype = 2) +
  theme_classic()

# Clustering con K-means
kmeans.result <- kmeans(X_scaled, centers = 5, nstart = 25)

# Visualización
fviz_cluster(kmeans.result, data = X_scaled, geom = "point") +
  ggtitle("K-means (k = 5)", subtitle = "") +
  theme_minimal()

# Comparación de los resultados con los datos reales
table(Cluster_Kmeans = kmeans.result$cluster, Clase_Real = df_filtrado$Clase)


# Clustering jerárquico aglomerativo -------------------------------------------

# Cálculo de la matriz de distancias
dist_matrix <- dist(X_scaled)

# Clustering jerárquico
hclust_model <- hclust(dist_matrix, method = "ward.D2")

# Visualización del dendograma
fviz_dend(hclust_model,
          k = 5,
          cex = 0.5,
          rect = TRUE,
          main = "Dendrograma de Clustering Jerárquico (Ward)",
          xlab = "Índice de Observaciones",
          ylab = "Distancia") + theme_classic()

# Visualización del heatmap
heatmap(as.matrix(dist_matrix),
        symm = TRUE, 
        distfun = function(x) as.dist(x),
        hclustfun = function(x) hclust(x, method = "ward.D2"),
        labRow = FALSE, labCol = FALSE)

# Comparación de los resultados con los datos reales
df_filtrado$Cluster_hclust <- cutree(hclust_model, k = 5)
table(Cluster_hclust = df_filtrado$Cluster_hclust, Clase_Real = df_filtrado$Clase)


# ##############################################################################
# Métodos de aprendizaje supervisado
# ##############################################################################

# SVM (Support Vector Machine)---------------------------------------------------

# Preparación de los datos
df_filtrado$Clase <- as.factor(df_filtrado$Clase)

# Partición de datos: 80% entrenamiento, 20% prueba
trainIndex <- createDataPartition(df_filtrado$Clase, p = 0.8, list = FALSE)
train_set <- df_filtrado[trainIndex, ]
test_set  <- df_filtrado[-trainIndex, ]

# Configuración: validación cruzada 5-fold
ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)

# Entrenamiento SVM con kernel lineal (adecuado para alta dimensionalidad)
svm_model <- train(
  Clase ~ ., 
  data = train_set %>% select(-ID), # Excluir ID (no es variable predictora)
  method = "svmLinear",
  trControl = ctrl,
  preProcess = c("center", "scale") # Normalización Z-score
)

# Predicciones SVM en conjunto de prueba
predictions <- predict(svm_model, newdata = test_set)

# Random Forest ----------------------------------------------------------------
# Entrenamiento Random Forest: robusto para datos de alta dimensionalidad
set.seed(123) # Para reproducibilidad
rf_model <- train(
  Clase ~ ., 
  data = train_set[, colnames(train_set) != "ID"], # Excluir ID
  method = "rf",
  trControl = ctrl, # Misma configuración que SVM
  ntree = 150 # Número de árboles (balance tiempo/precisión)
)

# Predicciones Random Forest
predictions_rf <- predict(rf_model, test_set)

# KNN --------------------------------------------------------------------------
# Entrenamiento K-Nearest Neighbors: clasificación por similitud
knn_model <- train(
  Clase ~ .,
  data = train_set[, colnames(train_set) != "ID"], # Excluir ID
  method = "knn",
  trControl = ctrl, # Misma configuración
  tuneGrid = expand.grid(k = c(3, 5, 7)) # Valores de k probados
)

# Predicciones KNN
predictions_knn <- predict(knn_model, test_set)

# Cálculo de métricas ----------------------------------------------------------
# Función para calcular métricas multiclase (5 clases)
get_metrics <- function(pred, actual) {
  cm <- confusionMatrix(pred, actual)
  # Promedio macro (todas las clases igual importancia)
  c(
    Accuracy = cm$overall["Accuracy"],
    Precision = mean(cm$byClass[, "Precision"], na.rm = TRUE),
    Sensitivity = mean(cm$byClass[, "Sensitivity"], na.rm = TRUE),
    Specificity = mean(cm$byClass[, "Specificity"], na.rm = TRUE),
    F1 = mean(cm$byClass[, "F1"], na.rm = TRUE)
  )
}

# Tabla comparativa de los 3 modelos
results <- data.frame(
  Modelo = c("SVM", "Random Forest", "KNN"),
  rbind(
    get_metrics(predictions, test_set$Clase),
    get_metrics(predictions_rf, test_set$Clase),
    get_metrics(predictions_knn, test_set$Clase)
  )
)

# Resultados -------------------------------------------------------------------
print("Resultados comparación modelos supervisados")
print(results)


################################################################################


