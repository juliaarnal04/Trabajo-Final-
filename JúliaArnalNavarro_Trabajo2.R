# JúliaArnalNavarro_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")
BiocManager::install("RCurl", lib = Sys.getenv("R_LIBS_USER")); library(RCurl)
install.packages("reshape2")  
library(reshape2)
library(ggplot2)

datos <- read.csv("C:/Users/julia/Downloads/datos_biomed (1).csv", header = TRUE)


# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)
ls()
head(datos)
summary(datos)
dim(datos)
str(datos)

# Hay 6 variables, ya que cuando hemos ejecutado dim nos ha mostrado: [1] 100 6. Donde 100 significa el numero de filas y 6 las variables
	# Estas variables son: ID, tratamiento, glucosa, presion, colesterol y Tratamiento_factor
	
# En sumary (datos) dentro del tratamiento_factor encontramos: 
	# Farmaco A (36)
	# Farmaco B (31)
	# Placeno (33)
# Por lo tanto, hay 3 tratamientos. 

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)

datos_box <- datos[, c("Tratamiento", "Glucosa", "Presion", "Colesterol")]

datos_largo <- melt(datos_box, id.vars = "Tratamiento",
                    variable.name = "Variable",
                    value.name = "Valor")
boxplot(Valor ~ interaction(Variable, Tratamiento),
        data = datos_largo,
        col = c("lightblue", "lightgreen", "lightpink"),
        las = 2,
        main = "Boxplots por Tratamiento y Variable",
        xlab = "Variable y Tratamiento",
        ylab = "Valor")


# 4. Realiza un violin plot (investiga qué es). (1 pt)

# Un violin plot es una gráfica que combina un boxplot y un kernel density plot.
# Muestra la distribución de los datos y otro tipo de información como la mediana, los cuartiles y posibles valores atípicos.

library(ggplot2)

datos <- read.csv("C:/Users/julia/Downloads/datos_biomed (1).csv", header = TRUE)

datos_largo <- reshape(
  datos,
  varying = c("Glucosa", "Presion", "Colesterol"),
  v.names = "Valor",
  timevar = "Variable",
  times = c("Glucosa", "Presion", "Colesterol"),
  idvar = "ID",
  direction = "long"
)

ggplot(datos_largo, aes(x = Tratamiento, y = Valor, fill = Tratamiento)) +
  geom_violin(trim = FALSE, alpha = 0.6) +     
  geom_boxplot(width = 0.1, fill = "white") +  
  facet_wrap(~Variable, scales = "free_y") +   
  theme_minimal() +
  labs(title = "Violin plots por tratamiento y variable",
       x = "Tratamiento",
       y = "Valor") +
  scale_fill_manual(values = c("FarmacoA" = "lightblue",
                               "FarmacoB" = "lightgreen",
                               "Placebo" = "lightpink"))



# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)

plot(datos$Glucosa, datos$Presion,
     col = c("blue", "green", "pink")[as.factor(datos$Tratamiento)], 
     pch = 19,  
     xlab = "Glucosa",
     ylab = "Presion",
     main = "Gráfico de dispersión: Glucosa vs Presion")

legend("bottomright", 
       legend = c("FarmacoA", "FarmacoB", "Placebo"), 
       col = c("blue", "green", "pink"), 
       pch = 19)


# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)

ggplot(datos, aes(x = Presion, y = Colesterol, color = Tratamiento)) +
  geom_point(size = 3) + 
  facet_grid(. ~ Tratamiento) + 
  theme_minimal() +
  labs(title = "Colesterol vs Presion por Tratamiento",
       x = "Presión",
       y = "Colesterol") +
  scale_color_manual(values = c("FarmacoA" = "blue",
                                "FarmacoB" = "green",
                                "Placebo" = "pink"))

# Un facet grid es un tipo de gráfico (ggplot2)
# Facet grid divide un gráfico en múltiples subgráficos (facetas) según una o más variables. 
# En este caso podemos observar como subgraficos los diferentes tratamamientos comparados con las variables de presion y colesterol. 

# 7. Realiza un histogramas para cada variable. (0.5 pts)

hist(datos$Glucosa, 
     main = "Histograma de Glucosa", 
     xlab = "Glucosa", 
     col = "lightblue", 
     border = "black")

hist(datos$Presion, 
     main = "Histograma de Presión", 
     xlab = "Presión", 
     col = "lightgreen", 
     border = "black")

hist(datos$Colesterol, 
     main = "Histograma de Colesterol", 
     xlab = "Colesterol", 
     col = "lightpink", 
     border = "black")


# 8. Crea un factor a partir del tratamiento. Investifa factor(). (1 pt)

datos$Tratamiento_factor <- factor(datos$Tratamiento)
levels(datos$Tratamiento_factor)

# Factor () se usa para representar variables categóricas (en este caso, trataminto) haciendo así categorías discretas (no números continuos)

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)

media_Glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = mean)
media_Glucosa

sd_Glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = sd)
sd_Glucosa

# La medica de glucosa es: 110.47, 105.58, 103.04
# La desviacion estandar de glucosa es: 13.58, 12.15, 17.18 
# Ambos datos son respectivamente del farmacoA, farmacoB y placebo 

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)

Placebo <- subset(datos, Tratamiento == "Placebo")
FarmacoA <- subset(datos, Tratamiento == "FarmacoA")
FarmacoB <- subset(datos, Tratamiento == "FarmacoB")

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)

shapiro.test(Placebo$Glucosa)
shapiro.test(FarmacoA$Glucosa)
shapiro.test(FarmacoB$Glucosa)

shapiro.test(Placebo$Presion)
shapiro.test(FarmacoA$Presion)
shapiro.test(FarmacoB$Presion)

shapiro.test(Placebo$Colesterol)
shapiro.test(FarmacoA$Colesterol)
shapiro.test(FarmacoB$Colesterol)

# Podemos observar que todos los datos siguen una distribución normal ya que p-value > 0.05
# Como siguen la distribución normal, hacemos ahora la comparativa; usaremos ANOVA

anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos)
summary(anova_glucosa)

anova_presion <- aov(Presion ~ Tratamiento, data = datos)
summary(anova_presion)

anova_colesterol <- aov(Colesterol ~ Tratamiento, data = datos)
summary(anova_colesterol)

# En la glucosa: Pr(>f)= 0.1. 
	# Como es un valor mayor que 0.05, indica que no hay diferencias muy grandes entre los diversos tratamientos
	# Hay una cantidad similar de glucosa en Placebo, FarmacoA y FarmacoB

# En la presion (Pr(>F) = 1.2e-05) y en el colesterol (Pr(>F) = 0.00429) -> tenemos valores menores que 0.05
	# Esto nos idica que hay diferencias significativas entre los tratamientos.
	# Por lo tanto, al menos un tratamiento tiene una media de presion y de colester0ol bastante diferente de los otros (no tiene porque ser el mismo)

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)

anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos)
summary(anova_glucosa)





