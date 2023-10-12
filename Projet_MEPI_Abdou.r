
library(sensitivity)




### Mod?le dont la sensibilit? doit ?tre analys?e dans le cadre du projet MODE-MPI 2023-2024

### Le mod?le est ici d?fini sous forme de fonction pour faciliter vos analyses de sensibilit? (AS)
### La fonction renvoie les sorties ponctuelles qui sont ? analyser dans l'AS
parametre=ValNominale

modAppli <- function(parametre){  

  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (en jours)
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de param?tres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4)

  # boucle des sc?narios de l'?chantillonnage de l'AS
  for (i in 1:nrow(parametre)) { 

    # STRUCTURE & PARAMETRES DU MODELE

    # Récuperation des para
    K = parametre[i,1];    # Capacité de charge de l'environnement
    sr = parametre[i,2];   # Taux de survie
    m1 = parametre[i,3];   # Taux de mortalité classe 1
    m2 = parametre[i,4];   # Taux de mortalité classe 1
    m3 = parametre[i,5];   # Taux de mortalité classe 3
    f2 = parametre[i,6];   # Taux de fécondité classe 2
    f3 = parametre[i,7];   # Taux de fécondité classe 3
    portee = parametre[i,8]; # Portée (nombre d'individus produits par génération)
    t1 = parametre[i,9];   # Taux de transmission classe 1
    t2 = parametre[i,10];  # Taux de transmission classe 2
    
    # Transitions entre les états de santé
    trans = parametre[i,11];  # Taux de transition de XX
    lat = parametre[i,12];    # Taux de latence
    rec = parametre[i,13];    # Taux de récupération
    loss = parametre[i,14];   # Taux de perte
    madd = parametre[i,15];   # Taux d'addition de nouveaux individus XX

    # INITIALISATION
    MAT <- array(0, dim=c(4,4,temps)); # nb indiv par classe d'?ge en ligne (derni?re ligne = pop tot), ?tat de sant? en colonne, pas de temps (dimension 3)
    nvinf <- array(0, dim=c(temps));
    # conditions initiales (la population est ? sa structure d'?quilibre, calcul?e par ailleurs)
    MAT[1,1,1] <- 27; # xx
    MAT[2,1,1] <- 23; # xx
    MAT[3,1,1] <- 36; # xx
    MAT[3,3,1] <- 1;  # xx
    # effectifs par ?tat de sant?
    MAT[4,1,1] <- sum(MAT[1:3,1,1]); 
    MAT[4,2,1] <- sum(MAT[1:3,2,1]); 
    MAT[4,3,1] <- sum(MAT[1:3,3,1]); 
    MAT[4,4,1] <- sum(MAT[1:3,4,1]);

    # SIMULATIONS
    # boucle du temps
    for (t in 1:(temps-1)) { 
     # classe d'?ge xx
      # RQ : les naissances sont XX, les nouveaux n?s ?tant dans l'?tat XX
      N <- sum(MAT[4,,t]);	# taille de la pop en t
	MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t]      + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
	MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			  + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
	MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec)  		  + lat*MAT[1,2,t]; 
	MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) 		  + rec*MAT[1,3,t]; 
     # classe d'?ge xx
	MAT[2,1,t+1] <- MAT[1,1,t]*t1	+ MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
	MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
	MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-rec)		+ lat*MAT[2,2,t];
	MAT[2,4,t+1] <- MAT[1,4,t]*t1	+ MAT[2,4,t]*(1-m2-t2-loss)			+ rec*MAT[2,3,t];
     # classe d'?ge xx
	MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
	MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
	MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec)			+ lat*MAT[3,2,t];
	MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss)			+ rec*MAT[3,3,t];
     # calcul des effectifs par ?tat de sant?
	MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); 
	MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); 
	MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); 
	MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
	nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N

    }# fin boucle temps

    # sorties ponctuelles ? analyser
    # XX
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    # xx
    sortie2 <- nvinf[temps]
    # xx
    sortie3 <- max(MAT[4,3,1:temps])
    # xx
    sortie4 <- sum(nvinf[1:365])
    
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
  }# fin boucle sc?narios AS
  return(sorties)
} # fin fonction du mod?le

# END


# Le non des parametres
nom = c(
  "K",
  "sr",
  "m1",
  "m2",
  "m3",
  "f2",
  "f3",
  "portee",
  "t1",
  "t2",
  "trans",
  "lat",
  "rec",
  "loss",
  "madd"
)

# Les valeurs nominales de base
ValNominale = c(
  110,
  0.5,
  0.0014,
  0.00029,
  0.0019,
  0.0019,
  0.0082,
  5,
  1 / 365,
  1 / 365,
  0.3,
  1 / 5,
  1 / 20,
  1 / 100,
  0.001
)


# Utilisation de la fonction modAppli avec les valeurs nominales
ValNominale = matrix(ValNominale, nrow = 1)
colnames(ValNominale) = nom
modAppli(parametre = ValNominale)



# Methode OAT -------------------------------------------------------------

# Création d'une matrix pour les 10 valeurs des parametres
# Ajout 10 lignes vides
vides <- matrix(NA, nrow = 10, ncol = ncol(ValNominale))
colnames(vides) <- colnames(ValNominale)
matrix_para <- rbind(ValNominale, vides)
colnames(matrix_para) = colnames(ValNominale)

# Valeur de variation des parametres
increment = c(
  1,
  0.05,
  0.0001,
  0.00001,
  0.0002,
  0.0002,
  0.0003,
  0.5,
  0.0002,
  0.0002,
  0.03,
  0.01,
  0.0005,
  0.0005,
  0.0004
)
# Boucle pour ajouter les valeur dans la matrix de variation des parametres
for (i in 2:nrow(matrix_para)) {
  matrix_para[i, 1] <- matrix_para[i - 1, 1] + increment[1]
  matrix_para[i, 2] <- matrix_para[i - 1, 2] + increment[2]
  matrix_para[i, 3] <- matrix_para[i - 1, 3] + increment[3]
  matrix_para[i, 4] <- matrix_para[i - 1, 4] + increment[4]
  matrix_para[i, 5] <- matrix_para[i - 1, 5] + increment[5]
  matrix_para[i, 6] <- matrix_para[i - 1, 6] + increment[6]
  matrix_para[i, 7] <- matrix_para[i - 1, 7] + increment[7]
  matrix_para[i, 8] <- matrix_para[i - 1, 8] + increment[8]
  matrix_para[i, 9] <- matrix_para[i - 1, 9] + increment[9]
  matrix_para[i, 10] <- matrix_para[i - 1, 10] + increment[10]
  matrix_para[i, 11] <- matrix_para[i - 1, 11] + increment[11]
  matrix_para[i, 12] <- matrix_para[i - 1, 12] + increment[12]
  matrix_para[i, 13] <- matrix_para[i - 1, 13] + increment[13]
  matrix_para[i, 14] <- matrix_para[i - 1, 14] + increment[14]
  matrix_para[i, 15] <- matrix_para[i - 1, 15] + increment[15]
}
matrix_para # Matrix contenant les valeurs nominale en 1ere ligne et les 10 valuers de chaque parametre (lignes 2:11)




# Fonction pour tester toutes les variations de parametre selectionne sachant que les autes son fixe: OAT
oat = function(matrices, parametres) {
  mat = matrices # recuperer la matrix globale
  
  # Creer les matrix vide a remplir
  sensibilite = matrix(NA, nrow = 1, ncol = 4)
  elasticite = matrix(NA, nrow = 1, ncol = 4)
  
  # Avec les valeur nominale
  val <- mat[1,]
  val <- as.matrix(val, nrow = 1)
  val <- t(val)
  resultat = as.matrix(modAppli(parametre = val))
  
  # Boucle pour caculer Y pour chaque variation du parametres concernées sachant que les autres sont fixes
  for (i in 2:nrow(mat)) {
    ligne = mat[1, ]
    ligne = as.matrix(ligne, nrow = 1)
    ligne = t(ligne)
    
    para = mat[i, parametres] # recuperer la nouvelle valeur du parametres
    
    dp = para - ligne[1, parametres] # calculer la variation du parametres: Xi - X0
    
    ligne[1, parametres] = para
    
    res = modAppli(parametre = ligne)
    
    resultat = rbind(resultat, res)
    
    
    # Calculez la sensibilité
    snblt <- (res - resultat[1, ]) / dp
    sensibilite = rbind(sensibilite, snblt)
    
    
    # Calculez l'élasticité
    X = mat[1, parametres]
    Y = resultat[1, ]
    rapport = X / Y
    
    elast <- sensibilite[i, ] * rapport
    elasticite = rbind(elasticite, elast)
    
  }
  
  nom_para = colnames(mat)[parametres]
  
  
  return(
    list(
      Parametres = nom_para,
      Resultat = resultat,
      Sensibilite = sensibilite[-1, ],
      Elasticite = elasticite[-1, ]
    )
  )
} # fin de la fonction


oat(matrices = parametres, parametres = 2)



# Créez un dataframe vide pour stocker les résultats
data_E <-
  data.frame(
    P = character(0),
    E1 = numeric(0),
    E2 = numeric(0),
    E3 = numeric(0),
    E4 = numeric(0)
  )

# boucle pour enrregistrer les return de oat dans un df
for (t in 1:15) {
  test <- oat(matrices = matrix_para, parametres = t)
  
  # Créez un dataframe avec les résultats pour le paramètre t
  data_elas <- data.frame(
    P = as.factor(rep(test$Parametres, 10)),
    E1 = as.numeric(test$Elasticite[, 1]),
    E2 = as.numeric(test$Elasticite[, 2]),
    E3 = as.numeric(test$Elasticite[, 3]),# Est-ce une coquille ? Devrait-il être E3 ou E4 ?OUI C'ETAIT UNE COQUILLE
    E4 = as.numeric(test$Elasticite[, 4])
  )
  
  # Ajoutez les résultats au dataframe principal
  data_E <- rbind(data_E, data_elas)
}

# Convertir les elasticité en valeur absolu
data_E[,2]= abs(data_E[,2])
data_E[,3]= abs(data_E[,3])
data_E[,4]= abs(data_E[,4])
data_E[,5]= abs(data_E[,5])

plot(data_E$P, data_E$E1, xlab = "Parametres", ylab = "Sensibilité relative de Sortie 1")
plot(data_E$P, data_E$E2, xlab = "Parametres", ylab = "Sensibilité relative de Sortie 2")
plot(data_E$P, data_E$E3, xlab = "Parametres", ylab = "Sensibilité relative de Sortie 3")
plot(data_E$P, data_E$E4, xlab = "Parametres", ylab = "Sensibilité relative de Sortie 4")


# QUESTION D:

#les parametres rec, loss, trans et lat sont les plus influents sur la sortie 1
#les parametres K et loss, sont les plus influents sur la sortie 2
#les parametres rec, trans, K et lat sont les plus influents sur la sortie 3
#les parametres K, loss et trans sont les plus influents sur la sortie 4

# Limites possibles:
# L'analyse de sensibilité OAT ne permet pas de voire les effets des interactions entre les parametres
# L'analyse de sensibilité OAT ne donne qu'une information valable uniquement autour des valeurs nominales appliquée. Elle ne renseigne en rien de
      #l’influence des paramètres d’entrée sur l’ensemble de leurs domaines de variabilité.
# L'analyse de sensibilité OAT n'est pas trés adapté a notre modele non lineaire car l'OAT est une méthode locale dont sa forme la plus simple donne des
# dérivée partielle du modèle par rapport aux parametres d'entrée







# Methode MORRIS ----------------------------------------------------------

library(sensitivity)

sorti_morris <- morris(
  model = modAppli,
  factors = nom, # le nom des parametres
  r = 5,
  design = list(
    type = "oat",n=100,
    levels = 6,
    grid.jump = 3
  ),
  binf = parametres[1, ], # le minimum des para
  bsup = parametres[11, ] # le maximum des para
)


# Vérification
sorti_morris$X[1,]
modAppli(parametre = t(as.matrix(sorti_morris$X[1,],nrow=1)))
sorti_morris$y[1,]
sorti_morris$X
sorti_morris$y


sorti_morris$de
sorti_morris$factors
sorti_morris$ee


plot(sorti_morris)
plot(sorti_morris, y_col = 1)
plot(sorti_morris, y_col = 2)
plot(x, y_col = 3)
plot(x, y_col = 4)




# Reccuper la sortie de Morris
mu <- apply(sorti_morris$ee, 3, function(M){
  apply(M, 2, mean)
})
mu.star <- apply(abs(sorti_morris$ee), 3, function(M){
  apply(M, 2, mean)
})
sigma <- apply(sorti_morris$ee, 3, function(M){
  apply(M, 2, sd)
})

mu.star= as.data.frame(mu.star)
sigma= as.data.frame(sigma)
row.names(sigma)

# Sortie 1
plot(sorti_morris, y_col = 1)
plot(0.1, 0.1, type="n",ylim=c(0, 0.005),xlim=c(0, 0.06),xlab="mu.star",ylab="sigma")
polygon(x = c(0,0,100),y = c(0,100,0),col = "blue",border = NA)
polygon(x = c(0,0,0.1),y = c(0,100,0),col = "yellow",border = NA)
polygon(x = c(0.1,100,0.1),y = c(0,0.1,0.1),col = "orange",border = NA)
points(mu.star$ycol1, sigma$ycol1)
text(mu.star$ycol1, sigma$ycol1, labels = row.names(sigma), pos = 3)

# Sortie 2
plot(sorti_morris, y_col = 2)
plot(0.1, 0.1, type="n",ylim=c(0, 0.040),xlim=c(0, 0.30),xlab="mu.star",ylab="sigma")
polygon(x = c(0,0,100),y = c(0,100,0),col = "blue",border = NA)
polygon(x = c(0,0,0.1),y = c(0,100,0),col = "yellow",border = NA)
polygon(x = c(0.1,100,0.1),y = c(0,0.01,0.01),col = "orange",border = NA)
points(mu.star$ycol2, sigma$ycol2)
text(mu.star$ycol2, sigma$ycol2, labels = row.names(sigma), pos = 3)

# Sortie 3
plot(sorti_morris, y_col = 3)
plot(0.1, 0.1, type="n",ylim=c(0, 1.4),xlim=c(0, 12),xlab="mu.star",ylab="sigma")
polygon(x = c(0,0,100),y = c(0,100,0),col = "blue",border = NA)
polygon(x = c(0,0,0.1),y = c(0,100,0),col = "yellow",border = NA)
polygon(x = c(0.1,900,0.1),y = c(0,0,0.1),col = "orange",border = NA)
points(mu.star$ycol3, sigma$ycol3)
text(mu.star$ycol3, sigma$ycol3, labels = row.names(sigma), pos = 3)


# Sortie 4
plot(sorti_morris, y_col = 4)
plot(0.1, 0.1, type="n",ylim=c(0, 12),xlim=c(0, 90),xlab="mu.star",ylab="sigma")
polygon(x = c(0,0,200),y = c(0,100,0),col = "blue",border = NA)
polygon(x = c(0,0,0.1),y = c(0,5000,0),col = "yellow",border = NA)
polygon(x = c(0.1,5000,0.1),y = c(0,0,0.1),col = "orange",border = NA)
points(mu.star$ycol4, sigma$ycol4)
text(mu.star$ycol4, sigma$ycol4, labels = row.names(sigma), pos = 3)


#- entrées ayant des effets négligeables (en jaune),
#– entrées ayant des effets linéaires et sans interaction (orange),
#– entrées ayant des effets non linéaires et/ou avec interactions (sans distinction de ces deux en bleu)




# Méthode FAST ------------------------------------------------------------


# Bornes des intervalles d’incertitude (lois uniformes)

borne_para <- apply(matrix_para,
                     2,
                     function(x) {
                       list(mean = mean(x), sd = sd(x))
                     })


borne_para <- apply( cbind(matrix_para[1,],
                           matrix_para[11,]),
                     1,
                     function(x){list(min=x[1],max=x[2])} )

# QUESTION b: Cas a 100 scenario
start_time100 <- Sys.time()

fast100 <- fast99(model=modAppli,
                      factors=nom,
                      n=100,
                      q=rep("qunif",15),
                      q.arg=borne_para)


end_time100 <- Sys.time() 
execution_time <- end_time100 - start_time100
execution_time # un petit doute ça ne dure que 0.780349 secs


plot(fast100) 
# Loss et trans ont des effet plus influents que les autres parametres. Faible effet des interactions


# les sortie
fast100$X
fast100$y
fast100$V


echant_100= as.data.frame(fast100$X)# Recuperer les valeurs echantillonnées
plot(echant_100$m1)
nrow(echant_100)

s_100=modAppli(parametre = echant_100) # utilisée les valeur echantionnée par fast99
hist(s_100[,1])
hist(s_100[,2])
hist(s_100[,3])
hist(s_100[,4])

library(dplyr)
library(tidyverse)

hist(s_100[,1], freq = FALSE, main = "")
density(s_100[,1]) |> lines(col = "blue")
title(main = "Histogramme et densité")
is.na(colSums(s_100[,1]))
drop_na(s_100)
str(s_100)


# Cas a 1000 scenario
start_time1000 <- Sys.time()

fast1000 <- fast99(model=modAppli,
                  factors=nom,
                  n=1000,
                  q=rep("qnorm",15),
                  q.arg=borne_para)


end_time1000 <- Sys.time() 
execution_time <- end_time1000 - start_time1000
execution_time # 

tell()

plot(fast1000) 
# 


# les sortie
fast1000$X
fast1000$y
fast1000$V


echant_1000= fast1000$X # Recuperer les valeurs echantillonnées
plot(echant_1000$K)
nrow(echant_1000)

s_1000=modAppli(parametre = echant_1000) # utilisée les valeur echantionnée par fast99
plot(s_1000[,1])
plot(s_1000[,2])
plot(s_1000[,3])
plot(s_1000[,4])

stock_1000= modAppli(echant_1000)
hist(stock_1000[,1])
hist(stock_1000[,2])
hist(stock_1000[,3])
hist(stock_1000[,4])

qunif(p = 0.5,min = 1,max = 4)


# Scénario ----------------------------------------------------------------

# si les individus infectieux atteigne 50% de la pop, les réduire de 50%
# Modiifcation de la fonctions
modAppli_2 <- function(parametre){  
  
  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (en jours)
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de param?tres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4)
  
  # boucle des sc?narios de l'?chantillonnage de l'AS
  for (i in 1:nrow(parametre)) { 
    
    # STRUCTURE & PARAMETRES DU MODELE
    
    # Récuperation des para
    K = parametre[i,1];    # Capacité de charge de l'environnement
    sr = parametre[i,2];   # Taux de survie
    m1 = parametre[i,3];   # Taux de mortalité classe 1
    m2 = parametre[i,4];   # Taux de mortalité classe 1
    m3 = parametre[i,5];   # Taux de mortalité classe 3
    f2 = parametre[i,6];   # Taux de fécondité classe 2
    f3 = parametre[i,7];   # Taux de fécondité classe 3
    portee = parametre[i,8]; # Portée (nombre d'individus produits par génération)
    t1 = parametre[i,9];   # Taux de transmission classe 1
    t2 = parametre[i,10];  # Taux de transmission classe 2
    
    # Transitions entre les états de santé
    trans = parametre[i,11];  # Taux de transition de XX
    lat = parametre[i,12];    # Taux de latence
    rec = parametre[i,13];    # Taux de récupération
    loss = parametre[i,14];   # Taux de perte
    madd = parametre[i,15];   # Taux d'addition de nouveaux individus XX
    eff = parametre[i,16];    # Efficacité du test
    seuil = parametre[i,17];  # Seuil de mortalité
    
    
    
    
    # INITIALISATION
    MAT <- array(0, dim=c(4,5,temps)); # nb indiv par classe d'?ge en ligne (derni?re ligne = pop tot), ?tat de sant? en colonne, pas de temps (dimension 3)
    nvinf <- array(0, dim=c(temps));
    
    # conditions initiales (la population est ? sa structure d'?quilibre, calcul?e par ailleurs)
    MAT[1,1,1] <- 27; # xx
    MAT[2,1,1] <- 23; # xx
    MAT[3,1,1] <- 36; # xx
    MAT[3,3,1] <- 1;  # xx
    # effectifs par ?tat de sant?
    MAT[4,1,1] <- sum(MAT[1:3,1,1]); 
    MAT[4,2,1] <- sum(MAT[1:3,2,1]); 
    MAT[4,3,1] <- sum(MAT[1:3,3,1]); 
    MAT[4,4,1] <- sum(MAT[1:3,4,1]);
    MAT[4,5,1] <- sum(MAT[1:3,5,1]);
    
    
    # SIMULATIONS
    # boucle du temps
    for (t in 1:(temps-1)) { 
      
          
          #  si mortalite en dessous du seui
          if (MAT[4,3,t]*madd < seuil*sum(MAT[4,,t])) { 
            
            
          # classe d'?ge xx
          # RQ : les naissances sont XX, les nouveaux n?s ?tant dans l'?tat XX
          N <- sum(MAT[4,,t]);	# taille de la pop en t
          MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t]      + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
          MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			  + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
          MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec)  		  + lat*MAT[1,2,t]; 
          MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) 		  + rec*MAT[1,3,t]; 
          # classe d'?ge xx
          MAT[2,1,t+1] <- MAT[1,1,t]*t1	+ MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
          MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
          MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-rec)		+ lat*MAT[2,2,t];
          MAT[2,4,t+1] <- MAT[1,4,t]*t1	+ MAT[2,4,t]*(1-m2-t2-loss)			+ rec*MAT[2,3,t];
          # classe d'?ge xx
          MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
          MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
          MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec)			+ lat*MAT[3,2,t];
          MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss)			+ rec*MAT[3,3,t];
          # calcul des effectifs par ?tat de sant?
          MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); 
          MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); 
          MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); 
          MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
          nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N
          
        } 
          else {   # Quarantaine si la mortalité lié a la mort depasse le suil 
        
        # classe d'?ge xx
        # RQ : les naissances sont XX, les nouveaux n?s ?tant dans l'?tat XX
        N <- sum(MAT[4,,t]);	# taille de la pop en t
        MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t]      + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
        MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			  + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
        MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-eff -rec)  		  + lat*MAT[1,2,t]; 
        MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-madd-t1-rec) + eff*MAT[1,3,t];
        MAT[1,5,t+1] <- MAT[1,5,t]*(1-m1-t1-loss) 		  + rec*(MAT[1,3,t] + MAT[1,4,t]); 
        # classe d'?ge xx
        MAT[2,1,t+1] <- MAT[1,1,t]*t1	+ MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
        MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
        MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-eff -rec)		+ lat*MAT[2,2,t];
        MAT[2,4,t+1] <- MAT[1,4,t]*t1 + MAT[2,4,t]*(1-m1-madd-t2-rec) + eff*MAT[2,3,t];
        MAT[2,5,t+1] <- MAT[1,5,t]*t1	+ MAT[2,5,t]*(1-m2-t2-loss)			+ rec*(MAT[2,3,t]+MAT[2,4,t]);
        # classe d'?ge xx
        MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
        MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
        MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-eff -rec)			+ lat*MAT[3,2,t];
        MAT[3,4,t+1] <- MAT[2,4,t]*t1 + MAT[3,4,t]*(1-m1-madd-rec) + eff*MAT[3,3,t];
        MAT[3,5,t+1] <- MAT[2,5,t]*t2	+ MAT[3,5,t]*(1-m3-loss)			+ rec*(MAT[3,3,t]+MAT[3,4,t]);
        # calcul des effectifs par ?tat de sant?
        MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); 
        MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); 
        MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); 
        MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
        MAT[4,5,t+1] <- sum(MAT[1:3,5,t+1]);
        nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N
        
      }# fin boucle quarantaine
          
    } # fin du boucle temps
        
    # sorties ponctuelles ? analyser
    # XX
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    # xx
    sortie2 <- nvinf[temps]
    # xx
    sortie3 <- max(MAT[4,3,1:temps])
    # xx
    sortie4 <- sum(nvinf[1:365])
    
    sortie5 <- data.frame( S= MAT[4,1,1:temps], L =MAT[4,2,1:temps], I= MAT[4,3,1:temps] , Q= MAT[4,4,1:temps] ,R=MAT[4,5,1:temps])
    
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
    
  }# fin boucle sc?narios AS
  return(sortie5)
} # fin fonction du mod?le


  
  
ValNominale = c(100, 0.5, 0.0014, 0.00029, 0.0019, 0.0019, 0.0082, 
                5, 1/365, 1/365, 0.3, 1/5, 1/20, 1/100, 0.01, 0.7, 0.0005)
  
modAppli_2(parametre = matrix(ValNominale, nrow=1))

