### Mod?le dont la sensibilit? doit ?tre analys?e dans le cadre du projet MODE-MPI 2023-2024

### Le mod?le est ici d?fini sous forme de fonction pour faciliter vos analyses de sensibilit? (AS)
### La fonction renvoie les sorties ponctuelles qui sont ? analyser dans l'AS

modAppli <- function(parametre){  

  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (en jours)
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de param?tres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4)
  List_mat = list()
  # boucle des sc?narios de l'?chantillonnage de l'AS
  for (i in 1:nrow(parametre)) { 

    # STRUCTURE & PARAMETRES DU MODELE

    # XX
    K = parametre[i,1];		# xx
    sr = parametre[i,2];	# xx
    m1 = parametre[i,3];	# xx
    m2 = parametre[i,4];	# xx
    m3 = parametre[i,5];	# xx
    f2 = parametre[i,6];	# xx
    f3 = parametre[i,7];	# xx
    portee = parametre[i,8];	# xx
    t1 = parametre[i,9];	# xx
    t2 = parametre[i,10];	# xx

    # XX
    trans = parametre[i,11]; # xx
    lat = parametre[i,12];	# xx
    rec = parametre[i,13];	# xx
    loss = parametre[i,14];	# xx
    madd = parametre[i,15];	# xx

    # INITIALISATION
    
    MAT <- array(0, dim=c(4,4,temps)); # nb indiv par classe d'?ge en ligne (derni?re ligne = pop tot), ?tat de sant? en colonne, pas de temps (dimension 3)
    nvinf <- array(0, dim=c(temps));
    
    # conditions initiales (la population est ? sa structure d'?quilibre, calcul?e par ailleurs)
    MAT[1,1,1] <- 27; # xx
    MAT[2,1,1] <- 23; # xx
    MAT[3,1,1] <- 36; # xx
    MAT[3,3,1] <- 1;  # xx
    # effectifs par ?tat de sant?
    MAT[4,1,1] <- sum(MAT[1:3,1,1]); MAT[4,2,1] <- sum(MAT[1:3,2,1]); MAT[4,3,1] <- sum(MAT[1:3,3,1]); MAT[4,4,1] <- sum(MAT[1:3,4,1]);

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
	MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
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
    
    List_mat[[i]] = MAT
    
    
  } # fin boucle sc?narios AS
  
  ensemble = list(sorties, List_mat)
  
  return(ensemble)
  
} # fin fonction du mod?le

# END


# Premi?re fonction graphique (partie 1), illustration de la premi?re simulation.

  # Arguments : 
      # classe --> dans [1:4] 1 = classe 1, 2 = classe 2, 3 = classe 3, 4 = population totale
      # Titre --> Cha?ne de caract?re

plot_initial = function(classe, titre) { # Faire en ggplot2
  
  plot(1:temps, predictions[classe,1,], type = "l", col = "yellow", cex = 2, 
       xlab = "temps (jours)", ylab = "Effectif", main = titre)
  lines(predictions[classe,2,], col = "orange", cex = 2)
  lines(predictions[classe,3,], col = "red", cex = 2)
  lines(predictions[classe,4,], col = "green", cex = 2)
  lines(predictions[classe,4,]+predictions[classe,3,]+predictions[classe,2,])
  
    if (classe == 1){
      legend(x = "topright", legend = c("Susceptibles",
                                       "Latents", 
                                       "Infectieux", 
                                       "Immunis?s", 
                                       "Population totale"), 
             col = c("yellow", "orange", "red", "green", "black"),lty = 1,
            cex = 0.5)
  }
}


    ### Partie 2 : Analyse de sensibilit? OAT

# Fonction de variation de chaque param?tre.

  # Arguments : 
      # parametres --> Valeurs des param?tres initiaux (matrice 1x15)
      # increment --> variations de chaque param?tres par it?ration (vecteur de longueur 15)
  
  # Sorties : 
      # parametres --> matrice 10 x 15 contenant toutes les valeurs de param?tres 

var_para <- function(parametres, increment){
  
  vides <- matrix(NA, nrow = 10, ncol = ncol(parametres))
  matrix_para <- rbind(parametres, vides)
  colnames(matrix_para) = nom
  
  for (i in 2:nrow(matrix_para)) {
    for (j in 1:length(increment)) {
    matrix_para[i, j] <- matrix_para[i - 1, j] + increment[j]
    }
  }
  return(matrix_para)
}

# Fonction graphique pour OAT chaque param?tre.

  # Arguments : 

# Fonction pour tester toutes les variations de parametre selectionne sachant que les autes son fixe: OAT
oat = function(matrices, parametres) {
  mat = matrices # recuperer la matrix globale
  
  # Creer les matrix vide a remplir
  sensibilite = matrix(NA, nrow = 1, ncol = 4)
  elasticite = matrix(NA, nrow = 1, ncol = 4)
  
  # Avec les valeur nominale
  
  val <- t(as.matrix(mat[1,], nrow = 1))
  
  resultat = as.matrix(modAppli(parametre = val)[[1]])
  
  # Boucle pour caculer Y pour chaque variation du parametres concernées sachant que les autres sont fixes
  for (i in 2:nrow(mat)) {

    ligne = t(as.matrix(mat[1, ], nrow = 1))
    
    para = mat[i, parametres] # recuperer la nouvelle valeur du parametres
    
    dp = para - ligne[1, parametres] # calculer la variation du parametres: Xi - X0
    
    ligne[1, parametres] = para
    
    res = modAppli(parametre = ligne)[[1]]
    
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




      ### Partie 3 : Analyse de sensibilit? par la m?thode Morris

# ----

# Fonction d'extraction des sorties MORRIS \mu, \mu* et \sigma

sorties_morris = function(vecteur_ee){
  mu <- apply(vecteur_ee, 3, function(M){
    apply(M, 2, mean)
  })
  mu.star <- apply(abs(vecteur_ee), 3, function(M){
    apply(M, 2, mean)
  })
  sigma <- apply(vecteur_ee, 3, function(M){
    apply(M, 2, sd)
  })
}
