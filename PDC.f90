PROGRAM Ecoulement_insta

IMPLICIT NONE

!-----------------------------------------------------------------------------------------
!						       Déclaration des variables
!-----------------------------------------------------------------------------------------
INTEGER, PARAMETER:: FICH=125
INTEGER, PARAMETER:: FICH2=126
DOUBLE PRECISION, PARAMETER:: Pi=4*ATAN(1.D0)
INTEGER:: i, j,m, Mt, Nt, n, nb, compteur0, compteur2, compteur4, compteur6,compteur8,compteur10
DOUBLE PRECISION:: longueur,Rayon, epaisseur2, epaisseur3, tsauvegarde, compteur, temps, dt, duree,dx
DOUBLE PRECISION:: Tamb, rho,rhop, Ceau, Cper, Tprep, hair, heau, vitesse, attente
DOUBLE PRECISION:: lambda1, lambda2, lambda3, S1, S2, V1, V2, Slat1, Slat2, Slat3
DOUBLE PRECISION:: F1, F2, F3, F4, F5, F6, F7, compteur3, tconsigne,tu, tempsf, trefroi, Tparoi, Rci, Rcp
DOUBLE PRECISION:: tu1, tu2, tu3, tu4, tu5, duree2, duree3, duree4, duree5, compteur5,compteur7,compteur9,compteur11
DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: T1, T2

!-----------------------------------------------------------------------------------------
!						       Importation des données du problème
!-----------------------------------------------------------------------------------------
CALL load()
PRINT*,'Chargement des données.................OK'
!-----------------------------------------------------------------------------------------
!						       Initialisation des variables
!-----------------------------------------------------------------------------------------
ALLOCATE(T1(1:Mt))
ALLOCATE(T2(1:Mt))

!Etat refroidi: toute la conduite est à Tambiant (K)
T1(:)=Tamb
T2(:)=Tamb
dt=duree/Nt
Tparoi=30
nb=0
temps=0
compteur0=0
compteur=0
compteur2=0
compteur3=0
compteur4=0
compteur5=0
compteur6=0
compteur7=0
compteur8=0
compteur9=0
compteur10=0
compteur11=0
attente=0  !Par défaut indique que la température seuil n'a pas été atteinte

OPEN(FICH2,FILE="result_PDC.txt", ACTION="WRITE", STATUS="UNKNOWN");
  WRITE(FICH2,*)
CLOSE(FICH2)

!-----------------------------------------------------------------------------------------
!						               Discrétisation spatiale
!-----------------------------------------------------------------------------------------

CALL Discretisation()
PRINT*,'Discrétisation.........................OK'

!Contrainte sur le pas de temps
PRINT*, "Il faut:", vitesse*S2*dt, "<<", V1
PRINT*, "Sinon l'hypothèse que la température dans l'écoulement reste constante"
PRINT*, "en un temps dt n'est pas vérifiée."
!-----------------------------------------------------------------------------------------
!						                   Calcul
!-----------------------------------------------------------------------------------------
CALL Calcul()
PRINT*,'Calcul.................................OK'
!Temps d'attente pour atteindre la valeur consigne
PRINT*, "Pour obtenir une température (pour le dernier profil) en sortie de", Tconsigne,"°C on devra attendre", attente  ,"s"
PRINT*, "(Si 0s est affiché c'est que la température n'est pas atteinte durant ce temps de simulation)."
PRINT*, "L'ECS est refroidi à ", Trefroi,"°C après un temps de refroidissement de ", Tempsf, "s soit ", Tempsf/60, "minutes"
PRINT*, "(si 0s est affiché c'est qu'on n'atteint pas la temperature considérée comme froide pendant cette simulatiton)"
!-----------------------------------------------------------------------------------------
!						                 Exportation
!-----------------------------------------------------------------------------------------
CALL export()
PRINT*,'Exportation............................OK'

CONTAINS

SUBROUTINE Calcul()
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------
    !						Discretisation temporelle [Explicite]
    !-----------------------------------------------------------------------------------------

    DO WHILE (temps<duree) !Condition d'arrêt

      DO m=1,Mt

          !Flux gauche
          IF (m>1) THEN
            F1=S1*lambda1*(T1(m-1)-T1(m))/dx
            F3=S2*lambda2*(T2(m-1)-T2(m))/dx
          ELSE
            F1=2*S1*lambda1*(Tprep-T1(m))/dx
            F3=2*S2*lambda2*(Tparoi-T2(m))/dx
          END IF

          !Flux droit
          IF (m<Mt) THEN
            F2=S1*lambda1*(T1(m+1)-T1(m))/dx
            F4=S2*lambda2*(T2(m+1)-T2(m))/dx
          ELSE
            F2=(Tamb-T1(Mt))/(dx/(2*lambda1*S1)+1/(hair*S1))  !Attention convection
            F4=(Tamb-T2(Mt))/(dx/(2*lambda2*S2)+1/(hair*S2))
          END IF

          !Flux conducto-convectifs
          F6=(T1(m)-T2(m))/(epaisseur2/(2*lambda2*Slat1)+rayon/(rayon*heau*Slat1+lambda1*Slat1))

          !Flux avec la résistance cylindrique correspondant à l'isolant
          Rcp=log((rayon+epaisseur2)/(rayon+epaisseur2/2))/(2*pi*lambda2*dx)
          Rci=log((rayon+epaisseur2+epaisseur3)/(rayon+epaisseur2))/(2*pi*lambda3*dx)
          F7=(T2(m)-Tamb)/(Rcp+Rci+1/(hair*Slat3))

          !Flux advectif
          IF (m==1) THEN
              F5=vitesse*rho*S1*Ceau*(Tprep-T1(1))
          ELSE
              F5=vitesse*rho*S1*Ceau*(T1(m-1)-T1(m))
          END IF

          !Bilan sur le PER
          T2(m)=dt*(F3+F4+F6-F7)/(rhop*Cper*V2)+T2(m)
          !Bilan sur l'écoulement
          T1(m)=dt*(F1+F2+F5-F6)/(rho*Ceau*V1)+T1(m)

      END DO
      !--------------FIN DE LA BOUCLE D'ESPACE----------------------
      !Condition pour obtenir le temps de refroidissement
      IF(T1(m)<Trefroi .AND. compteur0==0) THEN
        tempsf=temps
        compteur0=1
      END IF

      !Conditions pour définir le profil de puisage
      IF (temps>duree2 .AND. temps<duree3) THEN
        vitesse=0.9
        CALL profil(tu2,compteur4,compteur5)
      ELSE IF (temps>duree3 .AND. temps<duree4) THEN
        vitesse=0.9
        CALL profil(tu3,compteur6,compteur7)
      ELSE IF (temps>duree4 .AND. temps<duree5) THEN
        vitesse=0.9
        CALL profil(tu4,compteur8,compteur9)
      ELSE IF (temps>duree5) THEN
        vitesse=0.9
        CALL profil(tu5,compteur10,compteur11)
      ELSE
        CALL profil(tu1,compteur2,compteur3)
      END IF

      temps=temps+dt
      compteur=compteur+dt

      !-------------------------------------------------------------------------------------
      !						Ecriture de certain résultats dans un fichier
      !-------------------------------------------------------------------------------------
      IF (compteur>tsauvegarde-dt .AND. compteur<=tsauvegarde+dt) THEN
        nb=nb+1 ! nombre de données sauvegardées
        CALL export()
        compteur=0
        !PRINT*, "F1=",F1,"F3=",F3,"F2=",F2,"F4=",F4,"F5=",F5,"F6=",F6
      END IF

    END DO
    !--------------FIN DE LA BOUCLE DE TEMPS---------------------
END SUBROUTINE

SUBROUTINE profil(tu,compt,compt2)
  IMPLICIT NONE
  INTEGER:: compt
  DOUBLE PRECISION:: tu, compt2

  !Condition pour récupérer le temps d'attente pour atteindre Tconsigne
  IF (T1(Mt)>Tconsigne .AND. compt==0) THEN
    attente=temps
    compt=1
  END IF

  !Condition pour compter le temps d'utilisation
  IF (T1(Mt)>Tconsigne) THEN
    compt2=compt2+dt
  END IF

  !Condition pour couper l'ecoulement à Tutilisation1
  IF (compt2>tu) THEN
    vitesse=0
  END IF
END SUBROUTINE

SUBROUTINE load()
	IMPLICIT NONE

	OPEN(FICH,FILE="donnees_PDC_profil.f90", ACTION="READ");

	! ------------------------------------------------------------------------------------
	! Lecture des données globales
	! ------------------------------------------------------------------------------------

  !Géométrie et  maillage
	READ(FICH,*) Mt
	READ(FICH,*) Rayon
  READ(FICH,*) longueur
  READ(FICH,*) epaisseur2
  READ(FICH,*) epaisseur3

  !Discrétisation spatiale
  READ(FICH,*) Nt
  READ(FICH,*) duree
  READ(FICH,*) tsauvegarde

  !Propriété de l'ECS en écoulement
  READ(FICH,*) vitesse
	READ(FICH,*) rho
  READ(FICH,*) Ceau
  READ(FICH,*) Cper
  READ(FICH,*) lambda1
  READ(FICH,*) heau

  !Propriété du PER
  READ(FICH,*) lambda2
  READ(FICH,*) rhop

  !Propriétés de l'isolant et de l'air
  READ(FICH,*) lambda3
  READ(FICH,*) hair
  READ(FICH,*) Tamb
  READ(FICH,*) Tprep

  !Température consigne Temps d'utilisation
  READ(FICH,*)Tconsigne
  READ(FICH,*)Tu1
  READ(FICH,*)Tu2
  READ(FICH,*)Tu3
  READ(FICH,*)Tu4
  READ(FICH,*)Tu5
  READ(FICH,*)duree2
  READ(FICH,*)duree3
  READ(FICH,*)duree4
  READ(FICH,*)duree5
  READ(FICH,*)Trefroi

	CLOSE(FICH)
END SUBROUTINE load

SUBROUTINE Discretisation
  IMPLICIT NONE

  dx=longueur/Mt
  S1=pi*rayon**2
  S2=pi*((epaisseur2+rayon)**2-rayon**2)
  Slat1=2*pi*rayon*dx
  Slat2=2*pi*(rayon+epaisseur2)*dx
  Slat3=2*pi*(rayon+epaisseur2+epaisseur3)*dx
  V1=S1*dx
  V2=S2*dx
END SUBROUTINE Discretisation

SUBROUTINE export()
	IMPLICIT NONE

	OPEN(FICH2,FILE="result_PDC.txt", ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND");

  ! ------------------------------------------------------------------------------------
  ! Ecriture des données
  ! ------------------------------------------------------------------------------------
    !WRITE(FICH2,*) 'temps=', temps
    !WRITE(FICH2,*) 'temps=', temps, 'Valeur n°', nb, 'au niveau de l épaisseur'
    WRITE(FICH2,*) temps,T1(Mt),T2(Mt)

  CLOSE(FICH2)
END SUBROUTINE export

END PROGRAM
