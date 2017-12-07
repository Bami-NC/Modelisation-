PROGRAM Ecoulement_insta

IMPLICIT NONE

!-----------------------------------------------------------------------------------------
!						       Déclaration des variables
!-----------------------------------------------------------------------------------------
INTEGER, PARAMETER:: FICH=125
INTEGER, PARAMETER:: FICH2=126
DOUBLE PRECISION, PARAMETER:: Pi=4*ATAN(1.D0)
INTEGER:: i, j,m, Mt, Nt, n, nb, compteur2
DOUBLE PRECISION:: longueur,Rayon, epaisseur2, epaisseur3, tsauvegarde, compteur, temps, dt, duree,dx
DOUBLE PRECISION:: Tamb, rho,rhop, Ceau, Cper, Tprep, hair, heau, vitesse, attente
DOUBLE PRECISION:: lambda1, lambda2, lambda3, S1, S2, V1, V2, Slat1, Slat2, Slat3
DOUBLE PRECISION:: F1, F2, F3, F4, F5, F6, F7
DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: T1, T2

!-----------------------------------------------------------------------------------------
!						       Importation des données du problème
!-----------------------------------------------------------------------------------------
CALL load()
PRINT*,'Chargement des données.................OK'
!-----------------------------------------------------------------------------------------
!						       Initialisation des variables
!-----------------------------------------------------------------------------------------
ALLOCATE(T1(1:Mt)) !On stocke les temperatures dans une matrice
ALLOCATE(T2(1:Mt)) !On stocke les temperatures dans une matrice

!Etat refroidi: toute la conduite est à Tambiant (K)
T1(:)=Tamb
T2(:)=Tamb
dt=duree/Nt
nb=0
temps=0
compteur=0
compteur2=0
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
PRINT*, "Il faut:", vitesse*S1*dt, "<<", V2
PRINT*, "Sinon l'hypothèse que la température dans l'écoulement reste constante"
PRINT*, "en un temps dt n'est pas vérifiée."
!-----------------------------------------------------------------------------------------
!						                   Calcul
!-----------------------------------------------------------------------------------------
CALL Calcul()
PRINT*,'Calcul.................................OK'
!Temps d'attente pour atteindre la valeur consigne
PRINT*, "Pour obtenir une température en sortie de 28°C on devra attendre", anint (attente ) ,"s"
PRINT*, "(Si 0s est affiché c'est que la température n'est pas atteinte durant ce temps de simulation)."
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

          IF (m>1) THEN
            F1=S1*lambda1*(T1(m-1)-T1(m))/dx
            F3=S2*lambda2*(T2(m-1)-T2(m))/dx
          ELSE
            F1=2*S1*lambda1*(Tprep-T1(1))/dx
            F3=2*S2*lambda2*(Tprep-T2(1))/dx
          END IF

          IF (m<Mt) THEN
            F2=S1*lambda1*(T1(m+1)-T1(m))/dx
            F4=S2*lambda2*(T2(m+1)-T2(m))/dx
          ELSE
            F2=2*S1*lambda1*(Tamb-T1(Mt))/dx
            F4=2*S2*lambda2*(Tamb-T2(Mt))/dx
          END IF

          !Flux conducto-convectifs
          F6=(T1(m)-T2(m))/(epaisseur2/(2*lambda2*Slat1)+rayon/(rayon*heau*Slat1+lambda1*Slat1))

          !Flux avec la résistance cylindrique correspondant à l'isolant
          F7=(T2(m)-Tamb)/(1/(hair*Slat3)+log((rayon+epaisseur2+epaisseur3)/(rayon+epaisseur2))/(2*pi*lambda3*dx))

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

      IF (T1(Mt)>28 .AND. compteur2==0) THEN
          attente=temps
          compteur2=1
          vitesse=0 !dès qu'on atteint 28 on coupe l'ECS pour observer le refroidissement
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

SUBROUTINE load()
	IMPLICIT NONE

	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	OPEN(FICH,FILE="donnees_PDC.f90", ACTION="READ");

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
    DO i=1,Mt
      WRITE(FICH2,*) i*dx,",",T1(i),",",T2(i),","
    END DO

  CLOSE(FICH2)
END SUBROUTINE export

END PROGRAM
