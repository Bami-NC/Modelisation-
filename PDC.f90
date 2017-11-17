PROGRAM Ecoulement_insta

IMPLICIT NONE

!-----------------------------------------------------------------------------------------
!						       Déclaration des variables
!-----------------------------------------------------------------------------------------
INTEGER, PARAMETER:: FICH=125
INTEGER, PARAMETER:: FICH2=126
DOUBLE PRECISION, PARAMETER:: Pi=3.14
INTEGER:: i, j,m, Mt, Nt, n, nb
DOUBLE PRECISION:: longueur,Rayon, epaisseur2, epaisseur3, tsauvegarde, compteur, temps, dt, duree,dx
DOUBLE PRECISION:: Tamb, rho, Ceau, Tprep, hair, heau, vitesse
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
dt=duree/Nt

!Conversion °C->K
Tamb=273.15+Tamb
Tprep=273.15+Tprep

ALLOCATE(T1(1:Mt)) !On stocke les temperatures dans une matrice
ALLOCATE(T2(1:Mt)) !On stocke les temperatures dans une matrice

!Etat refroidi: toute la conduite est à Tambiant (K)
T1(:)=Tamb
T2(:)=Tamb
nb=0
temps=0
compteur=0

OPEN(FICH2,FILE="result_PDC.f90", ACTION="WRITE", STATUS="UNKNOWN");
  WRITE(FICH2,*)
CLOSE(FICH2)

!-----------------------------------------------------------------------------------------
!						       Discrétisation spatiale
!-----------------------------------------------------------------------------------------

CALL Discretisation()
PRINT*,'Discrétisation.........................OK'
!-----------------------------------------------------------------------------------------
!						       Calcul
!-----------------------------------------------------------------------------------------

CALL Calcul()
PRINT*,'Calcul.................................OK'
!-----------------------------------------------------------------------------------------
!						       Exportation
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
            F1=S1*lambda1*(Tprep-T1(1))/dx
            F3=S2*lambda2*(Tprep-T2(1))/dx
          END IF

          IF (m<Mt) THEN
            F2=S1*lambda1*(T1(m+1)-T1(m))/dx
            F4=S2*lambda2*(T2(m+1)-T2(m))/dx
          ELSE
            F2=S1*lambda1*(Tamb-T1(Mt))/dx
            F4=S2*lambda2*(Tamb-T2(Mt))/dx
          END IF

          !Flux conducto-convectifs
          F6=(T2(m)-T1(m))/(rayon/lambda1+epaisseur2/2*lambda2+1/heau*Slat1)

          !Flux avec la résistance cylindrique correspondant à l'isolant
          F7=(Tamb-T2(m))/(1/hair*Slat3+log((rayon+epaisseur2+epaisseur3)/(rayon+epaisseur2))/2*pi*lambda3*dx)

          IF (m>1 .AND. m<Mt) THEN
              F5=vitesse*rho*S1*Ceau*(T1(m-1)-T1(m))
          ELSE IF (m==1) THEN
              F5=vitesse*rho*S1*Ceau*(Tprep-T1(1))
          ELSE IF (m==Mt) THEN
              F5=vitesse*rho*S1*Ceau*(T1(Mt)-Tamb)
          END IF

          !Bilan sur le PER()
          T2(m)=dt*(F3+F4+F6+F7)/(rho*Ceau*V2)+T2(m)
          !Bilan sur l'écoulement
          T1(m)=dt*(F1+F2+F5+F6)/(rho*Ceau*V1)+T1(m)
      END DO

      !-------------------------------------------------
      temps=temps+dt
      compteur=compteur+dt

      !-------------------------------------------------------------------------------------
      !						Ecriture de certain résultats dans un fichier
      !-------------------------------------------------------------------------------------
      IF (compteur>tsauvegarde-dt .AND. compteur<=tsauvegarde+dt) THEN
        nb=nb+1 ! nombre de données sauvegardées
        !Conversion en Celsius
        DO i=1,M
          T1(i)=
          T2(i)
        END DO
        CALL export()
        compteur=0
      END IF

    END DO
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
  READ(FICH,*) lambda1
  READ(FICH,*) heau

  !Propriété du PER
  READ(FICH,*) lambda2

  !Propriétés de l'isolant et de l'air
  READ(FICH,*) lambda3
  READ(FICH,*) hair
  READ(FICH,*) Tamb
  READ(FICH,*) Tprep

	CLOSE(FICH)
END SUBROUTINE load

SUBROUTINE Discretisation
  IMPLICIT NONE

  dx=longueur/M
  S1=pi*rayon**2
  S2=pi*((epaisseur2+rayon)**2-rayon**2)
  Slat1=2*pi*rayon*dx
  Slat3=2*pi*(rayon+epaisseur2+epaisseur3)*dx
  V1=S1*dx
  V2=S2*dx
END SUBROUTINE Discretisation

SUBROUTINE export()
	IMPLICIT NONE

	OPEN(FICH2,FILE="result_PDC.f90", ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND");

  ! ------------------------------------------------------------------------------------
  ! Ecriture des données
  ! ------------------------------------------------------------------------------------

    WRITE(FICH2,*) 'temps=', temps, 'Valeur n°', nb, 'au niveau de l épaisseur'
    DO i=1,M
      WRITE(FICH2,*) T2(i)
    END DO
    WRITE(FICH2,*) '---------------------------------------------------------'

    DO i=1,M
      WRITE(FICH2,*) T1(i)
    END DO
    WRITE(FICH2,*) '---------------------------------------------------------'

  CLOSE(FICH2)
END SUBROUTINE export

END PROGRAM
