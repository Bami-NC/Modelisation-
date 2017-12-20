PROGRAM Ecoulement_instaMCP

IMPLICIT NONE

!-----------------------------------------------------------------------------------------
!						       Déclaration des variables
!-----------------------------------------------------------------------------------------
INTEGER, PARAMETER:: FICH=148
INTEGER, PARAMETER:: FICH2=128
DOUBLE PRECISION, PARAMETER:: Pi=4*ATAN(1.D0)
INTEGER:: i, j,m, Mt, Nt, n, nb, compteur0, compteur2, compteur4, compteur6, compteur8, compteur10
DOUBLE PRECISION:: longueur,Rayon, epaisseur2, epaisseur3, epaisseurMCP, tsauvegarde, compteur, temps, dt, duree,dx
DOUBLE PRECISION:: Tamb, rho,rhop, Ceau, Cper, Tprep, hair, heau, vitesse, attente, Tparoi
DOUBLE PRECISION:: rhoMCP, CS, CL,TF, LF, HG
DOUBLE PRECISION:: lambda1, lambda2, lambda3, lambdaMCP, S1, S2, S3, V1, V2,V3, Slat1, Slat2, Slat3,Slat4
DOUBLE PRECISION:: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10 , tempsf, Tconsigne, Trefroi, Rcp, Rcm, Rci
DOUBLE PRECISION:: tu1, tu2, tu3, tu4, tu5, duree2, duree3, duree4, duree5, compteur3, compteur5,compteur7,compteur9,compteur11
DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: T1, T2, T3,h, Y

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
ALLOCATE(T3(1:Mt))
ALLOCATE(h(1:Mt))
ALLOCATE(Y(1:Mt))
!Etat refroidi: toute la conduite est à Tambiant (K)
T1(:)=Tamb
T2(:)=Tamb
T3(:)=Tamb
dt=duree/Nt
Tparoi=30 !Température de la paroi du préparateur
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

OPEN(FICH2,FILE="result_PDC2.txt", ACTION="WRITE", STATUS="UNKNOWN");
  WRITE(FICH2,*)
CLOSE(FICH2)

!-----------------------------------------------------------------------------------------
!						               Discrétisation spatiale
!-----------------------------------------------------------------------------------------

CALL Discretisation()
PRINT*,'Discrétisation.........................OK'

!Contrainte sur le pas de temps
PRINT*, "Il faut:", vitesse*S3*dt, "<<", V1
PRINT*, "Sinon l'hypothèse que la température dans l'écoulement reste constante"
PRINT*, "en un temps dt n'est pas vérifiée."
!-----------------------------------------------------------------------------------------
!						                   Calcul
!-----------------------------------------------------------------------------------------
CALL Calcul()
PRINT*,'Calcul.................................OK'
!Temps d'attente pour atteindre la valeur consigne et temps de refroidissement
PRINT*, "Pour obtenir une température en sortie (pour le dernier profil) de",Tconsigne,"°C on devra attendre", attente ,"s"
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

    !------       Condition initiale et calcul de h0       ---------------------
    h(:) = calcH(Tamb)		! on travaille en massique
    !-----------------------------------------------------------------------------------------
    !						Discretisation temporelle [Explicite]
    !-----------------------------------------------------------------------------------------

    DO WHILE (temps<duree) !Condition d'arrêt

      DO m=1,Mt
          !Flux gauches
          IF (m>1) THEN
            F1=S1*lambda1*(T1(m-1)-T1(m))/dx
            F3=S2*lambda2*(T2(m-1)-T2(m))/dx
            F9=S3*lambdaMCP*(T3(m-1)-T3(m))/dx
          ELSE
            F1=2*S1*lambda1*(Tprep-T1(m))/dx    !en contact avec le préparateur
            F3=2*S2*lambda2*(Tparoi-T2(m))/dx   !en contact avec la paroi du préparateur
            F9=2*S3*lambdaMCP*(Tparoi-T3(m))/dx
          END IF
          !Flux droits
          IF (m<Mt) THEN
            F2=S1*lambda1*(T1(m+1)-T1(m))/dx
            F4=S2*lambda2*(T2(m+1)-T2(m))/dx
            F10=S3*lambdaMCP*(T2(m+1)-T3(m))/dx
          ELSE
            F2=(Tamb-T1(Mt))/(dx/(2*lambda1*S1)+1/(hair*S1))  !Attention convection
            F4=(Tamb-T2(Mt))/(dx/(2*lambda2*S2)+1/(hair*S2))
            F10=(Tamb-T3(Mt))/(dx/(2*lambdaMCP*S3)+1/(hair*S3))
          END IF

          !Flux conducto-convectifs
          F6=(T1(m)-T2(m))/(epaisseur2/(2*lambda2*Slat1)+rayon/(rayon*heau*Slat1+lambda1*Slat1))

          !Flux avec la résistance cylindrique correspondant à l'isolant
          Rcm=log((rayon+epaisseur2+epaisseurMCP)/(rayon+epaisseur2+epaisseurMCP/2))/(2*Pi*lambdaMCP*dx)
          Rci=log((rayon+epaisseur2+epaisseur3+epaisseurMCP)/(rayon+epaisseur2+epaisseurMCP))/(2*pi*lambda3*dx)
          F7=(T3(m)-Tamb)/(Rcm+Rci+1/(hair*Slat3))

          !Flux advectif
          IF (m==1) THEN
              F5=vitesse*rho*S1*Ceau*(Tprep-T1(1))
          ELSE
              F5=vitesse*rho*S1*Ceau*(T1(m-1)-T1(m))
          END IF

          !Flux conductif entre le PER et le MCP
          Rcp=log((rayon+epaisseur2)/(rayon+epaisseur2/2))/(2*Pi*lambda2*dx)
          Rcm=log((rayon+epaisseur2+epaisseurMCP/2)/(rayon+epaisseur2))/(2*Pi*lambda3*dx)
          F8=(T2(m)-T3(m))/(Rcp+Rcm)

          !Bilan sur le PER
          T2(m)=dt*(F3+F4+F6-F8)/(rhop*Cper*V2)+T2(m)
          !Bilan sur l'écoulement
          T1(m)=dt*(F1+F2+F5-F6)/(rho*Ceau*V1)+T1(m)
          !Bilan sur le MCP
          h(m)=dt*(F8+F9+F10-F7)/(rhoMCP*V3)+h(m)

          T3(m) = calcT(h(m))
		      Y(m) = calcY(h(m))

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

      !-------------------------------------------------------------------------
      !						Ecriture de certain résultats dans un fichier
      !-------------------------------------------------------------------------
      IF (compteur>tsauvegarde-dt .AND. compteur<=tsauvegarde+dt) THEN
        nb=nb+1 ! nombre de données sauvegardées
        CALL export()
        compteur=0
        !PRINT*, "F1=",F1,"F3=",F3,"F2=",F2,"F4=",F4,"F5=",F5,"F6=",F6
        !PRINT*, "F8=",F8,"F9=",F9,"F10=",F10
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

  !Condition pour couper l'ecoulement à TutilisationX
  IF (compt2>tu) THEN
    vitesse=0
  END IF
END SUBROUTINE

FUNCTION calcH(T)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: T
	DOUBLE PRECISION :: calcH
  !Dans le cas de la paraffine on définit un cp moyen entre 28-32
  IF (T<32 .AND. T>28) THEN
    CS=26000
    CL=26000
  ELSE
    CS=3500
    CL=2500
  END IF

	IF (T<TF) THEN
		! SOLIDE
		calcH = CS*(T-TF)
	ELSEIF (T>=TF) THEN
		! LIQUIDE
		calcH = LF+CL*(T-TF)
	ELSE ! T=TF
		! FUSION
		WRITE(*,*)"ERREUR T=TF"
		STOP
	ENDIF
END FUNCTION

FUNCTION calcT(h)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: h
	DOUBLE PRECISION :: calcT

	IF (h<0) THEN
		! SOLIDE
		calcT = TF + h/CS
	ELSEIF (h>=LF) THEN
		! LIQUIDE
		calcT = TF+(h-LF)/CL
	ELSE
		! EQUILIBRE
		calcT = TF
	ENDIF
END FUNCTION

FUNCTION calcY(h)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: h
	DOUBLE PRECISION :: calcY

	IF (h<0) THEN
		! SOLIDE
		calcY = 0
	ELSEIF (h>=LF) THEN
		! LIQUIDE
		calcY = 1
	ELSE
		! EQUILIBRE
		calcY = h/LF
	ENDIF
END FUNCTION

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

  !Température consigne et Temps d'utilisation
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

  !Propriétés du MCP
  READ(FICH,*)epaisseurMCP
  READ(FICH,*)lambdaMCP
	READ(FICH,*)rhoMCP
	READ(FICH,*)CS
	READ(FICH,*)CL
	READ(FICH,*)TF
	READ(FICH,*)LF

	CLOSE(FICH)
END SUBROUTINE load

SUBROUTINE Discretisation
  IMPLICIT NONE

  dx=longueur/Mt
  S1=pi*rayon**2
  S2=pi*((epaisseur2+rayon)**2-rayon**2)
  S3=pi*((epaisseur2+rayon+epaisseurMCP)**2-(rayon+epaisseurMCP)**2)
  Slat1=2*pi*rayon*dx
  Slat2=2*pi*(rayon+epaisseur2)*dx
  Slat4=2*pi*(rayon+epaisseur2+epaisseurMCP)*dx
  Slat3=2*pi*(rayon+epaisseur2+epaisseur3+epaisseurMCP)*dx
  V1=S1*dx
  V2=S2*dx
  V3=S3*dx
END SUBROUTINE Discretisation

SUBROUTINE export()
	IMPLICIT NONE

	OPEN(FICH2,FILE="result_PDC2.txt", ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND");

  ! ------------------------------------------------------------------------------------
  ! Ecriture des données
  ! ------------------------------------------------------------------------------------

    !WRITE(FICH2,*) 'temps=', temps
    !DO i=1,Mt
    !  WRITE(FICH2,*) i*dx, T1(i), T2(i), T3(i)
    !END DO

    WRITE(FICH2,*) temps, T1(Mt), T2(Mt), T3(Mt)


  CLOSE(FICH2)
END SUBROUTINE export

END PROGRAM
