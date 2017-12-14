PROGRAM Ecoulement_instaMCP

IMPLICIT NONE

!-----------------------------------------------------------------------------------------
!						       Déclaration des variables
!-----------------------------------------------------------------------------------------
INTEGER, PARAMETER:: FICH=127
INTEGER, PARAMETER:: FICH2=150
DOUBLE PRECISION, PARAMETER:: Pi=4*ATAN(1.D0)
INTEGER:: i, j, m, Mt, Nt, n, nb, compteur2, compteur4
DOUBLE PRECISION:: longueur,Rayon, epaisseur2, epaisseur3, epaisseurMCP, tsauvegarde, compteur, temps, dt, duree,dx
DOUBLE PRECISION:: Tamb, rho,rhop, Ceau, Cper, Tprep, hair, heau, vitesse, attente
DOUBLE PRECISION:: rhoMCP, CS, CL,TF, LF, HG , moitie
DOUBLE PRECISION:: lambda1, lambda2, lambda3, lambdaMCP, S1, S2, S3, V1, V2,V3, Slat1, Slat2, Slat3,Slat4
DOUBLE PRECISION:: F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, compteur3, tconsigne,tu, tempsf, trefroi
DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: T1, T2, T3, h, Y

!-----------------------------------------------------------------------------------------
!						       Importation des données du problème
!-----------------------------------------------------------------------------------------
CALL load()
PRINT*,'Chargement des données.................OK'
!-----------------------------------------------------------------------------------------
!						       Initialisation des variables
!-----------------------------------------------------------------------------------------
moitie=Mt/2
ALLOCATE(T1(1:Mt)) !On stocke les temperatures dans une matrice
ALLOCATE(T2(1:Mt)) !On stocke les temperatures dans une matrice
ALLOCATE(T3(nint(moitie):Mt))
ALLOCATE(h(nint(moitie):Mt))
ALLOCATE(Y(nint(moitie):Mt))
!Etat refroidi: toute la conduite est à Tambiant (K)
T1(:)=Tamb
T2(:)=Tamb
T3(:)=Tamb
dt=duree/Nt
nb=0
temps=0
compteur=0
compteur2=0
compteur3=0
compteur4=0
attente=0  !Par défaut indique que la température seuil n'a pas été atteinte

OPEN(FICH2,FILE="result_PDC2_aval.txt", ACTION="WRITE", STATUS="UNKNOWN");
  WRITE(FICH2,*)
CLOSE(FICH2)

OPEN(FICH2+1,FILE="resultsortie_PDC2.txt", ACTION="WRITE", STATUS="UNKNOWN");
  WRITE(FICH2+1,*)
CLOSE(FICH2+1)
!-----------------------------------------------------------------------------------------
!						               Discrétisation spatiale
!-----------------------------------------------------------------------------------------

CALL Discretisation()
PRINT*,'Discrétisation.........................OK'

!Contrainte sur le pas de temps
PRINT*, "Il faut:", vitesse*S1*dt, "<<", V3
PRINT*, "Sinon l'hypothèse que la température dans l'écoulement reste constante"
PRINT*, "en un temps dt n'est pas vérifiée."
!-----------------------------------------------------------------------------------------
!						                   Calcul
!-----------------------------------------------------------------------------------------
CALL Calcul()
PRINT*,'Calcul.................................OK'
!Temps d'attente pour atteindre la valeur consigne
PRINT*, "Pour obtenir une température en sortie de", Tconsigne,"°C on devra attendre", attente  ,"s"
PRINT*, "(Si 0s est affiché c'est que la température n'est pas atteinte durant ce temps de simulation)."
PRINT*, "L'ECS est refroidi à ", Trefroi,"°C après un temps de refroidissement de ", Tempsf, "s soit ", Tempsf/60, "minutes"
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
      !Sans MCP
      DO m=1,nint(moitie)

        Slat3=2*pi*(rayon+epaisseur2+epaisseur3)*dx

        IF (m>1) THEN
          F1=S1*lambda1*(T1(m-1)-T1(m))/dx
          F3=S2*lambda2*(T2(m-1)-T2(m))/dx
        ELSE
          F1=2*S1*lambda1*(Tprep-T1(1))/dx
          F3=2*S2*lambda2*(Tprep-T2(1))/dx
        END IF

        F2=S1*lambda1*(T1(m+1)-T1(m))/dx
        F4=S2*lambda2*(T2(m+1)-T2(m))/dx

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

      !Avec MCP
      DO m=nint(moitie)+1,Mt
        Slat3=2*pi*(rayon+epaisseur2+epaisseur3+epaisseurMCP)*dx

        F1=S1*lambda1*(T1(m-1)-T1(m))/dx
        F3=S2*lambda2*(T2(m-1)-T2(m))/dx
        F9=S3*lambdaMCP*(T3(m-1)-T3(m))/dx

        IF (m<Mt) THEN
          F2=S1*lambda1*(T1(m+1)-T1(m))/dx
          F4=S2*lambda2*(T2(m+1)-T2(m))/dx
          F10=S3*lambdaMCP*(T2(m+1)-T3(m))/dx
        ELSE
          F2=2*S1*lambda1*(Tamb-T1(m))/dx
          F4=2*S2*lambda2*(Tamb-T2(m))/dx
          F10=2*S3*lambdaMCP*(Tamb-T3(m))/dx
        END IF

        !Flux conducto-convectifs
        F6=(T1(m)-T2(m))/(epaisseur2/(2*lambda2*Slat1)+rayon/(rayon*heau*Slat1+lambda1*Slat1))

        !Flux avec la résistance cylindrique correspondant à l'isolant
        F7=log((rayon+epaisseur2+epaisseur3+epaisseurMCP)/(rayon+epaisseur2+epaisseurMCP))/(2*pi*lambdaMCP*dx)
        F7=(T3(m)-Tamb)/(1/(hair*Slat3)+F7)

        F5=vitesse*rho*S1*Ceau*(T1(m-1)-T1(m))

        F8=Slat2*(T2(m)-T3(m))/(epaisseur2/(2*lambda2)+epaisseurMCP/(2*lambdaMCP))

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
      !Condition pour récupérer le temps d'attente pour atteindre Tconsigne
      IF (T1(Mt)>Tconsigne .AND. compteur2==0) THEN
          attente=temps
          compteur2=1
      END IF

      !Condition pour compter le temps d'utilisation
      IF (T1(Mt)>Tconsigne) THEN
          compteur3=compteur3+dt
      END IF

      !Condition pour couper l'ecoulement à Tutilisation
      IF (compteur3>tu) THEN
        vitesse=0
      END IF

      !Condition pour connaitre le temps de refroidissement
      IF(T1(Mt)<Trefroi .AND. compteur4==0 .AND. compteur2==1) THEN
        compteur4=1
        tempsf=temps-attente-tu
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
        !PRINT*, "F8=",F8,"F9=",F9,"F10=",F10
      END IF

    END DO
    !--------------FIN DE LA BOUCLE DE TEMPS---------------------
END SUBROUTINE

FUNCTION calcH(T)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: T
	DOUBLE PRECISION :: calcH

  !Dans le cas de la paraffine on définit un cp moyen entre 28-32
  IF (T<32 .AND. T>28) THEN
    CS=20000
    CL=20000
  ELSE
    CS=3000
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

  !Température consigne Temps d'utilisation
  READ(FICH,*)Tconsigne
  READ(FICH,*)Tu
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

	OPEN(FICH2,FILE="result_PDC2_aval.txt", ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND");

  ! ------------------------------------------------------------------------------------
  ! Ecriture des données pour obtenir T(x)
  ! ------------------------------------------------------------------------------------

    WRITE(FICH2,*) 'temps=', temps
    DO i=1,nint(moitie)
      WRITE(FICH2,*) i*dx, T1(i), T2(i)
    END DO
    DO i=nint(moitie),Mt
      WRITE(FICH2,*) i*dx, T1(i), T2(i), T3(i)
    END DO
  CLOSE(FICH2)

  OPEN(FICH2+1,FILE="resultsortie_PDC2.txt", ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND");

  ! ------------------------------------------------------------------------------------
  ! Ecriture des données pour obtenir T(t)
  ! ------------------------------------------------------------------------------------
      WRITE(FICH2+1,*) temps, T1(Mt), T2(Mt), T3(Mt)

  CLOSE(FICH2+1)
END SUBROUTINE export

END PROGRAM
