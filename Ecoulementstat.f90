PROGRAM Ecoulement_insta

IMPLICIT NONE

!-----Déclaration des variables------
INTEGER, PARAMETER:: FICH=120
INTEGER, PARAMETER:: FICH2=121
DOUBLE PRECISION, PARAMETER:: Pi=3.14
INTEGER:: i, j, M, Nt, n, nb
DOUBLE PRECISION:: Rayon,T0, debit,dt,duree, dr,dx, lambda1, lambda2, S1, S2, V1, V2, rho, c, To
DOUBLE PRECISION:: F1, F2, F3, F4, F5, F6, F7, h, hf, longueur, epaisseur, tsauvegarde, Tinf, compteur, temps
DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: T1, T2


!-----------------------------------------------------------------------------------------
!						Chargement des parametres depuis le fichier
!-----------------------------------------------------------------------------------------

CALL load()
PRINT*,'Chargement des données.................OK'

dt=duree/Nt
!-----Initialisation des variables--------

ALLOCATE(T1(1:M))  !On stocke les temperatures dans une matrice
ALLOCATE(T2(1:M)) !On stocke les temperatures dans une matrice

T1(:)=T0
T2(:)=T0
nb=0
temps=0
compteur=0

OPEN(FICH2,FILE="result_ecoulement.f90", ACTION="WRITE", STATUS="UNKNOWN");
  WRITE(FICH2,*)
CLOSE(FICH2)
!-----------------------------------------------------------------------------------------
!						       Discrétisation spatiale
!-----------------------------------------------------------------------------------------

CALL Discretisation()
PRINT*,'Discrétisation.........................OK'


CALL Calcul()

!-----------------------------------------------------------------------------------------
!						       Exportation
!-----------------------------------------------------------------------------------------

CALL export()
PRINT*,'Exportation.........................OK'

CONTAINS

SUBROUTINE Calcul()
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------
    !						Discretisation temporelle [Explicite]
    !-----------------------------------------------------------------------------------------

    DO WHILE (temps<duree) !Condition d'arrêt

      DO j=1,M
          IF (j>1) THEN
            F1=S1*lambda1*(T1(j-1)-T1(j))/dx
            F3=S2*lambda2*(T2(j-1)-T2(j))/dx
          ELSE
            F1=S1*lambda1*(T0-T1(1))/dx
            F3=S2*lambda2*(T0-T2(1))/dx
          END IF

          IF (j<M) THEN
            F2=S1*lambda1*(T1(j+1)-T1(j))/dx
            F4=S2*lambda2*(T2(j+1)-T2(j))/dx
          ELSE
            F2=S1*lambda1*(T0-T1(M))/dx
            F4=S2*lambda2*(T0-T2(M))/dx
          END IF

          !Attention conducto-convection
          F5=(Tinf-T2(j))/(epaisseur/2*lambda2+1/h*S1)

          F6=(T2(j)-T1(j))/(epaisseur/2*lambda2+1/hf*S2)

          IF (j>1 .AND. j<M) THEN
            IF (debit>0) THEN
                F7=abs(debit)*c*(T1(j-1)-T1(j))
            ELSE IF (debit<0) THEN
                F7=abs(debit)*c*(T1(j+1)-T1(j))
            ELSE
                F7=0
            END IF
          ELSE IF (j==1) THEN
            IF (debit>0) THEN
                F7=abs(debit)*c*(T0-T1(1))
            ELSE IF (debit<0) THEN
                F7=abs(debit)*c*(T1(2)-T1(1))
            ELSE
                F7=0
            END IF
          ELSE IF (j==M) THEN
            IF (debit>0) THEN
                F7=abs(debit)*c*(T1(M-1)-T1(M))
            ELSE IF (debit<0) THEN
                F7=abs(debit)*c*(T0-T1(M))
            END IF
          END IF

          !Bilan sur l'epaisseur
          T2(j)=dt*(F3+F4+F6+F5)/(rho*c*V2)+T2(j)
          !Bilan sur le tube
          T1(j)=dt*(F1+F2+F7)/(rho*c*V1)+T1(j)
    END DO

!----------------- AUTRE METHODE ---------------------

  !     DO j=2,M-1
  !       F1=S1*lambda1*(T1(j-1)-T1(j))/dx
  !       F2=S1*lambda1*(T1(j+1)-T1(j))/dx
  !       F3=S2*lambda2*(T2(j-1)-T2(j))/dx
  !       F4=S2*lambda2*(T2(j+1)-T2(j))/dx
  !       F5=h*(T2(j)-Tinf)*2*pi*(rayon+epaisseur)*longueur
  !       F6=hf*(T2(j)-T1(j))*2*pi*rayon*longueur
  !     IF (debit>0) THEN
  !       F7=abs(debit)*c*(T1(j-1)-T1(j))
  !     ELSE IF (debit<0) THEN
  !       F7=abs(debit)*c*(T1(j+1)-T1(j))
  !     END IF
  !
  !     !Bilan sur l'epaisseur
  !     T2(j)=dt*(F3+F4+F6)/(rho*c*V2)+T2(j)
  !     !Bilan sur le tube
  !     T1(j)=dt*(F1+F2+F7)/(rho*c*V1)+T1(j)
  !
  !     END DO
  !
  !     !Pour j=1
  !     F1=S1*lambda1*(T0-T1(1))/dx
  !     F2=S1*lambda1*(T1(2)-T1(1))/dx
  !     F3=S2*lambda2*(T0-T2(1))/dx
  !     F4=S2*lambda2*(T2(2)-T2(1))/dx
  !     F5=h*(T2(1)-Tinf)*2*pi*(rayon+epaisseur)*longueur
  !     F6=hf*(T2(1)-T1(1))*2*pi*rayon*longueur
  !   IF (debit>0) THEN
  !     F7=abs(debit)*c*(T0-T1(1))
  !   ELSE IF (debit<0) THEN
  !     F7=abs(debit)*c*(T1(2)-T1(1))
  !   END IF
  !
  !   !Bilan en j=1 sur l'épaisseur
  !   T2(1)=dt*(F3+F4+F6)/(rho*c*V2)+T2(1)
  !   !Bilan sur le tube
  !   T1(1)=dt*(F1+F2+F7)/(rho*c*V1)+T1(1)
  !
  !   !Pour j=M
  !   F1=S1*lambda1*(T1(M-1)-T1(M))/dx
  !   F2=S1*lambda1*(T0-T1(M))/dx
  !   F3=S2*lambda2*(T2(M-1)-T2(M))/dx
  !   F4=S2*lambda2*(T0-T2(M))/dx
  !   F5=h*(T2(M)-Tinf)*2*pi*(rayon+epaisseur)*longueur
  !   F6=hf*(T2(M)-T1(M))*2*pi*rayon*longueur
  ! IF (debit>0) THEN
  !   F7=abs(debit)*c*(T1(M-1)-T1(M))
  ! ELSE IF (debit<0) THEN
  !   F7=abs(debit)*c*(T0-T1(M))
  ! ELSE
  !   F7=0
  ! END IF
  !
  ! !Bilan en j=M
  ! T2(M)=dt*(F3+F4+F6)/(rho*c*V2)+T2(M)
  ! !Bilan sur le tube
  ! T1(M)=dt*(F1+F2+F7)/(rho*c*V1)+T1(M)

!-------------------------------------------------
  temps=temps+dt
  compteur=compteur+dt

  !-------------------------------------------------------------------------------------
  !						Ecriture de certain résultats dans un fichier
  !-------------------------------------------------------------------------------------
  IF (compteur>tsauvegarde-dt .AND. compteur<=tsauvegarde+dt) THEN
    nb=nb+1 ! nombre de données sauvegardées
    CALL export()
    compteur=0
  END IF

END DO

END SUBROUTINE

! Chargement des données depuis le fichier donnees_mur.f90
SUBROUTINE load()
	IMPLICIT NONE

	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	OPEN(FICH,FILE="donnees_ecoulement.f90", ACTION="READ");

	! ------------------------------------------------------------------------------------
	! Lecture des données globales
	! ------------------------------------------------------------------------------------

	READ(FICH,*) M
	READ(FICH,*) Rayon
  READ(FICH,*) debit
	READ(FICH,*) rho
  READ(FICH,*) C
  READ(FICH,*) lambda1
  READ(FICH,*) lambda2
  READ(FICH,*) To
  READ(FICH,*) Tinf
  READ(FICH,*) h
  READ(FICH,*) hf
  READ(FICH,*) longueur
  READ(FICH,*) epaisseur
  READ(FICH,*) Nt
  READ(FICH,*) duree
  READ(FICH,*) tsauvegarde

	CLOSE(FICH)
END SUBROUTINE load

SUBROUTINE Discretisation
  IMPLICIT NONE

  dx=longueur/M
  dr=Rayon/M

  S1=pi*rayon**2
  S2=pi*((epaisseur+rayon)**2-rayon**2)
  V1=S1*dx
  V2=S2*dx

END SUBROUTINE Discretisation

! Exportation des données dans le fichier result.f90
SUBROUTINE export()
	IMPLICIT NONE


	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	OPEN(FICH2,FILE="result_ecoulement.f90", ACTION="WRITE", STATUS="UNKNOWN", position="append");


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
