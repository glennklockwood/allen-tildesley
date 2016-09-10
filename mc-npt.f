********************************************************************************
** FICHE F.12.  CONSTANT-NPT MONTE CARLO ALGORITHM.                           **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

        PROGRAM MCNPT

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** MONTE CARLO SIMULATION IN THE CONSTANT-NPT ENSEMBLE.          **
C    **                                                               **
C    ** THIS PROGRAM TAKES A CONFIGURATION OF LENNARD JONES ATOMS     **
C    ** AND PERFORMS A MONTE CARLO SIMULATION AT CONSTANT NPT. THE    **
C    ** BOX IS IN UNITS OF SIGMA THE LENNARD JONES DIAMETER.          **
C    ** THERE ARE NO LOOKUP TABLES INCLUDED.                          **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** MCDONALD, CHEM. PHYS. LETT. 3, 241, 1969.                     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                             NUMBER OF MOLECULES     **
C    ** REAL    RX(N),RY(N),RZ(N)             POSITIONS               **
C    ** REAL    VOL                           VOLUME                  **
C    ** REAL    BOX                           BOX LENGTH              **
C    ** REAL    DENS                          REDUCED DENSITY         **
C    ** REAL    TEMP                          REDUCED TEMPERATURE     **
C    ** REAL    SIGMA                         REDUCED LJ DIAMETER     **
C    ** REAL    DRMAX                         MAXIMUM DISPLACEMENT    **
C    ** REAL    V                             THE POTENTIAL ENERGY    **
C    ** REAL    W                             THE VIRIAL              **
C    ** REAL    PRESUR                        REQUIRED PRESSURE       **
C    ** REAL    DBOXMX                        MAX CHANGE IN BOX       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** CONDUCTS MONTE CARLO SIMULATION AT CONSTANT PRESSURE FOR A    **
C    ** SPECIFIED NUMBER OF CYCLES FROM A GIVEN INITIAL CONFIGURATION.**
C    **                                                               **
C    ** UNITS:                                                        **
C    **                                                               **
C    ** THIS PROGRAM USES THE USUAL REDUCED LJ UNITS. IN PARTICULAR   **
C    ** THE BOX LENGTH IS IN UNITS OF SIGMA.                          **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE ENERGY ( RXI, RYI, RZI, I, RCUT, BOX, V12, V6,     **
C    **    :                W12, W6 )                                 **
C    **    CALCULATES THE ENERGY AND VIRIAL FOR ATOM I IN THE FLUID   **
C    ** REAL FUNCTION RANF ( DUMMY )                                  **
C    **    RETURNS A UNIFORM RANDOM NUMBER BETWEEN ZERO AND ONE       **
C    ** SUBROUTINE READCN ( CNFILE, BOX )                             **
C    **    READS IN CONFIGURATION AND BOX VARIABLES                   **
C    ** SUBROUTINE SUMUP ( RCUT, RMIN, OVRLAP, BOX, V12, V6, W12, W6 )**
C    **    CALCULATES POTENTIAL AND VIRIAL FOR A CONFIGURATION        **
C    ** SUBROUTINE WRITCN ( CNFILE, BOX )                             **
C    **    WRITES OUT CONFIGURATION AND BOX VARIABLES                 **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        REAL        RX(N), RY(N), RZ(N)

        INTEGER     STEP, NSTEP, IPRINT, IRATIO, IRATB, I
        REAL        ACV, ACP, ACD, ACM, ACATMA, ACBOXA, NORM
        REAL        ACVSQ, ACPSQ, ACDSQ
        REAL        AVV, AVP, AVD
        REAL        FLV, FLP, FLD
        REAL        DENS, TEMP, RCUT, RMIN, PRESUR, BOX, VOL, PRES, VN
        REAL        BOXINV, BOXNEW, RATBOX, RAT12, RAT6, DVOL, DPV
        REAL        DELTHB, DRMAX, DBOXMX, BETA, RANF, DUMMY, RATIO
        REAL        RRBOX, BRATIO, DELTVB, RCUTN
        REAL        RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW
        REAL        V12OLD, V6OLD, V12NEW, V6NEW, VS
        REAL        W12OLD, W6OLD, W12NEW, W6NEW, WS, PS
        REAL        DELV12, DELV6, DELW12, DELW6, DELTV
        REAL        V6, W6, V12, W12, VLRC, VLRCN, WLRC, WLRCN
        REAL        SR3, SR9, VLRC6, WLRC6, VLRC12, WLRC12, PI
        CHARACTER   TITLE*80, CNFILE*30
        LOGICAL     OVRLAP

        PARAMETER ( PI = 3.1415927 )

C    *******************************************************************

        WRITE(*,'(1H1,'' **** PROGRAM MCNPT ****               '')')
        WRITE(*,'(//'' MONTE CARLO IN A CONSTANT-NPT ENSEMBLE  '')')
        WRITE(*,'(  '' LENNARD-JONES ATOMS                     '')')

C    ** BASIC SIMULATION PARAMETERS **

        WRITE(*,'('' ENTER RUN TITLE                           '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER CONFIGURATION FILENAME              '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' ENTER NUMBER OF CYCLES                    '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER INTERVAL BETWEEN PRINTS IN CYCLES   '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER INTERVAL FOR UPDATE OF MAXIMUM      '')')
        WRITE(*,'('' DISPLACEMENT OF ATOMS IN CYCLES           '')')
        READ (*,*) IRATIO
        WRITE(*,'('' ENTER INTERVAL FOR UPDATE OF MAXIMUM      '')')
        WRITE(*,'('' DISPLACEMENT OF THE BOX IN CYCLES         '')')
        READ (*,*) IRATB
        WRITE(*,'(/'' ENTER THE FOLLOWING IN L-J REDUCED UNITS  ''/)')
        WRITE(*,'('' ENTER POTENTIAL CUTOFF                    '')')
        READ (*,*) RCUT
        WRITE(*,'('' ENTER DESIRED PRESSURE                    '')')
        READ (*,*) PRESUR
        WRITE(*,'('' ENTER DESIRED TEMPERATURE                 '')')
        READ (*,*) TEMP

        WRITE(*,'(//1X,A)') TITLE
        WRITE(*,'('' CONFIGURATION FILENAME            '',A)') CNFILE
        WRITE(*,'('' NUMBER OF CYCLES                = '',I6)') NSTEP
        WRITE(*,'('' PRINT INTERVAL                  = '',I6)') IPRINT
        WRITE(*,'('' RATIO UPDATE INTERVAL FOR ATOMS = '',I6)') IRATIO
        WRITE(*,'('' RATIO UPDATE INTERVAL FOR BOX   = '',I6)') IRATB
        WRITE(*,'('' POTENTIAL CUTOFF                = '',F10.5)')RCUT
        WRITE(*,'('' DESIRED PRES.                   = '',F10.5)')PRESUR
        WRITE(*,'('' DESIRED TEMP                    = '',F10.5)')TEMP

C    ** READCN READS IN INITIAL CONFIGURATION **
C    ** ORIGIN SHOULD BE AT CENTRE OF BOX     **

        CALL READCN ( CNFILE, BOX )

C    ** SET DEPENDENT PARAMETERS **

        VOL    = BOX ** 3
        BOXINV = 1.0 / BOX
        DENS   = REAL ( N ) / VOL

        IF ( RCUT .GT. ( 0.5 * BOX ) ) STOP 'CUT-OFF TOO LARGE'

        DBOXMX = BOX / 40.0
        DRMAX  = 0.15
        RMIN   = 0.70
        BETA   = 1.0 / TEMP

        WRITE(*,'('' INITIAL DENSITY                 = '',F10.5)') DENS

C    ** CALCULATE LONG-RANGE CORRECTIONS FOR LJ POTENTIAL.    **
C    ** 6 IS FOR ATTRACTIVE CONTRIBUTIONS 12 IS FOR REPULSIVE **

        SR3    =    ( 1.0 / RCUT ) ** 3
        SR9    =    SR3 ** 3
        VLRC12 =    8.0 * PI * DENS * REAL ( N ) * SR9 / 9.0
        VLRC6  = -  8.0 * PI * DENS * REAL ( N ) * SR3 / 3.0
        WLRC12 =    4.0 * VLRC12
        WLRC6  =    2.0 * VLRC6
        VLRC   =    VLRC12 + VLRC6
        WLRC   =    WLRC12 + WLRC6

C    ** ZERO ACCUMULATORS **

        ACM    = 0.0
        ACATMA = 0.0
        ACBOXA = 0.0

        ACV = 0.0
        ACP = 0.0
        ACD = 0.0

        ACVSQ = 0.0
        ACPSQ = 0.0
        ACDSQ = 0.0

        FLV = 0.0
        FLP = 0.0
        FLD = 0.0

C    ** CALCULATE INITIAL ENERGY AND VIRIAL **

        CALL SUMUP ( RCUT, RMIN, OVRLAP, BOX, V12, V6, W12, W6 )

        IF ( OVRLAP ) STOP 'OVERLAP IN INITIAL CONFIGURATION'

C    ** CALCULATE THE INITIAL ENERGY AND VIRIAL **

        VS = ( V12 + V6 + VLRC ) / REAL ( N )
        WS = ( W12 + W6 + WLRC ) / REAL ( N )
        PS = DENS * TEMP + ( W12 + W6 + WLRC ) / VOL

C    ** ADD LONG RANGE CORRECTIONS **
C    ** INTO THE ENERGY AND VIRIAL **

        V12 = V12 + VLRC12
        V6  = V6  + VLRC6
        W12 = W12 + WLRC12
        W6  = W6  + WLRC6

        WRITE(*,'('' INITIAL V/N                     = '',F10.6)') VS
        WRITE(*,'('' INITIAL W/N                     = '',F10.6)') WS
        WRITE(*,'('' INITIAL P                       = '',F10.6)') PS
        WRITE(*,'(//'' **** START OF MARKOV CHAIN ****'')')
        WRITE(*,10001)

C    *******************************************************************
C    ** MAIN LOOP STARTS                                              **
C    *******************************************************************

        DO 100 STEP = 1, NSTEP

C       ** LOOP OVER ATOMS STARTS **

           DO 97 I = 1, N

              RXIOLD = RX(I)
              RYIOLD = RY(I)
              RZIOLD = RZ(I)

C          ** CALCULATE V FOR AN ATOM IN OLD STATE **

              CALL ENERGY ( RXIOLD, RYIOLD, RZIOLD, I, RCUT, BOX,
     :                      V12OLD, V6OLD, W12OLD, W6OLD )

C          ** MOVE ATOM I AND PICKUP CENTRAL IMAGE **

              RXINEW = RXIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX
              RYINEW = RYIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX
              RZINEW = RZIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX

              RXINEW = RXINEW - ANINT ( RXINEW * BOXINV ) * BOX
              RYINEW = RYINEW - ANINT ( RYINEW * BOXINV ) * BOX
              RZINEW = RZINEW - ANINT ( RZINEW * BOXINV ) * BOX

C          ** CALCULATE V FOR ATOM IN NEW STATE **

              CALL ENERGY ( RXINEW, RYINEW, RZINEW, I, RCUT, BOX,
     :                      V12NEW, V6NEW, W12NEW, W6NEW )

C          ** CHECK FOR ACCEPTANCE **

              DELV12 = V12NEW - V12OLD
              DELV6  = V6NEW  - V6OLD
              DELW12 = W12NEW - W12OLD
              DELW6  = W6NEW  - W6OLD
              DELTV  = DELV12 + DELV6
              DELTVB = BETA * DELTV

              IF ( DELTVB .LT. 75.0 ) THEN

                 IF ( DELTV .LE. 0.0 ) THEN

                    V12    = V12 + DELV12
                    V6     = V6  + DELV6
                    W12    = W12 + DELW12
                    W6     = W6  + DELW6
                    RX(I)  = RXINEW
                    RY(I)  = RYINEW
                    RZ(I)  = RZINEW
                    ACATMA = ACATMA + 1.0

                 ELSEIF ( EXP ( - DELTVB ) .GT. RANF ( DUMMY ) ) THEN

                    V12    = V12 + DELV12
                    V6     = V6  + DELV6
                    W12    = W12 + DELW12
                    W6     = W6  + DELW6
                    RX(I)  = RXINEW
                    RY(I)  = RYINEW
                    RZ(I)  = RZINEW
                    ACATMA = ACATMA + 1.0

                 ENDIF

              ENDIF

              VN   = ( V12 + V6 ) / REAL ( N )
              PRES = DENS * TEMP +  ( W12 + W6 )  / VOL

C       ** INCREMENT ACCUMULATORS **

              ACM = ACM + 1.0
              ACV = ACV + VN
              ACP = ACP + PRES
              ACD = ACD + DENS

              ACVSQ = ACVSQ + VN ** 2
              ACPSQ = ACPSQ + PRES ** 2
              ACDSQ = ACDSQ + DENS ** 2

97         CONTINUE

C       ** ENDS LOOP OVER ATOMS IN ONE CYCLE  **

C       ** ATTEMPT A BOX MOVE **

           BOXNEW = BOX + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DBOXMX
           RATBOX = BOX / BOXNEW
           RRBOX  = 1.0 / RATBOX
           RCUTN  = RCUT * RRBOX

C       ** CALCULATE SCALING PARAMETERS **

           RAT6   = RATBOX ** 6
           RAT12  = RAT6 * RAT6

C       ** SCALE ENERGY, AND VIRIAL INCLUDING LRC **

           V12NEW = V12  * RAT12
           V6NEW  = V6   * RAT6
           W12NEW = W12  * RAT12
           W6NEW  = W6   * RAT6

C       ** CALCULATE CHANGE IN ENERGY AND VOLUME **

           DELTV  = V12NEW + V6NEW - V12 - V6
           DPV    = PRESUR * ( BOXNEW ** 3 - VOL )
           DVOL   = 3.0 * TEMP * REAL ( N ) * ALOG ( RATBOX )
           DELTHB = BETA * ( DELTV + DPV + DVOL )

C       ** CHECK FOR ACCEPTANCE **

           IF ( DELTHB .LT. 75.0 ) THEN

              IF ( DELTHB .LE. 0.0 ) THEN

                 V12    = V12NEW
                 V6     = V6NEW
                 W12    = W12NEW
                 W6     = W6NEW

                 DO 98 I = 1, N

                    RX(I) = RX(I) * RRBOX
                    RY(I) = RY(I) * RRBOX
                    RZ(I) = RZ(I) * RRBOX

98               CONTINUE

                 BOX    = BOXNEW
                 RCUT   = RCUTN
                 ACBOXA = ACBOXA + 1.0

              ELSEIF ( EXP ( - DELTHB ) .GT. RANF ( DUMMY ) ) THEN

                 V12    = V12NEW
                 V6     = V6NEW
                 W12    = W12NEW
                 W6     = W6NEW

                 DO 99 I = 1, N

                    RX(I) = RX(I) * RRBOX
                    RY(I) = RY(I) * RRBOX
                    RZ(I) = RZ(I) * RRBOX

99               CONTINUE

                 BOX    = BOXNEW
                 RCUTN  = RCUT
                 ACBOXA = ACBOXA + 1.0

              ENDIF

           ENDIF

           BOXINV = 1.0 / BOX
           VOL    = BOX ** 3
           DENS   = REAL ( N ) / VOL

C       ** CALCULATE ENERGY AND PRESSURE **

           VN   = ( V12 + V6 ) / REAL ( N )
           PRES = DENS * TEMP + ( W12 + W6 ) / VOL

C       ** INCREMENT ACCUMULATORS **

           ACM = ACM + 1.0
           ACV = ACV + VN
           ACP = ACP + PRES
           ACD = ACD + DENS

           ACVSQ = ACVSQ + VN ** 2
           ACPSQ = ACPSQ + PRES ** 2
           ACDSQ = ACDSQ + DENS ** 2

C       ** ENDS ATTEMPTED BOX MOVE **

C       ** PERFORM PERIODIC OPERATIONS **

           IF ( MOD ( STEP, IRATIO ) .EQ. 0 ) THEN

C          ** ADJUST MAXIMUM DISPLACEMENT FOR ATOMS **

              RATIO = ACATMA / REAL ( N * IRATIO )

              IF ( RATIO .GT. 0.5 ) THEN

                 DRMAX = DRMAX * 1.05

              ELSE

                 DRMAX = DRMAX * 0.95

              ENDIF

              ACATMA = 0.0

           ENDIF

           IF ( MOD ( STEP, IRATB ) .EQ. 0 ) THEN

C          ** ADJUST MAXIMUM DISPLACEMENT FOR THE BOX **

              BRATIO = ACBOXA / REAL ( IRATB )

              IF ( BRATIO .GT. 0.5 ) THEN

                 DBOXMX = DBOXMX * 1.05

              ELSE

                 DBOXMX = DBOXMX * 0.95

              ENDIF

              ACBOXA = 0.0

           ENDIF

           IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN

C          ** OPTIONALLY PRINT INFORMATION **

              WRITE(*,'(1X,I8,5(2X,F10.5))')
     :                 STEP, VN, PRES, DENS, RATIO, BRATIO

           ENDIF

100     CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(/1X,''**** END OF MARKOV CHAIN **** ''//)')

C    ** WRITE OUT FINAL CONFIGURATION AND BOXLENGTH **

        CALL WRITCN ( CNFILE, BOX )

C    ** WRITE OUT FINAL AVERAGES **

        NORM = REAL ( ACM )
        AVV  = ACV / NORM
        AVP  = ACP / NORM
        AVD  = ACD / NORM

        ACVSQ = ( ACVSQ / NORM ) - AVV ** 2
        ACPSQ = ( ACPSQ / NORM ) - AVP ** 2
        ACDSQ = ( ACDSQ / NORM ) - AVD ** 2

        IF ( ACVSQ .GT. 0.0 ) FLV = SQRT ( ACVSQ )
        IF ( ACPSQ .GT. 0.0 ) FLP = SQRT ( ACPSQ )
        IF ( ACDSQ .GT. 0.0 ) FLD = SQRT ( ACDSQ )

        WRITE(*, 10002)
        WRITE(*,'('' AVERAGES'',3(2X,F10.5))') AVV, AVP, AVD
        WRITE(*,'('' FLUCTS  '',3(2X,F10.5))') FLV, FLP, FLD

        STOP

10001   FORMAT(//1X,'  CYCLE   ..POTENT..  ..PRESSURE.. ..DENSITY..',
     :              ' ..RATIO..  ..BRATIO..')
10002   FORMAT(//1X,'  CYCLE   ..POTENT..  ..PRESSURE.. ..DENSITY..')

        END



        SUBROUTINE SUMUP ( RCUT, RMIN, OVRLAP, BOX, V12, V6, W12, W6 )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** CALCULATES THE TOTAL POTENTIAL ENERGY FOR A CONFIGURATION.    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS                 **
C    ** REAL    RX(N(,RY(N),RZ(N) THE POSITIONS OF THE ATOMS          **
C    ** REAL    V                 THE POTENTIAL ENERGY                **
C    ** REAL    W                 THE VIRIAL                          **
C    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE SUBROUTINE RETURNS THE TOTAL POTENTIAL ENERGY AT THE      **
C    ** BEGINNING AND END OF THE RUN.                                 **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        RMIN, RCUT, V12, V6, W12, W6, BOX
        LOGICAL     OVRLAP

        REAL        RCUTSQ, RMINSQ, VIJ12, VIJ6
        REAL        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ
        REAL        SR2, SR6, RIJSQ, BOXINV
        INTEGER     I, J

C    *******************************************************************

        OVRLAP = .FALSE.
        RCUTSQ = RCUT * RCUT
        RMINSQ = RMIN * RMIN
        BOXINV = 1.0 / BOX

        V12    = 0.0
        V6     = 0.0
        W12    = 0.0
        W6     = 0.0

C    ** LOOP OVER ALL THE PAIRS IN THE LIQUID **

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)

           DO 99 J = I + 1, N

              RXIJ  = RXI - RX(J)
              RYIJ  = RYI - RY(J)
              RZIJ  = RZI - RZ(J)

              RXIJ  = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
              RYIJ  = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
              RZIJ  = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

              IF ( RIJSQ .LT. RMINSQ ) THEN

                 OVRLAP = .TRUE.
                 RETURN

              ELSEIF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2   = 1.0 / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ12 = SR6 *  SR6
                 VIJ6  = - SR6
                 V12   = V12 + VIJ12
                 V6    = V6  + VIJ6
                 W12   = W12 + VIJ12
                 W6    = W6  + VIJ6 * 0.5

              ENDIF

99         CONTINUE

100     CONTINUE

        V12 = 4.0 * V12
        V6  = 4.0 * V6
        W12 = 48.0 * W12 / 3.0
        W6  = 48.0 * W6  / 3.0

        RETURN
        END



        SUBROUTINE ENERGY ( RXI, RYI, RZI, I, RCUT, BOX,
     :                      V12, V6, W12, W6 )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** CALCULATES THE POTENTIAL ENERGY OF I WITH ALL OTHER ATOMS     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER I                 THE ATOM OF INTEREST                **
C    ** INTEGER N                 THE NUMBER OF ATOMS                 **
C    ** REAL    RX(N),RY(N),RZ(N) THE ATOM POSITIONS                  **
C    ** REAL    RXI,RYI,RZI       THE COORDINATES OF ATOM I           **
C    ** REAL    V                 THE POTENTIAL ENERGY OF ATOM I      **
C    ** REAL    W                 THE VIRIAL OF ATOM I                **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL MOVE OF ATOM I. IT IS CALLED BEFORE AND        **
C    ** AFTER THE RANDOM DISPLACEMENT OF I.                           **
C    *******************************************************************


        INTEGER     N
        PARAMETER ( N = 108 )
        REAL        RX(N), RY(N), RZ(N)
        REAL        RCUT, BOX, RXI, RYI, RZI, V12, V6, W12, W6
        INTEGER     I

        REAL        RCUTSQ, VIJ12, VIJ6
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, SR2, SR6, BOXINV
        INTEGER     J

C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        BOXINV = 1.0 / BOX

        V12    = 0.0
        V6     = 0.0
        W12    = 0.0
        W6     = 0.0

C    ** LOOP OVER ALL MOLECULES EXCEPT I  **

        DO 100 J = 1, N

        IF ( I .NE. J ) THEN

           RXIJ  = RXI - RX(J)
           RYIJ  = RYI - RY(J)
           RZIJ  = RZI - RZ(J)

           RXIJ  = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
           RYIJ  = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
           RZIJ  = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX

           RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

           IF ( RIJSQ .LT. RCUTSQ ) THEN

              SR2   = 1.0 / RIJSQ
              SR6   = SR2 * SR2 * SR2
              VIJ12 = SR6 * SR6
              VIJ6  = - SR6
              V12   = V12 + VIJ12
              V6    = V6  + VIJ6
              W12   = W12 + VIJ12
              W6    = W6  + VIJ6 * 0.5

           ENDIF

        ENDIF

100     CONTINUE

        V12 = 4.0 * V12
        V6  = 4.0 * V6
        W12 = 48.0 * W12 / 3.0
        W6  = 48.0 * W6  / 3.0

        RETURN
        END



        SUBROUTINE READCN ( CNFILE, BOX )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN THE CONFIGURATION FROM UNIT 10          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N)
        REAL        BOX

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )
        INTEGER     NN

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD',
     :         FORM = 'UNFORMATTED'                        )

        READ ( CNUNIT ) NN, BOX
        IF ( NN .NE. N ) STOP 'N ERROR IN READCN'
        READ ( CNUNIT ) RX, RY, RZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE, BOX )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 10          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N)
        REAL        BOX

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
     :         FORM = 'UNFORMATTED'                        )

        WRITE ( CNUNIT ) N, BOX
        WRITE ( CNUNIT ) RX, RY, RZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        REAL FUNCTION RANF ( DUMMY )

C    *******************************************************************
C    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
C    **                                                               **
C    **                 ***************                               **
C    **                 **  WARNING  **                               **
C    **                 ***************                               **
C    **                                                               **
C    ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
C    ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
C    *******************************************************************

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
        REAL        DUMMY
        SAVE        SEED
        DATA        SEED / 0 /

C    *******************************************************************

        SEED = MOD ( SEED * L + C, M )
        RANF = REAL ( SEED ) / M

        RETURN
        END



