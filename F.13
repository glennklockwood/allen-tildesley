********************************************************************************
** FICHE F.13.  THE HEART OF A CONSTANT MU VT MONTE CARLO PROGRAM             **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** ATTEMPTED CREATIONS AND DESTRUCTIONS IN GRAND CANONICAL MC.   **
C    **                                                               **
C    ** THESE ROUTINES ALLOW FOR A TRIAL DESTRUCTION OR CREATION IN A **
C    ** GRAND CANONICAL MONTE CARLO PROGRAM.                          **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF ATOMS BEFORE THE TRIAL  **
C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING THE TRIAL  **
C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS ALLOWED   **
C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF ATOM     **
C    ** REAL    RX(NMAX) ETC.       POSITIONS OF CURRENT ATOMS        **
C    ** REAL    V                   POTENTIAL ENERGY + LRC            **
C    ** REAL    W                   VIRIAL + LRC                      **
C    ** REAL    DELTV               CHANGE IN ENERGY                  **
C    ** REAL    DELTW               CHANGE IN VIRIAL                  **
C    ** REAL    TEMP                REDUCED TEMPERATURE               **
C    ** REAL    Z                   ABSOLUTE ACTIVITY COEFFICIENT     **
C    ** REAL    SIGMA               LENNARD JONES DIAMETER            **
C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
C    ** REAL    RMIN                REDUCED MINIMUM SEPARATION        **
C    ** LOGICAL OVRLAP              TRUE FOR SUBSTANTIAL ATOM OVERLAP **
C    ** LOGICAL CREATE              TRUE FOR AN ACCEPTED CREATION     **
C    ** LOGICAL GHOST               TRUE FOR AN ACCEPTED DESTRUCTION  **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE IN ( TEMP, Z, SIGMA, RCUT, N, V, W, CREATE )       **
C    **    PERFORMS A TRIAL CREATION                                  **
C    ** SUBROUTINE OUT ( TEMP, Z, SIGMA, RCUT, N, V, W, GHOST )       **
C    **    PERFORMS A TRIAL DESTRUCTION                               **
C    ** SUBROUTINE POTIN ( RXNEW, RYNEW, RZNEW, N, SIGMA, RCUT, RMIN, **
C    ** :                  DELTV, DELTW, OVRLAP )                     **
C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON CREATION         **
C    ** SUBROUTINE POTOUT ( IPULL, N, SIGMA, RCUT, DELTV, DELTW )     **
C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON DESTRUCTION      **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** REAL FUNCTION RANF ( DUMMY )  (GIVEN IN F.11)                 **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON ZERO TO ONE            **
C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
C    **    UPDATES LOCATE AFTER ADDITION (GIVEN IN F.14)              **
C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
C    **    UPDATES LOCATE AFTER REMOVAL (GIVEN IN F.14)               **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** ROUTINES IN AND OUT SHOULD BE CALLED WITH EQUAL PROBABILITY   **
C    ** IN A GRAND CANONICAL MONTE CARLO SIMULATION. IF A TRIAL       **
C    ** CREATION IS ACCEPTED THEN CREATE IS SET TO TRUE. IF A TRIAL   **
C    ** DESTRUCTION IS ACCEPTED THEN GHOST IS SET TO TRUE. THE        **
C    ** ROUTINES ARE WRITTEN FOR LENNARD-JONES ATOMS. THE BOX IS OF   **
C    ** UNIT LENGTH, ALL DISTANCES ARE SCALED TO THE BOX LENGTH.      **
C    ** TRIAL INPUTS WHICH RESULT IN A SEPARATION OF LESS THAN        **
C    ** 0.5*SIGMA ARE REJECTED. THE LONG-RANGE CORRECTIONS ARE        **
C    ** INCLUDED IN V AND W. ALL ACCUMULATORS ARE UPDATED IN THE MAIN **
C    ** PART OF THE PROGRAM WHICH IS NOT GIVEN HERE.                  **
C    *******************************************************************

        SUBROUTINE IN ( TEMP, Z, SIGMA, RCUT, N, V, W, CREATE )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / LOCATE

C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL CREATION                           **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP            TEMPERATURE                           **
C    ** REAL    Z               ABSOLUTE ACTIVITY                     **
C    ** REAL    SIGMA           LENNARD-JONES DIAMETER                **
C    ** REAL    RCUT            CUT-OFF DISTANCE                      **
C    ** REAL    V               POTENTIAL ENERGY                      **
C    ** REAL    W               VIRIAL                                **
C    ** INTEGER N               NUMBER OF ATOMS BEFORE TRIAL CREATION **
C    ** LOGICAL CREATE          TRUE FOR A SUCCESSFUL CREATION        **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 500 )

        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        REAL        TEMP, Z, SIGMA, RCUT, V, W
        INTEGER     LOCATE(NMAX)
        INTEGER     N
        LOGICAL     CREATE

        REAL        BETA, RXNEW, RYNEW, RZNEW, DELTV, DELTW, DELTCB
        REAL        RANF, DUMMY, RMIN
        INTEGER     NTRIAL
        LOGICAL     OVRLAP

C    *******************************************************************

        CREATE = .FALSE.
        BETA   = 1.0 / TEMP
        RMIN   = 0.5 * SIGMA
        NTRIAL = N + 1

        IF ( NTRIAL .GE. NMAX ) STOP 'MAXIMUM NUMBER OF ATOMS IN BOX'

C    ** GENERATE THE POSITION OF THE TRIAL ATOM **

        RXNEW  = RANF ( DUMMY ) - 0.5
        RYNEW  = RANF ( DUMMY ) - 0.5
        RZNEW  = RANF ( DUMMY ) - 0.5

C    ** CALCULATE ENERGY CHANGE ON ADDITION **

        CALL POTIN ( RXNEW, RYNEW, RZNEW, N, SIGMA, RCUT, RMIN, DELTV,
     :               DELTW, OVRLAP )

C    ** CHECK FOR ACCEPTANCE **

        IF ( .NOT. OVRLAP ) THEN

           DELTCB = BETA * DELTV - LOG ( Z / REAL ( NTRIAL ) )

           IF ( DELTCB .LT. 75.0 ) THEN

               IF ( DELTCB .LE. 0.0 ) THEN

                 CREATE = .TRUE.

                 CALL ADD ( RXNEW, RYNEW, RZNEW, N )

                 V    = V + DELTV
                 W    = W + DELTW
                 N    = NTRIAL

              ELSE IF ( EXP ( - DELTCB ) .GT. RANF ( DUMMY ) ) THEN

                 CREATE = .TRUE.

                 CALL ADD ( RXNEW, RYNEW, RZNEW, N )

                 V    = V + DELTV
                 W    = W + DELTW
                 N    = NTRIAL

              ENDIF

           ENDIF

        ENDIF

        RETURN
        END



        SUBROUTINE OUT ( TEMP, Z, SIGMA, RCUT, N, V, W, GHOST )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / LOCATE

C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL DESTRUCTION                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP         TEMPERATURE                              **
C    ** REAL    Z            ABSOLUTE ACTIVITY                        **
C    ** REAL    SIGMA        LENNARD-JONES DIAMETER                   **
C    ** REAL    RCUT         CUT-OFF DISTANCE                         **
C    ** REAL    V            POTENTIAL ENERGY                         **
C    ** REAL    W            VIRIAL                                   **
C    ** INTEGER N            NUMBER OF ATOMS BEFORE TRIAL DESTRUCTION **
C    ** LOGICAL GHOST        TRUE FOR A SUCCESSFUL DESTRUCTION        **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 500 )

        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        REAL        TEMP, Z, SIGMA, RCUT, V, W
        INTEGER     LOCATE(NMAX)
        INTEGER     N
        LOGICAL     GHOST

        REAL        BETA, DELTV, DELTW, DELTDB, RANF, DUMMY
        INTEGER     NTRIAL, NLOC, IPULL
        LOGICAL     OVRLAP

C    *******************************************************************

        GHOST  = .FALSE.
        BETA   = 1.0 / TEMP
        NTRIAL = N - 1

        IF ( NTRIAL .EQ. 1 ) STOP 'ONLY ONE ATOM REMAINS'

C    ** PICK A RANDOM ELEMENT FROM THE ACTIVE PART OF LOCATE **

        NLOC  = INT ( REAL ( NTRIAL ) * RANF ( DUMMY ) ) + 1
        IPULL = LOCATE(NLOC)

C    ** CALCULATE ENERGY CHANGE ON REMOVAL OF ATOM IPULL **

        CALL POTOUT ( IPULL, N, SIGMA, RCUT, DELTV, DELTW )

C    ** CHECK FOR ACCEPTANCE **

        DELTDB = BETA * DELTV - LOG ( REAL ( N ) / Z )

        IF ( DELTDB .LT. 75.0 ) THEN

           IF ( DELTDB .LT. 0.0 ) THEN

              GHOST = .TRUE.

              CALL REMOVE ( NLOC, IPULL, N )

              V = V + DELTV
              W = W + DELTW
              N = NTRIAL

           ELSE IF ( EXP( -DELTDB ) .GT. RANF ( DUMMY ) ) THEN

              GHOST = .TRUE.

              CALL REMOVE ( NLOC, IPULL, N )

              V = V + DELTV
              W = W + DELTW
              N = NTRIAL

           ENDIF

        ENDIF

        RETURN
        END



        SUBROUTINE POTIN ( RXI, RYI, RZI, N, SIGMA, RCUT, RMIN, DELTV,
     :                     DELTW, OVRLAP )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / LOCATE

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE ON ADDING AN ATOM.        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE ADDITION **
C    ** REAL    RXI,RYI,RZI       THE COORDINATES OF THE ADDED ATOM   **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    ** REAL    RMIN              MINIMUM ALLOWED APPROACH OF ATOMS   **
C    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL ADDITION OF AN ATOM TO THE FLUID. THE LONG     **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 500 )
        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        REAL        RCUT, RMIN, SIGMA, RXI, RYI, RZI, DELTV, DELTW
        INTEGER     N, NTRIAL, LOCATE(NMAX)
        LOGICAL     OVRLAP

        REAL        RCUTSQ, RMINSQ, SIGSQ, SR2, SR6, SR3, SR9
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************

        OVRLAP = .FALSE.
        RCUTSQ = RCUT * RCUT
        RMINSQ = RMIN * RMIN
        SIGSQ  = SIGMA * SIGMA

C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SIGCUB = SIGSQ * SIGMA
        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
        WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0

C    ** LOOP OVER ALL ATOMS  **

        DO 100 J = 1, N

C       ** PICK ACTIVE ATOMS FROM THE ARRAY LOCATE **

           JIN   = LOCATE(J)

           RXIJ  = RXI - RX(JIN)
           RYIJ  = RYI - RY(JIN)
           RZIJ  = RZI - RZ(JIN)

           RXIJ  = RXIJ - ANINT ( RXIJ )
           RYIJ  = RYIJ - ANINT ( RYIJ )
           RZIJ  = RZIJ - ANINT ( RZIJ )

           RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

           IF ( RIJSQ .LT. RMINSQ) THEN

              OVRLAP = .TRUE.
              RETURN

           ELSEIF ( RIJSQ .LT. RCUTSQ ) THEN

              SR2   = SIGSQ / RIJSQ
              SR6   = SR2 * SR2 * SR2
              VIJ   = SR6 * ( SR6 - 1.0 )
              WIJ   = SR6 * ( SR6 - 0.5 )
              DELTV = DELTV + VIJ
              DELTW = DELTW + WIJ

           ENDIF

100     CONTINUE

        DELTV = 4.0  * DELTV
        DELTW = 48.0 * DELTW / 3.0

C    ** ADD CHANGE IN LONG RANGE CORRECTION **

        DELTV = DELTV + ( 2.0 * REAL ( N ) + 1.0 ) * VLRC0
        DELTW = DELTW + ( 2.0 * REAL ( N ) + 1.0 ) * WLRC0

        RETURN
        END



        SUBROUTINE POTOUT ( IPULL, N, SIGMA, RCUT, DELTV, DELTW )

        COMMON / BLOCK1 / RX, RY, RZ
        COMMON / BLOCK2 / LOCATE

C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 500 )
        REAL        RX(NMAX), RY(NMAX), RZ(NMAX)
        REAL        RCUT, SIGMA, DELTV, DELTW
        INTEGER     N, IPULL, LOCATE(NMAX)

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )

C     ******************************************************************

        RCUTSQ = RCUT * RCUT
        SIGSQ  = SIGMA * SIGMA
        SIGCUB = SIGSQ * SIGMA

C    ** CALCULATE LONG RANGE CORRECTIONS **
C    ** NOTE: SPECIFIC TO LENNARD-JONES  **

        SR3    = ( SIGMA / RCUT ) ** 3
        SR9    = SR3 ** 3
        VLRC0  = ( 8.0  / 9.0 ) * PI * SIGCUB * (     SR9 - 3.0*SR3 )
        WLRC0  = ( 16.0 / 9.0 ) * PI * SIGCUB * ( 2.0*SR9 - 3.0*SR3 )

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0

        RXI    = RX(IPULL)
        RYI    = RY(IPULL)
        RZI    = RZ(IPULL)

C    ** LOOP OVER ALL ATOMS  EXCEPT IPULL **

        DO 100 J = 1, N

C       ** PICK ACTIVE ATOMS FROM LOCATE **

           JIN = LOCATE(J)

           IF ( JIN .NE. IPULL ) THEN

              RXIJ  = RXI - RX(JIN)
              RYIJ  = RYI - RY(JIN)
              RZIJ  = RZI - RZ(JIN)

              RXIJ  = RXIJ - ANINT ( RXIJ )
              RYIJ  = RYIJ - ANINT ( RYIJ )
              RZIJ  = RZIJ - ANINT ( RZIJ )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

              IF ( RIJSQ .LT. RCUTSQ ) THEN

                 SR2   = SIGSQ / RIJSQ
                 SR6   = SR2 * SR2 * SR2
                 VIJ   = SR6 * ( SR6 - 1.0 )
                 WIJ   = SR6 * ( SR6 - 0.5 )
                 DELTV = DELTV + VIJ
                 DELTW = DELTW + WIJ

              ENDIF

           ENDIF

100     CONTINUE

        DELTV =  4.0 * DELTV
        DELTW = 48.0 * DELTW / 3.0

C    ** ADD CHANGE IN LONG RANGE CORRECTION **

        DELTV = DELTV + ( 2.0 * REAL ( N ) - 1.0 ) * VLRC0
        DELTW = DELTW + ( 2.0 * REAL ( N ) - 1.0 ) * WLRC0

C    ** CHANGE SIGN OF DELTV AND DELTW FOR A REMOVAL **

        DELTV = - DELTV
        DELTW = - DELTW

        RETURN
        END



