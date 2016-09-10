********************************************************************************
** FICHE F.16.  HARD DUMB-BELL MONTE CARLO PROGRAM                            **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        PROGRAM MCBELL

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ

C    *******************************************************************
C    ** CONSTANT-NVT MONTE CARLO PROGRAM FOR HARD DUMB-BELLS.         **
C    **                                                               **
C    ** THE BOX IS OF UNIT LENGTH, -0.5 TO +0.5. THERE ARE NO LOOKUP  **
C    ** TABLES INCLUDED.                                              **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                   NUMBER OF MOLECULES               **
C    ** INTEGER NATOM               NUMBER OF ATOMS PER MOLECULE      **
C    ** INTEGER NSTEP               MAXIMUM NUMBER OF CYCLES          **
C    ** INTEGER IPRINT              PRINT INTERVAL                    **
C    ** INTEGR  ISAVE               SAVE INTERVAL                     **
C    ** INTEGER IRATIO              MAX DISPLACEMENT UPDATE INTERVAL  **
C    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
C    ** REAL    EX(N),EY(N),EZ(N)   ORIENTATIONS, UNIT AXIAL VECTOR   **
C    ** REAL    DAB(NATOM)          DISTANCE FROM COM TO NATOM        **
C    ** REAL    D                   REDUCED BOND LENGTH (D/SIGMA)     **
C    ** REAL    DENS                REDUCED DENSITY                   **
C    ** REAL    SIGMA               HARD SPHERE DIAMETER              **
C    ** REAL    DRMAX               REDUCED MAXIMUM DISPLACEMENT      **
C    ** REAL    DOTMIN              CONTROLS ANGULAR DISPLACEMENT     **
C    ** LOGICAL OVRLAP              TRUE IF DUMBBELLS OVERLAP         **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE CHECK ( SIGMA, DAB, OVRLAP )                       **
C    **    CHECKS FOR OVERLAPS IN A FLUID OF HARD DUMBELLS            **
C    ** SUBROUTINE ORIEN ( EXIOLD, EYIOLD, EZIOLD, DOTMIN, EXINEW,    **
C    **    :               EYINEW, EZINEW )                           **
C    **    PRODUCES A TRIAL RANDOM ORIENTATION FOR A MOLECULE         **
C    ** REAL FUNCTION RANF( DUMMY )                                   **
C    **    RETURNS A UNIFORM RANDOM NUMBER BETWEEN ZERO AND ONE       **
C    ** SUBROUTINE READCN ( CNFILE )                                  **
C    **    READS IN A CONFIGURATION                                   **
C    ** SUBROUTINE TEST ( RXI, RYI, RZI, I, EXI, EYI, EZI, SIGMA,     **
C    **    :                DAB, OVRLAP )                             **
C    **    CHECKS FOR OVERLAPS AFTER THE DISPLACEMENT OF MOLECULE I   **
C    ** SUBROUTINE WRITCN ( CNFILE )                                  **
C    **    WRITES OUT A CONFIGURATION                                 **
C    *******************************************************************

        INTEGER     N, NATOM
        PARAMETER ( N = 108, NATOM = 2 )

        REAL        RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)
        REAL        DAB(NATOM), DRMAX, DOTMIN, DENS, D, SIGMA, RATIO
        REAL        RXIOLD, RYIOLD, RZIOLD, RXINEW, RYINEW, RZINEW
        REAL        EXIOLD, EYIOLD, EZIOLD, EXINEW, EYINEW, EZINEW
        REAL        RANF, DUMMY, ACM, ACMMVA
        INTEGER     STEP, I, NSTEP, IRATIO, IPRINT, ISAVE
        LOGICAL     OVRLAP
        CHARACTER   TITLE*80, CNFILE*80

C    *******************************************************************

C    ** READ INPUT DATA **

        WRITE(*,'(1H1,'' **** PROGRAM MCBELL ****                 '')')
        WRITE(*,'(/   '' CONSTANT-NVT MONTE CARLO                 '')')
        WRITE(*,'(    '' FOR HARD DUMBELLS                       ''/)')
        WRITE(*,'('' ENTER THE RUN TITLE                          '')')
        READ (*,'(A)') TITLE
        WRITE(*,'('' ENTER NUMBER OF CYCLES                       '')')
        READ (*,*) NSTEP
        WRITE(*,'('' ENTER NUMBER OF CYCLES BETWEEN OUTPUT        '')')
        READ (*,*) IPRINT
        WRITE(*,'('' ENTER NUMBER OF CYCLES BETWEEN DATA SAVES    '')')
        READ (*,*) ISAVE
        WRITE(*,'('' ENTER INTERVAL FOR UPDATE OF MAX. DISPL.     '')')
        READ (*,*) IRATIO
        WRITE(*,'('' ENTER THE CONFIGURATION FILE NAME            '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'(/'' ENTER THE FOLLOWING IN LENNARD-JONES UNITS '',/)')
        WRITE(*,'('' ENTER THE DENSITY                            '')')
        READ (*,*) DENS
        WRITE(*,'('' ENTER THE MAXIMUM DISPLACEMENT               '')')
        READ (*,*) DRMAX
        WRITE(*,'('' ENTER THE REDUCED BOND LENGTH                '')')
        READ (*,*) D

C    ** WRITE INPUT DATA **

        WRITE(*,'(       //1X                    ,A     )') TITLE
        WRITE(*,'('' NUMBER OF ATOMS           '',I10   )') N
        WRITE(*,'('' NUMBER OF CYCLES          '',I10   )') NSTEP
        WRITE(*,'('' OUTPUT FREQUENCY          '',I10   )') IPRINT
        WRITE(*,'('' SAVE FREQUENCY            '',I10   )') ISAVE
        WRITE(*,'('' RATIO UPDATE FREQUENCY    '',I10   )') IRATIO
        WRITE(*,'('' CONFIGURATION FILE  NAME  '',A     )') CNFILE
        WRITE(*,'('' DENSITY                   '',F10.5 )') DENS
        WRITE(*,'('' MAX. DISPLACEMENT         '',F10.5 )') DRMAX
        WRITE(*,'('' BOND LENGTH               '',F10.5 )') D

C    ** SET DEPENDENT VARIABLES **

        SIGMA  = ( DENS / REAL ( N ) ) ** ( 1.0 / 3.0 )
        DAB(1) = D * SIGMA / 2.0
        DAB(2) = - DAB(1)
        DRMAX  = DRMAX * SIGMA
        DOTMIN = 0.2

C    ** WRITE OUT SOME USEFUL INFORMATION **

        WRITE( *, '( '' NUMBER OF MOLECULES   =  '', I10   )' )  N
        WRITE( *, '( '' NUMBER OF ATOMS       =  '', I10   )' )  NATOM
        WRITE( *, '( '' SIGMA  / BOX          =  '', F10.5 )' )  SIGMA
        WRITE( *, '( '' DAB(1) / BOX          =  '', F10.5 )' )  DAB(1)
        WRITE( *, '( '' DAB(2) / BOX          =  '', F10.5 )' )  DAB(2)
        WRITE( *, '( '' DRMAX  / BOX          =  '', F10.5 )' )  DRMAX
        WRITE( *, '( '' DOTMIN                =  '', F10.5 )' )  DOTMIN

C    ** READ IN INITIAL CONFIGURATION **

        CALL READCN ( CNFILE )

C    ** CHECK FOR OVERLAPS IN INITIAL CONFIGURATION **

        CALL CHECK ( SIGMA, DAB, OVRLAP )

        IF ( OVRLAP ) STOP 'OVERLAP IN INITIAL CONFIGURATION'

C    ** ZERO ACCUMULATORS **

        ACM    = 0.0
        ACMMVA = 0.0

        WRITE( *, '(//'' START OF MARKOV CHAIN                ''//)')
        WRITE( *, '( ''     ACM    RATIO     DRMAX     DOTMIN  '')')

C    *******************************************************************
C    ** LOOP OVER CYCLES BEGINS                                       **
C    *******************************************************************

        DO 100 STEP = 1, NSTEP

C       ** LOOP OVER MOLECULES **

           DO 99 I = 1, N

              RXIOLD = RX(I)
              RYIOLD = RY(I)
              RZIOLD = RZ(I)
              EXIOLD = EX(I)
              EYIOLD = EY(I)
              EZIOLD = EZ(I)

C          ** MOVE I AND PICKUP THE CENTRAL IMAGE **

              RXINEW = RXIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX
              RYINEW = RYIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX
              RZINEW = RZIOLD + ( 2.0 * RANF ( DUMMY ) - 1.0 ) * DRMAX

              RXINEW = RXINEW - ANINT ( RXINEW )
              RYINEW = RYINEW - ANINT ( RYINEW )
              RZINEW = RZINEW - ANINT ( RZINEW )

C          ** CHANGE THE ORIENTATION OF MOLECULE I **

              CALL ORIEN ( EXIOLD, EYIOLD, EZIOLD, DOTMIN,
     :                     EXINEW, EYINEW, EZINEW         )

C          ** CHECK FOR ACCEPTANCE **

              CALL TEST ( RXINEW, RYINEW, RZINEW, I,
     :                    EXINEW, EYINEW, EZINEW, SIGMA, DAB, OVRLAP )

              IF ( .NOT. OVRLAP ) THEN

C             ** ACCEPT MOVE **

                 RX(I) = RXINEW
                 RY(I) = RYINEW
                 RZ(I) = RZINEW
                 EX(I) = EXINEW
                 EY(I) = EYINEW
                 EZ(I) = EZINEW
                 ACMMVA = ACMMVA + 1.0

              ENDIF

              ACM = ACM + 1.0

99         CONTINUE

C       ****************************************************************
C       ** LOOP OVER MOLECULES COMPLETE                               **
C       ****************************************************************

C       ** PERFORM PERIODIC OPERATIONS  **

C       ** CHANGE MAXIMUM DISPLACEMENT **

           IF ( MOD ( STEP, IRATIO ) .EQ. 0 ) THEN

              RATIO = ACMMVA / REAL ( N * IRATIO )

              IF ( RATIO .GT. 0.5 ) THEN

                 DRMAX = DRMAX * 1.05
                 DOTMIN = DOTMIN * 1.025

              ELSE

                 DRMAX = DRMAX * 0.95
                 DOTMIN = DOTMIN * 0.975

              ENDIF

              ACMMVA = 0.0

           ENDIF

C       ** WRITE OUT RUNTIME INFORMATION **

           IF ( MOD ( STEP, IPRINT ) .EQ. 0 ) THEN

              WRITE(*,'(I8,3F10.4)') INT(ACM), RATIO, DRMAX, DOTMIN

           ENDIF

C       ** WRITE OUT THE CONFIGURATION AT INTERVALS **

           IF ( MOD ( STEP, ISAVE ) .EQ. 0 ) THEN

              CALL WRITCN ( CNFILE )

              CALL CHECK ( SIGMA, DAB, OVRLAP )

              IF ( OVRLAP ) STOP 'OVERLAP DURING THE RUN'

           ENDIF

100     CONTINUE

C    *******************************************************************
C    ** ENDS THE LOOP OVER CYCLES                                     **
C    *******************************************************************

C    ** CHECKS FOR OVRLAPS IN THE FINAL CONFIGURATION  **

        CALL CHECK ( SIGMA, DAB, OVRLAP )

        IF ( OVRLAP ) STOP 'OVERLAP IN FINAL CONFIGURATION'

C    ** WRITE OUT THE FINAL CONFIGURATION FROM THE RUN **

        CALL WRITCN ( CNFILE )

        STOP
        END



        SUBROUTINE ORIEN ( EXIOLD, EYIOLD, EZIOLD, DOTMIN,
     :                     EXINEW, EYINEW, EZINEW         )

C    *******************************************************************
C    ** FINDS A TRIAL RANDOM ORIENTATION OF A LINEAR MOLECULE.        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL      EXIOLD,EYIOLD,EZIOLD  OLD AXIAL VECTOR FOR I        **
C    ** REAL      EYINEW,EYINEW,EZINEW  NEW AXIAL VECTOR FOR I        **
C    ** REAL      DOT                   DOT PRODUCT OF OLD AND NEW    **
C    **                                 AXIAL VECTORS                 **
C    ** REAL      DOTMIN                MINIMUM ALLOWED DOT PRODUCT   **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE METHOD USE A REJECTION TECHNIQUE TO CREATE A TRIAL        **
C    ** ORIENTATION OF MOLECULE I SUBJECT TO THE CONSTRAINT THAT      **
C    ** THE COSINE OF THE ANGLE BETWEEN THE OLD AND NEW AXIAL         **
C    ** VECTORS, DOT, IS GREATER THAN ( 1.0 - DOTMIN ).               **
C    *******************************************************************

        REAL    EXIOLD, EYIOLD, EZIOLD, EXINEW, EYINEW, EZINEW, DOTMIN

        REAL    DOT, XI1, XI2, XI, XISQ
        REAL    RANF, DUMMY

C    *******************************************************************

C    ** INITIALISE DOT **

        DOT  = 0.0

C    ** ITERATIVE LOOP **

1000    IF ( ( 1.0 - DOT ) .GE. DOTMIN ) THEN

C       ** INITIALISE XISQ **

           XISQ = 1.0

C       ** INNER ITERATIVE LOOP **

2000       IF ( XISQ .GE. 1.0 ) THEN

              XI1  = RANF ( DUMMY ) * 2.0 - 1.0
              XI2  = RANF ( DUMMY ) * 2.0 - 1.0
              XISQ = XI1 * XI1 + XI2 * XI2

              GOTO 2000

           ENDIF

           XI     = SQRT ( 1.0 - XISQ )
           EXINEW = 2.0 * XI1 * XI
           EYINEW = 2.0 * XI2 * XI
           EZINEW = 1.0 - 2.0 * XISQ
           DOT    = EXINEW * EXIOLD + EYINEW * EYIOLD + EZINEW * EZIOLD

           GOTO 1000

        ENDIF

        RETURN
        END



        SUBROUTINE TEST ( RXI, RYI, RZI, I, EXI, EYI, EZI, SIGMA,
     :                    DAB, OVRLAP )

        COMMON / BLOCK1 / RX, RY,RZ, EX, EY, EZ

C    *******************************************************************
C    ** CHECKS FOR OVERLAP OF I WITH ALL OTHER MOLECULES.             **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER I                  THE MOLECULE OF INTEREST           **
C    ** INTEGER N                  NUMBER OF MOLECULES                **
C    ** INTEGER NATOM              NUMBER OF ATOMS PER MOLECULE       **
C    ** REAL    RXI,RYI,RZI        POSITION OF MOLECULE I             **
C    ** REAL    EXI,EYI,EZI,       ORIENTATION OF MOLECULE I          **
C    ** REAL    RX(N),RY(N),RZ(N)  MOLECULAR POSITIONS                **
C    ** REAL    EX(N),EY(N),EZ(N)  MOLECULAR ORIENTATIONS             **
C    ** REAL    DAB(NATOM)         POSITION OF ATOMS IN A MOLECULE    **
C    ** REAL    SIGMA              REDUCED ATOM DIAMETER              **
C    ** LOGICAL OVRLAP             TRUE IF MOLECULE I OVERLAPS        **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** CALLED AFTER A TRIAL DISPLACEMENT OF MOLECULE I TO ESTABLISH  **
C    ** WHETHER THERE IS AN OVERLAP IN THE TRIAL CONFIGURATION.       **
C    *******************************************************************

        INTEGER     NATOM, N
        PARAMETER ( NATOM = 2, N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        EX(N), EY(N), EZ(N)
        REAL        RXI, RYI, RZI, EXI, EYI, EZI
        REAL        SIGMA, DAB(NATOM)
        INTEGER     I
        LOGICAL     OVRLAP

        REAL        RXIJ, RYIJ, RZIJ, EXJ, EYJ, EZJ
        REAL        RXAB, RYAB, RZAB, DABI, SIGSQ, RABSQ
        INTEGER     J, IA, JB

C    *******************************************************************

        OVRLAP = .FALSE.
        SIGSQ  = SIGMA * SIGMA

C    ** LOOPS OVER MOLECULES EXCEPT I **

        DO 100 J = 1, N

           IF ( J .NE. I ) THEN

              EXJ = EX(J)
              EYJ = EY(J)
              EZJ = EZ(J)

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT ( RXIJ )
              RYIJ = RYIJ - ANINT ( RYIJ )
              RZIJ = RZIJ - ANINT ( RZIJ )

C          ** LOOPS OVER ATOMS **

              DO 99 IA = 1, NATOM

                 DABI = DAB(IA)

                 DO 98 JB = 1, NATOM

                    RXAB = RXIJ + EXI * DABI + EXJ * DAB(JB)
                    RYAB = RYIJ + EYI * DABI + EYJ * DAB(JB)
                    RZAB = RZIJ + EZI * DABI + EZJ * DAB(JB)

                    RABSQ = RXAB * RXAB + RYAB * RYAB + RZAB * RZAB

                    IF ( RABSQ .LT. SIGSQ ) THEN

                       OVRLAP = .TRUE.
                       RETURN

                    ENDIF

98               CONTINUE

99            CONTINUE

           ENDIF

100     CONTINUE

        RETURN
        END



        SUBROUTINE CHECK ( SIGMA, DAB, OVRLAP )

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ

C    *******************************************************************
C    ** ROUTINE TO CHECK FOR OVERLAPS IN A FLUID OF HARD DUMBBELLS    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 NUMBER OF MOLECULES                 **
C    ** INTEGER NATOM             NUMBER OF ATOMS PER MOLECULE        **
C    ** REAL    RX(N),RY(N),RZ(N) MOLECULAR POSITIONS                 **
C    ** REAL    EX(N),EY(N),EZ(N) MOLECULAR ORIENTATIONS              **
C    ** REAL    DAB(NATOM)        POSITION OF ATOMS IN A MOLECULE     **
C    ** REAL    SIGMA             REDUCED ATOM DIAMETER               **
C    ** LOGICAL OVRLAP            TRUE IF TWO DUMBELLS OVERLAP        **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** CALLED AT INTERVALS DURING THE RUN TO CHECK FOR OVERLAPS. IF  **
C    ** OVRLAP IS RETURNED WITH A TRUE VALUE THEN THERE IS AN ERROR   **
C    ** IN THE PROGRAM AND THE EXECUTION IS STOPPED.                  **
C    *******************************************************************

        INTEGER     NATOM, N
        PARAMETER ( NATOM = 2, N = 108 )

        REAL        RX(N), RY(N), RZ(N)
        REAL        EX(N), EY(N), EZ(N)
        REAL        SIGMA, DAB(NATOM)
        LOGICAL     OVRLAP

        REAL        RXI, RYI, RZI, RXIJ, RYIJ, RZIJ, EXI, EYI, EZI
        REAL        EXJ, EYJ, EZJ, RXAB, RYAB, RZAB, DABI, SIGSQ, RABSQ
        INTEGER     I, J, IA, JB

C    *******************************************************************

        OVRLAP = .FALSE.
        SIGSQ  = SIGMA * SIGMA

C    ** LOOPS OVER MOLECULES **

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           EXI = EX(I)
           EYI = EY(I)
           EZI = EZ(I)

           DO 99 J = I + 1, N

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT ( RXIJ )
              RYIJ = RYIJ - ANINT ( RYIJ )
              RZIJ = RZIJ - ANINT ( RZIJ )

              EXJ  = EX(J)
              EYJ  = EY(J)
              EZJ  = EZ(J)

C          ** LOOPS OVER ATOMS **

              DO 98 IA = 1, NATOM

                 DABI = DAB(IA)

                 DO 97 JB = 1, NATOM

                    RXAB = RXIJ + EXI * DABI + EXJ * DAB(JB)
                    RYAB = RYIJ + EYI * DABI + EYJ * DAB(JB)
                    RZAB = RZIJ + EZI * DABI + EZJ * DAB(JB)

                    RABSQ = RXAB * RXAB + RYAB * RYAB + RZAB * RZAB

                    IF ( RABSQ .LT. SIGSQ ) THEN

                       OVRLAP = .TRUE.
                       RETURN

                    ENDIF

97               CONTINUE

98            CONTINUE

99         CONTINUE

100     CONTINUE

        RETURN
        END



        SUBROUTINE READCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN THE CONFIGURATION FROM UNIT 10          **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )
        CHARACTER   CNFILE*(*)
        REAL        RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)

        INTEGER     CNUNIT
        PARAMETER ( CNUNIT = 10 )
        INTEGER     NN

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'OLD',
     :         FORM = 'UNFORMATTED'                        )

        READ ( CNUNIT ) NN
        IF ( NN .NE. N ) STOP 'N ERROR IN READCN'
        READ ( CNUNIT ) RX, RY, RZ
        READ ( CNUNIT ) EX, EY, EZ

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END



        SUBROUTINE WRITCN ( CNFILE )

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ

C    *******************************************************************
C    ** SUBROUTINE TO WRITE OUT THE CONFIGURATION TO UNIT 10          **
C    *******************************************************************

        INTEGER      N
        PARAMETER (  N = 108 )
        CHARACTER    CNFILE*(*)
        REAL         RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)

        INTEGER      CNUNIT
        PARAMETER (  CNUNIT = 10 )

C   ********************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE, STATUS = 'UNKNOWN',
     :         FORM = 'UNFORMATTED'                        )

        WRITE ( CNUNIT ) N
        WRITE ( CNUNIT ) RX, RY, RZ
        WRITE ( CNUNIT ) EX, EY, EZ

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



