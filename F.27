********************************************************************************
** FICHE F.27.  PROGRAM TO COMPUTE TIME CORRELATION FUNCTIONS                 **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************



        PROGRAM TCORR

        COMMON / BLOCK1 / STORX, STORY, STORZ
        COMMON / BLOCK2 / VX, VY, VZ
        COMMON / BLOCK3 / VACF, ANORM

C    *******************************************************************
C    ** CALCULATION OF TIME CORRELATION FUNCTIONS.                    **
C    **                                                               **
C    ** THIS PROGRAM ANALYZES DATA TO CALCULATE A TIME CORRELATION    **
C    ** FUNCTION IN ONE SWEEP (LOW MEMORY REQUIREMENT). IN THIS       **
C    ** EXAMPLE THE VELOCITY AUTO-CORRELATION FUNCTION IS CALCULATED. **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER  N                  NUMBER OF ATOMS                   **
C    ** INTEGER  NSTEP              NUMBER OF STEPS ON THE TAPE       **
C    ** INTEGER  IOR                INTERVAL FOR TIME ORIGINS         **
C    ** INTEGER  NT                 CORRELATION LENGTH, INCLUDING T=0 **
C    ** INTEGER  NTIMOR             NUMBER OF TIME ORIGINS            **
C    ** INTEGER  NLABEL             LABEL FOR STEP (1,2,3...NSTEP)    **
C    ** REAL     VX(N),VY(N),VZ(N)  VELOCITIES                        **
C    ** REAL     VACF(NT)           THE CORRELATION FUNCTION          **
C    ** NSTEP AND NT SHOULD BE MULTIPLES OF IOR.                      **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE STORE ( J1 )                                       **
C    **    ROUTINE TO STORE THE DATA FOR CORRELATION                  **
C    ** SUBROUTINE CORR ( J1, J2, IT )                                **
C    **    ROUTINE TO CORRELATE THE STORED TIME ORIGINS               **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** DATA IN FILE DFILE ON FORTRAN UNIT DUNIT.                     **
C    ** RESULTS IN FILE RFILE ON FORTRAN UNIT RUNIT.                  **
C    *******************************************************************

        INTEGER     N, NSTEP, IOR, NT, NDIM, DUNIT, RUNIT, NTIMOR
        INTEGER     FULLUP
        PARAMETER ( N = 256, NSTEP = 1000, IOR = 4, NT = 200  )
        PARAMETER ( DUNIT = 10, RUNIT = 11                    )
        PARAMETER ( NDIM = NT / IOR + 1, NTIMOR = NSTEP / IOR )
        PARAMETER ( FULLUP = NDIM - 1                         )

        REAL        VX(N), VY(N), VZ(N)
        REAL        STORX(NDIM,N), STORY(NDIM,N), STORZ(NDIM,N)
        REAL        VACF(NT), ANORM(NT)
        INTEGER     S(NTIMOR), TM(NTIMOR)
        INTEGER     TS, TSS, L, NINCOR, K, JA, IB, IN, IA, JO, I
        INTEGER     NLABEL
        CHARACTER   DFILE * 30
        CHARACTER   RFILE * 30

C    *******************************************************************

        WRITE(*,'('' **** PROGRAM TCORR ****                       '')')
        WRITE(*,'('' CALCULATION OF TIME CORRELATION FUNCTIONS     '')')

C    ** READ IN FILE NAMES **

        WRITE(*,'('' ENTER DATA FILE NAME                          '')')
        READ (*,'(A)') DFILE
        WRITE (*,'('' ENTER RESULTS FILE NAME                      '')')
        READ (*,'(A)') RFILE

C    ** INITIALIZE COUNTERS **

        NINCOR = FULLUP
        JA = 1
        IA = 1
        IB = 1

C    ** ZERO ARRAYS **

        DO 5 I = 1, NT

           VACF(I)  = 0.0
           ANORM(I) = 0.0

5       CONTINUE

C    ** OPEN DATA FILE AND RESULTS FILE **

        OPEN ( UNIT = DUNIT, FILE = DFILE, ACCESS = 'SEQUENTIAL',
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        OPEN ( UNIT = RUNIT, FILE = RFILE, STATUS = 'NEW' )

C   ** CALCULATION BEGINS **

        DO 40 L = 1, NTIMOR

           JA   = JA + 1
           S(L) = JA - 1

           READ ( DUNIT ) NLABEL, VX, VY, VZ

           TM(L) = NLABEL

C       ** STORE STEP AS A TIME ORIGIN **

           CALL STORE ( JA )

C       ** CORRELATE THE ORIGINS IN STORE **

           DO 10 IN = IA, L

              TSS = TM(L) - TM(IN)
              TS  = TSS + 1
              JO  = S(IN) + 1
              CALL CORR ( JO, JA, TS )

10         CONTINUE

C       ** READ IN DATA BETWEEN TIME ORIGINS. THIS CAN  **
C       ** BE CONVENIENTLY STORED IN ELEMENT 1 OF THE   **
C       ** ARRAYS STORX ETC. AND CAN THEN BE CORRELATED **
C       ** WITH THE TIME ORIGINS.                       **

           DO 30 K = 1, IOR - 1

              READ ( DUNIT ) NLABEL, VX, VY, VZ

              CALL STORE ( 1 )

              DO 20 IN = IA, L

                 TSS = NLABEL - TM(IN)
                 TS  = TSS + 1
                 JO  = S(IN) + 1
                 CALL CORR ( JO, 1, TS )

20            CONTINUE

30         CONTINUE

           IF ( L .GE. FULLUP ) THEN

              IF ( L .EQ. NINCOR ) THEN

                 NINCOR = NINCOR + FULLUP
                 JA     = 1

              ENDIF

              IA = IA + 1

           ENDIF

40      CONTINUE

        CLOSE ( UNIT = DUNIT )

C    ** NORMALISE CORRELATION FUNCTIONS **

        VACF(1) = VACF(1) / ANORM(1) / REAL ( N )

        DO 50 I = 2, NT

           VACF(I) = VACF(I) / ANORM(I) / REAL ( N ) / VACF(1)

50      CONTINUE

        WRITE ( RUNIT, '('' VELOCITY ACF '')')
        WRITE ( RUNIT, '(I6,E15.6)') ( I, VACF(I), I = 1, NT )

        CLOSE ( RUNIT )

        STOP
        END



        SUBROUTINE STORE ( J1 )

        COMMON/ BLOCK1 / STORX, STORY, STORZ
        COMMON/ BLOCK2 / VX, VY, VZ

C    *******************************************************************
C    ** SUBROUTINE TO STORE TIME ORIGINS                              **
C    *******************************************************************

        INTEGER     J1
        INTEGER     N, NT, IOR, NDIM
        PARAMETER ( N = 256, NT = 200, IOR = 4 )
        PARAMETER ( NDIM = NT / IOR + 1        )

        REAL        STORX(NDIM,N), STORY(NDIM,N), STORZ(NDIM,N)
        REAL        VX(N), VY(N), VZ(N)
        INTEGER     I


        DO 10 I = 1, N

           STORX(J1,I) = VX(I)
           STORY(J1,I) = VY(I)
           STORZ(J1,I) = VZ(I)

10      CONTINUE

        RETURN
        END



        SUBROUTINE CORR ( J1, J2, IT )

        COMMON/ BLOCK1 / STORX, STORY, STORZ
        COMMON/ BLOCK3 / VACF, ANORM

C    *******************************************************************
C    ** SUBROUTINE TO CORRELATE TIME ORIGINS                          **
C    *******************************************************************

        INTEGER     J1, J2, IT
        INTEGER     N, NT, IOR, NDIM
        PARAMETER ( N = 256, NT = 200, IOR = 4 )
        PARAMETER ( NDIM = NT / IOR + 1        )

        REAL        STORX(NDIM,N), STORY(NDIM,N), STORZ(NDIM,N)
        REAL        VACF(NT), ANORM(NT)
        INTEGER     I

C    *******************************************************************

        DO 10 I = 1, N

           VACF(IT) = VACF(IT) + STORX(J1,I) * STORX(J2,I)
     :                         + STORY(J1,I) * STORY(J2,I)
     :                         + STORZ(J1,I) * STORZ(J2,I)

10      CONTINUE

        ANORM(IT) = ANORM(IT) + 1.0

        RETURN
        END



