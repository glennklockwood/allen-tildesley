********************************************************************************
** FICHE F.35.  THE VORONOI CONSTRUCTION IN 2D AND 3D.                        **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** TWO SEPARATE PARTS: TWO AND THREE DIMENSIONAL VERSIONS.       **
C    *******************************************************************



C    *******************************************************************
C    ** FICHE F.35  -  PART A                                         **
C    ** THE VORONOI CONSTRUCTION IN 2D.                               **
C    *******************************************************************



        PROGRAM VORON2

        COMMON / BLOCK1 / RX, RY

C    *******************************************************************
C    ** CONSTRUCTION OF THE VORONOI POLYGON IN 2D.                    **
C    **                                                               **
C    ** THIS PROGRAM TAKES IN A CONFIGURATION IN A SQUARE BOX WITH    **
C    ** CONVENTIONAL PERIODIC BOUNDARY CONDITIONS AND FOR EACH ATOM   **
C    ** OBTAINS THE SURROUNDING VORONOI POLYGON, DEFINED AS THAT      **
C    ** REGION OF SPACE CLOSER TO THE CHOSEN ATOM THAN TO ANY OTHER.  **
C    ** NEIGHBOURING POLYGONS DEFINE NEIGHBOURING ATOMS.              **
C    ** THE PROGRAM IS SLOW BUT ESSENTIALLY FOOLPROOF.                **
C    ** WE USE THE MINIMUM IMAGE CONVENTION AND SET A CUTOFF BEYOND   **
C    ** WHICH ATOMS ARE ASSUMED NOT TO BE NEIGHBOURS: BOTH OF THESE   **
C    ** MEASURES ARE DANGEROUS FOR SMALL AND/OR RANDOM SYSTEMS.       **
C    ** WE DELIBERATELY DO NOT USE PREVIOUSLY-FOUND NEIGHBOURS IN     **
C    ** CONSTRUCTING NEIGHBOUR LISTS, SO THAT AN INDEPENDENT CHECK    **
C    ** MAY BE MADE AT THE END.                                       **
C    ** HERE WE SIMPLY PRINT OUT THE GEOMETRICAL INFORMATION AT THE   **
C    ** END.  THE OUTPUT IS QUITE LENGTHY.  IN PRACTICE, IT WOULD     **
C    ** PROBABLY BE ANALYZED DIRECTLY WITHOUT PRINTING OUT.           **
C    ** NB: BEWARE DEGENERATE CONFIGURATIONS, I.E. ONES IN WHICH MORE **
C    ** THAN THREE VORONOI DOMAINS SHARE A VERTEX. THE SQUARE LATTICE **
C    ** IS AN EXAMPLE.                                                **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                        NUMBER OF ATOMS              **
C    ** REAL    RX(N),RY(N)              POSITIONS                    **
C    ** REAL    PX(MAXCAN),PY(MAXCAN)    CANDIDATE RELATIVE POSITIONS **
C    ** REAL    PS(MAXCAN)               SQUARED RELATIVE DISTANCES   **
C    ** INTEGER NVER                     NUMBER OF VERTICES FOUND     **
C    ** INTEGER NEDGE                    NUMBER OF EDGES FOUND        **
C    ** INTEGER VERTS(MAXCAN)            VERTICES FOR EACH CANDIDATE  **
C    **                                  = 0 IF NOT A NEIGHBOUR       **
C    **                                  = 2 ( 1 EDGE ) IF NEIGHBOUR  **
C    ** REAL    RXVER(MAXVER)            VERTEX RELATIVE X-COORD      **
C    ** REAL    RYVER(MAXVER)            VERTEX RELATIVE Y-COORD      **
C    ** INTEGER IVER(MAXVER)             ATOMIC INDICES TAGGING       **
C    ** INTEGER JVER(MAXVER)             .. EACH VERTEX OF POLYGON    **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE, N, BOX )                          **
C    **    READS IN CONFIGURATION, NUMBER OF ATOMS, BOX SIZE          **
C    ** SUBROUTINE SORT ( MAXCAN, PX, PY, PS, TAG, NCAN )             **
C    **    SORTS NEIGHBOUR DETAILS INTO ASCENDING DISTANCE ORDER      **
C    ** SUBROUTINE WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE,          **
C    **     PX, PY, PS, VERTS, RXVER, RYVER, IVER, JVER )             **
C    **    CARRIES OUT THE VORONOI CONSTRUCTION                       **
C    *******************************************************************

        INTEGER     MAXN, MAXCAN, MAXVER
        PARAMETER ( MAXN = 108, MAXCAN = 50, MAXVER = 50 )

        REAL        RX(MAXN), RY(MAXN)

        REAL        PX(MAXCAN), PY(MAXCAN), PS(MAXCAN)
        INTEGER     TAG(MAXCAN), VERTS(MAXCAN)

        REAL        RXVER(MAXVER), RYVER(MAXVER)
        INTEGER     IVER(MAXVER), JVER(MAXVER)
        INTEGER     NABLST(MAXVER,MAXN), NNAB(MAXN), INAB, JNAB

        INTEGER     NCAN, NVER, NCOORD, NEDGE
        INTEGER     I, J, CAN, VER, N
        REAL        BOX, BOXINV, RCUT, RCUTSQ, COORD
        REAL        RXJ, RYJ, RZJ, RXIJ, RYIJ, RZIJ, RIJSQ
        CHARACTER   CNFILE*30
        LOGICAL     OK

C    *******************************************************************

        WRITE(*,'(1H1,'' **** PROGRAM VORON2 ****                  '')')
        WRITE(*,'(//1X,''VORONOI CONSTRUCTION IN 2D                '')')

C    ** BASIC PARAMETERS **

        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE

C    ** READCN MUST READ IN INITIAL CONFIGURATION  **

        CALL READCN ( CNFILE, N, BOX )

        WRITE(*,'(1X,I5,''-ATOM CONFIGURATION'')') N
        WRITE(*,'('' BOX LENGTH = '',F10.5)') BOX
        WRITE(*,'('' ENTER NEIGHBOUR CUTOFF IN SAME UNITS '')')
        READ (*,*) RCUT
        WRITE(*,'('' NEIGHBOUR CUTOFF = '',F10.5)') RCUT

        RCUTSQ = RCUT ** 2
        BOXINV = 1.0 / BOX

C    ** ZERO ACCUMULATORS **

        DO 100 J = 1, N

           NNAB(J) = 0

           DO 90 INAB = 1, NVER

              NABLST(INAB,J) = 0

90         CONTINUE

100     CONTINUE

C    *******************************************************************
C    ** MAIN LOOP STARTS                                              **
C    *******************************************************************

        DO 1000 J = 1, N

           IF ( MOD ( J, 2 ) .EQ. 0 ) THEN

              WRITE(*,'(///1X,''RESULTS FOR ATOM '',I5)') J

           ELSE

              WRITE(*,'(1H1,''RESULTS FOR ATOM '',I5)') J

           ENDIF

           RXJ = RX(J)
           RYJ = RY(J)
           CAN = 0

C       ** SELECT CANDIDATES **

           DO 500 I = 1, N

              IF ( I .NE. J ) THEN

                 RXIJ = RX(I) - RXJ
                 RYIJ = RY(I) - RYJ
                 RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
                 RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
                 RIJSQ  = RXIJ ** 2 + RYIJ ** 2

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    CAN = CAN + 1

                    IF ( CAN .GT. MAXCAN ) THEN

                       WRITE(*,'('' TOO MANY CANDIDATES '')')
                       STOP

                    ENDIF

                    PX(CAN)  = RXIJ
                    PY(CAN)  = RYIJ
                    PS(CAN)  = RIJSQ
                    TAG(CAN) = I

                 ENDIF

              ENDIF

500        CONTINUE

C       ** CANDIDATES HAVE BEEN SELECTED **

           NCAN = CAN

C       ** SORT INTO INCREASING DISTANCE ORDER **
C       ** THIS SHOULD IMPROVE EFFICIENCY      **

           CALL SORT ( MAXCAN, PX, PY, PS, TAG, NCAN )

C       ** PERFORM VORONOI CONSTRUCTION **

           CALL WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE,
     :                 PX, PY, PS, VERTS,
     :                 RXVER, RYVER, IVER, JVER )

C       ** WRITE OUT RESULTS **

           WRITE(*,'(/1X,''NUMBER OF NEIGHBOURS '',I5)') NEDGE
           WRITE(*,'(/1X,''NEIGHBOUR LIST '')')
           WRITE(*,10001)

           DO 800 CAN = 1, NCAN

              IF ( VERTS(CAN) .NE. 0 ) THEN

                 PS(CAN) = SQRT ( PS(CAN) )
                 WRITE(*,'(1X,I5,3X,I5,3X,2F12.5,3X,F12.5)')
     :              TAG(CAN), VERTS(CAN), PX(CAN), PY(CAN), PS(CAN)
                 NNAB(J) = NNAB(J) + 1
                 NABLST(NNAB(J),J) = TAG(CAN)

              ENDIF

800        CONTINUE

           WRITE(*,'(/1X,''NUMBER OF VERTICES   '',I5)') NVER
           WRITE(*,'(/1X,''VERTEX LIST '')')
           WRITE(*,10002)

           DO 900 VER = 1, NVER

              WRITE(*,'(1X,2I5,3X,2F12.5)')
     :           TAG(IVER(VER)), TAG(JVER(VER)),
     :           RXVER(VER), RYVER(VER)

900        CONTINUE

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(1H1,''FINAL SUMMARY'')')
        WRITE(*,10003)

        NCOORD = 0

        DO 2000 J = 1, N

           NCOORD = NCOORD + NNAB(J)

           WRITE(*,'(1X,I5,3X,I5,3X,30I3)') J, NNAB(J),
     :       ( NABLST(INAB,J), INAB = 1, NNAB(J) )

C       ** CHECK THAT IF I IS A NEIGHBOUR OF J **
C       ** THEN J IS ALSO A NEIGHBOUR OF I     **

           DO 1500 INAB = 1, NNAB(J)

              I = NABLST(INAB,J)

              OK = .FALSE.
              JNAB = 0

1200          IF ( ( .NOT. OK ) .AND. ( JNAB .LE. NNAB(I) ) ) THEN

                 OK = ( J .EQ. NABLST(JNAB,I) )
                 JNAB = JNAB + 1
                 GOTO 1200

              ENDIF

              IF ( .NOT. OK ) THEN

                 WRITE(*,'(1X,I3,'' IS NOT A NEIGHBOUR OF '',I3)') J, I

              ENDIF

1500       CONTINUE

2000    CONTINUE

        COORD = REAL ( NCOORD ) / REAL ( N )

        WRITE(*,'(/1X,'' AVERAGE COORDINATION NUMBER = '',F10.5)') COORD

        STOP

10001   FORMAT(/1X,'ATOM ',3X,'EDGE ',
     :         /1X,'INDEX',3X,'VERTS',3X,
     :         '      RELATIVE POSITION   ',3X,'  DISTANCE  ')
10002   FORMAT(/1X,'   INDICES         RELATIVE POSITION ')
10003   FORMAT(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ')

        END



        SUBROUTINE READCN ( CNFILE, N, BOX )

        COMMON / BLOCK1 / RX, RY

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION                   **
C    *******************************************************************

        INTEGER     MAXN
        PARAMETER ( MAXN = 108 )

        REAL        RX(MAXN), RY(MAXN), BOX
        INTEGER     N

        CHARACTER   CNFILE*(*)

        INTEGER     CNUNIT, I
        PARAMETER ( CNUNIT = 10 )

C    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) N, BOX
        IF ( N .GT. MAXN ) STOP ' N TOO LARGE '
        READ ( CNUNIT ) ( RX(I), I = 1, N ), ( RY(I), I = 1, N )

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END




        SUBROUTINE WORK ( MAXCAN, MAXV, NN, NV, NE, RX, RY, RS, VERTS,
     :                    VX, VY, IV, JV )

C    *******************************************************************
C    ** ROUTINE TO PERFORM VORONOI ANALYSIS                           **
C    **                                                               **
C    ** WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                **
C    ** I.E. THE EDGES OF THE POLYGON GO THROUGH THE POINTS.          **
C    *******************************************************************

        INTEGER     MAXCAN, NN, MAXV, NV, NE
        INTEGER     VERTS(MAXCAN)
        REAL        RX(MAXCAN), RY(MAXCAN), RS(MAXCAN)
        REAL        VX(MAXV), VY(MAXV)
        INTEGER     IV(MAXV), JV(MAXV)

        LOGICAL     OK
        INTEGER     I, J, L, NN1, N, V
        REAL        AI, BI, CI, AJ, BJ, CJ, DET, DETINV
        REAL        VXIJ, VYIJ
        REAL        TOL
        PARAMETER ( TOL = 1.E-6 )

C    *******************************************************************

C    ** IF THERE ARE LESS THAN 3 POINTS GIVEN **
C    ** WE CANNOT CONSTRUCT A POLYGON         **

        IF ( NN .LT. 3 ) THEN

           WRITE(*,'('' LESS THAN 3 POINTS GIVEN TO WORK '',I5)') NN
           STOP

        ENDIF

        NN1 = NN - 1
        V = 0

C    ** WE AIM TO EXAMINE EACH POSSIBLE VERTEX  **
C    ** DEFINED BY THE INTERSECTION OF 2 EDGES  **
C    ** EACH EDGE IS DEFINED BY RX,RY,RS.       **

        DO 400 I = 1, NN1

           AI =  RX(I)
           BI =  RY(I)
           CI = -RS(I)

           DO 300 J = I + 1, NN

              AJ =  RX(J)
              BJ =  RY(J)
              CJ = -RS(J)

              DET = AI * BJ - AJ * BI

              IF ( ABS ( DET ) .GT. TOL ) THEN

C             ** THE EDGES INTERSECT **

                 DETINV = 1.0 / DET

                 VXIJ = ( BI * CJ - BJ * CI ) * DETINV
                 VYIJ = ( AJ * CI - AI * CJ ) * DETINV

C             ** NOW WE TAKE SHOTS AT THE VERTEX **
C             ** USING THE REMAINING EDGES ..... **

                 OK = .TRUE.
                 L  = 1

100              IF ( OK .AND. ( L .LE. NN ) ) THEN

                    IF ( ( L .NE. I ) .AND. ( L .NE. J ) ) THEN

                       OK = ( RX(L) * VXIJ + RY(L) * VYIJ ) .LE. RS(L)

                    ENDIF

                    L = L + 1
                    GOTO 100

                 ENDIF

C             ** IF THE VERTEX MADE IT      **
C             ** ADD IT TO THE HALL OF FAME **
C             ** CONVERT TO CORRECT SCALE   **

                 IF ( OK ) THEN

                    V = V + 1
                    IF ( V .GT. MAXV ) STOP 'TOO MANY VERTICES'
                    IV(V)  = I
                    JV(V)  = J
                    VX(V) = 0.5 * VXIJ
                    VY(V) = 0.5 * VYIJ

                 ENDIF

              ENDIF

300        CONTINUE

400     CONTINUE

C    ** THE SURVIVING VERTICES DEFINE THE VORONOI POLYGON **

        NV = V

        IF ( NV .LT. 3 ) THEN

           WRITE(*,'('' LESS THAN 3 VERTICES FOUND IN WORK '',I5)') NV
           STOP

        ENDIF

C    ** IDENTIFY NEIGHBOURING POINTS **

        DO 500 N = 1, NN

           VERTS(N) = 0

500     CONTINUE

        DO 600 V = 1, NV

           VERTS(IV(V)) = VERTS(IV(V)) + 1
           VERTS(JV(V)) = VERTS(JV(V)) + 1

600     CONTINUE

C    ** POINTS WITH NONZERO VERTS ARE NEIGHBOURS **
C    ** IF NONZERO, VERTS SHOULD BE EQUAL TO 2   **

C    ** CHECK RESULT AND COUNT EDGES **

        OK = .TRUE.
        NE = 0

        DO 700 N = 1, NN

           IF ( VERTS(N) .GT. 0 ) THEN

              NE = NE + 1

              IF ( VERTS(N) .NE. 2 ) THEN

                 OK = .FALSE.

              ENDIF

           ENDIF

700     CONTINUE

        IF ( .NOT. OK ) THEN

           WRITE (*,'('' **** VERTEX ERROR: DEGENERACY ? **** '')')

        ENDIF

        IF ( NE .NE. NV ) THEN

           WRITE(*,'('' **** EDGE   ERROR: DEGENERACY ? ****  '')')

        ENDIF

        RETURN
        END



        SUBROUTINE SORT ( MAXCAN, RX, RY, RS, TAG, NN )

C    *******************************************************************
C    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
C    **                                                               **
C    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
C    *******************************************************************

        INTEGER MAXCAN, NN
        REAL    RX(MAXCAN), RY(MAXCAN), RS(MAXCAN)
        INTEGER TAG(MAXCAN)

        LOGICAL CHANGE
        INTEGER I, ITOP, I1, TAGI
        REAL    RXI, RYI, RSI

C    *******************************************************************

        CHANGE = .TRUE.
        ITOP = NN - 1

1000    IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

           CHANGE = .FALSE.

           DO 100 I = 1, ITOP

              I1 = I + 1

              IF ( RS(I) .GT. RS(I1) ) THEN

                 RXI = RX(I)
                 RYI = RY(I)
                 RSI = RS(I)
                 TAGI = TAG(I)

                 RX(I) = RX(I1)
                 RY(I) = RY(I1)
                 RS(I) = RS(I1)
                 TAG(I) = TAG(I1)

                 RX(I1) = RXI
                 RY(I1) = RYI
                 RS(I1) = RSI
                 TAG(I1) = TAGI

                 CHANGE = .TRUE.

              ENDIF

100        CONTINUE

           ITOP = ITOP - 1
           GOTO 1000

        ENDIF

        RETURN
        END



C    *******************************************************************
C    ** FICHE F.35  - PART B                                          **
C    ** THE VORONOI CONSTRUCTION IN 3D.                               **
C    *******************************************************************



        PROGRAM VORON3

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** CONSTRUCTION OF VORONOI POLYHEDRON IN 3D.                     **
C    **                                                               **
C    ** THIS PROGRAM TAKES IN A CONFIGURATION IN A CUBIC BOX WITH     **
C    ** CONVENTIONAL PERIODIC BOUNDARY CONDITIONS AND FOR EACH ATOM   **
C    ** OBTAINS THE SURROUNDING VORONOI POLYHEDRON, DEFINED AS THAT   **
C    ** REGION OF SPACE CLOSER TO THE CHOSEN ATOM THAN TO ANY OTHER.  **
C    ** NEIGHBOURING POLYHEDRA DEFINE NEIGHBOURING ATOMS.             **
C    ** THIS PROGRAM IS SLOW BUT ESSENTIALLY FOOLPROOF.               **
C    ** WE USE THE MINIMUM IMAGE CONVENTION AND SET A CUTOFF BEYOND   **
C    ** WHICH ATOMS ARE ASSUMED NOT TO BE NEIGHBOURS: BOTH OF THESE   **
C    ** MEASURES ARE DANGEROUS FOR SMALL AND/OR RANDOM SYSTEMS.       **
C    ** WE DELIBERATELY DO NOT USE PREVIOUSLY-FOUND NEIGHBOURS IN     **
C    ** CONSTRUCTING NEIGHBOUR LISTS, SO THAT AN INDEPENDENT CHECK    **
C    ** MAY BE MADE AT THE END.                                       **
C    ** HERE WE SIMPLY PRINT OUT THE GEOMETRICAL INFORMATION AT THE   **
C    ** END.  THE OUTPUT IS QUITE LENGTHY.  IN PRACTICE, IT WOULD     **
C    ** PROBABLY BE ANALYZED DIRECTLY WITHOUT PRINTING IT OUT.        **
C    ** NB: BEWARE DEGENERATE CONFIGURATIONS, I.E. ONES IN WHICH MORE **
C    ** THAN FOUR VORONOI DOMAINS SHARE A VERTEX. THE SIMPLE CUBIC    **
C    ** AND FACE-CENTRED CUBIC LATTICES ARE EXAMPLES.                 **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                       NUMBER OF MOLECULES           **
C    ** REAL    RX(N),RY(N),RZ(N)       POSITIONS                     **
C    ** REAL    PX(MAXCAN), ETC.        CANDIDATE RELATIVE POSITIONS  **
C    ** REAL    PS(MAXCAN)              SQUARED RELATIVE DISTANCES    **
C    ** INTEGER NVER                    NUMBER OF VERTICES FOUND      **
C    ** INTEGER NEDGE                   NUMBER OF EDGES FOUND         **
C    ** INTEGER NFACE                   NUMBER OF FACES FOUND         **
C    ** INTEGER EDGES(MAXCAN)           EDGES PER FACE FOR CANDIDATES **
C    **                                 = 0 FOR NON-NEIGHBOURS        **
C    ** REAL    RXVER(MAXVER)           VERTEX RELATIVE X-COORD       **
C    ** REAL    RYVER(MAXVER)           VERTEX RELATIVE Y-COORD       **
C    ** REAL    RZVER(MAXVER)           VERTEX RELATIVE Z-COORD       **
C    ** INTEGER IVER(MAXVER)            ATOM INDICES DEFINING         **
C    ** INTEGER JVER(MAXVER)            .. VERTICES IN VORONOI        **
C    ** INTEGER KVER(MAXVER)            .. POLYHEDRON.                **
C    **                                                               **
C    ** ROUTINES REFERENCED:                                          **
C    **                                                               **
C    ** SUBROUTINE READCN ( CNFILE, N, BOX )                          **
C    **    READS CONFIGURATION, NUMBER OF ATOMS, BOX SIZE.            **
C    ** SUBROUTINE SORT ( MAXCAN, PX, PY, PZ, PS, TAG, NCAN )         **
C    **    SORTS NEIGHBOUR DETAILS INTO ASCENDING DISTANCE ORDER      **
C    ** SUBROUTINE WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, NFACE,   **
C    **    PX, PY, PZ, PS, EDGES,                                     **
C    **    RXVER, RYVER, RZVER, IVER, JVER, KVER )                    **
C    **    CARRIES OUT VORONOI CONSTRUCTION                           **
C    *******************************************************************

        INTEGER     MAXN, MAXCAN, MAXVER
        PARAMETER ( MAXN = 108, MAXCAN = 50, MAXVER = 100 )

        REAL        RX(MAXN), RY(MAXN), RZ(MAXN)
        REAL        PX(MAXCAN), PY(MAXCAN), PZ(MAXCAN), PS(MAXCAN)
        INTEGER     TAG(MAXCAN), EDGES(MAXCAN)
        REAL        RXVER(MAXVER), RYVER(MAXVER), RZVER(MAXVER)
        INTEGER     IVER(MAXVER), JVER(MAXVER), KVER(MAXVER)
        INTEGER     NABLST(MAXVER,MAXN), NNAB(MAXN), INAB, JNAB

        INTEGER     NCAN, NVER, NEDGE, NFACE, NCOORD
        INTEGER     I, CAN, VER, J, N
        REAL        BOX, BOXINV, RCUT, RCUTSQ, COORD
        REAL        RXJ, RYJ, RZJ, RXIJ, RYIJ, RZIJ, RIJSQ
        CHARACTER   CNFILE*30
        LOGICAL     OK

C    *******************************************************************

        WRITE(*,'(1H1,'' **** PROGRAM VORON3 ****                  '')')
        WRITE(*,'(//1X,''VORONOI CONSTRUCTION IN 3D                '')')

C    ** BASIC PARAMETERS **

        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE

C    ** READCN MUST READ IN INITIAL CONFIGURATION  **

        CALL READCN ( CNFILE, N, BOX )

        WRITE(*,'(1X,I5,''-ATOM CONFIGURATION'')') N
        WRITE(*,'('' BOX LENGTH = '',F10.5)') BOX
        WRITE(*,'('' ENTER NEIGHBOUR CUTOFF IN SAME UNITS '')')
        READ (*,*) RCUT
        WRITE(*,'('' NEIGHBOUR CUTOFF = '',F10.5)') RCUT

        RCUTSQ = RCUT**2
        BOXINV = 1.0 / BOX

C    ** ZERO ACCUMULATORS **

        DO 100 J = 1, N

           NNAB(J) = 0

           DO 90 INAB = 1, NVER

              NABLST(INAB,J) = 0

90         CONTINUE

100     CONTINUE

C    *******************************************************************
C    ** MAIN LOOP STARTS                                              **
C    *******************************************************************

        DO 1000 J = 1, N

           WRITE(*,'(1H1,'' RESULTS FOR ATOM '',I5)') J

           RXJ = RX(J)
           RYJ = RY(J)
           RZJ = RZ(J)
           CAN = 0

C       ** SELECT CANDIDATES **

           DO 500 I = 1, N

              IF ( I .NE. J ) THEN

                 RXIJ = RX(I) - RXJ
                 RYIJ = RY(I) - RYJ
                 RZIJ = RZ(I) - RZJ
                 RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
                 RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
                 RZIJ = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX
                 RIJSQ  = RXIJ**2 + RYIJ**2 + RZIJ**2

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    CAN = CAN + 1

                    IF ( CAN .GT. MAXCAN ) THEN

                       WRITE(*,'('' TOO MANY CANDIDATES '')')
                       STOP

                    ENDIF

                    PX(CAN)  = RXIJ
                    PY(CAN)  = RYIJ
                    PZ(CAN)  = RZIJ
                    PS(CAN)  = RIJSQ
                    TAG(CAN) = I

                 ENDIF

              ENDIF

500        CONTINUE

C       ** CANDIDATES HAVE BEEN SELECTED **

           NCAN = CAN

C       ** SORT INTO ASCENDING ORDER OF DISTANCE **
C       ** THIS SHOULD IMPROVE EFFICIENCY        **

           CALL SORT ( MAXCAN, PX, PY, PZ, PS, TAG, NCAN )

C       ** PERFORM VORONOI ANALYSIS **

           CALL WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, NFACE,
     :                 PX, PY, PZ, PS, EDGES,
     :                 RXVER, RYVER, RZVER, IVER, JVER, KVER )

C       ** WRITE OUT RESULTS **

           WRITE(*,'(/1X,''NUMBER OF NEIGHBOURS '',I5)') NFACE
           WRITE(*,'(/1X,''NEIGHBOUR LIST '')')
           WRITE(*,10001)

           DO 800 CAN = 1, NCAN

              IF (EDGES(CAN) .NE. 0) THEN

                 PS(CAN) = SQRT ( PS(CAN) )
                 WRITE(*,'(1X,I5,3X,I5,3X,3F12.5,3X,F12.5)')
     :              TAG(CAN), EDGES(CAN),
     :              PX(CAN), PY(CAN), PZ(CAN), PS(CAN)

                 NNAB(J) = NNAB(J) + 1
                 NABLST(NNAB(J),J) = TAG(CAN)

              ENDIF

800        CONTINUE

           WRITE(*,'(/1X,''NUMBER OF EDGES '',I5)') NEDGE

           WRITE(*,'(/1X,''NUMBER OF VERTICES '',I5)') NVER
           WRITE(*,'(/1X,''VERTEX LIST '')')
           WRITE(*,10002)

           DO 900 VER = 1, NVER

              WRITE(*,'(1X,3I5,3X,3F12.5)')
     :        TAG(IVER(VER)), TAG(JVER(VER)), TAG(KVER(VER)),
     :        RXVER(VER), RYVER(VER), RZVER(VER)

900        CONTINUE

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        WRITE(*,'(1H1,''FINAL SUMMARY'')')
        WRITE(*,10003)

        NCOORD = 0

        DO 2000 J = 1, N

           NCOORD = NCOORD + NNAB(J)

           WRITE(*,'(1X,I5,3X,I5,3X,30I3)') J, NNAB(J),
     :        ( NABLST(INAB,J), INAB = 1, NNAB(J) )

C       ** CHECK THAT IF I IS A NEIGHBOUR OF J **
C       ** THEN J IS ALSO A NEIGHBOUR OF I     **

           DO 1500 INAB = 1, NNAB(J)

              I = NABLST(INAB,J)

              OK = .FALSE.
              JNAB = 0

1200          IF ( ( .NOT. OK ) .AND. ( JNAB .LE. NNAB(I) ) ) THEN

                 OK = ( J .EQ. NABLST(JNAB,I) )
                 JNAB = JNAB + 1
                 GOTO 1200

              ENDIF

              IF ( .NOT. OK ) THEN

                 WRITE(*,'(1X,I3,'' IS NOT A NEIGHBOUR OF '',I3)') J, I

              ENDIF

1500       CONTINUE

2000    CONTINUE

        COORD = REAL ( NCOORD ) / REAL ( N )

        WRITE(*,'(/1X,'' AVERAGE COORDINATION NUMBER = '',F10.5)') COORD

        STOP

10001   FORMAT(/1X,'ATOM ',3X,'FACE ',
     :         /1X,'INDEX',3X,'EDGES',3X,
     :         '            RELATIVE POSITION         ',3X,
     :         '  DISTANCE')
10002   FORMAT(/1X,'      INDICES           RELATIVE POSITION')
10003   FORMAT(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ')

        END



        SUBROUTINE READCN ( CNFILE, N, BOX )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION FROM UNIT 10      **
C    *******************************************************************

        INTEGER     MAXN
        PARAMETER ( MAXN = 108 )

        REAL        RX(MAXN), RY(MAXN), RZ(MAXN), BOX

        INTEGER     CNUNIT, N, I
        PARAMETER ( CNUNIT = 10 )

        CHARACTER   CNFILE*(*)


C    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) N, BOX
        IF ( N .GT. MAXN ) STOP ' N TOO LARGE '
        READ ( CNUNIT ) ( RX(I), I = 1, N ),
     :                  ( RY(I), I = 1, N ),
     :                  ( RZ(I), I = 1, N )

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END




        SUBROUTINE WORK ( MAXCAN, MAXV, NN, NV, NE, NF,
     :                    RX, RY, RZ, RS, EDGES,
     :                    VX, VY, VZ, IV, JV, KV )

C    *******************************************************************
C    ** ROUTINE TO PERFORM VORONOI ANALYSIS                           **
C    **                                                               **
C    ** WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                **
C    ** I.E. THE FACES OF THE POLYHEDRON GO THROUGH THE POINTS.       **
C    *******************************************************************

        INTEGER     MAXCAN, NN, MAXV, NV, NE, NF
        INTEGER     EDGES(MAXCAN)
        REAL        RX(MAXCAN), RY(MAXCAN), RZ(MAXCAN), RS(MAXCAN)
        REAL        VX(MAXV), VY(MAXV), VZ(MAXV)
        INTEGER     IV(MAXV), JV(MAXV), KV(MAXV)

        LOGICAL     OK
        INTEGER     I, J, K, L, NN1, NN2, N, V
        REAL        AI, BI, CI, DI, AJ, BJ, CJ, DJ, AK, BK, CK, DK
        REAL        AB, BC, CA, DA, DB, DC, DET, DETINV
        REAL        VXIJK, VYIJK, VZIJK
        REAL        TOL
        PARAMETER ( TOL = 1.E-6 )

C    *******************************************************************

C    ** IF THERE ARE LESS THAN 4 POINTS GIVEN **
C    ** WE CANNOT CONSTRUCT A POLYHEDRON      **

        IF ( NN .LT. 4 ) THEN

           WRITE(*,'('' LESS THAN 4 POINTS GIVEN TO WORK '',I5)') NN
           STOP

        ENDIF

        NN1 = NN - 1
        NN2 = NN - 2
        V = 0

C    ** WE AIM TO EXAMINE EACH POSSIBLE VERTEX  **
C    ** DEFINED BY THE INTERSECTION OF 3 PLANES **
C    ** EACH PLANE IS SPECIFIED BY RX,RY,RZ,RS  **

        DO 400 I = 1, NN2

           AI =  RX(I)
           BI =  RY(I)
           CI =  RZ(I)
           DI = -RS(I)

           DO 300 J = I + 1, NN1

              AJ =  RX(J)
              BJ =  RY(J)
              CJ =  RZ(J)
              DJ = -RS(J)

              AB = AI * BJ - AJ * BI
              BC = BI * CJ - BJ * CI
              CA = CI * AJ - CJ * AI
              DA = DI * AJ - DJ * AI
              DB = DI * BJ - DJ * BI
              DC = DI * CJ - DJ * CI

              DO 200 K = J + 1, NN

                 AK =  RX(K)
                 BK =  RY(K)
                 CK =  RZ(K)
                 DK = -RS(K)

                 DET = AK * BC + BK * CA + CK * AB

                 IF ( ABS ( DET ) .GT. TOL ) THEN

C                ** THE PLANES INTERSECT **

                    DETINV = 1.0 / DET

                    VXIJK = ( - DK * BC + BK * DC - CK * DB ) * DETINV
                    VYIJK = ( - AK * DC - DK * CA + CK * DA ) * DETINV
                    VZIJK = (   AK * DB - BK * DA - DK * AB ) * DETINV

C                ** NOW WE TAKE SHOTS AT THE VERTEX **
C                ** USING THE REMAINING PLANES .... **

                    OK = .TRUE.
                    L  = 1

100                 IF ( OK .AND. ( L .LE. NN ) ) THEN

                       IF ( ( L .NE. I ) .AND.
     :                      ( L .NE. J ) .AND.
     :                      ( L .NE. K )       ) THEN

                          OK = ( ( RX(L) * VXIJK +
     :                             RY(L) * VYIJK +
     :                             RZ(L) * VZIJK  ) .LE. RS(L) )

                       ENDIF

                       L = L + 1
                       GOTO 100

                    ENDIF

C                ** IF THE VERTEX MADE IT      **
C                ** ADD IT TO THE HALL OF FAME **
C                ** CONVERT TO CORRECT SCALE   **

                    IF ( OK ) THEN

                       V = V + 1

                       IF ( V .GT. MAXV ) STOP 'TOO MANY VERTICES'

                       IV(V)  = I
                       JV(V)  = J
                       KV(V)  = K
                       VX(V) = 0.5 * VXIJK
                       VY(V) = 0.5 * VYIJK
                       VZ(V) = 0.5 * VZIJK

                    ENDIF

                 ENDIF

200           CONTINUE

300        CONTINUE

400     CONTINUE

        NV = V

        IF ( NV .LT. 4 ) THEN

           WRITE(*,'('' LESS THAN 4 VERTICES FOUND IN WORK '',I5)') NV
           STOP

        ENDIF

C    ** IDENTIFY NEIGHBOURING POINTS **

        DO 500 N = 1, NN

           EDGES(N) = 0

500     CONTINUE

        DO 600 V = 1, NV

           EDGES(IV(V)) = EDGES(IV(V)) + 1
           EDGES(JV(V)) = EDGES(JV(V)) + 1
           EDGES(KV(V)) = EDGES(KV(V)) + 1

600     CONTINUE

C    ** POINTS WITH NONZERO EDGES ARE NEIGHBOURS **

C    ** CHECK EULER RELATION **

        NF = 0
        NE = 0

        DO 700 N = 1, NN

           IF ( EDGES(N) .GT. 0 ) NF = NF + 1
           NE = NE + EDGES(N)

700     CONTINUE

        IF ( MOD ( NE, 2 ) .NE. 0 ) THEN

           WRITE(*,'('' NONINTEGER NUMBER OF EDGES'',I5)') NE
           STOP

        ENDIF

        NE = NE / 2

        IF ( ( NV - NE + NF ) .NE. 2 ) THEN

           WRITE(*,'('' **** EULER ERROR: DEGENERACY ? **** '')')

        ENDIF

        RETURN
        END



        SUBROUTINE SORT ( MAXCAN, RX, RY, RZ, RS, TAG, NN )

C    *******************************************************************
C    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
C    **                                                               **
C    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
C    *******************************************************************

        INTEGER MAXCAN, NN
        REAL    RX(MAXCAN), RY(MAXCAN), RZ(MAXCAN), RS(MAXCAN)
        INTEGER TAG(MAXCAN)

        LOGICAL CHANGE
        INTEGER I, ITOP, I1, TAGI
        REAL    RXI, RYI, RZI, RSI

C    *******************************************************************

        CHANGE = .TRUE.
        ITOP = NN - 1

1000    IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

           CHANGE = .FALSE.

           DO 100 I = 1, ITOP

              I1 = I + 1

              IF ( RS(I) .GT. RS(I1) ) THEN

                 RXI = RX(I)
                 RYI = RY(I)
                 RZI = RZ(I)
                 RSI = RS(I)
                 TAGI = TAG(I)

                 RX(I) = RX(I1)
                 RY(I) = RY(I1)
                 RZ(I) = RZ(I1)
                 RS(I) = RS(I1)
                 TAG(I) = TAG(I1)

                 RX(I1) = RXI
                 RY(I1) = RYI
                 RZ(I1) = RZI
                 RS(I1) = RSI
                 TAG(I1) = TAGI

                 CHANGE = .TRUE.

              ENDIF

100        CONTINUE

           ITOP = ITOP - 1
           GOTO 1000

        ENDIF

        RETURN
        END



