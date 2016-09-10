********************************************************************************
** FICHE F.20.  ROUTINES TO CONSTRUCT AND USE CELL LINKED-LIST METHOD         **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** CONSTRUCTION OF CELL LINKED-LISTS AND USE IN FORCE ROUTINE.   **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** QUENTREC AND BROT, J. COMPUT. PHYS. 13, 430, 1975.            **
C    ** HOCKNEY AND EASTWOOD, COMPUTER SIMULATION USING PARTICLES,    **
C    **    MCGRAW HILL, 1981.                                         **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE MAPS                                               **
C    **    SETS UP MAP OF CELL STRUCTURE FOR USE IN FORCE             **
C    ** SUBROUTINE LINKS ( RCUT )                                     **
C    **    SETS UP HEAD OF CHAIN ARRAY AND LINKED LIST                **
C    ** SUBROUTINE FORCE ( SIGMA, RCUT, V, W )                        **
C    **    CALCULATES FORCES USING A LINKED LIST                      **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SUBROUTINE MAPS IS CALLED ONCE AT THE START OF A SIMULATION   **
C    ** TO ESTABLISH CELL NEIGHBOUR IDENTITIES.  AT EACH TIMESTEP,    **
C    ** SUBROUTINE LINKS IS CALLED TO SET UP THE LINKED LIST AND THIS **
C    ** IS IMMEDIATELY USED BY SUBROUTINE FORCE.                      **
C    *******************************************************************



        SUBROUTINE MAPS

        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** ROUTINE TO SET UP A LIST OF NEIGHBOURING CELLS                **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
C    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE SETS UP A LIST OF THE THIRTEEN NEIGHBOURING   **
C    ** CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX. THE      **
C    ** EFFECTS OF THE PERIODIC BOUNDARY CONDITIONS ARE INCLUDED.     **
C    ** THE SUBROUTINE IS CALLED ONCE AT THE BEGINNING OF THE         **
C    ** SIMULATION AND THE MAP IS USED IN THE FORCE SUBROUTINE        **
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ
        PARAMETER ( N = 1372, M = 5, NCELL = M * M * M )
        PARAMETER ( MAPSIZ = 13 * NCELL )

        INTEGER     LIST(N), HEAD(NCELL), MAP(MAPSIZ)
        INTEGER     IX, IY, IZ, IMAP, ICELL

C    *******************************************************************

C    ** STATEMENT FUNCTION TO GIVE CELL INDEX **

        ICELL ( IX, IY, IZ) = 1 + MOD ( IX - 1 + M, M )
     :                          + MOD ( IY - 1 + M, M ) * M
     :                          + MOD ( IZ - 1 + M, M ) * M * M

C    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

        DO 50 IZ = 1, M

           DO 40 IY = 1, M

              DO 30 IX = 1, M

                 IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 13

                 MAP( IMAP + 1  ) = ICELL( IX + 1, IY    , IZ     )
                 MAP( IMAP + 2  ) = ICELL( IX + 1, IY + 1, IZ     )
                 MAP( IMAP + 3  ) = ICELL( IX    , IY + 1, IZ     )
                 MAP( IMAP + 4  ) = ICELL( IX - 1, IY + 1, IZ     )
                 MAP( IMAP + 5  ) = ICELL( IX + 1, IY    , IZ - 1 )
                 MAP( IMAP + 6  ) = ICELL( IX + 1, IY + 1, IZ - 1 )
                 MAP( IMAP + 7  ) = ICELL( IX    , IY + 1, IZ - 1 )
                 MAP( IMAP + 8  ) = ICELL( IX - 1, IY + 1, IZ - 1 )
                 MAP( IMAP + 9  ) = ICELL( IX + 1, IY    , IZ + 1 )
                 MAP( IMAP + 10 ) = ICELL( IX + 1, IY + 1, IZ + 1 )
                 MAP( IMAP + 11 ) = ICELL( IX    , IY + 1, IZ + 1 )
                 MAP( IMAP + 12 ) = ICELL( IX - 1, IY + 1, IZ + 1 )
                 MAP( IMAP + 13 ) = ICELL( IX    , IY    , IZ + 1 )

30            CONTINUE

40         CONTINUE

50      CONTINUE

        RETURN
        END



        SUBROUTINE LINKS ( RCUT )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** ROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS    **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF ATOMS                    **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **
C    ** INTEGER LIST(N)            LINKED LIST OF ATOMS               **
C    ** INTEGER HEAD(NCELL)        HEAD OF CHAIN FOR EACH CELL        **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** EACH ATOM IS SORTED INTO ONE OF THE M**3 SMALL CELLS.         **
C    ** THE FIRST ATOM IN EACH CELL IS PLACED IN THE HEAD ARRAY.      **
C    ** SUBSEQUENT ATOMS ARE PLACED IN THE LINKED LIST ARRAY.         **
C    ** ATOM COORDINATES ARE ASSUMED TO BE BETWEEN -0.5 AND +0.5.     **
C    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ
        PARAMETER ( N = 1372, M = 5, NCELL = M * M * M )
        PARAMETER ( MAPSIZ = 13 * NCELL )

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        INTEGER     HEAD(NCELL), LIST(N), MAP(MAPSIZ)

        REAL        CELLI, RCUT, CELL
        INTEGER     ICELL, I

C    *******************************************************************

C    ** ZERO HEAD OF CHAIN ARRAY **

        DO 10 ICELL = 1, NCELL

           HEAD(ICELL) = 0

10      CONTINUE

        CELLI = REAL ( M )
        CELL  = 1.0 / CELLI

        IF ( CELL .LT. RCUT ) THEN

           STOP ' CELL SIZE TOO SMALL FOR CUTOFF '

        ENDIF

C    ** SORT ALL ATOMS **

        DO 20 I = 1, N

           ICELL = 1 + INT ( ( RX(I) + 0.5 ) * CELLI )
     :               + INT ( ( RY(I) + 0.5 ) * CELLI ) * M
     :               + INT ( ( RZ(I) + 0.5 ) * CELLI ) * M * M

           LIST(I)     = HEAD(ICELL)
           HEAD(ICELL) = I

20      CONTINUE

        RETURN
        END



        SUBROUTINE FORCE ( SIGMA, RCUT, V, W )

        COMMON / BLOCK1 / RX, RY, RZ, VX, VY, VZ, FX, FY, FZ
        COMMON / BLOCK2 / LIST, HEAD, MAP

C    *******************************************************************
C    ** ROUTINE TO COMPUTE FORCES AND POTENTIAL USING A LINK LIST     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  NUMBER OF ATOMS                    **
C    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
C    ** INTEGER NCELL              NUMBER OF SMALL CELLS (M**3)       **
C    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
C    ** INTEGER LIST(N)            THE LINKED LIST                    **
C    ** INTEGER HEAD(NCELL)        THE HEAD OF CHAIN ARRAY            **
C    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    ** REAL    FX(N),FY(N),FZ(N)  FORCES                             **
C    ** REAL    SIGMA              THE LJ LENGTH PARAMETER            **
C    ** REAL    RCUT               THE CUT-OFF DISTANCE               **
C    ** REAL    V                  THE POTENTIAL ENERGY               **
C    ** REAL    W                  THE VIRIAL                         **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** FORCE IS CALLED IN AN MD PROGRAM TO CALCULATE THE FORCE ON    **
C    ** EACH ATOM. THE ROUTINE IS WRITTEN FOR A LIQUID OF LENNARD     **
C    ** JONES ATOMS. SUBROUTINE FORCE REQUIRES A LINKED LIST, SET UP  **
C    ** USING SUBROUTINE LINKS, AND THE MAP OF THE SMALL CELLS SET UP **
C    ** USING SUBROUTINE MAPS.                                        **
C    *******************************************************************

        INTEGER     N, M, NCELL, MAPSIZ
        PARAMETER ( N = 1372, M = 5, NCELL = M * M * M )
        PARAMETER ( MAPSIZ = 13 * NCELL)

        REAL        RX(N), RY(N), RZ(N)
        REAL        VX(N), VY(N), VZ(N)
        REAL        FX(N), FY(N), FZ(N)
        INTEGER     HEAD(NCELL), LIST(N), MAP(MAPSIZ)

        REAL        RCUT, SIGMA, V, W

        REAL        RXI, RYI, RZI, FXIJ, FYIJ, FZIJ, FIJ, RCUTSQ
        REAL        SIGSQ, FXI, FYI, FZI, SR2, SR6, VIJ, WIJ
        REAL        RIJSQ, RXIJ, RYIJ, RZIJ
        INTEGER     ICELL, JCELL0, JCELL, I, J, NABOR

C    *******************************************************************

        SIGSQ  = SIGMA ** 2
        RCUTSQ = RCUT ** 2

C    ** ZERO FORCES AND POTENTIAL **

        DO 10 I = 1, N

           FX(I) = 0.0
           FY(I) = 0.0
           FZ(I) = 0.0

10      CONTINUE

        V = 0.0
        W = 0.0

C    ** LOOP OVER ALL CELLS **

        DO 5000 ICELL = 1, NCELL

           I = HEAD(ICELL)

C       ** LOOP OVER ALL MOLECULES IN THE CELL **

1000       IF ( I .GT. 0 ) THEN

              RXI = RX(I)
              RYI = RY(I)
              RZI = RZ(I)
              FXI = FX(I)
              FYI = FY(I)
              FZI = FZ(I)

C          ** LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL **

              J = LIST(I)

2000          IF ( J .GT. 0 ) THEN

                 RXIJ  = RXI - RX(J)
                 RYIJ  = RYI - RY(J)
                 RZIJ  = RZI - RZ(J)

                 RXIJ  = RXIJ - ANINT ( RXIJ )
                 RYIJ  = RYIJ - ANINT ( RYIJ )
                 RZIJ  = RZIJ - ANINT ( RZIJ )
                 RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    SR2   = SIGSQ / RIJSQ
                    SR6   = SR2 * SR2 * SR2
                    VIJ   = SR6 * ( SR6 - 1.0 )
                    WIJ   = SR6 * ( SR6 - 0.5 )
                    V     = V + VIJ
                    W     = W + WIJ
                    FIJ   = WIJ / RIJSQ
                    FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ
                    FXI   = FXI + FXIJ
                    FYI   = FYI + FYIJ
                    FZI   = FZI + FZIJ
                    FX(J) = FX(J) - FXIJ
                    FY(J) = FY(J) - FYIJ
                    FZ(J) = FZ(J) - FZIJ

                 ENDIF

                 J = LIST(J)

                 GO TO 2000

              ENDIF

C          ** LOOP OVER NEIGHBOURING CELLS **

              JCELL0 = 13 * (ICELL - 1)

              DO 4000 NABOR = 1, 13

                 JCELL = MAP ( JCELL0 + NABOR )

C             ** LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS **

                 J = HEAD(JCELL)

3000             IF ( J .NE. 0 ) THEN

                    RXIJ  = RXI - RX(J)
                    RYIJ  = RYI - RY(J)
                    RZIJ  = RZI - RZ(J)

                    RXIJ  = RXIJ - ANINT ( RXIJ )
                    RYIJ  = RYIJ - ANINT ( RYIJ )
                    RZIJ  = RZIJ - ANINT ( RZIJ )
                    RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                    IF ( RIJSQ. LT. RCUTSQ ) THEN

                       SR2   = SIGSQ / RIJSQ
                       SR6   = SR2 * SR2 * SR2
                       VIJ   = SR6 * ( SR6 - 1.0 )
                       WIJ   = SR6 * ( SR6 - 0.5 )
                       V     = V + VIJ
                       W     = W + WIJ
                       FIJ   = WIJ / RIJSQ
                       FXIJ  = FIJ * RXIJ
                       FYIJ  = FIJ * RYIJ
                       FZIJ  = FIJ * RZIJ
                       FXI   = FXI + FXIJ
                       FYI   = FYI + FYIJ
                       FZI   = FZI + FZIJ
                       FX(J) = FX(J) - FXIJ
                       FY(J) = FY(J) - FYIJ
                       FZ(J) = FZ(J) - FZIJ

                    ENDIF

                    J = LIST(J)

                    GO TO 3000

                 ENDIF

4000          CONTINUE

              FX(I) = FXI
              FY(I) = FYI
              FZ(I) = FZI

              I = LIST(I)

              GO TO 1000

           ENDIF

5000    CONTINUE

C    ** INCORPORATE ENERGY FACTORS **

        V = 4.0  * V
        W = 48.0 * W / 3.0

        DO 50 I = 1, N

           FX(I) = 48.0 * FX(I)
           FY(I) = 48.0 * FY(I)
           FZ(I) = 48.0 * FZ(I)

50      CONTINUE

        RETURN
        END



