********************************************************************************
** FICHE F.1.   PERIODIC BOUNDARY CONDITIONS IN VARIOUS GEOMETRIES            **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************


C    *******************************************************************
C    ** THREE ROUTINES ILLUSTRATING THE IMPLEMENTATION OF PERIODIC    **
C    ** BOUNDARY CONDITIONS FOR SIMULATIONS IN DIFFERENT GEOMETRIES.  **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** ADAMS DJ, CCP5 QUARTERLY, 10, 30, 1983.                       **
C    ** SMITH W, CCP5 QUARTERLY, 10, 37, 1983.                        **
C    ** TALBOT J, PRIVATE COMMUNICATION, 1987.                        **
C    **                                                               **
C    ** SUPPLIED ROUTINES:                                            **
C    **                                                               **
C    ** SUBROUTINE TOBOUN ( I )                                       **
C    **    IMPLEMENTS THE PERIODIC BOUNDARY CONDITIONS FOR MOLECULE   **
C    **    I IN A TRUNCATED OCTAHEDRAL BOX                            **
C    ** SUBROUTINE RDBOUN ( I )                                       **
C    **    IMPLEMENTS THE PERIODIC BOUNDARY CONDITIONS FOR MOLECULE   **
C    **    I IN A RHOMBIC DODECAHEDRAL BOX                            **
C    ** SUBROUTINE RHBOUN ( I )                                       **
C    **    IMPLEMENTS THE PERIODIC BOUNDARY CONDITIONS FOR MOLECULE   **
C    **    I IN A TWO DIMENSIONAL RHOMBOIDAL BOX.                     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                  THE NUMBER OF ATOMS                **
C    ** INTEGER I                  A PARTICULAR ATOM                  **
C    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THESE SUBROUTINES ARE CALLED AFTER THE ATOMS HAVE BEEN MOVED  **
C    ** IN AN MC OR MD CODE. IF AN ATOM HAS MOVED OUT OF THE BOX, IT  **
C    ** WILL BE REPLACED IN THE CENTRAL BOX ACCORDING TO THE PERIODIC **
C    ** BOUNDARY CONDITIONS. THIS CODE HAS BEEN WRITTEN AS SEPARATE   **
C    ** SUBROUTINES FOR CONVENIENCE. IN YOUR PARTICULAR APPLICATION   **
C    ** IT MAY NOT BE APPROPRIATE TO CALL IT AS A SUBROUTINE.         **
C    ** SIMILAR CODE MAY BE USED TO CALCULATE THE MINIMUM IMAGE       **
C    ** DISTANCE BETWEEN A PAIR I, J. IN THIS CASE RX(I) SHOULD       **
C    ** BE REPLACED BY RXIJ, RY(I) BY RYIJ, AND RZ(I) BY RZIJ.        **
C    *******************************************************************



        SUBROUTINE TOBOUN ( I )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** PERIODIC BOUNDARY CONDITIONS FOR A TRUNCATED OCTAHEDRON       **
C    **                                                               **
C    ** THE BOX IS CENTRED AT THE ORIGIN. THE AXES PASS THROUGH THE   **
C    ** CENTRES OF THE SIX SQUARE FACES OF THE TRUNCATED OCTAHEDRON   **
C    ** (SEE F1G. 1.10(A)). THE CONTAINING CUBE IS OF UNIT LENGTH     **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     I
        REAL        RX(N), RY(N), RZ(N)
        REAL        CORR, R75
        PARAMETER ( R75 = 4.0 / 3.0 )

C    *******************************************************************

        RX(I) = RX(I) - ANINT ( RX(I) )
        RY(I) = RY(I) - ANINT ( RY(I) )
        RZ(I) = RZ(I) - ANINT ( RZ(I) )

        CORR = 0.5 * AINT ( R75 * ( ABS ( RX(I) ) +
     :                              ABS ( RY(I) ) +
     :                              ABS ( RZ(I) ) ) )

        RX(I) = RX(I) - SIGN ( CORR, RX(I) )
        RY(I) = RY(I) - SIGN ( CORR, RY(I) )
        RZ(I) = RZ(I) - SIGN ( CORR, RZ(I) )

        RETURN
        END



        SUBROUTINE RDBOUN ( I )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** PERIODIC BOUNDARY CONDITIONS FOR A RHOMBIC DODECAHEDRON       **
C    **                                                               **
C    ** THE BOX IS CENTRED AT THE ORIGIN. THE X AND Y AXES JOIN THE   **
C    ** CENTRES OF OPPOSITE FACES OF THE DODECAHEDRON. THE Z AXIS     **
C    ** JOINS OPPOSITE VERTICES OF THE RHOMBIC DODECAHEDRON (SEE FIG. **
C    ** 1.10(B)). THE DIAGONAL OF THE RHOMBIC FACE IS OF UNIT LENGTH  **
C    ** AND THE SIDE OF THE CONTAINING CUBE IS SQRT(2.0).             **
C    ** NOTE THAT THE X AND Y AXES PASS THROUGH THE CUBE EDGES, WHILE **
C    ** THE Z AXIS PASSES THROUGH THE CUBE FACES.                     **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    RT2                SQRT(2.0) TO MACHINE ACCURACY      **
C    ** REAL    RRT2               1.0/SQRT(2.0)                      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     I
        REAL        RX(N), RY(N), RZ(N)

        REAL        RT2, RRT2, CORR
        PARAMETER ( RT2 = 1.4142136, RRT2 = 1.0 / RT2 )

C    *******************************************************************

        RX(I) = RX(I) - ANINT ( RX(I) )
        RY(I) = RY(I) - ANINT ( RY(I) )
        RZ(I) = RZ(I) - RT2 * ANINT ( RRT2 * RZ(I) )

        CORR = 0.5 * AINT ( ( ABS ( RX(I) ) +
     :                        ABS ( RY(I) ) +
     :                  RT2 * ABS ( RZ(I) ) )

        RX(I) = RX(I) - SIGN ( CORR, RX(I) )
        RY(I) = RY(I) - SIGN ( CORR, RY(I) )
        RZ(I) = RZ(I) - SIGN ( CORR, RZ(I) ) * RT2

        RETURN
        END



        SUBROUTINE RHBOUN ( I )

        COMMON / BLOCK1 / RX, RY, RZ

C    *******************************************************************
C    ** PERIODIC BOUNDARY CONDITIONS FOR A RHOMBIC BOX.               **
C    **                                                               **
C    ** PERIODIC CORRECTIONS ARE APPLIED IN TWO DIMENSIONS X, Y.      **
C    ** IN MOST APPLICATIONS THE MOLECULES WILL BE CONFINED IN THE    **
C    ** Z DIRECTION BY REAL WALLS RATHER THAN BY PERIODIC BOUNDARIES. **
C    ** THE BOX IS CENTRED AT THE ORIGIN. THE X AXIS LIES ALONG THE   **
C    ** SIDE OF THE RHOMBUS, WHICH IS OF UNIT LENGTH.                 **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    RT3                SQRT(3.0) TO MACHINE ACCURACY      **
C    ** REAL    RT32               SQRT(3.0)/2.0                      **
C    ** REAL    RRT3               1.0/SQRT(3.0)                      **
C    ** REAL    RRT32              2.0/SQRT(3.0)                      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 108 )

        INTEGER     I
        REAL        RX(N), RY(N), RZ(N)

        REAL        RT3, RRT3, RT32, RRT32
        PARAMETER ( RT3 = 1.7320508, RRT3 = 1.0 / RT3 )
        PARAMETER ( RT32 = RT3 / 2.0, RRT32 = 1.0 / RT32 )

C    *******************************************************************

        RX(I) = RX(I) - ANINT ( RX(I) - RRT3 * RY(I) )
     :                - ANINT ( RRT32 * RY(I) ) * 0.5
        RY(I) = RY(I) - ANINT ( RRT32 * RY(I) ) * RT32

        RETURN
        END



