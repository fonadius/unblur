      module napack_sort2

      contains
C
C      ________________________________________________________
C     |                                                        |
C     |            SORT AN ARRAY IN INCREASING ORDER           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY OF NUMBERS                       |
C     |                                                        |
C     |         W     --WORKING ARRAY (LENGTH  AT LEAST N)     |
C     |                                                        |
C     |         N     --NUMBER OF ARRAY ELEMENTS TO SORT       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --ORIGINAL ARRAY                         |
C     |                                                        |
C     |         Y     --INDICES OF X GIVING INCREASING ORDER   |
C     |________________________________________________________|
C
      pure subroutine sort2(X,Y,W,N)
        real(kind=8),   intent(inout)   ::  x(*)
        real(kind=8),   intent(inout)   ::  y(*)
        real(kind=8),   intent(inout)   ::  w(*)
        integer,    intent(in)      ::  n
        !
      REAL S,T
      INTEGER I,J,K,L,M,P,Q
      I = 1
10    K = I
20    J = I
      Y(I) = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      W(K) = I
      GOTO 10
30    IF ( K .EQ. 1 ) RETURN
      W(K) = N + 1
40    M = 1
      L = 1
50    I = L
      IF ( I .GT. N ) GOTO 120
      P = Y(I)
      S = X(P)
      J = W(I)
      K = J
      IF ( J .GT. N ) GOTO 100
      Q = Y(J)
      T = X(Q)
      L = W(J)
      Y(I) = L
60    IF ( S .GT. T ) GOTO 70
      W(M) = P
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      P = Y(I)
      S = X(P)
      GOTO 60
70    W(M)= Q
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      Q = Y(J)
      T = X(Q)
      GOTO 60
80    W(M) = Q
      K = M + L - J
      I = J - M
90    M = M + 1
      IF ( M .EQ. K ) GOTO 50
      W(M) = Y(M+I)
      GOTO 90
100   Y(I) = J
      L = J
110   W(M) = P
      K = M + K - I
      I = I - M
      GOTO 90
120   I = 1
130   K = I
      J = Y(I)
140   Y(I) = W(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      W(K) = I
      IF ( I .LE. N ) GOTO 130
      IF ( K .GT. 1 ) GOTO 40
      RETURN
      end subroutine

      end module napack_sort2
