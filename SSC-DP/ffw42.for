
c       This subroutine calculates the Gamma function impulse response,
c       sampled pulse response of duration D (Nash model) or the
c       averaged pulse response of duration D averaged over the
c	duration D (Nash model).
c
c
        Subroutine ffw42(h,d,m,an,ak,nopt,nacc,maxs,nprint,ier)
	  implicit none
c
	  integer maxw,idebug,nopt,m,nacc,maxs,nprint,ier
        parameter(maxw=3000)
	  integer iunito,i,j,k,nsamp,ifail,ntotal,iere,nt
        real h(200),z(maxw),s(maxw)
	  real an,ak,d,delt,t,const,SS,tbyk,x
!        Common/syspar/iuniti,iunito,idebug
c
        idebug=0
	  iunito=101

	  if(nopt.lt.1.or.nopt.gt.3)nopt = 1
        if(idebug.ne.0)then
        Write(iunito,1000)
1000    Format(/,15x,'Debug information',/,15x,30('-'))
        write(iunito,1002)d,m,an,ak,nopt,nacc,maxs
1002    Format(/,15x,'Subroutine ffw42 called',/,15x,'D',t21,'=',
     1  e13.4,/,15x,'m',t21,'=',i4,/,15x,'an',t21,'=',e13.4,/,15x,
     2  'ak',t21,'=',e13.4,/,15x,'nopt',t21,'=',i1,/,15x,
     3  'nacc',t21,'=',i3,/,15x,'maxs',t21,'=',i5)
        endif
c
        do 2 i=1,m
         h(i)=0.0
2       continue
        do 3 i=1,maxs
          z(i)=0.0
          s(i)=0.0
3       continue
c
        if(nacc.lt.1)nacc=1
        if(d.lt.0.0)d=1.0
        delt=d/nacc
        nsamp = nacc*m
        if((nsamp+nacc+1).gt.maxs)then
        write(iunito,998)maxs
998     Format(/,20x,'Dimensional error in ffw42 - maxs = ',i5)
        ier = 999
        go to 999
        endif
        t = 0
        tbyk = delt/ak
c
	call gamma(an,x,ifail)
        if (ifail.ne.0) then
        write(iunito,1003) ifail
1003    format(5x,'Error detected in gamma - ier =',i3)
        ier = 998
        goto 999
        endif
c
        z(1) = 0.0
        const = 1.0/(x*ak)
        if(nopt.eq.1)then
        if(abs(an-1.0).le.0.00005)z(1) = const
        do 5 j = 2,nsamp
        t = t+tbyk
        z(j) = const*exp(-t)*t**(an-1.0)
5       continue
        go to 100
        endif
c
        ntotal = nsamp+nacc+1
        do 10 j = 2,ntotal
        t = t+tbyk
        iere = 0
        call incgam(t,an,x,ss,iere)
        if(iere.ne.0)then
        write(iunito,996)iere
996     Format(/,20x,'Error detected in incgam - ier =',i3)
        ier = 997
        go to 999
        endif
        s(j) = ss
10      continue
        do 15 j = 2,nacc+1
        z(j) = s(j)/d
15      continue
        do 16 j = nacc+2,nsamp
        z(j) = (s(j)-s(j-nacc))/d
16      continue
        if(nopt.eq.2)go to 100
        h(1) = 0.0
        do 50 j = 1,nacc+1
50      h(1) = h(1)+z(j)
        h(1) = h(1)/nacc
        nt = nacc+1
        do 20 j = 2,m
        h(j) = 0.0
        do 25 k = nt+1,nt+nacc
        h(j) = h(j)+z(k)
25      continue
        nt = nt+nacc
        h(j) = h(j)/nacc
20      Continue
        go to 60
100     h(1) = z(1)
        k = nacc+1
        do 30 j = 2,m
        if(k.gt.nsamp)go to 60
        h(j) = z(k)
        k = k+nacc
30      Continue
60      if(nprint.ne.0)then
        write(iunito,520)d,delt,an,ak,m
520     Format(/,10x,'TABLE',/,10x,20('-'),//
     1  /,10x,'Method : NASH MODEL',
     2  /,10x,'Data : Pulsed / Sampled output and pulsed input',
     3  /,10x,'Duration',t30,'=',f10.2,
     4  /,10x,'Sampling interval',t30,'=',f10.2,
     5  /,10x,'Model parameters',/,15x,'n',t30,'=',f10.2,
     6  /,15x,'k',t30,'=',f10.2,
     7  /,10x,'No. of ordinates',t30,'=',i3,//,
     8  /,10x,'Catchment : ',//)
        if(nopt.eq.1)write(iunito,500)
        if(nopt.eq.2)write(iunito,505)
        if(nopt.eq.3)write(iunito,510)
500     Format(10x,'The sampled impulse response function',
     1  /,10x,38('-'))
505     Format(10x,'The sampled pulse response function',
     1  /,10x,36('-'))
510     Format(10x,'The averaged pulse response function',
     1  /,10x,37('-'))
        write(iunito,525)(h(i),i=1,m)
525     Format(10x,6f10.3)
        endif
999     return
        end
c
c
C================================================================
      SUBROUTINE INCGAM(X,P,GP,D6,IFAULT)
C
C     ADAPTATION OF MOORE'S ALG AS 187 IN APPL STAT, FOR COMPUTING
C     INCOMPLETE GAMMA INTEGRAL & ITS DERIVATIVES. THIS ADAPTATION
C     ONLY USES THOSE PORTIONS PERTAINING TO THE INCOMPLETE GAMMA
C     INTEGRAL AND DOES NOT INCLUDE THE DERIVATIVES.
C
C     X - ARGUEMENT OF INC GAMMA FN
C     P - SHAPE FACTOR OF INC GAMMA FN
C     GP- GAMMA FN OF P
C     D6- INC GAMMA FN OF X, FOR SHAPE FACTOR P
C     IFAULT- ERROR CODE
C
C
C     THE FOLLOWING CHANGES HAVE BEEN MADE IN THIS ADAPTATION:
C     1) AS 187 RESULTS WERE STORED IN VECTOR D(6); RESULTS HERE
C        FOR INCOMPLETE GAMMA INTEGRAL ARE STORED IN SCALAR D6;
C     2) EXPRESSIONS FOR VARIABLE "F" HAVE BEEN REPLACED WITH
C        EXACT EQUIVALENT EXPRESSIONS;
C     3)
C
C==============================================================
      DIMENSION PN(6)
      DATA E,OFLO,TMAX,ZERO/1.0E-6,1.0E30,100.0,1.0E-30/
      IFAULT=0
      PM1=P-1.0
C     XLOG=ALOG(X)
      IF (X.GT.1.0.AND.X.GE.P) GOTO 30
C     F=EXP(P*XLOG-GP1LOG-X)
      F=(X**P)*EXP(-X)/(P*GP)
      TMAXP=TMAX-P
      C=1.0
      S=1.0
      A=P
    1 A=A+1.0
      C=C*X/A
      S=S+C
      IF (A.GT.TMAXP) GOTO 1001
      IF (C.GT.E*S) GOTO 1
      D6=S*F
      RETURN
C
C        CONTINUED FRACTION EXPANSION
C
   30 F=(X**P)*EXP(-X)/GP
C
      A=PM1
      B=X+1.0-A
      TERM=0.0
      PN(1)=1.0
      PN(2)=X
      PN(3)=X+1.0
      PN(4)=X*B
      S0=PN(3)/PN(4)
   32 A=A-1.0
      B=B+2.0
      TERM=TERM+1.0
      AN=A*TERM
      PN(5)=B*PN(3)+AN*PN(1)
      PN(6)=B*PN(4)+AN*PN(2)
C
      IF (ABS(PN(6)).LT.ZERO) GOTO 35
      S=PN(5)/PN(6)
      C=ABS(S-S0)
      IF (C*P.GT.E) GOTO 34
      IF (C.LE.E*S) GOTO 42
C
   34 S0=S
   35 DO 36 I=1,4
      I2=I+2
      PN(I)=PN(I2)
   36 CONTINUE
C
      IF (TERM.GT.TMAX) GOTO 1001
      IF (ABS(PN(5)).LT.OFLO) GOTO 32
      DO 41 I=1,4
      PN(I)=PN(I)/OFLO
   41 CONTINUE
      GOTO 32
C
   42 D6=1.0-F*S
      RETURN
C
C         SET FAULT INDICATOR
C
 1001 IFAULT=1
      RETURN
      END
c
c
	SUBROUTINE GAMMA(X,GAM,IER)
C- Calculates the Gamma function for positive arguments.
C- M.Bruen, Hydrology Dept., UCG, 1983
C- Cf. Hart,J.F. et al., "Computer Approximations", Wiley, 1968
C-
C- The calling arguments are
C-
C- 	X	The argument of the Gamma function.
C-		 It must be a positive number less than 34.8
C-	GAM	The resulting value of the Gamma function.
C-	IER	An error parameter which indicates if this
C-		 routine is successful, IER=0, or fails,
C-		 IER > 0
C-
	IER=999
	IF(X.LT.0.0)RETURN
	IF(X.GT.34.8)RETURN
C- THIS GIVES A GAMMA WHICH IS JUST LESS THAN THE MAXIMUM
C- WHICH CAN BE STORED IN THE DEC-20
	IER=0
	IF(X.LE.20.0) GOTO 10
C- STIRLINGS FORM
	Y=1.0/(X*X)
	P=(0.77783067E-3*Y-0.277765545E-2)*Y+0.8333333309E-1
	P=P/X
	GAM=(X-0.5)*ALOG(X)-X+0.9189385+P
	GAM=EXP(GAM)
	RETURN
10	Y=AINT(X)
	N=Y-2.
	Y=X-Y
	GAM=(((.1082985985E-1*Y-.3427052255E-2)*Y+.77549276E-1)*Y)
	GAM=(((GAM+.8017824769E-1)*Y+.4121029027)*Y+.4227663678)*Y
	GAM=GAM+1.000000199
	T1=1.0
	YP2=Y+2.0
	IF(N)40,70,60
40	CONTINUE
C- NEGATIVE N
	N=IABS(N)
	DO 45 I=1,N
45	T1=T1*(YP2-I)
	T1=1.0/T1
	GOTO 70
60	CONTINUE
C- POSITIVE N
	N=N-1
	DO 65 I=0,N
65	T1=T1*(YP2+I)
70	GAM=GAM*T1
	RETURN
	END
