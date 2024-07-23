C*** THIS IS A DYNAMIC RELAXATION                                       DRE00010
C*** PROGRAM                                                            DRE00020
c
c    version II of Flat prepared 4th October 1989
c
      include 'flat2.com'
      OPEN (2,STATUS='OLD',FILE='input.dat')                                DRE00150
      OPEN (3,STATUS='unknown',FILE='output.dat')
      OPEN (4,STATUS='unknown',FILE='data.dat')
      open (50,status = 'unknown',file = 'stress.dat')
      write(6,*) 'input number of required iterations'
      READ(5,*)NMAX
      WRITE(6,*) 'INPUT MASS ADDITION'
      READ(5,*)ASS
      WRITE(3,14)ASS
  14  FORMAT(' ADDED MASS COMPONENT',E15.7,//)
      WRITE(6,*) 'INPUT DENSITY'
      READ(5,*)DEN
      WRITE(3,13)DEN
 13   FORMAT(' DENSITY OF MATERIAL ',E15.7,//)
      WRITE(6,*)'OUTPUT /SLACK LENGTHS T/F'
      READ(5,*)ANSWER
      KT=0                                                              DRE00180
      DT=0.1                                                            DRE00190
      DAMP=0.0D0                                                        DRE00200
      ENGL=0.0D0
      ENG=0.0D0                                                         DRE00210
      CALL DATA                                                         DRE00230
      CALL LENG                                                         DRE00240
      IF(NE.GT.0) CALL AREAT
      IF(NE.GT.0) CALL FORM             
      WRITE(*,89)                                                       DRE00260
  89  FORMAT(///,' INITIALISE LOAD,MASS,FOR FIRST TIME')         
      IF (NSTR.GT.0)CALL AVENOR 
      IF (DEN.GT.0.0)CALL UDL
      CALL RESI
      CALL LOAD                                                         DRE00280
      CALL MASS                                                         DRE00290
      CALL INIT                                                         DRE00300
      CALL PRINT                                                        DRE00310
      WRITE(3,67)AF,B                                                   DRE00320
 67   FORMAT (' AF=',F10.3,'      B=',F10.3)                            DRE00330
      DO 100 KT=1,KMAX                                                  DRE00340
      DO 200 NT=1,NMAX                                                  DRE00350
      IT=((KT-1)*NMAX) +NT                                              DRE00360
C     WRITE(*,452) IT                                                   DRE00370
 452  FORMAT(/,' INITIALISE LOAD AND RES IT=',I5)                       DRE00380
      CALL RESI                                                         DRE00390
      CALL LOAD                                                         DRE00400
C     WRITE(*,4444) IT                                                  DRE00410
 4444 FORMAT(' INITIALISE STEP IT=',I5)                                 DRE00420
      CALL STEP                                                         DRE00430
      if(ENG.LT.1.0E-10) goto 300
      IF(DAMP.GT.0.001) GOTO 200                                        DRE00440
      IF(ENG.LT.ENGL) CALL PEAK                                         DRE00450
      ENGL=ENG                                                          DRE00460
 200  CONTINUE                                                          DRE00470
      CALL PRINT                                                        DRE00480
 100  CONTINUE                                                          DRE00500
 300  CALL RSTART                                                       DRE00510
      CALL REACTION
      STOP                                                              DRE00520
      END                                                               DRE00530
      SUBROUTINE DATA    
      include 'flat2.com'
C****  DATA INPUT ROUTINE                                               DRE00670
      READ (2,678)TITEL                                                 DRE00680
 678  FORMAT(A32)                                                       DRE00690
      WRITE(3,1)TITEL                                                   DRE00700
 1    FORMAT(      /   /,' DATA AND RESULTS  ',/,1X,A80)                DRE00710
      READ(2,3345)DT,DAMP
 3345 FORMAT(2F10.5)
      WRITE(3,2323)DT,DAMP
 2323 FORMAT(///,1X,'DT=',F10.5,5X,'DAMP='F10.5)
      READ(2,2)NN,NE,NB,NM,NL,NDIME,NSTR,NBE,KMAX,INMAX,ILN           
      read(2,2) nq
 2    FORMAT(BZ,2I6,9I5)                                                   DRE00730
      WRITE(3,3)NN,NE,NB,NM,NL,NDIME,NSTR,NBE,KMAX,NMAX                 DRE00740
 3    FORMAT(/,' NO OF NODES=',I5,/,' NO OF ELEMENTS=',I5,/,            DRE00750
     1' NO OF RESTRAINED BOUNDARY NODES=',I5,/,' NO OF MATERIAL TYPES=',DRE00760
     2I5,/,' NO OF LINKS=',I5,/,' NO OF DIMENSIONS=',I5,/,' NO OF STRS='DRE00770
     3,I5,/,' NO BEAMS=',I5,/,' NO OF ITGRPS=',I5,/,' NO OF ITBTWP=',I5)DRE00780
      WRITE(3,567)ILN                                                   DRE00790
 567  FORMAT(' O IF NO LOADS=',I5)                                      DRE00800
      DO 44 L=1,NN                                                      DRE00810
 44   READ(2,4) N,(CORD(N,I),I=1,NDIME)                                 DRE00820
 4    FORMAT(BZ,I5,3F15.5)                                              DRE00830
      DO 414 N=1,NN                                                     DRE00840
 414  WRITE(3,5) N,(CORD(N,I),I=1,NDIME )                               DRE00850
 5    FORMAT(' NODE=',I5,10X,'COORDS=',3(2X ,F10.4))                    DRE00860
      IF(NE.EQ.0) GOTO 56                                               DRE00870
      DO 5441 L=1,NE                                                    DRE00880
 5441 READ(2,6)  N,(NOD(L,I),I=1,4),(TLEN(L,I),I=1,3)                   DRE00890
 6    FORMAT(BZ,5I5,3F15.8)                                             DRE00900
      DO 9871 N=1,NE                                                    DRE00910
 9871 WRITE(3,7) N,(NOD(N,I),I=1,4),(TLEN(N,I),I=1,3)                   DRE00920
 7    FORMAT(' ELEMENT=',I3,2X,'NODE1=',I3,2X,'NODE2=',I3,2X,'NODE3='   DRE00930
     1,I3,2X,'MAT NO=',I5,'   INITIAL LENGTHS:',3(3X,F15.8))            DRE00940
  56  IF(NL.EQ.0) GOTO 57                                               DRE00950
      DO 491 L=1,NL                                                     DRE00960
 491  READ(2,76) N,(NODL(L,I),I=1,3), RLEN(L)                           DRE00970
 76   FORMAT(BZ,4I5,F15.8)                                              DRE00980
      WRITE(3,77)(N,(NODL(N,I),I=1,3),RLEN(N),N=1,NL)                   DRE00990
 77   FORMAT(' LINK=',I3,2X,'NODE1=',I3,2X,'NODE2=',I3,2X,              DRE01000
     1'MAT NO=',I5,'  INT LENGTH=',F15.8)                               DRE01010
 57   CONTINUE                                                          DRE01020
      DO 577 L=1,NM                                                     DRE01030
 577  READ(2,11) N,MAT(N),(AMAT(N,I),I=1,5)                             DRE01040
 11   FORMAT(BZ,2I5,5e12.5)                                             DRE01050
      WRITE(3,12)                                                       DRE01060
 12   FORMAT( / ,' MATERIAL DATA')                                      DRE01070
      DO 1111 N=1,NM                                                    DRE01080
 1111 WRITE(3,13)N,MAT(N),(AMAT(N,I),I=1,5)                             DRE01090
 13   FORMAT('  MATERIAL NUMBER',I2,3X,I5,5X,5e12.5)  
      DO 88 N=1,NB                                                      DRE01110
  88  READ(2,8)(NCO(N,I),I=1,(NDIME+1))                                 DRE01120
  8   FORMAT(BZ,4I5)                                                    DRE01130
      WRITE(3,9)                                                        DRE01140
 9    FORMAT(/,' BOUNDARY CONDITIONS')                                  DRE01150
      DO 100 N=1,NB                                                     DRE01160
 100  WRITE(3,10)(NCO(N,I),I=1,(NDIME+1))                               DRE01170
 10   FORMAT(  ' NODE=',I5,5X,3(5X,I2))                                 DRE01180
C*** READ IN THE APPLIED LOADS PUT IN F ARRAY                           DRE01190
C*** REMEMBER MUST HAVE LAST LINE FOR LAST NODE EVEN IF NO APPLIED LOAD DRE01200
      WRITE(3,999)                                                      DRE01210
 999  FORMAT(///, ' LOADING CONDITIONS'/)                               DRE01220
      IF(ILN.EQ.0)  GOTO 628                                            DRE01230
      ILN=0                                                             DRE01240
 777  ILN=ILN+1                                                         DRE01250
      I=ILN                                                             DRE01260
      READ(2,16)NC(I),(ALOAD(I,II),II=1,NDIME)                          DRE01270
 16   FORMAT(I5,3F15.5)                                                 DRE01280
      WRITE(3,17)NC(I),(ALOAD(I,II),II=1,NDIME)                         DRE01290
 17   FORMAT(1H ,'NODE=',I5,5X,3(4X,F15.5))                             DRE01300
      IF(NC(I).LT.NN)GOTO 777  
 628  CONTINUE                                                          DRE01320
C***** READ IN G STRINGS
      IF(NSTR.EQ.0) GOTO 727
      CALL GSDATA
      CALL AVENOR
 727  CONTINUE
      RETURN                                                            DRE01330
      END                                                               DRE01340
      SUBROUTINE MASS                                                   DRE01350
C ******* CALCULATES MASSES                                             DRE01360
C ** USES LINKS ONLY AND THEN TAKES LARGEST/10. VALUE FOR THE ELEMENTS  DRE01370
      include 'flat2.com'
      DO 50 N=1,NN                                                      DRE01500
      DO 50 K=1,NDIME                                                   DRE01510
 50   BB(N,K)=0.0D0                                                     DRE01520
C *********  CHECK STIFFNESS OF LINKS                                   DRE01530
      IF(NL.EQ.0)GOTO 7788                                              DRE01550
      DO 52 L=1,NL                                                      DRE01560
      N=NODL(L,3)                                                       DRE01570
      N1= NODL(L,1)                                                     DRE01580
      N2= NODL(L,2)                                                     DRE01590
      E=AMAT(N,2)                                                       DRE01600
      A=AMAT(N,1)
      MMM=MAT(N)
      EA=A*E 
      IF(MMM.EQ.5)EA=AMAT(N,3)
      DL=0.0                                                            DRE01620
      DO 51 K=1,NDIME                                                   DRE01630
      D(K) = CORD(N2,K) -CORD(N1,K)                                     DRE01640
  51  DL= DL+(D(K)*D(K))                                                DRE01650
      DL=DSQRT(DL)                                                      DRE01660
      DO 52 I=1,2                                                       DRE01670
      NODE=NODL(L,I)                                                    DRE01680
      DO 52 K=1,NDIME                                                   DRE01690
      XRT=((EA)*D(K)**2/(DL**3)) 
 52   BB(NODE,K) = BB(NODE,K) + XRT                                     DRE01720
 7788 CONTINUE                                                          DRE01730
      IF(NE.LT.1) GOTO 567                                              DRE01750
C ***** CALCULATE MASSES                                                DRE01770
      DO 56 N=1,NN                                                      DRE01790
      DO 56 K=1,NDIME                                                   DRE01800
      BB(N,K) = ((DT/0.5)**2)*((BB(N,K))/2.0D0)              
C      BB(N,K) =BB(N,K) +   1.0     
      BB(N,K)=ASS
  56  CONTINUE         
 567  CONTINUE
C******  FIXED NODES                                                    DRE01840
C     CALL PRINT                                                        DRE01850
      DO 1715   N=1,NB                                                  DRE01860
      NODE=NCO(N,1)                                                     DRE01870
      DO 1715 K=2,NDIME+1                                               DRE01880
      IF(NCO(N,K).EQ.0) GOTO 1715                                       DRE01890
      BB(NODE,K-1)=10.0D0 ** 25                                         DRE01900
 1715 CONTINUE                                                          DRE01910
      RETURN                                                            DRE01920
      END                                                               DRE01930
      SUBROUTINE STEP
C****** ONE ITERATION                                                   DRE01950
      include 'flat2.com'
      ENG =0.0D0                                                        DRE02080
      DO 26 N=1,NN                                                      DRE02090
      DO 26 K=1,NDIME                                                   DRE02100
      V(N,K) =(V(N,K) *AF)+(R(N,K)*DT*B/BB(N,K))                        DRE02110
      IF(DABS(V(N,K)).LT.10.0D-25) V(N,K) =0.0                          DRE02120
C     WRITE(3,788)N,K,ENG,BB(N,K),V(N,K)                                DRE02130
C788  FORMAT(1X,2(2X,I5),3(3X,F15.4))                                   DRE02140
      ENG= ENG+ ((V(N,K)**2)*BB(N,K))                                   DRE02150
 26   CORD(N,K) =CORD(N,K) +(DT*V(N,K))                                 DRE02160
      WRITE(*,789)IT,ENG                                                DRE02170
  789 FORMAT('  ENERGY AT ITERATION',I5,'    EQUALS ',E20.8)            DRE02180
      RETURN                                                            DRE02190
      END                                                               DRE02200
      SUBROUTINE INIT                                                   DRE02210
C **** INITIAL ITERATION STEP                                           DRE02220
      include 'flat2.com'
      AF= (1.0D0-(DAMP/2.0D0))/(1.0D0+(DAMP/2.0D0))                     DRE02350
      B= 1.0D0/ (1.0D0+(DAMP/2.0D0))                                    DRE02360
      DO 50 N=1,NN                                                      DRE02370
      DO 50 K=1,NDIME                                                   DRE02380
      V(N,K) =B*R(N,K)*DT/(BB(N,K)*(1.0D0 + AF))                        DRE02390
 50   CORD(N,K) =CORD(N,K) +DT*V(N,K)                                   DRE02400
      RETURN                                                            DRE02410
      END                                                               DRE02420
      SUBROUTINE LOAD                                                   DRE02430
C***** CALCULATES LOADS AND ADDS TO RESIDUALS                           DRE02440
      include 'flat2.com'
      IF(ILN.EQ.0) GOTO 5717                                            DRE02630
C******** THESE ARE APPLIED POINT LOADS                                 DRE02640
      DO 65 J=1,ILN                                                     DRE02650
      N  = NC(J)                                                        DRE02660
      DO 65 K=1,NDIME                                                   DRE02670
 65   R(N,K)=R(N,K)+ALOAD(J,K)                                          DRE02680
 5717 CONTINUE                                                          DRE02690
      IF (DEN.LT.0.0)GOTO 5678
      DO 495 N=1,NN
      DO 485 K=1,NDIME
  485 R(N,K)=R(N,K)+RR(N,K)
  495 CONTINUE
 5678 CONTINUE
      RETURN                                                            DRE02730
      END                                                               DRE02740
      SUBROUTINE SELF                                                   DRE02430
C***** CALCULATES LOADS AND ADDS TO RESIDUALS                           DRE02440
      include 'flat2.com'
      CALL AREAT
C****** DO FOR EACH ELEMENT                                             DRE05900
      DO 57 L=1,NE                                                      DRE05910
      WTT=AREA(L)*DEN/3.0D0
      N1=NOD(L,1)
      N2=NOD(L,2)
      N3=NOD(L,3)
      R(N1,3)=R(N1,3)-WTT
      R(N2,3)=R(N2,3)-WTT
      R(N3,3)=R(N3,3)-WTT
  57  CONTINUE
      RETURN                                                            DRE02730
      END                                                               DRE02740
      SUBROUTINE RESI                                                   DRE02750
      include 'flat2.com'
      DO 63 N=1,NN
      DO 63 K=1,NDIME
 63   R(N,K)=0.0D0
C****** REMEMBER G STRINGS MUST BE DONE FIRST
      IF (NSTR.GT.0)CALL GSRES
C****  LINKS IN HERE                                                    DRE02890
      IF(NL.GT.0)CALL LINKS 
C**** TRIANGULAR ELEMENTS IN HERE                                       DRE02950
      IF(NE.GT.0)CALL TRIANG
      RETURN                                                            DRE02980
      END                                                               DRE02990
      SUBROUTINE LENG                                                   DRE03000
C ******** CALCULATES THE INITIAL LINK LENGTHS                          DRE03010
      include 'flat2.com'
      IF(NL.EQ.0) GOTO 71                                               DRE03140
      DO 53 L=1,NL                                                      DRE03150
      IF(RLEN(L).LT.0.001) GOTO 31                                      DRE03160
C**** HER LINK LENGTH DIRECTLY FROM INPUT LIST                          DRE03170
      ODL(L)=RLEN(L)                                                    DRE03180
      GOTO 53                                                           DRE03190
  31  CONTINUE                                                          DRE03200
      MT=NODL(L,3)                                                      DRE03210
      NLT=MAT(MT)                                                       DRE03220
      IF(NLT.EQ.3)GOTO 761                                              DRE03230
C****** HERE LINK CALCULATED FROM INITIAL GEOMETRY                      DRE03240
      N1=NODL(L,1)                                                      DRE03250
      N2=NODL(L,2)                                                      DRE03260
      DL=0.0D0                                                          DRE03270
      DO 52 K=1,NDIME                                                   DRE03280
      D(K) = CORD(N2,K) -CORD(N1,K)                                     DRE03290
 52   DL=DL + D(K) ** 2                                                 DRE03300
      ODL(L) = DSQRT(DL)                                                DRE03310
      GOTO 53                                                           DRE03320
C ***** HERE LINK READ IN VALUE                                         DRE03330
 761  ODL(L) = AMAT(MT,4)                                               DRE03340
  53  CONTINUE                                                          DRE03350
  71  CONTINUE                                                          DRE03360
      IF(NE.EQ.0) GOTO 4455                                             DRE03370
C ***** TRIANGULAR ELEMENTS IN HERE                                     DRE03380
C    CALCULATE INITIAL LENGTHS                                          DRE03390
      DO 76 L=1,NE                                                      DRE03400
      MT=NOD(L,4)                                                       DRE03410
      NLT=MAT(MT)                                                       DRE03420
C **** CONSIDER EACH NODE                                               DRE03430
      NX=1                                                              DRE03440
      DO 772 NSS=1,3                                                    DRE03450
      NS=4-NSS                                                          DRE03460
      IF(TLEN(L,NS).LT.0.0001)GOTO 661                                  DRE03470
C*** HERE INITIAL LENGTH FROM LIST ON FILE                              DRE03480
      OLS(L,NS) = TLEN(L,NS)                                            DRE03490
      GOTO 771                                                          DRE03500
 661  CONTINUE                                                          DRE03510
      IF(NLT.NE.3) GOTO 66                                              DRE03520
C***** HERE INITIAL LENGTH FROM GROUP                                   DRE03530
      OLS(L,NS) =AMAT(MT,NS+2)                                          DRE03540
      GOTO 771                                                          DRE03550
  66  CONTINUE                                                          DRE03560
      N1=NOD(L,NS)                                                      DRE03570
      N2=NOD(L,NX)                                                      DRE03580
      DIF=0.0D0                                                         DRE03590
      DO 72 K=1,NDIME                                                   DRE03600
 72   DIF= DIF+(CORD(N2,K)-CORD(N1,K))**2                               DRE03610
      OLS(L,NS)=DSQRT(DIF)                                              DRE03620
 771  CONTINUE                                                          DRE03630
      NX=NS                                                             DRE03640
 772  CONTINUE                                                          DRE03650
C     WRITE(3,567)L,OLS(L,1),OLS(L,2),OLS(L,3)                          DRE03660
 567  FORMAT(' SLACK LENGTHS USED IN PROG',I5,3(5X,F10.2))              DRE03670
 76   CONTINUE                                                          DRE03680
 4455 CONTINUE                                                          DRE03690
      RETURN                                                            DRE03700
      END                                                               DRE03710
      SUBROUTINE PEAK                                                   DRE03720
      include 'flat2.com'
      WRITE(3,789)IT,ENG                                                DRE03850
      WRITE(*,789)IT,ENG                                                DRE03860
  789 FORMAT('  ENERGY AT PEAK   ITERATION',I5,'    EQUALS ',E20.8)     DRE03870
      ENG=0.0D0                                                         DRE03880
C     WRITE(*,45)                                                       DRE03890
C 45  FORMAT('   PEAK CALLED ')                                         DRE03900
      DO 26 N=1,NN                                                      DRE03910
      DO 26  K=1,NDIME                                                  DRE03920
      DASH= BB(N,K)                                                     DRE03930
C     WRITE(3,6789)N,K,BB(N,K),DASH                                     DRE03940
C6789 FORMAT(5X,2(I5,2X),2(F20.3,5X))                                   DRE03950
      CORD(N,K)=                                                        DRE03960
     =CORD(N,K)-(V(N,K)*DT*1.5)+((DT*DT*R(N,K))/(2.0D0*DASH))           DRE03970
      V(N,K)=0.0D0
  26  CONTINUE                                                          DRE03980
      IF(NE.LT.1) GOTO 518                                              DRE03990
      CALL FORM                                                         
  518 CONTINUE                                                          DRE04010
      IF(NSTR.EQ.0)GOTO 561
      CALL AVENOR
 561  CONTINUE
      IF (DEN.GT.0.0)CALL UDL
      CALL RESI                                                         DRE04030
      CALL LOAD                                                         DRE04040
      CALL MASS
      CALL INIT                                                         DRE04050
      RETURN                                                            DRE04060
      END                                                               DRE04070
      SUBROUTINE LINKS                                                  DRE04080
C ****** CALCULATES THE RESIDUAL TERMS FOR LINKS                        DRE04090
      include 'flat2.com'
      DO 52 L=1,NL                                                      DRE04220
      N1= NODL(L,1)                                                     DRE04230
      N2=NODL(L,2)                                                      DRE04240
      MT=NODL(L,3)                                                      DRE04250
C********* LINK OR CABLE TYPE                                           DRE04260
      NLT=MAT(MT)                                                       DRE04270
      DL=0.0D0                                                          DRE04280
C **************** CALCULATE LENGTHS                                    DRE04290
      DO 53 K=1,NDIME                                                   DRE04300
      D(K) =CORD(N2,K) -CORD(N1,K)                                      DRE04310
  53  DL=DL+ D(K)**2                                                    DRE04320
      DL=DSQRT(DL)                                                      DRE04330
      IF(NLT.EQ.1) GOTO 7631                                            DRE04340
      IF(NLT.EQ.2) GOTO 7641                                            DRE04350
      IF(NLT.EQ.5) GOTO 7888                                            DRE04360
C ******* NLT IS 0 OR 3 HERE                                            DRE04370
C ************   NOW A CABLE HERE                                       DRE04380
      E=AMAT(MT,2)                                                      DRE04390
      A= AMAT(MT,1)                                                     DRE04400
      TT(L) =(DL-ODL(L))*E*A/ODL(L)                                     DRE04410
C *************  ADD ON ANY PRESTRESS                                   DRE04420
      TT(L) = TT(L) + AMAT(MT,5)                                        DRE04430
      IF(TT(L).LT.0.0000000001) TT(L) = 0.0D0                           DRE04440
      GOTO 7456                                                         DRE04450
 7641 CONTINUE                                                          DRE04460
C*********** HERE CABLE OF FIXED TENSION                                DRE04470
      TT(L) = AMAT(MT,5)                                                DRE04480
      GOTO 7456                                                         DRE04490
 7888 CONTINUE                                                          DRE04500
C***** HERE A BAR OF CONSTANT FORCE DENSITY                             DRE04510
      TT(L)=AMAT(MT,5)*DL                                               DRE04520
      GOTO 7456                                                         DRE04530
 7631 CONTINUE                                                          DRE04540
C ************** HERE A BAR                                             DRE04550
      E=AMAT(MT,2)                                                      DRE04560
      A= AMAT(MT,1)                                                     DRE04570
      TT(L) =(DL-ODL(L))*E*A/ODL(L)                                     DRE04580
 7456 CONTINUE                                                          DRE04590
C ********  ADD TO RES VECTOR                                           DRE04600
C     WRITE (3,4561)L,TT(L)                                             DRE04610
 4561 FORMAT (' FROM LINKS     NO=',I5,'    TENSION',F20.5)             DRE04620
      DO 102 K=1,NDIME                                                  DRE04630
      DR=D(K) *TT(L)/DL                                                 DRE04640
      R(N1,K)=R(N1,K)+DR                                                DRE04650
 102  R(N2,K)=R(N2,K)-DR                                                DRE04660
 52   CONTINUE                                                          DRE04670
 716  CONTINUE                                                          DRE04680
      RETURN                                                            DRE04690
      END                                                               DRE04700
      SUBROUTINE TRIANG                                                 DRE04710
C***** CALCULATES RESIDUALS FOR TRIANGLES                               DRE04720
      include 'flat2.com'
c     WRITE(6,*)'TRINAG CALLED'
      CALL AREAT
      DO 57 L=1,NE                                                      DRE04850
      MT=NOD(L,4)                                                       DRE04860
C     WRITE(3,897)L,NE                                                  DRE04870
C 897 FORMAT(' TRIA CHECK INITIAL AND C',I5,5X,I5)                      DRE04880
      NLT=MAT(MT)                                                       DRE04890
      PRE=AMAT(MT,4)
      PRE2=AMAT(MT,5)
      THIC=AMAT(MT,3)
C     WRITE(3,567)L,OLS(L,1),OLS(L,2),OLS(L,3)                          DRE03660
 567  FORMAT(' ORIG SLACK LENGTHS USED ',I5,3(5X,F10.2))    
C     WRITE(3,5672)L,DTL(L,1),DTL(L,2),DTL(L,3)      
 5672 FORMAT(' TRIANG- LENT OF TRI SIDES',I5,3(5X,F10.2))               DRE05090
C***** CALCULATE TRIANGLE ANGLES                                        DRE05100
      SUM=DTL(L,1)**2 + DTL(L,2)**2 + DTL(L,3)**2    
      ANG(2)=DACOS((SUM-2.0D0*DTL(L,2)**2)/(2.0D0*DTL(L,1)*DTL(L,3))) 
      ANG(3)=DACOS((SUM-2.0D0*DTL(L,3)**2)/(2.0D0*DTL(L,1)*DTL(L,2)))
      ANG(1)=DACOS((SUM-2.0D0*DTL(L,1)**2)/(2.0D0*DTL(L,2)*DTL(L,3))) 
C     WRITE(3,9911)ANG(1),ANG(2),ANG(3)
 9911 FORMAT ('    ANGLES',3F10.5)
      IF(NLT.EQ.2) GOTO 712                                             DRE05150
      GOTO 812                                                          DRE05160
C****  NOW SPECIFIC TYPES                                               DRE05170
 712  CONTINUE                                                          DRE05180
C****** HERE MEMBRANE IS OF FIXED TENSION                               DRE05190
C**** FORMFINDING     CONSIDER EACH SIDE                                DRE05200
      NX=1                                                              DRE05210
      DO 76 NSS=1,3                                                     DRE05220
      NS=4-NSS                                                          DRE05230
      T=PRE*THIC*DTL(L,NS)/(2.0D0*DTAN(ANG(NS))) 
c      IF (NS.EQ.1) T=T+(AREA(L)*(PRE2-PRE)/DTL(L,1))
C     WRITE(6,56)T
   56 FORMAT (' T=',E15.7)
      N1= NOD(L,NS)                                                     DRE05250
      N2= NOD(L,NX)                                                     DRE05260
      NX= NS                                                            DRE05270
C***** EACH DIMENSION                                                   DRE05280
      DO 78 K=1,NDIME                                                   DRE05290
      DR = DTT(L,NS,K)*T/DTL(L,NS) 
C     WRITE(3,5671)N1,N2,K,T,DR                                         DRE05310
 5671 FORMAT(' RES CALC' ,3(4X,I5),3(F20.5))                            DRE05320
      R(N1,K)=R(N1,K) + DR                                              DRE05330
  78  R(N2,K)=R(N2,K)-DR                                                DRE05340
  76  CONTINUE                                                          DRE05350
      GOTO 57                                                           DRE05360
 812  CONTINUE                                                          DRE05370
C***** HERE TRIANGLE                                                    DRE05380
C**** OK ALREADY CALCULATED                                             DRE05390
C**** CALCULATED EXTENSION FOR EACH SIDE                                DRE05400
      DO 27 NS=1,3                                                      DRE05410
  27  EX(NS) =-OLS(L,NS) +DTL(L,NS)
C     WRITE(3,566)L,OLS(L,1),OLS(L,2),OLS(L,3)                          DRE03660
 566  FORMAT(' SLACK LENGTHS USED IN PROG',I5,3(5X,F10.2))              DRE03670
C      WRITE(3,564)L,DTL(L,1),DTL(L,2),DTL(L,3) 
 564  FORMAT(' NEW LENGTHS USED IN PROG',I5,3(5X,F10.2))   
C      WRITE(3,4111)L,EX(1),EX(2),EX(3)                 
 4111 FORMAT(' TRIA - SIDE EXT',I5,3(1X,F15.5)) 
C***** MULTIPY EX BY OK TO GET EDGE FORCES                              DRE05450
C******* DO FOR EACH SIDE                                               DRE05460
      NX = 1                                                            DRE05470
      DO 271  NSS=1,3                                                   DRE05480
      T=0.0D0                                                           DRE05490
      NS=4-NSS                                                          DRE05500
      DO 291  II= 1,3                                                   DRE05510
 291  T=OK(L,NS,II)*EX(II) + T                                          DRE05520
C***** GET OUT NOS OF SIDES                                             DRE05530
      N1= NOD(L,NS)                                                     DRE05540
      N2= NOD(L,NX)                                                     DRE05550
C***** HERE PRETENSION MUST BE CONSIDERED                               DRE05560
C***** ADD INTO RESIDUALS                                               DRE05600
      DO 781 K=1,NDIME                                                  DRE05610
      DR =DTT(L,NS,K) * T/DTL(L,NS)
C     WRITE(3,45) DR,T, DTL(L,NS),DTT(L,NS,K)
C45   FORMAT(' IN TRI RES COMP',4 (2X,F15.6))                           DRE05640
      R(N1,K) = R(N1,K)+ DR                                             DRE05650
      R(N2,K) = R(N2,K) -DR                                             DRE05660
  781 CONTINUE                                                          DRE05670
      NX=NS                                                             DRE05680
  271 CONTINUE                                                          DRE05690
  57  CONTINUE                                                          DRE05700
  516 CONTINUE                                                          DRE05710
      RETURN                                                            DRE05720
      END                                                               DRE05730
      SUBROUTINE FORM                                                   DRE05740
C****** FORMS STIFFNESSES FOR THE TRIANGLES                             DRE05750
      include 'flat2.com'
      DIMENSION AX(3),BX(3),CX(3),ANAX(3)                               DRE05880
      DIMENSION        RINT(3,3)                                        DRE05890
C ***** DO FOR EACH ELEMENT                                             DRE05900
      DO 57 L=1,NE                                                      DRE05910
C***** INITIAL LENGTHS ALREADY CALCULATED IN LENGTH                
      WRITE(3,567)L,DTL(L,1),DTL(L,2),DTL(L,3) 
 567  FORMAT(' FORM LEN OF TRI SIDES',I5,3(5X,F10.2))        
C***** HERE USE CURRENT LENGTHS TO CALCULATE THE ANGLES OF TRIANGLE     DRE06080
      SUM=DTL(L,1)**2 + DTL(L,2)**2 + DTL(L,3)**2 
      ANG(2)=DACOS((SUM-2.0D0*DTL(L,2)**2)/(2.0D0*DTL(L,1)*DTL(L,3))) 
      ANG(3)=DACOS((SUM-2.0D0*DTL(L,3)**2)/(2.0D0*DTL(L,1)*DTL(L,2)))  
      ANG(1)=DACOS((SUM-2.0D0*DTL(L,1)**2)/(2.0D0*DTL(L,2)*DTL(L,3)))  
      WRITE(3,9999) L,ANG(1),ANG(2),ANG(3)                              DRE06130
 9999 FORMAT(' FORM ANGLES',I5,3(5X,F10.5))                             DRE06140
C**** CALCULATE TERMS                                                   DRE06150
C***** TO CALCULATE THE INITIAL ANGLES                                  DRE06160
      PIE=DATAN(1.0D0)*4.0D0                                            DRE06170
      ANAX(1) = 0.0D0                                                   DRE06180
      ANAX(2) = PIE -ANG(3)                                             DRE06190
      ANAX(3) =PIE+     ANG(2)                                        
      WRITE(3,7676) ANAX(1),ANAX(2),ANAX(3)                             DRE06210
 7676 FORMAT (' ANG 1,2,3 :',3(5X,F20.5))                               DRE06220
      DO 7111 J=1,3                                                     DRE06230
      AX(J) =  DCOS(ANAX(J)) **2                                        DRE06240
      BX(J) = DSIN(ANAX(J)) **2                                         DRE06250
      CX(J) = DSIN(ANAX(J))*DCOS(ANAX(J))                               DRE06260
      WRITE(3,5676) AX(1),BX(2),CX(3)                                   DRE06270
 5676 FORMAT (' A,B,C X :',3(5X,F20.5))                                 DRE06280
 7111 CONTINUE                                                          DRE06290
      BOT =AX(1) *(BX(2) *CX(3) -CX(2)*BX(3))                           DRE06300
     1     -BX(1)*(AX(2)*CX(3)-CX(2)*AX(3))                             DRE06310
     2     +CX(1)*(AX(2)*BX(3)-BX(2)*AX(3))                             DRE06320
      WRITE (3,4421) BOT                                                DRE06330
 4421 FORMAT( '    BOT=',F20.5)                                         DRE06340
      G(1,1)=((BX(2)*CX(3))-(BX(3)*CX(2)))/(BOT*OLS(L,1))               DRE06350
      G(1,2)=((BX(3)*CX(1))-(BX(1)*CX(3)))/(BOT*OLS(L,2))               DRE06360
      G(1,3)=((BX(1)*CX(2))-(BX(2)*CX(1)))/(BOT*OLS(L,3))               DRE06370
      G(2,1)=((AX(3)*CX(2))-(AX(2)*CX(3)))/(BOT*OLS(L,1))               DRE06380
      G(2,2)=((AX(1)*CX(3))-(AX(3)*CX(1)))/(BOT*OLS(L,2))               DRE06390
      G(2,3)=((AX(2)*CX(1))-(AX(1)*CX(2)))/(BOT*OLS(L,3))               DRE06400
      G(3,1)=((AX(2)*BX(3))-(AX(3)*BX(2)))/(BOT*OLS(L,1))               DRE06410
      G(3,2)=((AX(3)*BX(1))-(AX(1)*BX(3)))/(BOT*OLS(L,2))               DRE06420
      G(3,3)=((AX(1)*BX(2))-(AX(2)*BX(1)))/(BOT*OLS(L,3))               DRE06430
      DO 7933  J=1,3                                                    DRE06440
 7933    WRITE(3,7322)J,(G(J,I),I=1,3)                                  DRE06450
 7322  FORMAT(' ROW OF G ',I5,3(3X,F20.1))                              DRE06460
C****** FORM THE D MATRIX                                               DRE06470
      MT= NOD(L,4)                                                      DRE06480
      THIC= AMAT(MT,3)
      E=AMAT(MT,2)                                                      DRE06500
      VV=AMAT(MT,1)                                                     DRE06510
      COM= E/(1.0D0-VV*VV)                                              DRE06520
      DD(1,1)=COM                                                       DRE06530
      DD(1,2)=COM*VV                                                    DRE06540
      DD(1,3)=0.0D0                                                     DRE06550
      DD(2,1)=COM*VV                                                    DRE06560
      DD(2,2)=COM                                                       DRE06570
      DD(2,3)=0.0D0                                                     DRE06580
      DD(3,1)=0.0D0                                                     DRE06590
      DD(3,2)=0.0D0                                                     DRE06600
      DD(3,3)=COM*(1.0D0-VV)/2.0D0                                      DRE06610
      DO 933  J=1,3                                                     DRE06620
  933    WRITE(3,1322)J,(DD(J,I),I=1,3)                                 DRE06630
 1322  FORMAT(' ROW OF D ',I5,3(3X,F20.1))                              DRE06640
C**** CHECK BUCKLING                                                    DRE06650
      IF(MAT(MT).GT.0) GOTO 5688                                        DRE06670
      CALL BUCK(L,E)                                                    DRE06680
 5688 CONTINUE                                                          DRE06690
      DO 60 J=1,3                                                       DRE06700
      DO 60 I=1,3                                                       DRE06710
      RINT(J,I) = 0.0D0                                                 DRE06720
      DO 60 K=1,3                                                       DRE06730
   60 RINT(J,I)= RINT(J,I) +DD(J,K) *G(K,I)                             DRE06740
C**** CALCULATE AREA FROM NEW LENGHTS S FORMULA ZIENK                   DRE06750
      DO 61 J=1,3                                                       DRE06780
      DO 61 I=1,3                                                       DRE06790
      OK(L,J,I) =0.0D0                                                  DRE06800
      DO 61 K=1,3                                                       DRE06810
  61  OK(L,J,I) = OK(L,J,I) +G(K,J) *RINT(K,I)*AREA(L)*THIC             DRE06820
      WRITE(3,422)L,AREA(L)                                             DRE06830
 422  FORMAT('  FORM ELEMENT OK ',I4,5X,'AREA',F20.6)                   DRE06840
      DO 233  J=1,3                                                     DRE06850
  233    WRITE(3,322)J,(OK(L,J,I),I=1,3)                                DRE06860
 322  FORMAT(' ROW ',I5,3(3X,F20.1))                                    DRE06870
  57  CONTINUE                                                          DRE06880
      RETURN                                                            DRE06890
      END                                                               DRE06900
      SUBROUTINE BUCK(L,E)                                              DRE06910
C***** CALCULATES CURRENT ANGLES IN THE TRIANGLES                       DRE06920
      include 'flat2.com'
      DIMENSION VS(3),ST(3),RINT(3,3)                                   DRE07050
C**** CALCULATES ANGLE OF MAXIMUM STRESS                                DRE07060
C**** FIRST CALCULATE STRESSES                                          DRE07070
C**** CALCULATE D.G.DEL                                                 DRE07080
C     USE ORIGINAL D                                                    DRE07090
C*** CALCULATE DEL                                                      DRE07100
      DO 27 NS = 1,3                                                    DRE07110
 27   EX(NS) = OLS(L,NS)-DTL(L,NS)
      DO 28 J=1,3                                                       DRE07130
      VS(NS) = 0.0D0                                                    DRE07140
      DO 28 I=1,3                                                       DRE07150
  28  VS(J)=VS(J) +G(J,I) *EX(I)                                        DRE07160
      DO 29 J=1,3                                                       DRE07170
      ST(J) =0.0D0                                                      DRE07180
      DO 29 I=1,3                                                       DRE07190
  29  ST(J) = ST(J)+DD(J,I) *VS(J)                                      DRE07200
      CCX=(ST(1)+ST(2))/2.0D0                                           DRE07210
      BBX=(ST(1)-ST(2))/2.0D0                                           DRE07220
      AAX=DSQRT((BBX**2)+(ST(3)**2))                                    DRE07230
      SMAX=CCX+AAX                                                      DRE07240
      SMIN =CCX-AAX                                                     DRE07250
      DIF=ST(2) -SMIN                                                   DRE07260
      IF(DABS(DIF).GT.0.0000001) GOTO 8923                              DRE07270
      BANG=0.0D0                                                        DRE07280
      GOTO 567                                                          DRE07290
 8923 CONTINUE                                                          DRE07300
C***** HERE ELEMENT HAS ANGLE NE 0.0                                    DRE07310
      BANG = DATAN(ST(3)/DIF)                                           DRE07320
      STRPER =0.0D0                                                     DRE07330
C***** CHECK STRESS                                                     DRE07340
      NBCX= 1                                                           DRE07350
      IF(SMAX.GT.STRPER) GOTO 831                                       DRE07360
      NBCX=0                                                            DRE07370
      DD(3,3) =0.0D0                                                    DRE07380
      DD(1,1) =0.0D0                                                    DRE07390
      DD(1,2) = 0.0D0                                                   DRE07400
      DD(2,1)= 0.0D0                                                    DRE07410
      DD(2,2)=E                                                         DRE07420
 831  IF(SMIN.GT.STRPER) GOTO 859                                       DRE07430
      NBCX=0                                                            DRE07440
      DD(3,3) =0.0D0                                                    DRE07450
      DD(2,1) =0.0D0                                                    DRE07460
      DD(2,2) =0.0D0                                                    DRE07470
      DD(1,2) =0.0D0                                                    DRE07480
      IF(SMAX.GT.STRPER)DD(1,1)=E                                       DRE07490
 859  IF(NBCX.EQ.1) GOTO 567                                            DRE07500
C**** HERE CALCULATE T                                                  DRE07510
      CAN=DCOS(BANG)                                                    DRE07520
      SAN=DSIN(BANG)                                                    DRE07530
      TR(1,1)=CAN*CAN                                                   DRE07540
      TR(1,2)=SAN*SAN                                                   DRE07550
      TR(1,3)=-2.0D0*SAN*CAN                                            DRE07560
      TR(2,1)=SAN*SAN                                                   DRE07570
      TR(2,2)=CAN*CAN                                                   DRE07580
      TR(3,1)= CAN*SAN                                                  DRE07590
      TR(3,2)= -SAN*CAN                                                 DRE07600
      TR(3,3)= CAN*CAN-SAN*SAN                                          DRE07610
C***** HERE TRANSFORM DD                                                DRE07620
      DO 55 J=1,3                                                       DRE07630
      DO 55 I=1,3                                                       DRE07640
      RINT(J,I)=0.0D0                                                   DRE07650
      DO 55 K=1,3                                                       DRE07660
 55   RINT(J,I) = RINT(J,I)+DD(J,K)*TR(I,K)                             DRE07670
      DO 58 J=1,3                                                       DRE07680
      DO 58 I=1,3                                                       DRE07690
      DD(J,I)=0.0D0                                                     DRE07700
      DO 58 K=1,3                                                       DRE07710
 58   DD(J,I) = DD(J,I)+TR(J,K)*RINT(K,I)                               DRE07720
C**** NOW DD CONTAINS TERMS RELATIVE TO THE LOCAL AXES FOR G            DRE07730
      RETURN                                                            DRE07740
 567  CONTINUE                                                          DRE07750
C****** HERE LOCAL AXES X' IS ALONG SIDE ONE OF THE ELEMENT             DRE07760
      RETURN                                                            DRE07770
      END                                                               DRE07780
      SUBROUTINE RSTART                                                 DRE07790
      include 'flat2.com'
      DIMENSION AX(3),BX(3),CX(3),ANAX(3)
      DIMENSION VS(3,3),ST(3),RINT(3,3)                                 DRE07930
C****  DATA INPUT ROUTINE                                               DRE07940
      WRITE(4,678)TITEL                                                 DRE07950
 678  FORMAT(A32)                                                       DRE07960
      WRITE(3,1)TITEL                                                   DRE07970
 1    FORMAT(      /   /,' DATA AND RESULTS  ',/,1X,A80)                DRE07980
      WRITE(4,345)DT,DAMP
      WRITE(3,345)DT,DAMP
 345  FORMAT(2F10.5)
      WRITE(4,2)NN,NE,NB,NM,NL,NDIME,NSTR,NBE,KMAX,NMAX,ILN             DRE07990
      write(4,2) nq
 2    FORMAT(BZ,2I6,9I5)                                                   DRE08000
      WRITE(3,3)NN,NE,NB,NM,NL,NDIME,NSTR,NBE,KMAX,NMAX           
 3    FORMAT(/,' NO OF NODES=',I5,/,' NO OF ELEMENTS=',I5,/,            DRE08020
     1' NO OF RESTRAINED BOUNDARY NODES=',I5,/,' NO OF MATERIAL TYPES=',DRE08030
     2I5,/,' NO OF LINKS=',I5,/,' NO OF DIMENSIONS=',I5,/,' NO OF STRS='DRE08040
     3,I5,/,' NO BEAMS=',I5,/,' NO OF ITGRPS=',I5,/,' NO OF ITBTWP=',I5)DRE08050
      WRITE(3,267)ILN                                                   DRE08060
 267  FORMAT(/,' O IF NO LOADS=',I5)                                    DRE08070
      write(3,268) nq
 268  format(' ','NO OF QUADS =',1X,I5) 
      DO 44 N=1,NN                                                      DRE08080
  44  WRITE(4,4) N,(CORD(N,I),I=1,NDIME)                                DRE08090
 4    FORMAT(BZ,I5,3F15.5)                                              DRE08100
      DO 213 N=1,NN                                                     DRE08110
  213 WRITE(3,5) N,(CORD(N,I),I=1,NDIME)                                DRE08120
 5    FORMAT(' NODE=',I5, 5X,F15.10,2(F15.10))                          DRE08130
      IF(NE.EQ.0) GOTO 56                                               DRE08140
       DO 2157 L=1,NE
C***** CALCULATE SLACK LENGTHS OF TRIANGLE SIDES IF FF
      MT = NOD(L,4)
      NLT = MAT(MT)
      IF (NLT.NE.2) GOTO 7223
      EE = AMAT(MT,2)
      VV = AMAT (MT,1)
      THIC = AMAT(MT,3)
      PRE= AMAT(MT,4)
      FACT=EE/(EE+(PRE*(1.0D0-VV)))
      NX=1
      DO 7177 NSS=1,3
      NS= 4-NSS
      N1=NOD(L,NS)
      N2=NOD(L,NX)
      NX=NS
      DTL(L,NS)=0.0D0
      DO 572 K=1,NDIME
      DTT(L,NS,K)=CORD(N2,K)-CORD(N1,K)+10.0D-25
      DTL(L,NS) =DTL(L,NS) + DTT(L,NS,K)**2
 572  CONTINUE
      DTL(L,NS)=DSQRT(DTL(L,NS))
      OLS(L,NS)=DTL(L,NS)*FACT
 7177   CONTINUE
 7223  WRITE(4,6)L,(NOD(L,I),I=1,4),(OLS(L,I),I=1,3) 
 6    FORMAT(BZ,5I5,3F15.8)                                             DRE08160
 2157 CONTINUE
      DO 7777 N=1,NE                                                    DRE08170
 7777 WRITE(3,7) N,(NOD(N,I),I=1,4),(OLS(N,I),I=1,3)                    DRE08180
 7    FORMAT(' ELEMENT=',I3,2X,'NODE1=',I3,2X,'NODE2=',I3,2X,'NODE3='   DRE08190
     1,I3,2X,'MAT NO=',I5,'INITIAL LENGTHS:',3(3X,F15.8))               DRE08200
  56  IF(NL.EQ.0) GOTO 57                                               DRE08210
      WRITE(4,76)(N,(NODL(N,I),I=1,3),ODL(N),N=1,NL)                    DRE08220
 76   FORMAT(BZ,4I5,F15.8)                                              DRE08230
      WRITE(3,77)(N,(NODL(N,I),I=1,3),ODL(N),N=1,NL)                    DRE08240
 77   FORMAT(' LINK=',I3,2X,'NODE1=',I3,2X,'NODE2=',I3,2X,              DRE08250
     1'MAT NO=',I5,'  INITIAL LENGTH=',F15.8)                           DRE08260
 57   CONTINUE                                                          DRE08270
      DO 577 N=1,NM                                                     DRE08280
 577  WRITE(4,11)N,MAT(N),(AMAT(N,I),I=1,5)                             DRE08290
 11   FORMAT(BZ,2I5,5e12.5) 
      WRITE(3,12)                                                       DRE08310
 12   FORMAT(1H0,'MATERIAL DATA')                                       DRE08320
      DO 1313 N=1,NM                                                    DRE08330
 1313 WRITE(3,13) N,MAT(N),(AMAT(N,I),I=1,5)                            DRE08340
 13   FORMAT('  MATERIAL NUMBER',I2,3X,I5,3X,5f12.5)
      DO 88 N=1,NB                                                      DRE08360
 88   WRITE(4,8)(NCO(N,I),I=1,(NDIME+1))                                DRE08370
  8   FORMAT(BZ,4I5)                                                    DRE08380
      WRITE(3,9)                                                        DRE08390
 9    FORMAT(/,' BOUNDARY CONDITIONS')                                  DRE08400
      DO 100 N=1,NB                                                     DRE08410
 100  WRITE(3,10)(NCO(N,I),I=1,(NDIME+1))                               DRE08420
 10   FORMAT(/,' NODE=',I5,5X,3(5X,I2))                                 DRE08430
      WRITE(3,999)                                                      DRE08440
 999  FORMAT(///, ' LOADING CONDITIONS'/)                               DRE08450
      IF(ILN.EQ.0)GOTO 8989                                             DRE08460
      DO 567 I=1,ILN                                                    DRE08470
      WRITE(4,16)NC(I),(ALOAD(I,II),II=1,NDIME)                         DRE08480
 16   FORMAT(I5,3F15.5)                                                 DRE08490
      WRITE(3,17)NC(I),(ALOAD(I,II),II=1,NDIME)                         DRE08500
 17   FORMAT(1H ,'NODE=',I5,5X,3(4X,F15.5))                             DRE08510
 567  CONTINUE                                                          DRE08520
8989  CONTINUE
      IF (NSTR.EQ.0)GOTO 8844
C**** WRITES STRING DATA                                             
      WRITE(3,234)NSTR
 234  FORMAT(/,'NO OF STRINGS=',I5)
      DO 58 J=1,NSTR                                                    
      WRITE(4,781)NGS(J)                                                
  781 FORMAT(I5)                                                        
      WRITE(3,782)J,NGS(J)                                             
  782 FORMAT(' STRING NO=',I5,'   NO OF LINKS=',I5,/,                   
     .' STRING NODE:')                                                 
      WRITE(4,784)(NSNOD(J,JD),JD=1,(NGS(J)+1))                        
  784 FORMAT (8I5)                                                     
      WRITE(3,785) (NSNOD(J,JD),JD=1,(NGS(J)+1))                      
  785 FORMAT(1X,8I5)
   58 CONTINUE                                                          
 8844 CONTINUE
      DO 94 N=1,NN                                                     
  94  WRITE(3,466) N,(CORD(N,I),I=1,NDIME)                              DRE08540
 466  FORMAT(BZ,I5,3F15.5)                                              DRE08550
      DO 266 N=1,NN                                                     DRE08560
 266  WRITE(3,233) N,(R(N,I),I=1,NDIME)                                 DRE08570
 233  FORMAT('  NODE = ',I5,5X,'RESIDUALS:',3(4X,E15.7))                DRE08580
      DO 6000 N=1,NN                                                    DRE08590
 6000 WRITE(3,633)N,(V(N,I),I=1,NDIME)                                  DRE08600
 633  FORMAT('  NODE = ',I5,5X,'VELOCITIES:',3(4X,F10.5))               DRE08610
      DO 899 N=1,NN                                                     DRE08620
 899  WRITE(3,833) N,(BB(N,I),I=1,NDIME)                                DRE08630
 833  FORMAT('  NODE = ',I5,5X,'MASSES:',3(4X,E15.7))                   DRE08640
      IF (NL.LT.1) GOTO 6711                                            DRE08650
      DO 52 L=1,NL                                                      DRE08660
      N1= NODL(L,1)                                                     DRE08670
      N2=NODL(L,2)                                                      DRE08680
      MT=NODL(L,3)                                                      DRE08690
C********* LINK OR CABLE TYPE                                           DRE08700
      NLT=MAT(MT)                                                       DRE08710
      DL=0.0D0                                                          DRE08720
C **************** CALCULATE LENGTHS                                    DRE08730
      DO 53 K=1,NDIME                                                   DRE08740
      D(K) =CORD(N2,K) -CORD(N1,K)                                      DRE08750
  53  DL=DL+ D(K)**2                                                    DRE08760
      DL=DSQRT(DL)                                                      DRE08770
      IF(NLT.EQ.1) GOTO 7631                                            DRE08780
      IF(NLT.EQ.2) GOTO 7641                                            DRE08790
      IF(NLT.EQ.5) GOTO 5557                                            DRE08800
C ******* NLT IS 0 OR 3 HERE                                            DRE08810
C ************   NOW A CABLE HERE                                       DRE08820
      E=AMAT(MT,2)                                                      DRE08830
      A= AMAT(MT,1)                                                     DRE08840
      TT(L) =(DL-ODL(L))*E*A/ODL(L)                                     DRE08850
C *************  ADD ON ANY PRESTRESS                                   DRE08860
      TT(L) = TT(L) + AMAT(MT,5)                                        DRE08870
      IF(TT(L).LT.0.0000000001) TT(L) = 0.0D0                           DRE08880
      GOTO 7456                                                         DRE08890
 7641 CONTINUE                                                          DRE08900
C*********** HERE CABLE OF FIXED TENSION                                DRE08910
      TT(L) = AMAT(MT,5)                                                DRE08920
      GOTO 7456                                                         DRE08930
 5557 CONTINUE                                                          DRE08940
C**** HERE CONSTANT FORCE DENSITY                                       DRE08950
      TT(L) =AMAT(MT,5)*DL                                              DRE08960
      GOTO 7456                                                         DRE08970
 7631 CONTINUE                                                          DRE08980
C ************** HERE A BAR                                             DRE08990
      E=AMAT(MT,2)                                                      DRE09000
      A= AMAT(MT,1)                                                     DRE09010
      TT(L) =(DL-ODL(L))*E*A/ODL(L)                                     DRE09020
 7456 CONTINUE                                                          DRE09030
C *CAL FORCE NOW OUTPUT                                                 DRE09040
C      WRITE (3,4561)L,N1,N2,DL,TT(L)                                   
 4561 FORMAT(' LINK=',I5,' NODES:',2(2X,I5),'  CLEN=',F15.9,' TENSION=',DRE09060
     . F15.9)                                                           DRE09070
 52   CONTINUE                                                          DRE09080
 716  CONTINUE                                                          DRE09090
 6711 CONTINUE                                                          DRE09100
      IF(NE.LT.1) GOTO 3311                                             DRE09110
      CALL FORM
      DO 157 L=1,NE                                                     DRE09120
      MT=NOD(L,4)                                                       DRE09130
C     WRITE(3,897)L,NE                                                  DRE09140
C 897 FORMAT(' TRIA CHECK INITIAL AND C',I5,5X,I5)                      DRE09150
      NLT=MAT(MT)                                                       DRE09160
      PRE=AMAT(MT,4)                                                    DRE09170
      THIC=AMAT(MT,3)
C**** CONSIDER EACH NODE                                                DRE09190
      NX=1                                                              DRE09200
      DO 71 NSS=1,3                                                     DRE09210
      NS=4-NSS                                                          DRE09220
      N1=NOD(L,NS)                                                      DRE09230
      N2=NOD(L,NX)                                                      DRE09240
      NX=NS                                                             DRE09250
      DTL(L,NS)=0.0D0
      DO 72 K=1,NDIME                                                   DRE09270
      DTT(L,NS,K) =CORD(N2,K)-CORD(N1,K)+10.D-25                          DRE09280
      DTL(L,NS) =DTL(L,NS) + DTT(L,NS,K) **2 
C     WRITE(3,456) N1,N2,K,DTT(L,NS,K),DTL(L,NS) 
C456  FORMAT(' TRIAN W34',3(2X,I4),2(5X,F10.5))                         DRE09310
  72  CONTINUE                                                          DRE09320
      DTL(L,NS) = DSQRT(DTL(L,NS))  
  71  CONTINUE                                                          DRE09340
C     WRITE(3,567)L,DTL(1L,),DTL(L,2),DTL(L,3)   
C567  FORMAT(' TRIANG- LENT OF TRI SIDES',I5,3(5X,F10.2))               DRE09360
C***** CALCULATE TRIANGLE ANGLES                                        DRE09370
      SUM=DTL(L,1)**2 + DTL(L,2)**2 + DTL(L,3)**2               
      ANG(2)=DACOS((SUM-2.0D0*DTL(L,2)**2)/(2.0D0*DTL(L,1)*DTL(L,3)))   
      ANG(3)=DACOS((SUM-2.0D0*DTL(L,3)**2)/(2.0D0*DTL(L,1)*DTL(L,2)))   
      ANG(1)=DACOS((SUM-2.0D0*DTL(L,1)**2)/(2.0D0*DTL(L,2)*DTL(L,3)))    
C**** CALCULATED EXTENSION FOR EACH SIDE                                DRE09420
      DO 27 NS=1,3                                                      DRE09430
  27  EX(NS) =-OLS(L,NS) +DTL(L,NS)                                       DRE09440
C     WRITE(3,4111)L,EX(1),EX(2),EX(3)                                  DRE09450
C4111 FORMAT(' TRIA - SIDE EXTENSIONS',I5,3(5X,F20.5))                  DRE09460
C***** MULTIPY EX BY OK TO GET EDGE FORCES                              DRE09470
C******* DO FOR EACH SIDE                                               DRE09480
C      DO 1233  J=1,3     
C 1233    WRITE(3,122)J,(OK(L,J,I),I=1,3) 
 122  FORMAT(' ROW ',I5,3(3X,F20.10))                                   DRE09510
      NX = 1                                                            DRE09520
      DO 271  NSS=1,3                                                   DRE09530
      T=0.0D0                                                           DRE09540
      NS=4-NSS                                                          DRE09550
      DO 291  II= 1,3                                                   DRE09560
 291  T=OK(L,NS,II)*EX(II) + T                                          DRE09570
C***** GET OUT NOS OF SIDES                                             DRE09580
      N1= NOD(L,NS)                                                     DRE09590
      N2= NOD(L,NX)                                                     DRE09600
C***** HERE PRETENSION MUST BE CONSIDERED                               DRE09610
      WRITE(3,45)L,N1,N2,T, DTL(L,NS),OLS(L ,NS),EX(NS) 
 45   FORMAT(' ELE=',I5,' N1=',I5,' N2=',I5,'  TENSION=',F10.4,         DRE09630
     1'  EX LEN=',F10.4,' ORIG LEN=',F10.4,' EXT=',F15.9)               DRE09640
      NX=NS                                                             DRE09650
  271 CONTINUE                                                          DRE09660
C***** HERE USE CURRENT LENGTHS TO CALCULATE THE ANGLES OF TRIANGLE     DRE09670
      SUM=DTL(L,1)**2 + DTL(L,2)**2 + DTL(L,3)**2    
      ANG(2)=DACOS((SUM-2.0D0*DTL(L,2)**2)/(2.0D0*DTL(L,1)*DTL(L,3)))   
      ANG(3)=DACOS((SUM-2.0D0*DTL(L,3)**2)/(2.0D0*DTL(L,1)*DTL(L,2)))
      ANG(1)=DACOS((SUM-2.0D0*DTL(L,1)**2)/(2.0D0*DTL(L,2)*DTL(L,3)))
C     WRITE(3,9999) L,ANG(1),ANG(2),ANG(3)                              DRE09720
C9999 FORMAT(' FORM ANGLES',I5,3(5X,F10.5))                             DRE09730
C**** CALCULATE TERMS                                                   DRE09740
C***** TO CALCULATE THE INITIAL ANGLES                                  DRE09750
      PIE=DATAN(1.0D0)*4.0D0                                            DRE09760
      ANAX(1) = 0.0D0                                                   DRE09770
      ANAX(2) = PIE -ANG(3)                                             DRE09780
      ANAX(3) =     PIE+ANG(2)                                          DRE09790
C     WRITE(3,7676) ANAX(1),ANAX(2),ANAX(3)                             DRE09800
 7676 FORMAT (' ANG 1,2,3 :',3(5X,F20.5))                               DRE09810
      DO 7111 J=1,3                                                     DRE09820
      AX(J) =  DCOS(ANAX(J)) **2                                        DRE09830
      BX(J) = DSIN(ANAX(J)) **2                                         DRE09840
      CX(J) = DSIN(ANAX(J))*DCOS(ANAX(J))                               DRE09850
C     WRITE(3,5676) AX(1),BX(2),CX(3)                                   DRE09860
 5676 FORMAT (' A,B,C X :',3(5X,F20.5))                                 DRE09870
 7111 CONTINUE                                                          DRE09880
      BOT =AX(1) *(BX(2) *CX(3) -CX(2)*BX(3))                           DRE09890
     1     -BX(1)*(AX(2)*CX(3)-CX(2)*AX(3))                             DRE09900
     2     +CX(1)*(AX(2)*BX(3)-BX(2)*AX(3))                             DRE09910
C     WRITE (3,4421) BOT                                                DRE09920
 4421 FORMAT( '    BOT=',F20.5)                                         DRE09930
      G(1,1)=((BX(2)*CX(3))-(BX(3)*CX(2)))/(BOT*OLS(L,1))               DRE09940
      G(1,2)=((BX(3)*CX(1))-(BX(1)*CX(3)))/(BOT*OLS(L,2))               DRE09950
      G(1,3)=((BX(1)*CX(2))-(BX(2)*CX(1)))/(BOT*OLS(L,3))               DRE09960
      G(2,1)=((AX(3)*CX(2))-(AX(2)*CX(3)))/(BOT*OLS(L,1))               DRE09970
      G(2,2)=((AX(1)*CX(3))-(AX(3)*CX(1)))/(BOT*OLS(L,2))               DRE09980
      G(2,3)=((AX(2)*CX(1))-(AX(1)*CX(2)))/(BOT*OLS(L,3))               DRE09990
      G(3,1)=((AX(2)*BX(3))-(AX(3)*BX(2)))/(BOT*OLS(L,1))               DRE10000
      G(3,2)=((AX(3)*BX(1))-(AX(1)*BX(3)))/(BOT*OLS(L,2))               DRE10010
      G(3,3)=((AX(1)*BX(2))-(AX(2)*BX(1)))/(BOT*OLS(L,3))               DRE10020
C     DO 7933  J=1,3                                                    DRE10030
C7933    WRITE(3,7322)J,(G(J,I),I=1,3)                                  DRE10040
C7322  FORMAT(' ROW OF G ',I5,3(3X,F20.1))                              DRE10050
C****** FORM THE D MATRIX                                               DRE10060
      MT= NOD(L,4)                                                      DRE10070
      THIC= AMAT(MT,3)
      E=AMAT(MT,2)                                                      DRE10090
      VV=AMAT(MT,1)                                                     DRE10100
      COM= E/(1.0D0-VV*VV)                                              DRE10110
      DD(1,1)=COM                                                       DRE10120
      DD(1,2)=COM*VV                                                    DRE10130
      DD(1,3)=0.0D0                                                     DRE10140
      DD(2,1)=COM*VV                                                    DRE10150
      DD(2,2)=COM                                                       DRE10160
      DD(2,3)=0.0D0                                                     DRE10170
      DD(3,1)=0.0D0                                                     DRE10180
      DD(3,2)=0.0D0                                                     DRE10190
      DD(3,3)=COM*(1.0D0-VV)/2.0D0                                      DRE10200
C     DO 933  J=1,3                                                     DRE10210
C 933    WRITE(3,1322)J,(DD(J,I),I=1,3)                                 DRE10220
C1322  FORMAT(' ROW OF D ',I5,3(3X,F20.1))                              DRE10230
      DO 28 J=1,3                                                       DRE10240
      DO 28 K=1,3                                                       DRE10250
      VS(J,K ) = 0.0D0                                                  DRE10260
      DO 28 I=1,3                                                       DRE10270
  28  VS(J,K)=VS(J,K) +DD(J,I) *G(I,K)                                  DRE10280
      DO 29 J=1,3                                                       DRE10290
      ST(J) =0.0D0                                                      DRE10300
      DO 29 I=1,3                                                       DRE10310
  29  ST(J) = ST(J)+VS(J,I) *EX(I)                                      DRE10320
C      WRITE(3,4350)((VS(J,I),I=1,3),J=1,3) 
 4350 FORMAT('  ',3(2X ,F20.12))                                        DRE10340
      CCX=(ST(1)+ST(2))/2.0D0                                           DRE10350
      BBX=(ST(1)-ST(2))/2.0D0                                           DRE10360
      AAX=DSQRT((BBX**2)+(ST(3)**2))                                    DRE10370
      SMAX=CCX+AAX                                                      DRE10380
      SMIN =CCX-AAX                                                     DRE10390
      DIF=ST(2) -SMIN                                                   DRE10400
      IF(DABS(DIF).GT.0.0000001) GOTO 8923                              DRE10410
      BANG=0.0D0                                                        DRE10420
      GOTO 5688                                                         DRE10430
 8923 CONTINUE                                                          DRE10440
C***** HERE ELEMENT HAS ANGLE NE 0.0                                    DRE10450
      BANG = DATAN(ST(3)/DIF)                                           DRE10460
 5688 CONTINUE                                                          DRE10470
      WRITE(3,5411)ST(1),ST(2),ST(3),SMAX,SMIN,BANG                     DRE10480
 5411 FORMAT(' ST1=',F15.5,' ST2=',F15.5,'ST3=',F15.5,' SMAX=',F15.5,   DRE10490
     1' SMIN=',F15.5,' BANG=',F15.5)                                    DRE10500
      write(50,888)st(1),st(2),st(3),smax,smin
888   format(5e12.5)
 157  CONTINUE                                                          DRE10510
 3311 CONTINUE                                                          DRE10520
      RETURN                                                            DRE10530
      END                                                               DRE10540
      SUBROUTINE PRINT                                                  DRE10550
      include 'flat2.com' 
      DO 44 N=1,NN                                                      DRE10680
  44  WRITE(3,4) N,(CORD(N,I),I=1,NDIME)                                DRE10690
 4    FORMAT(BZ,I5,3F20.12)                                             DRE10700
      DO 2 N=1,NN                                                       DRE10710
 2    WRITE(3,233) N,(R(N,I),I=1,NDIME)                                 DRE10720
 233  FORMAT('  NODE = ',I5,5X,'RESIDUALS:',3(4X,E15.7))                DRE10730
      DO 6 N=1,NN                                                       DRE10740
  6   WRITE(3,633)N,(V(N,I),I=1,NDIME)                                  DRE10750
 633  FORMAT('  NODE = ',I5,5X,'VELOCITIES:',3(4X,E15.7))               DRE10760
      DO 8 N=1,NN                                                       DRE10770
  8   WRITE(3,833) N,(BB(N,I),I=1,NDIME)                                DRE10780
 833  FORMAT('  NODE = ',I5,5X,'MASSES:',3(4X,E15.7))                   DRE10790
      RETURN                                                            DRE10800
      END                                                               DRE10810
      SUBROUTINE NORMAL                                                 DRE13240
C**** CALCULATES NORMALS TO ELEMNTS                                     DRE13250
      include 'flat2.com'
C**** THIS SOULD ONLY BE CALLED AT THE PEAK                             DRE13380
      DO 20 J=1,NE                                                      DRE13390
      N1=NOD(J,1)                                                       DRE13400
      N2=NOD(J,2)                                                       DRE13410
      N3=NOD(J,3)                                                       DRE13420
      X21=CORD(N2,1)-CORD(N1,1)                                         DRE13430
      Y21=CORD(N2,2)-CORD(N1,2)                                         DRE13440
      Z21=CORD(N2,3)-CORD(N1,3)                                         DRE13450
      X31=CORD(N3,1)-CORD(N1,1)                                         DRE13460
      Y31=CORD(N3,2)-CORD(N1,2)                                         DRE13470
      Z31=CORD(N3,3)-CORD(N1,3)                                         DRE13480
      VX(J,1)=(Y21*Z31)-(Y31*Z21)                                       DRE13490
      VX(J,2)=(Z21*X31)-(Z31*X21)                                       DRE13500
      VX(J,3)=(X21*Y31)-(X31*Y21)                                       DRE13510
 20   CONTINUE      
C**** NORMALISE VECTOR
      DO 67 J=1,NE
      SUM=0.0D0
      DO 617 JJ=1,3
 617  SUM=SUM+VX(J,JJ)**2
      SUM=DSQRT(SUM)
      DO 627 JJ=1,3
 627  VX(J,JJ)=VX(J,JJ)/SUM
 67   CONTINUE     
      RETURN                                                            DRE13530
      END                                                               DRE13540
      SUBROUTINE PRESSU                                                 DRE13550
      include 'flat2.com'
C**** CALCULATES THE PRESSURE AT EACH NODE  DUE TO LOADING              DRE13680
      DO 571 J=1,NN                                                     DRE13690
      DO 571 JJ=1,3                                                     DRE13700
 571  PDIS(J,JJ) =0.0D0                                                 DRE13710
C**** LOADING AT EACH NODE FOR AN INTERNAL PRESSURE                     DRE13720
      DO 56 J=1,NE                                                      DRE13730
      MT=NOD(L,4)                                                       DRE13740
      PRES=AMAT(L,5)                                                    DRE13750
      DO 57 JD=1,3                                                      DRE13760
  57  P(JD)=V(J,JD)* PRES/6.0D0                                         DRE13770
C******ADD INTO APPROPRIATE LOCATIONS                                   DRE13780
      DO 671 JN=1,3                                                     DRE13790
      NODE = NOD(J,JN)                                                  DRE13800
      DO 671  JD=1,3                                                    DRE13810
  671 PDIS(NODE,JD)=PDIS(NODE,JD) +P(JD)                                DRE13820
  56  CONTINUE                                                          DRE13830
      RETURN                                                            DRE13840
      END                                                               DRE13850
      SUBROUTINE SLACKS                                                 DRE13860
      include 'flat2.com'
C***** CALCULATES SLACK LENGTHS OF GSTRINGS                             DRE13990
      DO 276 J=1,NSTR                                                   DRE14000
      IJ=NGS(J)                                                         DRE14010
      N1=NSNOD(J,1)                                                     DRE14020
      N2=NSNOD(J,IJ+1)                                                  DRE14030
      DL=0.0D0                                                          DRE14040
      DO 53 K=1,NDIME                                                   DRE14050
      D(K)=CORD(N2,K)-CORD(N1,K)                                        DRE14060
  53  DL=DL+D(K)**2                                                     DRE14070
      DL=DSQRT(DL)                                                      DRE14080
      RIJ=IJ                                                            DRE14090
      SLS(J)=DL/(RIJ*2.0D0)                                             DRE14100
      WRITE(3,45)J,SLS(J)
 45   FORMAT(//1X,'STRING=',I5,'LENGTH=',F20.5)  
 276  CONTINUE      
      RETURN                                                            DRE14110
      END                                                               DRE14120
      SUBROUTINE GSRES                                                  DRE14130
      include 'flat2.com'
C***** CALCULATE RESIDUALS FOR STRINGS                                  DRE14260
C**** REMEBER THIS MUST BE DONE FIRST                                   DRE14270
      DO 52 NS=1,NSTR                                                   DRE14280
      DO 52  L=1,NGS(NS)                                                DRE14290
      N1= NSNOD(NS,L)                                                   DRE14300
      N2= NSNOD(NS,L+1)                                                 DRE14310
      DL= 0.0D0                                                         DRE14320
      DO 53 K=1,NDIME                                                   DRE14330
      D(K) = CORD(N2,K)-CORD(N1,K)                                      DRE14340
   53 DL= DL +D(K)**2                                                   DRE14350
      DL = DSQRT(DL)                                                    DRE14360
      TTX=(DL-SLS(NS))*EES/SLS(NS)                                      DRE14390
      DO 102 K=1,NDIME                                                  DRE14400
      DR=D(K)*TTX/DL                                                    DRE14410
      R(N1,K)=R(N1,K) +DR                                               DRE14420
      R(N2,K)=R(N2,K) -DR                                               DRE14430
 102  CONTINUE                                                          DRE14440
  52  CONTINUE                                                          DRE14520
C***** ZERO NODES AT ENDS OF STRINGS                                    DRE14530
      DO 9952 NS=1,NSTR                                                 DRE14540
      N1= NSNOD(NS,1)                                                   DRE14550
      N2= NSNOD(NS,NGS(NS)+1)                                           DRE14560
      DO 9953 ND=1,NDIME                                                DRE14570
      R(N1,ND)=0.                                                       DRE14580
      R(N2,ND)=0.                                                       DRE14590
 9953 CONTINUE                                                          DRE14600
 9952 CONTINUE                                                          DRE14610
C******** SUBTRACT NORMAL COMPONRNTS
      DO 627  J=1,NN
      RJN=0.
      DO 625 I=1,NDIME
 625  RJN=RJN+SV(J,I)*R(J,I)
      DO 626  I=1,NDIME
      R(J,I)=R(J,I)-RJN*SV(J,I)
 626  CONTINUE
 627  CONTINUE
      RETURN
      END                                                               DRE14630
      SUBROUTINE AVENOR                                                 DRE14660
      include 'flat2.com'
C**** SHOULD BE CALLED AFTER PEAK IF STRINGS                            DRE14790
      CALL NORMAL                                                       DRE14800
      DO 27 J=1,NN
      DO 27 JX=1,NDIME
 27   SV(J,JX)=0.
      DO 57  J=1,NE
      DO 57  JX=1,3
      NODE=NOD(J,JX)
      DO 58 JU=1,NDIME
 58   SV(NODE,JU)=SV(NODE,JU)+VX(J,JU)
 57   CONTINUE
C**** NORMALISE THE VECTOR
      DO 53 J=1,NN
      SUM=0.0D0
      DO 89 K=1,NDIME
  89  SUM=SUM + (SV(J,K)**2)
      SUM=DSQRT(SUM)
      DO 899 K=1,NDIME
 899  SV(J,K)=SV(J,K)/SUM
 53   CONTINUE
      RETURN                                                            DRE15000
      END                                                               DRE15010
      SUBROUTINE GSDATA                                                 DRE15020
      include 'flat2.com'
      WRITE(*,*)'INPUT THE E VALUE FOR  G STRINGS'
      READ(*,*)EES    
      WRITE(3,345)EES
 345  FORMAT(' E VALUES FOR G STRINGS',E15.7)                        
C**** READS IN STRING DATA                                              DRE15150
      DO 58 J=1,NSTR                                                    DRE15160
      READ(2,781)NGS(J)                                                 DRE15170
  781 FORMAT(I5)                                                        DRE15180
      WRITE(3,782)J,NGS(J)                                              DRE15190
  782 FORMAT(' STRING NO=',I5,'   NO OF LINKS=',I5,/,                   DRE15200
     .' STRING NODE:')                                                  DRE15210
      READ(2,784)(NSNOD(J,JD),JD=1,(NGS(J)+1))                          DRE15220
  784 FORMAT (8I5)                                                       DRE15230
      WRITE(3,785) (NSNOD(J,JD),JD=1,(NGS(J)+1))                        DRE15240
  785 FORMAT(1X,8I5)
   58 CONTINUE                                                          DRE15250
C**** RESIDUAL TERMS FOR THE STRING                                     DRE15260
C**** CHECK ELEMENT NOS CONNECTED TO A STRING NODE                      DRE15270
      DO 126 J=1,NSTR                                                   DRE15280
      IJ=NGS(J)                                                         DRE15290
      DO 126 I=1,(IJ+1)                                                 DRE15300
      NODE=NSNOD(J,I)                                                   DRE15310
      IC=0                                                              DRE15320
      DO 128 N=1,NE                                                     DRE15330
      DO 123 IN=1,3                                                     DRE15340
      NX=NOD(N,IN)                                                      DRE15350
      IF(NX.NE.NODE) GOTO 123                                           DRE15360
      IC=IC+1                                                           DRE15370
      IECON(J,I,IC)=N                                                   DRE15380
      NECON(J,I) =IC                                                    DRE15390
  123 CONTINUE                                                          DRE15400
  128 CONTINUE                                                          DRE15410
  126 CONTINUE                                                          DRE15420
      DO 731 J=1,NSTR                                                   DRE15430
      WRITE(3,228)J,NGS(J)                                              DRE15440
  228 FORMAT(/,' STRING NO=',I5,5X,'NO OF LINKS IN THE STRING=',I5)     DRE15450
      NJS=NGS(J) +1                                                     DRE15460
      DO 732 I=1,NJS                                                    DRE15470
      WRITE(3,229) I,NSNOD(J,I)                                         DRE15480
  229 FORMAT(' NODE',I5,5X,' IN THE STRING IS=',I5)                     DRE15490
      WRITE(3,320)(IECON(J,I,IC),IC=1,NECON(J,I))                       DRE15500
  320 FORMAT(' CONNECTED TO ELEMENT NOS:',10(2X,I3))                    DRE15510
  732 CONTINUE                                                          DRE15520
  731 CONTINUE                                                          DRE15530
C****CALCULATE SLACK LENGTHS OF STRINGS
      CALL SLACKS
      RETURN                                                            DRE15540
      END                                                               DRE15550
      SUBROUTINE UDL                                                    DRE02430
C***** CALCULATES LOADS AND ADDS TO RESIDUALS                           DRE02440
      include 'flat2.com'
      DO 5 I=1,NN
      DO 5 J=1,NDIME
 5    RR(I,J)=0.0D0
      call NORMAL 
      CALL AREAT
C      WRITE(7,9)
C9     FORMAT(/' DIRECTION COSINES')
C      DO 10 I=1,NE
C10    WRITE(7,11)I,(VX(I,J),J=1,3)
C11    FORMAT(I5,3F15.5)
C****** DO FOR EACH ELEMENT                                             DRE05900
c*****  PRES AND PLANL ACT DOWNWARDS NEGATIVE
c****  SUCTION IS THEREFORE POSITIVE
      DO 57 L=1,NE                                                      DRE05910
      mt=NOD(L,4)
      planl=AMAT(mt,5)
      pres =AMAT(mt,5)
      wt(1)=VX(L,1)*pres*AREA(L)/3.0D0
      wt(2)=VX(L,2)*pres*AREA(L)/3.0D0
      wt(3)=(VX(L,3)*(pres+planl)-DEN)*AREA(L)/3.0D0
C      WRITE(7,14)(wt(JJ),JJ=1,3)
C14    FORMAT(3F10.5)
      do 40 i=1,3
      do 40 j=1,3
40    RR(NOD(L,i),j)=RR(NOD(L,i),j)+wt(j)
  57  CONTINUE
      RETURN                                                            DRE02730
      END                                                               DRE02740
      SUBROUTINE REACTION
      include 'flat2.com'
      DIMENSION REACT(100,3),TOTAL(3)
C..............................CALCULATE & OUTPUT REACTIONS
      DO 30 I=1,3
30    TOTAL(I)=0.0
      DO 60 I=1,NB
      DO 50 J=1,NDIME
      REACT(I,J)= -R(NCO(I,1),J)
50    TOTAL(J)=TOTAL(J)+REACT(I,J)
      WRITE(3,100)NCO(I,1),(REACT(I,J),J=1,3)
100   FORMAT(/' REACTION AT NODE ',I5,3X,'ARE ',3(3X,E15.7))
60    CONTINUE
      WRITE(3,110)
110   FORMAT(/' GLOBAL TOTAL REACTIONS')
      DO 120 I=1,3
120   WRITE(3,130)I,TOTAL(I)
130   FORMAT(/' DIRECTION ',3X,I5,3X,E15.7)
      RETURN
      END
      SUBROUTINE AREAT
      include 'flat2.com'
C****** DO FOR EACH ELEMENT                                             DRE05900
      DO 57 L=1,NE                                                      DRE05910
C**** HERE CURRENT SIDE LENGTHS CALCULATED                              DRE05930
      NX=1                                                              DRE05940
      DO 71 NSS=1,3                                                     DRE05950
      NS=4-NSS                                                          DRE05960
      N1=NOD(L,NS)                                                      DRE05970
      N2=NOD(L,NX)                                                      DRE05980
      NX=NS                                                             DRE05990
      DTL(L,NS)=0.0D0 
      DO 72 K=1,NDIME                                                   DRE06010
      DTT(L,NS,K) =CORD(N2,K)-CORD(N1,K) +10.D-25 
  72  DTL(L,NS) =DTL(L,NS) +(DTT(L,NS,K)**2) 
      DTL(L,NS) = DSQRT(DTL(L,NS))    
  71  CONTINUE                                                          DRE06050
C     WRITE(3,567)L,DTL(L,1),DTL(L,2),DTL(L,3) 
 567  FORMAT(' FORM LENGTHS IN AREAT',I5,3(5X,F10.2))                   DRE06070
      S=(DTL(L,1)+DTL(L,2)+DTL(L,3))/2.0D0                                
      AREA(L) =DSQRT(S*(S-DTL(L,1))*(S-DTL(L,2))*(S-DTL(L,3)))
C      WRITE(7,13)L,AREA(L)
C13    FORMAT(/' AREA ',I5,3X,F10.5)
 57    CONTINUE   
       RETURN
       END
