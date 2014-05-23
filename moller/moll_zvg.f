


C **********************************************************************
C *                                                                    *
C *              Espalhamento Moller com emissao de foton              *
C *             ------------------------------------------             *
C *                                                                    *
C *                 Roberval Walsh & Anibal Ramalho                    *
C *                ---------------------------------                   *
C *                                                                    *
C **********************************************************************


C                                                Ultima versao: 27/09/01






C Declaracao de variaveis

      IMPLICIT NONE
	 
      REAL*8 PROOTS,PLAMBDA,ALUM,ESIST
      INTEGER NX,NY
      INTEGER NEVENT
      INTEGER HVX,HVY
      REAL*8 PDHX,PDHY
C PREENCHER OS PARAMETROS ABAIXO COM OS VALORES ADEQUADOS
C #######################################################
C Energia c.m. (GeV)
      PARAMETER(PROOTS=0.5D3)
C Escala do fator de forma (GeV)
      PARAMETER(PLAMBDA=1.5D3)
C Luminosidade (pb^-1)
      PARAMETER(ALUM=50.D3)
C Erro sistematico
      PARAMETER(ESIST=2.D-2)
C Tamanho da malha
      PARAMETER(NX=8,NY=8)
C Numero de eventos
      PARAMETER(NEVENT=1000000)
C Escolha do par de acoplamentos para a malha
C ( 1 -> h1z, 2 -> h2z, 3 -> h3z, 4 -> h4z )
C ( 5 -> h1g, 6 -> h2g, 7 -> h3g, 8 -> h4g )
      PARAMETER(HVX=1,HVY=5)
C Tamanho dos acoplamentos (para a malha)
      PARAMETER(PDHX=4.1D-1,PDHY=1.5D-1)
C #######################################################
C Momentos iniciais e finais das particulas       
      REAL*8 PA(0:4),PB(0:4),P1(0:4),P2(0:4),P3(0:4)
      COMMON/MOMENTA/PA,PB,P1,P2,P3
C Vetor de polarizacao do foton
      COMPLEX*16 EPS(0:4)
      COMMON/FOTON/EPS
C Bilineares
      COMPLEX*16 ETA2A(0:4),ETA3B(0:4),ETA2B(0:4),ETA3A(0:4)
      COMMON/BILINEAR/ETA2A,ETA3B,ETA2B,ETA3A
C Helicidades      
      INTEGER HA,HB,H1,H2,H3
C Energia do centro de massa
      REAL*8 ROOTS,S
      COMMON/CME/ROOTS,S
C Limites de integracao
      REAL*8 EN1MIN,EN1MAX,CT1MIN,CT1MAX,CT2MIN,CT2MAX,
     #       PHI2MIN,PHI2MAX
      REAL*8 U1MIN,U1MAX,ETA1MIN,ETA1MAX,ETA2MIN,ETA2MAX
      COMMON/LIMITE1/EN1MIN,EN1MAX,CT1MIN,CT1MAX,CT2MIN,CT2MAX,
     #               PHI2MIN,PHI2MAX 
      COMMON/LIMITE2/U1MIN,U1MAX,ETA1MIN,ETA1MAX,ETA2MIN,ETA2MAX
      REAL*8 DPHI2,DU1,DETA1,DETA2
C Variaveis para cortes
      REAL*8 THETA,THETAEG,ENERGY
      COMMON/CORTES/THETA,THETAEG,ENERGY
C Peso e volume do espaco de fase e fluxo
      REAL*8 PSWT,PSVOL,FLUXO
C Parametros eletrofracos
      REAL*8 ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ,AMZ
      COMMON /ELECTROWEAK/ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      REAL*8 C1V,C1A,C2V,C2A
      COMMON/CVCA/C1V,C1A,C2V,C2A
      REAL*8 C1D,C1S,C2D,C2S
      COMMON/CVCA1/C1D,C1S,C2D,C2S
C Amplitude
      COMPLEX*16 AMP0
      COMPLEX*16 AMPZ(4),AMPG(4),AMP1(8)
      REAL*8 AMPSQ0
      REAL*8 AMPSQ1(NX,NY)
C Secao de choque total
      REAL*8 SIGMA,DSIGMA
      REAL*8 SIGMA1(NX,NY),DSIGMA1(NX,NY)
C Variaveis auxiliares
      INTEGER I,J
      REAL*8 PI
      PARAMETER (PI=3.14159265358979D0)
      INTEGER KSTATE
      REAL*8 SOMA01,SOMA02
      REAL*8 SOMA11(NX,NY),SOMA12(NX,NY)
      REAL*8 AMEDIA,FATOR
      REAL*8 RMOD2
      REAL*8 DOT3
      REAL*8 CT1,CT2,EN1,PT1
C Do gerador de numeros aleatorio RAN1 - Numerical Recipes
      INTEGER IDUM
      COMMON/RAN/IDUM
C Histogramas
      INTEGER NBIN,NDIM
      PARAMETER (NBIN=8,NDIM=NBIN+1)
C Escala de energia      
      REAL*8 ALAMBDA
      COMMON/ACOPL3/ALAMBDA
C  Storage Matrices
      REAL*8 DX,BX
      REAL*8 VX01(NDIM),VX02(NDIM),VX03(NDIM),VX04(NDIM)
      REAL*8 VY01(NBIN),VY02(NBIN),VY03(NBIN),VY04(NBIN)
      REAL*8 AREA01,AREA02,AREA03,AREA04
      INTEGER IX,IY
      INTEGER NHX,NHY
      REAL*8 DHX,DHY
      REAL*8 AH0X(NX),AH0Y(NY)
      REAL*8 VY1(NBIN,NX,NY),VY2(NBIN,NX,NY)
      REAL*8 VY3(NBIN,NX,NY),VY4(NBIN,NX,NY)
      REAL*8 VYAUX(NBIN)
      REAL*8 AREA1(NX,NY),AREA2(NX,NY),AREA3(NX,NY),AREA4(NX,NY)
      REAL*8 ERRO21(NBIN),ERRO22(NBIN),ERRO23(NBIN),ERRO24(NBIN)
      REAL*8 CHI21(NX,NY),CHI22(NX,NY),CHI23(NX,NY),CHI24(NX,NY)

C ----------------------------------------------------------------------
C Inicio do programa pricipal

      OPEN(14,FILE='X2_PT1.DAT',STATUS='UNKNOWN')
      OPEN(15,FILE='SM_PT1.DAT',STATUS='UNKNOWN')
      OPEN(16,FILE='XM_PT1.DAT',STATUS='UNKNOWN')
      OPEN(20,FILE='SIGTOT.DAT',STATUS='UNKNOWN')

C Dimensao da malha
      NHX=NX
      NHY=NY
C Tamanho do acoplamento
      DHX=PDHX
      DHY=PDHY

      CALL ACOPL(NHX,NHY,DHX,DHY,AH0X,AH0Y)

C Escala de energia
      ALAMBDA=PLAMBDA

C Inicializacao do gerador de numeros aleatorios       
      IDUM=-100
	 
C Energia de centro de massa       
      ROOTS=PROOTS
      S=ROOTS**2
      FLUXO=2.D0*S
	
C Acoplamentos e parametros eletrofracos
      ALPHA=1.D0/128.D0
      STW2=0.23124D0
      CTW2=1.D0-STW2
      AMZ=91.187D0
      AMZ2=AMZ*AMZ
      GZMZ=2.49D0*AMZ
      C1V=-1.D0
      C1A=0.D0
      C2V=(STW2-0.25D0)/DSQRT(STW2*CTW2)
      C2A=-1.D0/(4.D0*DSQRT(STW2*CTW2))
      C1D=C1V-C1A
      C1S=C1V+C1A
      C2D=C2V-C2A
      C2S=C2V+C2A

C Massas das particulas       
      PA(4)=0.D0
      PB(4)=0.D0
      P1(4)=0.D0
      P2(4)=0.D0
      P3(4)=0.D0
	 
C Momentos das particulas iniciais
      PA(0)=ROOTS/2.D0
      PA(1)=0.D0
      PA(2)=0.D0
      PA(3)=ROOTS/2.D0
      PB(0)=ROOTS/2.D0
      PB(1)=0.D0
      PB(2)=0.D0
      PB(3)=-ROOTS/2.D0

C Cortes
      ENERGY=10.D0
      THETA=5.D0*PI/180.D0
      THETAEG=1.D0*PI/180.D0
	 
C Limites de integracao
      EN1MIN=ENERGY
      EN1MAX=ROOTS/2.D0
      CT1MIN=-DCOS(THETA)
      CT1MAX=DCOS(THETA)
      CT2MIN=-DCOS(THETA)
      CT2MAX=DCOS(THETA)
      PHI2MIN=0.D0
      PHI2MAX=2.D0*PI

      U1MAX=-DLOG(2.D0*EN1MIN/ROOTS)
      U1MIN=-DLOG(2.D0*EN1MAX/ROOTS)
      ETA1MAX=0.5D0*DLOG((1.D0+CT1MAX)/(1.D0-CT1MAX))
      ETA1MIN=0.5D0*DLOG((1.D0+CT1MIN)/(1.D0-CT1MIN))
      ETA2MAX=0.5D0*DLOG((1.D0+CT2MAX)/(1.D0-CT2MAX))
      ETA2MIN=0.5D0*DLOG((1.D0+CT2MIN)/(1.D0-CT2MIN))
		  
C Volume do espaco de fase
      DU1=U1MAX-U1MIN
      DETA1=ETA1MAX-ETA1MIN
      DETA2=ETA2MAX-ETA2MIN
      DPHI2=PHI2MAX-PHI2MIN

      PSVOL=DU1*DETA1*DETA2*DPHI2

C Preparando os bins
      DO I=1,NBIN
       VY01(I)=0.D0
       VY02(I)=0.D0
       VY03(I)=0.D0
       VY04(I)=0.D0
       DO IX=1,NHX
        DO IY=1,NHY
      	  VY1(I,IX,IY)=0.D0
      	  VY2(I,IX,IY)=0.D0
      	  VY3(I,IX,IY)=0.D0
      	  VY4(I,IX,IY)=0.D0
        END DO
       END DO
      END DO
      DX=1.D0
      BX=DX/NBIN
      DO I=1,NDIM
       VX01(I)=-1.D0+(I-1)*(BX*2.D0)  ! CT1
       VX02(I)=VX01(I)                ! CT2
       VX03(I)=(I-1)*BX               ! EN1
       VX04(I)=VX03(I)                ! PT1
      END DO

C Inicializacao de variaveis
      AMPSQ0=0.D0
      SOMA01=0.D0
      SOMA02=0.D0
      AREA01=0.D0
      AREA02=0.D0
      AREA03=0.D0
      AREA04=0.D0
      AMP0=(0.D0,0.D0)
      DO IX=1,NHX
       DO IY=1,NHY
        AMPSQ1(IX,IY)=0.D0
        SOMA11(IX,IY)=0.D0
        SOMA12(IX,IY)=0.D0
        AREA1(IX,IY)=0.D0
        AREA2(IX,IY)=0.D0
        AREA3(IX,IY)=0.D0
        AREA4(IX,IY)=0.D0
       END DO
      END DO
      DO I=1,4
       AMPG(I)=(0.D0,0.D0)
       AMPZ(I)=(0.D0,0.D0)
      END DO
      DO I=1,8
       AMP1(I)=(0.D0,0.D0)
      END DO

      DO I=1,NEVENT

       CALL CINEMATICA(KSTATE,PSWT)
		  
       IF (KSTATE.EQ.1) THEN

        CT1=DOT3(PA,P1)/DSQRT(DOT3(PA,PA)*DOT3(P1,P1))
        CT2=DOT3(PA,P2)/DSQRT(DOT3(PA,PA)*DOT3(P2,P2))
        EN1=2.D0*P1(0)/ROOTS
        PT1=DSQRT(P1(1)**2+P1(2)**2)*2.D0/ROOTS
				  
        DO HA=-1,1,2
         DO HB=-1,1,2
          DO H1=-1,1,2
          CALL GEPS(P1,H1,EPS)
           DO H2=-1,1,2
            CALL UBGU(P2,H2,PA,HA,ETA2A)
            CALL UBGU(P2,H2,PB,HB,ETA2B)
            DO H3=-1,1,2
             CALL UBGU(P3,H3,PB,HB,ETA3B)
             CALL UBGU(P3,H3,PA,HA,ETA3A)
             CALL AMPLTD(HA,HB,H1,H2,H3,AMP0)
             AMPSQ0=AMPSQ0+RMOD2(AMP0)
             IF ((HVX.LE.4).OR.(HVY.LE.4)) THEN
              CALL AMPLTD1(HA,HB,H1,H2,H3,AMPZ)
             END IF
             IF ((HVX.GE.5).OR.(HVY.GE.5)) THEN
              CALL AMPLTD2(HA,HB,H1,H2,H3,AMPG)
             END IF
             DO J=1,4
              AMP1(J)=AMPZ(J)
              AMP1(J+4)=AMPG(J)
             END DO
             DO IX=1,NHX
              DO IY=1,NHY
               AMPSQ1(IX,IY)=AMPSQ1(IX,IY)+
     #         RMOD2(AMP0+(AH0X(IX))*AMP1(HVX)+
     #         (AH0Y(IY))*AMP1(HVY))
              END DO
             END DO
            END DO
           END DO
          END DO
         END DO
        END DO

        AMPSQ0=AMPSQ0*PSWT
        SOMA01=SOMA01+AMPSQ0
        SOMA02=SOMA02+AMPSQ0**2
C  Secoes de choque diferenciais - Modelo Padrao
        CALL STORE(VX01,NBIN,VY01,CT1,AMPSQ0)
        CALL STORE(VX02,NBIN,VY02,CT2,AMPSQ0)
        CALL STORE(VX03,NBIN,VY03,EN1,AMPSQ0)
        CALL STORE(VX04,NBIN,VY04,PT1,AMPSQ0)
        AMPSQ0=0.D0

        DO IX=1,NHX
         DO IY=1,NHY
          AMPSQ1(IX,IY)=AMPSQ1(IX,IY)*PSWT
          SOMA11(IX,IY)=SOMA11(IX,IY)+AMPSQ1(IX,IY)
          SOMA12(IX,IY)=SOMA12(IX,IY)+AMPSQ1(IX,IY)**2
         END DO
        END DO
C  Secoes de choque diferenciais - Modelo Estendido
        CALL STORE1(VX01,NBIN,VY1,CT1,AMPSQ1,NHX,NHY)
        CALL STORE1(VX02,NBIN,VY2,CT2,AMPSQ1,NHX,NHY)
        CALL STORE1(VX03,NBIN,VY3,EN1,AMPSQ1,NHX,NHY)
        CALL STORE1(VX04,NBIN,VY4,PT1,AMPSQ1,NHX,NHY)

        DO IX=1,NHX
         DO IY=1,NHY
          AMPSQ1(IX,IY)=0.D0
         END DO
        END DO

       END IF
      END DO
	
C Secao de choque total - Modelo Padrao -----------------------------
      AMEDIA=0.25D0
      SOMA01=SOMA01*PSVOL/NEVENT
      SOMA02=SOMA02*PSVOL**2/NEVENT
      DSIGMA=DSQRT((SOMA02-SOMA01**2)/NEVENT)
      FATOR=3.8939D08*(ALPHA*4.D0*PI)**3/FLUXO/(2.D0*PI)**5
C Fator de simetria
      FATOR=FATOR*0.5D0
      SIGMA=SOMA01*FATOR*AMEDIA
      DSIGMA=DSIGMA*FATOR*AMEDIA

C  Normalization factors for distributions
      CALL AREAD(NBIN,BX*2.D0,VY01,AREA01)  ! CT1
      CALL AREAD(NBIN,BX*2.D0,VY02,AREA02)  ! CT2
      CALL AREAD(NBIN,BX,VY03,AREA03)       ! EN1
      CALL AREAD(NBIN,BX,VY04,AREA04)       ! PT1
      AREA01=SIGMA/AREA01
      AREA02=SIGMA/AREA02
      AREA03=SIGMA/AREA03
      AREA04=SIGMA/AREA04
C  Prepare bins
      DO I=1,NBIN
       VY01(I)=AREA01*VY01(I)
       VY02(I)=AREA02*VY02(I)
       VY03(I)=AREA03*VY03(I)
       VY04(I)=AREA04*VY04(I)
      END DO
C -------------------------------------------------------------------
C Secao de choque total - Modelo estendido  -------------------------
      DO IX=1,NHX
       DO IY=1,NHY
        SOMA11(IX,IY)=SOMA11(IX,IY)*PSVOL/NEVENT
        SOMA12(IX,IY)=SOMA12(IX,IY)*PSVOL**2/NEVENT
        DSIGMA1(IX,IY)=DSQRT((SOMA12(IX,IY)-SOMA11(IX,IY)**2)/NEVENT)
        SIGMA1(IX,IY)=SOMA11(IX,IY)*FATOR*AMEDIA
        DSIGMA1(IX,IY)=DSIGMA1(IX,IY)*FATOR*AMEDIA
       END DO
      END DO

C  Normalization factors for distributions
      CALL AREAD1(NBIN,BX*2.D0,VY1,AREA1,NHX,NHY)  ! CT1
      CALL AREAD1(NBIN,BX*2.D0,VY2,AREA2,NHX,NHY)  ! CT2
      CALL AREAD1(NBIN,BX,VY3,AREA3,NHX,NHY)       ! EN1
      CALL AREAD1(NBIN,BX,VY4,AREA4,NHX,NHY)       ! PT1

      DO IX=1,NHX
       DO IY=1,NHY
        AREA1(IX,IY)=SIGMA1(IX,IY)/AREA1(IX,IY)
        AREA2(IX,IY)=SIGMA1(IX,IY)/AREA2(IX,IY)
        AREA3(IX,IY)=SIGMA1(IX,IY)/AREA3(IX,IY)
        AREA4(IX,IY)=SIGMA1(IX,IY)/AREA4(IX,IY)
C  Prepare bins
        DO I=1,NBIN
         VY1(I,IX,IY)=AREA1(IX,IY)*VY1(I,IX,IY)
         VY2(I,IX,IY)=AREA2(IX,IY)*VY2(I,IX,IY)
         VY3(I,IX,IY)=AREA3(IX,IY)*VY3(I,IX,IY)
         VY4(I,IX,IY)=AREA4(IX,IY)*VY4(I,IX,IY)
        END DO
       END DO
      END DO

      CALL FPLOT(15,NBIN,VX04,VY04)
      DO IX=1,NHX
       DO IY=1,NHY
        DO I=1,NBIN
         VYAUX(I)=VY4(I,IX,IY)
        END DO
        CALL FPLOT(16,NBIN,VX04,VYAUX)
       END DO
      END DO
C -------------------------------------------------------------------
C Calculo do chi^2
      DO I=1,NBIN
       VY01(I)=VY01(I)*ALUM*BX*2.D0
       VY02(I)=VY02(I)*ALUM*BX*2.D0
       VY03(I)=VY03(I)*ALUM*BX
       VY04(I)=VY04(I)*ALUM*BX
       ERRO21(I)=VY01(I)+(ESIST*VY01(I))**2
       ERRO22(I)=VY02(I)+(ESIST*VY02(I))**2
       ERRO23(I)=VY03(I)+(ESIST*VY03(I))**2
       ERRO24(I)=VY04(I)+(ESIST*VY04(I))**2
      END DO

      WRITE(20,*) SIGMA,DSIGMA
	   
      DO IX=1,NHX
       DO IY=1,NHY
        DO I=1,NBIN
         VY1(I,IX,IY)=VY1(I,IX,IY)*ALUM*BX*2.D0
         VY2(I,IX,IY)=VY2(I,IX,IY)*ALUM*BX*2.D0
         VY3(I,IX,IY)=VY3(I,IX,IY)*ALUM*BX
         VY4(I,IX,IY)=VY4(I,IX,IY)*ALUM*BX
         CHI21(IX,IY)=CHI21(IX,IY)+(VY1(I,IX,IY)-VY01(I))**2/ERRO21(I)
         CHI22(IX,IY)=CHI22(IX,IY)+(VY2(I,IX,IY)-VY02(I))**2/ERRO22(I)
         CHI23(IX,IY)=CHI23(IX,IY)+(VY3(I,IX,IY)-VY03(I))**2/ERRO23(I)
         CHI24(IX,IY)=CHI24(IX,IY)+(VY4(I,IX,IY)-VY04(I))**2/ERRO24(I)
        END DO
c        WRITE(11,*) AH0X(IX),AH0Y(IY),CHI21(IX,IY)
c        WRITE(12,*) AH0X(IX),AH0Y(IY),CHI22(IX,IY)
c        WRITE(13,*) AH0X(IX),AH0Y(IY),CHI23(IX,IY)
        WRITE(14,*) AH0X(IX),AH0Y(IY),CHI24(IX,IY)
        WRITE(20,*) SIGMA1(IX,IY),DSIGMA1(IX,IY)
       END DO
      END DO
	

      PRINT*,SIGMA,DSIGMA,NEVENT
	
      STOP
      END
C **********************************************************************

C **********************************************************************

C                             Cinematica

C **********************************************************************

      SUBROUTINE CINEMATICA(KSTATE,PSWT)

C                        Declaracao de variaveis
C                       -------------------------
   
      IMPLICIT NONE
	 
C Momentos iniciais, finais e internos
      REAL*8 PA(0:4),PB(0:4),P1(0:4),P2(0:4),P3(0:4)
      COMMON/MOMENTA/PA,PB,P1,P2,P3
C Qij=Pi-Pj
C Kij=Pi+Pj
      REAL*8 QA1(0:4),QB3(0:4),QB1(0:4),Q2A(0:4),K12(0:4),K13(0:4)
      REAL*8 QB2(0:4),Q3A(0:4),QA2(0:4),Q3B(0:4),QA3(0:4),Q2B(0:4)
      COMMON/AUXMOM/QA1,QB3,QB1,Q2A,K12,K13,QB2,Q3A,QA2,Q3B,QA3,Q2B

C Energia do centro de massa
      REAL*8 ROOTS,S
      COMMON/CME/ROOTS,S
C Variaveis sorteadas       
      REAL*8 X1,X2,X3,X4
C Variaveis de integracao       
      REAL*8 U1,ETA1,ETA2
      REAL*8 EN1,CT1,CT2,PHI2
C Limites de integracao
      REAL*8 EN1MIN,EN1MAX,CT1MIN,CT1MAX,CT2MIN,CT2MAX,
     #       PHI2MIN,PHI2MAX
      REAL*8 U1MIN,U1MAX,ETA1MIN,ETA1MAX,ETA2MIN,ETA2MAX
      COMMON/LIMITE1/EN1MIN,EN1MAX,CT1MIN,CT1MAX,CT2MIN,CT2MAX,
     #               PHI2MIN,PHI2MAX 
      COMMON/LIMITE2/U1MIN,U1MAX,ETA1MIN,ETA1MAX,ETA2MIN,ETA2MAX
C Cortes
      REAL*8 THETA,THETAEG,ENERGY
      COMMON/CORTES/THETA,THETAEG,ENERGY

      REAL*8 CT3,T12,T13
C Variaveis auxiliares
      REAL*8 ST1,ST2,CT12
      REAL*8 PI
      PARAMETER (PI=3.14159265358979D0)
      INTEGER KSTATE
C Peso do espaco de fase       
      REAL*8 PSWT
C Do gerador de numeros aleatorio RAN1 - Numerical Recipes
      REAL*8 RAN1
      INTEGER IDUM
      COMMON/RAN/IDUM
      INTEGER J
      REAL*8 DOT3
	

C                               Inicio
C                              --------

      KSTATE=1
      PSWT=0.D0

C Massas das particulas       
      PA(4)=0.D0
      PB(4)=0.D0
      P1(4)=0.D0
      P2(4)=0.D0
      P3(4)=0.D0
	 
C Sorteio
      X1=RAN1(IDUM)
      X2=RAN1(IDUM)
      X3=RAN1(IDUM)
      X4=RAN1(IDUM)
      U1=(U1MAX-U1MIN)*X1+U1MIN
      ETA1=(ETA1MAX-ETA1MIN)*X2+ETA1MIN
      ETA2=(ETA2MAX-ETA2MIN)*X3+ETA2MIN
      PHI2=(PHI2MAX-PHI2MIN)*X4+PHI2MIN

C  Mudanca de variaveis
      EN1=0.5D0*ROOTS*DEXP(-U1)
      CT1=DTANH(ETA1)
      CT2=DTANH(ETA2)

C Identidades
      ST1=DSQRT(1.D0-CT1**2)
      ST2=DSQRT(1.D0-CT2**2)
      CT12=ST1*ST2*DCOS(PHI2)+CT1*CT2

Calculo das componentes dos momentos finais
      P1(0)=EN1
      P1(1)=P1(0)*ST1
      P1(2)=0.D0
      P1(3)=P1(0)*CT1
	
      P2(0)=(2.D0*ROOTS*EN1-S)/(2.D0*(EN1*(1.D0-CT12)-ROOTS))
      P2(1)=P2(0)*ST2*DCOS(PHI2)
      P2(2)=P2(0)*ST2*DSIN(PHI2)
      P2(3)=P2(0)*CT2
	
      P3(0)=PA(0)+PB(0)-P1(0)-P2(0)
      P3(1)=PA(1)+PB(1)-P1(1)-P2(1)
      P3(2)=PA(2)+PB(2)-P1(2)-P2(2)
      P3(3)=PA(3)+PB(3)-P1(3)-P2(3)

C Calculo de variaveis para os cortes
      CT3=DOT3(P3,PA)/DSQRT(DOT3(P3,P3)*DOT3(PA,PA))
      T12=DOT3(P1,P2)/DSQRT(DOT3(P1,P1)*DOT3(P2,P2))
      T13=DOT3(P1,P3)/DSQRT(DOT3(P1,P1)*DOT3(P3,P3))


C Cortes
      IF (DABS(CT1).GT.DCOS(THETA)) KSTATE=-1
      IF (DABS(CT2).GT.DCOS(THETA)) KSTATE=-1
      IF (DABS(CT3).GT.DCOS(THETA)) KSTATE=-1
      IF (T12.GT.DCOS(THETAEG)) KSTATE=-1
      IF (T13.GT.DCOS(THETAEG)) KSTATE=-1
      IF (P1(0).LT.ENERGY) KSTATE=-1
      IF (P2(0).LT.ENERGY) KSTATE=-1
      IF (P3(0).LT.ENERGY) KSTATE=-1

C Peso do espaco de fase
      PSWT=(2.D0*PI/8.D0)*(P1(0)*P2(0)/DABS(P1(0)*(1.D0-CT12)-ROOTS))
     #     *ST1**2
     #     *EN1
     #     *ST2**2

C Calculo do momentos internos
      DO J=0,3
       QA1(J)=PA(J)-P1(J)
       QB3(J)=PB(J)-P3(J)
       QB1(J)=PB(J)-P1(J)
       Q2A(J)=P2(J)-PA(J)
       K12(J)=P1(J)+P2(J)
       K13(J)=P1(J)+P3(J)
       QB2(J)=PB(J)-P2(J)
       Q3A(J)=P3(J)-PA(J)
       QA2(J)=-Q2A(J)           
       Q3B(J)=-QB3(J)           
       QA3(J)=-Q3A(J)           
       Q2B(J)=-QB2(J)
       Q3B(J)=-QB3(J)
      END DO

C                                 Fim
C                                _____

      RETURN
      END

C **********************************************************************


C **********************************************************************

C              Calculo das amplitudes no modelo padrao

C **********************************************************************

      SUBROUTINE AMPLTD(HA,HB,H1,H2,H3,AMP0)

C                     Declaracao de variaveis
C                    -------------------------
      IMPLICIT NONE
	
C Amplitudes
      COMPLEX*16 AMP(8),AMP0

C Helicidades
      INTEGER HA,HB,H1,H2,H3

C Momentos iniciais, finais e internos
      REAL*8 PA(0:4),PB(0:4),P1(0:4),P2(0:4),P3(0:4)
      COMMON/MOMENTA/PA,PB,P1,P2,P3
      REAL*8 QA1(0:4),QB3(0:4),QB1(0:4),Q2A(0:4),K12(0:4),K13(0:4)
      REAL*8 QB2(0:4),Q3A(0:4),QA2(0:4),Q3B(0:4),QA3(0:4),Q2B(0:4)
      COMMON /AUXMOM/QA1,QB3,QB1,Q2A,K12,K13,QB2,Q3A,QA2,Q3B,QA3,Q2B

C Bilinears; vetor de polarizacao
      COMPLEX*16 ETA2A(0:4),ETA3B(0:4),ETA2B(0:4),ETA3A(0:4)
      COMMON /BILINEAR/ETA2A,ETA3B,ETA2B,ETA3A
      COMPLEX*16 EPS(0:4)
      COMMON /FOTON/EPS

C Parametros eletrofracos
      REAL*8 ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      COMMON /ELECTROWEAK/ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      REAL*8 C1D,C1S,C2D,C2S
      COMMON /CVCA1/C1D,C1S,C2D,C2S

C Diversos
      REAL*8 DK,DOT4
      COMPLEX*16 ZDOTR,ZDOTZ,ZZZR,DZ


C                               Inicio
C                              --------

c M1+M2 (NOVA)
      AMP(1)=1.D0/DOT4(QA1,QA1)*
     #       (ZDOTR(ETA3B,QA1)*ZDOTZ(ETA2A,EPS)-
     #        ZDOTZ(ETA3B,EPS)*ZDOTR(ETA2A,QA1)+
     #        ZDOTZ(ETA3B,ETA2A)*ZDOTR(EPS,QA1)+
     #        (0.D0,1.D0)*HA*ZZZR(ETA3B,ETA2A,EPS,QA1))*
     #       ((C1D*C1D*DK(HA,+1)+C1S*C1S*DK(HA,-1))*
     #        (C1D*DK(HB,+1)+C1S*DK(HB,-1))/DOT4(QB3,QB3)+
     #        (C1D*C2D*DK(HA,+1)+C1S*C2S*DK(HA,-1))*
     #        (C2D*DK(HB,+1)+C2S*DK(HB,-1))/DZ(DOT4(QB3,QB3)))


c M3+M4 (NOVA)
      AMP(2)=1.D0/DOT4(QB1,QB1)*
     #       (ZDOTR(ETA2A,QB1)*ZDOTZ(ETA3B,EPS)-
     #        ZDOTZ(ETA2A,EPS)*ZDOTR(ETA3B,QB1)+
     #        ZDOTZ(ETA2A,ETA3B)*ZDOTR(EPS,QB1)+
     #        (0.D0,1.D0)*HB*ZZZR(ETA2A,ETA3B,EPS,QB1))*
     #       ((C1D*C1D*DK(HB,+1)+C1S*C1S*DK(HB,-1))*
     #        (C1D*DK(HA,+1)+C1S*DK(HA,-1))/DOT4(Q2A,Q2A)+
     #        (C1D*C2D*DK(HB,+1)+C1S*C2S*DK(HB,-1))*
     #        (C2D*DK(HA,+1)+C2S*DK(HA,-1))/DZ(DOT4(Q2A,Q2A)))


c M5+M6 (NOVA)
      AMP(3)=1.D0/DOT4(K12,K12)*
     #       (ZDOTR(ETA3B,K12)*ZDOTZ(ETA2A,EPS)-
     #        ZDOTZ(ETA3B,EPS)*ZDOTR(ETA2A,K12)+
     #        ZDOTZ(ETA3B,ETA2A)*ZDOTR(EPS,K12)-
     #        (0.D0,1.D0)*HA*ZZZR(ETA3B,ETA2A,EPS,K12))*
     #       ((C1D*C1D*DK(HA,+1)+C1S*C1S*DK(HA,-1))*
     #        (C1D*DK(HB,+1)+C1S*DK(HB,-1))/DOT4(QB3,QB3)+
     #        (C1D*C2D*DK(HA,+1)+C1S*C2S*DK(HA,-1))*
     #        (C2D*DK(HB,+1)+C2S*DK(HB,-1))/DZ(DOT4(QB3,QB3)))


c M7+M8 (NOVA)
      AMP(4)=1.D0/DOT4(K13,K13)*
     #       (ZDOTR(ETA2A,K13)*ZDOTZ(ETA3B,EPS)-
     #        ZDOTZ(ETA2A,EPS)*ZDOTR(ETA3B,K13)+
     #        ZDOTZ(ETA2A,ETA3B)*ZDOTR(EPS,K13)-
     #        (0.D0,1.D0)*HB*ZZZR(ETA2A,ETA3B,EPS,K13))*
     #       ((C1D*C1D*DK(HB,+1)+C1S*C1S*DK(HB,-1))*
     #        (C1D*DK(HA,+1)+C1S*DK(HA,-1))/DOT4(Q2A,Q2A)+
     #        (C1D*C2D*DK(HB,+1)+C1S*C2S*DK(HB,-1))*
     #        (C2D*DK(HA,+1)+C2S*DK(HA,-1))/DZ(DOT4(Q2A,Q2A)))

	
c M9+M10 (NOVA)
      AMP(5)=1.D0/DOT4(QA1,QA1)*
     #       (ZDOTR(ETA2B,QA1)*ZDOTZ(ETA3A,EPS)-
     #        ZDOTZ(ETA2B,EPS)*ZDOTR(ETA3A,QA1)+
     #        ZDOTZ(ETA2B,ETA3A)*ZDOTR(EPS,QA1)+
     #        (0.D0,1.D0)*HA*ZZZR(ETA2B,ETA3A,EPS,QA1))*
     #       ((C1D*C1D*DK(HA,+1)+C1S*C1S*DK(HA,-1))*
     #        (C1D*DK(HB,+1)+C1S*DK(HB,-1))/DOT4(QB2,QB2)+
     #        (C1D*C2D*DK(HA,+1)+C1S*C2S*DK(HA,-1))*
     #        (C2D*DK(HB,+1)+C2S*DK(HB,-1))/DZ(DOT4(QB2,QB2)))


c M11+M12 (NOVA)
      AMP(6)=1.D0/DOT4(QB1,QB1)*
     #       (ZDOTR(ETA3A,QB1)*ZDOTZ(ETA2B,EPS)-
     #        ZDOTZ(ETA3A,EPS)*ZDOTR(ETA2B,QB1)+
     #        ZDOTZ(ETA3A,ETA2B)*ZDOTR(EPS,QB1)+
     #        (0.D0,1.D0)*HB*ZZZR(ETA3A,ETA2B,EPS,QB1))*
     #       ((C1D*C1D*DK(HB,+1)+C1S*C1S*DK(HB,-1))*
     #        (C1D*DK(HA,+1)+C1S*DK(HA,-1))/DOT4(Q3A,Q3A)+
     #        (C1D*C2D*DK(HB,+1)+C1S*C2S*DK(HB,-1))*
     #        (C2D*DK(HA,+1)+C2S*DK(HA,-1))/DZ(DOT4(Q3A,Q3A)))


c M13+M14 (NOVA)
      AMP(7)=1.D0/DOT4(K13,K13)*
     #       (ZDOTR(ETA2B,K13)*ZDOTZ(ETA3A,EPS)-
     #        ZDOTZ(ETA2B,EPS)*ZDOTR(ETA3A,K13)+
     #        ZDOTZ(ETA2B,ETA3A)*ZDOTR(EPS,K13)-
     #        (0.D0,1.D0)*HA*ZZZR(ETA2B,ETA3A,EPS,K13))*
     #       ((C1D*C1D*DK(HA,+1)+C1S*C1S*DK(HA,-1))*
     #        (C1D*DK(HB,+1)+C1S*DK(HB,-1))/DOT4(QB2,QB2)+
     #        (C1D*C2D*DK(HA,+1)+C1S*C2S*DK(HA,-1))*
     #        (C2D*DK(HB,+1)+C2S*DK(HB,-1))/DZ(DOT4(QB2,QB2)))


c M15+M16 (NOVA)
      AMP(8)=1.D0/DOT4(K12,K12)*
     #       (ZDOTR(ETA3A,K12)*ZDOTZ(ETA2B,EPS)-
     #        ZDOTZ(ETA3A,EPS)*ZDOTR(ETA2B,K12)+
     #        ZDOTZ(ETA3A,ETA2B)*ZDOTR(EPS,K12)-
     #        (0.D0,1.D0)*HB*ZZZR(ETA3A,ETA2B,EPS,K12))*
     #       ((C1D*C1D*DK(HB,+1)+C1S*C1S*DK(HB,-1))*
     #        (C1D*DK(HA,+1)+C1S*DK(HA,-1))/DOT4(Q3A,Q3A)+
     #        (C1D*C2D*DK(HB,+1)+C1S*C2S*DK(HB,-1))*
     #        (C2D*DK(HA,+1)+C2S*DK(HA,-1))/DZ(DOT4(Q3A,Q3A)))


c M1+M2+M3+M4+M5+M6+M7+M8-M9-M10-M11-M12-M13-M14-M15-M16
      AMP0=AMP(1)+AMP(2)+AMP(3)+AMP(4)-AMP(5)-AMP(6)-AMP(7)-AMP(8)


C                                 Fim
C                                _____

      RETURN
      END

C **********************************************************************

C **********************************************************************

C             Calculo das amplitudes com o vertice ZgZ

C **********************************************************************

      SUBROUTINE AMPLTD1(HA,HB,H1,H2,H3,AMP)

C                     Declaracao de variaveis
C                    -------------------------
      IMPLICIT NONE
	
C Amplitudes
      COMPLEX*16 AMP(4),AMP17(4),AMP18(4)

C Helicidades
      INTEGER HA,HB,H1,H2,H3

C Momentos iniciais, finais e internos
      REAL*8 PA(0:4),PB(0:4),P1(0:4),P2(0:4),P3(0:4)
      COMMON/MOMENTA/PA,PB,P1,P2,P3
      REAL*8 QA1(0:4),QB3(0:4),QB1(0:4),Q2A(0:4),K12(0:4),K13(0:4)
      REAL*8 QB2(0:4),Q3A(0:4),QA2(0:4),Q3B(0:4),QA3(0:4),Q2B(0:4)
      COMMON /AUXMOM/QA1,QB3,QB1,Q2A,K12,K13,QB2,Q3A,QA2,Q3B,QA3,Q2B

C Bilinears; vetor de polarizacao
      COMPLEX*16 ETA2A(0:4),ETA3B(0:4),ETA2B(0:4),ETA3A(0:4)
      COMMON /BILINEAR/ETA2A,ETA3B,ETA2B,ETA3A
      COMPLEX*16 EPS(0:4)
      COMMON /FOTON/EPS

C Parametros eletrofracos
      REAL*8 ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      COMMON /ELECTROWEAK/ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      REAL*8 C1D,C1S,C2D,C2S
      COMMON /CVCA1/C1D,C1S,C2D,C2S

C Diversos
      REAL*8 DK,DOT4
      COMPLEX*16 ZDOTR,ZDOTZ,ZZZR,ZZRR,DZ
      REAL*8 AUX1,AUX2
      REAL*8 ALAMBDA
      COMMON/ACOPL3/ALAMBDA

      INTEGER I
      REAL*8 DIP173,DIP183,DIP174,DIP184


C                               Inicio
C                              --------


c M21-M22
      AUX1=(C2D*DK(HA,+1)+C2S*DK(HA,-1))*
     #     (C2D*DK(HB,+1)+C2S*DK(HB,-1))*
     #     (DOT4(QB3,QB3)-DOT4(Q2A,Q2A))/AMZ2/
     #     DZ(DOT4(Q2A,Q2A))/DZ(DOT4(QB3,QB3))
	
      AUX2=(C2D*DK(HA,+1)+C2S*DK(HA,-1))*
     #     (C2D*DK(HB,+1)+C2S*DK(HB,-1))*
     #     (DOT4(QB2,QB2)-DOT4(Q3A,Q3A))/AMZ2/
     #     DZ(DOT4(Q3A,Q3A))/DZ(DOT4(QB2,QB2))
	
      DIP173=1.D0/(1.D0-(DOT4(Q2A,Q2A)+DOT4(QB3,QB3))/ALAMBDA**2)**3
      DIP183=1.D0/(1.D0-(DOT4(Q3A,Q3A)+DOT4(QB2,QB2))/ALAMBDA**2)**3
      DIP174=1.D0/(1.D0-(DOT4(Q2A,Q2A)+DOT4(QB3,QB3))/ALAMBDA**2)**4
      DIP184=1.D0/(1.D0-(DOT4(Q3A,Q3A)+DOT4(QB2,QB2))/ALAMBDA**2)**4

C H1Z
      AMP17(1)=AUX1*DIP173*
     #         (ZDOTZ(ETA2A,EPS)*ZDOTR(ETA3B,P1)-
     #          ZDOTZ(ETA3B,EPS)*ZDOTR(ETA2A,P1))
      AMP18(1)=AUX2*DIP183*
     #         (ZDOTZ(ETA3A,EPS)*ZDOTR(ETA2B,P1)-
     #          ZDOTZ(ETA2B,EPS)*ZDOTR(ETA3A,P1))

C H2Z
      AMP17(2)=AUX1/AMZ2*DIP174*
     #          (ZDOTR(ETA2A,QB3)*
     #           (DOT4(QB3,P1)*ZDOTZ(ETA3B,EPS)-
     #            ZDOTR(ETA3B,P1)*ZDOTR(EPS,QB3))-
     #           ZDOTR(ETA3B,Q2A)*
     #           (DOT4(Q2A,P1)*ZDOTZ(ETA2A,EPS)-
     #           ZDOTR(ETA2A,P1)*ZDOTR(EPS,Q2A)))

      AMP18(2)=AUX2/AMZ2*DIP184*
     #         (ZDOTR(ETA3A,QB2)*
     #          (DOT4(QB2,P1)*ZDOTZ(ETA2B,EPS)-
     #           ZDOTR(ETA2B,P1)*ZDOTR(EPS,QB2))-
     #          ZDOTR(ETA2B,Q3A)*
     #          (DOT4(Q3A,P1)*ZDOTZ(ETA3A,EPS)-
     #           ZDOTR(ETA3A,P1)*ZDOTR(EPS,Q3A)))

C H3Z
      AMP17(3)=AUX1*DIP173*ZZZR(ETA3B,ETA2A,EPS,P1)
      AMP18(3)=AUX2*DIP183*ZZZR(ETA2B,ETA3A,EPS,P1)

C H4Z
      AMP17(4)=AUX1/AMZ2*DIP174*
     #         (ZDOTR(ETA2A,QB3)*ZZRR(ETA3B,EPS,QB3,P1)-
     #          ZDOTR(ETA3B,Q2A)*ZZRR(ETA2A,EPS,Q2A,P1))
      AMP18(4)=AUX2/AMZ2*DIP184*
     #         (ZDOTR(ETA3A,QB2)*ZZRR(ETA2B,EPS,QB2,P1)-
     #          ZDOTR(ETA2B,Q3A)*ZZRR(ETA3A,EPS,Q3A,P1))

      DO I=1,4
        AMP(I)=AMP17(I)-AMP18(I)
      END DO

      RETURN
      END
C **********************************************************************
C **********************************************************************

C             Calculo das amplitudes com o vertice ZgZ

C **********************************************************************

      SUBROUTINE AMPLTD2(HA,HB,H1,H2,H3,AMP)

C                     Declaracao de variaveis
C                    -------------------------
      IMPLICIT NONE
	
C Amplitudes
      COMPLEX*16 AMP(4),AMP19(4),AMP20(4),AMP21(4),AMP22(4)

C Helicidades
      INTEGER HA,HB,H1,H2,H3

C Momentos iniciais, finais e internos
      REAL*8 PA(0:4),PB(0:4),P1(0:4),P2(0:4),P3(0:4)
      COMMON/MOMENTA/PA,PB,P1,P2,P3
      REAL*8 QA1(0:4),QB3(0:4),QB1(0:4),Q2A(0:4),K12(0:4),K13(0:4)
      REAL*8 QB2(0:4),Q3A(0:4),QA2(0:4),Q3B(0:4),QA3(0:4),Q2B(0:4)
      COMMON /AUXMOM/QA1,QB3,QB1,Q2A,K12,K13,QB2,Q3A,QA2,Q3B,QA3,Q2B

C Bilinears; vetor de polarizacao
      COMPLEX*16 ETA2A(0:4),ETA3B(0:4),ETA2B(0:4),ETA3A(0:4)
      COMMON /BILINEAR/ETA2A,ETA3B,ETA2B,ETA3A
      COMPLEX*16 EPS(0:4)
      COMMON /FOTON/EPS

C Parametros eletrofracos
      REAL*8 ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      COMMON /ELECTROWEAK/ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      REAL*8 C1D,C1S,C2D,C2S
      COMMON /CVCA1/C1D,C1S,C2D,C2S

C Diversos
      REAL*8 DK,DOT4
      COMPLEX*16 ZDOTR,ZDOTZ,ZZZR,ZZRR,DZ
      REAL*8 AUX19,AUX20,AUX21,AUX22
      REAL*8 ALAMBDA
      COMMON/ACOPL3/ALAMBDA

      INTEGER I
      REAL*8 DIP193,DIP203,DIP213,DIP223,DIP194,DIP204,DIP214,DIP224


C                               Inicio
C                              --------


      AUX19=-(C2D*DK(HA,+1)+C2S*DK(HA,-1))/AMZ2/DZ(DOT4(Q2A,Q2A))
      AUX20=-(C2D*DK(HB,+1)+C2S*DK(HB,-1))/AMZ2/DZ(DOT4(QB3,QB3))
      AUX21=-(C2D*DK(HA,+1)+C2S*DK(HA,-1))/AMZ2/DZ(DOT4(Q3A,Q3A))
      AUX22=-(C2D*DK(HB,+1)+C2S*DK(HB,-1))/AMZ2/DZ(DOT4(QB2,QB2))

      DIP193=1.D0/(1.D0-(DOT4(Q2A,Q2A)+DOT4(QB3,QB3))/ALAMBDA**2)**3
      DIP203=DIP193
      DIP213=1.D0/(1.D0-(DOT4(Q3A,Q3A)+DOT4(QB2,QB2))/ALAMBDA**2)**3
      DIP223=DIP213
      DIP194=1.D0/(1.D0-(DOT4(Q2A,Q2A)+DOT4(QB3,QB3))/ALAMBDA**2)**4
      DIP204=DIP194
      DIP214=1.D0/(1.D0-(DOT4(Q3A,Q3A)+DOT4(QB2,QB2))/ALAMBDA**2)**4
      DIP224=DIP214

C H1g
      AMP19(1)=AUX19*DIP193*
     #         (ZDOTZ(ETA2A,EPS)*ZDOTR(ETA3B,P1)-
     #          ZDOTZ(ETA3B,EPS)*ZDOTR(ETA2A,P1))
      AMP20(1)=-AUX20*DIP203*
     #         (ZDOTZ(ETA2A,EPS)*ZDOTR(ETA3B,P1)-
     #          ZDOTZ(ETA3B,EPS)*ZDOTR(ETA2A,P1))
      AMP21(1)=AUX21*DIP213*
     #         (ZDOTZ(ETA3A,EPS)*ZDOTR(ETA2B,P1)-
     #          ZDOTZ(ETA2B,EPS)*ZDOTR(ETA3A,P1))
      AMP22(1)=-AUX22*DIP223*
     #         (ZDOTZ(ETA3A,EPS)*ZDOTR(ETA2B,P1)-
     #          ZDOTZ(ETA2B,EPS)*ZDOTR(ETA3A,P1))

C H2g
      AMP19(2)=AUX19/AMZ2*DIP194*
     #         ZDOTR(ETA2A,QB3)*
     #          (DOT4(QB3,P1)*ZDOTZ(ETA3B,EPS)-
     #           ZDOTR(ETA3B,P1)*ZDOTR(EPS,QB3))
      AMP20(2)=AUX20/AMZ2*DIP204*
     #         ZDOTR(ETA3B,Q2A)*
     #          (DOT4(Q2A,P1)*ZDOTZ(ETA2A,EPS)-
     #           ZDOTR(ETA2A,P1)*ZDOTR(EPS,Q2A))
      AMP21(2)=AUX21/AMZ2*DIP214*
     #         ZDOTR(ETA3A,QB2)*
     #          (DOT4(QB2,P1)*ZDOTZ(ETA2B,EPS)-
     #           ZDOTR(ETA2B,P1)*ZDOTR(EPS,QB2))
      AMP22(2)=AUX22/AMZ2*DIP224*
     #         ZDOTR(ETA2B,Q3A)*
     #          (DOT4(Q3A,P1)*ZDOTZ(ETA3A,EPS)-
     #           ZDOTR(ETA3A,P1)*ZDOTR(EPS,Q3A))

C H3g
      AMP19(3)=AUX19*DIP193*ZZZR(ETA3B,ETA2A,EPS,P1)
      AMP20(3)=-AUX20*DIP203*ZZZR(ETA3B,ETA2A,EPS,P1)
      AMP21(3)=AUX21*DIP213*ZZZR(ETA2B,ETA3A,EPS,P1)
      AMP22(3)=-AUX22*DIP223*ZZZR(ETA2B,ETA3A,EPS,P1)

C H4g
      AMP19(4)=AUX19/AMZ2*DIP194*
     #         ZDOTR(ETA2A,QB3)*ZZRR(ETA3B,EPS,QB3,P1)
      AMP20(4)=AUX20/AMZ2*DIP204*
     #         ZDOTR(ETA3B,Q2A)*ZZRR(ETA2A,EPS,Q2A,P1)
      AMP21(4)=AUX21/AMZ2*DIP214*
     #         ZDOTR(ETA3A,QB2)*ZZRR(ETA2B,EPS,QB2,P1)
      AMP22(4)=AUX22/AMZ2*DIP224*
     #         ZDOTR(ETA2B,Q3A)*ZZRR(ETA3A,EPS,Q3A,P1)

      DO I=1,4
        AMP(I)=AMP19(I)+AMP20(I)-AMP21(I)-AMP22(I)
      END DO

      RETURN
      END
C **********************************************************************
C **********************************************************************
      SUBROUTINE UBGU(P2,IH2,P1,IH1,ZETA)

C  zeta^mu = ubar(p2,ih2).(gamma^mu).u(p1,ih1) , mu=0,1,2,3
	
      COMPLEX*16 EX1P,EX1M,EX2P,EX2M
      COMPLEX*16 UB(0:3),U(0:3),ZETA(0:4)
      REAL*8 A1,B1,A2,B2,P1(0:4),P2(0:4)
      REAL*8 PMOD1,PMOD2,CT1,CT2,ST1,ST2,CF1,CF2,SF1,SF2,PT1,PT2
      INTEGER IH1,IH2

      PT1=DSQRT(P1(1)*P1(1)+P1(2)*P1(2))
      IF (PT1.EQ.0.D0) THEN
         PMOD1=DABS(P1(3))         
         EX1P=(1.D0,0.D0)
         EX1M=(1.D0,0.D0)
         IF (P1(3).LT.0.D0) THEN
      	CT1=0.D0
      	ST1=1.D0
         ELSE
      	CT1=1.D0
      	ST1=0.D0
         END IF
      ELSE
         PMOD1=DSQRT(PT1*PT1+P1(3)*P1(3))
         CT1=DSQRT(0.5D0*(1.D0+P1(3)/PMOD1))
         ST1=DSQRT(1.D0-CT1*CT1)
         CF1=P1(1)/PT1
         SF1=P1(2)/PT1
         EX1P=CF1+(0.D0,1.D0)*SF1
         EX1M=CF1-(0.D0,1.D0)*SF1
      END IF
      A1=DSQRT(P1(0)+P1(4))
      B1=PMOD1/(A1*A1)
      IF (IH1.EQ.-1) THEN
         U(0)=-EX1M*ST1
         U(1)=CT1
         U(2)=B1*EX1M*ST1
         U(3)=-B1*CT1
      ELSE IF (IH1.EQ.1) THEN
         U(0)=CT1
         U(1)=EX1P*ST1
         U(2)=B1*CT1
         U(3)=B1*EX1P*ST1
      END IF
	
      PT2=DSQRT(P2(1)*P2(1)+P2(2)*P2(2))
      IF (PT2.EQ.0.D0) THEN
         PMOD2=DABS(P2(3))
         EX2P=(1.D0,0.D0)
         EX2M=(1.D0,0.D0)
         IF (P2(3).LT.0.D0) THEN
      	CT2=0.D0
      	ST2=1.D0
         ELSE
      	CT2=1.D0
      	ST2=0.D0
         END IF
      ELSE
         PMOD2=DSQRT(PT2*PT2+P2(3)*P2(3))
         CT2=DSQRT(0.5D0*(1.D0+P2(3)/PMOD2))
         ST2=DSQRT(1.D0-CT2*CT2)
         CF2=P2(1)/PT2
         SF2=P2(2)/PT2
         EX2P=CF2+(0.D0,1.D0)*SF2
         EX2M=CF2-(0.D0,1.D0)*SF2
      END IF
      A2=DSQRT(P2(0)+P2(4))
      B2=PMOD2/(A2*A2)
      IF (IH2.EQ.-1) THEN
         UB(0)=-EX2P*ST2
         UB(1)=CT2
         UB(2)=-B2*EX2P*ST2
         UB(3)=B2*CT2
      ELSE IF (IH2.EQ.1) THEN
         UB(0)=CT2
         UB(1)=EX2M*ST2
         UB(2)=-B2*CT2
         UB(3)=-B2*EX2M*ST2
      END IF      

      ZETA(0)=A1*A2*(UB(0)*U(0)+UB(1)*U(1)-UB(2)*U(2)-UB(3)*U(3))
      ZETA(1)=A1*A2*(UB(0)*U(3)+UB(1)*U(2)-UB(2)*U(1)-UB(3)*U(0))
      ZETA(2)=-UB(0)*U(3)+UB(1)*U(2)+UB(2)*U(1)-UB(3)*U(0)
      ZETA(2)=(0.D0,1.D0)*A1*A2*ZETA(2)
      ZETA(3)=A1*A2*(UB(0)*U(2)-UB(1)*U(3)-UB(2)*U(0)+UB(3)*U(1))

      RETURN
      END
C**********************************************************************
C  Vetor de polarizacao do foton real (complexo-conjugado)
      SUBROUTINE GEPS(P1,H1,EPS)

      IMPLICIT NONE
	
      REAL*8 P1(0:4)
      INTEGER H1

C  Photon's polarization vector
      COMPLEX*16 EPS1(0:4),EPS2(0:4)
      COMPLEX*16 EPS(0:4)

      REAL*8 Q0,Q1,Q01

      REAL*8 DOT3
      INTEGER I

      Q0=DSQRT(DOT3(P1,P1))
      Q1=DSQRT(P1(1)*P1(1)+P1(2)*P1(2))
      Q01=1./(Q0*Q1)

      EPS1(0)=0.
      EPS1(1)=Q01*(P1(1)*P1(3))
      EPS1(2)=Q01*(P1(2)*P1(3))
      EPS1(3)=Q01*(-Q1*Q1)

      EPS2(0)=0.
      EPS2(1)=(1./Q1)*(-P1(2))
      EPS2(2)=(1./Q1)*(P1(1))
      EPS2(3)=0.

      DO I=0,3
       IF (H1.EQ.+1) THEN
        EPS(I)=DSQRT(.5D0)*(-EPS1(I)-(0.,1.)*EPS2(I))
        EPS(I)=CONJG(EPS(I))
       ELSE    
        EPS(I)=DSQRT(.5D0)*(+EPS1(I)-(0.,1.)*EPS2(I))
        EPS(I)=CONJG(EPS(I))
       END IF         
      END DO

      RETURN
      END      
C **********************************************************************
C Quadrado do modulo de um numero complexo
      REAL*8 FUNCTION RMOD2(Z)
	
      IMPLICIT NONE

      COMPLEX*16 Z

      RMOD2=DREAL(Z)*DREAL(Z)+DIMAG(Z)*DIMAG(Z)

      RETURN
      END
C **********************************************************************
C Delta de Kronecker
      REAL*8 FUNCTION DK(I,J)
	
      IMPLICIT NONE
	
      INTEGER I,J
	
      DK=0.5D0*IABS(I+J)

      RETURN
      END
C **********************************************************************
C Produto escalar de dois vetores
      DOUBLE PRECISION FUNCTION DOT3(P1,P2)

      IMPLICIT NONE

      REAL*8 P1(0:4),P2(0:4)
	
      DOT3=P1(1)*P2(1)+P1(2)*P2(2)+P1(3)*P2(3)
	
      RETURN
      END
C***********************************************************************
C Produto escalar de dois 4-vetores reais
      DOUBLE PRECISION FUNCTION DOT4(P1,P2)

      IMPLICIT NONE
	
      REAL*8 P1(0:4),P2(0:4)
	
      DOT4=P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)
	
      RETURN
      END
C **********************************************************************
C Breit-Wigner do Z
      COMPLEX*16 FUNCTION DZ(Q2)
	
      IMPLICIT NONE
      REAL*8 ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      COMMON /ELECTROWEAK/ALPHA,CTW,CTW2,STW,STW2,AMZ2,GZMZ
      REAL*8 Q2
	
      DZ = (Q2-AMZ2+(0.D0,1.D0)*GZMZ)
	
      RETURN
      END      
C **********************************************************************
C Produto escalar de dois 4-vetores complexos
      COMPLEX*16 FUNCTION ZDOTZ(ZVA,ZVB)
	
      IMPLICIT NONE
      COMPLEX*16 ZVA(0:4),ZVB(0:4)

      ZDOTZ=ZVA(0)*ZVB(0)-ZVA(1)*ZVB(1)-ZVA(2)*ZVB(2)-ZVA(3)*ZVB(3)

      RETURN
      END
C **********************************************************************
C Produto escalar de um 4-vetor complexo com um 4-vetor real
      COMPLEX*16 FUNCTION ZDOTR(ZVA,VB)
	
      IMPLICIT NONE
      COMPLEX*16 ZVA(0:4)
      REAL*8 VB(0:4)

      ZDOTR=ZVA(0)*VB(0)-ZVA(1)*VB(1)-ZVA(2)*VB(2)-ZVA(3)*VB(3)

      RETURN
      END
C **********************************************************************
C Contracao de Levi-Civita com 2 4-vetores complexos e 2 4-vetores reais
      COMPLEX*16 FUNCTION ZZRR(ZVA,ZVB,VC,VD)

      IMPLICIT NONE
      COMPLEX*16 ZVA(0:4),ZVB(0:4)
      REAL*8 VC(0:4),VD(0:4)
      COMPLEX*16 Z1,Z2,Z3,Z4,Z5,Z6

      Z1=(ZVA(0)*ZVB(1)-ZVA(1)*ZVB(0))*(VC(2)*VD(3)-VC(3)*VD(2))
      Z2=(ZVA(0)*ZVB(3)-ZVA(3)*ZVB(0))*(VC(1)*VD(2)-VC(2)*VD(1))
      Z3=(ZVA(0)*ZVB(2)-ZVA(2)*ZVB(0))*(VC(3)*VD(1)-VC(1)*VD(3))
      Z4=(ZVA(1)*ZVB(2)-ZVA(2)*ZVB(1))*(VC(0)*VD(3)-VC(3)*VD(0))
      Z5=(ZVA(1)*ZVB(3)-ZVA(3)*ZVB(1))*(VC(2)*VD(0)-VC(0)*VD(2))
      Z6=(ZVA(2)*ZVB(3)-ZVA(3)*ZVB(2))*(VC(0)*VD(1)-VC(1)*VD(0))

      ZZRR=Z1+Z2+Z3+Z4+Z5+Z6

      RETURN
      END
C **********************************************************************
C Contracao de Levi-Civita com 3 4-vetores complexos e 1 4-vetor reaL
      COMPLEX*16 FUNCTION ZZZR(ZVA,ZVB,ZVC,VD)

      IMPLICIT NONE

      COMPLEX*16 ZVA(0:4),ZVB(0:4),ZVC(0:4)
      REAL*8 VD(0:4)
      COMPLEX*16 Z1,Z2,Z3,Z4,Z5,Z6

      Z1=(ZVA(0)*ZVB(1)-ZVA(1)*ZVB(0))*(ZVC(2)*VD(3)-ZVC(3)*VD(2))
      Z2=(ZVA(0)*ZVB(3)-ZVA(3)*ZVB(0))*(ZVC(1)*VD(2)-ZVC(2)*VD(1))
      Z3=(ZVA(0)*ZVB(2)-ZVA(2)*ZVB(0))*(ZVC(3)*VD(1)-ZVC(1)*VD(3))
      Z4=(ZVA(1)*ZVB(2)-ZVA(2)*ZVB(1))*(ZVC(0)*VD(3)-ZVC(3)*VD(0))
      Z5=(ZVA(1)*ZVB(3)-ZVA(3)*ZVB(1))*(ZVC(2)*VD(0)-ZVC(0)*VD(2))
      Z6=(ZVA(2)*ZVB(3)-ZVA(3)*ZVB(2))*(ZVC(0)*VD(1)-ZVC(1)*VD(0))

      ZZZR=Z1+Z2+Z3+Z4+Z5+Z6

      RETURN
      END
C **********************************************************************
C Gerador de numeros aleatorios
      REAL*8 FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 13]2Y_213.
C *************************************************************************
      SUBROUTINE STORE (VBOUND,NBIN,SUM,X,Y)

C  Given a monotonic (either increasing or decreasing) array VBOUND
C  of length NBIN+1 and a value X, this routine finds the bin index J 
C  corresponding to X and adds Y to the J-th component of an array SUM
C  of length NBIN
	 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VBOUND(NBIN+1),SUM(NBIN)
	
      N=NBIN+1
      JL=0
      JU=N+1
 10   IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((VBOUND(N).GT.VBOUND(1)).EQV.(X.GT.VBOUND(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      IF (J.EQ.0.OR.J.EQ.N) PRINT*,'OUT OF RANGE'
      SUM(J)=SUM(J)+Y

      RETURN
      END
C***********************************************************************
      SUBROUTINE FPLOT (KUNIT,NBIN,VBOUND,SUM)
	
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VBOUND(NBIN+1),SUM(NBIN)
	
      DO I=1,NBIN
       XMED=0.5D0*(VBOUND(I)+VBOUND(I+1))
       WRITE(KUNIT,*) XMED,SUM(I)
      END DO
	
      RETURN
      END      
C***********************************************************************
      SUBROUTINE AREAD(NBIN,BWIDTH,Y,AREA)
	
C Calculo da area definida por uma curva obtida a partir de pontos 
C "experimentais".


      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION Y(NBIN)

      AREA=3.D0*(Y(1)+Y(NBIN))/8.D0 + 7.D0*(Y(2)+Y(NBIN-1))/6.D0 + 
     #     23.D0*(Y(3)+Y(NBIN-2))/24.D0

      DO I=4,NBIN-3
      AREA=AREA+Y(I)
      END DO

      AREA=AREA*BWIDTH

      RETURN
      END
C***********************************************************************
      SUBROUTINE ACOPL(NHX,NHY,DHX,DHY,AHX,AHY)

      IMPLICIT NONE

      INTEGER NHX,NHY
      REAL*8 DHX,DHY
      REAL*8 AHX(NHX),AHY(NHY)
      INTEGER I

      DO I=1,NHX
        AHX(I)=0.D0
      END DO
      DO I=1,NHY
        AHY(I)=0.D0
      END DO

      DHX=DHX*2.D0
      DHY=DHY*2.D0

      DO I=1,NHX
        AHX(I)=-DHX/2.D0+(I-1)*DHX/(NHX-1)
      END DO
      DO I=1,NHY
        AHY(I)=-DHY/2.D0+(I-1)*DHY/(NHY-1)
      END DO

      RETURN
      END
C *************************************************************************
      SUBROUTINE STORE1(VBOUND,NBIN,SUM,X,Y,NX,NY)

C  Given a monotonic (either increasing or decreasing) array VBOUND
C  of length NBIN+1 and a value X, this routine finds the bin index J 
C  corresponding to X and adds Y to the J-th component of an array SUM
C  of length NBIN
	 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NX,NY,NBIN
      REAL*8 X
      REAL*8 VBOUND(NBIN+1),SUM(NBIN,NX,NY),Y(NX,NY)

      N=NBIN+1
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((VBOUND(N).GT.VBOUND(1)).EQV.(X.GT.VBOUND(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
        GO TO 10
      ENDIF
      J=JL
      IF (J.EQ.0.OR.J.EQ.N) PRINT*,'OUT OF RANGE'
      DO IX=1,NX
       DO IY=1,NY
        SUM(J,IX,IY)=SUM(J,IX,IY)+Y(IX,IY)
       END DO
      END DO

      RETURN
      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE AREAD1(NBIN,BWIDTH,Y,AREA,NX,NY)
	
C Calculo da area definida por uma curva obtida a partir de pontos 
C "experimentais".


      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION Y(NBIN,NX,NY),AREA(NX,NY)

      DO I=1,NX
       DO J=1,NY

        AREA(I,J)=3.D0*(Y(1,I,J)+Y(NBIN,I,J))/8.D0 +
     #            7.D0*(Y(2,I,J)+Y(NBIN-1,I,J))/6.D0 +
     #            23.D0*(Y(3,I,J)+Y(NBIN-2,I,J))/24.D0

        DO K=4,NBIN-3
         AREA(I,J)=AREA(I,J)+Y(K,I,J)
        END DO

        AREA(I,J)=AREA(I,J)*BWIDTH
	  
       END DO
      END DO

      RETURN
      END
C***********************************************************************
