!     ------------------------------------------------------------------
!     |     DATE           PROGRAMMER           DESCRIPTION OF CHANGE  |  
!     |     ====           ==========            ===================== |
!     |  10/09/2021      MÁRIO RUI ARRUDA                              |
!     |  IST LISBON                                                    |
!     ------------------------------------------------------------------

!     ---------UMAT ONLY FOR PLANE ELEMENST WITH TSAI-WU DAMAGE---------

! When using this UMAT in papers, always cite these works https://doi.org/10.1016/j.jcomc.2021.100122
! and always provide credit to the original authors

!     ------------------------------------------------------------------
!     ----------ABAQUS INPUT VARIABLES IN UMAT SUBROUTINE---------------
!     ------------------------------------------------------------------

subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, &
           ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, & 
           nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, &
           npt, layer, kspt, kstep, kinc)

IMPLICIT NONE
 

!     ------------------------------------------------------------------
!     -------------------ABAQUS DIMENSION VARIABLES---------------------
!     ------------------------------------------------------------------
 
!     -----------------ABAQUS CHARACTER VARIABLES----------------------
character(len=80) :: cmname 

!     -------------------ABAQUS INTEGER VARIABLES-----------------------  
integer :: ntens,ndi,nshr,nstatv,nprops,noel,npt,kspt,kstep,kinc,nprecd,layer

!     -------------------ABAQUS REAL VARIABLES--------------------------
real(kind=8) :: celent,sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt 

!     -------------------ABAQUS ARRAY VARIABLES-------------------------
real(kind=8) :: stress(ntens),statev(nstatv),ddsdde(ntens,ntens),&
                ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),&
                predef(1),dpred(1),props(nprops),coords(3),drot(3,3),&
                dfgrd0(3,3),dfgrd1(3,3),time(2)

					 
!     ------------------------------------------------------------------    
!     -----------------DECLARATION OF VARIABLES-------------------------
!     ------------------------------------------------------------------

integer :: i,j,k,l,m,n,kiter,ktotal

real(kind=8) :: E1,E2,G,v12,v21,Xt,Xc,Yt,Yc,St,Sl,TWold,Tw,Twf,Twm,aft,afc,amt,amc,alpha,&
								afto,afco,amto,amco,Gft,Gfc,Gmt,Gmc,eta,dmax,dfto,dfco,dmto,dmco,dso,&
								dvft,dvfc,dvmt,dvmc,dvfto,dvfco,dvmto,dvmco,dvso,seqft0,seqfc0,seqmt0,seqmc0,&
								df,dvf,dm,dvm,Dcoef,catLc,dft,dfc,dmt,dmc,ds,dvs,F1,F2,F11,F22,F12,Bxy,F66,&
								ueqft,ueqft0,ueqftu,ueqfc,ueqfc0,ueqfcu,ueqmt,ueqmt0,ueqmtu,ueqmc,ueqmc0,ueqmcu,&
                ueqftp,ueqfcp,ueqmtp,ueqmcp,pft,pfc,pmt,pmc,pms,ashear,ashearo,tau0,gamma0,&
                ueqg,ueqg0,ueqgu,ueqgp

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4

!     -------------INITIATIONS OF ARRAYS--------------------------------
real(kind=8) :: eij(ntens),eijbp(ntens),eijbn(ntens),eijo(ntens),eijbpo(ntens),eijbno(ntens),&
                sij(ntens),sijbp(ntens),sijbn(ntens),sijo(ntens),sijbpo(ntens),sijbno(ntens),&
                sije(ntens),sijebp(ntens),sijebn(ntens),sijeo(ntens),sijebpo(ntens),sijebno(ntens),&
                ddsddee(ntens,ntens)
                
!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3-------------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0, FOUR=4.d0,&
                           HALF=0.5d0,TOL0=1.d-8, TOLseq=0.1d0			
						
!     --------------ELASTIC AND MECHANICAL PROPERTIES-------------------

E1=props(1)     ! LONGITUDINAL ELASTIC MODULUS
E2=props(2)     ! TRANSVERSAL ELASTIC MODULUS
G=props(3)      ! SHEAR ELASTIC MODULUS
v12=props(4)    ! POISSON COEFICIENT
v21=v12*E2/E1   ! SYMMETRIC POISSON COEFICIENT

Xt=props(5)     ! LONGITUDINAL TENSILE RESISTENCE
Xc=props(6)     ! LONGITUDINAL COMPRESSIVE RESISTENCE
Yt=props(7)     ! TRANSVERSAL TENSILE RESISTENCE
Yc=props(8)     ! TRANSVERSAL COMPRESSIVE RESISTENCE
Sl=props(9)     ! LONGITUDINAL SHEAR RESISTENCE
St=props(10)    ! TRANSVERSAL SHEAR RESISTENCE
alpha=props(11)

Gft=props(12)   ! FRACTURE ENERGY FOR FIBRE TENSION
Gfc=props(13)   ! FRACTURE ENERGY FOR FIBRE COMPRESSION
Gmt=props(14)   ! FRACTURE ENERGY FOR MATRIX TENSION
Gmc=props(15)   ! FRACTURE ENERGY FOR MATRIX COMPRESSION

eta=props(16)   ! VISCOUS REGULARIZATION COEFICIENT
dmax=props(17)  ! MAXIMUM ALLOWED DAMAGE FOR CONVERGENCY

pft=props(18)   ! RESIDUAL LONGITUDINAL TENSILE RESISTENCE
pfc=props(19)   ! RESIDUAL LONGITUDINAL COMPRESSIVE RESISTENCE
pmt=props(20)   ! RESIDUAL TRANSVERSAL TENSILE RESISTENCE
pmc=props(21)   ! RESIDUAL TRANSVERSAL COMPRESSIVE RESISTENCE


!     ---------------STATE FIELD VARIABLES FOR ABAQUS-------------------
if (kinc == 1) then 
	do i=1, nstatv
		statev(i)=ZERO  ! VARIABLES INITIATION (if not=0, depending on compiler LINUX/WINDOWS)
  end do
end if

dfto=statev(1)      ! FIBRE TENSION DAMAGE FROM PREVIOUS INCREMENT
dfco=statev(2)      ! FIBRE COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dmto=statev(3)      ! MATRIX TENSION DAMAGE FROM PREVIOUS INCREMENT
dmco=statev(4)      ! MATRIX COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dso=statev(5)       ! SHEAR DAMAGE FROM PREVIOUS INCREMENT

afto=statev(6)      ! ALPHA TENSION FIBER CONSTANT
afco=statev(7)      ! ALPHA COMPRESSION FIBER CONSTANT
amto=statev(8)      ! ALPHA TENSION MATRIX CONSTANT
amco=statev(9)      ! ALPHA COMPRESSION MATRIX CONSTANT
TWold=statev(10)    ! TSAI-WU FAILURE FUNCTION

dvfto=statev(11)    ! VISCOUS FIBRE TENSION DAMAGE FROM PREVIOUS INCREMENT
dvfco=statev(12)    ! VISCOUS FIBRE COMRESSION DAMAGE FROM PREVIOUS INCREMENT
dvmto=statev(13)    ! VISCOUS MATRIX TENSION DAMAGE FROM PREVIOUS INCREMENT
dvmco=statev(14)    ! VISCOUS MATRIX COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dvso=statev(15)     ! VISCOUS SHEAR DAMAGE FROM PREVIOUS INCREMENT

seqft0=statev(16)   ! EQUIVALENTE FIBRE TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqfc0=statev(17)   ! EQUIVALENTE FIBRE COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqmt0=statev(18)   ! EQUIVALENTE MATRIX TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqmc0=statev(19)   ! EQUIVALENTE MATRIX COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT

ueqft0=statev(20)   ! EQUIVALENTE FIBRE TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqfc0=statev(21)   ! EQUIVALENTE FIBRE COMPRESSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqmt0=statev(22)   ! EQUIVALENTE MATRIX TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqmc0=statev(23)   ! EQUIVALENTE MATRIX COMPRESSON DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT

sijeo(1)=statev(24) ! EFFECTIVE S11 STRESS AT THE BEGGING OF THE INCREMENT
sijeo(2)=statev(25) ! EFFECTIVE S22 STRESS AT THE BEGGING OF THE INCREMENT
sijeo(3)=statev(26) ! EFFECTIVE S12 STRESS AT THE BEGGING OF THE INCREMENT

ashearo=statev(27)  ! ALPHA SHEAR CONSTANT
tau0=statev(28)     ! INITIAL SHEAR DAMAGE STRESS
gamma0=statev(29)   ! INITIAL SHEAR DAMAGE DISTORTION


!     ------------------------------------------------------------------
!     ---------FIELD PREDICTOR FOR VARIABLES AND STIFFNESS -------------
!     ------------------------------------------------------------------

!     --------------FINAL STRAIN INCREMENT------------------------------
do i=1, ntens
  eij(i)=stran(i)+dstran(i) ! TOTAL STRAIN IN THE BEGGIND OF THE INCREMENT
  eijo(i)=stran(i)          ! STRAIN FROM PREVIOUS INCREMENT
  sijo(i)=stress(i)         ! STRESS FROM PREVIOUS INCREMENT
end do

! !     ----STABILITY CONDITION DURING EXPLICIT MATERIAL INTEGRATION------
do i=1,ntens
	if (abs(dstran(i)) > (Xt/E1)) then
		pnewdt=0.5
	end if
end do

!     -------VERIFICATION OF FIBRE AND MATRIX TENSION/COMPRESSION-------

df=ONE-(ONE-dfto)*(ONE-dfco)
dvf=ONE-(ONE-dvfto)*(ONE-dvfco)
dm=ONE-(ONE-dmto)*(ONE-dmco)
dvm=ONE-(ONE-dvmto)*(ONE-dvmco)

ds=dso
aft=afto
afc=afco
amt=amto
amc=amco
ashear=ashearo

!     -----SECANT DAMAGED MATRIX AT THE BEGINNING OF THE INCREMENT------
call calc_stiff_assembler(ntens,ddsdde,Dcoef,E1,E2,G,df,dm,ds,v12,v21)

!     -----SECANT EFFECTIVE MATRIX AT THE BEGINNING OF THE INCREMENT----
call calc_efect_stiff_assembler(ntens,ddsddee,Dcoef,E1,E2,G,df,dm,ds,v12,v21)

!    -----------MULTIPLICATION OF SECANT MATRIX WITH STRAIN-------------
call calc_s_ke_multiply(ndi,ntens,sij,eij,ddsdde)
call calc_s_ke_multiply(ndi,ntens,sije,eij,ddsddee)


!     ------------------------------------------------------------------
!     ------------------CREATING MACAULY BRACKETS-----------------------
!     ------------------------------------------------------------------
call calc_macaulay_bracket(ndi,ntens,eij,eijbp,eijbn)        ! BRACKETS FOR STRAIN
call calc_macaulay_bracket(ndi,ntens,sij,sijbp,sijbn)        ! BRACKETS FOR STRESS
call calc_macaulay_bracket(ndi,ntens,sije,sijebp,sijebn)     ! BRACKETS FOR EFECTIVE STRESS
call calc_macaulay_bracket(ndi,ntens,eijo,eijbpo,eijbno)     ! BRACKETS FOR PREVIOUS STRAIN
call calc_macaulay_bracket(ndi,ntens,sijo,sijbpo,sijbno)     ! BRACKETS FOR PREVIOUS STRESS
call calc_macaulay_bracket(ndi,ntens,sijeo,sijebpo,sijebno)  ! BRACKETS FOR PREVIOUS EFECTIVE STRESS


!     ------------------------------------------------------------------
!     --------------FAILURE CRITERIA VERIFICATION-----------------------
!     ------------------------------------------------------------------

!     -------------CARACTERISTIC LENGTH FROM ABAQUS---------------------
catLc=celent
F11=ONE/(Xt*Xc)
F22=ONE/(Yt*Yc)
F1=ONE/Xt-ONE/Xc
F2=ONE/Yt-ONE/Yc
F66=ONE/(St*Sl)
Bxy=min(Yt/(ONE-v21),Xt/(ONE-v12))
F12=ONE/(TWO*Bxy**2)*(ONE-Bxy*(F1+F2)-Bxy**2*(F11+F22))
if (F12**2 > (F11*F22)) then
  F12=-HALF*sqrt(F11*F22)
end if

!     --------VERIFICATION OF FIBRE TENSION/COMPRESSION FAILURE---------
TW=F1*sije(1)+F2*sije(2)+TWO*F12*sije(1)*sije(2)+F11*sije(1)**2+F22*sije(2)**2+F66*sije(3)**2

if (sije(1) >= ZERO .and. aft < TOLseq) then
	if (TW >= ONE) then
    aft=(abs(sijebp(1)/Xt)+alpha*abs(sije(3)/St))
    if (aft >= TOLseq) then
      call calc_activation(catLc,F11,F66,F1,sijeo(1),sijeo(3),eijbpo(1),eijo(3),sijbpo(1),sijo(3),ueqft0,seqft0,alpha)
      statev(20)=ueqft0
      statev(16)=seqft0
      aft=max(afto,min(aft,ONE)) !!!! THIS CONDITION alpha<1.0 and delta_alpha>0
      statev(6)=aft
    end if
	end if
else if (sije(1) < ZERO .and. afc < TOLseq) then
	if (TW >= ONE) then
    afc=(abs(sijebn(1)/Xc)+alpha*abs(sije(3)/St))
    if (afc >= TOLseq) then
      call calc_activation(catLc,F11,F66,F1,sijeo(1),sijeo(3),eijbno(1),eijo(3),sijbno(1),sijo(3),ueqfc0,seqfc0,alpha)
      statev(21)=ueqfc0
      statev(17)=seqfc0
      afc=max(afco,min(afc,ONE)) !!!! THIS CONDITION alpha<1.0 and delta_alpha>0
      statev(7)=afc
    end if
	end if
end if

!     --------VERIFICATION OF MATRIX TENSION/COMPRESSION FAILURE---------
if (sije(2) >= ZERO .and. amt < TOLseq) then
	if (TW >= ONE) then
    amt=(abs(sijebp(2)/Yt)+alpha*abs(sije(3)/St))
    if (amt >= TOLseq) then
      call calc_activation(catLc,F22,F66,F2,sijeo(2),sijeo(3),eijbpo(2),eijo(3),sijbpo(2),sijo(3),ueqmt0,seqmt0,alpha)
      statev(22)=ueqmt0
      statev(18)=seqmt0
      amt=max(amto,min(amt,ONE)) !!!! THIS CONDITION alpha<1.0 and delta_alpha>0
      statev(8)=amt
    end if
	end if
else if (sije(2) < ZERO .and. amc < TOLseq) then
	if (TW >= ONE) then
    amc=(abs(sijebn(2)/Yc)+alpha*abs(sije(3)/St))
    if (amc >= TOLseq) then
      call calc_activation(catLc,F22,F66,F2,sijeo(2),sijeo(3),eijbno(2),eijo(3),sijbno(2),sijo(3),ueqmc0,seqmc0,alpha)
		  statev(23)=ueqmc0
		  statev(19)=seqmc0
      amc=max(amco,min(amc,ONE)) !!!! THIS CONDITION alpha<1.0 and delta_alpha>0
      statev(9)=amc
    end if
	end if
end if

if (TW >= ONE) then
  TW=1.001d0
end if
statev(10)=TW


!     ------------------------------------------------------------------
!     -----------------HASHIN DAMAGE EVOLUTION--------------------------
!     ------------------------------------------------------------------

dft=dfto;dfc=dfco;dmt=dmto;dmc=dmco;ds=dso
dvft=dvfto;dvfc=dvfco;dvmt=dvmto;dvmc=dvmco;dvs=dvso

!     -------------FIBER TENSION DAMAGE EVOLUTION-----------------------
if (TW >= ONE .and. aft >= TOLseq) then
	ueqft=sqrt(eijbp(1)**2+alpha*eij(3)**2)*catLc
  aft=(abs(sijebp(1)/Xt)+alpha*abs(sije(3)/St))
  aft=max(afto,min(aft,ONE))
  call calc_damage_lin(ueqft,ueqft0,ueqftu,eta,dft,dfto,dvft,dvfto,dmax,dtime,pft,ueqftp,seqft0,Gft,aft)
  !call calc_damage_exp(ueqft,ueqft0,ueqftu,eta,dft,dfto,dvft,dvfto,dmax,dtime,pft,ueqftp,seqft0,Gft,aft)
  statev(1)=dft
	statev(11)=dvft
end if

!     -------------FIBER COMPRESSION DAMAGE EVOLUTION-------------------
if (TW >= ONE .and. afc >= TOLseq) then
	ueqfc=sqrt(eijbn(1)**2+alpha*eij(3)**2)*catLc
  afc=(abs(sijebn(1)/Xc)+alpha*abs(sije(3)/St))
  afc=max(afco,min(afc,ONE))
  call calc_damage_lin(ueqfc,ueqfc0,ueqfcu,eta,dfc,dfco,dvfc,dvfco,dmax,dtime,pfc,ueqfcp,seqfc0,Gfc,afc)
  !call calc_damage_exp(ueqfc,ueqfc0,ueqfcu,eta,dfc,dfco,dvfc,dvfco,dmax,dtime,pfc,ueqfcp,seqfc0,Gfc,afc)
  statev(2)=dfc
  statev(12)=dvfc
end if

!     -------------MATRIX TENSION DAMAGE EVOLUTION----------------------
if (TW >= ONE .and. amt >= TOLseq) then
	ueqmt=sqrt(eijbp(2)**2+alpha*eij(3)**2)*catLc
  amt=(abs(sijebp(2)/Yt)+alpha*abs(sije(3)/St))
  amt=max(amto,min(amt,ONE))
  call calc_damage_lin(ueqmt,ueqmt0,ueqmtu,eta,dmt,dmto,dvmt,dvmto,dmax,dtime,pmt,ueqmtp,seqmt0,Gmt,amt)
  !call calc_damage_exp(ueqmt,ueqmt0,ueqmtu,eta,dmt,dmto,dvmt,dvmto,dmax,dtime,pmt,ueqmtp,seqmt0,Gmt,amt)
  statev(3)=dmt
	statev(13)=dvmt
end if

!     -----------MATRIX COMPRESSION DAMAGE EVOLUTION--------------------
if (TW >= ONE .and. amc >= TOLseq) then
	ueqmc=sqrt(eijbn(2)**2+alpha*eij(3)**2)*catLc
  amc=(abs(sijebn(2)/Yc)+alpha*abs(sije(3)/St))
  amc=max(amco,min(amc,ONE))
  call calc_damage_lin(ueqmc,ueqmc0,ueqmcu,eta,dmc,dmco,dvmc,dvmco,dmax,dtime,pmc,ueqmcp,seqmc0,Gmc,amc)
  !call calc_damage_exp(ueqmc,ueqmc0,ueqmcu,eta,dmc,dmco,dvmc,dvmco,dmax,dtime,pmc,ueqmcp,seqmc0,Gmc,amc)
  statev(4)=dmc
	statev(14)=dvmc
end if

ds=ONE-sqrt((ONE-dft)*(ONE-dfc)*(ONE-dmt)*(ONE-dmc))
dvs=ONE-sqrt((ONE-dvft)*(ONE-dvfc)*(ONE-dvmt)*(ONE-dvmc))
statev(5)=ds
statev(15)=dvs

statev(6)=aft
statev(7)=afc
statev(8)=amt
statev(9)=amc


!     ------------------------------------------------------------------
!     -----------------FINAL DAMAGE BEHAVIOUR---------------------------
!     ------------------------------------------------------------------

df=ONE-(ONE-dft)*(ONE-dfc)
dvf=ONE-(ONE-dvft)*(ONE-dvfc)
dm=ONE-(ONE-dmt)*(ONE-dmc)
dvm=ONE-(ONE-dvmt)*(ONE-dvmc)


!     ------------------------------------------------------------------
!     -------FIELD CORRECTOR FOR STRESS AND SECANT STIFFNESS -----------
!     ------------------------------------------------------------------

!     ---------VISCOUS SECANT DAMAGED MATRIX IN THE ITERATION-----------
call calc_stiff_assembler(ntens,ddsdde,Dcoef,E1,E2,G,dvf,dvm,dvs,v12,v21)

!     --------VISCOUS EFFECTIVE SECANT MATRIX IN THE ITERATION----------
call calc_efect_stiff_assembler(ntens,ddsddee,Dcoef,E1,E2,G,dvf,dvm,dvs,v12,v21)

!    ---------------------VISCOUS STRESS CORRECTOR----------------------
call calc_s_ke_multiply(ndi,ntens,stress,eij,ddsdde)
call calc_s_ke_multiply(ndi,ntens,sije,eij,ddsddee)

statev(24)=sije(1)
statev(25)=sije(2)
statev(26)=sije(3)

!    -----------------UPDATED DEFORMATION ENERGY------------------------
call calc_half_s_e_multiply(ndi,ntens,sse,stress,sijo,dstran)

return
end


!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||USER SUBROUTINES||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


!     __________________________________________________________________
!     ___________SUBROUTINE FOR MACAULY BRACKET OPERATION_______________
!     __________________________________________________________________

subroutine calc_macaulay_bracket(ndi,ntens,mat,mat_bp,mat_bn)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i
integer, intent(in):: ndi,ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: mat(ntens)
real(kind=8), intent(inout) :: mat_bp(ntens),mat_bn(ntens)

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2---------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
	mat_bp(i)=HALF*(mat(i)+abs(mat(i)))  ! BRACKETS FOR DAMAGE
	mat_bn(i)=HALF*(-mat(i)+abs(mat(i))) ! BRACKETS FOR DAMAGE
end do

return
end subroutine


!     __________________________________________________________________
!     _________SUBROUTINE FOR STRAIN STIFFNESS MULTIPLICATION___________
!     __________________________________________________________________

subroutine calc_s_ke_multiply(ndi,ntens,mat_si,mat_ej,mat_Cij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: mat_Cij(ntens,ntens), mat_ej(ntens)
real(kind=8), intent(inout) :: mat_si(ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
  mat_si(i)=ZERO
  do j=1,ndi ! NORMAL STRAIN
	  mat_si(i)=mat_si(i)+mat_Cij(i,j)*mat_ej(j)
  enddo
enddo
do i=ndi+1,ntens ! SHEAR STRAIN
  mat_si(i)=mat_Cij(i,i)*mat_ej(i) 
enddo

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FOR DEFORMATION ENERGY__________________
!     __________________________________________________________________

subroutine calc_half_s_e_multiply(ndi,ntens,energy,sij,sij_old,deij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: sij(ntens), sij_old(ntens), deij(ntens)
real(kind=8), intent(inout) :: energy

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ntens
  energy=energy+HALF*(sij_old(i)+sij(i))*deij(i)
end do

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE STIFFNESS ASSEMBLER____________________
!     __________________________________________________________________

subroutine calc_stiff_assembler(ntens,Cij,Dvar,E1,E2,G,df,dm,ds,v12,v21)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens

!     ----------------------REAL MATRIX VARIABLES----------------------- 
real(kind=8), intent(in)  :: E1,E2,G,df,dm,ds,v12,v21
real(kind=8), intent(inout) :: Cij(ntens,ntens), Dvar

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

Dvar=ONE-(ONE-df)*(ONE-dm)*v12*v21    
Cij(1,1)=(ONE-df)*E1*ONE/Dvar
Cij(2,2)=(ONE-dm)*E2*ONE/Dvar
Cij(3,3)=(ONE-ds)*G    
Cij(1,2)=(ONE-df)*(ONE-dm)*v21*E1*ONE/Dvar
Cij(2,1)=(ONE-df)*(ONE-dm)*v12*E2*ONE/Dvar
Cij(1,3)=ZERO;Cij(3,1)=ZERO
Cij(2,3)=ZERO;Cij(3,2)=ZERO

return
end subroutine


!     __________________________________________________________________
!     ____________SUBROUTINE EFECTIVE STIFFNESS ASSEMBLER_______________
!     __________________________________________________________________

subroutine calc_efect_stiff_assembler(ntens,Cij,Dvar,E1,E2,G,df,dm,ds,v12,v21)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens

!     ---------------------REAL MATRIX VARIABLES------------------------ 
real(kind=8), intent(in)  :: E1,E2,G,df,dm,ds,v12,v21
real(kind=8), intent(inout) :: Cij(ntens,ntens), Dvar

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

Dvar=ONE-(ONE-df)*(ONE-dm)*v12*v21    
Cij(1,1)=E1*ONE/Dvar
Cij(2,2)=E2*ONE/Dvar
Cij(3,3)=G    
Cij(1,2)=(ONE-dm)*v21*E1*ONE/Dvar
Cij(2,1)=(ONE-df)*v12*E2*ONE/Dvar
Cij(1,3)=ZERO;Cij(3,1)=ZERO
Cij(2,3)=ZERO;Cij(3,2)=ZERO

return
end subroutine


!     __________________________________________________________________
!     ______________SUBROUTINE FOR ACTIVATION FUNCTION__________________
!     __________________________________________________________________

subroutine calc_activation(catLc,Fii,Fjj,Fi,siie,sije,eii,eij,sii,sij,ueq0,seq0,alpha)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8) :: a,b,c,root
real(kind=8), intent(in)  :: catLc,Fii,Fjj,Fi,siie,sije,eii,eij,sii,sij,alpha
real(kind=8), intent(inout) :: ueq0,seq0

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

a=Fii*siie**2+alpha*Fjj*sije**2
b=Fi*siie
c=-ONE
root=abs((-b+sqrt(b**2-FOUR*a*c))/(TWO*a))
ueq0=root*sqrt(eii**2+alpha*eij**2)*catLc
seq0=root**2*(sii*eii+alpha*sij*eij)/(ueq0/catLc)

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE FOR DAMAGE EVOLUTION___________________
!     __________________________________________________________________

subroutine calc_damage_lin(ueq,ueq0,uequ,eta,d,dold,dv,dvold,dmax,dtime,pres,ueqp,seq0,G,al)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8), intent(in)  :: ueq,ueq0,eta,dmax,dvold,dold,dtime,pres,seq0,G,al
real(kind=8), intent(inout) :: d,dv,ueqp,uequ

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

uequ=ueq0+TWO*G/seq0
ueqp=uequ-pres*(uequ-ueq0)
if (ueq > ueq0 .and. ueq<ueqp) then
  d=uequ*(ueq-ueq0)/(ueq*(uequ-ueq0))
  d=d*al
	dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
else if (ueq>=ueqp) then
  d=ONE-pres*ueq0/ueq
  d=d*al
  dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
end if

d=min(d,dmax)
d=max(d,dold)
dv=min(dv,dmax)
dv=max(dv,dvold)

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE FOR DAMAGE EVOLUTION___________________
!     __________________________________________________________________

subroutine calc_damage_exp(ueq,ueq0,uequ,eta,d,dold,dv,dvold,dmax,dtime,pres,ueqp,seq0,G,al)

IMPLICIT NONE

!     ---------------------REAL MATRIX VARIABLES------------------------
real(kind=8), intent(in)  :: ueq,ueq0,eta,dmax,dvold,dold,dtime,pres,seq0,G,al
real(kind=8), intent(inout) :: d,dv,ueqp,uequ

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

if (pres==ZERO) then
  if (ueq > ueq0) then
    d=ONE-ueq0/ueq*exp(-seq0/G*(ueq-ueq0))
    d=d*al
    dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
  end if
else
  ueqp=-G/seq0*log(pres)+ueq0
  if (ueq > ueq0 .and. ueq<ueqp) then
    d=ONE-ueq0/ueq*exp(-seq0/G*(ueq-ueq0))
    d=d*al
    dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
  else if (ueq>=ueqp) then
    d=ONE-pres*ueq0/ueq
    d=d*al
    dv=eta/(eta+dtime)*dvold+dtime/(eta+dtime)*d
  end if
end if

d=min(d,dmax)
d=max(d,dold)
dv=min(dv,dmax)
dv=max(dv,dvold)

return
end subroutine