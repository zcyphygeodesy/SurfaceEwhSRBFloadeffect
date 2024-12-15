!  SurfaceEwhSRBFloadeffect.f90 
!
!  FUNCTIONS:
!  SurfaceEwhSRBFloadeffect - Entry point of console application.
!
      program SurfaceEwhSRBFloadeffect
      implicit none
	character*800::calcpntfl,ewhgridfl,loadlovfl
      real*8 para(20)
!---------------------------------------------------------------------
	write(calcpntfl,*) "calcpnt.txt" !The space calculation point file
	write(ewhgridfl,*) "soilewh20180131.dat"!The residual EWH variation grid file
	write(loadlovfl,*) "Love_load_cm.dat"!The load love number file
      para(1)=1!累积逼近次数cumulative approach times
      !选择法方程解算方法 Select the method of the solution of normal equation
      para(2)=1!1-LU分解,2-Cholesky分解,3-最小二乘QR分解,4,5-最小范数奇异值分解
      !输入主SRBF参数 input the main SRBF parameters
      !0-the radial multipole kernel function, 1-Poisson wavelet kernel function
      para(3)=0!0-径向多级子核函数，1-Poisson小波核函数
      para(4)=0!多级次数 the order number m
      para(5)=15;para(6)=900!SRBF最小最大阶数 Minimum and maximum degree
      para(7)=30.d0!Bjerhammar球面相对地面的平均深度 the Bjerhammar sphere burial depth (km)
      para(8)=150.d0!SRBF中心作用距离(km) the action distance of SRBF center
      para(9)=1800!主SRBF的Reuter等级K the Reuter network level K for main SRBF
      !输入累积逼近SRBF参数（与主SRBF参数相同）
      !input the SRBF parameters for cumulative approach
      para(10)=1; para(11)=0; para(12)=45; para(13)=1800
      para(14)=20.d0; para(15)=90.d0; para(16)=1800
      call EwhSRBFLoadeffectpnt(calcpntfl,ewhgridfl,loadlovfl,para)
      end
!
!******************************************************************
!
      subroutine Equsolve(BB,xx,nn,BL,knd,bf)
!调用mkl_lapack95_ilp64.lib解大型方程组BB.xx=BL
!knd=1 LU分解,2 Cholesky分解,3 最小二乘QR分解,4最小范数奇异值分解
!knd=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square qr decomposition,
!   =4 Minimum norm singular value decomposition)
!bf(8)-解的性质
!---------------------------------------------------------------
      USE OMP_LIB
      implicit none
	integer::nn,knd,nm,inf,lwk,rk,astat(20)
	real*8::BB(nn,nn),BL(nn),xx(nn),bf(8),rnd
      real*8,allocatable::s(:),wk(:)
	integer,allocatable::iwk(:),ipv(:)
	integer::i,ki
!-----------------------------------------------------------------------
      nm=nn;bf=0.d0;rk=0;xx=BL
 	allocate(s(nn), stat=astat(1))
 	allocate(wk(nn*nn), stat=astat(2))
 	allocate(iwk(nn*nn), stat=astat(3))
 	allocate(ipv(nn), stat=astat(4))
	if (sum(astat(1:4)) /= 0) then
        bf(1)=1.d0;return
      endif
      lwk=-1;s=0.d0
      if(knd==3)call dgels('No transpose',nm,nm,1,BB,nm,xx,nm,wk,lwk,inf)
      if(knd==4)call dgelsd(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,iwk,inf)
      if(knd==5)call dgelss(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,inf)
      nm=nn;lwk = nm**nm
      if(knd==1)call dgesv(nm,1,BB,nm,ipv,xx,nm,inf)
      if(knd==2)call dposv('Upper',nm,1,BB,nm,xx,nm,inf)
      if(knd==3)call dgels('No transpose',nm,nm,1,BB,nm,xx,nm,wk,lwk,inf)
      if(knd==4)call dgelsd(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,iwk,inf)
      if(knd==5)call dgelss(nm,nm,1,BB,nm,xx,nm,s,-1.d0,rk,wk,lwk,inf)
      if( inf >0 ) then
        bf(2)=1.d0;goto 2001!计算失败-,请调整参数重新计算
      endif
      bf(3)=rk
2001	deallocate(s,wk,iwk,ipv)
	return
      end
