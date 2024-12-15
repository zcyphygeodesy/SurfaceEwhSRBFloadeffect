      subroutine EwhSRBFLoadeffectpnt(calcpntfl,ewhgridfl,loadlovfl,para)
      implicit none
	character*800::calcpntfl,ewhgridfl,loadlovfl,str,astr
	character*80000::line
	integer ki,kj,ni,nj,i,j,k,n,m,kk,mk,sn,mn,kln,astat(5)
	integer nn,mm,nn1,mm1,krbf,krbf1,order,order1,minN,minN1,maxN,maxN1,Kt,Kt1,lvl,lvl1
      real*8 rec(8000),flv(40000,3),GRS(6),BLH(3),rln(3),r0,sta(6,4),rr,tmp,val,bf(8),para(20)
	real*8 hd(6),rhd(6),RAD,ae,fk,dpth,dr,nta,br,dlat,blat,dpth1,dr1,nta1,br1,dlat1,blat1
      real*8 sinf,tt,cosa,sina,mr,u15(15),ff(15),fw(15),sf(15),rst(15)
      real*8 GMr,gr,NFD(5),dln(2),rlnk(5),fi,fi1,la,la1,dla,dfi2,dnla2
	integer nlon,nlat,knd,nk,hgt,itern,NF,nd,edgn,NF1,nd1,edgn1,ks,kp,obsn,pp
	integer::status=0
	real*8,allocatable::ewh(:,:),obs(:,:),BPB(:,:),BPL(:),BB(:),xx(:),B15(:,:),chd(:),RBFi(:,:)
	real*8,allocatable::rlatlon(:,:),mpn(:,:),mdp(:,:),mp2(:,:),RBF(:,:),sr(:),dl(:),lon(:,:),RBFn(:,:)
	integer,allocatable::nln(:),nrd(:,:),node(:),enode(:),gpnt(:,:)!格网、未知数序号,每个观测量有效节点序号
	real*8,allocatable::BPB1(:,:),BPL1(:),BB1(:),x1(:),xy(:,:),B151(:,:)
	real*8,allocatable::rlatlon1(:,:),mpn1(:,:),mdp1(:,:),mp21(:,:),RBF1(:,:),sr1(:),dl1(:),lon1(:,:)
	integer,allocatable::nln1(:),nrd1(:,:),node1(:),enode1(:),gpnt1(:,:)!格网、未知数序号,每个观测量有效节点序号
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0;ae=GRS(2)
    	RAD=datan(1.d0)/45.d0;BLH(3)=0.d0;mr=36.d2/RAD;fk=3.d0/5.517d0*GRS(1)
      !单位变换
      u15(1)=1.d2;u15(2)=1.d3;u15(3)=1.d8;u15(4)=1.d8;u15(5)=mr;u15(6)=mr
      u15(7)=mr;u15(8)=mr;u15(9)=1.d3;u15(10)=1.d3;u15(11)=1.d3;u15(12)=1.d3
      u15(13)=1.d12;u15(14:15)=1.d9!mE,E
	!读负荷勒夫数
      flv=0.d0
      open(unit=8,file=loadlovfl,status="old",iostat=status)
      if(status/=0) goto 902
      do i=1,6
        read(8,'(a)') line
      enddo
      n=0
	do while(.not.eof(8))
        n=n+1
	  read(8,*,end=903)i,(flv(n,j),j=1,3)
      enddo
903   close(8)
      para(1)=1!itern-累积逼近次数
      itern=nint(para(1));knd=nint(para(2))!itern-累积逼近次数；knd-法方程解算方法
      krbf=nint(para(3));order=nint(para(4));minN=nint(para(5));maxN=nint(para(6))
      dpth=para(7)*1.d3;dr=para(8)*1.d3/ae/RAD!dr球面角距°
      krbf1=nint(para(10));order1=nint(para(11));minN1=nint(para(12));maxN1=nint(para(13))
      dpth1=para(14)*1.d3;dr1=para(15)*1.d3/ae/RAD!dr球面角距°
      lvl=nint(para(9));lvl1=nint(para(16))!Reuter等级K
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!test
      !打开地面等效水高格网文件
      open(unit=8,file=ewhgridfl,status="old",iostat=status)
      if(status/=0) goto 902
      read(8,'(a)') line
      call PickReclong(line,kln,rec,sn)
      if(sn<6)then
         close(8);goto 902
      endif
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/dble(nlon)
	hd(6)=(hd(4)-hd(3))/dble(nlat)
 	allocate(ewh(nlat,nlon), stat=astat(1))
 	allocate(obs(nlat*nlon*2,5), stat=astat(2))
 	allocate(chd(nlat*nlon*2), stat=astat(3))
	if (sum(astat(1:3)) /= 0) goto 902
      k=0;BLH(3)=0.d0
 	do i=1,nlat
	   read(8,*,end=905)(ewh(i,j),j=1,nlon)
         BLH(1)=hd(3)+(real(i)-0.5d0)*hd(6)
         do j=1,nlon
            BLH(2)=hd(1)+(real(j)-0.5d0)*hd(5)
            if(ewh(i,j)<9.d3)then
              k=k+1;call BLH_RLAT(GRS,BLH,rln)
              obs(k,1:3)=rln(1:3);obs(k,4)=ewh(i,j)*1.d-2;obs(k,5)=obs(k,4)
            endif
         enddo
      enddo
905   close(8)
      obsn=k
   !计算格网节点/未知数个数Kt
      !nln(nn)-平行圈方向格网数nn=maxi-mini+1
      !sr(nn)平行圈方向格网面积与赤道格网面积之差的百分比
      !dl(nn)平行圈方向经度间隔°Kt格网总点数-节点数、未知数
      !lon(nn,mm)格网中心经度,mm为平行圈方向最多格网数
      rhd(1:4)=hd(1:4);BLH(2)=(hd(1)+hd(2))/2.d0;BLH(3)=0.d0!!!!!!!目标格网范围用球坐标表示
      BLH(1)=hd(3);call BLH_RLAT(GRS,BLH,rln);rhd(3)=rln(2);r0=rln(1)
      BLH(1)=hd(4);call BLH_RLAT(GRS,BLH,rln);rhd(4)=rln(2);r0=(r0+rln(1))/2.d0!!!!!!!!!
      br=r0-dpth;nta=br/r0!补偿深度dpth和宽度参数nta
      dlat=180.d0/dble(lvl);nd=nint(dr/dlat+0.5d0)!dlat格网间隔,积分半径对应的格网数
      nn=nint((rhd(4)-rhd(3))/dlat+0.5);mm=nint((rhd(2)-rhd(1))/dlat+0.5)!mm平行圈方向最大格网数
      allocate(nln(nn),sr(nn),dl(nn),nrd(nn,mm),gpnt(nn,mm),lon(nn,mm),rlatlon(2*(nn+mm),2),enode(2*(nn+mm)))
      call ReuterGrid(rhd,lvl,Kt,blat,nn,mm,nln,sr,dl,nrd,lon)!Kt节点数/未知数个数
      gpnt=0!计算格网中测点数，修正Reuter格网节点数Kt,序号nrd
      call Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
      NF=nint(dr*3600)!影响半径等分,间隔1″,NF+1→[0,dr]
      allocate(RBF(NF+1,11),BB(kt),BPL(kt),node(Kt),xx(kt))
 	allocate(mpn(maxN-minN+1,NF+1), stat=astat(1))
 	allocate(mdp(maxN-minN+1,NF+1), stat=astat(2))!!!!!!!放到最后计算开空间
 	allocate(mp2(maxN-minN+1,NF+1), stat=astat(3))!!!!!!!!!!
 	allocate(BPB(Kt,Kt), stat=astat(4))
 	allocate(B15(Kt,15), stat=astat(5))
	if (sum(astat(1:5)) /= 0) goto 601
   !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
      call LegPn02(mpn,mdp,mp2,minN,maxN,NF,dr)
   !由初始补偿深度dpth,计算SRFB曲线
      call SRBF11all(RBF,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
      call Stat1d(obs(1:obsn,4),obsn,sta(1,1:4))!trn-1迭代次数
      write(*, 104) (sta(1,j),j=1,4)
104   format('    The source EWH (cm):   meam',f10.3,'     STD',f10.3,'     Min',f10.3,'     Max',f10.3)
      write(*, *)"     Please wait......"
   !组建观测方程、法方程并求解!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      BPB=0.d0;BPL=0.d0;rlnk(1)=br
      do k=1,obsn!-等效水高观测点数
        rln(1:3)=obs(k,1:3);val=obs(k,4);rr=rln(1)
        call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
        mk=0;BB=0.d0;node=0!BB-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
        do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
          if(i<1.or.i>nn)goto 1001!lon(i,j),blat第一行平行圈格网中心地心纬度°
          rlnk(2)=blat+(i-1.d0)*dlat
          do j=kj-nd,kj+nd
            if(j<1.or.j>nln(i)) goto 1002 
            if(nrd(i,j)<1) goto 1002 !nrd(i,j)未知数序号
            rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
            if(dln(2)>dr)goto 1002
            call RBFvalue(RBF(:,1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
            mk=mk+1;node(mk)=nrd(i,j)
            BB(nrd(i,j))=rr*tmp*dexp((minN-1.0)*dlog(r0/rr))
1002        continue
          enddo
1001      continue    
        enddo
        do i=1,mk
          ki=node(i)
          BPL(ki)=BPL(ki)+BB(ki)*val
          do j=1,i
            kj=node(j)
            BPB(ki,kj)=BPB(ki,kj)+BB(ki)*BB(kj)
	    enddo
        enddo
      enddo
      tmp=0.d0
	do i=1,Kt
	   do j=1,i-1
	      BPB(j,i)=BPB(i,j)
         enddo
         tmp=tmp+BPB(i,i)**2/dble(Kt)
      enddo
      tmp=dsqrt(tmp)
      !以Reuter格网四周节点未知数为零组成观测方程，抑制边缘效应。
      !节点序号数组enode
      do i=1,edgn!edgn-Reuter格网四周节点数
         ki=enode(i); BPB(ki,ki)=BPB(ki,ki)+tmp/dsqrt(dble(obsn))
      enddo
	do i=1,Kt
         BPB(i,i)=BPB(i,i)+tmp*2.d-4
      enddo
      xx=0.d0!knd=0LU分解,1Cholesky分解,2最小二乘QR分解,3,4最小范数奇异值分解
      call Equsolve(BPB,xx,Kt,BPL,knd,bf)
      kk=0!计算残差等效水高
      do k=1,obsn
        rln(1:3)=obs(k,1:3);rr=rln(1)
        call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
        mk=0;BB=0.d0;node=0!BB-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
        do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
          if(i<1.or.i>nn)goto 2001!lon(i,j),blat第一行平行圈格网中心地心纬度°
          rlnk(2)=blat+(i-1.d0)*dlat
          do j=kj-nd,kj+nd
            if(j<1.or.j>nln(i)) goto 2002 
            if(nrd(i,j)<1) goto 2002 !nrd(i,j)未知数序号
            rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
            if(dln(2)>dr)goto 2002
            call RBFvalue(RBF(:,1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
            mk=mk+1;node(mk)=nrd(i,j)
            BB(nrd(i,j))=rr*tmp*dexp((minN-1.0)*dlog(r0/rr))
2002        continue
          enddo
2001      continue    
        enddo
        val=0.d0
        do i=1,mk
          pp=node(i)
          val=val+BB(pp)*xx(pp)
        enddo
        obs(k,4)=obs(k,5)-val
        if(ki>=nd/2.and.ki<=nn-nd/2)then
          if(kj>=nd/2.and.kj<=nln(ki)-nd/2)then
            kk=kk+1;chd(kk)=obs(k,4)
          endif
        endif
      enddo
      call Stat1d(chd(1:kk),kk,sta(2,1:4))!;obs(k,5)=obs(k,4)
      write(*, 105) (sta(2,j),j=1,4)
105   format('  1st residual EWH (cm):   meam',f10.3,'     STD',f10.3,'     Min',f10.3,'     Max',f10.3)
      if(itern==0)goto 1111
!第二个SRBF用于累积逼近!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      br1=r0-dpth1;nta1=br1/r0!补偿深度dpth和宽度参数nta
      dlat1=180.d0/dble(lvl1);nd1=nint(dr1/dlat1+0.5d0)!dlat格网间隔,积分半径对应的格网数
      nn1=nint((rhd(4)-rhd(3))/dlat1+0.5);mm1=nint((rhd(2)-rhd(1))/dlat1+0.5)!mm平行圈方向最大格网数
      allocate(nln1(nn1),sr1(nn1),dl1(nn1),nrd1(nn1,mm1),gpnt1(nn1,mm1),lon1(nn1,mm1),rlatlon1(2*(nn1+mm1),2),enode1(2*(nn1+mm1)))
      call ReuterGrid(rhd,lvl1,Kt1,blat1,nn1,mm1,nln1,sr1,dl1,nrd1,lon1)!Kt节点数/未知数个数
      gpnt1=0!计算格网中测点数，修正Reuter格网节点数Kt1,序号nrd1
      call Edgnode(enode1,rlatlon1,lvl1,edgn1,lon1,blat1,nln1,gpnt1,nn1,mm1)
      NF1=nint(dr1*3600)!影响半径等分,间隔1″,NF+1→[0,dr]
      allocate(RBF1(NF1+1,11),BB1(Kt1),BPL1(Kt1),node1(Kt1),x1(Kt1),xy(Kt1,itern))
 	allocate(mpn1(maxN1-minN1+1,NF1+1), stat=astat(1))
 	allocate(mdp1(maxN1-minN1+1,NF1+1), stat=astat(2))!!!!!!!放到最后计算开空间
 	allocate(mp21(maxN1-minN1+1,NF1+1), stat=astat(3))!!!!!!!!!!
 	allocate(BPB1(Kt1,Kt1), stat=astat(4))
 	allocate(B151(Kt1,15), stat=astat(5))
	if (sum(astat(1:5)) /= 0) then
         deallocate(nln1,sr1,dl1,nrd1,gpnt1,lon1,rlatlon1,enode1,RBF1,BB1,BPL1,node1,x1,xy)
         goto 1111
      endif
   !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
      call LegPn02(mpn1,mdp1,mp21,minN1,maxN1,NF1,dr1)
   !由初始补偿深度dpth,计算SRFB曲线
      call SRBF11all(RBF1,flv,order1,krbf1,mpn1,mdp1,mp21,minN1,maxN1,NF1,nta1)
      do kp=1,itern
        obs(1:obsn,5)=obs(1:obsn,4)
        BPB1=0.d0;BPL1=0.d0;rlnk(1)=br1
        do k=1,obsn!-监测量
          rln(1:3)=obs(k,1:3);val=obs(k,4);rr=rln(1)
          call RtGridij(rln,ki,kj,blat1,lvl1,nn1,mm1,nln1,dl1,lon1)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB1=0.d0;node1=0!BB-RBF系数,node1(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd1,ki+nd1!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
            if(i<1.or.i>nn1)goto 5001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat1+(i-1.d0)*dlat1
            do j=kj-nd1,kj+nd1
              if(j<1.or.j>nln1(i)) goto 5002 
              if(nrd1(i,j)<1) goto 5002 !nrd(i,j)未知数序号
              rlnk(3)=lon1(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr1)goto 5002
              call RBFvalue(RBF1(:,1),NF1,dr1,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              mk=mk+1;node1(mk)=nrd1(i,j)
              BB1(nrd1(i,j))=rr*tmp*dexp((minN1-1.0)*dlog(r0/rr))
5002          continue
            enddo
5001        continue    
          enddo
          do i=1,mk
            ki=node1(i)
            BPL1(ki)=BPL1(ki)+BB1(ki)*val
            do j=1,i
              kj=node1(j)
              BPB1(ki,kj)=BPB1(ki,kj)+BB1(ki)*BB1(kj)
	      enddo
          enddo
        enddo
        tmp=0.d0
  	  do i=1,Kt1
	     do j=1,i-1
	        BPB1(j,i)=BPB1(i,j)
           enddo
           tmp=tmp+BPB1(i,i)**2/dble(Kt1)
        enddo
        tmp=dsqrt(tmp)
        !以Reuter格网四周节点未知数为零组成观测方程，抑制边缘效应。
        !节点序号数组enode
        do i=1,edgn1!edgn-Reuter格网四周节点数
           ki=enode1(i); BPB1(ki,ki)=BPB1(ki,ki)+tmp/dsqrt(dble(obsn))
        enddo
	  do i=1,Kt
           BPB1(i,i)=BPB1(i,i)+tmp*2.d-4
        enddo
        x1=0.d0!knd=0LU分解,1Cholesky分解,2最小二乘QR分解,3,4最小范数奇异值分解
        call Equsolve(BPB1,x1,Kt1,BPL1,knd,bf)!4、5测试通过
        xy(1:Kt1,kp)=x1(1:Kt1);kk=0!计算残差等效水高
        do k=1,obsn
          rln(1:3)=obs(k,1:3);rr=rln(1)
          call RtGridij(rln,ki,kj,blat1,lvl1,nn1,mm1,nln1,dl1,lon1)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB1=0.d0;node1=0!BB-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd1,ki+nd1!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
            if(i<1.or.i>nn1)goto 6001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat1+(i-1.d0)*dlat1
            do j=kj-nd1,kj+nd1
              if(j<1.or.j>nln(i)) goto 6002 
              if(nrd1(i,j)<1) goto 6002 !nrd(i,j)未知数序号
              rlnk(3)=lon1(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr1)goto 6002
              call RBFvalue(RBF1(:,1),NF1,dr1,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
              mk=mk+1;node1(mk)=nrd1(i,j)
              BB1(nrd1(i,j))=rr*tmp*dexp((minN1-1.0)*dlog(r0/rr))
6002          continue
            enddo
6001        continue    
          enddo
          val=0.d0
          do i=1,mk
            pp=node1(i)
            val=val+BB1(pp)*x1(pp)
          enddo
          obs(k,4)=obs(k,5)-val
          if(ki>=nd1/2.and.ki<=nn-nd1/2)then
            if(kj>=nd1/2.and.kj<=nln1(ki)-nd1/2)then
              kk=kk+1;chd(kk)=obs(k,4)
            endif
          endif
        enddo
        call Stat1d(chd(1:kk),kk,sta(kp+2,1:4))!;obs(k,5)=obs(k,4)
        write(*, 106) kp+1,(sta(kp+2,j),j=1,4)
106     format('  ',i1,'st residual EWH (cm):   meam',f10.3,'     STD',f10.3,'     Min',f10.3,'     Max',f10.3)
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1111  kk=0;open(unit=8,file=calcpntfl,status="old",iostat=status)
      if(status/=0)goto 909
      open(unit=10,file="reslt.txt",status="replace")
      read(8,'(a)') line
      write(10,101)trim(line)
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickReclong(line,kln,rec,sn)
         if(sn<4)goto 906
         BLH(1)=rec(3);BLH(2)=rec(2);BLH(3)=rec(4)
         call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
         call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
         ff(1)=rr;ff(2)=fk/gr/rr;ff(3)=fk/rr/rr;ff(4)=ff(3);ff(5)=ff(2)/rr;ff(6)=ff(5)
         ff(7)=-ff(5);ff(8)=ff(2);ff(9)=ff(2);ff(10)=ff(3)/rr;ff(11)=ff(5)/rr
         rst(1:15)=0.d0
         call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
         mk=0;B15=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
         do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
           if(i<1.or.i>nn)goto 3001!lon(i,j),blat第一行平行圈格网中心地心纬度°
           rlnk(2)=blat+(i-1.d0)*dlat
           do j=kj-nd,kj+nd
             if(j<1.or.j>nln(i)) goto 3002 
             if(nrd(i,j)<1) goto 3002 !nrd(i,j)未知数序号
             rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
             if(dln(2)>dr)goto 3002
             mk=mk+1;node(mk)=nrd(i,j);sf=0.d0
             do k=1,11!!!!!!!!!!!!!!!!!!!!!!!!!11
               call RBFvalue(RBF(:,k),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
               sf(k)=ff(k)*tmp*dexp((minN-1.0)*dlog(r0/rr))
             enddo
             tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
             fi=rln(2)*RAD;fi1=rlnk(2)*RAD;la=rln(3)*RAD;la1=rln(3)*RAD;dla=(rlnk(3)-rln(3))*RAD
	       cosa=(dcos(fi)*dsin(fi1)-dsin(fi)*dcos(fi1)*dcos(dla))/sinf
	       sina=dcos(fi1)*dsin(dla)/sinf
             dfi2=-dsin(fi)*dsin(fi1)-dcos(fi)*dcos(fi1)*dcos(dla)+tt*cosa**2;dfi2=dfi2/sinf
             dnla2=-dcos(fi1)*dsin(dla)+tt*dcos(fi)*sinf**2;dnla2=dnla2/sinf
             fw(1:5)=sf(1:5);fw(6)=fw(5);fw(7)=sf(6);fw(8)=sf(6);fw(9)=sf(7);fw(10)=sf(7)
             fw(11)=sf(8);fw(12)=sf(9);fw(13)=sf(10);fw(14)=sf(11);fw(15)=-sf(11)/((dcos(fi))**2)
             fw(5)=fw(5)*cosa;fw(6)=fw(6)*sina;fw(7)=fw(7)*cosa;fw(8)=fw(8)*sina
             fw(9)=fw(9)*cosa;fw(10)=fw(10)*sina;fw(14)=fw(14)*dfi2;fw(15)=fw(15)*dnla2
             B15(nrd(i,j),1:15)=fw(1:15)
3002         continue
           enddo
3001       continue    
         enddo
         ff=0.d0
         do i=1,mk
           ki=node(i)
           do k=1,15
             ff(k)=ff(k)+B15(ki,k)*xx(ki)
           enddo
         enddo
         rst(1:15)=rst(1:15)+ff(1:15)*u15(1:15)
         do kp=1,itern
            call RtGridij(rln,ki,kj,blat1,lvl1,nn1,mm1,nln1,dl1,lon1)!计算rln在格网中的序号ki>0,kj>0
            x1(1:Kt1)=xy(1:Kt1,kp);mk=0;B151=0.d0;node1=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
            do i=ki-nd1,ki+nd1!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
              if(i<1.or.i>nn1)goto 7001!lon(i,j),blat第一行平行圈格网中心地心纬度°
              rlnk(2)=blat1+(i-1.d0)*dlat1
              do j=kj-nd1,kj+nd1
                if(j<1.or.j>nln1(i)) goto 7002 
                if(nrd1(i,j)<1) goto 7002 !nrd(i,j)未知数序号
                rlnk(3)=lon1(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
                if(dln(2)>dr1)goto 7002
                mk=mk+1;node1(mk)=nrd1(i,j);sf=0.d0
                do k=1,11!!!!!!!!!!!!!!!!!!!!!!!!!11
                  call RBFvalue(RBF1(:,k),NF1,dr1,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
                  sf(k)=ff(k)*tmp*dexp((minN-1.0)*dlog(r0/rr))
                enddo
                tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
                fi=rln(2)*RAD;fi1=rlnk(2)*RAD;la=rln(3)*RAD;la1=rln(3)*RAD;dla=(rlnk(3)-rln(3))*RAD
	          cosa=(dcos(fi)*dsin(fi1)-dsin(fi)*dcos(fi1)*dcos(dla))/sinf
	          sina=dcos(fi1)*dsin(dla)/sinf
                dfi2=-dsin(fi)*dsin(fi1)-dcos(fi)*dcos(fi1)*dcos(dla)+tt*cosa**2;dfi2=dfi2/sinf
                dnla2=-dcos(fi1)*dsin(dla)+tt*dcos(fi)*sinf**2;dnla2=dnla2/sinf
                fw(1:5)=sf(1:5);fw(6)=fw(5);fw(7)=sf(6);fw(8)=sf(6);fw(9)=sf(7);fw(10)=sf(7)
                fw(11)=sf(8);fw(12)=sf(9);fw(13)=sf(10);fw(14)=sf(11);fw(15)=-sf(11)/((dcos(fi))**2)
                fw(5)=fw(5)*cosa;fw(6)=fw(6)*sina;fw(7)=fw(7)*cosa;fw(8)=fw(8)*sina
                fw(9)=fw(9)*cosa;fw(10)=fw(10)*sina;fw(14)=fw(14)*dfi2;fw(15)=fw(15)*dnla2
                B151(nrd1(i,j),1:15)=fw(1:15)
7002            continue
              enddo
7001          continue    
            enddo
            ff=0.d0
            do i=1,mk
              ki=node1(i)
              do k=1,15
                ff(k)=ff(k)+B151(ki,k)*x1(ki)
              enddo
            enddo
            rst(1:15)=rst(1:15)+ff(1:15)*u15(1:15)
         enddo
         kk=kk+1
         !输出等效水高,高程异常,地面重力,扰动重力,地倾斜,垂线偏差,水平,径向,正常高,扰动梯度,水平梯度
         write(10,'(a,13F12.4,2ES13.4)')trim(line),(rst(i),i=1,15)
906      continue
	enddo
 	close(8)
      close(10)
      allocate(RBFi(NF+1,11),RBFn(maxN-minN+1,11))
      call SRBFone(RBFi,RBFn,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
      open(unit=10,file="SRBFspc.txt",status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     * " ewh,hgtanomaly,terrgrav,gravdisturb,tilt,vdeflect,hdisplace,radial,ortheight,gradient,vertgrad" 
      tmp=dr/dble(NF)*RAD*ae*1.d-3
 	do i=1,NF
	   write(10,'(F12.3,11F13.5)')-(NF-i+1.0)*tmp,(RBFi(NF-i+2,j),j=1,11)
      enddo
 	do i=1,NF+1
	   write(10,'(F12.3,11F13.5)')(i-1.0)*tmp,(RBFi(i,j),j=1,11)
      enddo
      close(10)
      open(unit=10,file="SRBFdgr.txt",status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     * " ewh,hgtanomaly,terrgrav,gravdisturb,tilt,vdeflect,hdisplace,radial,ortheight,gradient,vertgrad" 
 	do i=minN,maxN
	   write(10,'(I8,11F13.5)')i,(RBFn(i-minN+1,j),j=1,11)
      enddo
      close(10)
      deallocate(RBFi,RBFn)
909   if(itern>0.1)deallocate(mpn1,mdp1,mp21,BPB1,B151)
      if(itern>0.1)deallocate(nln1,sr1,dl1,nrd1,gpnt1,lon1,rlatlon1,enode1,RBF1,BB1,BPL1,node1,x1,xy)
      deallocate(mpn,mdp,mp2,BPB,B15)
601   deallocate(nln,sr,dl,nrd,gpnt,lon,rlatlon,enode,RBF,BB,BPL,node,xx)
904   deallocate(ewh,obs,chd)
902	continue
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      write (*,*)'  The program also outputs the SRBF spatial curve file SRBFspc.rbf and spectral curve file
     & SRBFdgr.rbf of 11 kinds of geodetic variations into the current directory.'
      pause
101   format(a,40F12.4)
103   format(4F12.6,2F12.8)
      end
