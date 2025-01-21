! **********************************CCD Image preprocessing*************************************
!program name:03match.f90
! function:1）读取gaia星表和卫星历表
!          2) 根据星表及卫星历表位置及图像检测形象的信息，程序自动确认底片旋转角度及视场中目标星坐标
!          3) 自动完成匹配并归算目标星位置
!          4) 输出目标星的时间、位置、O-C、曝光时间、图像名等信息
!    input:0）由detect程序得到的*.reg文件
!          1) 配置文件信息（adias.cfg）
!            1.1) 望远镜标签：3tele_label=ss156
!            1.2) 望远镜焦距：3tele_focal=15600.0
!            1.3) CCD比例尺：3ccd_scale=0.048
!            1.4) 底片模型参数：3modeltype=6
!            1.5) 底片旋转角度初始值：3plate_angle=0    
!            1.6) 目标星历表文件：3obj_ephfile=F:\Adias\obs\EPH_N1_200608_DE431.DAT
!            1.7) gaia星表：3gaia_catfile=F:\Adias\obs\CAT_N1_200608.DAT
!            1.8) gaia星表最小星等：3gaia_minmag=3
!            1.9) gaia星表最小星等：3gaia_maxmag=16
!            1.10) 参考星匹配限：3match_limit=5.0
!   output:1）参考星信息文件：*.ref.reg
!          2) 目标星的详细信息（时间、位置、O-C，图像sig0,曝光时间，图像名称）：object.out
!  子程序：1）find_obj_baseangle(fitsfile,field_angle,pscale,fl,limit_match,eph_T,eph_ra,eph_de,n_eph,&
!         & gaia_ra,gaia_de,gaia_pr,gaia_pd,n_gaia,modeltype,&
!         & year,month,obj_T,hh,mm,ss,exptime,obj_x0,obj_y0,obj_obsra,obj_obsde,obj_ephra,obj_ephde,&
!        &n_match1,sig1,ref_x0,ref_y0,ref_ra0,ref_de0,ref_mag0,par0,nostar)fits_header(fitsfile,bitpix,naxis1,naxis2,bscale,bzero,gain,exptime,year,month,day,hh,mm,ss)
!                      主要子程序，基于固定底片旋转角度+星表+目标历表，自动寻找目标，并归算其观测位置，并输出参考星信息
!          2) fits_header(fitsfile,naxis1,naxis2,bscale,bzero,gain,exptime,year,month,day,hh,mm,ss)
!                       读取fits头文件信息
!          3) CAL_RL(RA01,DEC01,RA02,DEC02,RL)
!                      计算两颗星大圆弧上距离（用于检验是否匹配）
!          4) xy2rade(par,nm,x,y,ra0,de0,ra,de)
!                      由量度坐标转换为赤道坐标
!          5) sol_par(x,y,ra,de,n,ra0,de0,par,nm,sig0,nused,vksi,vyit,p,IFLAG1)   
!                      由参考星的量度坐标、星表位置等信息归算底片常数par
!          6) rade2ky(ra,de,ra0,de0,ksi,yit)
!                      由赤道坐标转为理想坐标
!          7) LEAST(A,Y,X,N1,K,NE,ERROR,SITA2,nrefused,rmag,p,vksi,vyit,IFLAG1) 
!                      最小二乘法解线性方程中待解参数
!          8) MUL(A1,A2,A12,K1,K2,K3,K4) 
!                      两矩阵相乘
!          9) MATINV (N,NUP,A)
!                      矩阵求逆
!          10) am2hms(ra,de,hh,mm,ss,sign,ddd,dmm,dss) 
!                      将赤经赤纬度，变为时分秒，度分秒
!          11) sla_CALDJ (IY, IM, ID, DJM, J)
!                      通过年月日计算简约儒略日，J标记年月日格式是否正确
!          12) sla_CLDJ (IY, IM, ID, DJM, J)
!                      两矩阵相乘
!          13) sort_T(t,line,n) ENLGR(X,Y,N,T,Z)
!                      按照行中某一列对全行进行排序
!          14) ENLGR(X,Y,N,T,Z) 未用
!                      对N个数据的时间、赤经（赤纬）数组，内插T时刻时赤经（赤纬）位置
!          15) interT(x,y,n,t,z) 本程序用
!                      对N个数据的时间、赤经（赤纬）数组，内插T时刻时赤经（赤纬）位置
! by zhy 2019.09.01
! V2.2 add the enhance function by zhy 2019.09.24
! V3.0  Chinese comment by zhy 2020.01.30
!20200212 确定旋转角度部分，增加粗匹配和细匹配，粗10度步长，2*limit_match匹配限，匹配星数量>6，sig<1则开始细匹配，步长1，匹配限恢复
!20200213  
!1）根据当前幅图像目标历表位置，从星表提取15角分为半径大小视场的参考星，用于去匹配，加快速度；
!2）对于视场检测星象，使用前100颗进行匹配，认为目标出现在前80颗（加快速度）；
!3）粗匹配，3度步长，匹配限1.5*limit_match，当匹配星数>6，0<sig<1,从此角度前3度进入细匹配，共匹配6度，此时1度步长，1倍匹配限；
!4）如果第一幅匹配失败，自动读取下一幅进行匹配，最多使用3幅匹配，若3幅均失败，自动使用前一天匹配角度进行匹配。
!5）目前对于156tele从2006-2011年，共6年，10个月，65天，只有2008年9月无法匹配得出结果（原因待查）；
!6) 程序匹配+残差归算，花费时间共计：10分钟内。
!20200216：注意，观测数据请尽量在同一个月，历表一定要在同月
!20200216 增加时间修正delta_T
!20200216 归算背景时曝光时间归一化（背景和星象均已减去bias)(修改！！！)
!20200217 停下计算背景星等的测光工作，集中关注星象定心精度。
!20200217 目前find_obj程序中增加，提取目标星亮度，00ref文件中为使用参考星和归算目标星的星等.
!20200222 find程序中，第一次匹配使用1.5倍limit_match，粗匹配成功后，重新计算par,使用细匹配，1/3*limit_match
!                    设定如果参考星数>20,只使用前20颗星参与底片归算，最小二乘时，也会删除误差大的星，目的只使用亮星匹配
!         配置文件中确定角度可直接给出，如果不确定填写0.方便参考星少无法匹配时用。
!20210814处理天王星卫星数据，156观测数据，不改头文件观测时间标记的字符，目前将读取头文件的子程序中关键字改为‘DATE-STA=’（原来乔老师写程序修改头文件为了满足astrometica软件使用）
!!20220724修改，待检测，为什么当obj_obsra/de为0时也输出。600行
      program match
      implicit none
      integer*4 i,j,k,m,n,ns,nr,nref,ic,ic1,kk,i11,i2,j1,j2,iii,nostar
      integer*4 n_eph,n_gaia,n_det,n_match,n_objpos,jjj,n_obj_det
      integer*4 year,month,day,hh,mm,modeltype,refused1,n_match1,obj_j,n_match0,jss
      integer*4 naxis1,naxis2,bitpix,iflag2,i1,s_loop,n_lose,s_loop0,n_stop
      parameter(m=900000)
      character*300 tempcha,objfitsfile,regfile,refregfile,fitsfilepath,ephsource
      character*300 objccdviewfile,ephfile,gaiacatfile,objoutfile,line(m),ephpath(10),gaiacatpath(10)
      character*100 cc,fitspath,newoutfile0,refoutfile,tele_label,month_a
      integer*4 rahj,ramj,decdj,decmj,obj,obj_total,fsize
      real*8 rasj,decsj
      character*1 pmj,obj_label
      real*8 eph_T(m),eph_ra(m),eph_de(m),min_mag,max_mag,limit_match
      real*8 gaia_ra(m),gaia_de(m),gaia_pr(m),gaia_pd(m),gaia_mag(m),cat_ra,cat_de,cat_pr,cat_pd,cat_mag
      real*8 gaia_ra0(m),gaia_de0(m),gaia_mag0(m),fitsdata(m)
      real*8 obj_x,obj_y,obj_ephra,obj_ephde,obj_T,obj_obsra,obj_obsde,obj_resra,obj_resde
      real*8 par(30),par1(30),bscale,bzero,gain,ss,exptime,flux,ff,fff
      real*8 det_x(m),det_y(m),det_ra(m),det_de(m),det_flux(m),x0,y0,ra0,de0
      real*8 ref_x(m),ref_y(m),ref_ra(m),ref_de(m),ref_mag(m),det_xn(m),det_yn(m)
      real*8 calpm_epoch,rl,pi,q,fl,pscale,djm,field_angle,angle_step,field_angle0,sig00,field_angle_bak
      real*8 sig0,VKSI1(m),VYIT1(m),weight1(m),p,T(m),obj_x0,obj_y0,limit_match0,delta_t,obj_flux,snr0
      real*8 ref_flux(m),bkgd,bkgd_flux,bkgd_mag,std_mag,bkgd0_mag(m),bkgdsigma,arcsec_pix
      integer*4 n_day,n_day_lose,n_confirm,n_gaia0,n_new,npx

!01--读配置文件       
       open(1,file='D:\00adias\exe\fitspath.in',status='old')
!       open(1,file='D:\00adias\exe\200809fitspath.in',status='old')
       n_day=0
       n_day_lose=0
77     read(1,'(a100)',end=777) fitspath
       fitspath=trim(fitspath)
       n_day=n_day+1  !用以记录共处理几天数据
    write(*,'(i3,a,a,a)')  n_day,'-','本日fits文件：',trim(fitspath)
!   按照配置文件内容逐日处理   
     open(11,file='D:\00adias\exe\adias.cfg')
!     open(11,file='D:\00adias\exe\200809adias.cfg')
      do i=1,m
        read(11,'(a100)',end=1000) tempcha

         ic=index(tempcha,'3tele_label=') 
	      if(ic.gt.0) then
  	      ic1=index(tempcha,'=')
	        read(tempcha(ic1+1:ic1+50),*) tele_label
          endif 
          
        ic=index(tempcha,'3tele_focal=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+10),*) fl
        endif

        ic=index(tempcha,'3ccd_scale=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) pscale
        endif

        ic=index(tempcha,'3ccd_fieldsize=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) fsize
        endif
        
        ic=index(tempcha,'3modeltype=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) modeltype
        endif
    
        ic=index(tempcha,'3plate_angle=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) field_angle
        endif          

        ic=index(tempcha,'3obj_ephsource=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+100),*) ephsource
        endif        

        ic=index(tempcha,'3gaia_minmag=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) min_mag
        endif

        ic=index(tempcha,'3gaia_maxmag')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) max_mag
        endif
        
        ic=index(tempcha,'3match_limit=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) limit_match
        endif
        
        ic=index(tempcha,'3delta_T=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) delta_t
        endif
        
        ic=index(tempcha,'4specified-output=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+100),*) newoutfile0
        endif 
        
        ic=index(tempcha,'3obj_total=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+30),*) obj_total !待处理目标数
      
        endif     
 !      根据待处理目标数将各文件路径放入ephpath和gaiacatpath数组中进行存储
        if(obj_total.gt.0.and.obj_total.lt.10)then
          do obj=1,obj_total
            write(obj_label,'(i1)')obj
            ic=index(tempcha,'3obj_ephfile'//trim(obj_label)//'=')
            if(ic.gt.0) then
              ic1=index(tempcha,'=')
              read(tempcha(ic1+1:ic1+100),*) ephpath(obj)
            endif
           
            ic=index(tempcha,'3gaia_catfile'//trim(obj_label)//'=')
            if(ic.gt.0) then
              ic1=index(tempcha,'=')
              read(tempcha(ic1+1:ic1+100),*) gaiacatpath(obj)
            endif          
          enddo 
          endif

      enddo   !结束读配置文件  
     
1000  close(11)
      
    !  limit_match=8 
      pi=4*datan(1.0d0)
      q=pi/180.0d0
!多目标判断
      obj=0    
9999  obj=obj+1 !用于标记目标数
      write(*,'(a,i2,a,i2)')'当前处理目标：',obj,'/',obj_total
      gaiacatfile=gaiacatpath(obj)
      ephfile=ephpath(obj)
      
      print*,'Gaia星表路径：',gaiacatpath(obj)
      print*,'目标历表路径：',ephpath(obj)

!!02--************read GAIA catalog-------------
      open(4,file=gaiacatfile,status='old')
      k=1
      do j=1,60
          read(4,*)
      enddo
      k=1
      do j=1,m
          read(4,'(a300)',end=9951) tempcha
          if(tempcha(1:4).eq.'    '.or.tempcha(1:4).eq.'#END') goto 311
!gaia exp:21 17 38.5772643376 -15 30 01.349582279 319.41072917498 -15.50039051142    -2.110    -3.630 14.8833    
          read(tempcha,'(39x,2f16.11,2f10.3,f8.4)') cat_ra,cat_de,cat_pr,cat_pd,cat_mag
          
          if(cat_mag.ge.min_mag.and.cat_mag.le.max_mag) then  
            gaia_ra(k)=cat_ra
            gaia_de(k)=cat_de
            gaia_pr(k)=cat_pr
            gaia_pd(k)=cat_pd
            gaia_mag(k)=cat_mag
            k=k+1   
          endif  
      enddo
311   n_gaia=k-1 !参考星表星数
9951  close(4)   
!!031--************read satellite ephemeris（IMCCE）-------------
    if(trim(ephsource).eq.'IMCCE') then
	   write(*,*)'---------------读取IMCCE历表---------------'
      eph_T=0
      eph_ra=0
      eph_de=0
      open(2,file=ephfile,status='old')
      do i=1,10
        read(2,*)
      enddo
      k=1
	  do i=1,m
110     read(2,'(a300)',end=111) tempcha
      !  print*,tempcha
        if(tempcha(1:3).eq.'---') goto 110
        read(tempcha,*) year,month,day,hh,mm,ss,eph_ra(k),eph_de(k)
!all the data must be in the same month
        eph_T(k)=day+(hh*3600.0+mm*60.0+ss)/86400d0
        eph_ra(k)=eph_ra(k)*15d0
        k=k+1
      enddo
111 	n_eph=k-1
      close(2) 
    endif
 !读取jpl历表结束============ 
!!032--************read satellite ephemeris(JPL)-------------  
    if(trim(ephsource).eq.'JPL') then
	write(*,*)'---------------读取JPL历表---------------'
      eph_T=0
      eph_ra=0
      eph_de=0
      open(2,file=ephfile,status='old')
      do i=1,40
        read(2,*)
      enddo
      k=1
	  do i=1,m
210     read(2,'(a300)',end=211) tempcha
      !  print*,tempcha
        if(tempcha(1:5).eq.'$$EOE') goto 211
      read(tempcha,220)year,month_a,day,hh,mm,jss,&
	  &rahj,ramj,rasj,pmj,decdj,decmj,decsj
220   FORMAT(X,I4,X,A3,X,I2.2,X,I2.2,X,I2.2,X,I2.2,5X,I2.2,X,I2.2,X,F9.6,X,A1,I2.2,X,I2.2,X,F8.5)
      IF (month_a .EQ. 'Jan') month= 1
      IF (month_a .EQ. 'Feb') month= 2
      IF (month_a .EQ. 'Mar') month= 3
      IF (month_a .EQ. 'Apr') month= 4
      IF (month_a .EQ. 'May') month= 5
      IF (month_a .EQ. 'Jun') month= 6
      IF (month_a .EQ. 'Jul') month= 7
      IF (month_a .EQ. 'Aug') month= 8
      IF (month_a .EQ. 'Sep') month= 9
      IF (month_a .EQ. 'Oct') month= 10
      IF (month_a .EQ. 'Nov') month= 11
      IF (month_a .EQ. 'Dec') month= 12
!all the data must be in the same month
!内插时按照天的小数内插，所以以上月份判断实际上没有用。
!程序后续改进，可进行不分月份内插。
		if(pmj=='-')then
			decdj=(-1)*decdj
			decmj=(-1)*decmj
			decsj=(-1)*decsj 
		end if	  

	  eph_T(k)=day+(hh*3600.0+mm*60.0+jss)/86400d0
	  eph_ra(k)=rahj*15.0d0+ramj*15.0d0/60.0d0+rasj*15.0d0/3600.0d0
	  eph_de(k)=decdj+decmj/60.0d0+decsj/3600.0d0  
   !   print*,eph_T(k),eph_ra(k),eph_de(k)

        k=k+1
      enddo
211 	n_eph=k-1
      close(2)  
      endif
 !读取jpl历表结束============   
      
      write(*,'(a,i5)')'   读取星表恒星数：',n_gaia
      write(*,'(a,i5)')'   读取历表位置数：',n_eph
      write(*,*)       '   观测使用望远镜：',trim(tele_label)
 !     pause!0815
    !  print*,field_angle
    !  pause
 !     if(field_angle.gt.0) goto 222 !若此角度设置大于0，则认为直接使用此角度
!03--使用第一幅图像确定当日底片旋转角度（请确定前三幅图像有效）
      if(field_angle>0)  goto 100
!       field_angle=345.0  !一般旋转角在355左右，设定此值，方便快速确定角度
      n_confirm=1
!031--处理第一幅图像，计算旋转角（起始角度为配置文件读取，步长为3，当匹配上时，步长变为1所有星计算匹配成功星数）          
      open(11,file=trim(fitspath)//'\'//'fits.lst') 
333   read(11,*,end=999) objfitsfile
		objfitsfile=trim(objfitsfile)
      ref_x=0.0
      ref_y=0.0
      ref_ra=0.0
      ref_de=0.0
!      field_angle=0
      angle_step=2.0 !5
      limit_match0=limit_match*1.5 
   !   print*,'匹配极限：',limit_match,limit_match!0815
      if(trim(tele_label)=='ss156')then
        call fits_header_156(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
      !  print*,bzero,gain,exptime,year,month,day,hh,mm,ss!0815
        elseif(trim(tele_label)=='ss156_2014')then
        call fits_header_156_2014(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='ss156_2016')then
        call fits_header_156_2016(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='km100')then
        call fits_header_km100(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='km100B')then
        call fits_header_km100B(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='lj240')then
        call fits_header_lj240(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        else
            write(*,*)'请检查fits头文件格式！'
        endif
        
!all the data must be in the same month  

        obj_T=day+(hh*3600+mm*60+ss)/86400d0
        
        call sla_CALDJ (year,month,day,djm,jjj)
        djm=djm+(hh*3600+mm*60+ss)/86400d0
        calpm_epoch=year+((month-1)*30+day)/365d0!use to compute the proper motion
!interpolation the epheris to get the c-data   
        call ENLGR(eph_T,eph_ra,n_eph,obj_T,obj_ephra)!单位为度
        call ENLGR(eph_T,eph_de,n_eph,obj_T,obj_ephde)
        !call interT(eph_T,eph_ra,n_eph,obj_T,obj_ephra)
        !call interT(eph_T,eph_de,n_eph,obj_T,obj_ephde)
!!************对参考星表计算自行-------------
!        print*,'本图目标时间及radec位置：',obj_T,obj_ephra,obj_ephde!0815
        n=1
        do k=1,n_gaia 
            if(gaia_ra(k).gt.(obj_ephra-fsize/60.0).and.gaia_ra(k).lt.(obj_ephra+fsize/60.0)&
                &.and.gaia_de(k).gt.(obj_ephde-fsize/60.0).and.gaia_de(k).lt.(obj_ephde+fsize/60.0))then !人为约束半度以内星表
              gaia_de0(n)=gaia_de(k)+(calpm_epoch-2016.0)*gaia_pd(k)/1000d0/3600d0  !注意自行的单位
              gaia_ra0(n)=gaia_ra(k)+(calpm_epoch-2016.0)*gaia_pr(k)/1000d0/3600d0/dcosd(gaia_de(n))
              gaia_mag0(n)=gaia_mag(k)
              n=n+1
        endif
        enddo   
        n_gaia0=n-1
        write(*,'(a,i5)')  '   用于匹配恒星数：',n_gaia0
            n_match0=0
            sig00=0.0
            s_loop=0
            s_loop0=0
            obj_x0=0.0d0
            obj_y0=0.0d0
            ref_x=0.0d0
            ref_y=0.0d0
            ref_ra=0.0d0
            ref_de=0.0d0
            n_stop=0
        write(*,'(a,a)')'   使用此图确定旋转角度:',trim(objfitsfile)    
123     if(field_angle.lt.0) field_angle=field_angle+360
  !      print*,'field_angle',field_angle!0815
  !      print*,'子函数参数：',objfitsfile,tele_label,field_angle,pscale,fl,limit_match0
        call find_obj_baseangle(objfitsfile,tele_label,field_angle,pscale,fl,limit_match0,obj_ephra,obj_ephde,&
        !call find_obj_baseangle02(objfitsfile,field_angle,pscale,fl,limit_match0,obj_ephra,obj_ephde,&
         & gaia_ra0,gaia_de0,gaia_mag0,n_gaia0,modeltype,&
         & obj_x,obj_y,obj_flux,snr0,obj_obsra,obj_obsde,n_match1,sig0,ref_x,ref_y,ref_flux,ref_ra,ref_de,ref_mag,par1,nostar)
   !  print*,'子函数参数：',objfitsfile,tele_label,field_angle,pscale,fl,limit_match0,&!obj_ephra,obj_ephde,&
        !call find_obj_baseangle02(objfitsfile,field_angle,pscale,fl,limit_match0,obj_ephra,obj_ephde,&
         !& gaia_ra0,gaia_de0,gaia_mag0,n_gaia0,modeltype,&
    !     & obj_x,obj_y,obj_flux,snr0,obj_obsra,obj_obsde,n_match1,sig0!,ref_x,ref_y,ref_flux,ref_ra,ref_de,ref_mag,par1,nostar
!       write(*,'(a,4f10.2)')'目标星信息',obj_x,obj_y,obj_flux,snr0!0815
        if(s_loop==0)then
            field_angle0=field_angle
            n_match0=n_match1
            sig00=sig0
            s_loop0=s_loop
            obj_x0=obj_x
            obj_y0=obj_y
        endif
        if(n_match1.gt.n_match0.and.abs(sig0).lt.0.3.and.abs(sig0).gt.0.001)then
            field_angle0=field_angle
            n_match0=n_match1
            sig00=sig0
            s_loop0=s_loop
            obj_x0=obj_x
            obj_y0=obj_y
            write(*,'(a,f7.2,i7,f7.3)')'   field_angle/n_match1/sig0:',field_angle0,n_match0,sig00!大于前一天匹配数量或满足上面要求才输出
        endif
!032--基于初始旋转角度完成匹配后，判断其匹配sigma和匹配数量
!032--若sig0>0.5,或者match<6,增加旋转角度，循环加1，重新调用子程序，直到匹配结果满足条件
!032--以此确定旋转角度，后续图像均使用此角度。
!n_stop用来记录匹配上后求精确角度的范围，避免浪费时间，设定20度，认为如果匹配星超过6颗则认为开始匹配上
      if(s_loop .lt.((360/angle_step)+3).and.n_stop.le.6)then!保证循环360度且需要粗匹配后的精细匹配时有6次
       ! if(sig0 .gt.0.5.or.n_match1.lt.6) then
          if(abs(sig00) .lt.1.0.and.abs(sig00).gt.0.001.and.n_match0.gt.6) then !确定粗匹配成功，调小匹配限和步长，进行细匹配
              n_stop=n_stop+1
              if(n_stop==1)then
                  sig00=0
                  n_match0=0
                field_angle=field_angle-3 !如果匹配满足以上条件，将角度调小3度，进行细匹配，这里减3，后面加2，相当于减1
                angle_step=1.0
               limit_match0=limit_match 
              endif
          endif
          field_angle=field_angle+angle_step 
            if(field_angle.gt.360.0)then
                field_angle=field_angle-360.0
            endif
         !  print*,'sloop-field_angle',s_loop,field_angle
            s_loop=s_loop+1
         !   print*,trim(objfitsfile),field_angle,pscale,fl,limit_match,obj_x,obj_y,n_match1,sig0
            goto 123
        !endif
      endif

!033--输出确认旋转角度的信息
        write(*,'(a)')          '================Confirm the angle=================='     
        write(*,'(a,a)')        'This image:',trim(objfitsfile)
        if(n_match0.lt.3) write(*,'(a,i3)')'参考星少(<3) n_match=',n_match0
        if(abs(sig00).gt.0.5)  write(*,'(a,f7.3)')'匹配误差大(>0.5) sig0=',sig00
        write(*,'(a,2f10.3)')   'The object position:',obj_x0,obj_y0
        write(*,'(a,i3,f10.3)') 'Matched stars/sigma:',n_match0,sig00 
        write(*,'(a,f7.3)') '          The angle:',field_angle0
        
        if(nostar==1.or.abs(sig00)<0.001.or.abs(sig00)>0.3.or.n_match0<3)then
            if(n_confirm.le.3)then
                write(*,'(i5,a)')  n_confirm,'--上图确定角度失败，使用下一幅确定角度...'
                n_confirm=n_confirm+1
                goto 333!若第一幅图的reg文件内没有星象，直接读取第二幅
            else
                write(*,'(a)')' 本日图像旋转角度确定失败！！请检查：'
                write(*,'(a)')'         1) 星表、历表是否正确；2）本日前5幅图像是否正常！' 
                write(*,'(a)')' ===========使用前一天旋转角度进行今日匹配============='
                field_angle0=field_angle_bak
               ! angle_step=5.0
               ! goto 333
               ! if(n_confirm.gt.6)then
               !    close(11)!关闭今日图像列表文件
                   n_day_lose=n_day_lose+1
                !   goto 77
               ! endif
            endif
        endif     
        write(*,*)
        
        write(*,'(a)')          '================Processing ......=================='
      close(11)!关闭lst文件，正式处理时重新打开
      
!04--!!!!!!!======================利用以上旋转角度计算本日所有数据。
            field_angle=field_angle0
!041--建立结果输出文件（object.out）
!==============输出所有参考星，方便长跨度参考星数量归算及背景星等
100 write(*,'(a)')'       使用配置文件中已设置的旋转角度进行匹配。'
      write(obj_label,'(i1)')obj
      refoutfile=trim(newoutfile0)//'\'//'00ref_all_'//trim(obj_label)//'.out'
      open(18,file=refoutfile,status='unknown',access='append')     

222  objoutfile=trim(fitspath)//'\object_'//trim(obj_label)//'.out'
!042--!按照列表逐幅进行处理，自动检测目标
      open(7,file=objoutfile,status='replace')        
      open(11,file=trim(fitspath)//'\'//'fits.lst')  
      n_objpos=0!统计本程序参与匹配的图像数
      n_lose=0!统计因图像质量问题（）无法匹配的图像数目
      
!↓↓↓↓↓↓↓↓开始按顺序处理图像
       do i=1,m
555     read(11,*,end=9991) objfitsfile
        n_objpos=n_objpos+1
		objfitsfile=trim(objfitsfile)
      if(trim(tele_label)=='ss156')then
        call fits_header_156(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='ss156_2014')then
        call fits_header_156_2014(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='ss156_2016')then
        call fits_header_156_2016(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='km100')then
        call fits_header_km100(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='km100B')then
        call fits_header_km100B(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='lj240')then
        call fits_header_lj240(objfitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        else
            write(*,*)'请检查fits头文件格式！'
        endif
!all the data must be in the same month  
        obj_T=day+(hh*3600+mm*60+ss+delta_t)/86400d0
        call sla_CALDJ (year,month,day,djm,jjj)
        djm=djm+(hh*3600+mm*60+ss)/86400d0
        calpm_epoch=year+((month-1)*30+day)/365d0!use to compute the proper motion
!interpolation the epheris to get the c-data   
        call ENLGR(eph_T,eph_ra,n_eph,obj_T,obj_ephra)
        call ENLGR(eph_T,eph_de,n_eph,obj_T,obj_ephde)
        !call interT(eph_T,eph_ra,n_eph,obj_T,obj_ephra)
        !call interT(eph_T,eph_de,n_eph,obj_T,obj_ephde)
!!************对参考星表计算自行-------------
        n=1
        do k=1,n_gaia 
            if(gaia_ra(k).gt.(obj_ephra-fsize/60.0).and.gaia_ra(k).lt.(obj_ephra+fsize/60.0)&
                &.and.gaia_de(k).gt.(obj_ephde-fsize/60.0).and.gaia_de(k).lt.(obj_ephde+fsize/60.0))then
              gaia_de0(n)=gaia_de(k)+(calpm_epoch-2016)*gaia_pd(k)/1000d0/3600d0  !注意自行的单位
              gaia_ra0(n)=gaia_ra(k)+(calpm_epoch-2016)*gaia_pr(k)/1000d0/3600d0/dcosd(gaia_de0(n))
              gaia_mag0(n)=gaia_mag(k)
              n=n+1
        endif
        enddo   
        n_gaia0=n-1
        write(*,'(a,i5)')  '   用于匹配恒星数：',n_gaia0             
        
        ref_x=0.0
        ref_y=0.0
        ref_ra=0.0
        ref_de=0.0
       ! field_angle=3.0
        !基于之前计算角度，进行自动寻找目标及匹配（目标归算使用迭代par）
        call find_obj_baseangle(objfitsfile,tele_label,field_angle,pscale,fl,limit_match,obj_ephra,obj_ephde,&
         & gaia_ra0,gaia_de0,gaia_mag0,n_gaia0,modeltype,&
         & obj_x,obj_y,obj_flux,snr0,obj_obsra,obj_obsde,n_match1,sig0,ref_x,ref_y,ref_flux,ref_ra,ref_de,ref_mag,par1,nostar)
        
        if(nostar==1) then !若图像reg文件中无星象或者少于3颗，继续读下一幅
             n_lose=n_lose+1
             goto 555
        endif
!!!!!!!!!!
!!              write(*,'(a,i2,a,a)') 'QQQ  ', obj, '     ',obj_label
              write(obj_label,'(i1)') obj  !Qiao 强制赋值
!!              write(*,'(a,i2,a,a)') 'qqq  ', obj, '     ',obj_label
!!             pause
        
!043--将匹配信息存入reg文件         
        open(6,file=trim(objfitsfile)//trim(obj_label)//'.ref.reg',status='replace') 
           write(6,'(a)')'global color=red font="helvetica 10 normal" selec&
           &t=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source' 
           write(6,'(a)') 'physical'
        do n=1,n_match1
			  write(6,'(a7,2f11.3,2f6.1,a5,2f15.8,f8.3,2f20.10,f15.3)')&
                  &    'ellipse',ref_x(n),ref_y(n),&   
                  &             15.0d0,15.0d0,'  #  ',ref_ra(n),ref_de(n),ref_mag(n)
        enddo
        close(6)
!输出不同目标的参考星        
!        open(16,file=trim(objfitsfile)//trim(obj_label)//'.ref.out',status='replace') 
!        do n=1,n_match1
!           x0=(ref_x(n)-obj_x)*pscale/fl!以目标为中心位置，计算AD
!           y0=(ref_y(n)-obj_y)*pscale/fl 
!           call xy2rade(par1(1:6),6,x0,y0,obj_ephra*q,obj_ephde*q,ra0,de0)
!         !   print*,ra0/q,de0/q
!         !   print*,'00'
!			  write(16,'(2f11.3,2f15.8,f8.3,2f11.3,2f15.8)')&
!                  &  ref_x(n),ref_y(n),ref_ra(n),ref_de(n),ref_mag(n),ra0/q,de0/q,(ra0/q-ref_ra(n))*3600.0,(de0/q-ref_de(n))*3600.0
!         !     write(*,'(2f11.3,2f15.8,f8.3,2f11.3,2f15.8)')&
!         !         &  ref_x(n),ref_y(n),ref_ra(n),ref_de(n),ref_mag(n),ra0/q,de0/q,(ra0/q-ref_ra(n))*3600.0,(de0/q-ref_de(n))*3600.0
!        enddo
!        close(16)
!将所有参考星输入到指定文件，并计算图像背景星等
!背景星等=每平方角秒所含像素灰度总和参与计算星等
!        call fits_data(objfitsfile,fitsdata,exptime,npx) 
!        call cal_b(fitsdata,npx,bkgd,bkgdsigma)
!        !156tele:0.635arcsec/pix
!        arcsec_pix=pscale/fl/q*3600 !每像素角秒值
!        bkgd_flux=(bkgd)*(1.0/arcsec_pix/arcsec_pix)
!        !bkgd_flux=(bkgd-2150.0-2100)*(1.0/arcsec_pix/arcsec_pix)
!        !bkgd_flux=(bkgd-2150.0)*(1.0/arcsec_pix/arcsec_pix)/exptime        
!      !  ref_flux=ref_flux/exptime
        
!        call cal_star_mag(n_match1,ref_flux,ref_mag,obj_flux,bkgd0_mag,bkgd_mag,n_new,std_mag)!目标星星等
!        do n=1,n_match1            
!			  write(18,'(i5,i3,f10.6,2f10.3,i5,2f20.10,5f15.3,i7,f10.3)')&
!                  &    year,month,obj_T,bkgdsigma,exptime,n_match1,ref_ra(n),ref_de(n),ref_flux(n),ref_mag(n),bkgd_flux,bkgd0_mag(n),bkgd_mag,n_new,std_mag
!        ! write(*,'(i5,i3,f10.6,2f10.3,2f20.10,5f15.3,i7,f10.3)')&
!         !         &    year,month,obj_T,bkgdsigma,exptime,ref_ra(n),ref_de(n),ref_flux(n),ref_mag(n),bkgd_flux,bkgd0_mag(n),bkgd_mag,n_new,std_mag
!   
!        enddo
!        close(18)
        
!044--屏幕输出该幅图像匹配的相关信息     
        write(*,'(i3.3,a,a)')n_objpos,'-',trim(objfitsfile)
        write(*,'(a,i4)')         '              Matched  stars：',n_match1 
      !  write(*,'(a,2f10.3)')         'Background magtitude/std_mag：',bkgd_mag,std_mag
        write(*,'(a,f7.3)')       '        Plate rotation angle：',field_angle
        write(*,'(a,2f10.3)')       '       Object position x/y：',obj_x,obj_y
        if(modeltype==6)then
            write(*,'(a,i2,a,6f10.4)')'       Plate constants-',modeltype,'-par:',par1(1:modeltype)
        elseif(modeltype==12)then
            write(*,'(a,i2,a,12f10.4)')'      Plate constants-',modeltype,'-par:',par1(1:modeltype)
        elseif(modeltype==20)then
            write(*,'(a,i2,a,20f10.4)')'      Plate constants-',modeltype,'-par:',par1(1:modeltype)
        elseif(modeltype==30)then
            write(*,'(a,i2,a,30f10.4)')'      Plate constants-',modeltype,'-par:',par1(1:modeltype)
        endif
        write(*,'(a,f7.3)')       ' Plate model reduction sigma：',sig0
        if(sig0<0.0001)write(*,'(a)')'  ========This image-Matching Failed!' 
        write(*,*)
!045--对于sig0<0.15的图像结果计算其O-C残差（resra已*cos(de)）
       ! if(sig0.lt.1.0.and.sig0.gt.0.0001)then
         if(abs(sig0).lt.1.0.and.abs(sig0).gt.0.0001)then!.and.snr0.gt.10)then  !20211217增加去重复星像后，不考虑信噪比，全部目标均可用。使用目标星的snr限制目标星的数量，而非以历表为标准。!Qiao

        obj_resde=(obj_obsde-obj_ephde)*3600.0d0
        obj_resra=(obj_obsra-obj_ephra)*3600.0d0*dcosd(obj_ephde)
!        write(*,'(a,4f7.3)')'目标位置',obj_obsra,obj_obsde,obj_ephra,obj_ephde
!         write(*,'(a,2f7.3)')'目标残差',obj_resra,obj_resde
!046--将观测位置变换格式，输出到object.out文件（未按时间排序）
        call am2hms(obj_obsra*60.0d0,obj_obsde*60.0d0,i1,i2,ff,cc,j1,j2,fff)
        if(abs(obj_resra).ge.100.0) obj_resra=100.0 !20220724修改，待检测，为什么当obj_obsra/de为0时也输出。
        if(abs(obj_resde).ge.100.0) obj_resde=100.0
        write(7,'(i6,1x,i2.2,f10.6,2i3,f8.4,1x,a1,i2.2,i3,f7.3,4f12.7,3f10.4,2i3,f6.2,f9.2,a,a,f7.2)')&
     &        year,month,obj_T,i1,i2,ff,cc,j1,j2,fff,obj_obsra,obj_obsde,obj_ephra,obj_ephde,&
     &        obj_resra,obj_resde,sig0,hh,mm,ss,exptime,' ',trim(objfitsfile),field_angle  
       else
           n_lose=n_lose+1
       endif
    enddo
!↑↑↑↑↑↑处理完图像
!complete one image match and o-c calculation,output *.ref.reg, output:one line data in object.out file.
9991 close(11)!objectoutfile
     close(7)!object.out
!05--输出匹配总结信息
      write(*,'(a)')        '================03Match-Program execution summary=================='
      write(*,'(a,i2,a,i2)')        '目前处理目标数：',obj,'/',obj_total
      write(*,'(a,i3)')     '   Total number of data in match-pro：',n_objpos
      write(*,'(a,i3)')     'The number of sig0<1.0 in match-pro：',n_objpos-n_lose
      write(*,'(a)')        'Output the reference stars''s infor to the files named by *.ref.reg.'
      write(*,'(a)')        'Output the object''s obs-data to the file named by object.out.'
      write(*,'(a)')        '==================================================================='
            
!--排序，对object.out文件中的信息按时间排序
      iii=0
      open(9,file=objoutfile)
     do i=1,m
        read(9,'(a300)',end=101) line(i)
        read(line(i),'(10x,f10.6,138x)') T(i)
        iii=iii+1
     enddo 
101  close(9)
     call sort_T(T,line,iii)  
   !   objoutfile=trim(fitspath)//'\object_sort.out'
      open(15,file=objoutfile,status='replace')
     do i=1,iii
        write(15,'(a300)') line(i)
     enddo 
     close(15)    
    if(obj.lt.obj_total)     goto 9999 !若小于总目标数，跳转，重新赋值参考星和历表文件路径
    write(*,'(a,i3)') '   本日匹配共归算天然卫星目标个数',obj
    write(*,'(a)')    '================================================================='
    field_angle_bak=field_angle0!将前一天旋转角作为经验值
!处理下一天数据
    goto 77
777 close(1)    
    write(*,'(a,i3)') '   本目标匹配结束，共匹配天数：',n_day
    write(*,'(a,i3)') '   角度匹配失败天数（请确认）：',n_day_lose
    write(*,'(a)')    '================================================================='  
    stop
999 close(11)
    write(*,*)'本日图像不足3幅！'
    end program match
!*********************************************************************************
    
    
!=========================subtoutine==========================  
     subroutine find_obj_baseangle(fitsfile,tele_label,field_angle,pscale,fl,limit_match00,obj_ephra,obj_ephde,&
         & gaia_ra0,gaia_de0,gaia_mag0,n_gaia0,modeltype,&
         & obj_x0,obj_y0,obj_flux0,snr0,obj_obsra,obj_obsde,n_match1,sig1,ref_x0,ref_y0,ref_flux0,ref_ra0,ref_de0,ref_mag0,par0,nostar)
     implicit none
      integer*4 year,month,day,hh,mm,modeltype,refused1,n_match1,obj_j
      integer*4 j,k,m,n,nostar,iflag_par,n_gaia0,irematch 
      integer*4 n_eph,n_gaia,n_det,n_match,jjj,selectflag
      integer*4 naxis1,naxis2,bitpix,iflag2,i1,s_loop,n_lose,n_match_use
      character*300 fitsfile,regfile
      character*100 tele_label
      parameter(m=900000)
      real*8 eph_T(m),eph_ra(m),eph_de(m),min_mag,max_mag,limit_match,snr0
      real*8 gaia_ra(m),gaia_de(m),gaia_pr(m),gaia_pd(m),gaia_mag(m),gaia_mag0(m)
      real*8 obj_x,obj_y,obj_ephra,obj_ephde,obj_T,obj_x0,obj_y0,obj_xn,obj_yn
      real*8 par00(30),par(30),par1(30),par0(30),bscale,bzero,gain,ss,exptime,flux,ff,fff
      real*8 det_x(m),det_y(m),det_ra(m),det_de(m),det_flux(m)
      real*8 ref_x(m),ref_y(m),ref_flux(m),ref_ra(m),ref_de(m),ref_mag(m),det_xn(m),det_yn(m),snr(m)
      real*8 ref_x0(m),ref_y0(m),ref_flux0(m),ref_ra0(m),ref_de0(m),gaia_ra0(m),gaia_de0(m),ref_mag0(m)
      real*8 calpm_epoch,rl,pi,q,fl,pscale,djm,field_angle,sig1,obj_obsra,obj_obsde
      real*8 sig0,VKSI1(m),VYIT1(m),weight1(m),p,T(m),ra00,de00,obj_flux0,limit_match00
      real*8 x_center,y_center,ra_center,de_center,pre_x,pre_y,xi,eta     
      
      nostar=0
      pi=4*datan(1.0d0)
      q=pi/180.0d0
      par=0
      par00(1)=1 
      par00(2)=0
      par00(3)=0
      par00(4)=0
      par00(5)=1
      par00(6)=0

      par(1)=dcos(field_angle*q)*par00(1)+dsin(field_angle*q)*par00(4)  
      par(2)=dcos(field_angle*q)*par00(2)+dsin(field_angle*q)*par00(5)  
      par(3)=0d0
      par(4)=-dsin(field_angle*q)*par00(1)+dcos(field_angle*q)*par00(4)  
      par(5)=-dsin(field_angle*q)*par00(2)+dcos(field_angle*q)*par00(5)  
      par(6)=0d0

      ref_x0=0
      ref_y0=0
      ref_ra0=0
      ref_de0=0
      det_x=0
      det_y=0
      det_flux=0

      if(trim(tele_label)=='ss156')then
        call fits_header_156(fitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='ss156_2014')then
        call fits_header_156_2014(fitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='ss156_2016')then
        call fits_header_156_2016(fitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='km100')then
        call fits_header_km100(fitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='km100B')then
        call fits_header_km100B(fitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        elseif(trim(tele_label)=='lj240')then
        call fits_header_lj240(fitsfile,naxis1,naxis2,bscale,&
      &  bzero,gain,exptime,year,month,day,hh,mm,ss)
        else
            write(*,*)'请检查fits头文件格式！'
        endif

!read the detection *.reg file to get xy
        regfile=trim(fitsfile)//'.reg'!5        
!open the regfile and read the stars-position to match the catalog 
         n_det=1
        det_x=0.0
        det_y=0.0        
        open(5,file=regfile,status='old')
        read(5,*)
        read(5,*)

        read(5,'(7x,2f11.3,17x,f21.4,f10.2)',end=1111) det_x(1),det_y(1),det_flux(1),snr(1) 
        do j=2,m
          read(5,'(7x,2f11.3,17x,f21.4,f10.2)',end=888) det_x(j),det_y(j),det_flux(j),snr(j)
          n_det=n_det+1
        enddo
      !逐个将探测目标作为目标星 
888   close(5)
      if(n_det.lt.3)then
          write(*,*) trim(fitsfile),'检测星少于3颗'
          nostar=1
          return
      endif
      
      n_match1=0
      sig1=0
      iflag_par=0
      j=0   
199   j=j+1
      obj_x=det_x(j)
      obj_y=det_y(j)
      par1=par
      selectflag=0 !记录底片归算次数
      irematch=0
     !归算视场中心位置
666  if(selectflag<1)then 
         limit_match=limit_match00!*2.0
     else
         limit_match=limit_match00/1.5
     endif

        det_xn=0.0
        det_yn=0.0
        det_ra=0.0
        det_de=0.0
        ref_x=0.0
        ref_y=0.0
        ref_ra=0.0
        ref_de=0.0
        n_match=0
        sig0=0
        
        do k=1,n_det  
            if(det_x(k).gt.0.and.det_y(k).gt.0)then
            !!使用目标星作为中心位置
               det_xn(k)=(det_x(k)-obj_x)*pscale/fl
               det_yn(k)=(det_y(k)-obj_y)*pscale/fl
               if(iflag_par==0)then
                 call xy2rade(par1(1:6),6,det_xn(k),det_yn(k),obj_ephra*q,obj_ephde*q,det_ra(k),det_de(k))
               else
                 call xy2rade(par1(1:modeltype),modeltype,det_xn(k),det_yn(k),obj_ephra*q,obj_ephde*q,det_ra(k),det_de(k))
               endif         
              do n=1,n_gaia0
                call CAL_RL(det_ra(k),det_de(k),gaia_ra0(n)*q,gaia_de0(n)*q,rl)  
                if(rl/q*3600.le.limit_match) then
                   n_match=n_match+1 
                   ref_x(n_match)=det_xn(k)
                   ref_y(n_match)=det_yn(k)
                   ref_flux(n_match)=det_flux(k)
                   ref_ra(n_match)=gaia_ra0(n)
                   ref_de(n_match)=gaia_de0(n)
                   ref_mag(n_match)=gaia_mag0(n) 
                endif            
              enddo 
            endif
        enddo
    !试验该星不是目标星，换下一颗
!use these ref stars to recompute the field constant 

        iflag_par=0
        if(j==1) n_match1=n_match
        if(j.gt.1.and.n_match.gt.3.and.(j.le.n_det.or.selectflag.gt.0))then
	        call sol_par(ref_x,ref_y,ref_ra,ref_de,n_match,obj_ephra,obj_ephde,par1,&
        &    modeltype,sig0,refused1,vksi1,vyit1,weight1,IFLAG2)     
            sig0=sig0/q*3600.0
           !查找本幅图像中匹配最多的为目标星n_match1为最多星数。
            if (n_match.gt.n_match1.and.sig0.lt.1.and.sig0.gt.0.001)then
                ref_x0=0
                ref_y0=0
                ref_ra0=0
                ref_de0=0
                ref_mag0=0
                par0=0
                
                n_match1=n_match
                sig1=sig0
                obj_x0=det_x(j)
                obj_y0=det_y(j)
                obj_flux0=det_flux(j)
                snr0=snr(j)
                ref_x0(1:n_match)=ref_x(1:n_match)*fl/pscale+obj_x0
                ref_y0(1:n_match)=ref_y(1:n_match)*fl/pscale+obj_y0
                ref_flux0(1:n_match)=ref_flux(1:n_match)
                ref_ra0(1:n_match)=ref_ra(1:n_match)
                ref_de0(1:n_match)=ref_de(1:n_match)
                ref_mag0(1:n_match)=ref_mag(1:n_match)
                par0=par1
            endif
        endif
    !检验下一颗是否是目标星                       
        if(j.le.n_det.and.selectflag.lt.1) goto 199
       
        det_xn=0.0
        det_yn=0.0
        ref_x=0.0
        ref_y=0.0   
!使用目标星xy和rade为中心，以及所有参考星xy和rade，进行重新底片常数归算，         
        ref_x=(ref_x0-obj_x0)*pscale/fl
        ref_y=(ref_y0-obj_y0)*pscale/fl 
	        call sol_par(ref_x,ref_y,ref_ra0,ref_de0,n_match1,obj_ephra,obj_ephde,par0,&
        &    modeltype,sig0,refused1,vksi1,vyit1,weight1,IFLAG2)                           
            sig1=sig0/q*3600.0
!其后，使用此参数解算参考星中心位置天球座标(x_cneter,y_center,ra_center,de_center)            
        x_center=sum(ref_x0)/n_match1
        y_center=sum(ref_y0)/n_match1
        obj_xn=(x_center-obj_x0)*pscale/fl
        obj_yn=(y_center-obj_y0)*pscale/fl 
        call xy2rade(par0(1:modeltype),modeltype,obj_xn,obj_yn,obj_ephra*q,obj_ephde*q,ra_center,de_center) 
!        write(*,'(a,4f11.3)')'参考星星座中心xyrade',x_center,y_center,ra_center/q,de_center/q
!        pause
!最后，使用中心位置xy和rade为中心，所有参考星xy和rade，进行重新底片常数归算，以及解算目标星位置。
        ref_x=(ref_x0-x_center)*pscale/fl
        ref_y=(ref_y0-y_center)*pscale/fl         
	        call sol_par(ref_x,ref_y,ref_ra0,ref_de0,n_match1,ra_center/q,de_center/q,par0,&
        &    modeltype,sig0,refused1,vksi1,vyit1,weight1,IFLAG2)                           
            sig1=sig0/q*3600.0
!计算目标星赤经赤纬位置    (不适用直接查找的结果，修正后， 使用后面算法，通过预报结果查找目标星的方法)    
!        obj_xn=(obj_x0-x_center)*pscale/fl
!        obj_yn=(obj_y0-y_center)*pscale/fl             
!         call xy2rade(par0(1:modeltype),modeltype,obj_xn,obj_yn,ra_center,de_center,obj_obsra,obj_obsde)
!        write(*,'(a,4f11.3)')'目标星中心xyrade',obj_x0,obj_y0,obj_obsra/q,obj_obsde/q    
!        pause
!计算出底片常数后，利用目标星预报赤经赤纬位置，反算其在底片上的xy位置
!利用此位置在所有检测星中去对比进而找出目标星，此方法对于多目标靠的很近时有效。
            pre_x=0
            pre_y=0
 !        write(*,'(a,4f11.3)')'参考星星座中心xyrade',x_center,y_center,ra_center/q,de_center/q
        call rade2ky(obj_ephra*q,obj_ephde*q,ra_center,de_center, xi, eta)
        call xieta2xy(xi,eta,par0,6,pre_x,pre_y)
        pre_x=pre_x*fl/pscale+x_center
        pre_y=pre_y*fl/pscale+y_center
!         write(*,'(a,4f11.3)')'目标星中心xyrade',pre_x,pre_y,obj_ephra,obj_ephde
        write(*,'(a,2f11.3)')'预报xy： ',pre_x,pre_y
!        pause

!  从xyfile中挑选天然卫星的实测位置
        call select_obj_xy(pre_x,pre_y,regfile,obj_x,obj_y,selectflag)   !先运行这个，准备xy初值，供cal_target_xy使用
!        call selectxy(sx,sy,tarfile,objx,objy,selectflag)  !然后运行这个
        write(*,'(a,2f11.3,i4)')'实测xy： ',obj_x,obj_y,selectflag
        if(selectflag.eq.0) goto 9962
        if(selectflag.eq.1) irematch=irematch+1
        if(selectflag.eq.1.and.irematch.lt.2)    goto 666
!计算目标星赤经赤纬位置        
        obj_xn=(obj_x-x_center)*pscale/fl
        obj_yn=(obj_y-y_center)*pscale/fl             
         call xy2rade(par0(1:modeltype),modeltype,obj_xn,obj_yn,ra_center,de_center,obj_obsra,obj_obsde)
         obj_obsra=obj_obsra/q
         obj_obsde=obj_obsde/q  !子函数输出目标的观测位置为degree为单位
!         write(*,'(a,4f11.3)')'参考星星座中心xyrade',x_center,y_center,ra_center/q,de_center/q
!         write(*,'(a,4f11.3)')'目标星中心xyrade',obj_x,obj_y,obj_obsra/q,obj_obsde/q
          
        return
9962    write(*,*)'此图像未找到与预报目标位置相近的目标' 
        
        return
!        print*,'pipei',n_match,sig0,sig1
!        write(*,'(a,a)')'This image:',trim(fitsfile)
!        if(n_match1.lt.6) write(*,'(a,i3)')'参考星少！n_match=',n_match1
!        if(sig1.gt.0.5)  write(*,'(a,f7.3)')'匹配误差大！sig0=',sig1
!        write(*,'(a,2f10.3)')'目标星位置',obj_x0,obj_y0
!        write(*,'(a,i3,f10.3)')'匹配星数/底片模型误差',n_match1,sig1
1111  nostar=1
      write(*,'(a,a)')'This image has no stars detected',fitsfile
    end subroutine find_obj_baseangle
 
    
      subroutine select_obj_xy(sx,sy,xyfile,objx,objy,selectflag)
!         function:select the x/y of the target from all the detected stars
!            input:sx,sy,regfile
!           output:objx,objy:检测星中的目标星的位置
!                   selectflag:遴选目标成功与否的标志    
      implicit none
      integer i,m,selectflag
      parameter(m=10000)
      real*8 sx,sy,objx,objy,x,y
      character*300 xyfile

      selectflag=0
      objx=0
      objy=0

      open(66,file=xyfile)
      read(66,*)
      read(66,*)
      do i=1,m
        read(66,'(7x,2f11.3)',end=999) x,y
        if(dabs(x-999d0).ge.1d-6.and.dabs(y-999d0).ge.1d-6) then
          if(dsqrt((x-sx)**2+(y-sy)**2).le.3d0) then
            objx=x
            objy=y
            selectflag=1
            goto 999
          endif
        endif
      enddo
999   close(66)          

      end subroutine select_obj_xy

      subroutine xieta2xy(xi,eta,par,nm,sx,sy)
      implicit none
      integer nm
      real*8 xi,eta,par(nm),sx,sy,a,b,c,d,e,f

      if(nm.eq.6) then

        a=par(1)
        b=par(2)
        c=par(3)
        d=par(4)
        e=par(5)
        f=par(6)

        sx=(e*xi-b*eta+b*f-c*e)/(a*e-b*d)
        sy=(d*xi-a*eta+a*f-c*d)/(b*d-a*e)
    
      else
        write(*,'(a)') ' 模型不是6参数模型！！！退出！'
        stop
      endif

      end subroutine xieta2xy
!
      subroutine rade2xieta(ra,de,ra0,de0,xi,eta)
!  
!  根据目标的赤经赤纬和节点的赤经赤纬，计算目标的理想坐标
!
!     输入：
!        RA,DEC           d           目标的赤经赤纬
!        RA0,DEC0         d           节点的赤经赤纬
!     输出：
!        XI,ETA           d           目标的理想坐标
!
!
      implicit none
      real*8 ra,de,ra0,de0,xi,eta,a,b,c
      a=dcos(de)*dsin(ra-ra0)
      b=dcos(de0)*dsin(de)-dsin(de0)*dcos(de)*dcos(ra-ra0)
      c=dsin(de0)*dsin(de)+dcos(de0)*dcos(de)*dcos(ra-ra0)     
      
      xi=a/c
      eta=b/c

      end subroutine rade2xieta

      SUBROUTINE CAL_RL(RA01,DEC01,RA02,DEC02,RL)
!       function:Calculate the great arc length of two equatorial coordinates
!          input: RA01,DEC01,RA02,DEC02
!         output:RL（unit arcdegree）
    	IMPLICIT NONE
    	REAL*8 RA01,DEC01,RA02,DEC02,RL

    	RL=DACOS(DSIN(DEC01)*DSIN(DEC02)+DCOS(DEC01)*DCOS(DEC02)*&
         &  DCOS(RA01-RA02))

    	RETURN
    	END SUBROUTINE CAL_RL

      subroutine xy2rade(par,nm,x,y,ra0,de0,ra,de)
         !fucntion:xy-->rade
         !input:par(6),nm=6,
         !      xy:xy(has corrected by x0y0(center)x=(x-x0)*pscale/fl)
         !      ra0,de0:the ra0,de0 of the image center position(arcdegree)
         !output：ra/de (correspoding to the center)
       implicit none
       integer nm
       real*8 par(nm),x,y,ra0,de0,ra,de,ksi,yit,pi,q,tand1,tand2
       
       if(nm.eq.6) then
    	  ksi=par(1)*x+par(2)*y+par(3)
	      yit=par(4)*x+par(5)*y+par(6)
       elseif(nm.eq.12) then
	     ksi=par(1)*x+par(2)*y+par(3)+par(4)*x**2+par(5)*x*y+par(6)*y**2
	     yit=par(7)*x+par(8)*y+par(9)+par(10)*x**2+par(11)*x*y+par(12)*y**2
       elseif(nm.eq.20) then
	     ksi=par(1)*x+par(2)*y+par(3)+par(4)*x**2+par(5)*x*y+par(6)*y**2+&
            & par(7)*x**3+par(8)*x**2*y+par(9)*x*y**2+par(10)*y**3
	     yit=par(11)*x+par(12)*y+par(13)+par(14)*x**2+par(15)*x*y+&
            & par(16)*y**2+par(17)*x**3+par(18)*x**2*y+par(19)*x*y**2+par(20)*y**3
       elseif(nm.eq.30) then
	     ksi=par(1)*x+par(2)*y+par(3)+par(4)*x**2+par(5)*x*y+par(6)*y**2+&
            & par(7)*x**3+par(8)*x**2*y+par(9)*x*y**2+par(10)*y**3&
            & +par(11)*y**4+par(12)*x**3*y+par(13)*x**2*y**2+par(14)*x*y**3+par(15)*y**4
	     yit=par(16)*x+par(17)*y+par(18)+par(19)*x**2+par(20)*x*y+&
            & par(21)*y**2+par(22)*x**3+par(23)*x**2*y+par(24)*x*y**2+par(25)*y**3&
            & +par(26)*y**4+par(27)*x**3*y+par(28)*x**2*y**2+par(29)*x*y**3+par(30)*y**4
       endif
       tand1=ksi/DCOS(de0)/(1.0d0-yit*DTAN(de0))
       ra=ra0+DATAN(tand1)
	    
       tand2=(yit+DTAN(de0))*DCOS(DATAN(tand1))/(1.0d0 &
             &	-yit*DTAN(de0)) 
       de=DATAN(tand2)

      end subroutine xy2rade     
      
	  subroutine sol_par(x,y,ra,de,n,ra0,de0,par,nm,sig0,nused,vksi,vyit,p,IFLAG1)
!    function:解算底片常数
!       input:x:星象量度坐标x
!             y:星象量度坐标y
!             ra:星象星表位置
!             de:星象星表位置
!              n:使用星象数量
!            ra0:底片中心位置（用于将xy归算至理想坐标）
!            de0:底片中心位置（用于将xy归算至理想坐标）
!             nm:使用底片常数模型常数数量
!     output:par:常数解算结果
!            sig0:底片模型参数解算误差（arcsec）
!            nused：求解最终参数使用的观测数据数量
!            vksi：ksi归算残差矩阵N1
!            vyit:yit归算残差矩阵N1
!            p：使用1与未使用0数据显示矩阵N1
!            IFLAG1：如果迭代次数超过30次，标记为1，认为虽然星够，但是迭代不满足收敛条件
	   implicit none
	   integer n,nm,i,nused,IFLAG1
	   real*8 x(n),y(n),ra(n),de(n),par(nm),sig0
	   real*8 a(2*n,nm),b(nm,1),c(2*n,1),ra0,de0,pi,q,w(3,3)
       real*8 ksi(n),yit(n),sigpar(nm),rmag(n),p(2*n),vksi(n),vyit(n)

	   pi=4d0*datan(1.0d0)
	   q=pi/180.0d0

	   a=0.0d0
	   b=0.0d0
	   c=0.0d0
!yu----rade2ky,星表位置转为理想坐标
!      call ra0de02w(ra0*q,de0*q,w)
!	   do i=1,n
!	      call radew2ky(ra(i)*q,de(i)*q,w,ksi(i),yit(i))  
!	   enddo
    
!zhang---rade2ky
!首先将赤经赤纬转为理想坐标，以便求出理想坐标与量度坐标关系
 	   do i=1,n
	      call rade2ky(ra(i)*q,de(i)*q,ra0*q,de0*q,ksi(i),yit(i))  
       enddo      
!建立量度坐标、理想坐标矩阵a,c
!可解参数4-6-12-20-30
       if(nm.eq.4) then
         do i=1,n
	       a(2*i-1,1)=x(i)
	       a(2*i-1,2)=y(i)
	       a(2*i-1,3)=1.0
	       a(2*i-1,4)=0.0

	       a(2*i,1)=y(i)
	       a(2*i,2)=-x(i)
	       a(2*i,3)=0.0
	       a(2*i,4)=1.0
      
	       c(2*i-1,1)=ksi(i) 
	       c(2*i,1)=yit(i) 
	     enddo
       elseif(nm.eq.6) then
	     do i=1,n
	       a(2*i-1,1)=x(i)
	       a(2*i-1,2)=y(i)
	       a(2*i-1,3)=1.0

	      a(2*i,4)=x(i)
          a(2*i,5)=y(i)
	      a(2*i,6)=1.0
      
	      c(2*i-1,1)=ksi(i) 
	      c(2*i,1)=yit(i)
	    enddo
       elseif(nm.eq.12) then
	    do i=1,n
	      a(2*i-1,1)=x(i)
	      a(2*i-1,2)=y(i)
	      a(2*i-1,3)=1.0
	      a(2*i-1,4)=x(i)**2
	      a(2*i-1,5)=x(i)*y(i)
	      a(2*i-1,6)=y(i)**2

	      a(2*i,7)=x(i)
	      a(2*i,8)=y(i)
	      a(2*i,9)=1.0
	      a(2*i,10)=x(i)**2
	      a(2*i,11)=x(i)*y(i)
	      a(2*i,12)=y(i)**2
      
	      c(2*i-1,1)=ksi(i) 
	      c(2*i,1)=yit(i) 
	    enddo
       elseif(nm.eq.20) then
	    do i=1,n
	      a(2*i-1,1)=x(i)
	      a(2*i-1,2)=y(i)
	      a(2*i-1,3)=1.0
	      a(2*i-1,4)=x(i)**2
	      a(2*i-1,5)=x(i)*y(i)
	      a(2*i-1,6)=y(i)**2
	      a(2*i-1,7)=x(i)**3
	      a(2*i-1,8)=x(i)**2*y(i)
	      a(2*i-1,9)=x(i)*y(i)**2
	      a(2*i-1,10)=y(i)**3

	      a(2*i,11)=x(i)
	      a(2*i,12)=y(i)
	      a(2*i,13)=1.0
	      a(2*i,14)=x(i)**2
	      a(2*i,15)=x(i)*y(i)
	      a(2*i,16)=y(i)**2
	      a(2*i,17)=x(i)**3
	      a(2*i,18)=x(i)**2*y(i)
	      a(2*i,19)=x(i)*y(i)**2
	      a(2*i,20)=y(i)**3
      
	      c(2*i-1,1)=ksi(i) 
	      c(2*i,1)=yit(i) 
	    enddo
       elseif(nm.eq.30) then
	    do i=1,n
	      a(2*i-1,1)=x(i)
	      a(2*i-1,2)=y(i)
	      a(2*i-1,3)=1.0
	      a(2*i-1,4)=x(i)**2
	      a(2*i-1,5)=x(i)*y(i)
	      a(2*i-1,6)=y(i)**2
	      a(2*i-1,7)=x(i)**3
	      a(2*i-1,8)=x(i)**2*y(i)
	      a(2*i-1,9)=x(i)*y(i)**2
	      a(2*i-1,10)=y(i)**3
	      a(2*i-1,11)=x(i)**4
	      a(2*i-1,12)=x(i)**3*y(i)
	      a(2*i-1,13)=x(i)**2*y(i)**2
	      a(2*i-1,14)=x(i)*y(i)**3
	      a(2*i-1,15)=y(i)**4

	      a(2*i,16)=x(i)
	      a(2*i,17)=y(i)
	      a(2*i,18)=1.0
	      a(2*i,19)=x(i)**2
	      a(2*i,20)=x(i)*y(i)
	      a(2*i,21)=y(i)**2
	      a(2*i,22)=x(i)**3
	      a(2*i,23)=x(i)**2*y(i)
	      a(2*i,24)=x(i)*y(i)**2
	      a(2*i,25)=y(i)**3
	      a(2*i,26)=x(i)**4
	      a(2*i,27)=x(i)**3*y(i)
	      a(2*i,28)=x(i)**2*y(i)**2
	      a(2*i,29)=x(i)*y(i)**3
	      a(2*i,30)=y(i)**4
      
	      c(2*i-1,1)=ksi(i) 
	      c(2*i,1)=yit(i) 
	    enddo
       endif	  
!使用最小二乘算法，求解b系数（nm大小）
	   call least(a,c,b,2*n,nm,2*n,sigpar,sig0,nused,rmag,p,vksi,vyit,IFLAG1)

       do i=1,nm
		par(i)=b(i,1)
	   enddo
    end subroutine sol_par
  
 
    subroutine rade2ky(ra,de,ra0,de0,ksi,yit) 
!    function:ra/de to ksi/yit
!       input:ref star infor:ra/de(in arc)
!             center ra0/de0(in arc)
!    output:ref star ksi/yit
!    by zhang 20190823
    implicit none 
    real*8 ra,de,ra0,de0,ksi,yit
    real*8 ksi01,ksi02,yit01,yit02
    
    ksi01=dcos(de)*dsin(ra-ra0)
    ksi02=dsin(de)*dsin(de0)+dcos(de)*dcos(de0)*dcos(ra-ra0)
    ksi=ksi01/ksi02
    yit01=dsin(de)*dcos(de0)-dcos(de)*dsin(de0)*dcos(ra-ra0)
    yit02=dsin(de)*dsin(de0)+dcos(de)*dcos(de0)*dcos(ra-ra0)
    yit=yit01/yit02
    end subroutine rade2ky

	SUBROUTINE LEAST(A,Y,X,N1,K,NE,ERROR,SITA2,nrefused,rmag,p,vksi,vyit,IFLAG1)
!    function:calculate the matrix X by least square
!       input:A:观测资料矩阵
!           Y:等号右侧矩阵
!           X:待解参数
!           N1:观测资料数量（方程数量）
!           K:待解参数数量
!           NE: 观测资料数量（方程数量）   
!    output:X:求解底片常数数值
!           ERROR:各个参数解算误差
!           SITA2:方程组结算误差（sqrt(sum(res**2)/nep-k)）??除号下是否应该直接是nep-1
!           nrefused：求解最终参数使用的观测数据数量
!           rmag：定义了，实际未用，应该是想要考虑亮度归算残差的
!           p：使用1与未使用0数据显示矩阵N1
!           vksi：ksi归算残差矩阵N1
!           vyit:yit归算残差矩阵N1
!           IFLAG1：如果迭代次数超过30次，标记为1，认为虽然星够，但是迭代不满足收敛条件
	 implicit real*8(a-h,o-z)
 	 REAL*8 A(N1,K),X(K,1),Y(N1,1),AT(K,NE),ATA(K,K),ATAAT(K,NE)
	 REAL*8 V(NE),V2(NE),SITA2,ERROR(K),P(NE),rmag(ne/2),vksi(ne/2),vyit(ne/2)
	 INTEGER IFLAG1

      parameter(SFA=2.6,SFB=0.01)
 
      IFLAG1=0
      pi=4*datan(1.0d0)
	  q=pi/180.0d0
      SITA2=0.0D0
      SITA2P=0.0D0
	  P=1.0
	  AT=0.0
	  NBIG=0
	  vksi=0.0d0
	  vyit=0.0d0
	  SITA20=0.0D0

	  ILOOP=0
666	  NEP=NE-NBIG*2
	  SITA2P=SITA2
	  ILOOP=ILOOP+1
	  SITA2=0.0

	  DO I=1,K
 		DO J=1,NE
			AT(I,J)=A(J,I)
	    ENDDO
	  ENDDO

	  DO I=1,NE
		DO J=1,K
			AT(J,I)=AT(J,I)*P(I)
		ENDDO
      ENDDO
!ATAX=ATY(最小二乘解超定方程组)
  	  CALL MUL(AT,A,ATA,K,NE,NE,K)
	  CALL MATINV(K,K,ATA)
      CALL MUL(ATA,AT,ATAAT,K,K,K,NE)
      CALL MUL(ATAAT,Y,X,K,NE,NE,1)
!此处x即待解参数，K个
      v=0.0
 	  DO 10 I=1,NE
		DO 11 J=1,K
			V(I)=V(I)+A(I,J)*X(J,1)
11		    CONTINUE
		    V(I)=V(I)-Y(I,1)
		    V2(I)=V(I)**2
		    SITA2=SITA2+V2(I)*P(I)
10    CONTINUE
!只要剩余星象数量不少于待解参数，就计算其std，继续迭代
	        IF(NEP.NE.K) THEN 
		       SITA2=DSQRT(SITA2/(NEP-K))
	        ELSE
		       SITA2=0.0
	        ENDIF
            IF(DABS(SITA2-SITA20).GE.SFB*SITA20.AND.ILOOP.LE.30) THEN
		       SITA20=SITA2
		       NBIG=0
		       P=1.0
		       DO I=2,NE,2
		         IF(DABS(V(I)).GE.SFA*SITA2.or.DABS(V(I-1)).GE.SFA*SITA2) THEN !若残差过大，则删除此星
				   P(I-1)=0.0
				   P(I)=0.0
				   NBIG=NBIG+1
		         ENDIF
		       ENDDO
		       GOTO 666
            ELSEIF(ILOOP.GT.30) THEN
	           IFLAG1=1
            ENDIF
	        do i=2,ne,2
		      vksi(i/2)=v(i-1)/q*3600.0
		      vyit(i/2)=v(i)/q*3600.0
	        enddo   
	        nrefused=ne/2-nbig
            DO 30 I=1,K
		       ERROR(I)=SITA2*SQRT(ATA(I,I))
30             CONTINUE
     END SUBROUTINE LEAST


	SUBROUTINE MUL(A1,A2,A12,K1,K2,K3,K4)
!    function:Two matrix multiplication
!    input:A1:matrix 1
!          A2:matrix 2        
!          K1/K2:A(K1,K2)
!          K3/K4:A(K3,K4)
!    output:A12:A12=A1*A2,A12(K1,K4)
    
	 implicit real*8(a-h,o-z)
      REAL*8 A1(K1,K2),A2(K3,K4),A12(K1,K4)
	  DO 30 I=1,K1
        DO 30 J=1,K4
           A12(I,J)=0.0
 30        CONTINUE
           DO 10 I=1,K1
		      DO 10 J=1,K4
			     L=1
 20			     A12(I,J)=A12(I,J)+A1(I,L)*A2(L,J)
                 L=L+1
			     IF(L.LE.K2) GOTO 20
 10	             CONTINUE
	 END SUBROUTINE MUL

	 SUBROUTINE MATINV (N,NUP,A)
!    function:Calculate the inverse matrix
!    input:N,NUP,A
!    output:A:new matrix
     
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(N,N),B(1000),X(1000),EX(1000),IK(1000),JK(1000)
  10    DET=1.0D0
  11    DO 100 K=1,N
        AMAX=0.0D0
  21    DO 30 I=K,N
        DO 30 J=K,N
  23    IF (DABS(AMAX)-DABS(A(I,J))) 24,24,30
  24    AMAX=A(I,J)
        IK(K)=I
        JK(K)=J
  30    CONTINUE
  31    IF (AMAX) 41,32,41
  32    DET=0.0D0
        GOTO 140
  41    I=IK(K)
        IF (I-K) 21,51,43
  43    DO 50 J=1,N
        SAVE=A(K,J)
        A(K,J)=A(I,J)
  50    A(I,J)=-SAVE
  51    J=JK(K)
        IF (J-K) 21,61,53
  53    DO 60 I=1,N
        SAVE=A(I,K)
        A(I,K)=A(I,J)
  60    A(I,J)=-SAVE
  61    DO 70 I=1,N
        IF (I-K) 63,70,63
  63    A(I,K)=-A(I,K)/AMAX
  70    CONTINUE
  71    DO 80 I=1,N
        DO 80 J=1,N
        IF (I-K) 74,80,74
  74    IF (J-K) 75,80,75
  75    A(I,J)=A(I,J)+A(I,K)*A(K,J)
  80    CONTINUE
  81    DO 90 J=1,N
        IF (J-K) 83,90,83
  83    A(K,J)=A(K,J)/AMAX
  90    CONTINUE
        A(K,K)=1.0D0/AMAX
 100    DET=DET*AMAX
 101    DO 130 L=1,N
        K=N-L+1
        J=IK(K)
        IF (J-K) 111,111,105
 105    DO 110 I=1,N
        SAVE=A(I,K)
        A(I,K)=-A(I,J)
 110    A(I,J)=SAVE
 111    I=JK(K)
        IF (I-K) 130,130,113
 113    DO 120 J=1,N
        SAVE=A(K,J)
        A(K,J)=-A(I,J)
 120    A(I,J)=SAVE
 130    CONTINUE
 200	FORMAT (1X,'A(I,J)=',3F10.5,' B(I)=',F10.5)
        DO 138 J=1,N
        X(J)=0.D0
        DO 138 I=1,N
 138    X(J)=X(J)+A(I,J)*B(I)
 140	RETURN
    END SUBROUTINE MATINV

	subroutine am2hms(ra,de,hh,mm,ss,sign,ddd,dmm,dss)
!    function:arcminite -->hh,mm,ss
!       input:c1:ra(in degeree) de(in degree)
!      output:hh,mm,ss:new form of ra
!             flag:de's sign
!             ddd,dmm,dss:new form of de
!c     角分转换至 时分秒   和    °′″
	implicit real*8(a-h,o-z)
	real*8 ra,de,de1,ss,dss
    integer*4 hh,mm,ddd,dmm
	character*100 sign

	hh=int(ra/60.0/15.0)!分
	mm=int(ra-hh*60.0*15.0)/15.0!秒
	ss=(ra-hh*60.0*15.0-mm*15)*60.0/15.0

	if(de.ge.0) sign='+'
	if(de.lt.0) sign='-'

	de1=abs(de)
	ddd=int(de1/60.0)
	dmm=int(de1-ddd*60.0)
	dss=(de1-ddd*60.0-dmm)*60.0

	end
    
    
    SUBROUTINE sla_CALDJ (IY, IM, ID, DJM, J)
!*     - - - - - -
!*      C A L D J
!*     - - - - - -
!*  Gregorian Calendar to Modified Julian Date
!*
!*  (Includes century default feature:  use sla_CLDJ for years
!*   before 100AD.)
!*
!*  Given:
!*     IY,IM,ID     int    year, month, day in Gregorian calendar
!*
!*  Returned:
!*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
!*     J            int    status:
!*                           0 = OK
!*                           1 = bad year   (MJD not computed)
!*                           2 = bad month  (MJD not computed)
!*                           3 = bad day    (MJD computed)
!*
!*  Acceptable years are 00-49, interpreted as 2000-2049,
!*                       50-99,     "       "  1950-1999,
!*                       100 upwards, interpreted literally.
!*
!*  Called:  sla_CLDJ
!*  P.T.Wallace   Starlink   November 1985
!*  Copyright (C) 1995 Rutherford Appleton Laboratory
!*  License:
!*    This program is free software; you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation; either version 2 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program (see SLA_CONDITIONS); if not, write to the
!*    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
!*    Boston, MA  02110-1301  USA

      IMPLICIT NONE
      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J
      INTEGER NY
!*  Default century if appropriate
      IF (IY.GE.0.AND.IY.LE.49) THEN
         NY=IY+2000
      ELSE IF (IY.GE.50.AND.IY.LE.99) THEN
         NY=IY+1900
      ELSE
         NY=IY
      END IF
!*  Modified Julian Date
      CALL sla_CLDJ(NY,IM,ID,DJM,J)

      END SUBROUTINE sla_CALDJ



      SUBROUTINE sla_CLDJ (IY, IM, ID, DJM, J)
!*+
!*     - - - - -
!*      C L D J
!*     - - - - -
!*
!*  Gregorian Calendar to Modified Julian Date
!*
!*  Given:
!*     IY,IM,ID     int    year, month, day in Gregorian calendar
!*
!*  Returned:
!*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
!*     J            int    status:
!*                           0 = OK
!*                           1 = bad year   (MJD not computed)
!*                           2 = bad month  (MJD not computed)
!*                           3 = bad day    (MJD computed)
!*
!*  The year must be -4699 (i.e. 4700BC) or later.
!*
!*  The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
!*
!*  Last revision:   27 July 2004
!*
!*  Copyright P.T.Wallace.  All rights reserved.
!*
!*  License:
!*    This program is free software; you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation; either version 2 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program (see SLA_CONDITIONS); if not, write to the
!*    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
!*    Boston, MA  02110-1301  USA
      IMPLICIT NONE

      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J

!*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31,28,31,30,31,30,31,31,30,31,30,31 /
!*  Preset status.
      J = 0
!*  Validate year.
      IF ( IY .LT. -4699 ) THEN
         J = 1
      ELSE
!*     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN
!*        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 )  MTAB(2) = 28
!*        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J=3
!*        Modified Julian Date.
            DJM = DBLE ( ( 1461 * ( IY - (12-IM)/10 + 4712 ) ) / 4&
     &               + ( 306 * MOD ( IM+9, 12 ) + 5 ) / 10&
     &               - ( 3 * ( ( IY - (12-IM)/10 + 4900 ) / 100 ) ) / 4&
     &               + ID - 2399904 )

!*        Bad month.
         ELSE
            J=2
         END IF
      END IF
    END SUBROUTINE sla_CLDJ

    subroutine interT(x,y,n,t,z)
   ! function: 
     implicit none
     integer i,n
     real*8 x(n),y(n),t,z

     i=1
188  if(t .lt. x(i))then
        z=y(i-1)+(y(i)-y(i-1))/(x(i)-x(i-1))*(t-x(i-1))  
        return
     endif
     i=i+1
     goto 188
    end subroutine interT
    
	SUBROUTINE ENLGR(X,Y,N,T,Z)
	DIMENSION X(N),Y(N)
	DOUBLE PRECISION X,Y,T,Z,S
	Z=0.0
	IF (N.LE.0) RETURN
	IF (N.EQ.1) THEN
	  Z=Y(1)
	  RETURN
	END IF
	IF (N.EQ.2) THEN
	  Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
	  RETURN
	END IF
	I=1
10	IF (X(I).LT.T) THEN
	  I=I+1
	  IF (I.LE.N) GOTO 10
	END IF
	K=I-4
	IF (K.LT.1) K=1
	M=I+3
	IF (M.GT.N) M=N

	DO 30 I=K,M
	  S=1.0
	  DO 20 J=K,M
	    IF (J.NE.I) THEN
	      S=S*(T-X(J))/(X(I)-X(J))
	    END IF
20	  CONTINUE
	  Z=Z+S*Y(I)
30	CONTINUE
	RETURN
    END
    
    subroutine sort_T(t,line,n)
    implicit real*8(a-h,o-z)
    integer i,j,k,n
	real*8 t(n),tempt
	character*300 line(n),templine

	do i=1,n-1
		k=i
		do j=i+1,n
			if(t(j).le.t(k)) then
				k=j
			endif
		enddo
		tempt=t(i)
		t(i)=t(k)
		t(k)=tempt

		templine=line(i)
		line(i)=line(k)
		line(k)=templine
	enddo

    end subroutine sort_T
    
    subroutine cal_star_mag(n0,ref_flux,ref_mag,star_flux,star0_mag,star_mag,n_new,std_mag)
    !star0_mag: 每颗参考星对应的星等
    !star_mag:使用多颗去野值平均后的星等
    implicit none
    integer*4 n0,i,nn,n_new,iloop
    real*8 ref_flux(n0),ref_mag(n0),star0_mag(n0),star00_mag(n0),star_flux,star_mag
    real*8 ref_flux1(n0),ref_mag1(n0),star1_mag(n0),ocmag(n0)
    real*8 mean_mag,std_mag

    do i=1,n0
       star0_mag(i)=ref_mag(i)-2.5*log(star_flux/ref_flux(i))
    enddo
    !对以上计算的星等序列进行整理归算。
    mean_mag=sum(star0_mag)/n0
    do i=1,n0
        ocmag(i)=(star0_mag(i)-mean_mag)**2
    enddo
    std_mag=sqrt(sum(ocmag)/(n0-1))
    !print*,'计算星等',star0_mag(1:n0)
    !print*,'pingjunzhu',mean_mag
    !print*,'原始星归算星等残差',std_mag
    nn=n0
    star00_mag=star0_mag
    iloop=0
533 n_new=1
    iloop=iloop+1
    star1_mag=0.0
    ocmag=0
    do i=1,nn
        if((abs(star00_mag(i)-mean_mag).lt.2.0*std_mag))then
            star1_mag(n_new)=star00_mag(i)
            n_new=n_new+1
        endif
    enddo
    n_new=n_new-1
    mean_mag=sum(star1_mag)/n_new
    do i=1,n_new
        ocmag(i)=(star1_mag(i)-mean_mag)**2
    enddo
    std_mag=sqrt(sum(ocmag)/(n_new-1))
    if(nn.gt.n_new)then
        iloop=iloop+1
        nn=n_new
        star00_mag=star1_mag
        goto 533
    endif
   ! print*,'mean,std',mean_mag,std_mag
    star_mag=sum(star1_mag)/n_new
    !print*,'原始星-使用星-循环次数',n0,n_new,iloop
    
    end subroutine cal_star_mag