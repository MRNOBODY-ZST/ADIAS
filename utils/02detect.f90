! **********************************CCD Image preprocessing*************************************
!program name:02detect.f90
! function:1�����������󣬲�����������λ��
!          2) ����*.reg�ļ����������������λ�ú�������Ϣ
!    input:0��(ע�����ͼ��ͬ������fits_headerʱ��ע��ʱ�䡣Ŀǰ�Ѿ��޸ģ������׼ͷ�������������ļ��ع⿪ʼ����ʱ�̵ĸ�ʽ)
!          1) �����ļ���Ϣ��adias.cfg��
!             1.1�����������ֵϵ����2bkgd_threshold=5.0
!             1.2���������ֵ��2snr_threshold=5.0 
!             1.3���������Ķ�λ�ķ�����2pos_method=2��1-�����أ�2-���������أ�3-���������أ�
!   output:1����������ļ���*.reg
    
!  �ӳ���1��fits_header(fitsfile,bitpix,naxis1,naxis2,bscale,bzero,gain,exptime,year,month,day,hh,mm,ss)
!                      ��ȡfitsͼ��ͷ�ļ���Ϣ��ע�ⲻͬͼ��ͷ�ļ�ʱ����Ϣ�Ĳ�ͬ��
!          2) fits_data(fitsfile,fitsdata,exptime,npx)
!                      ��ȡfitsͼ�����ݣ������ӶԷ��������صĶ�ȡ��
!          3) detect_xzj(fitsfile,bkgd_threshold,sumi,starx,stary,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma)
!                      ��ͨ��������+�����ض���������
!          4) detect_guass(fitsfile,bkgd_threshold,sumi,starxx,staryy,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma)
!                      ��ͨ��������+��˹����������
!          5) cal_b(mdat,npx,bkgd,bkgdsigma)
!                      ����fitsͼ�񱳾�ֵ����sigma���ѿ۳�����Ӱ�죩
!          6) sorting_7terms(x01,x02,x03,x04,x05,x06,x07,n)
!                      ��x01Ϊ�ο��������н�������
!          7) seq2mat(seqnum,naxis1,naxis2,nli,nar) 
!                      һά�����Ϊ��ά����
!          8) mat2seq(nli,nar,naxis1,naxis2,seq) 
!                      ��ά�����Ϊһά����
!          9)g��˹����ӳ��򣨴����䣩Ŀǰ���μ����˹λ�ÿ��ã����Ȳ��������ء���ε������Ż�
!ע�⣺����֮ǰ��pre������ʹ���˳�������ͼ��ķ�����Ԥ��������Ľ���������ڲ��Ŀ��ʹ��
! by zhy 2019.09.02
! V2.2 add the enhance function by zhy 2019.09.24
! V3.0  Chinese comment by zhy 2020.01.29
!20200212�޸�����ȹ�ʽ��ʹ�������=��������-����ֵ/�����
!20200216�޸�detect_xzj�ӳ��򣬽��������Ϊԭʼͼ�����ȣ�û�м������������ݣ�������match���ּ���ҹ���
!20200216�������Ϊԭʼͼ��-����������bias����ֵ
      program detect
       implicit real*8(a-h,o-z)
       parameter(m=1000)
       integer pos_method
	   real*8 bkgd_threshold,snr_threshold
	   integer*4 ic,ic1,nstar,i,n_fits,n_outstar,flatflag,superflag
       character*200 tempcha,fitspath,fitsfile0,fitsfile,xyfile
       real*8 gain,bkgd,bkgdsigma
	   real*8 sumi(50000),starx(50000),stary(50000),snr(50000)
       integer*4 star_id(50000),star_pix(50000),overflag(50000)
       integer*4 n_day
       character*100 tele_label

!01--�������ļ�       
!      open(1,file='fitspath.in',status='old')
        open(1,file='D:\00adias\exe\fitspath.in',status='old')
       n_day=0
77     read(1,'(a100)',end=777) fitspath
       fitspath=trim(fitspath)
       n_day=n_day+1  !���Լ�¼������������
    write(*,'(i3,a,a,a)')  n_day,'-','����fits�ļ���',trim(fitspath)
    
!   ����fitspath.in�������ļ����ݿ�ʼ���մ���ͼ��   
!     open(11,file='F:\Adias_jpl\02detect\Console1\adias.cfg')
     open(11,file='D:\00adias\exe\adias.cfg')
      do i=1,m
	      read(11,'(a100)',end=1000) tempcha

        ic=index(tempcha,'1flatflag=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+1),*) flatflag
        endif

        ic=index(tempcha,'1superflag=')
        if(ic.gt.0) then
          ic1=index(tempcha,'=')
          read(tempcha(ic1+1:ic1+1),*) superflag
        endif
          
	      ic=index(tempcha,'2bkgd_threshold=') 
	      if(ic.gt.0) then
  	      ic1=index(tempcha,'=')
	        read(tempcha(ic1+1:ic1+10),*) bkgd_threshold
          endif
          
	      ic=index(tempcha,'2snr_threshold=') 
	      if(ic.gt.0) then
  	      ic1=index(tempcha,'=')
	        read(tempcha(ic1+1:ic1+10),*) snr_threshold
          endif   
        
	      ic=index(tempcha,'2pos_method=') 
	      if(ic.gt.0) then
  	      ic1=index(tempcha,'=')
	        read(tempcha(ic1+1:ic1+10),*) pos_method
          endif          

      enddo
1000  close(11)		
!      bkgd_threshold=3.0    
!      snr_threshold=10.0
!02--�����ȡͼ��*_n.fit��������������
      open(11,file=trim(fitspath)//'\'//'fits.lst')
      do n_fits=1,m   
        read(11,'(a100)',end=999) fitsfile0	  
           write(*,'(i3.3,2a)') n_fits,'-',trim(fitsfile0)
           xyfile=trim(fitsfile0)//'.reg'
           open(15,file=xyfile)
           write(15,'(a)')'global color=green font="helvetica 10 normal" selec&
           &t=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source' 
           write(15,'(a)') 'physical'
!���ȥ���ظ�����ʱ�õĵȼ۱�
!      open(12,file='getfile.bat',status='replace')
 !     write(12,'(2a)') 'cd '
!      write(12,'(a)') 'del F:\00Adias\02detect\02detect\fort.333'
!      write(12,'(a)') 'del F:\00Adias\02detect\02detect\fort.334'
!      write(12,'(a)') 'del F:\00Adias\02detect\02detect\fort.335'
!      write(12,'(a)') 'del F:\00Adias\02detect\02detect\fort.336'
!      close(12)
!      call system('getfile.bat')
!==========��ͨ����+�����ط�������������
           call detect_xzj_yy(flatflag,superflag,fitsfile0,bkgd_threshold,pos_method,sumi,starx,stary,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma) 
          ! call detect_xzj_zz(flatflag,superflag,fitsfile0,bkgd_threshold,pos_method,sumi,starx,stary,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma) 
!==========��ͨ����+��˹��������������(�ݲ�����)
           !call detect_gauss(fitsfile,bkgd_threshold,sumi,starx,stary,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma)
           n_outstar=0!ͳ������snr����������
           do i=1,nstar
             if(snr(i).gt.snr_threshold .and.star_pix(i).gt.5)then
                write(15,'(a7,2f11.3,2f6.1,a5,f21.4,f10.2,3i5,2f15.3)')& 
                &     'ellipse',starx(i),stary(i),real(10),real(10),&
                &     '  #  ',sumi(i),snr(i),star_id(i),star_pix(i),overflag(i),bkgd,bkgdsigma
                n_outstar=n_outstar+1
             endif
           enddo  
           write(*,'(a,2f10.3)')  '           Image-bkgd,sigma: ',bkgd,bkgdsigma          
           write(*,'(a,i4)'),     'The stars��Satisfy the SNR and pix_num condition��: ',n_outstar
           write(*,*)
           close(15)            
      enddo
999   close(11)
!03--��Ļ�����Ϣ
      write(*,'(a)')   '================02Detect-Program execution summary=================='
      write(*,'(a,i4)')'    ���������(*.reg)�������ͼ��',n_fits-1
      write(*,'(a)')   '    Output the stars'' information to the files named by *.reg.'
      write(*,'(a)')   '===================================================================='    
!������һ������
    goto 77
777 close(1)
    write(*,'(a,i3)') '�������������������ͼ��������',n_day
    write(*,'(a)')    '================================================================='
    end program detect   
!*********************************************************************************
    
!=========================subtoutine==========================               
      subroutine detect_xzj_yy(flatflag,superflag,fitsfile0,bkgd_threshold,pos_method,sumi_real,starx,stary,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma)
!       function:detect stars and output the stars center
!          input:fitsfile:fits file name
!                bkgd_threshold:fitsdata-bkgd_threshold*bkgdsigma
!         output:sumi:the stars' flux
!                starx:the x-position of star center
!                stary:the x-position of star center
!                snr:the stars' snr
!                star_id:the id of the star( already sorted by brightness)
!                star_pix:The number of stars contained in star
!                overflag:Whether the stars overflow
!                nstar:Number of detected stars
!                bkgd:the image bkgd
!                bkgdsigma:the image bkgd
!Note:���Ŀ������=naxis1*naxis2/pi/6**2  �����������=pi*15**2, minpix=3
      implicit none
      character*200 fitsfile,fitsfile0 
      integer pos_method
      integer*4 i,j,naxis1,naxis2,npx,bitpix,s,seq0,seq1,k,flatflag,superflag
      integer*4 maxobj,minpix,maxpix,nli,nar,nstar,pos_nstar,npx0
      real*8 bzero,bscale,gain,bkgd,bkgdsigma,bkgd_threshold,bkgd0,bkgdsigma0
      parameter(s=9000000)
      real*8 mdat(s),mdat01(s),abox(3000,3000),lty(4),mdat0(s),abox0(3000,3000)
	  real*8 sumx(50000),sumy(50000),maxflux
      real*8 starx(50000),stary(50000),sumi(50000),pi,exptime,ss,sumi_real(50000)
      integer*4 star_id(50000),star_pix(50000),overflag(50000)
	  real*8 snr(50000),snr01(50000)  
      integer*4 idpix(s),numpix(s),objpos(5000,5000)
      integer*4 numobj,ic,ic1
      
      integer idbox(3000,3000),idbox2(3000,3000),nr,ni1,ni2,nj1,nj2,kread,iflag,i1,j1,neid,eid1(50000),eid2(50000),tempeid1,tempeid2
      character*200  tempcha1(50000),tempchap 
      
!      nstar=0
      pi=4*datan(1.0d0)
      bkgd=0.0d0
      bkgdsigma=0.0d0
      idbox2=0
      idbox=0
      if((flatflag+superflag).lt.0.5)then
        fitsfile=fitsfile0
        write(*,'(a,i5)') 'ʹ��ԭʼͼ�����ͼ���⣡'
      else
        ic=index(fitsfile0,'.fit') 
        if(ic.gt.0) then
           ic1=index(fitsfile0,'.fit') 	    
 	       fitsfile0=trim(fitsfile0)
	       fitsfile=fitsfile0(1:ic1-1)//'_n.fit'
        endif
      endif
!��ȡfitsͷ�ļ���������Ϣ����ȷ�����Ŀ����maxobj���������������maxpix�����maxflux          
!���Ŀ�����������ӳ���ȫ�ǰ뾶Ϊ6��������ex: 1024*1056/3.1415/6/6=9561��ʵ����û����ô���
!����������������ǣ�����뾶Ϊ15������������������һ��û����ô����������뾶��5����
      call fits_header_main(fitsfile,bitpix,naxis1,naxis2,bscale,bzero,gain,exptime)

      if(dabs(gain).le.1e-9) gain=1.0    
      maxflux=2**bitpix-1
!      maxflux=62000  !tong yu
      maxobj=int(naxis1*naxis2/pi/6**2)        
      minpix=3
      maxpix=int(pi*15**2)
      call fits_data(fitsfile,mdat,exptime,npx) 
      call cal_b(mdat,npx,bkgd,bkgdsigma)   !the whole image 
     ! call fits_data(fitsfile0,mdat0,exptime,npx0) !Ϊ�˵���ԭʼͼ�����ȣ�Ŀǰ����
     ! call cal_b(mdat0,npx0,bkgd0,bkgdsigma0)   !Ϊ�˵���ԭʼͼ�����ȣ�Ŀǰ����      
      
!��ͨ����ǰ�����ݽ��д���
!each pixes ADU-��background��bkgd_threshold*bkgdsigma����bkgd_thresholdһ����5����ͼ��Ԥ�������������ǿ����sigma��С
!��δԤ����sigma��������Ϊ2.6����3
      do i=1,npx
        mdat01(i)=mdat(i)-(bkgd+bkgd_threshold*bkgdsigma)
        if(mdat01(i).le.0.0d0) mdat01(i)=0.0d0
      enddo 
      abox=0.0
      do i=1,npx
        call seq2mat(i,naxis1,naxis2,nli,nar)    
        abox(nli,nar)=mdat01(i)
        abox0(nli,nar)=mdat(i)
      enddo
!set the image edge pixel data =0
      do i=1,naxis1
        abox(i,1)=0.0d0
        abox(i,naxis2)=0.0d0
      enddo
      do j=1,naxis2
        abox(1,j)=0.0d0
        abox(naxis1,j)=0.0d0
      enddo
!!!!!������Ҫʱ�������	
!    write(*,'(a,3i5,i8)'),'maxobj/minpix/maxpix/npx:',maxobj,minpix,maxpix,npx
      
!��������ͨ��⣬��ÿ���ر��
!lty1=���ϣ�lty2=���ϣ�lty3=���ϣ�lty4=��
!Ŀ��ţ�numobj�����ر��idpix(numobj)=numobj,Ŀ�꺬������numpix(numobj)
!�������£���������������ؼ�⣬��ĳ����>0ʱ���ֱ��������жϣ�
!1)��4�ھ�<0���½�Ŀ��numobj+1����Ŀ��������numpix(numobj)+1�������ر��idpix(seq)=numobj
!2��lty1>0��1�ı�Ÿ���ǰ���أ�Ŀ��������+1
!3��lty2>0��2�ı�Ÿ���ǰ���أ�Ŀ��������+1
!4��lty3>0��3�ı�Ÿ���ǰ���أ�Ŀ��������+1
!5��lty4>0��4�ı�Ÿ���ǰ���أ�Ŀ��������+1
!maxobj�趨���Ŀ�������=ȫ����/�׾��������������������Ų������ʾ
!maxpixĿ�꺬���������=2���׾�Բ����,�������������ʾ
!����ͨ����һ�μ��ɣ�����Ҫ�ظ���ţ���������ظ����
!���ԭ����ɨ����һ��3���أ���������ͬ���ٿ����ͬ�У���������ͬ������4������
      numobj=0
      numpix=0
      idpix=0
      objpos=0
      do i=2,naxis1-1  !��
        do j=2,naxis2-1  !��
          !  print*,'i,j',i,j,abox(i,j)
          lty(1)=abox(i-1,j+1)
          lty(2)=abox(i-1,j)
          lty(3)=abox(i-1,j-1)
          lty(4)=abox(i,j-1)
          if(abox(i,j).gt.1d-9) then
            call mat2seq(i,j,naxis1,naxis2,seq0)
!����ֻ��һ�����أ������0
            if(dabs(lty(1)+lty(2)+lty(3)+lty(4)).le.1d-9) then
              numobj=numobj+1
              idbox(i,j)=numobj
              idpix(seq0)=numobj !���ر��
              numpix(numobj)=numpix(numobj)+1  !Ŀ���ۼ�����
              objpos(idpix(seq0),numpix(numobj))=seq0
              if(numobj.gt.maxobj) then
                !write(*,*)'Detect stars more than maxobj:',numobj,maxobj
                !write(*,*)'Warning:Image quality problem!'
                !pause
                !goto 953
              endif
              if(numpix(numobj).gt.maxpix) then
                !write(*,*)'Pixels of star more than maxpix:',numpix(numobj),maxpix
                !write(*,*)'Warning:please check the radius value (maybe too small)!'
                !pause
                !goto 952
              endif                  
! ����
            elseif(lty(1).gt.1d-9) then
              call mat2seq(i-1,j+1,naxis1,naxis2,seq1)            
	          idbox(i,j)=idbox(i-1,j+1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
                !write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
                !write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                !goto 952
              endif    
! ����
            elseif(lty(2).gt.1d-9) then
	          idbox(i,j)=idbox(i-1,j)
              call mat2seq(i-1,j,naxis1,naxis2,seq1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
               ! write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
                !write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                !goto 952
              endif     
! ����         
            elseif(lty(3).gt.1d-9) then
	          idbox(i,j)=idbox(i-1,j-1)
              call mat2seq(i-1,j-1,naxis1,naxis2,seq1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
               ! write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
               ! write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                !goto 952
              endif     
! ���      
            elseif(lty(4).gt.1d-9) then
	          idbox(i,j)=idbox(i,j-1)
              call mat2seq(i,j-1,naxis1,naxis2,seq1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
               ! write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
               ! write(*,*)'Warning:please check the radius1 value (maybe too small)!'
               !goto 952
              endif
            endif
          endif
952     enddo  !952 'this image includes super star(>pi*r^2)'
      enddo
! ========================������ͨ����
953 continue   
    
!      write(*,*) 'numobj',numobj
      nr=250
      idbox2=idbox
      do i=2,naxis2-1
	  do j=naxis1-1,2,-1

        if(abox(i,j).gt.1d-9.and.abox(i,j-1).gt.1d-9) then
	      if(idbox2(i,j).eq.idbox2(i,j-1)) then
			  goto 9522
	      else 
              ni1=i-nr
              if(ni1.le.1) ni1=1
              ni2=i+nr
              if(ni2.ge.naxis2) ni2=naxis2
              nj1=j-nr
              if(nj1.le.1) nj1=1
              nj2=j+nr
              if(nj2.ge.naxis1) nj2=naxis1
              
              k=1
              do i1=ni1,ni2
                do j1=nj1,nj2
                
                  if((idbox2(i1,j1).eq.idbox2(i,j-1)).and.(i1.ne.i.and.j1.ne.j-1)) then
                    idbox2(i1,j1)=idbox2(i,j)
                    k=k+1
                  endif
                enddo
              enddo       
              idbox2(i,j-1)=idbox2(i,j)
              
            endif	          
          endif	        
9522    enddo
      enddo

!  �����ȼ۱�2017.9.8
      open(333,file='fort.333',status='replace')
      open(334,file='fort.334',status='replace')
      open(335,file='fort.335',status='replace')
      open(336,file='fort.336',status='replace')

      do i=2,naxis2-1
        do j=2,naxis1-1
          if(abox(i,j).gt.1d-9) then
            
            if(idbox2(i,j).ne.idbox2(i,j-1).and.idbox2(i,j-1).gt.0) then
              if(idbox2(i,j).lt.idbox2(i,j-1)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i,j-1)
              else
                write(333,'(2i10)') idbox2(i,j-1),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i,j+1).and.idbox2(i,j+1).gt.0) then
              if(idbox2(i,j).lt.idbox2(i,j+1)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i,j+1)
              else
                write(333,'(2i10)') idbox2(i,j+1),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i-1,j-1).and.idbox2(i-1,j-1).gt.0) then
              if(idbox2(i,j).lt.idbox2(i-1,j-1)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i-1,j-1)
              else
                write(333,'(2i10)') idbox2(i-1,j-1),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i-1,j).and.idbox2(i-1,j).gt.0) then
              if(idbox2(i,j).lt.idbox2(i-1,j)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i-1,j)
              else
                write(333,'(2i10)') idbox2(i-1,j),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i-1,j+1).and.idbox2(i-1,j+1).gt.0) then
              if(idbox2(i,j).lt.idbox2(i-1,j+1)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i-1,j+1)
              else
                write(333,'(2i10)') idbox2(i-1,j+1),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i+1,j-1).and.idbox2(i+1,j-1).gt.0) then
              if(idbox2(i,j).lt.idbox2(i+1,j-1)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i+1,j-1)
              else
                write(333,'(2i10)') idbox2(i+1,j-1),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i+1,j).and.idbox2(i+1,j).gt.0) then
              if(idbox2(i,j).lt.idbox2(i+1,j)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i+1,j)
              else
                write(333,'(2i10)') idbox2(i+1,j),idbox2(i,j)
              endif 
            endif

            if(idbox2(i,j).ne.idbox2(i+1,j+1).and.idbox2(i+1,j+1).gt.0) then
              if(idbox2(i,j).lt.idbox2(i+1,j+1)) then
                write(333,'(2i10)') idbox2(i,j),idbox2(i+1,j+1)
              else
                write(333,'(2i10)') idbox2(i+1,j+1),idbox2(i,j)
              endif 
            endif
          endif
        enddo
      enddo

      rewind(333)

!  ȥ�ظ�1
      k=1
      kread=0
      read(333,'(a20)',end=93) tempcha1(k)
      kread=kread+1

      do i=2,naxis1*naxis2
        read(333,'(a20)',end=93) tempchap
        iflag=0
        do j=1,k
          if(tempchap.eq.tempcha1(j)) iflag=1
        enddo 

        if(iflag.eq.0) then
          k=k+1
          tempcha1(k)=tempchap
        endif
      enddo
93    close(333)

      if(kread.eq.0) then
        close(334)
        close(335)
        close(336)
        neid=0
        goto 88
      endif

      do i=1,k
        write(334,'(a20)')  tempcha1(i)
      enddo

      rewind(334)
      
      k=1
      read(334,*,end=92) eid1(k),eid2(k)
      write(335,'(2i10)') eid1(k),eid2(k)
      do i=2,naxis1*naxis2
        read(334,*,end=92) tempeid1,tempeid2
        do j=1,k
          if(tempeid1.eq.eid2(j)) then
            tempeid1=eid1(j)
            goto 91
          elseif(tempeid2.eq.eid2(j)) then
            tempeid2=tempeid1
            tempeid1=eid1(j)
            goto 91
          endif
        enddo
91      k=k+1
        eid1(k)=tempeid1
        eid2(k)=tempeid2
        write(335,'(2i10)') eid1(k),eid2(k)
      enddo
92    close(334)

!  ȥ�ظ�2
      rewind(335)
      k=1
      read(335,'(a20)',end=932) tempcha1(k)

      do i=2,naxis1*naxis2
        read(335,'(a20)',end=932) tempchap
        iflag=0
        do j=1,k
          if(tempchap.eq.tempcha1(j)) iflag=1
        enddo 

        if(iflag.eq.0) then
          k=k+1
          tempcha1(k)=tempchap
        endif
      enddo
932    close(335)

      do i=1,k
        write(336,'(a20)')  tempcha1(i)
      enddo

      rewind(336)
      do i=1,naxis1*naxis2
        read(336,*,end=90) eid1(i),eid2(i)
      enddo
90    neid=i-1
      close(336)
!      write(*,*) neid,eid2(1:neid)
88    continue
      
      	do i=2,naxis2-1
        do j=2,naxis1-1  
          if(abox(i,j).gt.1d-9) then
            do k=1,neid
              if(idbox2(i,j).eq.eid2(k)) idbox2(i,j)=eid1(k)
            enddo
          endif
        enddo
      enddo
    
    
!! ===================XZJ
!!compute the center of the object
! starx(k):x-value
! stary(k):y-value
! stari(k):flux-value��mdat01---has subtracted (bkgd+bkgd_threshold*bkgdsigma),not final value������
!ע�⣬�˴�ʹ�����ػҶ�ֵΪmdat01(i)=mdat(i)-(bkgd+bkgd_threshold*bkgdsigma)
      k=1
      sumx=0.0
      sumy=0.0 
      sumi=0.0
      sumi_real=0.0
      starx=0.0
      stary=0.0  
      numpix=0

!�߽�����������ı߸��۳�10������
      do i=10,naxis2-10
        do j=10,naxis1-10
          if(idbox2(i,j).gt.0) then !.and.abox0(i,j).lt.maxflux) then !20220718������̨1m��������--�ع������Ȼ������������
              sumx(idbox2(i,j))=sumx(idbox2(i,j))+j*abox(i,j)**pos_method
              sumy(idbox2(i,j))=sumy(idbox2(i,j))+i*abox(i,j)**pos_method
              sumi(idbox2(i,j))=sumi(idbox2(i,j))+abox(i,j)**pos_method
              sumi_real(idbox2(i,j))=sumi_real(idbox2(i,j))+abox(i,j)
              numpix(idbox2(i,j))=numpix(idbox2(i,j))+1
          endif  
          if(abox0(i,j).ge.maxflux) overflag(idbox2(i,j))=1
        enddo
      enddo
            
!!!!!!!!!!!!!!!!2021-8-31      
      k=1
      
      do i=1,numobj
        if(numpix(i).ge.minpix.and.numpix(i).le.maxpix.and.overflag(i).lt.1) then

	      starx(k)=sumx(i)/sumi(i)
	      stary(k)=sumy(i)/sumi(i)
          sumi_real(k)=sumi_real(i)
          star_pix(k)=numpix(i)
          star_id(k)=i
          snr(k)=sumi_real(i)/dsqrt(sumi_real(i)+numpix(i)*bkgdsigma**2)

          k=k+1
            
            
!          do j=1,numpix(i)
!            call seq2mat(objpos(i,j),naxis1,naxis2,nli,nar)  
!            if(mdat01(objpos(i,j)).gt.maxflux) overflag(k)=1
!            sumx(k)=sumx(k)+nar*mdat01(objpos(i,j))**pos_method
!            sumy(k)=sumy(k)+nli*mdat01(objpos(i,j))**pos_method
!            sumi(k)=sumi(k)+mdat01(objpos(i,j))**pos_method 
!           ! sumi_real(k)=sumi_real(k)+(mdat0(objpos(i,j))-bkgd0)!�����������
!           sumi_real(k)=sumi_real(k)+mdat(objpos(i,j))
!          enddo
!maxim software:snr=sumi/bkgdsima/sqrt(pi*radius1**2)
!yuyong:tempsnrflux=sumvalue/dsqrt(sumvalue+numpix1*(bkgdsigma**2))
!tang:snr=mdat(i)-bkgd/bkgdsigma
          !snr(k)=sumi_real(k)/(bkgdsigma0*numpix(i)) 
          !snr(k)=(sumi_real(k)-bkgd*numpix(i))/(bkgdsigma*numpix(i)) ��zhangֱ�Ӷ��壨�۳���biasӰ�죩
          ! sumi_real(k)=sumi_real(k)-numpix(i)*bkgd
          ! snr(k)=(sumi_real(k)-bkgd*numpix(i))/bkgdsigma 
!          snr(k)=sumi_real(k)/sqrt(sumi_real(k)+numpix(i)*(bkgdsigma**2)) !yuyong ��������������������ɷֲ���
!          starx(k)=sumx(k)/sumi(k)
!          stary(k)=sumy(k)/sumi(k)
!          star_id(k)=i
!          star_pix(k)=numpix(i)
!          k=k+1
        endif
      enddo 
      nstar=k-1
!�����ض��������Ľ���
!��������������
      call sorting_7terms(sumi_real,starx,stary,snr,star_id,star_pix,overflag,nstar)
      
      call system('del fort.333')
      call system('del fort.334')
      call system('del fort.335')
      call system('del fort.336')
      
    end subroutine detect_xzj_yy

    
      subroutine fits_header_main(fitsfile,bitpix,naxis1,naxis2,bscale,bzero,gain,exptime)
!         function:get the fits header information
!            input:fitsfile(fits file name)
!           output: bitpix:image 
!                   naxis1/2:image size
!                   bscale:
!                   bzero:
!                   gain:
!                   exptime:
     implicit real*8(a-h,o-z)
       parameter(maxheadrec=10)   
       integer*4 naxis1,naxis2,nmaxre,nbyf,npxf,i,bitpix
       parameter (nmaxre=15000,nbyf=2880,npxf=1440)
       character*200 fitsfile
       real*8 bscale,bzero,gain,exptime
       character *2880 rhead

      open(2,file=fitsfile,status='old',access='direct',recl=nbyf)    
      gain=0.0d0
      headerr=0
      do i=1,maxheadrec
          read(2,rec=i) rhead

        ic=index(rhead,'BITPIX  =')
        if(ic.ne.0) then
          read(rhead(ic+9:ic+29),'(i21)') bitpix
        endif           
          
          ic=index(rhead,'NAXIS1  =')
          if(ic.ne.0) then
              read(rhead(ic+9:ic+29),'(i21)') naxis1
          endif

          ic=index(rhead,'NAXIS2  =')
          if(ic.ne.0) then
              read(rhead(ic+9:ic+29),'(i21)') naxis2
          endif


          ic=index(rhead,'BSCALE  =')
          if(ic.ne.0) then
              read(rhead(ic+9:ic+29),*) bscale
          endif

          ic=index(rhead,'BZERO   =')
          if(ic.ne.0) then
              read(rhead(ic+9:ic+29),*) bzero
          endif

          ic=index(rhead,'GAIN    =')
          if(ic.ne.0) then
              read(rhead(ic+9:ic+29),'(f21.0)') gain
          endif

          ic=index(rhead,'EXPTIME =')
          if(ic.ne.0) then
              read(rhead(ic+9:ic+29),*) exptime
          endif

          ic=index(rhead,'END            ')
          if(ic.ne.0) then 
              goto 95
          endif
      enddo
95    close(2)
    end subroutine fits_header_main
    
      subroutine detect_gauss(fitsfile,bkgd_threshold,tele_label,sumi,starxx,staryy,snr,star_id,star_pix,overflag,nstar,bkgd,bkgdsigma)
!       function:detect stars and output the stars center(by gauss)
!          input:fitsfile:fits file name
!                bkgd_threshold:fitsdata-bkgd_threshold*bkgdsigma
!         output:sumi:the stars' flux
!                starxx:the x-position of star center
!                staryy:the x-position of star center
!                snr:the stars' snr
!                star_id:the id of the star( already sorted by brightness)
!                star_pix:The number of stars contained in star
!                overflag:Whether the stars overflow
!                nstar:Number of detected stars
!                bkgd:the image bkgd
!                bkgdsigma:the image bkgd
!Note:���Ŀ������=naxis1*naxis2/pi/6**2  �����������=pi*15**2, minpix=3
      
      implicit none
      character*200 fitsfile  
      character*100 tele_label
      integer*4 i,j,naxis1,naxis2,npx,bitpix,s,seq0,seq1,k,z,t,iloop
      integer*4 maxobj,minpix,maxpix,nli,nar,nstar,pos_nstar,p
      real*8 bzero,bscale,gain,bkgd,bkgdsigma,bkgd_threshold
      parameter(s=9000000)
      real*8 mdat(s),mdat01(s),abox(3000,3000),abox5(3000,3000),lty(4),x(5,1),a(400,5),b(400,1)
	  real*8 sumx(50000),sumy(50000),maxflux,sumi_real(50000),starx0(50000),stary0(50000)
      real*8 starx(50000),stary(50000),sumi(50000),pi,exptime,ss,h,starxx(50000),staryy(50000)
      integer*4 star_id(50000),star_pix(50000),overflag(50000)
	  real*8 pos_x(50000),pos_y(50000),pos_i(50000),snr(50000),snr01(50000)  
      integer*4 idpix(s),numpix(s),objpos(5000,5000)
      integer year,month,day,hh,mm,numobj
 	  integer nused,IFLAG1,ne
	  real*8 sigpar(5),rmag(500),sig0,sx,sy,si,b_noise
      
     pi=4*datan(1.0d0)
      bkgd=0
      bkgdsigma=0
!��ȡfitsͷ�ļ���������Ϣ����ȷ�����Ŀ����maxobj���������������maxpix�����maxflux          
!���Ŀ�����������ӳ���ȫ�ǰ뾶Ϊ6��������ex: 1024*1056/3.1415/6/6=9561��ʵ����û����ô���
!����������������ǣ�����뾶Ϊ15������������������һ��û����ô����������뾶��5����
      call fits_header_main(fitsfile,bitpix,naxis1,naxis2,bscale,bzero,gain,exptime)
      if(dabs(gain).le.1e-9) gain=1.0    
      maxflux=2**bitpix-1
      maxobj=int(naxis1*naxis2/pi/6**2)        
      minpix=3
      maxpix=int(pi*15**2)
      call fits_data(fitsfile,mdat,exptime,npx) 
      call cal_b(mdat,npx,bkgd,bkgdsigma)   !the whole image 
      
!��ͨ����ǰ�����ݽ��д���
!each pixes-��background��bkgd_threshold*bkgdsigma����bkgd_thresholdһ����5����ͼ���������sigma��С
      do i=1,npx
        mdat01(i)=mdat(i)-(bkgd+bkgd_threshold*bkgdsigma)
        if(mdat01(i).le.0.0d0) mdat01(i)=0.0d0
      enddo 
      abox=0.0
      do i=1,npx
        call seq2mat(i,naxis1,naxis2,nli,nar)    
        abox(nli,nar)=mdat01(i)
      enddo
!set the image edge pixel data =0
      do i=1,naxis2
        abox(i,1)=0.0d0
        abox(i,naxis1)=0.0d0
      enddo
      do j=1,naxis1
        abox(1,j)=0.0d0
        abox(naxis2,j)=0.0d0
      enddo
!!!!!��Ҫʱ�������	
!    write(*,'(a,3i5,i8)'),'maxobj/minpix/maxpix/npx:',maxobj,minpix,maxpix,npx
      
!��������ͨ��⣬��ÿ���ر��
!lty1=���ϣ�lty2=���ϣ�lty3=���ϣ�lty4=��
!Ŀ��ţ�numobj�����ر��idpix(numobj)=numobj,Ŀ�꺬������numpix(numobj)
!�������£���������������ؼ�⣬��ĳ����>0ʱ���ֱ��������жϣ�
!1)��4�ھ�<0���½�Ŀ��numobj+1����Ŀ��������numpix(numobj)+1�������ر��idpix(seq)=numobj
!2��lty1>0��1�ı�Ÿ���ǰ���أ�Ŀ��������+1
!3��lty2>0��2�ı�Ÿ���ǰ���أ�Ŀ��������+1
!4��lty3>0��3�ı�Ÿ���ǰ���أ�Ŀ��������+1
!5��lty4>0��4�ı�Ÿ���ǰ���أ�Ŀ��������+1
!maxobj�趨���Ŀ�������=ȫ����/�׾��������������������Ų������ʾ
!maxpixĿ�꺬���������=2���׾�Բ����,�������������ʾ
!����ͨ����һ�μ��ɣ�����Ҫ�ظ���ţ���������ظ����
!���ԭ����ɨ����һ��3���أ���������ͬ���ٿ����ͬ�У���������ͬ������4������
      numobj=0
      numpix=0
      idpix=0
      objpos=0
      do i=2,naxis2-1
        do j=2,naxis1-1
          lty(1)=abox(i-1,j+1)
          lty(2)=abox(i-1,j)
          lty(3)=abox(i-1,j-1)
          lty(4)=abox(i,j-1)
          if(abox(i,j).gt.1d-9) then
            call mat2seq(i,j,naxis1,naxis2,seq0)
!����ֻ��һ�����أ������0
            if(dabs(lty(1)+lty(2)+lty(3)+lty(4)).le.1d-9) then
              numobj=numobj+1
              idpix(seq0)=numobj !���ر��
              numpix(numobj)=numpix(numobj)+1  !Ŀ���ۼ�����
              objpos(idpix(seq0),numpix(numobj))=seq0
              if(numobj.gt.maxobj) then
                write(*,*)'Detect stars more than maxobj:',numobj,maxobj
                write(*,*)'Warning:Image quality problem!'
                goto 953
              endif
              if(numpix(numobj).gt.maxpix) then
                write(*,*)'Pixels of star more than maxpix:',numpix(numobj),maxpix
                write(*,*)'Warning:please check the radius value (maybe too small)!'
                goto 952
              endif                  
! ����
            elseif(lty(1).gt.1d-9) then
              call mat2seq(i-1,j+1,naxis1,naxis2,seq1)            
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
                write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
                write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                goto 952
              endif    
! ����
            elseif(lty(2).gt.1d-9) then
              call mat2seq(i-1,j,naxis1,naxis2,seq1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
                write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
                write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                goto 952
              endif     
! ����         
            elseif(lty(3).gt.1d-9) then
              call mat2seq(i-1,j-1,naxis1,naxis2,seq1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
                write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
                write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                goto 952
              endif     
! ���      
            elseif(lty(4).gt.1d-9) then
              call mat2seq(i,j-1,naxis1,naxis2,seq1)
              idpix(seq0)=idpix(seq1)
              numpix(idpix(seq0))=numpix(idpix(seq0))+1
              objpos(idpix(seq0),numpix(idpix(seq0)))=seq0
              if(numpix(idpix(seq0)).gt.maxpix) then
                write(*,*)'Pixels of star more than maxpix:',numpix(idpix(seq0)),maxpix
                write(*,*)'Warning:please check the radius1 value (maybe too small)!'
                goto 952
              endif
            endif
          endif
952     enddo  !952 'this image includes super star(>pi*r^2)'
      enddo
! ========================������ͨ����
953 continue   
!! ===================XZJ
!!compute the center of the object
! starx(k):x-value
! stary(k):y-value
! stari(k):flux-value��mdat01---has subtracted (bkgd+bkgd_threshold*bkgdsigma),not final value������
!ע�⣬�˴�ʹ�����ػҶ�ֵΪmdat01(i)=mdat(i)-(bkgd+bkgd_threshold*bkgdsigma)  
      k=1
      sumx=0.0
      sumy=0.0 
      sumi=0.0
      starx=0.0
      stary=0.0       
      do i=1,numobj
        if(numpix(i).ge.minpix.and.numpix(i).le.maxpix) then
          do j=1,numpix(i)
            call seq2mat(objpos(i,j),naxis1,naxis2,nli,nar)  
            if(mdat01(objpos(i,j)).gt.maxflux) overflag(k)=1
            sumx(k)=sumx(k)+nli*mdat01(objpos(i,j))**2
            sumy(k)=sumy(k)+nar*mdat01(objpos(i,j))**2
            sumi(k)=sumi(k)+mdat01(objpos(i,j))**2 
            sumi_real(k)=sumi_real(k)+mdat(objpos(i,j))!�����������
            !�����������ظ�˹
             ! a(j,1)=nli**2
             ! a(j,2)=nar**2
             ! a(j,3)=nli
             ! a(j,4)=nar
             ! a(j,5)=1.0
             !call mat2seq(nli,nar,naxis1,naxis2,seq0)
             !if(mdat(seq0).le.0) mdat(seq0)=1.0
             !b(j,1)=log(mdat(seq0))
          enddo
!maxim software:snr=sumi/bkgdsima/sqrt(pi*radius1**2)
!yuyong:tempsnrflux=sumvalue/dsqrt(sumvalue+numpix1*(bkgdsigma**2))
!tang:snr=mdat(i)-bkgd/bkgdsigma
          !  snr01(k)=sumi(k)/bkgdsigma/sqrt(numpix(i)) 
          !snr(k)=sumi(k)/sqrt(sumi(k)+numpix(i)*(bkgdsigma**2))
          snr(k)=sumi_real(k)/sqrt(sumi_real(k)+numpix(i)*(bkgdsigma**2)) !yuyong            
          starx(k)=sumx(k)/sumi(k)!��¼�����ؽ��
          stary(k)=sumy(k)/sumi(k)
          star_id(k)=i
          star_pix(k)=numpix(i)
      !=========����������Ϊ��ֵ��6Ϊ�뾶����ԲȦ���˹
          t=0
!ģ������
          !do p=1,size(abox(:,1))
          !    do z=1,size(abox(1,:))
          !      call RANDOM_NUMBER(b_noise)
          !      abox(p,z)=b_noise/10.0
          !    enddo
          !    enddo
          !
          !abox(8,10)=abox(8,10)+2.0
          !abox(9,10)=abox(9,10)+7.0
          !abox(10,10)=abox(10,10)+10.0
          !abox(11,10)=abox(11,10)+7.0
          !abox(12,10)=abox(12,10)+2.0
          !abox(10,8)=abox(10,8)+2.0
          !abox(10,9)=abox(10,9)+7.0
          !abox(10,10)=abox(10,10)+10.0
          !abox(10,11)=abox(10,11)+7.0
          !abox(10,12)=abox(10,12)+2.0
          !abox(9,9)=abox(9,9)+5
          !abox(11,9)=abox(11,9)+5
          !abox(9,11)=abox(9,11)+5
          !abox(11,11)=abox(11,11)+5
          !
          !starx(k)=10
          !stary(k)=10
          !sx=0
          !sy=0
          !si=0
100       if(starx(k).gt.6 .and.stary(k).gt.6 .and.starx(k).lt.(naxis1-6) .and.stary(k).lt.(naxis2-6))then
          do p=1,11
             nli=int(starx(k))-6+p
             do z=1,11
                nar=int(stary(k))-6+z
        !  !guass parameters  
                t=t+1
                a(t,1)=nli**2
                a(t,2)=nar**2
                a(t,3)=nli
	            a(t,4)=nar
                a(t,5)=1.0
                !call mat2seq(nli,nar,naxis1,naxis2,seq0)
                !if(mdat01(seq0).le.0) mdat01(seq0)=1.0
                !b(t,1)=log(mdat01(seq0)) !ʹ�ü�����ǰ
                if(abox(nli,nar).le.0) abox(nli,nar)=1.0
                b(t,1)=log(abox(nli,nar))!ʹ�ü�������
                !print*,a(t,1:5),b(t,1)
                
               !t sx=sx+nli*abox(nli,nar)!**2
                !t sy=sy+nar*abox(nli,nar)!**2
                !t si=si+abox(nli,nar)!**2 
            enddo
          enddo
          !t starx(k)=sx/si
          !t stary(k)=sy/si
!ʵ�飬�����xzj���Ϊ���Ļ�Ȧ���и�˹��λ�����Ȳ��ȶ���0.07-0.09��zhj0.02����
!����ʹ�ú�������ͬ�������ؼ��㣬����������Ȳ��ü�����֮ǰ������mdat0
          call least(a(1:t,:),b(1:t,:),x,t,5,t,sigpar,sig0,nused,IFLAG1)
          !============��Ȧ������λ�á�
          !call least(a(1:numpix(i),:),b(1:numpix(i),:),x,numpix(i),5,numpix(i),sigpar,sig0,nused,IFLAG1)
          !LEAST(A,Y,X,N1,K,NE,ERROR,SITA2,nrefused,IFLAG1)
            starxx(k)=-1*x(3,1)/(2*x(1,1))
            staryy(k)=-1*x(4,1)/(2*x(2,1))
            !ǰ��������ģ�0.001����
           ! if(abs(starx(k)-starxx(k)).gt.0.001 .or.abs(stary(k)-staryy(k)).gt.0.001)then
           ! starx(k)=starxx(k)
           ! stary(k)=staryy(k)
           ! iloop=iloop+1
           ! goto 100
           ! endif
           ! print*,'iloop',iloop
            write(*,'(a,i7)')'numpix(i)',numpix(i)
            write(*,'(a,2f10.3)')'xzj:x-y',starx(k),stary(k)
            write(*,'(a,2f10.3)')'guass:x-y',starxx(k),staryy(k)
            write(*,'(a,2f10.3)')'xzj-guass:x-y',starx(k)-starxx(k),stary(k)-staryy(k)
            write(*,'(a,f7.3,i5)')'sig--nused',sig0,nused
       !     pause
       endif
        k=k+1
       endif
      enddo 
      nstar=k-1
       call sorting_7terms(sumi,starxx,staryy,snr,star_id,star_pix,overflag,nstar)   
!note: guass
      !z=ax^2+by^2+cx+dy+f  ������X��abcdf��
      !A=��x^2,y^2,x,y,1�� ������������
      !B=z  ���ػҶ�ֵ
      !sigma_x^2=-a
      !sigma_y^2=-b
      !x0=-c/2a
      !y0=-d/2b
    end subroutine detect_gauss    
         
      subroutine fits_data(fitsfile,fitsdata,exptime,npx)
!         function:get the fitsdata to fitsdata(npx)
!            input:fitsfile(fits file name)
!           output: fitsdata(һά���飺naxis1*naxis2,��*bscale+bzero)
!                   exptime ������ͼ����ع�ʱ�䣩        
       implicit real*8(a-h,o-z)
       parameter(maxheadrec=10)
       character*200 fitsfile
       integer *4 nmaxre,nbyf,npxf,i
       parameter (nmaxre=15000,nbyf=2880,npxf=1440)
       integer *1 recdat(nbyf),recsw(nbyf)
       integer *2 recpx(npxf)
       equivalence (recsw,recpx)
       integer *4 npx,ntotpx,nrec,s
       integer *4 naxis1,naxis2,seqnum,nli,nar,bitpix
       character*2880 rhead
       real*8 bzero,bscale,gain,exptime
       parameter(s=9000000)
       real*8 fitsdata(s)
!  ��ȡͼ��ͷ�ļ������Ϣ
       open(33,file=fitsfile,status='old',access='direct',recl=nbyf)    
       do i=1,maxheadrec
         read(33,rec=i) rhead

         ic=index(rhead,'BITPIX  =')
         if(ic.ne.0) then
           read(rhead(ic+9:ic+29),'(i21)') bitpix
         endif    

         ic=index(rhead,'NAXIS1  =')
         if(ic.ne.0) then
           read(rhead(ic+9:ic+29),'(i21)') naxis1
         endif

         ic=index(rhead,'NAXIS2  =')
         if(ic.ne.0) then
           read(rhead(ic+9:ic+29),'(i21)') naxis2
         endif

         ic=index(rhead,'BSCALE  =')
         if(ic.ne.0) then
           read(rhead(ic+9:ic+29),*) bscale
         endif

         ic=index(rhead,'BZERO   =')
         if(ic.ne.0) then
           read(rhead(ic+9:ic+29),*) bzero
         endif

         ic=index(rhead,'EXPTIME =')
         if(ic.ne.0) then
           read(rhead(ic+9:ic+29),*) exptime
         endif 

         ic=index(rhead,'END            ')
         if(ic.ne.0) then 
           goto 93
         endif
       enddo
93    nrec=i+1
      npx =1
      ntotpx = naxis1*naxis2

!  ���㹲�ж��ٸ����ݿ�nbyf=2880,npxf=1440
      numdatarec=int(ntotpx*bitpix/8/nbyf)
      do i=1,numdatarec
        read(33,rec=nrec) recdat
        do k = 2,2880,2             
          recsw(k-1)=recdat(k)
          recsw(k)=recdat(k-1)
        enddo      

        do j = 1,npxf
          fitsdata(npx) = recpx(j)*bscale+bzero       
          if(fitsdata(npx).lt.0) fitsdata(npx)=fitsdata(npx)+65536  
          npx = npx + 1
        enddo
		nrec=nrec+1
      enddo
!   ��ȡʣ�����أ�2019.10.11����         
      read(33,rec=nrec) recdat
      do k = 2,2*(ntotpx-1440*(numdatarec)),2							
          recsw(k-1)=recdat(k)
          recsw(k)=recdat(k-1)
      enddo      

      do j = 1,ntotpx-1440*(numdatarec)
          fitsdata(npx) = recpx(j)*BScale+bzero     
          if(fitsdata(npx).lt.0) fitsdata(npx)=fitsdata(npx)+65536  
          npx = npx + 1
      enddo
      npx=npx-1  
	  close(33)
!һά�����Ϊ��ά���󣨿�ʵ�������ά����Ŀǰδ�ã�
!      fitsdata=0.0d0
!      do i=1,npx
!        call seq2mat(i,naxis1,naxis2,nli,nar)    
!        fitsdata(nli,nar)=mdat(i)
!      enddo
      end subroutine fits_data

     subroutine cal_b(mdat,npx,bkgd,bkgdsigma)
!       function:calculate the background and bkgdsigma of fits
!          input:mdat��fitsdata
!                npx:the size of the mdat
!         output:bkgd:the average of fitsdata
!                bkgdsigma:Background relief
!ע�⣺����ֵΪ�۳�>2.6sigma������sigma���첻���Ϊ������sigmaֵ
!      �����ͱ�����ͼ���������Ӱ�죬��Ϊ����ֵ�ͱ������
        implicit real*8(a-h,o-z)
        integer*4 npx
        real*8 mdat(npx),bkgd,bkgdsigma
        parameter(sfa=2.6d0)
        integer k,nstep

        nstep=1
        k=1
        sumvalue=0.0d0
        do i=1,npx,nstep
          sumvalue=sumvalue+mdat(i)
          k=k+1
        enddo
        avervalue=sumvalue/real(k-1)
        sumvalue2=0.0d0
        do i=1,npx,nstep
          sumvalue2=sumvalue2+(mdat(i)-avervalue)**2
        enddo
        sigma=dsqrt(sumvalue2/(k-1-1))
!  ��������ֵ-ƽ��ֵ����С��2.6sigma���뱳��ֵ�������򲻼���
10      sumvalue=0.0d0
        k=1
        do i=1,npx,nstep
          temp=mdat(i)-avervalue 
          if(dabs(temp).le.sfa*sigma) then
            sumvalue=sumvalue+mdat(i)
            k=k+1
          endif
        enddo
        avervalue2=sumvalue/(k-1)
        sumvalue2=0.0d0
        do i=1,npx,nstep
          temp=mdat(i)-avervalue 
          if(dabs(temp).le.sfa*sigma) then
            sumvalue2=sumvalue2+(mdat(i)-avervalue2)**2
          endif
        enddo
        sigma2=dsqrt(sumvalue2/(k-1-1))

        if(sigma2.le.1d-6) then
 !        print*,'ͼ���Ϊƽ̹(bkgd,sigma):',avervalue,sigma
          goto 20
        endif

        if(dabs(sigma2-sigma).ge.0.01*sigma2) then
!         print*,'ͼ��ƽ̹(bkgd,sigma):',avervalue,sigma
          avervalue=avervalue2
          sigma=sigma2
          goto 10
        endif
! ��Ϊ������߲��첻�����һ�μ���ɴ�������ͼƽ����ƽ̹�ȡ�
! һ��������5������
20      bkgd=avervalue
        bkgdsigma=sigma
        end subroutine cal_b



      subroutine sorting_7terms(x01,x02,x03,x04,x05,x06,x07,n)
!       function:according to he x01,reorder the sequencex01,x02,x03,x04,x05,x06,x07
!          input:x01,x02,x03,x04,x05,x06,x07:old sequence  
!                N:the size of sequence
!         output:x01,x02,x03,x04,x05,x06,x07:new sequence
      implicit real*8(a-h,o-z)
      real*8 x01(n),x02(n),x03(n),x04(n),tempx
      integer*4 x05(n),x06(n),x07(n),tempi

      do i=1,n-1
        k=i
        do j=i+1,n
          if(x01(j).ge.x01(k)) then
            k=j
          endif
        enddo
        tempx=x01(i)
        x01(i)=x01(k)
        x01(k)=tempx

        tempx=x02(i)
        x02(i)=x02(k)
        x02(k)=tempx

        tempx=x03(i)
        x03(i)=x03(k)
        x03(k)=tempx

        tempx=x04(i)
        x04(i)=x04(k)
        x04(k)=tempx

        tempi=x05(i)
        x05(i)=x05(k)
        x05(k)=tempi
    
        tempi=x06(i)
        x06(i)=x06(k)
        x06(k)=tempi   
            
        tempi=x07(i)
        x07(i)=x07(k)
        x07(k)=tempi 
      enddo
      end subroutine sorting_7terms

      subroutine seq2mat(seqnum,naxis1,naxis2,nli,nar) 
!       function:Calculate the column Numbers in the matrix
!          input:seqnum:sequence number
!                naxis1,naxis2:The size of the matrix
!         output:nli,nar��column Numbers in the matrix(naxis1,naxis2)  
! nli--n_line;nar-n_arrager
        implicit real*8(a-h,o-z)
        integer seqnum,naxis1,naxis2,nli,nar,dd
    
        if(mod(seqnum,naxis1).eq.0) then
          nli=int(seqnum/naxis1)
        else
          nli=int(seqnum/naxis1)+1
        endif
        nar=seqnum-(nli-1)*naxis1
        end subroutine seq2mat

      subroutine mat2seq(nli,nar,naxis1,naxis2,seqnum)
!       function:Calculate the sequence Numbers in the matrix
!          input:nli,nar��column Numbers in the matrix(naxis1,naxis2)
!                naxis1,naxis2:The size of the matrix
!         output:seqnum:sequence number
! nli--n_line;nar-n_arrager
        implicit real*8(a-h,o-z)
        integer seqnum,naxis1,naxis2,nli,nar
        
        seqnum=naxis1*(nli-1)+nar
    end subroutine mat2seq
    
	SUBROUTINE LEAST(A,Y,X,N1,K,NE,ERROR,SITA2,nrefused,IFLAG1)
!    function:calculate the matrix X by least square
!       input:A:�۲����Ͼ���
!           Y:�Ⱥ��Ҳ����
!           X:�������
!           N1:�۲���������������������
!           K:�����������
!           NE: �۲���������������������   
!    output:X:����Ƭ������ֵ
!           ERROR:���������������
!           SITA2:�����������sqrt(sum(res**2)/nep-k)��??�������Ƿ�Ӧ��ֱ����nep-1
!           nrefused��������ղ���ʹ�õĹ۲���������
!           rmag�������ˣ�ʵ��δ�ã�Ӧ������Ҫ�������ȹ���в��
!           p��ʹ��1��δʹ��0������ʾ����N1
!           vksi��ksi����в����N1
!           vyit:yit����в����N1
!           IFLAG1�����������������30�Σ����Ϊ1����Ϊ��Ȼ�ǹ������ǵ�����������������
	 implicit real*8(a-h,o-z)
 	 REAL*8 A(N1,K),X(K,1),Y(N1,1),AT(K,NE),ATA(K,K),ATAAT(K,NE)
	 REAL*8 V(NE),V2(NE),SITA2,ERROR(K),P(NE),rmag(ne/2),vksi(ne/2),vyit(ne/2)
	 INTEGER IFLAG1
!a(400,5),b(400,1),x(5,1)

      parameter(SFA=2.6,SFB=0.01)
 
      IFLAG1=0
      pi=4*datan(1.0d0)
	  q=pi/180.0d0
      SITA2=0.0D0
      SITA2P=0.0D0
	  P=1.0
	  AT=0.0
	  NBIG=0
	  SITA20=0.0D0
!write(*,'(a,5f11.2)')'inleast a=',a(1:n1,1:5)
!pause
	  ILOOP=0
666	  NEP=NE-NBIG
	  SITA2P=SITA2
	  ILOOP=ILOOP+1
	  SITA2=0.0
	  DO I=1,K
 		DO J=1,NE
			AT(I,J)=A(J,I)
            !write(*,'(a10,2i6,2f10.1)') 'in least',J,i,A(j,i)
	    ENDDO
      ENDDO
      
	  DO I=1,NE
		DO J=1,K
			AT(J,I)=AT(J,I)*P(I)
		ENDDO
      ENDDO
!ATAX=ATY(��С���˽ⳬ��������)
  	  CALL MUL(AT,A,ATA,K,NE,NE,K)
	  CALL MATINV(K,K,ATA)
      CALL MUL(ATA,AT,ATAAT,K,K,K,NE)
      CALL MUL(ATAAT,Y,X,K,NE,NE,1)
   !   print*,'x',x
    !xx=-1*x(3,1)/(2*x(1,1))!������λ��
    !yy=-1*x(4,1)/(2*x(2,1))
    !print*,'xx-yy',xx,yy
    !pause     
!�˴�x�����������K��
      v=0.0
 	  DO 10 I=1,NE
		DO 11 J=1,K
			V(I)=V(I)+A(I,J)*X(J,1)
11		    CONTINUE
		    V(I)=V(I)-Y(I,1)
		    V2(I)=V(I)**2
		    SITA2=SITA2+V2(I)*P(I)
10    CONTINUE
!ֻҪʣ���������������ڴ���������ͼ�����std����������
	        IF(NEP.NE.K) THEN 
               ! SITA2=DSQRT(SITA2/(NEP-1))!zz
		       SITA2=DSQRT(SITA2/(NEP-K))
	        ELSE
		       SITA2=0.0
	        ENDIF
            IF(DABS(SITA2-SITA20).GE.SFB*SITA20.AND.ILOOP.LE.30) THEN
		       SITA20=SITA2
		       NBIG=0
		       P=1.0
		       DO I=2,NE,2
		         IF(DABS(V(I)).GE.SFA*SITA2) THEN !���в������ɾ������
				   P(I)=0.0
				   NBIG=NBIG+1
		         ENDIF
		       ENDDO
		       GOTO 666
            ELSEIF(ILOOP.GT.30) THEN
	           IFLAG1=1
            ENDIF
!            print*,'iloop=',iloop
!	        do i=2,ne,2
!		      vksi(i/2)=v(i-1)/q*3600.0
!		      vyit(i/2)=v(i)/q*3600.0
!	        enddo   
	        nrefused=ne-nbig
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
     
      
      subroutine least_square(a,b,x,m,n)
      !fuction:Solving parameter sequence by least squares:AX=B-->ATAX=ATB
      !m ������������У�n��������
 !     use mkl
      integer m,n,i,L
      real*8 a(m,n),b(n,1),x(m,1)
      real*8 c(n,m),d(m,m),e(1,m)      
      !real*8 f(3,3),g(1,3),h(3,3)
      !data f/1,1,1,2,2,2,3,3,3/
      !data g/2,2,2/
      
      if (m<n)then
          c=transpose(a)
          d=matmul(a,c)
          e=matmul(c,b)
          call AGAUS(d,e,m,x)  
        !  print*,'e',size(e)
         ! pause
 !         print*,'agaus finish'
        !  pause
      else if(m==n)then
          call AGAUS(a,b,m,x)
      else
          write(*,'(a,i5,i5)')'In least_square, m<n, No roots! m,n:',m,n
      endif
    end subroutine least_square

    !������Է�����
    SUBROUTINE AGAUS(A,B,N,X) 
	DIMENSION A(N,N),X(N),B(N),JS(N) 
	DOUBLE PRECISION A,B,X,T 
	L=1 
	DO 50 K=1,N-1 
	  D=0.0 
	  DO 210 I=K,N 
	  DO 210 J=K,N 
	    IF (ABS(A(I,J)).GT.D) THEN 
	      D=ABS(A(I,J)) 
	      JS(K)=J 
	      IS=I 
	    END IF 
210   CONTINUE 
    !  print*,D
     ! pause
	  IF (D+1.0.EQ.1.0) THEN 
	    L=0 
	  ELSE 
	    IF (JS(K).NE.K) THEN 
	      DO 220 I=1,N 
	        T=A(I,K) 
	        A(I,K)=A(I,JS(K)) 
	        A(I,JS(K))=T 
220	      CONTINUE 
	    END IF 
	    IF (IS.NE.K) THEN 
	      DO 230 J=K,N 
	        T=A(K,J) 
	        A(K,J)=A(IS,J) 
	        A(IS,J)=T 
230	      CONTINUE 
	      T=B(K) 
	      B(K)=B(IS) 
	      B(IS)=T 
	    END IF 
	  END IF 
	  IF (L.EQ.0) THEN 
	   ! WRITE(*,100) 
	    RETURN 
	  END IF 
	  DO 10 J=K+1,N 
	    A(K,J)=A(K,J)/A(K,K) 
10	  CONTINUE 
	  B(K)=B(K)/A(K,K) 
	  DO 30 I=K+1,N 
	    DO 20 J=K+1,N 
	      A(I,J)=A(I,J)-A(I,K)*A(K,J) 
20	    CONTINUE 
	    B(I)=B(I)-A(I,K)*B(K) 
30	  CONTINUE 
50  CONTINUE 
    !print*,ABS(A(N,N))
    !pause
	IF (ABS(A(N,N))+1.0.EQ.1.0) THEN 
	  L=0 
      !WRITE(*,100)
	  !write(*,*) 
	  RETURN 
	END IF 
	X(N)=B(N)/A(N,N) 
	DO 70 I=N-1,1,-1 
	  T=0.0 
	  DO 60 J=I+1,N 
	    T=T+A(I,J)*X(J) 
60	  CONTINUE 
	  X(I)=B(I)-T 
70	CONTINUE 
100 FORMAT(1X,' guass FAIL ') 
	JS(N)=N 
	DO 150 K=N,1,-1 
	  IF (JS(K).NE.K) THEN 
	    T=X(K) 
	    X(K)=X(JS(K)) 
	    X(JS(K))=T 
	  END IF 
150	CONTINUE 
	RETURN 
    END     
    