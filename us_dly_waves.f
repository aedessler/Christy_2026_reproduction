rm fort.*
rm source.f
ln -s ushcn-v2-stations.txt fort.10
ln -s tmistndly_dec_mar_25a.txt fort.21
ln -s tmxstndly_may_sep_25.txt fort.22
ln -s tmistndly_dec_mar_adj.txt fort.26
ln -s tmxstndly_may_sep_adj.txt fort.27
ln -s ushcn_stn_order.txt fort.30
ln -s usa48_halfdeg.txt fort.35
ln -s usreg_half.txt fort.36
cat >source.f<<'eof'
  
      program us
      parameter (nyb = 1899, nye = 2025, nx = 116, ny = 50)
      integer ival(nyb:nye,0:12,31),miss(10)
      integer nyval(1218,nyb:nye)
      dimension xsort(20000),xsval(49,nyb:nye),xmed(1218,nyb:nye)
      dimension xeval(1218,nyb:nye),xseval(49,nyb:nye)
      dimension xval(200),xlat(1218),xlon(1218),ddd(0:10)
      dimension dist(1218,nx,ny),xg(nyb:nye,nx,ny)
      dimension xsreg(0:9,nyb:nye),dtot(nyb:nye)
      integer ns(500000),nsval(49,nyb:nye),ntot(49)
      integer nst(1218),ntval(200),lastdy(nyb:nye,0:12)
      integer monthe(0:12),jval(nyb:nye,200),kval(200),ktval(200)
      integer nval(nyb:nye,1218),nus48(nx,ny),nreg(9,nx,ny)
      integer nusreg(nx,ny)
      character*2 ast(49)
      character*5 astr(0:9)
      data monthe/31,31,28,31,30,31,30,31,31,30,31,30,31/
      data ast/'AL','AZ','AR','CA','CO','CN','DE','FL',
     + 'GA','ID','IL','IN','IA','KS','KE','LA','ME','MD',
     + 'MA','MI','MN','MS','MO','MT','NB','NV','NH','NJ',
     + 'NM','NY','NC','ND','OH','OK','OR','PE','RI','SC',
     + 'SD','TN','TX','UT','VE','VA','WA','WV','WI','WY','US'/
      data astr/'CONUS','PacSW','PacNW','4Crns','NPlns',
     +          'SPlns','UMidW','OhVly','SoEst','NoEst'/ 

      do iyr = nyb,nye
        do imo = 0,12
          lastdy(iyr,imo) = monthe(imo)
        end do
      end do
      do iyr = 1904,nye,4   ! leap year
        lastdy(iyr,2) = 29
      end do
      do iy = ny,1,-1
c       read (35,*),(nus48(ix,iy),ix=1,nx)
      end do

      do iy = ny,1,-1
        read (36,*) (nusreg(ix,iy),ix=1,nx)
c       print *,iy
        print 600,iy,(nusreg(ix,iy),ix=1,nx)
      end do 
  600 format (i4,1x,130i1)

      do ix = 1,nx
        print 404,(nusreg(ix,iy),iy=1,ny)
      end do

      print *,' Enter 0 for no adj, 1 for adjusted '
      read *,iadj
 
 
      do istn = 1,1218
        read (30,*),inum,istnnum
        ns(istn) = inum
        nst(istn) = istnnum
      end do 
      print *,' Enter 1 for Min, 2 for Max'
      read *,ivar
      if(ivar .eq. 1) then
        iunit = 21 + iadj*5
        nmon1 = 0
        nmon = 4
        jtot = 122
      else
        iunit = 22 + iadj*5
        nmon1 = 5
        nmon = 5
        jtot = 153
      end if
c
c     read perentile
c
      print *,' Enter percentile as integer '
      read *,nperc
      print *,' Enter number of days in run'
      read *,nrun
c**********
  400 format(4i7,2f8.2)
  401 format(2i7,5f8.2)
  402 format(52f5.1)
  404 format(52i2)
      do istn = 1,1218
        read (10,*) id,xlt,xln
        xlat(ns(istn)) = xlt
        xlon(ns(istn)) = xln
c       print 400,id,istn,ns(istn),id,xlt,xln
      end do
c
c     calcualte station distance from each grid point
c     half degree
c
      r = 6335.44
      pi180 = 3.1415926535 / 180.
      do iy = 1,ny
        xlt = iy*0.5 + 24.75
        do ix = 1,nx
          xln = ix*0.5 - 124.75
          do istn = 1,1218
               if( nusreg(ix,iy) .gt. 0) then
                  clat1 = cos(xlt*pi180)
                  clat2 = cos(xlat(ns(istn))*pi180)
                  slat1 = sin(xlt*pi180)
                  slat2 = sin(xlat(ns(istn))*pi180)
                  clon1 = cos(xln*pi180)
                  clon2 = cos(xlon(ns(istn))*pi180)
                  slon1 = sin(xln*pi180)
                  slon2 = sin(xlon(ns(istn))*pi180)
                  slon21 = sin((xlon(ns(istn))-xln)*pi180)
                  clon21 = cos((xlon(ns(istn))-xln)*pi180)
c                 print 130,clat1,clat2,slat1,slat2,clon1,
c    + clon2,slon1,slon2,slon21,clon21
                  x1 = clat1*slon21*clat1*slon21
                  x2 = clat2*slat1 - slat2*clat1*clon21
                  x2 = x2 * x2
                  x3 = slat2 * slat1
                  x4 = clat2 * clat1 * clon21
                  x5 = sqrt(x1 + x2) / (x3 + x4)
                  dist(ns(istn),ix,iy) = atan(x5) * r
c                 if( ix .eq. 90) then
c                   if( iy .eq. 35 ) then
c                   print 401,ns(istn),istn,dist(ns(istn),ix,iy),
c    +     xlat(ns(istn)),xlon(ns(istn)),xlt,xln
c                   end if
c                 end if
           
                end if

          end do  ! istn
        
        end do  ! ix

      end do  ! iy

c       print *,id,xlat,xlon
c           inum = inum + 1
c           nd(id) = inum
c           nid(inum) = id
c           asts(inum) = ast
c           astns(inum) = astn
c           xlt(inum) = xlat
c           xln(inum) = xlon
c           print 129,inum,nd(id),nid(inum),asts(inum)
c           do ix = 1,nx
c             do iy = 1,ny
c               print *,' ix,iy,nusreg ',ix,iy,nusreg(ix,iy)
c               if( nus48(ix,iy) .gt. 0) then
c                 clat1 = cos(xt(iy)*pi180)
c                 clat2 = cos(xlat*pi180)
c                 slat1 = sin(xt(iy)*pi180)
c                 slat2 = sin(xlat*pi180)
c                 clon1 = cos(xl(ix)*pi180)
c                 clon2 = cos(xlon*pi180)
c                 slon1 = sin(xl(iy)*pi180)
c                 slon2 = sin(xlon*pi180)
c                 slon21 = sin((xlon-xl(ix))*pi180)
c                 clon21 = cos((xlon-xl(ix))*pi180)
c                 print 130,clat1,clat2,slat1,slat2,clon1,
c    + clon2,slon1,slon2,slon21,clon21
c                 x1 = clat1*slon21*clat1*slon21
c                 x2 = clat2*slat1 - slat2*clat1*clon21
c                 x2 = x2 * x2
c                 x3 = slat2 * slat1
c                 x4 = clat2 * clat1 * clon21
c                 x5 = sqrt(x1 + x2) / (x3 + x4)
c                 xd(inum,ix,iy) = atan(x5) * r
c                 if( inum .eq. 50) then
c                   print *,ix,iy,inum,xd(inum,ix,iy)
c                 end if
c               end if
c             end do
c           end do
c           xel(inum) = xele
c           num = inum
c     end do

c***************
c
c     adjustment manual
c
      print *,' Enter 1 for adjustment'
      read *,nadj
      if( nadj .eq. 1) then
        print *,' Enter value of adjustment year1 year2, stninc'
        read *,madj,myear1,myear2,madjinc
      end if
      xlonw = -130
      xlone = -50
      xlatn = 50
      xlats = 20

      do istn = 1,1218  ! 1218
      
      do myr = nyb,nye

        do imo = nmon1,nmon1+nmon-1   !
        
        read (iunit,*) nstn,iyr,mmo,
     +    (ival(myr,mmo,id),id=1,31)

        if( nadj .eq. 1) then
        if( xlon(ns(istn)).gt.xlonw.and.xlon(ns(istn)).lt.xlone) then
        if( xlat(ns(istn)).gt.xlats.and.xlat(ns(istn)).lt.xlatn) then
c       print *,ns(istn),xlon(ns(istn)),xlat(ns(istn))
        if( myr .ge. myear1 .and. myr .le. myear2) then
          if( mod(istn,madjinc) .eq. 0) then
          do id = 1,31
            if( ival(myr,mmo,id) .gt. -90) then
              ival(myr,mmo,id) = ival(myr,mmo,id) + madj
            end if
          end do
          end if
        end if
        end if
        end if
        end if
c       if( myr .eq. 2021) then
c         do id = 1,31
c           ival(myr,mmo,id) = -999
c         end do
c       end if
c
c     December will be month zero of year
c
        if( ivar .eq. 1) then
        if(mmo .eq. 12) then
          do id = 1,31
            ival(iyr+1,0,id) = ival(myr,12,id)
          end do
        end if
        end if

        end do   ! imo
        nst(istn) = nstn

      end do   ! myr

c     print *,'  '
c
c     all data for station istn is in ival
c
c     calculate threshhold for each day of season
c     place data in single string for the season
c

      do iyr = nyb,nye
        do jdy = 1,jtot
          jval(iyr,jdy) = -999
        end do

        icnt = 0
        do imo = nmon1,nmon1+nmon-1
          do idy = 1,lastdy(iyr,imo)
            icnt = icnt + 1
            jval(iyr,icnt) = ival(iyr,imo,idy)
          end do
        end do
c       print 202,iyr,(jval(iyr,ic),ic=1,icnt)
      end do
      
        jtot1 = jtot
        if( lastdy(iyr,2) .eq. 28) jtot1 = jtot - 1
        do idy = 1,jtot1

          icnt = 0
          jd1 = idy-3
          jd2 = idy+3
          if( jd1 .lt. 1) jd1 = 1
          if( jd2 .gt. jtot1) jd2 = jtot1

          do myr = nyb,nye

            do jday = jd1,jd2
            if( jval(myr,jday) .gt. -990) then
              icnt = icnt + 1
              xsort(icnt) = jval(myr,jday)
            end if
            end do
          end do  ! myr

c           if( istn .eq. 1) then
c           print 202,idy,(nint(xsort(ic)),ic=1,icnt)
c           end if

          if( icnt .gt. 200) then
            call sort(icnt,xsort(1:icnt))
c           print 202,idy,icnt,(nint(xsort(ix)),ix=1,icnt)
            ntval(idy) = xsort(int(icnt*nperc*.01 + 1) )
c           if( istn .lt. 10) then
c             print *,nst(istn),idy,ntval(idy),icnt
c           end if
            
          else
            ntval(idy) = -999
          end if
        end do  ! idy
c
c     count runs per year
c
     
      do iyr = nyb,nye
        iycnt = 0
        do idy = 1,jtot
          if( jval(iyr,idy) .gt. -900) then
            iycnt = iycnt + 1
          end if
        end do
        do idy = 1,jtot
          kval(idy) = 0
          ktval(idy) = 0
          if( ntval(idy) .gt. -900) then

          if( ivar .eq. 1) then    ! tmin
            if( jval(iyr,idy) .le. ntval(idy)) then
              kval(idy) = 1
            end if
          else   ! tmax
            if( jval(iyr,idy) .ge. ntval(idy)) then
             kval(idy) = 1
c            if( iyr .eq. 1936 .and. nst(istn) .eq. 32356) then
c            print 620,nst(istn),iyr,idy,jval(iyr,idy),ntval(idy)
c            end if
            end if
          end if
         
          end if

        end do ! idy

        icnt = 0    ! days in run
        jcnt = 0    ! index of runs
        do idy = 2,jtot
          if( kval(idy) .gt. 0) then
            if( kval(idy-1) .eq. 0) then  ! start
              jcnt = jcnt + 1  ! start
              icnt =  1
            end if
          end if
          if( kval(idy) .gt. 0) then
            if( kval(idy-1) .gt. 0) then
              icnt = icnt + 1
              if( idy-1 .eq.1) jcnt = jcnt + 1
            end if
          end if
          if (kval(idy-1) .gt. 0) then  ! end
            if( kval(idy) .eq. 0) then
              ktval(jcnt) = icnt 
              if( iyr .gt. 1984) then
                if( icnt .ge. nrun) then
                  do jdy = idy  - icnt, idy
              print 620,nst(istn),iyr,jdy,icnt,
     +               jval(iyr,jdy),ntval(jdy)
                  end do
                end if
              end if
            end if
          end if

        end do
        nval(iyr,istn) = 0
        do ikt = 1,jcnt
c         print *,nst(istn),iyr,ikt,ktval(ikt)
          if( ktval(ikt) .ge. nrun) then
            nval(iyr,istn) = nval(iyr,istn) + ktval(ikt)
          end if
        end do
        if( iycnt .lt. nint(jtot *0.7) ) then
            jcnt = 1
            nval(iyr,istn) = -99
        end if
c       print *,nst(istn),iyr,nval(iyr,istn),iycnt
      end do ! iyr


      end do ! istn
  620 format(i7,i5,4i4)

      do istate = 1,48
        do iyr = nyb,nye
        xden = 0.
        xsval(istate,iyr) = 0 
        do istn = 1,1218
          jstate = nst(istn)/10000
          if( jstate .eq. istate) then
          if( nval(iyr,istn) .gt. -90) then
            xsval(istate,iyr) = xsval(istate,iyr) + nval(iyr,istn)
            xden = xden + 1
          end if
          end if
        end do
        if( xden .gt. 0) then
          xsval(istate,iyr) = xsval(istate,iyr)/xden
        else
          xsval(istate,iyr) = -99
        end if
        end do  ! iyr
      end do

c     print *,'  '


c    
c     area weighting
c
      rad = 115

c
c     complete grid
c
      ddtot = 0. 
      do iy = 1,ny
        xlt = iy*.5 + 24.75
        do ix = 1,nx
          if( nusreg(ix,iy) .gt. 0) then
          ddtot = ddtot +  cos(xlt*pi180)
          end if
        end do
      end do 
c
c     eliminate boxes of stations
c
      do istn = 1,1218  
        xxlat = xlat(istn)
        xxlon = xlon(istn)
        if( xxlat .gt. 40 .and. xxlat .lt. 42.5) then  ! Nebraska
c       if( xxlon .lt. -97.5 .and. xxlon .gt. -100.5) then
        if( xxlon .lt. -96.0 .and. xxlon .gt. -102.0) then
          print *,' Irrigation',nst(istn)
          do iyr = nyb,nye
c           nval(iyr,ns(istn)) = -99
          end do
        end if
        end if
        if( xxlat .gt. 34 .and. xxlat .lt. 36.5) then  ! Arkansas
c       if( xxlon .lt. -90.5 .and. xxlon .gt. -91.5) then
        if( xxlon .lt. -90.0 .and. xxlon .gt. -92.0) then
          print *,' Irrigation',nst(istn)
          do iyr = nyb,nye
c           nval(iyr,ns(istn)) = -99
          end do
        end if
        end if
        if( xxlat .gt. 36 .and. xxlat .lt. 38) then  ! California 
c       if( xxlon .lt. -120 .and. xxlon .gt. -121.5) then
        if( xxlon .lt. -119.25 .and. xxlon .gt. -122.25) then
          print *,' Irrigation',nst(istn)
          do iyr = nyb,nye
c           nval(iyr,ns(istn)) = -99
          end do
        end if
        end if
      end do

      do iyr = nyb,nye

      do iy = 1,ny
       do ix = 1,nx

         if( nusreg(ix,iy) .gt. 0) then

         xg(iyr,ix,iy) = 0.
         dd = 0.
         dstn = 0.
         do istn = 1,1218   ! doing this in USHCN order
           if( dist(ns(istn),ix,iy) .lt. rad ) then ! assigned to COOP
             if( nval(iyr,ns(istn)) .gt. -90) then
               ddx = rad/dist(ns(istn),ix,iy)
               xg(iyr,ix,iy) =
     +         xg(iyr,ix,iy) + nval(iyr,ns(istn))*ddx*ddx
               dd = dd + ddx*ddx
               dstn = dstn + 1
             end if
           end if
         end do

         if( dstn .gt. 1.1) then
           xg(iyr,ix,iy) = xg(iyr,ix,iy) / dd
         elseif( dstn .gt. 0.1 .and. sqrt(dd) .gt. 1.5) then
           xg(iyr,ix,iy) = xg(iyr,ix,iy) / dd
         else
           xg(iyr,ix,iy) = -999
         end if

         else   ! not a grid
           xg(iyr,ix,iy) = -999
         end if

       end do
      end do

      if( iyr .eq. 1984) then
        do ix = 1,nx
c         print 402,(xg(iyr,ix,iy),iy=1,ny)
        end do
        
      end if
c
c     area average
c
      do ireg = 0,9
       xsreg(ireg,iyr) = 0
       ddd(ireg) = 0.
      end do

      do iy = 1,ny
        xlt = iy*0.5 + 24.75
        do ix = 1,nx
          if( xg(iyr,ix,iy) .gt. -90) then
c
c    us48
c
          xsreg(0,iyr) = xsreg(0,iyr) + xg(iyr,ix,iy)*cos(xlt*pi180)
          ddd(0) = ddd(0) + cos(xlt*pi180)
c
c    regions
c
          jreg = nusreg(ix,iy)   ! jreg ranges 1 to 9
          xsreg(jreg,iyr) = xsreg(jreg,iyr) +
     +                             xg(iyr,ix,iy) * cos(xlt*pi180)
          ddd(jreg) = ddd(jreg) + cos(xlt*pi180)
          
          end if
        end do
      end do
  
      do ireg = 0,9
        if( ddd(ireg) .gt. 0.1) then
          xsreg(ireg,iyr) = xsreg(ireg,iyr) / ddd(ireg)
        else
          xsreg(ireg,iyr)= -99
        end if
      end do
        xsval(49,iyr) = xsreg(0,iyr)
c       dtot(iyr) = ddd(0) / ddtot
        print *,iyr,ddd(0)/ddtot


      end do  ! iyr

      print *,'  '
      print *,' Fraction of Coverage '
      do iyr = nyb,nye
c       print 211,iyr,dtot(iyr)
      end do
      print *,'  '


      do jstn = 1,1218,32
        kstn = jstn + 31
        if( kstn .gt. 1218) kstn = 1218
        print 300,(nst(istn),istn=jstn,kstn)
        do iyr = nyb,nye
          print 301,iyr,(nval(iyr,istn),istn = jstn,kstn)
        end do
        print *,'   '
      end do



c    
      print *,'  '
      print *,' State average ',nperc,' percentile',nrun,' days'
      print 200,(ast(is),is=1,24)
      do iyr = nyb,nye
        print 101,iyr,(xsval(is,iyr),is=1,24)
      end do
      print 200,(ast(is),is=1,24)
      print *,'  '
      print 201,(ast(is),is=25,49)
      do iyr = nyb,nye
        print 102,(xsval(is,iyr),is=25,49)
      end do
      print 201,(ast(is),is=25,49)

c    
      print *,'  '
      print *,' Regional averages percentile',nrun,' days'
      print 210,(astr(is),is=1,9),astr(0)
      do iyr = nyb,nye
        print 211,iyr,(xsreg(is,iyr),is=1,9),xsreg(0,iyr)
      end do
     
      print *,nadj,madj,myear

      stop

      print *,'  '
      print *,' State average magnitude of '
      print 200,(ast(is),is=1,24)
      do iyr = nyb,nye
        print 103,iyr,(xseval(is,iyr),is=1,24)
      end do
      print *,'  '
      print 201,(ast(is),is=25,49)
      do iyr = nyb,nye
        print 104,(xseval(is,iyr),is=25,49)
      end do

  100 format(i5,i7,i5,6i7)
  101 format(i5,48f6.2)
  102 format(49f6.2)
  103 format(i5,48f6.1)
  104 format(49f6.1)
  200 format(5x,24(3x,a2,1x))
  201 format(25(3x,a2,1x))
  202 format(48i4)
  210 format(5x,10(1x,a5,1x))
  211 format(i5,10f7.2)
  300 format(5x,48i7)
  301 format(i5,48i7)
      stop
      end
c********************************************************************** c
      subroutine mdian1(x,n,xmed)
      dimension x(n)
      call sort(n,x)
      n2 = n/2
      if( 2*n2 .eq. n) then
        xmed = 0.5*(x(n2)+x(n2+1))
      else
        xmed = x(n2+1)
      endif
      return
      end
c*********************************************************************
      subroutine rankx(x,n,y)
      dimension x(n),y(n)
      call sort(n,x)
      do i = 1,n
        y(i) = x(i)
      end do
      return
      end
c*********************************************************************
c
      subroutine sort(n,ra)
      dimension ra(n)
      l = n/2 + 1
      ir = n
   10 continue
      if( l .gt. 1 ) then
        l = l - 1
        rra = ra(l)
      else
        rra = ra(ir)
        ra(ir) = ra(1)
        ir = ir - 1
        if( ir .eq. 1 ) then
          ra(1) = rra
          return
        end if
      end if
      i = l
      j = l + l
   20 if( j .le. ir ) then
        if( j.lt.ir) then
          if( ra(j) .lt. ra(j+1)) j = j + 1
        end if
        if( rra .lt. ra(j) ) then
          ra(i) = ra(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if
        go to 20
      end if
      ra(i) = rra
      go to 10
      end
c
      subroutine sort2(n,ra,rb)
      dimension ra(n),rb(n)
      l = n/2 + 1
      ir = n
   10 continue
      if( l .gt. 1 ) then
        l = l - 1
        rra = ra(l)
        rrb = rb(l)
      else
        rra = ra(ir)
        rrb = rb(ir)
        ra(ir) = ra(1)
        rb(ir) = rb(1)
        ir = ir - 1
        if( ir .eq. 1 ) then
          ra(1) = rra
          rb(1) = rrb
          return
        end if
      end if
      i = l
      j = l + l
   20 if( j .le. ir ) then
        if( j.lt.ir) then
          if( ra(j) .lt. ra(j+1)) j = j + 1
        end if
        if( rra .lt. ra(j) ) then
          ra(i) = ra(j)
          rb(i) = rb(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if
        go to 20
      end if
      ra(i) = rra
      rb(i) = rrb
      go to 10
      end
eof
pgf90 -o source.exe source.f
./source.exe
#mv fort.50 stnmiss_mi.txt
rm source.f source.exe
