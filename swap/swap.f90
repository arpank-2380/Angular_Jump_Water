  module swap_module
  use tools
  use io_module
  implicit none

  real, parameter :: mas=18.0d0  ! mass of a water molecule
  real, parameter :: mn=16.0d0   ! mass of oxygen
  real, parameter :: md=1.0d0    ! mass of hydrogen (deutoriam)
  integer, parameter :: nmode = 2            ! number of O-H vector per water molecule

!  real, allocatable :: cx(:,:),cy(:,:),cz(:,:)
  real, allocatable :: cmx(:,:),cmy(:,:),cmz(:,:), &
& snx(:,:),sny(:,:),snz(:,:), shx(:,:,:),shy(:,:,:),shz(:,:,:)
  real, allocatable :: rnh(:),rhn(:,:), dooa(:),doob(:), &
& rn(:,:,:),rhnp(:)   
  integer, allocatable ::  nhbr(:),itag(:,:),nswap(:),jtag(:,:)
  real  masinv, mddimension,  &
& ctheta, rhnx, rhny, rhnz, rmx, rmy, rmz, rnx, rny, rnz, &
&  rnhx, rnhy, rnhz 
  real time                  !just dummy variable  
  integer, allocatable ::  swaptime(:,:) 
  integer o, oa, ob, h, i, icount, iflag, im, &
&  irlxa, irlxb, itemp, itime, j, k, kk, ll, lnhbr, &
&  ltot, m,  mode, ndo1, ndo2, &
&  nonzero, nro, nroa, nrob, ntime, nzero, jend, ios, &
&  cmd_arg_count
  real a(3),b(3),c(3),v(3),rcom_a(3),rcom_b(3)
  real psi, phi, theta_act 
  character*2 symbol
  integer(kind=8), allocatable :: no_of_swap(:)
  integer(kind=8) itotal
  integer iblock

  contains

  subroutine swap 
  real run_start_swap, run_end_swap, theta, &
&       roo_ctheta_tmp, roo_ctheta, roo, roo_tmp, rho
  real tie_breaker_tmp, tie_breaker 
  integer tag_tmp 
 
  call cpu_time ( run_start_swap )
  call read_trajectory(traj_fmt,reading)

  cmx(:,:) = 0.0d0
  cmy(:,:) = 0.0d0
  cmz(:,:) = 0.0d0
  snx(:,:) = 0.0d0
  shx(:,:,:) = 0.0d0
  sny(:,:) = 0.0d0
  shy(:,:,:) = 0.0d0
  snz(:,:) = 0.0d0
  shz(:,:,:) = 0.0d0
  rnh(m) = 0.0d0
  rhn(:,:) = 0.0d0
  rn(:,:,:) = 0.0d0  
  rhnp(:) = 0.0d0
  itag(:,:) = 0
  jtag(:,:) = 0
  swaptime(:,:) = 0

  !write(*,*) "boxl=", boxl

   write(*,*) "Calculating distances..."

!       Calculate centre-of-mass coordinates
  do i=1,ntotal
     do j=1,nn
        cmx(i,j)=(cx(i,j)*mn+cx(i,nn+2*j-1)*md+cx(i,nn+2*j-0)*md)*masinv
        cmy(i,j)=(cy(i,j)*mn+cy(i,nn+2*j-1)*md+cy(i,nn+2*j-0)*md)*masinv
        cmz(i,j)=(cz(i,j)*mn+cz(i,nn+2*j-1)*md+cz(i,nn+2*j-0)*md)*masinv
     enddo
  enddo

!       Calculate site vectors
  do i=1,ntotal
     do  j=1,nn
         snx(i,j)=cx(i,j)-cmx(i,j)
         shx(i,j,1)=cx(i,nn+2*j-1)-cmx(i,j)
         shx(i,j,2)=cx(i,nn+2*j-0)-cmx(i,j)

         sny(i,j)=cy(i,j)-cmy(i,j)
         shy(i,j,1)=cy(i,nn+2*j-1)-cmy(i,j)
         shy(i,j,2)=cy(i,nn+2*j-0)-cmy(i,j)

         snz(i,j)=cz(i,j)-cmz(i,j)
         shz(i,j,1)=cz(i,nn+2*j-1)-cmz(i,j)
         shz(i,j,2)=cz(i,nn+2*j-0)-cmz(i,j)
     enddo
  enddo

  
  write(*,*) "Preparing neighbor list for each O-H mode..."
  do 100 i=1,ntotal
  do 200 j=1,nn

!	calculate intramolecular O-H distance
  do m=1,2
  rnhx=-snx(i,j)+shx(i,j,m)
  rnhy=-sny(i,j)+shy(i,j,m)
  rnhz=-snz(i,j)+shz(i,j,m)
  rnh(m)=sqrt(rnhx**2+rnhy**2+rnhz**2)
  enddo

!********************************************
  do 300 k=1,nn

!	  calculate intermolecular COM  distances
  rmx=cmx(i,k)-cmx(i,j)
  rmy=cmy(i,k)-cmy(i,j)
  rmz=cmz(i,k)-cmz(i,j)

!	Apply periodic boundary condition w.r.t. centre-of mass coordinates
  rmx=rmx-boxl*anint(rmx/boxl)
  rmy=rmy-boxl*anint(rmy/boxl)
  rmz=rmz-boxl*anint(rmz/boxl)

!	  calculate intermolecular H-O  distance
  do 16 m=1,2
  mode=2*(j-1)+m
  rhnx=rmx+snx(i,k)-shx(i,j,m)
  rhny=rmy+sny(i,k)-shy(i,j,m)
  rhnz=rmz+snz(i,k)-shz(i,j,m)
  rhn(k,mode)=sqrt(rhnx**2+rhny**2+rhnz**2)

!	Calculate O-O distances 
  rnx=rmx+snx(i,k)-snx(i,j)
  rny=rmy+sny(i,k)-sny(i,j)
  rnz=rmz+snz(i,k)-snz(i,j)
  rn(k,mode,i)=sqrt(rnx**2+rny**2+rnz**2)

!       elimination of self-interaction
    if (k.eq.j) then
     rhn(k,mode)=100.
     rn(k,mode,i)=100.
    endif
16	   continue
300	   continue
              

!       tagging of H-bonded neighbours:
  nhbr(:) = 0
  do 17 m=1,2
  mode=2*(j-1)+m
  ltot=0
  do  k=1,nn
!       tagging of H-bonded neighbours of an O atom using distance criteria only:
    if ((rn(k,mode,i).le.rn_cut).and.(rhn(k,mode).le.rhn_cut)) then
    ltot=ltot+1
    nhbr(ltot)=k
    end if
  end do

  tag_tmp = 0
  tie_breaker = 1000.0
  do 18 lnhbr=1,ltot
  kk=nhbr(lnhbr)
  rho=rhn(kk,mode)
  roo=rn(kk,mode,i)

!       calculation of H--O-O angle using properties of triangle
  roo_ctheta= (roo**2+rnh(m)**2-rho**2)/(2.*rnh(m)) 
  tie_breaker_tmp = roo - abs(roo_ctheta)
  ctheta = roo_ctheta/roo
  if ( ctheta .gt. 1.0d0) then
     ctheta = 1.0d0
  else if ( ctheta .lt. -1.0d0) then
     ctheta = -1.0d0
  endif 

  !using roo(1-costheta) to chose tagged Oxygen in case many satisfies distance-angle crieteria
  ! oxygen with smaller roo(1-costheta) is selected.
  if (theta .le. theta_cut ) then
     if (tie_breaker_tmp .lt. tie_breaker ) then
        tag_tmp = kk
        tie_breaker = tie_breaker_tmp
     endif
  endif
18 continue

! Calculation of H-O--O angle
  if (tag_tmp .ne. 0) then
!       application of angle criteria over atoms tagged by distance criteria:           
     itag(mode,i)=tag_tmp
  end if
17        continue
200       continue
100	  continue

 
                
write(*,*) "Identifying H-bond switching times for each O-H mode" 

!         finding the first nonzero-element of the neighbour list of each mode:
  do 500 mode=1,nh
   nonzero=1
  do while (itag(mode,nonzero).eq.0)
    nonzero=nonzero+1
  enddo  

    nswap(mode)=0
    itemp=itag(mode,nonzero)
    jtag(mode,0)=itag(mode,nonzero)
    swaptime(mode,0)=1

!         tagging a swap: when a particular mode changes its partner then it is a swap.
    do i=nonzero+1,ntotal
!         applying swap condition:partner change from 0 to 5(say) is not a swap but a bond making. similarly 
!         partner change from 5(say) to 0 is not a swap but a bond breaking. so are omitting those cases.
       if ((itemp.ne.itag(mode,i)).and.(itemp.ne.0).and.(itag(mode,i).ne.0)) then
          nswap(mode)=nswap(mode)+1
          j=nswap(mode)
          jtag(mode,j)=itag(mode,i) 
          swaptime(mode,j)=i 
          itemp=itag(mode,i)
       end if
    enddo
500 continue

!520 continue
!521 continue

write(*,*) "Calculating distances, angles, etc. for each H-bond switching within given time window."

!         total no. of swap will be different which is required for averaging: no_of_swap(i) basically keeps track 
!         the total number of swap considered at i th time-step.

!         calculations of parameters: o is the central O atom and oa & ob are the initial & final acceptor of the H-bond
   itotal=0
   do 600 mode=1,nh
   jend=nswap(mode)
   do 601 j=1,jend
!          tagging of oa & ob:
    oa=jtag(mode,j-1)
    ob=jtag(mode,j)
    im=mod(mode,2)
!          tagging of o:
     if (im.eq.1) then
     o=(mode-im)/2+1
     else
     o=mode/2
     endif
!           tagging of rotating Hydrogen:
    h=nn+mode
!           iflag is a 'flag variable'--when it is 1 then hydrogen bond swap will be considered, otherwise not
    iflag=1
!           'nzero' is the time when o changes its partner---it is taken as zero time though it is not the true zero time.
    nzero=swaptime(mode,j)
!           itime is the initial time & ntime is the final time up to which distance & angles are monitored.
!           initial time should not be less than the time when the previous swap have already been occured
!           similarly finaltime should be lesser than the time of next swap.
    itime=max(nzero-ndata,swaptime(mode,j-1),1)
    ntime=min(nzero+ndata-1,swaptime(mode,j+1)-1,ntotal)
!           irlxa & irlxb are the time of residence of oa & ob: it should be more than 'ndata' but using this condition 
!           total no. of swap may be very small. so depending on the situation we can use different limiting value of 
!           irlxa & irlxb:
    irlxa=swaptime(mode,j)-swaptime(mode,j-1)
    irlxb=swaptime(mode,j+1)-swaptime(mode,j)

    ! switchin event would be discarded if there is (i) a previous switch and (ii) after switch within "dangling_time_grid"  
    if((irlxb.le. 0)) then
     iflag = 0 
    else if((irlxa.le.dangling_time_grid).and.(irlxb.le.dangling_time_grid)) then
     iflag=0
    ! switchin event would be discarded if thre is a previous switch within "dangling_time_grid" and actually that partner
    ! is same as the partner what would be after switch. This event is after dangling comming back to ex-partner. 
    else if (irlxa.le.dangling_time_grid) then
     if ( ( itag(mode,itime-1) .eq.ob ) .or. ( itag(mode,itime-1) .eq. 0 ) ) then
        iflag=0
     endif
    ! switchin event would be discarded if thre is a after switch within "dangling_time_grid" and actually that partner
    ! is same as the partner before previous switch. Also, this event signifies after dangling comming back to ex-partner. 
    else if (irlxb.le.dangling_time_grid) then
      if ( ( itag(mode,ntime+1) .eq. oa ) .or. ( itag(mode,ntime+1).eq.0 ) ) then
        iflag=0
      endif
    endif


!           before & after swap oa & ob would be the neighbour of o. but to increase the no. of swap we have to relax the condition

    icount=0
    if(iflag.ne.0) then
      do 620 i=itime,nzero-1
        if (itag(mode,i).ne.oa) then
          icount=icount+1
        endif
        if (icount.gt.0) then
          iflag=0
          go to 620
        endif
620   continue
    endif

    icount=0
    if (iflag.ne.0) then
      do 621 i=nzero,ntime
        if (itag(mode,i).ne.ob) then
           icount=icount+1
        endif
        if (icount.gt.0) then
           iflag=0
           go to 621
        endif
621   continue
    endif             


    if (iflag.eq.0) then
        goto 602
    else 

! Before the swap: one additional point is taken for running average smoothing with 3 points

    do i=itime,nzero-1
       no_of_swap(i-nzero)=no_of_swap(i-nzero)+1
       rcom_a(1)=cmx(i,oa)-cmx(i,o)
       rcom_a(2)=cmy(i,oa)-cmy(i,o)
       rcom_a(3)=cmz(i,oa)-cmz(i,o)
       rcom_b(1)=cmx(i,ob)-cmx(i,o)
       rcom_b(2)=cmy(i,ob)-cmy(i,o)
       rcom_b(3)=cmz(i,ob)-cmz(i,o)
!	  applying periodic boundary condition w.r.t C.O.M 
       do kk=1,3
         rcom_a(kk)=rcom_a(kk)-boxl*anint(rcom_a(kk)/boxl)
         rcom_b(kk)=rcom_b(kk)-boxl*anint(rcom_b(kk)/boxl)
       enddo

!         shifting the origin of vector space at the COM of o:
       a(1)=snx(i,oa)+rcom_a(1)
       a(2)=sny(i,oa)+rcom_a(2)
       a(3)=snz(i,oa)+rcom_a(3)
       b(1)=snx(i,ob)+rcom_b(1)
       b(2)=sny(i,ob)+rcom_b(2)
       b(3)=snz(i,ob)+rcom_b(3)
       c(1)=snx(i,o)
       c(2)=sny(i,o)
       c(3)=snz(i,o)
        if (im.eq.1) then
        v(1)=shx(i,o,1)-c(1)
        v(2)=shy(i,o,1)-c(2)
        v(3)=shz(i,o,1)-c(3)
        else
        v(1)=shx(i,o,2)-c(1)
        v(2)=shy(i,o,2)-c(2)
        v(3)=shz(i,o,2)-c(3)
        endif

       call projection_on_plane(c,a,b,v,psi,phi,theta_act)
       write(52,'(I7,5X,6G14.6)')i-nzero, rn(oa,mode,i), rn(ob,mode,i), &
&                 rn(oa,2*ob,i), psi, phi, theta_act



! Calculating the number of H-bond accepted by O*,Oa & Ob before swap:
! No of H-bond accepted by oa= No. of modes whose neighbour is oa
! No of H-bond accepted by ob= No. of modes whose neighbour is ob
       nro=0
       nroa=0
       nrob=0
       do ll=1,nh
          if (itag(ll,i).eq.o) then  
             nro=nro+1
         
          else if (itag(ll,i).eq.oa) then
             nroa=nroa+1
          else if (itag(ll,i).eq.ob) then
             nrob=nrob+1
          end if  
       enddo
       

! Calculating the number of H-bond donated by O* before swap:
! if the neighbour of a particular mode at a particular time is 0, then it is not H-bonded (donated), otherwise 
! it will have one neighbour, i.e. no. of H-bond donated will be 1 at that time.
! ndo1 is the no. of H-bond donated by the rotating Hydrogen and ndo2 is that donated by other Hydrogen

       if (im.eq.1) then
         if (itag(mode,i).ne.0) then
           ndo1=1
         else
           ndo1=0
         endif

         if (itag(mode+1,i).ne.0) then
           ndo2=1
         else
           ndo2=0
         endif

      else
          if (itag(mode,i).ne.0) then
           ndo1=1
         else
           ndo1=0
         endif

         if (itag(mode+1,i).ne.0) then
           ndo2=1
         else
           ndo2=0
         endif
      endif

       
      write(53,'(I10, 3x, I8,4x,I8,5x,I8,7x,I8,8x,I8)') i-nzero,nro,nroa,nrob,ndo1,ndo2


   enddo

! After the swap: 1 extra point taken for smoothing
! similar analysis was done for the time after the swap

   do i=nzero,ntime
       no_of_swap(i-nzero+1)=no_of_swap(i-nzero+1)+1
       rcom_a(1)=cmx(i,oa)-cmx(i,o)
       rcom_a(2)=cmy(i,oa)-cmy(i,o)
       rcom_a(3)=cmz(i,oa)-cmz(i,o)
       rcom_b(1)=cmx(i,ob)-cmx(i,o)
       rcom_b(2)=cmy(i,ob)-cmy(i,o)
       rcom_b(3)=cmz(i,ob)-cmz(i,o)
!      applying periodic boundary condition w.r.t C.O.M
       do kk=1,3
         rcom_a(kk)=rcom_a(kk)-boxl*anint(rcom_a(kk)/boxl)
         rcom_b(kk)=rcom_b(kk)-boxl*anint(rcom_b(kk)/boxl)
       enddo


       a(1)=snx(i,oa)+rcom_a(1)
       a(2)=sny(i,oa)+rcom_a(2)
       a(3)=snz(i,oa)+rcom_a(3)
       b(1)=snx(i,ob)+rcom_b(1)
       b(2)=sny(i,ob)+rcom_b(2)
       b(3)=snz(i,ob)+rcom_b(3)
       c(1)=snx(i,o)
       c(2)=sny(i,o)
       c(3)=snz(i,o)
        if (im.eq.1) then
        v(1)=shx(i,o,1)-c(1)
        v(2)=shy(i,o,1)-c(2)
        v(3)=shz(i,o,1)-c(3)
        else
        v(1)=shx(i,o,2)-c(1)
        v(2)=shy(i,o,2)-c(2)
        v(3)=shz(i,o,2)-c(3)
        endif

        call projection_on_plane(c,a,b,v, psi,phi,theta_act)
        write(52,'(I7,5X,6G14.6)') i-nzero+1, rn(oa,mode,i), rn(ob,mode,i), &
&                  rn(oa,2*ob,i), psi, phi, theta_act

! Calculating the number of H-bond accepted by O*,Oa & Ob after swap:
       nro=0
       nroa=0
       nrob=0
       do ll=1,nh
          if (itag(ll,i).eq.o) then
             nro=nro+1

          else if (itag(ll,i).eq.oa) then
             nroa=nroa+1
          else if (itag(ll,i).eq.ob) then
             nrob=nrob+1
          end if
       enddo

!  Calculating the number of H-bond donated by O* after swap :
       if (im.eq.1) then
          if (itag(mode,i).ne.0) then
            ndo1=1
          else
            ndo1=0
          endif

          if (itag(mode+1,i).ne.0) then
            ndo2=1
          else
            ndo2=0
          endif

       else
           if (itag(mode,i).ne.0) then
            ndo1=1
          else
            ndo1=0
          endif

          if (itag(mode+1,i).ne.0) then
            ndo2=1
          else
            ndo2=0
          endif
       endif


        
       write(53,'(I10, 3x, I8,4x,I8,5x,I8,7x,I8,8x,I8)') i-nzero+1,nro,nroa,nrob,ndo1,ndo2

   enddo
   itotal=itotal+1 
   endif
602	  continue 
601	  continue
600	  continue
         
 
 write(*,*)'Maximum no of swaps considered in Block:',iblock, "is:", itotal
 call cpu_time( run_end_swap )
 write(*,*) "Swap calculation for block", iblock," finished in ", run_end_swap-run_start_swap, " s."
 write(*,*) "***************************************************************************************"

 return
 end subroutine swap       

 end module swap_module                                                                                                                                                
