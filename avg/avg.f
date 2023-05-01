	  allocatable :: oacrude(:),oasum(:),oaavg(:)
	  allocatable :: obcrude(:),obsum(:),obavg(:),
     x	                 oabcrd(:),oabsum(:),oabavg(:),
     x	                 phicrd(:),phisum(:),phiavg(:),
     x	                 psicrd(:),psisum(:),psiavg(:),
     x	                 tcrude(:),tsum(:),tavg(:),
     x	                 avndo1(:),avndo2(:),
     x	                 avnro(:),avnroa(:),avnrob(:)

          integer(kind=8), allocatable:: ndo1(:),ndo2(:),
     x    ndo1s(:),ndo2s(:), nro(:),nros(:),
     x    nroa(:),nroas(:),nrob(:),nrobs(:),ndata(:)

          character*40 arg1, arg2
          integer (kind = 8) norm_jump_wait_time(2) 
          double precision avg_jump_waiting_time(2)


          n_arg_count = command_argument_count()
          if ( n_arg_count .eq. 0) then
             write(*,*) "Expects 2 command line arguments."
             write(*,*) "Arg1: time_grid"
             write(*,*) "Arg2: time_lag"
             stop
          endif
          
          CALL GET_COMMAND_ARGUMENT(1,arg1)
          CALL GET_COMMAND_ARGUMENT(2,arg2)

          read(arg1,*) itime
          read(arg2,*) timestep
          !write(*,*) itime, timestep

	  allocate(oacrude(-itime:itime),oasum(-itime:itime),
     x             oaavg(-itime:itime),oabcrd(-itime:itime),
     x             oabsum(-itime:itime),oabavg(-itime:itime),
     x              obcrude(-itime:itime),
     x             obsum(-itime:itime),obavg(-itime:itime),
     x    phicrd(-itime:itime),phisum(-itime:itime),
     x    phiavg(-itime:itime),psicrd(-itime:itime),
     x    psisum(-itime:itime),psiavg(-itime:itime),
     x    tcrude(-itime:itime),tsum(-itime:itime), tavg(-itime:itime),
     x    ndo1(-itime:itime),ndo2(-itime:itime),
     x    ndo1s(-itime:itime),ndo2s(-itime:itime),
     x    avndo1(-itime:itime),avndo2(-itime:itime),
     x    nro(-itime:itime),nros(-itime:itime),avnro(-itime:itime),
     x    nroa(-itime:itime),nroas(-itime:itime),avnroa(-itime:itime),
     x    nrob(-itime:itime),nrobs(-itime:itime),avnrob(-itime:itime),
     x    ndata(-itime:itime))

	  open(51,file='dist_angle.out',status='old')
          do i=1,24
             read(51,'(A)') 
          enddo
	  open(52,file='hbond.out',status='old')
          do i=1,26
             read(52,'(A)')
          enddo
	  open(53,file='jumptime_dist.dat',status='old')
          do i=1,18
             read(53,'(A)')
          enddo
	  open(54,file='avg_dist_angle.dat')
	  open(55,file='avg_hbond.dat')

	  do i=-itime,itime
	     oasum(i)=0.0d0
	     obsum(i)=0.0d0
	     oabsum(i)=0.0d0
	     phisum(i)=0.0d0
	     psisum(i)=0.0d0
	     tsum(i)=0.0d0
	     ndo1s(i)=0
	     ndo2s(i)=0
	     nros(i)=0
	     nroas(i)=0
	     nrobs(i)=0
	  enddo 


	 call print_info_avg_dist_angle
	 call print_info_avg_hbond 

	  
100	        read(51,*,IOSTAT=iostat1) k1, oa,ob,oab,psi,phi,theta
	        if (iostat1.lt.0) then 
	          go to 101
	        endif
	        oacrude(k1)=oa
	        oasum(k1)=oasum(k1)+oacrude(k1)

	        obcrude(k1)=ob
	        obsum(k1)=obsum(k1)+obcrude(k1)

                oabcrd(k1)=oab
                oabsum(k1)=oabsum(k1)+oabcrd(k1)

	        psicrd(k1)=psi
	        psisum(k1)=psisum(k1)+psicrd(k1)
	
	        phicrd(k1)=phi
	        phisum(k1)=phisum(k1)+phicrd(k1)

	        tcrude(k1)=theta
	        tsum(k1)=tsum(k1)+tcrude(k1)
	       

                   go to 100

101	        continue

200	        read(52,*,IOSTAT=iostat2)k2,nr1,nr2,nr3,nd1,nd2
	        if (iostat2.lt.0) then 
	          go to 201
	        endif
	        ndo1(k2)=nd1
	        ndo2(k2)=nd2

                ndo1s(k2)=ndo1s(k2)+ndo1(k2)
	        ndo2s(k2)=ndo2s(k2)+ndo2(k2)

	       nro(k2)=nr1
	       nroa(k2)=nr2
	       nrob(k2)=nr3

	       nros(k2)=nros(k2)+nro(k2)
	       nroas(k2)=nroas(k2)+nroa(k2)
	       nrobs(k2)=nrobs(k2)+nrob(k2)

	         go to 200

201	     continue

	    do i=-itime,itime
	      read(53,*) k3, ndata(i)
	    enddo 


	  do i=-itime,itime
	    if (i.ne.0) then
	     oaavg(i)=oasum(i)/float(ndata(i))
	     obavg(i)=obsum(i)/float(ndata(i))
	     psiavg(i)=psisum(i)/float(ndata(i))
	     phiavg(i)=phisum(i)/float(ndata(i))
	     tavg(i)=tsum(i)/float(ndata(i))
	     oabavg(i)=oabsum(i)/float(ndata(i))
	     avndo1(i)=float(ndo1s(i))/float(ndata(i))
	     avndo2(i)=float(ndo2s(i))/float(ndata(i))
	     avnro(i)=float(nros(i))/float(ndata(i))
	     avnroa(i)=float(nroas(i))/float(ndata(i))
	     avnrob(i)=float(nrobs(i))/float(ndata(i))
	    endif
	  enddo

	  do i=-itime,itime
	     if (i.ne.0) then
	       if (i.lt.0) then
	         time=-(timestep/2.)-timestep*float(abs(i)-1)
	       else
                 time=(timestep/2.)+timestep*float(i-1)
	       end if
	       write(54,*) time, oaavg(i),obavg(i),oabavg(i),
     x                      psiavg(i)-60.0d0,phiavg(i),tavg(i)
	       write(55,*) time,avnro(i),avnroa(i),avnrob(i),
     x                     avndo1(i),avndo2(i)
	     endif
	  enddo


	!Calculation of jump waiting time
         avg_jump_waiting_time(:) = 0.0d0
         norm_jump_wait_time(:) = 0
         do i=-itime,itime
             if(i .lt. 0) then;
                iflag = 1
             else
                iflag = 2
             endif

             avg_jump_waiting_time(iflag) = avg_jump_waiting_time(iflag)
     x       + float(ndata(i))*float(abs(i))*timestep       
             norm_jump_wait_time(iflag) = norm_jump_wait_time(iflag) 
     x                                     + ndata(i)
         enddo

          !write(*,*) avg_jump_waiting_time
         avg_jump_waiting_time(:) = avg_jump_waiting_time(:)/
     x                           float(norm_jump_wait_time(:))
         

         !write(*,*) float(norm_jump_wait_time)

         write(54,'(A)')"#*********************************************"

         write(54,*) "#Average jump waiting time from left (fs) = ", 
     x               avg_jump_waiting_time(1)
         write(54,*) "#Average jump waiting time from right (fs) = ",
     x               avg_jump_waiting_time(2)
         write(54,*) "#Average jump waiting time (fs) = ",
     x              0.5*(avg_jump_waiting_time(1)+
     x                   avg_jump_waiting_time(2))

         write(54,'(A)')"#*********************************************"



	  stop
	  end

	subroutine print_info_avg_dist_angle
	write(54,'(A)') "#<=== Average from dist_angle.out ===>#"
	write(54,'(A)') "#Column-1: Time grid point"
	write(54,'(A)') "#Column-2: Distance between O*---Oa"
	write(54,'(A)') "#Column-3: Distance between O*---Ob"
	write(54,'(A)') "#Column-4: Distance between Oa---Ob"
	write(54,'(A)') "#Column-5: Angle psi - 60 degree"
	write(54,'(A)') "#Column-6: Angle phi "
	write(54,'(A)') "#Column-7: Angle theta"
	end subroutine print_info_avg_dist_angle

        subroutine print_info_avg_hbond
        write(55,'(A)') "#<===== Average from hbond.out =====>#"
        write(55,'(A)') "#Column-1: Time grid point"
        write(55,'(A)') "#Column-2: H-bond accepted by O*"
        write(55,'(A)') "#Column-3: H-bond accepted by Oa"
        write(55,'(A)') "#Column-4: H-bond accepted by Ob"
        write(55,'(A)') "#Column-5: H-bond donated by H*"
        write(55,'(A)') "#Column-6: H-bond donated by H**"
        end subroutine print_info_avg_hbond


