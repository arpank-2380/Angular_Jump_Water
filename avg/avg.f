	  allocatable :: oacrude(:),oasum(:),oaavg(:)
	  allocatable :: obcrude(:),obsum(:),obavg(:),
     x	                 oabcrd(:),oabsum(:),oabavg(:),
     x	                 phicrd(:),phisum(:),phiavg(:),
     x	                 psicrd(:),psisum(:),psiavg(:),
     x	                 tcrude(:),tsum(:),tavg(:),
     x	                 avndo1(:),avndo2(:),
     x	                 avnro(:),avnroa(:),avnrob(:)

         allocatable :: oasum_sq(:),oavar(:),
     x                  obsum_sq(:),obvar(:),
     x                  oabsum_sq(:),oabvar(:), 
     x                  phisum_sq(:),phivar(:),
     x                  psisum_sq(:),psivar(:),
     x                  tsum_sq(:),tvar(:),
     x                  varndo1(:),varndo2(:),
     x                  varnro(:),varnroa(:),varnrob(:)         

          integer(kind=8), allocatable:: ndo1(:),ndo2(:),
     x    ndo1s(:),ndo2s(:), nro(:),nros(:),
     x    nroa(:),nroas(:),nrob(:),nrobs(:),ndata(:),
     x    ndo1s_sq(:),ndo2s_sq(:), nros_sq(:),
     x    nroas_sq(:),nrobs_sq(:)     

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

         !variables for variance calculations 
          allocate(oasum_sq(-itime:itime), oavar(-itime:itime),
     x             oabsum_sq(-itime:itime), oabvar(-itime:itime),  
     x             obsum_sq(-itime:itime), obvar(-itime:itime),
     x             phisum_sq(-itime:itime), phivar(-itime:itime),
     x             psisum_sq(-itime:itime), psivar(-itime:itime),
     x             tsum_sq(-itime:itime), tvar(-itime:itime),
     x             ndo1s_sq(-itime:itime), ndo2s_sq(-itime:itime),
     x             varndo1(-itime:itime), varndo2(-itime:itime),
     x             nros_sq(-itime:itime), varnro(-itime:itime),
     x             nroas_sq(-itime:itime), varnroa(-itime:itime),
     x             nrobs_sq(-itime:itime), varnrob(-itime:itime))

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
          open(56,file='var_dist_angle.dat')
          open(57,file='var_hbond.dat')

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

             oasum_sq(i)=0.0d0
             obsum_sq(i)=0.0d0
             oabsum_sq(i)=0.0d0
             phisum_sq(i)=0.0d0
             psisum_sq(i)=0.0d0
             tsum_sq(i)=0.0d0
             ndo1s_sq(i)=0
             ndo2s_sq(i)=0
             nros_sq(i)=0
             nroas_sq(i)=0
             nrobs_sq(i)=0             
	  enddo 


	 call print_info_avg_dist_angle(54)
	 call print_info_avg_hbond(55)
         call print_info_avg_dist_angle(56)
         call print_info_avg_hbond(57)    

	  
100	        read(51,*,IOSTAT=iostat1) k1, oa,ob,oab,psi,phi,theta
	        if (iostat1.lt.0) then 
	          go to 101
	        endif
	        oacrude(k1)=oa
	        oasum(k1)=oasum(k1)+oacrude(k1)
                oasum_sq(k1)=oasum_sq(k1)+oacrude(k1)*oacrude(k1)

	        obcrude(k1)=ob
	        obsum(k1)=obsum(k1)+obcrude(k1)
                obsum_sq(k1)=obsum_sq(k1)+obcrude(k1)*obcrude(k1)

                oabcrd(k1)=oab
                oabsum(k1)=oabsum(k1)+oabcrd(k1)
                oabsum_sq(k1)=oabsum_sq(k1)+oabcrd(k1)*oabcrd(k1)

	        psicrd(k1)=psi
	        psisum(k1)=psisum(k1)+psicrd(k1)
                psisum_sq(k1)=psisum_sq(k1)+psicrd(k1)*psicrd(k1)
	
	        phicrd(k1)=phi
	        phisum(k1)=phisum(k1)+phicrd(k1)
                phisum_sq(k1)=phisum_sq(k1)+phicrd(k1)*phicrd(k1)

	        tcrude(k1)=theta
	        tsum(k1)=tsum(k1)+tcrude(k1)
                tsum_sq(k1)=tsum_sq(k1)+tcrude(k1)*tcrude(k1)
	       

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
                ndo1s_sq(k2)=ndo1s_sq(k2)+ndo1(k2)*ndo1(k2)
                ndo2s_sq(k2)=ndo2s_sq(k2)+ndo2(k2)*ndo2(k2)

	       nro(k2)=nr1
	       nroa(k2)=nr2
	       nrob(k2)=nr3

	       nros(k2)=nros(k2)+nro(k2)
	       nroas(k2)=nroas(k2)+nroa(k2)
	       nrobs(k2)=nrobs(k2)+nrob(k2)
               nros_sq(k2)=nros_sq(k2)+nro(k2)*nro(k2)
               nroas_sq(k2)=nroas_sq(k2)+nroa(k2)*nroa(k2)
               nrobs_sq(k2)=nrobs_sq(k2)+nrob(k2)*nrob(k2)

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

             oavar(i) = variance(oasum_sq(i),oasum(i),ndata(i))
             obvar(i) = variance(obsum_sq(i),obsum(i),ndata(i))
             psivar(i) = variance(psisum_sq(i),psisum(i),ndata(i))
             phivar(i) = variance(phisum_sq(i),phisum(i),ndata(i))
             tvar(i) = variance(tsum_sq(i),tsum(i),ndata(i))
             oabvar(i) = variance(oabsum_sq(i),oabsum(i),ndata(i))

             varndo1(i) = variance(float(ndo1s_sq(i)),float(ndo1s(i)),
     x                    ndata(i))
             varndo2(i) = variance(float(ndo2s_sq(i)),float(ndo2s(i)),
     x                    ndata(i))
             varnro(i) = variance(float(nros_sq(i)), float(nros(i)),
     x                    ndata(i))
             varnroa(i) = variance(float(nroas_sq(i)), float(nroas(i)),
     x                    ndata(i))
             varnrob(i) = variance(float(nrobs_sq(i)), float(nrobs(i)),
     x                    ndata(i))


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
               write(56,*) time, oavar(i),obvar(i),oabvar(i),
     x                      psivar(i),phivar(i),tvar(i)
               write(57,*) time,varnro(i),varnroa(i),varnrob(i),
     x                     varndo1(i),varndo2(i)
              
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

	subroutine print_info_avg_dist_angle(i_file)
	write(i_file,'(A)') "#<=== Average from dist_angle.out ===>#"
	write(i_file,'(A)') "#Column-1: Time grid point"
	write(i_file,'(A)') "#Column-2: Distance between O*---Oa"
	write(i_file,'(A)') "#Column-3: Distance between O*---Ob"
	write(i_file,'(A)') "#Column-4: Distance between Oa---Ob"
	write(i_file,'(A)') "#Column-5: Angle psi - 60 degree"
	write(i_file,'(A)') "#Column-6: Angle phi "
	write(i_file,'(A)') "#Column-7: Angle theta"
	end subroutine print_info_avg_dist_angle

        subroutine print_info_avg_hbond(i_file)
        write(i_file,'(A)') "#<===== Average from hbond.out =====>#"
        write(i_file,'(A)') "#Column-1: Time grid point"
        write(i_file,'(A)') "#Column-2: H-bond accepted by O*"
        write(i_file,'(A)') "#Column-3: H-bond accepted by Oa"
        write(i_file,'(A)') "#Column-4: H-bond accepted by Ob"
        write(i_file,'(A)') "#Column-5: H-bond donated by H*"
        write(i_file,'(A)') "#Column-6: H-bond donated by H**"
        end subroutine print_info_avg_hbond

        function variance(sum_x2, sum_x, n)
        real sum_var_sq, sum_var
        integer(kind=8) n
        variance = (sum_x2 - (sum_x**2) / n) / (n - 1)
        return
        end


