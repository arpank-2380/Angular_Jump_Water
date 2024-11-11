 module io_module
 implicit none

 save

 character*200 :: infil, traj_file
 integer  ::     nn, nh, nh2o, ntotal, ntotal2, natom, last_frame, first_frame, step_frame=1
 integer  :: start_read_frame=0
 integer  ::     ndata, dangling_time_grid=0, nblock
 integer (kind = 8) :: trajectory_line
 real :: boxl, boxl2, boxl3, time_lag, output_time_lag 
 real :: rhn_cut=2.45, rn_cut=3.5, theta_cut=45.0, max_mem=1.0
 real :: box_alpha, box_beta, box_gamma
 real, allocatable :: cx(:,:),cy(:,:),cz(:,:)
 character*6 :: traj_fmt
 logical :: file_exists, reading


 contains

 subroutine read_trajectory(traj_fmt,reading)
 character*6, intent(in) :: traj_fmt
 logical reading
 integer i, j, kk, ios, iframe, rframe, remainder_frame
 real x, y, z, time, read_start_time, read_end_time, read_frac
 character*6 symbol, junk1, junk2
 character*30 junk3
 integer junk4, junk5
 real junk6, junk7
 logical store_data
 
 call cpu_time( read_start_time )
 
 cx(:,:) = 0.0d0
 cy(:,:) = 0.0d0
 cz(:,:) = 0.0d0

 write(*,*) "Reading trajectory ..."
 iframe = 0
 rframe = 0
 select case ( traj_fmt )
          case ( "old" )
              if (.not. reading ) then
                 start_read_frame = 1
                 trajectory_line = 0
              endif
              reading = .true.
              do i=start_read_frame,last_frame
                 if ( iframe .eq. ntotal ) exit 
                 rframe = rframe + 1
                 store_data = .false. 
                 if (( i .ge. first_frame) .and. (i .le. last_frame) ) then
                    remainder_frame = mod(i-first_frame ,step_frame)
                    if ( remainder_frame .eq. 0 ) then
                       store_data = .true.
                       iframe = iframe + 1
                    endif
                 endif

                 do j=1,natom
                    trajectory_line=trajectory_line+1
                    read(51,*,iostat=ios)time,symbol,x,y,z
                    if (ios.ne.0) then
                        write(*,*) "Cannot read trajectory line =", trajectory_line
                        stop
                    endif
                    if ( store_data ) then 
                       kk=j/3
                       If (mod(j,3).eq.1) then
                          cx(iframe,kk+1)=x
                          cy(iframe,kk+1)=y
                          cz(iframe,kk+1)=z
                       else if (mod(j,3).eq.2) then
                          cx(iframe,nn+1+2*kk)=x
                          cy(iframe,nn+1+2*kk)=y
                          cz(iframe,nn+1+2*kk)=z
                       else
                          cx(iframe,nn+2*kk)=x
                          cy(iframe,nn+2*kk)=y
                          cz(iframe,nn+2*kk)=z
                       endif
                    endif
               enddo
              enddo

          case( "gro" )
              if ( .not. reading ) then
                 read(51,*) junk1, boxl, boxl2, boxl3, box_alpha, box_beta, box_gamma, junk2, junk4, junk5
                 trajectory_line=1
                 start_read_frame=1 
              endif
              reading = .true. 
              !write(*,*) "started reading frame:", start_read_frame, "at line ", trajectory_line+1
              do i=start_read_frame,last_frame
                 !write(*,*) iframe, rframe
                 if ( iframe .eq. ntotal ) exit
                 rframe = rframe + 1
                 store_data = .false.
                 if ( ( i .ge. first_frame ) .and. (i .le. last_frame) ) then
                    remainder_frame = mod (i-first_frame, step_frame)
                    if ( remainder_frame .eq. 0 ) then
                       store_data = .true.
                       iframe = iframe + 1
                    endif
                 endif

                 do j=1,natom
                    trajectory_line=trajectory_line+1
                    read(51,*,iostat=ios) junk1, junk4, symbol, junk2, junk3, junk5, &
                                      & x,y,z, junk6, junk7
                    if (ios.ne.0) then
                        write(*,*) "Cannot read trajectory line =", trajectory_line
                        stop
                    endif

                    if ( store_data ) then
                        kk=j/3
                        If (mod(j,3).eq.1) then
                           cx(iframe,kk+1)=x
                           cy(iframe,kk+1)=y
                           cz(iframe,kk+1)=z
                        else if (mod(j,3).eq.2) then
                           cx(iframe,nn+1+2*kk)=x
                           cy(iframe,nn+1+2*kk)=y
                           cz(iframe,nn+1+2*kk)=z
                        else
                           cx(iframe,nn+2*kk)=x
                           cy(iframe,nn+2*kk)=y
                           cz(iframe,nn+2*kk)=z
                        endif
                    endif
               enddo
               read(51,'(A)')
               trajectory_line = trajectory_line + 1
              enddo

          case ( "xyz" )
              if (.not. reading ) then
                 start_read_frame = 1
                 trajectory_line = 0
              endif
              reading = .true.
              do i=start_read_frame,last_frame
                 if ( iframe .eq. ntotal ) exit
                 rframe = rframe + 1
                 store_data = .false.
                 if (( i .ge. first_frame) .and. (i .le. last_frame) ) then
                    remainder_frame = mod(i-first_frame ,step_frame)
                    if ( remainder_frame .eq. 0 ) then
                       store_data = .true.
                       iframe = iframe + 1
                    endif
                 endif

                 read(51,'(A)')
                 read(51,'(A)')
                 trajectory_line=trajectory_line+2
                 do j=1,natom
                    trajectory_line=trajectory_line+1
                    read(51,*,iostat=ios)symbol,x,y,z
                    if (ios.ne.0) then
                        write(*,*) "Cannot read trajectory line =", trajectory_line
                        stop
                    endif
                    if ( store_data ) then
                       kk=j/3
                       If (mod(j,3).eq.1) then
                          cx(iframe,kk+1)=x
                          cy(iframe,kk+1)=y
                          cz(iframe,kk+1)=z
                       else if (mod(j,3).eq.2) then
                          cx(iframe,nn+1+2*kk)=x
                          cy(iframe,nn+1+2*kk)=y
                          cz(iframe,nn+1+2*kk)=z
                       else
                          cx(iframe,nn+2*kk)=x
                          cy(iframe,nn+2*kk)=y
                          cz(iframe,nn+2*kk)=z
                       endif
                    endif
                 enddo
              enddo



          case( "pdb" )
              if ( .not. reading ) then
                 read(51,*) junk1, boxl, boxl2, boxl3, box_alpha, box_beta, box_gamma, junk2, junk4, junk5
                 trajectory_line=1
                 start_read_frame=1 
              endif
              reading = .true. 
              !write(*,*) "started reading frame:", start_read_frame, "at line ", trajectory_line+1
              do i=start_read_frame,last_frame
                 !write(*,*) iframe, rframe
                 if ( iframe .eq. ntotal ) exit
                 rframe = rframe + 1
                 store_data = .false.
                 if ( ( i .ge. first_frame ) .and. (i .le. last_frame) ) then
                    remainder_frame = mod (i-first_frame, step_frame)
                    if ( remainder_frame .eq. 0 ) then
                       store_data = .true.
                       iframe = iframe + 1
                    endif
                 endif

                 do j=1,natom
                    trajectory_line=trajectory_line+1
                    read(51,'(30A)',advance='no') junk3
                    read(51,*,iostat=ios) x,y,z, junk6, junk7
                    if (ios.ne.0) then
                        write(*,*) "Cannot read trajectory line =", trajectory_line
                        stop
                    endif

                    if ( store_data ) then
                        kk=j/3
                        If (mod(j,3).eq.1) then
                           cx(iframe,kk+1)=x
                           cy(iframe,kk+1)=y
                           cz(iframe,kk+1)=z
                        else if (mod(j,3).eq.2) then
                           cx(iframe,nn+1+2*kk)=x
                           cy(iframe,nn+1+2*kk)=y
                           cz(iframe,nn+1+2*kk)=z
                        else
                           cx(iframe,nn+2*kk)=x
                           cy(iframe,nn+2*kk)=y
                           cz(iframe,nn+2*kk)=z
                        endif
                    endif
                 enddo
               read(51,'(A)')
               trajectory_line = trajectory_line + 1
              enddo
          case default
               write(*,'(A)') " Didn't understand specified trajectory_format '"//trim(traj_fmt)//"'."
               stop
       end select
       start_read_frame = start_read_frame + rframe 
       write(*,*) "Finished reading frame", start_read_frame-1
       write(*,*) "No of frames read = ", rframe
       write(*,*) "No of frames stored = ", iframe
   call cpu_time ( read_end_time )
   write(*,*) "Reading trajectory finished in", read_end_time-read_start_time, "s."
   write(*,*) "***************************************************************************************"        
 return
 end subroutine read_trajectory

 subroutine read_input_file
 ! input file read related 
 character infil_input*200, label*100
 integer posstart, posend, poseq, poscomment, ios, iblock, total_frames
 integer (kind = 8) tot_real_numbers
 real mem_req
 !infil= "swap.input" 
 

 INQUIRE(file=infil,exist=file_exists)
 if ( file_exists ) then
    open(50, file=infil, status='old')
 else
    write(*,'(A)') "Input file: "//trim(infil)//" doesn't exist."
    stop
 end if

 ios = 0
 do while (ios .eq. 0)
    read(50,'(A)',iostat=ios) infil_input
    if (ios.eq.0) then
        poscomment = scan(infil_input,'#!')
        if(poscomment.ne.0) then
           infil_input=infil_input(1:poscomment-1)
        endif
        posstart=verify(infil_input,'   ')
        
        if (posstart.eq.0) then
          cycle
        end if

        posend = scan(infil_input(posstart:),' = ')+posstart-2
        poseq = scan(infil_input,'=',.TRUE.)

        if (poseq.eq.0) then
          write(*,'(A)') "Error in "//infil//" '=' sign not found"
          stop
        else
          label = infil_input(posstart:posend)
          infil_input=adjustl(infil_input(poseq+1:))

          select case (label)
          case("trajectory")
              read(infil_input,'(A)',iostat=ios) traj_file
          case("no_water")
              read(infil_input,*,iostat=ios) nh2o    
          case("last_frame") 
              read(infil_input,*,iostat=ios) last_frame
          case("time_lag")
              read(infil_input,*,iostat=ios) time_lag     
          case("first_frame")
              read(infil_input,*,iostat=ios) first_frame 
          case("step_frame")
              read(infil_input,*,iostat=ios) step_frame
          case("rOH_cut") 
              read(infil_input,*,iostat=ios) rhn_cut
          case("rOO_cut")
             read(infil_input,*,iostat=ios) rn_cut             
          case ("angle_cut")
              read(infil_input,*,iostat=ios) theta_cut
          case("time_grid")
              read(infil_input,*,iostat=ios) ndata
          case("dangling_time_cut")
              read(infil_input,*,iostat=ios) dangling_time_grid
          case("box_length")
              read(infil_input,*,iostat=ios) boxl
          case("trajectory_format")
              read(infil_input,*,iostat=ios) traj_fmt
          case("memory")
              read(infil_input,*,iostat=ios) max_mem
          case default
              write(*,'(A)') "Keyword not understood" 
          end select
        endif 
    endif
 enddo    

 output_time_lag = time_lag * float(step_frame)
 total_frames=int((last_frame-first_frame)/step_frame)+1
 write(*,*) "**************************************************************************************"
 write(*,*) "Frames to be considered: ", total_frames
 tot_real_numbers=int8(nh2o)*2*int8(nh2o*total_frames)+int8(8*3*nh2o*total_frames)+ &
       &           int8(8*nh2o*total_frames)+int8(2*nh2o*nh2o)+int8(5*nh2o) + &
       &           int8(2*nh2o) + 4*int8(nh2o)*int8(total_frames)
 mem_req=float(tot_real_numbers*sizeof(boxl))/float(1024*1024*1024)
 write(*,*) "Memory bottleneck(GB):", mem_req
 nblock = int(mem_req/max_mem) + 1
 write(*,*) "Trajectory will be divided into", nblock, "blocks."

 if(mod(total_frames,nblock) .eq. 0 ) then
   ntotal = int(total_frames/nblock)
   ntotal2 = ntotal
 else
   ntotal = int(total_frames/nblock)+1
   ntotal2 = total_frames - ntotal*(nblock-1)
 end if
 
 write(*,*) "Block index         Number of frames"
 do iblock=1,nblock
    if (iblock .ne. nblock) then
       write(*,'(I6,20x,I10)') iblock, ntotal
    else
       write(*,'(I6,20x,I10)') iblock, ntotal2
    endif
 enddo

 tot_real_numbers=int8(nh2o)*2*int8(nh2o*ntotal)+int8(8*3*nh2o*ntotal)+ &
       &           int8(8*nh2o*ntotal)+int8(2*nh2o*nh2o)+int8(5*nh2o) + &
       &           int8(2*nh2o) + 4*int8(nh2o)*int8(ntotal)

 mem_req=float(tot_real_numbers*sizeof(boxl))/float(1024*1024*1024)
 write(*,*) "Memory bottleneck reduced to (GB):", mem_req
   

  end subroutine

  subroutine write_credit(file_unit)
  integer, intent(in) :: file_unit
     write(file_unit,'(A)') "###############################################################################"
     write(file_unit,'(A)') "#                   This output is produced by swap.f90                       #"
     write(file_unit,'(A)') "#                       written by Arpan Kundu at                             #"
     write(file_unit,'(A)') "#                 Indian Institute of Technology Kanpur                       #"
     write(file_unit,'(A)') "###############################################################################"
     !write(file_unit,'(A)') "#TRAJECTORY = "//trim(traj_file)// " ."
 
  end subroutine write_credit

  subroutine print_info_dist_angle
     write(52,'(A)') "#TRAJECTORY = "//trim(traj_file)// " ."
     write(52,'(A14,I10,A18,I12)') "#First frame =" , first_frame , ".     Last frame= ", last_frame 
     write(52,'(A16,F10.2)') "#time_lag(fs) = ", output_time_lag 
     write(52, '(A)') "#*******************************************************************************"
     write(52, '(A)') "#                    ============ Information ============                     *"
     write(52, '(A)') "#This file registers distance and angle coordinates for each swap considered by*"
     write(52, '(A)') "#swap.f90 program. When the time grid point changes the value from negative to *"
     write(52, '(A)') "#positive it signifies a new swap. Notice, depending on the choice of the      *"
     write(52, '(A)') "#  dangling_time_grid variable different swaps may have different time windows *"
     write(52, '(A)') "#because some events are faster and some are slower. Therefore each time grid  *"
     write(52, '(A)') "# would have different number of data-poits which is also stored at            *"
     write(52, '(A)') "# jumptime_dist.out (jumptime distribution) file.                              *"
     write(52, '(A)') "# Definition of O*, Oa, Ob and angle phi, psi and theta can be found in Laage -*"
     write(52, '(A)') "# Hynes original Science paper and also on Arpan's master's thesis             *"
     write(52, '(A)') "#*******************************************************************************"
     write(52, '(A)') "#                                                                               "  
     write(52, '(A)') "#----------------------------------------------------------------------------------------------"
     write(52, '(A)') "#Time-grid     R(O*-Oa)      R(O*-Ob)      R(Oa-Ob)     angle_psi     angle_phi    angle_theta "
     write(52, '(A)') "#----------------------------------------------------------------------------------------------"           

  end subroutine print_info_dist_angle

    subroutine print_info_hbond
     write(53,'(A)') "#TRAJECTORY = "//trim(traj_file)// " ."
     write(53,'(A14,I10,A18,I12)') "#First frame =" , first_frame , ".     Last frame= ", last_frame
     write(53,'(A16,F10.2)') "#time_lag(fs) = ", output_time_lag
     write(53, '(A)') "#*******************************************************************************"
     write(53, '(A)') "#                    ============ Information ============                     *"
     write(53, '(A)') "#This file registers the state of hydrogen bond receptor and donor for each    *" 
     write(53, '(A)') "#swap considered by swap.f90 program. State 1 means donated (or accepted) and 0*"
     write(53, '(A)') "# means otherwise                                                              *"
     write(53, '(A)') "#   When the time grid point changes the value from negative to positive it    *"
     write(53, '(A)') "#     signifies a new swap. Notice, depending on the choice of the             *"
     write(53, '(A)') "#  dangling_time_grid variable different swaps may have different time windows *"
     write(53, '(A)') "#because some events are faster and some are slower. Therefore each time grid  *"
     write(53, '(A)') "# would have different number of data-poits which is also stored at            *"
     write(53, '(A)') "# jumptime_dist.out (jumptime distribution) file.                              *"
     write(53, '(A)') "# Definition of O*, Oa, Ob and angle phi, psi and theta can be found in Laage -*"
     write(53, '(A)') "# Hynes original Science paper and also on Arpan's master's thesis             *"
     write(53, '(A)') "#*******************************************************************************"
     write(53, '(A)') "#                                                                               "
     write(53, '(A)') "#-----------------------------------------------------------------------------------------"
     write(53, '(A)') "#Time-grid     Acc(O*)        Acc(Oa)        Acc(Ob)      Donated(H*)      Donated(H**)   "
     write(53, '(A)') "#-----------------------------------------------------------------------------------------"

  end subroutine print_info_hbond


  subroutine print_info_datapoint
    write(54,'(A)') "#TRAJECTORY = "//trim(traj_file)// " ."
    write(54,'(A14,I10,A18,I12)') "#First frame =" , first_frame , ".     Last frame= ", last_frame
    write(54,'(A16,F10.2)') "#time_lag(fs) = ", output_time_lag
    write(54,'(A)') "#******************************************************************************"
    write(54,'(A)') "# This file containts the number of datapoints available for each time-grid   *"
    write(54,'(A)') "# point where the origin is the transition point of a H-bond switching event. *"
    write(54,'(A)') "# This information can be latter used to calculate jump statistics (e.g using *"
    write(54,'(A)') "# avg.f).                                                                     *"
    write(54,'(A)') "#******************************************************************************"
    write(54,'(A)') "#                                                                              "
    write(54,'(A)') "#-----------------------------------------------                               "
    write(54,'(A)') "#Time grid points         Number of data points                                "
    write(54,'(A)') "#-----------------------------------------------                               "
 end subroutine print_info_datapoint



  end module io_module
