program main
  use tools
  use io_module
  use swap_module
  implicit none
  real run_start_main, run_end_main

  call cpu_time( run_start_main )
  cmd_arg_count = command_argument_count()


  if ( cmd_arg_count .eq. 0) then
     infil= "swap.input"
  else
     CALL GET_COMMAND_ARGUMENT(1,infil)
  endif

  call write_credit(6)


  call read_input_file

  inquire (file=traj_file,exist=file_exists) 
  if ( file_exists ) then 
     open (unit=51,file=traj_file,status='old')
  else
     write(*,*) "Trajectory not found!"
     stop
  endif

  open(52,file='dist_angle.out')
  open(53,file='hbond.out')
  open(54,file='jumptime_dist.dat')
  
  call write_credit(52)
  call write_credit(53)
  call write_credit(54)        

  call print_info_dist_angle
  call print_info_hbond
  call print_info_datapoint
  

  masinv=1.0d0/mas

  natom=3*nh2o

  nn=nh2o
  nh=nmode*nh2o

  reading = .false.

  write(*,*) "Allocating memory....."

  allocate(&
&          cx(ntotal,natom),cy(ntotal,natom),cz(ntotal,natom), &
&          cmx(ntotal,natom),cmy(ntotal,natom),cmz(ntotal,natom),&
&          snx(ntotal,natom),sny(ntotal,natom),snz(ntotal,natom),&
&          shx(ntotal,natom,nmode),shy(ntotal,natom,nmode),&
&          shz(ntotal,natom,nmode) &
&          )
  allocate( no_of_swap(-1*ndata:ndata) )

  allocate ( rnh(nh),rhn(nn,nh),rn(nn,nh,ntotal), rhnp(nn), &           
&            nhbr(nh),itag(nh,ntotal) )

  allocate(nswap(nh),jtag(nh,0:ntotal),swaptime(nh,0:ntotal))

  write(*,*) "Memory allocation finished..."

  no_of_swap(:) = 0

  do iblock = 1, nblock
     if ( iblock .eq. 1 ) reading = .false.
     if ( iblock .eq. nblock ) ntotal = ntotal2
     call swap 
  enddo

 do i=-ndata,ndata
   write(54,'(I10,18x,I12)') i,no_of_swap(i)
 enddo
 
 
 close(51)
 close (52)
 close (53)
 close(54)

 write(*,*) "-----------------------------------------------------"
 write(*,*) " Total number of swap considered = ", no_of_swap(1)
 write(*,*) "-----------------------------------------------------" 
 call cpu_time( run_end_main )
 write(*,*) "Job finished in ", run_end_main-run_start_main, " s."

 stop
 end program main                                                                                                                                                      
