! Capt'n General testing program

    PROGRAM GENCAP
    implicit none
    character*300 :: modfile, filename
    character*100 :: outfile(7)
    double precision :: mx, jx, sigma_0, capped_sd, maxcapture, rho0, usun, u0, vesc
    double precision :: maxcap
    integer :: nq(7), nv(7), i, j, num_isotopes, spin_dependency, cpl
    character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
                        "c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]
	
    ! Choose velocity and momentum transfer powers in differential cross-section
    nq = [0,-1,1,2,0,0,0]
    nv = [0,0,0,0,-1,1,2]

    outfile = ['const.dat','qm1--.dat','q1---.dat','q2---.dat','vm1--.dat','v1---.dat','v2---.dat']

    ! Choose solar model file
    !modfile = "solarmodels/model_gs98_nohead.dat"
    !modfile = "solarmodels/struct_b16_agss09_nohead.dat"
    modfile = "solarmodels/struct_b16_agss09_nohead_Tsmoothed.dat" !temperature smoothed to not nonsense

    ! number of isotopes capt'n will loop over in the calculation: up to 29 isotopes
    num_isotopes = 29

    ! zero for spin_independent, one for spin_dependent
    spin_dependency = 0

    print*
    print*, "Initializing Capt'n General..."

    ! Initialise capture calculations
    rho0 = 0.4d0 !
    usun = 235.d0 !
    u0 = 235.d0 !
    vesc = 550.d0 !
    call captn_init(modfile,rho0,usun,u0,vesc)
    
    print*
    print*, "Running Capt'n General..."

    do j = 1,7
      print*
      print*, "Power: ", outfile(j)

      open(94,file = outfile(j))
      write(94,*) "Number of Isotopes: ", num_isotopes
      write(94,*) "Spin Dependency: ", spin_dependency
      write(94,*) "Power: ", outfile(j)
      write(94,*) "Sigma_0 | ", "DM Mass | ", "Captured Dark Matter"
      do i = 1,10
        mx = 5.d0 + dble(i)/5.
        sigma_0 = 10d0**(-45+dble(i)/2.)
        print*
        print*, "mx: ", mx, "sigma_0: ", sigma_0, "cm^2"

        print*, "Geometrical limit on capture rate: ", maxcap(mx), "s^-1"

        call captn_general(mx,sigma_0,num_isotopes,nq(j),nv(j),spin_dependency,capped_sd)
        print*, "Capture rate: ", capped_sd, "s^-1"

        write(94,*) sigma_0, mx, capped_sd
      end do
      close(94)
    end do

    print*
    print*, "Running Capt'n Oper..."
    
   call captn_init_oper()
   num_isotopes = 16
   jx = 0.5
   do cpl=1, 14
     filename = "oper_"//trim(cplConsts(cpl))//"_model.dat"
     open(55,file=filename)
     write(55,*) "DM_Mass | ", "  Captures | ", "  MaxCaptures"

     if (cpl==1) then
       call populate_array(1.65d-8, cpl, 0)
     else if (cpl==2) then
       call populate_array(0.d0, cpl-1, 0)
       call populate_array(1.65d-8, cpl+1, 0)
     else
       call populate_array(0.d0, cpl, 0)
       call populate_array(1.65d-8, cpl+1, 0)
     endif
     
     print*
     print*, "Running coupling constant: ", cplConsts(cpl)
     do i = 1,10
       mx = 5.d0 + dble(i)/5.
       call captn_oper(mx,jx,num_isotopes,capped_sd)
       maxcapture = maxcap(mx)
       write(55,*) mx, capped_sd, maxcapture
       print*, "mx: ", mx, "capped: ", capped_sd, "max_capture:", maxcapture
     end do
     close(55)
   end do

   print*
   print*, "Done!"

END PROGRAM GENCAP
