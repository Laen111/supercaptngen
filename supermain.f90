! Super Capt'n General testing program
! Uses `supercap.f90` source file
! call `supercaptn_init` to start
! then `populate_array_super` to populate the coupling array with desired parameters
! then `supercaptn` to perform the calculation
! can then call `populate_array_super` again to pick a new combination of couplings

    PROGRAM SUPERCAP
    implicit none
    character*300 :: filename ! writing filename
    integer :: cpl, i, j      ! loop indicies
    double precision, parameter :: pi=3.141592653, hbarc = 1.97d-14
    double precision :: dm_Density, ejecta_Mass, ISM_Density, dist_SN, energy_SN,age_SN     ! initialization parameters

    double precision :: dm_Mass, dm_Spin, dm_Vel ! WIMP dark matter particle parameters
    integer :: num_isotopes   ! number of isotopes summed over (8 is all of them)
    double precision :: coupleVal ! the value I will be assigning to the couplings when tested
    double precision :: dm_Scattered    ! (the output) number of scatters calculated at a given mass, spin, and velocity

    character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
                        "c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]
                        !   Strings of the (isoscalar) coupling names, just for writing convenience

    num_isotopes = 8  ! number of isotopes super capt'n will loop over in the calculation: up to 8 isotopes
    dm_Spin = 0.5     ! WIMP dark matter spin
    ! coupleVal = 1.65d-8
    ! coupleVal = sqrt(3.14*1.d-30)/(1.*1.973267d-44) ! get a coupling value in GeV^-2 for a specified xSection in cm^2
    coupleVal = sqrt(pi*1.d-30/(1./.4))/hbarc

    print*
    print*, "Initializing Super Capt'n..."

    dm_Density = 0.3 ! GeV cm^{-3}
    ejecta_Mass = 5. ! m_Sun
    ISM_Density = 3.73d-3 ! cm^{-3}
    dist_SN = 300. ! pc
    energy_SN = 8.d50 ! erg
    age_SN = 6.8d4 ! years
    call supercaptn_init(dm_Density, ejecta_Mass, ISM_Density, dist_SN, energy_SN, age_SN)

    print*
    print*, "Running Super Capt'n..."

   do cpl=1, 14
     ! one unique filename for each coupling constant
     filename = "Oper_"//trim(cplConsts(cpl))//"_Phi.dat"
     open(55,file=filename)
     write(55,*) "Coupling Constant value: ", coupleVal, "GeV^-2, or approximately", (coupleVal*hbarc)**2 * (1./.4)/pi, "cm^-2"
     write(55,*) "Dark Matter Mass ", " | ", " Dark Matter Velocity ", " | ", " Dark Matter Scattered"

     ! set the one coupling 'cpl' to a default value, skipping c2-0 (the 2nd index is not used)
     ! zero out the old coupling so that only the desired coupling is non-zero
     if (cpl==1) then
       call populate_array_super(coupleVal, 1, 0)
     else if (cpl==2) then
       call populate_array_super(0.d0, 1, 0)
       call populate_array_super(coupleVal, 3, 0)
     else
       call populate_array_super(0.d0, cpl, 0)
       call populate_array_super(coupleVal, cpl+1, 0)
     endif

     print*
     print*, "Running coupling constant: ", cplConsts(cpl)
     do i = 1,1!1
       dm_Mass = 1.d0!10**(dble(i-1)/5.)
       do j = 1,1001
        ! dm_Vel = 10**(dble(j-1)/10. + 1.) ! chris' notes test the range of velocities from 10^8 to 10^9 cm s^-1 for SI xSec
        ! dm_Vel = 4.32d3*(dble(j)*0.0005 + 1.)
        ! dm_Vel = 1.d3 * (dble(j-1)/100. + 4.)
        dm_Vel = 4320. + dble(j-1) * (4495.-4320.)/1000.
        call supercaptn(dm_Mass, dm_Spin, dm_Vel, num_isotopes, dm_Scattered)
        write(55,*) dm_Mass, dm_Vel, dm_Scattered
        ! Use this to check for negative scattering numbers (couplings 5, 7, 8, 13, 14)
        ! if ( dm_Scattered < 0. ) then
        !   print*, "Dark Matter Mass: ", dm_Mass, "Dark Matter Velocity: ", dm_Vel, "Scattered: ", dm_Scattered
        ! end if
        print*, "Dark Matter Mass: ", dm_Mass, "Dark Matter Velocity: ", dm_Vel, "Scattered: ", dm_Scattered
       end do
     end do
     close(55)
   end do

   print*
   print*, "Done!"

END PROGRAM SUPERCAP
