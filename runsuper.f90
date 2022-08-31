program RunSuper
    implicit none
    character*300 :: filename ! writing filename
    integer :: i ! loop indicies
    double precision, parameter :: pi=3.141592653, hbarc = 1.973269804d-14, mnuc=0.938272088
    double precision :: dm_Density, ejecta_Mass, ISM_Density, dist_SN, energy_SN,age_SN, stellar_wind_vel, stellar_mass_loss ! initialization parameters
    integer :: novaTypeSelection

    double precision :: dm_Mass, dm_Spin, dm_Vel ! WIMP dark matter particle parameters
    integer :: num_isotopes=11   ! number of isotopes summed over (11 is all of them)
    double precision :: coupleValArray(15,2) ! the value I will be assigning to the couplings when tested
    double precision :: mu ! reduced dm-nucleon mass
    double precision :: velInit, velFinal ! velocity range parameters
    integer :: velNum ! velocity range parameter
    double precision :: dm_Scattered, scatteringOpacityCond    ! (the output) number of scatters calculated at a given mass, spin, and velocity
    character (len=2) :: int_str


! Filename
    filename = "test.dat" ! the file you want to write to
! Velocity range
    velInit = 4300.00    ! starting velocity [km s^-1]
    velFinal = 4500.00   ! ending velocity [km s^-1]
    velNum = 10000       ! integer number of velocity steps
!-----------------------------------------------------------------------------------------------------------------------------------
! Dark matter parameters
    dm_Mass = 1.d0 ! dark matter mass
    dm_Spin = 0.5 ! WIMP dark matter spin

    mu = (dm_Mass*mnuc)/(dm_Mass+mnuc)  ! the DM-nucleon reduced mass, used for effective cross section approximation
    ! sqrt(4.*pi*XSec)/(mu*hbarc) gives you the coupling from an effective cross section in cm^2
    ! the isoscalar couplings you want to run
    coupleValArray(1,1) = sqrt(4.*pi*1.d-30)/(mu*hbarc) ! Spin Independent
    ! c2_0 is not allowed
    coupleValArray(3,1) = 0.d0
    coupleValArray(4,1) = 0.d0 ! Spin Dependent
    coupleValArray(5,1) = 0.d0
    coupleValArray(6,1) = 0.d0
    coupleValArray(7,1) = 0.d0
    coupleValArray(8,1) = 0.d0
    coupleValArray(9,1) = 0.d0
    coupleValArray(10,1) = 0.d0
    coupleValArray(11,1) = 0.d0
    coupleValArray(12,1) = 0.d0
    coupleValArray(13,1) = 0.d0
    coupleValArray(14,1) = 0.d0
    coupleValArray(15,1) = 0.d0

    ! the isovector couplings you want to run
    coupleValArray(1,2) = 0.d0
    ! c2_1 is not allowed
    coupleValArray(3,2) = 0.d0
    coupleValArray(4,2) = 0.d0
    coupleValArray(5,2) = 0.d0
    coupleValArray(6,2) = 0.d0
    coupleValArray(7,2) = 0.d0
    coupleValArray(8,2) = 0.d0
    coupleValArray(9,2) = 0.d0
    coupleValArray(10,2) = 0.d0
    coupleValArray(11,2) = 0.d0
    coupleValArray(12,2) = 0.d0
    coupleValArray(13,2) = 0.d0
    coupleValArray(14,2) = 0.d0
    coupleValArray(15,2) = 0.d0
    
!-----------------------------------------------------------------------------------------------------------------------------------
! Supernova initialization parameters
    dm_Density = 0.3 ! GeV cm^{-3}
    ejecta_Mass = 5. ! m_Sun
    ISM_Density = 3.73d-3 ! cm^{-3}
    dist_SN = 300. ! pc
    energy_SN = 8.d50 ! erg
    age_SN = 6.8d4 ! years
    novaTypeSelection = 2 ! 1:type Ia (old) 2:type II (new)
    stellar_wind_vel = 10. ! km s^{-1}
    stellar_mass_loss = 1.d-5 ! m_Sun yr^{-1}
    call supercaptn_init(dm_Density, ejecta_Mass, ISM_Density, dist_SN, energy_SN, age_SN, &
                            novaTypeSelection, stellar_wind_vel, stellar_mass_loss)

!-----------------------------------------------------------------------------------------------------------------------------------
! Main script follows
    ! second coupling never used, just set it to negative to make it obvious
    coupleValArray(2,1) = -1.
    coupleValArray(2,2) = -1.
    
    ! supercaptn_init zeros out the internal coupling array so that only the desired coupling(s) is(are) non-zero
    ! skipping c2-0 (the 2nd index is not used)
    do i = 1, 15
        if ( i.ne.2 ) then
            call populate_array_super(coupleValArray(i,1), i, 0)
            call populate_array_super(coupleValArray(i,2), i, 1)
        end if
    end do

    open(55,file=filename)
    ! writing header with all the values used
    write(55, *) "Velocity spectrum at Earth, using the NREO formalism [arXiv:1501.03729]"
    write(55, *) "Dark matter density [GeV cm^-3]: ", dm_Density, &
                    " | ISM density [cm^-3]: ", ISM_Density, &
                    " | SN ejecta mass [M_Sun]: ", ejecta_Mass, &
                    " | SN distance [pc]: ", dist_SN, &
                    " | SN energy [erg]: ", energy_SN, &
                    " | SN age [tropical years]: ", age_SN
    write(55,*) "Dark matter mass [GeV]: ", dm_Mass, &
                " | Dark matter spin [unitless]: ", dm_Spin
    write(55,*) "The NREO couplings in use: "
    do i = 1, 15
        if ( i.ne.2 ) then
            write(int_str,"(I0.2)") i
            write(55,*) "c"//trim(int_str)//"_0 [GeV^-2]: ", coupleValArray(i,1), &
                            "or approximately [cm^2]: ", (coupleValArray(i,1)*mu*hbarc)**2 / (4.*pi), &
                        " | c"//trim(int_str)//"_1 [GeV^-2]: ", coupleValArray(i,2), &
                            "or approximately [cm^2]: ", (coupleValArray(i,2)*mu*hbarc)**2 / (4.*pi)
        end if
    end do

    write(55,*)! writing a blank line
    write(55,*) "To extract the data into python numpy arrays using one line: ", &
                "velArray, scatterArray, opacityCond = numpy.loadtxt(filename, skiprows=21).T"
    write(55,*) "Dark Matter Velocity [km s^-1]", " | ", &
        " Dark Matter Scattered [(cm s^-1)^-1 (s)^-1 (cm^2)^-1]", " | ", &
        " Scattering Opacity Condition [unitless]"

    ! calculate all the scatterings and write to the file
    print*, "Starting calculation..."
    do i = 1,velNum+1   ! create velocity velNum points linearly spaced on the range velInit to velFinal
        dm_Vel = velInit + dble(i-1) * (velFinal-velInit)/dble(velNum)
        call supercaptn(dm_Mass, dm_Spin, dm_Vel, num_isotopes, dm_Scattered, scatteringOpacityCond)
        write(55,*) dm_Vel, dm_Scattered, scatteringOpacityCond
    end do

    close(55)
    print*, "Done!"
end program RunSuper