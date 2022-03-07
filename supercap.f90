!   Super Capt'n
!   Module to house everying specific to Super Capt'n
!   Neal Avis Kozar 2021
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.


module supermod
    implicit none
    double precision, parameter :: hbar=6.582d-25 !GeV*s
    double precision, parameter :: hbarc=1.97d-14 !GeV*cm
    double precision, parameter :: c0=1.d0 !
    double precision, parameter :: c1=2.99792458d10 !cm*s^-1
    double precision, parameter :: mnuc=0.938 !GeV
    double precision, parameter :: pi=3.141592653
    double precision, parameter  :: year= 3.154d7 !seconds
    double precision, parameter :: AtomicNumber_super(8) = (/ 1., 4., 12., 20., 24., 28., 32., 56. /) !the isotopes the catena paper uses
    character (len=4) :: isotopes_super(8) = [character(len=4) :: "H", "He4", "C12", "Ne20", "Mg24", "Si28", "S32", "Fe56"] !the isotopes in text form to match against the W functions
    double precision, parameter :: AtomicSpin_super(8) = (/ 0.5, & ! {}^{1}\text{H}
                                                            0., & ! {}^{4}\text{He}
                                                            0., & ! {}^{12}\text{C}
                                                            0., & ! {}^{20}\text{Ne}
                                                            0., & ! {}^{24}\text{Mg}
                                                            0., & ! {}^{28}\text{Si}
                                                            0., & ! {}^{32}\text{S}
                                                            0. & ! {}^{56}\text{Fe}
                                                            /) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    double precision, parameter :: MassFrac_super(8) = (/   0.493, & ! {}^{1}\text{H}
                                                            0.35, & ! {}^{4}\text{He}
                                                            0.015, & ! {}^{12}\text{C}
                                                            0.005, & ! {}^{20}\text{Ne}
                                                            0.005, & ! {}^{24}\text{Mg}
                                                            0.02, & ! {}^{28}\text{Si}
                                                            0.005, & ! {}^{32}\text{S}
                                                            0.007 & ! {}^{56}\text{Fe}
                                                            /) ! the isotopic abundances from Chris' notes, sourced from SNe sims: https://arxiv.org/abs/astro-ph/0112478
    double precision :: coupling_Array(14,2)
    double precision :: W_array_super(8,8,2,2,7)
    double precision :: yConverse_array_super(8)

    double precision :: mdm, rhoX, Mej, ISM, Dist, Esn,Age

    contains

    subroutine novaParameters(lam_FE, lam_ST, R_0, t_0)
    implicit none
    double precision :: lam_FE, Lam_ST, R_0, t_0
        ! from Chris' notes
        lam_FE = 4./7.
        lam_ST = 2./5.
        R_0 = ((3.*Mej) / (4*pi*ISM*1.27*mnuc))**(1./3.)
        t_0 = R_0**(7./4.) * ((Mej*ISM*1.27*mnuc) / (0.38*Esn**2))**(1./4.)  ! include missing units of c to get t_0 in sec

    end subroutine novaParameters

    ! This gives you the radius of the SNe shockwave front as a function of time
    ! NEEDS TO BE OUTPUT AS cm!
    function Rshock(t)
        implicit none
        double precision :: Rshock
        double precision, intent(in) :: t

        double precision :: lam_FE, lam_ST
        double precision :: R_0, t_0

        ! ! from Chris' notes
        ! lam_FE = 4./7.
        ! lam_ST = 2./5.
        ! R_0 = ((3*Mej) / (4*pi*ISM*1.27*mnuc))**(1./3.)
        ! t_0 = R_0**(7./4.) * ((Mej*ISM*1.27*mnuc) / (0.38*Esn**2))**(1./4.) / c0 ! include missing units of c to get t_0 in sec

        ! retrieve lam_FE, lam_ST, R_0, and t_0 parameters
        call novaParameters(lam_FE, lam_ST, R_0, t_0)

        Rshock = R_0 * ((t/t_0)**(-5.*lam_FE) + (t/t_0)**(-5.*lam_ST))**(-1./5.)
        ! Rshock = 1.

    end function Rshock

    ! This gives you the velocity of the SNe shockwave front as a function of time
    ! NEEDS TO BE OUTPUT AS cm s^{-1}!
    function Vshock(t)
        implicit none
        double precision :: Vshock
        double precision, intent(in) :: t

        double precision :: lam_FE, lam_ST
        double precision :: R_0, t_0

        ! ! from Chris' notes
        ! lam_FE = 4./7.
        ! lam_ST = 2./5.
        ! R_0 = ((3*Mej) / (4*pi*ISM*1.27*mnuc))**(1./3.)
        ! t_0 = R_0**(7./4.) * ((Mej*ISM*1.27*mnuc) / (0.38*Esn**2))**(1./4.) / c0 ! include missing units of c to get t_0 in sec

        ! retrieve lam_FE, lam_ST, R_0, and t_0 parameters
        call novaParameters(lam_FE, lam_ST, R_0, t_0)

        Vshock = R_0/t_0 * (Rshock(t)/R_0)**6 * (lam_FE*(t/t_0)**(-5.*lam_FE-1.) + lam_ST*(t/t_0)**(-5.*lam_ST-1.))
        ! Vshock = 1.

    end function Vshock

end module supermod


subroutine populate_array_super(val, couple, isospin)
    ! in the arxiv:1501.03729 paper, the non-zero values chosen were 1.65*10^-8 (represented as 1.65d-8 in the code)
    ! I was trying to directly edit 'couple' and 'isospin' to use in the array indices, but Fortran was throwing segfaults when doing this
    use supermod
    implicit none
    double precision, intent(in) :: val
    integer, intent(in) :: couple, isospin

    integer :: cpl, iso

    ! isospin can be 0 or 1
    if ((-1.lt.isospin).and.(isospin.lt.2)) then
        iso = isospin + 1 !fortran arrays start at 1
    else
        stop "Error: isospin can only be 0 or 1!"
    endif

    ! couple can be integer from 1 to 15, BUT 2 IS NOT ALLOWED!
    if (couple.lt.1) then
        stop "Error: you cannot pick a coupling constant lower than one!"
    else if (couple.eq.1) then
        cpl = couple
    else if (couple.eq.2) then
        stop "Error: you cannot use the second coupling constant!"
    else if (couple.gt.2) then
        cpl = couple - 1 !the coupling array doesn't have a slot for 2, so all constants other than one are shifted in row number
    else if (couple.gt.15) then
        stop "Error: you cannot pick a coupling constant past 15!"
    endif

    ! val is the value you want to populate with
    ! set the value picked in the slot chosen
    coupling_Array(cpl,iso) = val
end subroutine populate_array_super

! designed to calculate the scattering from supernovae on dark matter
! mx_in: dark matter mass in GeV
! jx_in: dark matter spin
! vel_in: dark matter final velocity in km s^{-1}
!
subroutine supercaptn(mx_in, jx_in, vel_in, niso, scattered)
    use supermod
    implicit none
    double precision, intent(in) :: mx_in, jx_in, vel_in
    integer, intent(in):: niso
    double precision, intent(out) :: scattered !this is the output

    integer eli, funcType, tau, taup, term_R, term_W, q_pow, w_pow ! loop indicies
    double precision :: a, J, j_chi, mu_T, vel, time, V_s, R_s    ! parameters
    double precision :: result, DsigmaDe    ! used to tally results

    ! parameters used in prefactor calculation
    integer :: q_functype, q_index
    double precision ::  RFuncConst, WFuncConst, prefactor_functype, prefactor_current
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
    double precision :: prefactor_array(niso,11,2)

    mdm = mx_in ! input in GeV
    j_chi = jx_in
    vel = vel_in * 1.d5/c1 ! convert km s^{-1} to cm s^{-1} to c
    time = age - Dist/vel ! time for DM to reach earth (traveling Dist to earth at upscattered velocity vel), in seconds
    R_s = Rshock(time) ! given in cm
    V_s = Vshock(time) ! given in cm s^{-1}

    write(*,*) "DM velocity =", vel, "Time =", time, "R_shock =", R_s, "V_shock =", V_s
    write(*,*) "DM kinetic energy =", 0.5*mdm*vel**2

    if (time .lt. 0.d0) then
      scattered = 0.d0
    else


    do eli = 1, niso
        do q_pow = 1, 11
            do w_pow = 1, 2
                prefactor_array(eli,q_pow,w_pow) = 0.d0
            end do
        end do
    end do

    ! First I set the entries in prefactor_array(niso,11,2)
    ! These are the constants that mulitply the corresonding integral evaluation
    do eli=1,niso !isotopeChosen, isotopeChosen
        ! I'll need mu_T to include in the prefactor when there is a v^2 term
        a = AtomicNumber_super(eli)
        mu_T = (mnuc*a*mdm)/(mnuc*a+mdm)

        ! the current response function type in order: M, S2, S1, P2, MP2, P1, D, S1D
        do funcType = 1,8

            ! contribution to q^2 count from sum over function types
            q_functype = 0
            prefactor_functype = 1.
            if ( functype.gt.3 ) then
                q_functype = 1
                prefactor_functype = 1./mnuc**2
            end if

            ! the first index on each response function
            do tau=1,2

                ! the second index on each response function
                do taup=1,2

                    ! the possible y-terms for each W function in order: y^0, y^1, y^2, y^3, y^4, y^5, y^6
                    do term_W = 1,7

                        WFuncConst = W_array_super(funcType,eli,tau,taup,term_W)

                        ! skip if the result gets multiplied by zero in the WFunction
                        if (WFuncConst.ne.0.) then

                            ! the possible terms for each R function in order: c, v2, q2, v2q2, q4, v2q4
                            do term_R = 1,6

                                ! pick out the appropriate term's constant from a given R function of tau, taup, and term_R
                                ! currently passes mnuc, and c0 - these are constants that could be shared to it through the shared module?
                                select case (funcType)
                                case (1)
                                    RFuncConst = RM(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array) !!!!!!!!!!!!!!! in the R functions the R term starts at zero, should change it to start at 1 like other Fortran things do for consistency
                                case (2)
                                    RFuncConst = RS2(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (3)
                                    RFuncConst = RS1(mnuc,c0,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (4)
                                    RFuncConst = RP2(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (5)
                                    RFuncConst = RMP2(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (6)
                                    RFuncConst = RP1(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (7)
                                    RFuncConst = RD(mnuc,tau,taup,term_R-1,j_chi,coupling_Array)
                                case (8)
                                    RFuncConst = RS1D(tau,taup,term_R-1,j_chi,coupling_Array)
                                case default
                                    stop "Um, I ran out of R functions to choose from?"
                                end select

                                ! skip if the result gets multiplied by zero in the RFunction
                                if (RFuncConst.ne.0.) then

                                    ! calculates the total number of q^2
                                    q_index = 1 + q_functype + term_W - 1 + floor((term_R-1.)/2.)
                                    prefactor_current = prefactor_functype * RFuncConst * WFuncConst * &
                                                           (c0* yConverse_array_super(eli))**(term_W-1)

                                    ! check if term_R is even (in my index convention this corresponds to it having a v^2 in the term)
                                    ! v^2 = w^2 - q^2/(2mu_T)^2
                                    if ( mod(term_R,2).eq.0 ) then
                                        ! this is the -q^2/(2mu_T)^2 contribution
                                        ! it has one extra q^2 contribution compared to the current W & R function contributions
                                        ! TRACK THIS LINE, SHOULD BE ONLY PLACE WHERE NEGATIVES APPEAR?
                                        prefactor_array(eli,q_index+1,1) = prefactor_array(eli,q_index+1,1) - prefactor_current * &
                                            (c0**2/(4.*mu_T**2)) ! The Rfunctions are programmed with the 1/c0^2 in their v_perp^2 term (so I need to un-correct it for the- q^2/(2*mu_T)^2, and leave it be for the w^2/c^2)



                                        ! this is the +w^2 contribution
                                        ! it has the same q^2 contribution, but has a v_perp^2 contribution
                                        prefactor_array(eli,q_index,2) = prefactor_array(eli,q_index,2) + prefactor_current

                                    else
                                        prefactor_array(eli,q_index,1) = prefactor_array(eli,q_index,1) + prefactor_current

                                    end if
                                end if
                            end do !term_R
                        end if
                    end do !term_W
                end do !taup
            end do !tau
        end do !functype
    end do !eli

    ! now with all the prefactors computed, any 0.d0 entries in prefactor_array means that we can skip that integral evaluation!
    scattered = 0.d0
    ! print*,'WARNING ONLY LOOKING AT H!!!!!'
    do eli=1,niso

        a = AtomicNumber_super(eli)
        J = AtomicSpin_super(eli)

        result = 0.d0

        write(*,*) "Element: ", isotopes_super(eli), " Maximum energy =", 2. * mdm * (A*mnuc*V_s/c0/(A*mnuc+mdm))**2

        ! condition on the maximal energy the DM can get from scattering off the SNe as given by Chris
        ! 1/2 * m_A * V_s^2 * 4 * m_A*m_x/(m_A+m_x)^2
        if (((0.5*mdm*vel**2 .lt. 2*A**2 * mnuc**2 * mdm/(A*mnuc+mdm)**2*V_s**2)) .and. (age .gt. Dist/vel)) then

            do w_pow=1,2

                do q_pow=1,11
                    if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                        ! the energy in the integral is given by delta function: E = (mdm * w^2)/2.
                        ! the momentum transfer is defined using the energy of moving DM: E = q^2/(mdm*2)
                        ! gives: q = mdm * w
                        ! means I can squeeze both the q and w powers onto w, and leave mdm just on q powers:
                        result = result + prefactor_array(eli,q_pow,w_pow) * mdm**(q_pow-1) * (V_s/c0)**(w_pow-1)*vel**(q_pow-1)
                    end if
                end do !q_pow
            end do !w_pow

            ! CHECK THIS FOR UNITS AGAIN, IT PROBABLY NEEDS TO BE CHANGED TO GET PHI(v) OUT OF IT
            DsigmaDe = result * (2. * mnuc*a*c0**2)/(V_s**2 * (2.*J+1.))
            if (eli .eq. 1) then
            open(22,file='dsigmade.dat',ACCESS="APPEND")
            write(22,*) V_s, DsigmaDe
            close(22)
            end if


            scattered = scattered + (Mej*MassFrac_super(eli))/(a*mnuc) * DsigmaDe
            if ( eli.eq.1 ) then
                scattered = scattered + 4./3.*pi*R_s**3*ISM*c0**2 * DsigmaDe
            end if
        end if
    end do !eli

    scattered = scattered * (rhoX*V_s*vel)/(4.*pi*Dist**2) !natural units
    scattered = scattered/hbarc**3 !recover units

end if !End condition on t > 0
    ! if (capped .gt. 1.d100) then
    !   print*,"Capt'n General says: Oh my, it looks like you are capturing an", &
    !     "infinite amount of dark matter in the Sun. Best to look into that."
    ! end if
end subroutine supercaptn

subroutine supercaptn_init(rhoX_in, Mej_in, ISM_in, Dist_in, Esn_in,Age_in)
    ! input velocities in km/s, not cm/s!!!
    use supermod
    implicit none
    double precision, intent(in) :: rhoX_in, Mej_in, ISM_in, Dist_in, Esn_in, Age_in

    integer :: i, j, k, l, m
    character (len=2) :: terms(7) = [character(len=2) :: "y0", "y1", "y2", "y3", "y4", "y5", "y6"]
    real :: WM, WS2, WS1, WP2, WMP2, WP1, WD, WS1D

    ! load values into module
    rhoX = rhoX_in*hbarc**3 ! loaded in GeV cm^{-3} -> GeV^4
    Mej = Mej_in * 1.98841d30 * 5.60958860d26! convert Solar Masses to kg to GeV c^{-2}
    ISM = ISM_in*hbarc**3 ! loaded in cm^{-3}
    Dist = Dist_in * 3.08567758149d18 ! convert pc to cm
    Dist = Dist/hbarc !GeV^-1
    Age = Age_in*year/hbar !GeV^-1
    ! this might need to be in GeV, check the Vshock and Rshock functions to see if units work out in ergs
    Esn = Esn_in * 624.151 ! convert ergs to GeV

    ! mass fraction is given by MassFrac_super, from SNe sims
    ! it doesn't vary with a radial coordinate, so it is only a fixed frac per isotope

    ! this array stores each of the constants of the W polynomials from paper arxiv:1501.03729's appendix individually
    ! array index m handles the 8 varients of the W functions in order [M, S", S', P", MP", P', Delta, S'Delta]
    ! index i handles the 8 isotopes [H, He4, C12, Ne20, Mg24, Si28, S32, Fe56]
    ! index j & k handle the two superscripts for each W function, each taking values of 0 and 1
    ! index L determines the power of each constant ranging from y^0 to y^6
    do m=1,8
        do i=1,8
            do j=1,2
                do k=1,2
                    do l=1,7
                        if (m.eq.1) then
                            W_array_super(m,i,j,k,l) = WM(j-1,k-1,isotopes_super(i),terms(l))
                        else if (m.eq.2) then
                            W_array_super(m,i,j,k,l) = WS2(j-1,k-1,isotopes_super(i),terms(l))
                        else if (m.eq.3) then
                            W_array_super(m,i,j,k,l) = WS1(j-1,k-1,isotopes_super(i),terms(l))
                        else if (m.eq.4) then
                            W_array_super(m,i,j,k,l) = WP2(j-1,k-1,isotopes_super(i),terms(l))
                        else if (m.eq.5) then
                            W_array_super(m,i,j,k,l) = WMP2(j-1,k-1,isotopes_super(i),terms(l))
                        else if (m.eq.6) then
                            W_array_super(m,i,j,k,l) = WP1(j-1,k-1,isotopes_super(i),terms(l))
                        else if (m.eq.7) then
                            W_array_super(m,i,j,k,l) = WD(j-1,k-1,isotopes_super(i),terms(l))
                        else
                            W_array_super(m,i,j,k,l) = WS1D(j-1,k-1,isotopes_super(i),terms(l))
                        end if
                    end do
                end do
            end do
        end do
    end do

    ! initiate the coupling_Array (full of the coupling constants) with all zeros
    ! populate_array will place the non-zero value into a chosen slot at runtime
    coupling_Array = reshape((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), (/14, 2/))

    ! yconv comes from arxiv:1501.03729 page 10, where yconv = (b/{2 hbar c})^2
    do i = 1, 8
        yConverse_array_super(i) = 264.114/(45.d0*AtomicNumber_super(i)**(-1./3.)-25.d0*AtomicNumber_super(i)**(-2./3.))
    end do
end subroutine supercaptn_init
