!   Super Capt'n
!   Module to house everying specific to Super Capt'n
!   Neal Avis Kozar 2021-22
!   all units of distance: GeV^-1
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: GeV^-1
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.
!   We've swapped to all natural units internally


module supermod
    implicit none
    double precision, parameter :: hbar=6.582119569d-25 !GeV*s, taken from PDG 2020
    double precision, parameter :: hbarc=1.973269804d-14 !GeV*cm, taken from PDG 2020
    double precision, parameter :: c0=1.d0
    double precision, parameter :: c1=2.99792458d10 !cm*s^-1, taken from PDG 2020
    double precision, parameter :: mnuc=0.938272088 !GeV, taken from PDG 2020
    ! double precision, parameter :: AMU=0.931494102 !GeV, taken from PDG 2020, could use instead of the proton mass?
    double precision, parameter :: pi=3.141592653 ! taken from PDG 2020
    double precision, parameter :: year=3.15569251d7 !seconds, taken as Tropical Year PDG 2020
    double precision, parameter :: kg_SolarM=1.98841d30 !kg*SolarMass^-1
    double precision, parameter :: GeV_kg=5.60958860d26 !GeV*kg^-1
    double precision, parameter :: cm_parsec=3.08567758149d18 !cm*pc^-1
    double precision, parameter :: GeV_ergs=624.151 !GeV*ergs^-1
    double precision, parameter :: AtomicNumber_super(11) = (/ 1., 4., 12., 14., 16., 20., 23., 24., &
                                                                28., 32., 56. /) !the isotopes the catena paper uses
    character (len=4) :: isotopes_super(11) = [character(len=4) :: "H", "He4", "C12", "N14", "O16", "Ne20", "Na23", "Mg24", &
                                                                    "Si28", "S32", "Fe56"] !the isotopes in text form to match against the W functions
    double precision, parameter :: AtomicSpin_super(11) = (/ 0.5, & ! {}^{1}\text{H}
                                                            0., & ! {}^{4}\text{He}
                                                            0., & ! {}^{12}\text{C}
                                                            1., & ! {}^{14}\text{N}
                                                            0., & ! {}^{16}\text{O}
                                                            0., & ! {}^{20}\text{Ne}
                                                            1.5, & ! {}^{23}\text{Na}
                                                            0., & ! {}^{24}\text{Mg}
                                                            0., & ! {}^{28}\text{Si}
                                                            0., & ! {}^{32}\text{S}
                                                            0. & ! {}^{56}\text{Fe}
                                                            /) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    double precision, parameter :: MassFrac_super(11) = (/   0.493, & ! {}^{1}\text{H}
                                                            0.35, & ! {}^{4}\text{He}
                                                            0.015, & ! {}^{12}\text{C}
                                                            0.004, & ! {}^{14}\text{N}
                                                            0.1, & ! {}^{16}\text{O}
                                                            0.005, & ! {}^{20}\text{Ne}
                                                            0.0004, & ! {}^{23}\text{Na}
                                                            0.005, & ! {}^{24}\text{Mg}
                                                            0.02, & ! {}^{28}\text{Si}
                                                            0.005, & ! {}^{32}\text{S}
                                                            0.007 & ! {}^{56}\text{Fe}
                                                            /) ! the isotopic abundances from Chris' notes, sourced from SNe sims: https://arxiv.org/abs/astro-ph/0112478
    double precision :: coupling_Array(14,2)
    double precision :: W_array_super(8,11,2,2,7)
    double precision :: yConverse_array_super(11)

    double precision :: mdm, rhoX, Mej, ISM, Dist, Esn, Age, V_w, Mdot
    integer :: novaType

    contains

    ! This gets you the parameters used in Rshock and Vshock
    subroutine novaParameters(lam_FE, lam_ST, R_0, t_0)
    implicit none
    double precision :: lam_FE, Lam_ST, R_0, t_0
        ! from Chris' notes
        if ( novaType == 1 ) then
            lam_FE = 4./7.
            lam_ST = 2./5.
            R_0 = ((3*Mej) / (4*pi*ISM*1.27*mnuc))**(1./3.)
            t_0 = R_0**(7./4.) * ((Mej*ISM*1.27*mnuc) / (0.38*Esn**2))**(1./4.)  ! include missing units of c to get t_0 in sec
        else if ( novaType == 2 ) then
            lam_FE = 6./7.
            lam_ST = 2./3.
            R_0 = Mej*V_w/Mdot
            t_0 = R_0**(7./6.) * ( 9*Mdot/V_w * (18*Mej)**(2)/(40*Esn)**(3) )**(1./6.)
        else
            stop "novaType can only be 1 (type Ia) or 2 (type II)."
        end if
    end subroutine novaParameters

    ! This gives you the radius of the SNe shockwave front as a function of time
    ! Returns in units of GeV^{-1}
    function Rshock(t)
        implicit none
        double precision :: Rshock
        double precision, intent(in) :: t

        double precision :: lam_FE, lam_ST
        double precision :: R_0, t_0

        ! retrieve lam_FE, lam_ST, R_0, and t_0 parameters
        call novaParameters(lam_FE, lam_ST, R_0, t_0)

        Rshock = R_0 * ((t/t_0)**(-5.*lam_FE) + (t/t_0)**(-5.*lam_ST))**(-1./5.)
        ! Rshock = 1.

    end function Rshock

    ! This gives you the velocity of the SNe shockwave front as a function of time
    ! Returns in units of c
    function Vshock(t)
        implicit none
        double precision :: Vshock
        double precision, intent(in) :: t

        double precision :: lam_FE, lam_ST
        double precision :: R_0, t_0

        ! retrieve lam_FE, lam_ST, R_0, and t_0 parameters
        call novaParameters(lam_FE, lam_ST, R_0, t_0)

        Vshock = R_0/t_0 * (Rshock(t)/R_0)**6 * (lam_FE*(t/t_0)**(-5.*lam_FE-1.) + lam_ST*(t/t_0)**(-5.*lam_ST-1.))
        ! Vshock = 1.

    end function Vshock

    ! adapted the catpn GFFI so that the limits run from E=0 to E = 2 m_x * ( (m_i*V_s)/(m_i + m_x) )**2, and Ekin is now 1/2 m_i V_s^2, momentum is now p = m_i * V_s
    function GFFI_H_oper(V_s, mq)
        double precision :: GFFI_H_oper
        double precision :: V_s
        integer :: mq
        double precision :: Emax!, p,  mu

        ! p = mnuc * V_s
        ! mu = mdm/mnuc
        Emax = 2*mdm * ( (mnuc*V_s)/(mnuc + mdm) )**2

        if (mq >= 0) then
            ! GFFI_H_oper = (p**2 * mu)**mq * Ekin * 1./(1.+mq) * (Emax/Ekin)**(mq+1)
            GFFI_H_oper = (2.d0 * mdm)**mq * Emax**(mq+1)/(mq + 1.d0)
        else
            stop "You cannot pass a negative integer to GFFI"
        end if

    end function GFFI_H_oper
    
    function GFFI_A_oper(V_s,A,mq)
        double precision :: GFFI_A_oper
        double precision :: V_s, A
        integer :: mq
        double precision :: Emax, mN, Ei!, p, mu
        double precision :: dgamic
        
        ! p = mN * V_s
        ! mu = mdm / (mN)
        mN = A*mnuc
        Emax = 2*mdm * ( (mN*V_s)/(mN + mdm) )**2
        ! Ei = 1./4.d0/mN/264.114*(45.d0*A**(-1./3.)-25.d0*A**(-2./3.))
        Ei = 9.39010d-4/mdm * (45.d0*A**(-1./3.)-25.d0*A**(-2./3.)) ! this is hbarc^2 / mdm * (1/b)^2, pulling b from page 10 of [arxiv:1501.03729]

        if (mq == 0) then
            GFFI_A_oper = Ei * (1.d0 - exp(-Emax/Ei))
        else if (mq > 0) then
            GFFI_A_oper = (2.d0 * mdm)**mq * Ei**(1 + mq) * (gamma(1.d0 + dble(mq)) - dgamic(1.+dble(mq), Emax/Ei))
        else
            stop "You cannot pass a negative integer to GFFI"
        end if

    end function GFFI_A_oper

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
subroutine supercaptn(mx_in, jx_in, vel_in, niso, scattered, scatteringCheck)
    use supermod
    implicit none
    double precision, intent(in) :: mx_in, jx_in, vel_in
    integer, intent(in):: niso
    double precision, intent(out) :: scattered, scatteringCheck !this is the output

    integer eli, funcType, tau, taup, term_R, term_W, q_pow, w_pow ! loop indicies
    double precision :: a, J, j_chi, mu_T, vel, time, V_s, R_s    ! parameters
    double precision :: result, DsigmaDe, N_i, sigma_i    ! used to tally results

    ! parameters used in prefactor calculation
    integer :: q_functype, q_index
    double precision ::  RFuncConst, WFuncConst, prefactor_functype, prefactor_current
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
    double precision :: prefactor_array(niso,11,2)

    mdm = mx_in ! input in GeV
    j_chi = jx_in
    vel = vel_in * 1.d5/c1 ! convert km s^{-1} to cm s^{-1} to c
    if ( vel .gt. 1. ) stop "The dark matter velocity cannot exceed the speed of light"
    time = age - Dist/vel ! the time (t=0 at SNe detonation) when the DM scattered, in GeV^{-1}

    R_s = Rshock(time) ! given in GeV^{-1}
    V_s = Vshock(time) ! given in c

    if (time .lt. 0.d0) then
        scattered = 0.d0
        scatteringCheck = 0.d0
    else
        do eli = 1, niso
            do q_pow = 1, 11
                do w_pow = 1, 2
                    prefactor_array(eli,q_pow,w_pow) = 0.d0
                end do
            end do
        end do

        ! First I set the entries in prefactor_array(niso,11,2)
        ! These are the constants that mulitply the corresonding q^{2n} w^{2m} terms, n:[0->10], and m:[0->1]
        do eli=1, niso
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
                                            prefactor_array(eli,q_index+1,1) = prefactor_array(eli,q_index+1,1) &
                                                - prefactor_current * (c0**2/(4.*mu_T**2)) ! The Rfunctions are programmed with the 1/c0^2 in their v_perp^2 term (so I need to un-correct it for the- q^2/(2*mu_T)^2, and leave it be for the w^2/c^2)

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

        ! now with all the prefactors computed, any 0.d0 entries in prefactor_array means that we can skip that q^{2n} w^{2m} term
        scatteringCheck = 0.d0
        scattered = 0.d0
        do eli=1,niso

            a = AtomicNumber_super(eli)
            J = AtomicSpin_super(eli)

            sigma_i = 0.d0
            N_i = MassFrac_super(eli)*Mej / (4*pi*mnuc*a*R_s**2)

            result = 0.d0

            ! condition on the maximal energy the DM can get from scattering off the SNe as given by Chris
            ! 1/2 * m_A * V_s^2 * 4 * m_A*m_x/(m_A+m_x)^2
            if (vel .lt. (2*A*mnuc*V_s)/(A*mnuc + mdm)) then

                do w_pow=1,2

                    do q_pow=1,11
                        if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                            ! the energy in the integral is given by delta function: E = (mdm * w^2)/2.
                            ! the momentum transfer is defined using the energy of moving DM: E = q^2/(mdm*2)
                            ! gives: q = mdm * w
                            if ( eli .eq. 1 ) then ! we're doing Hydrogen which doesn't get an exponential in its W function fits
                                result = result + prefactor_array(eli,q_pow,w_pow) *(mdm*vel)**(2*(q_pow-1)) *(V_s)**(2*(w_pow-1))
                                sigma_i = sigma_i + prefactor_array(eli,q_pow,w_pow)*(V_s)**(2*(w_pow-1))*GFFI_H_oper(V_s,q_pow-1)!(mdm*vel)**(2*(q_pow-1)) STUFF FROM GFFI
                            else
                                result = result + prefactor_array(eli,q_pow,w_pow) *(mdm*vel)**(2*(q_pow-1)) *(V_s)**(2*(w_pow-1)) &
                                    * exp(-2*yConverse_array_super(eli)*(mdm*vel)**2)
                                sigma_i = sigma_i + prefactor_array(eli,q_pow,w_pow)*(V_s)**(2*(w_pow-1))*GFFI_A_oper(V_s,a,q_pow-1)!(mdm*vel)**(2*(q_pow-1)) * exp(-2*yConverse_array_super(eli)*(mdm*vel)**2) STUFF FROM GFFI
                            end if
                        end if
                    end do !q_pow
                end do !w_pow

                DsigmaDe = result * (2. * mdm*c0**2)/(V_s**2 * (2*J+1))
                sigma_i = sigma_i * (2. * mdm*c0**2)/(V_s**2 * (2*J+1))

                scattered = scattered + (Mej*MassFrac_super(eli))/(a*mnuc) * DsigmaDe
                if ( eli.eq.1 ) then
                    scattered = scattered + (R_s*Mdot)/(V_w*mnuc) * DsigmaDe
                end if

                scatteringCheck = scatteringCheck + N_i*sigma_i

            end if
        end do !eli

        scattered = scattered * (rhoX*V_s*vel)/(4.*pi*Dist**2) !natural units
        if ( scatteringCheck > 1.d0 ) then
            scattered = scattered / (2*scatteringCheck)
            ! print*, "Scattering opacity check is greater than 1.0:", scatteringCheck
        end if
        scattered = scattered/hbarc**3 !recover units

    end if !End condition on t > 0

end subroutine supercaptn

subroutine supercaptn_init(rhoX_in, Mej_in, ISM_in, Dist_in, Esn_in, Age_in, novaTypeChosen, V_w_in, Mdot_in)
    ! input velocities in km/s, not cm/s!!!
    use supermod
    implicit none
    double precision, intent(in) :: rhoX_in, Mej_in, ISM_in, Dist_in, Esn_in, Age_in, V_w_in, Mdot_in
    integer, intent(in) :: novaTypeChosen

    integer :: i, j, k, l, m
    character (len=2) :: terms(7) = [character(len=2) :: "y0", "y1", "y2", "y3", "y4", "y5", "y6"]
    real :: WM, WS2, WS1, WP2, WMP2, WP1, WD, WS1D

    ! load values into module
    rhoX = rhoX_in*hbarc**3 ! loaded in GeV cm^{-3} -> GeV^{4}
    Mej = Mej_in * kg_SolarM * GeV_kg ! loaded in Solar Masses -> kg -> GeV
    ISM = ISM_in*hbarc**3 ! loaded in cm^{-3} -> GeV^{3}
    Dist = Dist_in * cm_parsec / hbarc ! loaded in pc -> cm -> GeV^{-1}
    Age = Age_in*year/hbar ! loaded in years -> seconds -> GeV^{-1}
    Esn = Esn_in * GeV_ergs ! loaded in ergs -> GeV
    V_w = V_w_in * 100000 / c1 ! loaded in km/s -> cm/s -> unitless
    Mdot = Mdot_in * kg_SolarM * GeV_kg / year * hbar ! loaded in MSun/yr -> kg/yr -> GeV/yr -> GeV/s -> GeV^{2}

    novaType = novaTypeChosen
    ! mass fraction is given by MassFrac_super, from SNe sims
    ! it doesn't vary with a radial coordinate, so it is only a fixed frac per isotope

    ! this array stores each of the constants of the W polynomials from paper arxiv:1501.03729's appendix individually
    ! array index m handles the 8 variants of the W functions in order [M, S", S', P", MP", P', Delta, S'Delta]
    ! index i handles the 11 isotopes [H, He4, C12, N14, O16, Ne20, Na23, Mg24, Si28, S32, Fe56]
    ! index j & k handle the two superscripts for each W function, each taking values of 0 and 1
    ! index L determines the power of each constant ranging from y^0 to y^6
    do m=1,8
        do i=1,11
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
    do i = 1, 11
        yConverse_array_super(i) = 264.114/(45.d0*AtomicNumber_super(i)**(-1./3.)-25.d0*AtomicNumber_super(i)**(-2./3.))
    end do
end subroutine supercaptn_init
