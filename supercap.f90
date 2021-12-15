!   Capt'n Oper
!   Module to house everying specific to captn operator
!   Neal Avis Kozar 2020
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.


module supermod
    implicit none
    double precision, parameter :: hbar=6.582d-25 !GeV*s
    double precision, parameter :: pi=3.141592653, NAvo=6.0221409d23
    double precision, parameter :: c0=2.99792458d10, mnuc=0.938
    
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

    integer :: q_shared
    logical :: w_shared

    double precision :: mdm, vesc_shared, a_shared, mu, muplus
    
    ! values shared by module, set by supercaptn_init
    double precision :: usun, u0, rho0, vesc_halo
    double precision :: Mej, ISM, Dist, Esn

    contains

    ! !   this is the function f_sun(u) in 1504.04378 eqn 2.2 divided by u
    ! !velocity distribution,
    ! function vdist_over_u(u)
    !     double precision :: u, vdist_over_u, normfact
    !     vdist_over_u = (3./2.)**(3./2.)*4.*rho0*u/sqrt(pi)/mdm/u0**3 &
    !     *exp(-3.*(usun**2+u**2)/(2.*u0**2))*sinh(3.*u*usun/u0**2)/(3.*u*usun/u0**2)
    !     !normfact = .5*erf(sqrt(3./2.)*(vesc_halo-usun)/u0) + &
    !     !.5*erf(sqrt(3./2.)*(vesc_halo+usun)/u0)+ u0/(sqrt(6.*pi)*usun) &
    !     !*(exp(-3.*(usun+vesc_halo)/2./u0**2)-exp(-3.*(usun-vesc_halo)/2./u0**2))
    !     normfact = 1.
    !     !print*,normfact
    !     vdist_over_u = vdist_over_u/normfact
    ! end function vdist_over_u

    ! ! having removed the scaling momentum, are the units off here? I'm looking at the p/c0 in particular
    ! function GFFI_H_oper(w,vesc,mq)
    !     double precision :: p, mu,w,vesc,u,muplus,GFFI_H_oper,G
    !     integer mq
    !     p = mdm*w
    !     mu = mdm/mnuc
    !     muplus = (1.+mu)/2.
    !     u = sqrt(w**2-vesc**2)
    !     if (mq .ne. -1) then
    !         G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*1./(1.+mq)*((mu/muplus**2)**(mq+1)-(u**2/w**2)**(mq+1))
    !     else
    !         G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*log(mu/muplus**2*w**2/(u)**2)
    !     endif
    !     GFFI_H_oper = G
    ! end function GFFI_H_oper
    
    ! function GFFI_A_oper(w,vesc,A,mq)
    !     double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
    !     double precision :: dgamic,GFFI_A_oper
    !     integer :: mq
    !     p = mdm*w
    !     mu = mdm/mnuc/A
    !     muplus = (1.+mu)/2.
    !     u = sqrt(w**2-vesc**2)
    !     mN = A*mnuc
    !     Ei = 1./4.d0/mN/264.114*(45.d0*A**(-1./3.)-25.d0*A**(-2./3.))
    !     B = .5*mdm*w**2/Ei/c0**2
    !     if (mq .eq. 0) then
    !         GFFI_A_oper = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
    !     else
    !         GFFI_A_oper = ((p)/c0)**(2*mq)*Ei*c0**2/(B*mu)**mq*(dgamic(1.+dble(mq),B*u**2/w**2) &
    !             - dgamic(1.+dble(mq),B*mu/muplus**2))
    !     end if
    ! end function GFFI_A_oper

    ! This gives you the radius of the SNe shockwave front as a function of time
    ! What are the units?
    function Rshock(t)
        implicit none
        double precision :: Rshock
        double precision, intent(in) :: t
        double precision :: lam_FE, lam_ST
        double precision :: R_0, t_0
        
        ! from Chris' notes
        lam_FE = 4./7.
        lam_ST = 2./5.
        R_0 = ((3*Mej) / (4*pi*ISM*1.27*mnuc))**(1./3.)
        t_0 = R_0**(7./4.) * ((Mej*ISM*1.27*mnuc) / (0.38*Esn**2))**(1./4.)

        ! Rshock = R_0 * ((t/t_0)**(-5.*lam_FE) + (t/t_0)**(-5.*lam_ST))**(-1./5.)
        Rshock = 1.

    end function Rshock

    ! This gives you the velocity of the SNe shockwave front as a function of time
    ! What are the units?
    function Vshock(t)
        implicit none
        double precision :: Vshock
        double precision, intent(in) :: t
        double precision :: lam_FE, lam_ST
        double precision :: R_0, t_0
        
        ! from Chris' notes
        lam_FE = 4./7.
        lam_ST = 2./5.
        R_0 = ((3*Mej) / (4*pi*ISM*1.27*mnuc))**(1./3.)
        t_0 = R_0**(7./4.) * ((Mej*ISM*1.27*mnuc) / (0.38*Esn**2))**(1./4.)

        ! Vshock = R_0/t_0 * (Rshock(t)/R_0)**6 * (lam_FE*(t/t_0)**(-5.*lam_FE-1.) + lam_ST*(t/t_0)**(-5.*lam_ST-1.))
        Vshock = 1.
    
    end function Vshock

end module supermod


subroutine populate_array_super(val, couple, isospin)
    ! in the arxiv:1501.03729 paper, the non-zero values chosen were 1.65*10^-8 (represented as 1.65d-8 in the code)
    ! I was trying to directly edit 'couple' and 'isospin' to use in the array indices, but Fortran was throwing segfaults when doing this
    use supermod
    implicit none
    integer :: couple, isospin
    double precision :: val
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

! ! this is the integral over R in eqn 2.3 in arxiv:1501.03729
! ! note that Omega there is expanded and broken into terms of the form: const. * q^2n * exp{E_R/E_i}
! ! I've doen this so that I can tap into the GFFI functions in eqn 2.9 of arxiv:1504.04378
! ! THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
! function integrand_super(u, foveru)
!     use supermod
!     double precision :: u, w, integrand_oper, foveru !int, vesc
!     external foveru

!     ! vesc = tab_vesc(ri_for_omega)

!     w = sqrt(u**2+vesc_shared_arr(rindex_shared)**2)

!     !Switch depending on whether we are capturing on Hydrogen or not
!     if (a_shared .gt. 2.d0) then
!         integrand_oper = foveru(u)*GFFI_A_oper(w,vesc_shared_arr(rindex_shared),a_shared,q_shared)
!     else
!         integrand_oper = foveru(u)*GFFI_H_oper(w,vesc_shared_arr(rindex_shared),q_shared)
!     end if
!     if (w_shared) then
!         integrand_oper = integrand_oper * w**2
!     end if

! end function integrand_super

! designed to calculate the scattering from supernovae on dark matter
subroutine supercaptn(mx_in, jx_in, w, niso, scattered)
    use supermod
    implicit none
    double precision, intent(in) :: mx_in, jx_in, w
    integer, intent(in):: niso
    double precision, intent(out) :: scattered !this is the output

    integer eli, funcType, tau, taup, term_R, term_W, q_pow, w_pow ! loop indicies
    double precision :: a, J, j_chi, mu_T
    double precision :: result, DsigmaDe

    integer :: q_functype, q_index
    double precision ::  RFuncConst, WFuncConst, prefactor_functype, prefactor_current
    double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2 !R functions stored in their own source files
    double precision :: prefactor_array(niso,11,2)

    mdm = mx_in
    j_chi = jx_in

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
                                                            yConverse_array_super(eli)**(term_W-1)

                                    ! check if term_R is even (in my index convention this corresponds to it having a v^2 in the term)
                                    ! v^2 = w^2 - q^2/(2mu_T)^2
                                    if ( mod(term_R,2).eq.0 ) then
                                        ! this is the -q^2/(2mu_T)^2 contribution
                                        ! it has one extra q^2 contribution compared to the current W & R function contributions
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
    ! umin = 0.d0
    ! do ri=1,nlines
    !     vesc = tab_vesc(ri)
    !     rindex_shared = ri !make accessible via the module, used to get vesc value (and calculate w velocity) in integrand
    !     vesc_shared_arr(ri) = vesc !make accessible via the module

    do eli=1,niso !isotopeChosen, isotopeChosen
        ! int_res = 0.d0
        a = AtomicNumber_super(eli)
        ! a_shared = a !make accessible via the module

        ! mu = mdm/(mnuc*a)
        ! muplus = (1.+mu)/2.
        ! muminus = (mu-1.d0)/2.

        J = AtomicSpin_super(eli)

        ! DO WE STILL NEED TO DO THIS CHOP FOR SNe?
        ! Chop the top of the integral off at the smaller of the halo escape velocity or the minimum velocity required for capture.
        ! umax = min(vesc * sqrt(mu)/abs(muminus), vesc_halo)

        result = 0.d0
        do w_pow=1,2
            ! ! toggles whether we integrate with the w^2 term on
            ! w_shared = .false.
            ! if(w_pow.eq.2) then
            !     w_shared = .true.
            ! end if

            do q_pow=1,11
                if ( prefactor_array(eli,q_pow,w_pow).ne.0. ) then
                    ! q_shared = q_pow - 1
                    ! !Call integrator
                    ! call dsntdqagse(integrand_oper,vdist_over_u,umin,umax, &
                    !     epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

                    ! the energy in the integral is given by delta function: E = (mdm * w^2)/2.
                    ! the momentum transfer is defined using the energy of moving DM: E = q^2/(mdm*2)
                    ! gives: q = mdm * w
                    ! means I can squeeze both the q and w powers onto w, and leave mdm just on q powers:
                    result = result + prefactor_array(eli,q_pow,w_pow) * mdm**(q_pow-1) * w**(w_pow+q_pow-2)
                end if
            end do !q_pow
        end do !w_pow

        ! CHECK THIS FOR UNITS AGAIN, IT PROBABLY NEEDS TO BE CHANGED TO GET PHI(v) OUT OF IT
        ! factor_final = (2*mnuc*a)/(2*J+1) !* NAvo*tab_starrho(ri)*tab_mfr_oper(ri,eli)/(mnuc*a) * &
            ! tab_r(ri)**2*tab_dr(ri) * (hbar*c0)**2

        DsigmaDe = result * (2*mnuc*a)/(2*J+1)

        scattered = scattered + (Mej*MassFrac_super(eli))/(a*mnuc) * DsigmaDe
        if ( eli.eq.1 ) then
            scattered = scattered + (4.*pi*Rshock(Dist/w)*ISM)/3. * DsigmaDe
        end if

    end do !eli
    ! end do !ri

    scattered = scattered * (rho0*Vshock(Dist/w))/(4*pi*w*Dist**2)
    ! capped = 4.d0*pi*Rsun**3*capped

    ! if (capped .gt. 1.d100) then
    !   print*,"Capt'n General says: Oh my, it looks like you are capturing an", &
    !     "infinite amount of dark matter in the Sun. Best to look into that."
    ! end if
end subroutine supercaptn

subroutine supercaptn_init(rho0_in, usun_in, u0_in, vesc_in, Mej_in, ISM_in, Dist_in,Esn_in)
    ! input velocities in km/s, not cm/s!!!
    use supermod
    implicit none
    integer :: i, j, k, l, m
    character (len=2) :: terms(7) = [character(len=2) :: "y0", "y1", "y2", "y3", "y4", "y5", "y6"]
    real :: WM, WS2, WS1, WP2, WMP2, WP1, WD, WS1D
    
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in
    double precision, intent(in) :: Mej_in, ISM_in, Dist_in, Esn_in
    
    ! load values into module
    usun = usun_in*1.d5
    u0 =  u0_in*1.d5
    rho0 = rho0_in
    vesc_halo = vesc_in*1.d5
    Mej = Mej_in*1.98841d30 ! convert Solar Masses to kg
    ISM = ISM_in ! loaded in cm^-3
    Dist = Dist_in ! loaded in cm
    Esn = Esn_in !loaded in ergs

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

    do i = 1, 8
        yConverse_array_super(i) = 264.114/(45.d0*AtomicNumber_super(i)**(-1./3.)-25.d0*AtomicNumber_super(i)**(-2./3.))
    end do
end subroutine supercaptn_init
