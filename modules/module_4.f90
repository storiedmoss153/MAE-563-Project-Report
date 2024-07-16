module module_4
    use universal_module
    use module_3
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    contains

    function M4_test_mach(p1,pt3,eta_n) result(M_prime)
        ! Test Mach number, M' calcualtion
        implicit none

        real(dp), intent(in) :: p1, pt3, eta_n
        real(dp) :: var1, M_prime

        var1 = 1 - (p1 / pt3)**((gamma_3 - 1) / gamma_3)
        M_prime = (2 * eta_n * var1 / ((gamma_3 - 1) * (1 - eta_n * var1)))**0.5

    end function M4_test_mach

    function M4_nozzle_choke_test(Tt3,p1,pt3,eta_n) result(nozzle_conditions)
        ! Test for nozzle choking, outputs total temperature, static pressure,
        ! and Mach at state e
        implicit none

        real(dp), intent(in) :: Tt3, p1, pt3, eta_n
        real(dp) :: Me, pe, Tte, M_prime, nozzle_conditions(3)

        Tte = Tt3
        M_prime = M4_test_mach(p1,pt3,eta_n)

        if (M_prime < 1) then
            Me = M_prime
            pe = p1
        else
            Me = 1
            pe = pt3 * (1 - ((gamma_3 - 1) / (gamma_3 + 1)) / eta_n)**& 
                       (gamma_3 / (gamma_3 - 1))
        end if

        nozzle_conditions = [Me,pe,Tte]
    end function M4_nozzle_choke_test

    function M4_total_to_static_relations(Tte,pe,Me) result(Temp_n_pres)
        ! Total to static relations for state e
        implicit none

        real(dp), intent(in) :: Tte, pe, Me
        real(dp) :: Tte_Te, pte_pe, Temp_n_pres(2)

        Tte_Te = 1 + 0.5 * (gamma_3 - 1) * Me**2
        pte_pe = Tte_Te**(gamma_3 / (gamma_3 - 1))

        Temp_n_pres = [Tte / Tte_Te, pe * pte_pe]
    end function M4_total_to_static_relations

    function M4_exit_mass_flux(Te,pe,Ve,Area_e) result(me_dot)
        ! Exit mass flux from nozzle
        implicit none

        real(dp), intent(in) :: Te, pe, Ve, Area_e
        real(dp) :: me_dot

        me_dot = pe * Ve * Area_e / (R_3 * Te)
    end function M4_exit_mass_flux

    function M4_all_function(Tt3,p1,pt3,eta_n,delta_s13,Area_e) &
             result(M4_output_matrix)
        ! Module 4 function that combines all previous functions
        ! and outputs required data for state e
        implicit none

        real(dp), intent(in) :: Tt3, p1, pt3, eta_n, delta_s13, Area_e
        real(dp) :: Me, pe, pte, Te, Tte, ae, cpe, Ve, me_dot
        real(dp) :: delta_s3e, delta_s1e
        real(dp) :: nozzle_choke(3), T2S_relations(2)
        real(dp) :: M4_output_matrix(11)

        nozzle_choke = M4_nozzle_choke_test(Tt3,p1,pt3,eta_n)
        Me = nozzle_choke(1)
        pe = nozzle_choke(2)
        Tte = nozzle_choke(3)

        T2S_relations = M4_total_to_static_relations(Tte,pe,Me)
        Te = T2S_relations(1)
        pte = T2S_relations(2)

        cpe = M3_pressure_coefficient(Te)
        ae = M1_speed_of_sound(Te,gamma_3,R_3)
        Ve = M1_flow_speed(Me,ae)
        me_dot = M4_exit_mass_flux(Te,pe,Ve,Area_e) * 1e3

        delta_s3e = M3_entropy_increase(Tt3,Tte,pt3,pte,cpe)
        delta_s1e = delta_s13 + delta_s3e

        M4_output_matrix = [Me,pte,pe,Tte,Te,cpe,ae,Ve,me_dot, &
                            delta_s3e,delta_s1e]
    end function M4_all_function
end module module_4