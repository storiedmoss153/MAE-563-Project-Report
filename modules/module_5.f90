module module_5
    use universal_module
    use module_4
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    contains

    function M5_nozzle_external_efficiency(p1,pt3,eta_n) result(eta_n_ext)
        ! Nozzle external efficiency calculations
        implicit none

        real(dp), intent(in) :: p1,pt3,eta_n
        real(dp) :: M_prime, eta_n_ext

        M_prime = M4_test_mach(p1,pt3,eta_n)

        if (M_prime < 1) then
            eta_n_ext = 1
        else
            eta_n_ext = M_prime**(-0.3)
        end if
    end function M5_nozzle_external_efficiency

    function M5_total_to_static_relations(eta_n_ext,p1,pte,Tte) &
             result(output_array)
        ! Total to static relations for state 4
        implicit none

        real(dp), intent(in) :: eta_n_ext, p1, pte, Tte
        real(dp) :: Tt4, T4_Tt4, p4, pt4_p4, M4, output_array(5)

        Tt4 = Tte
        T4_Tt4 = 1 - eta_n_ext * (1 - (p1 / pte)**((gamma_3 - 1) / gamma_3))
        M4 = (2 * (T4_Tt4**(-1) - 1) / (gamma_3 - 1))**0.5
        p4 = p1
        pt4_p4 = (1 + 0.5 * (gamma_3 - 1) * M4**2)**(gamma_3 / (gamma_3 - 1))

        output_array = [Tt4,Tt4 * T4_Tt4,p4,p4 * pt4_p4,M4]
    end function M5_total_to_static_relations

    function M5_entropy_increase(Tt4,T4,Tte,pt4,pte) result(delta_se4)
        ! Entropy increase from state e to 4
        implicit none

        real(dp), intent(in) :: Tt4, T4, Tte, pt4, pte
        real(dp) :: cp, delta_se4

        cp =  M3_pressure_coefficient(T4)
        delta_se4 = cp * log(Tt4 / Tte) - R_3 * log(pt4 / pte)
    end function M5_entropy_increase

    function M5_all_function(p1,pt3,pte,Tte,eta_n,delta_s1e) &
             result(M5_output_array)
        ! Module 5 function that combines all previous functions
        ! and outputs required data for state 4
        implicit none

        real(dp), intent(in) :: p1, pt3, pte, Tte, eta_n, delta_s1e
        real(dp) :: eta_n_ext, T2S_relations(5)
        real(dp) :: M4, pt4, p4, Tt4, T4, cp4, a4, V4, delta_se4, delta_s14
        real(dp) :: M5_output_array(10)

        eta_n_ext = M5_nozzle_external_efficiency(p1,pt3,eta_n)
        
        T2S_relations = M5_total_to_static_relations(eta_n_ext,p1,pte,Tte)
        Tt4 = T2S_relations(1)
        T4 = T2S_relations(2)
        p4 = T2S_relations(3)
        pt4 = T2S_relations(4)
        M4 = T2S_relations(5)

        cp4 = M3_pressure_coefficient(T4)

        a4 = M1_speed_of_sound(T4,gamma_3,R_3)
        V4 = M1_flow_speed(M4,a4)

        delta_se4 = M5_entropy_increase(Tt4,T4,Tte,pt4,pte)
        delta_s14 = delta_s1e + delta_se4

        M5_output_array = [M4,pt4,p4,Tt4,T4,cp4,a4,V4,delta_se4,delta_s14]
    end function M5_all_function
end module module_5