module module_2
    use universal_module
    use module_1
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    contains

    function M2_total_to_static_relations(Tt1,p1,M1,M2,eta_d) result(out)
        ! Total to static relaitons for state 2
        implicit none

        real(dp), intent(in) :: Tt1, p1, M1, M2, eta_d
        real(dp) :: Tt2, T2, pt2, p2, Tt2_T2, out(4)

        Tt2 = Tt1
        Tt2_T2 = (1 + (gamma_1 - 1) / 2 * M2**2)
        T2 = Tt2 / Tt2_T2
        pt2 = p1 * (1 + eta_d * (gamma_1 - 1) / 2 * M1**2)**(gamma_1 / &
                                                            (gamma_1 - 1))
        p2 = pt2 / (Tt2_T2**(gamma_1 / (gamma_1-1)))

        out = [T2,Tt2,p2,pt2]
    end function M2_total_to_static_relations

    function M2_entropy_increase(Tt1,Tt2,pt1,pt2) result(delta_s)
        ! Entropy increase from state 1 to 2
        implicit none

        real(dp), intent(in) :: Tt1, Tt2, pt1, pt2
        real(dp) :: cp, delta_s

        cp = M1_pressure_coefficient()
        delta_s = cp * log(Tt2 / Tt1) - R_1 * log(pt2 / pt1)

    end function M2_entropy_increase
    
    function M2_all_function(Tt1,p1,pt1,M1,M2,eta_d) result(M2_output_matrix)
        ! Module 2 function that combines all previous functions
        ! and outputs required data for state 2
        implicit none

        real(dp), intent(in) :: Tt1, p1, pt1, M1, M2, eta_d
        real(dp) :: T2S_relations(4)
        real(dp) :: M2_output_matrix(8)
        real(dp) :: T2, Tt2, p2, pt2, delta_s12, cp2, a2, V2

        T2S_relations = M2_total_to_static_relations(Tt1,p1,M1,M2,eta_d)
        T2 = T2S_relations(1)
        Tt2 = T2S_relations(2)
        p2 = T2S_relations(3)
        pt2 = T2S_relations(4)

        delta_s12 = M2_entropy_increase(Tt2,Tt2,pt1,pt2)

        cp2 = M1_pressure_coefficient()
        a2 = M1_speed_of_sound(T2,gamma_1,R_1)
        V2 = M1_flow_speed(M2,a2)

        M2_output_matrix = [Tt2,T2,pt2,p2,cp2,a2,V2,delta_s12]
    end function M2_all_function
end module module_2