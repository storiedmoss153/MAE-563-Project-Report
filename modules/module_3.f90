module module_3
    use universal_module
    use module_2
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp), parameter :: R_3 = 286.9                                  ! J/kg-K
    real(dp), parameter :: gamma_3 = 1.30
    real(dp), parameter :: cp_a = 986
    real(dp), parameter :: cp_b = 0.179

    contains

    function M3_pressure_coefficient(T) result(cp)
        ! Pressure Coefficient for states 3, e, & 4
        implicit none

        real(dp), intent(in) :: T
        real(dp) :: cp

        cp = cp_a + (cp_b * T)                                          ! J/kg-K
    end function M3_pressure_coefficient

    function M3_thermal_choking(Tt2,M2) result(Tt3_choked)
        ! Calculations for thermally choked T3 value
        implicit none

        real(dp), intent(in) :: Tt2, M2
        real(dp) :: Tt3_choked

        Tt3_choked = Tt2 * ((1 + gamma_3 * M2**2)**2 / (2 * (gamma_3 + 1) * &
         M2**2) * (1 + 0.5 * (gamma_3 - 1) * M2**2)**(-1))

    end function M3_thermal_choking

    function M3_combustor_conditions(Tt2,Tt3_max,M2) result(out_matrix)
        ! Total temperature and mach for state 3
        implicit none

        real(dp), intent(in) :: Tt2, Tt3_max, M2
        real(dp), allocatable :: Mach(:)
        real(dp) :: Tt3, Tt3_choked, out_matrix(2)
        real(dp) :: a_Mod3, b_Mod3, c_Mod3, c_mach, M3
        integer :: i

        Tt3_choked = M3_thermal_choking(Tt2,M2)

        if (Tt3_choked < Tt3_max) then
            M3 = 1
            Tt3 = Tt3_choked
        else
            Tt3 = Tt3_max
            c_mach = M2**2 * Tt3 * (1 + 0.5 * (gamma_3 - 1) * M2**2) &
                     / (Tt2 * (1 + gamma_3 * M2**2)**2)
            a_Mod3 = c_mach * gamma_3**2 - (gamma_3 - 1) / 2
            b_Mod3 = 2 * c_mach * gamma_3 - 1
            c_Mod3 = c_mach
            Mach = sqrt(quadratic_equation_rr(a_Mod3,b_Mod3,c_Mod3))

            do i = 1,size(Mach)
                if (M2 < 1 .and. Mach(i) <= 1) then
                    M3 = mach(i)
                else if (M2 > 1 .and. Mach(i) >= 1) then
                    M3 = Mach(i)
                else
                    M3 = Mach(i)
                end if
            end do
        end if

        out_matrix = [M3,Tt3]
    end function M3_combustor_conditions

    function M3_heat_added(Tt2,Tt3) result(q23)
        ! Heat per unit mass added in the combustor
        implicit none

        real(dp), intent(in) :: Tt2, Tt3
        real(dp) :: q23

        q23 = cp_a * (Tt3 - Tt2) + 0.5 * cp_b * (Tt3**2 - Tt2**2)
    end function M3_heat_added

    function M3_total_to_static_relations(Tt3,p2,M3) result(Temp_n_pres)
        ! Total to static relations for state 3
        implicit none

        real(dp), intent(in) :: Tt3, p2, M3
        real(dp) :: p3, pt3, T3, Tt3_T3
        real(dp) :: Temp_n_pres(3)

        p3 = p2
        Tt3_T3 = 1 + 0.5 * (gamma_3 - 1) * M3**2
        T3 = Tt3 / Tt3_T3
        pt3 = p3 * Tt3_T3**(gamma_3 / (gamma_3 - 1))
        Temp_n_pres = [T3,p3,pt3]
    end function M3_total_to_static_relations

    function M3_entropy_increase(Tt2,Tt3,pt2,pt3,cp3) result(delta_s23)
        ! Entropy increase from state 2 to 3, and state 3 to e
        implicit none

        real(dp), intent(in) :: Tt2, Tt3, pt2, pt3, cp3
        real(dp) :: delta_s23

        delta_s23 = cp3 * log(Tt3 / Tt2) - R_3 * log(pt3 / pt2)
    end function M3_entropy_increase

    function M3_all_function(Tt2,Tt3_max,p2,pt2,M2,delta_s12,Tt3_in,cp3_in) &
         result(M3_output_matrix)
        ! Module 3 function that combines all previous functions
        ! and outputs required data for state 3
        implicit none

        real(dp), intent(in) :: Tt2, Tt3_max, p2, pt2, M2, delta_s12
        real(dp), optional :: Tt3_in, cp3_in
        real(dp) :: M3, q23, T3, Tt3, p3, cp3, pt3, a3, V3, delta_s23
        real(dp) :: delta_s13
        real(dp) :: comb_conditions(2), T2S_relaitons(3)
        real(dp) :: M3_output_matrix(9)

        comb_conditions = M3_combustor_conditions(Tt2,Tt3_max,M2)
        M3 = comb_conditions(1)

        if(present(Tt3_in)) then
            Tt3 = Tt3_in
        else
            Tt3 = comb_conditions(2)
        end if

        q23 = M3_heat_added(Tt2,Tt3) / 1e6

        T2S_relaitons = M3_total_to_static_relations(Tt3,p2,M3)
        T3 = T2S_relaitons(1)
        p3 = T2S_relaitons(2)
        pt3 = T2S_relaitons(3)

        if(present(cp3_in)) then
            cp3 = cp3_in
        else
            cp3 = M3_pressure_coefficient(T3)
        end if
        a3 = M1_speed_of_sound(T3,gamma_3,R_3)
        V3 = M1_flow_speed(M3,a3)

        delta_s23 = M3_entropy_increase(Tt2,Tt3,pt2,pt3,cp3)
        delta_s13 = delta_s12 + delta_s23

        M3_output_matrix = [Tt3,M3,q23,T3,p3,pt3,cp3,delta_s23,delta_s13]
    end function M3_all_function
end module module_3