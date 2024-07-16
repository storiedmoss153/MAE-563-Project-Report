module module_1
    use universal_module
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp), parameter :: Ts = 288                                     ! K
    real(dp), parameter :: ps = 101.3                                   ! kPa
    real(dp), parameter :: R_1 = 286.9                                  ! J/kg-K
    real(dp), parameter :: gamma_1 = 1.40
    real(dp), parameter :: z_star = 8404                                ! m

    contains

    function M1_pressure_coefficient() result(cp)
        ! Pressure Coefficient for states 1 & 2
        implicit none

        real(dp) :: cp

        cp = gamma_1 * R_1 / (gamma_1 - 1)                              ! J/kg-K
    end function M1_pressure_coefficient

    function M1_isentropic_atmos_model(z) result(Temp_n_pres)
        ! Isentropic model of the atmosphere
        implicit none

        real(dp), intent(in) :: z
        real(dp) :: T1_Ts, T1, p1, Temp_n_pres(2)

        if (z <= 7958._dp) then
            T1_Ts = (1 - (gamma_1 - 1) * (z / z_star) / gamma_1)
            T1 = Ts * T1_Ts
            p1 = ps * T1_Ts**(gamma_1 / (gamma_1 - 1))
        else
            T1 = 210
            p1 = 33.6 * exp(-(z - 7958) / 6605)
        end if

        Temp_n_pres = [T1,p1]
    end function M1_isentropic_atmos_model

    function M1_total_to_static_relations(M1,T1,p1) result(out)
        ! Total to static relations for state 1
        implicit none

        real(dp), intent(in) :: M1, T1, p1
        real(dp) :: Tt1_T1, pt1_p1, out(2)
        
        Tt1_T1 = 1 + (gamma_1 - 1) / 2 * M1**2
        pt1_p1 = Tt1_T1**(gamma_1 / (gamma_1 - 1))

        out = [Tt1_T1 * T1,pt1_p1 * p1]
    end function M1_total_to_static_relations

    function M1_speed_of_sound(T,gamma,R) result(a)
        ! Speed of sound calculations for all states
        implicit none

        real(dp), intent(in) :: T, gamma, R
        real(dp) :: a

        a = sqrt(gamma * R * T)
    end function M1_speed_of_sound

    function M1_flow_speed(M,a) result(V)
        ! Air flow speed calculations for all states
        implicit none

        real(dp), intent(in) :: M, a
        real(dp) :: V
        
        V = M * a
    end function M1_flow_speed

    function M1_all_function(z,M1,cp1_in) result(M1_output_matrix)
        ! Module 1 function that combines all previous functions
        ! and outputs required data for state 1
        implicit none

        real(dp), intent(in) :: z, M1
        real(dp), optional :: cp1_in
        real(dp) :: Temp_n_pres(2), T2S_relations(2)
        real(dp) :: M1_output_matrix(7)
        real(dp) :: T1, p1, Tt1, pt1, cp1, a1, V1

        Temp_n_pres = M1_isentropic_atmos_model(z)
        T1 = Temp_n_pres(1)
        p1 = Temp_n_pres(2)
        
        T2S_relations = M1_total_to_static_relations(M1,T1,p1)
        Tt1 = T2S_relations(1)
        pt1 = T2S_relations(2)
        if(present(cp1_in)) then
            cp1 = cp1_in
        else
            cp1 = M1_pressure_coefficient()
        end if

        a1 = M1_speed_of_sound(T1,gamma_1,R_1)
        V1 = M1_flow_speed(M1,a1)

        M1_output_matrix = [T1,p1,Tt1,pt1,cp1,a1,V1]
    end function M1_all_function
end module module_1