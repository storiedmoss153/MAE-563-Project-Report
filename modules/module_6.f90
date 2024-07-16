module module_6
    use universal_module
    use module_5
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(dp), parameter :: g_acc = 9.8                                   ! m/s^2

    contains

    function M6_mass_and_fuel_calcualtions(q23,qf,me_dot) result(mass_array)
        ! Calculations for the inlet mass flux (m1_dot)
        ! and fuel-air mass ratio (f)
        implicit none

        real(dp), intent(in) :: q23, qf, me_dot
        real(dp) :: mi_dot, mf_dot, f, mass_array(3)

        mi_dot = me_dot / (1 + q23 / qf)
        mf_dot = me_dot - mi_dot
        f = mf_dot / mi_dot
        mass_array = [mi_dot,mf_dot,f]
    end function M6_mass_and_fuel_calcualtions

    function M6_thrust_calculations(mi_dot,mf_dot,f,Ve,V1,pe,p1,Area_e) &
             result(thrust_array)
        ! Calculations for jet thrust, pressure thrust, total thrust,
        ! thrust specific fuel consumption (TSFC), and specific impulse (Isp)
        implicit none

        real(dp), intent(in) :: mi_dot, mf_dot, f, Ve, V1, pe, p1, Area_e
        real(dp) :: jet_thrust, pressure_thrust, total_thrust, TSFC, Isp
        real(dp) :: thrust_array(5)

        jet_thrust = mi_dot * ((1 + f) * Ve - V1)
        pressure_thrust = 1e3 * (pe - p1) * Area_e
        total_thrust = jet_thrust + pressure_thrust
        TSFC = (mf_dot * 3600) / total_thrust
        Isp = total_thrust / (mf_dot * g_acc)

        thrust_array = [jet_thrust,pressure_thrust,total_thrust,TSFC,Isp]
    end function M6_thrust_calculations

    function M6_efficiency_calculaitons(Ve,V1,q23,pe,total_thrust,p_inf, &
                                        Area_e,me_dot,mi_dot) &
             result(efficiency_array)
        ! Calculations for the fully expanded equivalent exit velocity (Veq),
        ! thermal efficiency (eta_n), propulsive efficiency (eta_p),
        ! overall efficiency (eta_o), and propulsive power (P)
        implicit none

        real(dp), intent(in) :: Ve, V1, q23, pe, total_thrust, p_inf
        real(dp), intent(in) :: Area_e, me_dot, mi_dot
        real(dp) :: Veq, eta_th, eta_p, eta_o, P
        real(dp) :: efficiency_array(5)

        Veq = Ve + 1000 * (pe - p_inf) * Area_e / me_dot
        eta_th = ((me_dot * 0.5 * Veq**2) - (mi_dot * 0.5 * V1**2)) / &
                 (mi_dot * q23 * 1e6)
        eta_p = 2 / (1 + Veq / V1)
        eta_o = eta_th * eta_p
        P = total_thrust * V1

        efficiency_array = [Veq, eta_th, eta_p, eta_o, P]
    end function M6_efficiency_calculaitons
    
    function M6_all_function(q23,qf,me_dot,Ve,V1,pe,p1,p_inf,Area_e) &
             result(M6_output_array)
        ! Module 6 function that combines all previous functions
        ! and outputs required data for all performance perameters
        implicit none

        real(dp), intent(in) :: q23, qf, me_dot, Ve, V1, pe, p1, p_inf, Area_e
        real(dp) :: mi_dot, mf_dot, f
        real(dp) :: jet_thrust, pressure_thrust, total_thrust, TSFC, Isp
        real(dp) :: Veq, eta_th, eta_p, eta_o, prop_power
        real(dp) :: mass_array(3), thrust_array(5), efficiency_array(5)
        real(dp) :: M6_output_array(13)

        mass_array = M6_mass_and_fuel_calcualtions(q23,qf,me_dot)
        mi_dot = mass_array(1)
        mf_dot = mass_array(2)
        f = mass_array(3)

        thrust_array = M6_thrust_calculations(mi_dot,mf_dot,f,Ve,V1, &
                                              pe,p1,Area_e)
        jet_thrust = thrust_array(1)
        pressure_thrust = thrust_array(2)
        total_thrust = thrust_array(3)
        TSFC = thrust_array(4)
        Isp = thrust_array(5)

        efficiency_array = M6_efficiency_calculaitons(Ve,V1,q23,pe, &
                           total_thrust,p_inf,Area_e,me_dot,mi_dot)
        Veq = efficiency_array(1)
        eta_th = efficiency_array(2)
        eta_p = efficiency_array(3)
        eta_o = efficiency_array(4)
        prop_power = efficiency_array(5) / 1e6

        M6_output_array = [mf_dot,mi_dot,f,jet_thrust, &
                           pressure_thrust,total_thrust, &
                           Veq,TSFC,Isp,eta_th,eta_p,eta_o,prop_power]
    end function M6_all_function
end module module_6