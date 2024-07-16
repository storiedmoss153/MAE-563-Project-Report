module auxiliary_module
    use module_6
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    contains

    subroutine read_initial_parameters(filename,parameters)
        ! Subrouting to read initial parameters found in the file
        ! './data_files/input_variables_case#.csv
        implicit none
    
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: parameters(:)
        real(dp) :: z, M1, eta_d, M2, qf, Tt3_max, eta_n, Area_e
        character(len=20) :: variable_name
        integer :: ios, unit
    
        open(newunit=unit, file=filename,status='old',action='read',iostat=ios)
    
        read(unit,*) variable_name, z
        read(unit,*) variable_name, M1
        read(unit,*) variable_name, eta_d
        read(unit,*) variable_name, M2
        read(unit,*) variable_name, qf
        read(unit,*) variable_name, Tt3_max
        read(unit,*) variable_name, eta_n
        read(unit,*) variable_name, Area_e

        close(unit)
        rewind(unit)
    
        parameters = [z, M1, eta_d, M2, qf, Tt3_max, eta_n, Area_e]
    end subroutine read_initial_parameters

    subroutine write_state_data_to_file(filename,state_array)
        ! subrouting to write an array to a file
        implicit none

        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: state_array(:)
        integer, allocatable :: n_array(:)
        integer :: i, n, unit

        n_array = shape(state_array)
        n = n_array(1)

        open (newunit=unit, file=filename, status='replace', action='write')

        do i = 1, n
            write(unit, '(I0, ",", I0 ",", F20.10)') i, i, state_array(i)
            ! write(unit, '(F20.10)') state_array(i)
        end do

        close(unit)
    end subroutine write_state_data_to_file

    function all_all_module(z,M1,eta_d,M2,qf,Tt3_max,eta_n,Area_e) &
             result(all_output_array)
        implicit none

        real(dp), intent(in) :: z, M1, eta_d, M2, qf, Tt3_max, eta_n, Area_e
        real(dp) :: dmyary_1(7), dmyary_2(8), dmyary_3(9)
        real(dp) :: dmyary_4(11), dmyary_5(10), dmyary_6(13)
        real(dp) :: all_output_array(3)

        dmyary_1 = M1_all_function(z,M1)
        dmyary_2 = M2_all_function(dmyary_1(3),dmyary_1(2),dmyary_1(4), &
                                   M1,M2,eta_d)
        dmyary_3 = M3_all_function(dmyary_2(1),Tt3_max,dmyary_2(4), &
                                   dmyary_2(3),M2,dmyary_2(8))
        dmyary_4 = M4_all_function(dmyary_3(1),dmyary_1(2),dmyary_3(6), &
                                   eta_n,dmyary_3(9),Area_e)
        dmyary_5 = M5_all_function(dmyary_1(2),dmyary_3(6),dmyary_4(2), &
                                   dmyary_4(4),eta_n,dmyary_4(11))
        dmyary_6 = M6_all_function(dmyary_3(3),qf,dmyary_4(9),dmyary_4(8), &
                                   dmyary_1(7),dmyary_4(3),dmyary_1(2), &
                                   dmyary_1(2),Area_e)
        all_output_array = [dmyary_6(12),dmyary_6(6),dmyary_6(8)]
    end function all_all_module
end module auxiliary_module