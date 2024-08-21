module universal_module
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    public pi

    real(dp), parameter :: pi = 4.0*atan(1.0)

    contains

    subroutine linspace(from, to, array)
        implicit none
        real(dp), intent(in) :: from, to
        real(dp), intent(out) :: array(:)
        real(dp) :: range
        integer :: n, i
        n = size(array)
        range = to - from
    
        if (n == 0) return
    
        if (n == 1) then
            array(1) = from
            return
        end if
        
        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end subroutine linspace

    subroutine write_2_arrays_to_file(filename, array1, array2)
        implicit none

        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: array1(:), array2(:)
        integer :: i, n, unit
    
        ! Ensure both arrays have the same size
        n = size(array1)
        if (n /= size(array2)) then
            print *, 'Error: Arrays must have the same size.'
            stop
        end if
    
        ! Open the file for writing
        open (newunit=unit, file=filename, status='replace', action='write')
    
        ! Write the data to the file
        do i = 1, n
            write(unit, '(*(F20.10, ","))') array1(i), array2(i)
        end do
    
        ! Close the file
        close(unit)
    end subroutine write_2_arrays_to_file

    subroutine write_3_arrays_to_file(filename, array1, array2, array3)
        implicit none

        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: array1(:), array2(:), array3(:)
        integer :: i, n, unit
    
        ! Ensure both arrays have the same size
        n = size(array1)
        if (n /= size(array2)) then
            print *, 'Error: Arrays must have the same size.'
            stop
        end if
    
        ! Open the file for writing
        open (newunit=unit, file=filename, status='replace', action='write')
    
        ! Write the data to the file
        do i = 1, n
            write(unit, '(*(F20.10, ","))') &
                  array1(i), array2(i), array3(i)
        end do
    
        ! Close the file
        close(unit)
    end subroutine write_3_arrays_to_file

    subroutine write_4_arrays_to_file(filename, array1, array2, array3, array4)
        implicit none

        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: array1(:), array2(:), array3(:), array4(:)
        integer :: i, n, unit
    
        ! Ensure both arrays have the same size
        n = size(array1)
        if (n /= size(array2)) then
            print *, 'Error: Arrays must have the same size.'
            stop
        end if
    
        ! Open the file for writing
        open (newunit=unit, file=filename, status='replace', action='write')
    
        ! Write the data to the file
        do i = 1, n
            write(unit, '(*(F20.10, ","))') &
                  array1(i), array2(i), array3(i), array4(i)
        end do
    
        ! Close the file
        close(unit)
    end subroutine write_4_arrays_to_file

    function count_lines(filename,nheadder) result(nlines)
        implicit none

        character(len=*), intent(in) :: filename
        integer :: nlines, nheadder, i, iounit
        real :: dummy1, dummy2, dummy3

        nlines = 0

        open (newunit=iounit, file=filename, status='old', action='read')
        do i = 1, nheadder
            read (iounit,*)
        end do
        do
            read (iounit,*, end=10) dummy1, dummy2, dummy3
            nlines = nlines + 1
        end do

10      close (iounit)
        rewind (iounit)
    end function count_lines

    function split_data(filename,nheadder,ncolumn) result(data_matrix)
        implicit none

        character(len=*), intent(in) :: filename
        real, allocatable :: data_matrix(:,:)
        integer :: nheadder, ncolumn, nlines, iounit, i

        nlines = count_lines(filename,nheadder)

        allocate(data_matrix(nlines,ncolumn))

        open (newunit=iounit, file=filename, status='old', action='read')
        if (nheadder /= 0) then
            do i = 1, nheadder
                read (iounit,*)
            end do
        end if
        do i = 1,nlines
            read (iounit,*) data_matrix(i,:)
        end do

        close (iounit)
        rewind (iounit)
    end function split_data

    function quadratic_equation_rr(a,b,c) result(x)
        implicit none

        real(dp), intent(in) :: a, b, c
        real(dp) :: discriminant
        real(dp), allocatable :: x(:)

        discriminant = b**2 - 4*a*c

        if (discriminant > 0) then
            allocate(x(2))
            x(1) = (-b + sqrt(discriminant)) / (2 * a)
            x(2) = (-b - sqrt(discriminant)) / (2 * a)
        else if (discriminant == 0) then
            allocate(x(1))
            x = -b / (2 * a)
        else
            allocate(x(0))
        end if
    end function quadratic_equation_rr
end module universal_module