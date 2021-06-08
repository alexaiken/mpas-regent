
MODULE module
CONTAINS

subroutine subroutineA(a) bind(c, name="subroutineA_")
    integer, intent(inout) :: a
    a = 7
end subroutine subroutineA

subroutine subroutineB(b) bind(c, name="subroutineB_")
    real(8), dimension(3), intent(inout) :: b
    do i = 1, 3
      b(i) = i + 11.2
    enddo
end subroutine subroutineB

END MODULE module
