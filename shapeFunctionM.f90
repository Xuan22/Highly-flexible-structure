!! Xuan 02/27/2016
module shapeFunctionM
	use types
	use str_constants, only: a, b
	use global_var, only: i, j, k
	implicit none
contains
!! ---------------------------------------------------------------------------------------------------------------------
	subroutine shapeFuncM(xi, yj, Nm, Dm_mat)
		implicit none
		real(dp), intent(in) :: xi, yj
		real(dp), dimension(5,56), intent(out) :: Nm                              !! dimensions
		real(dp), dimension(11,56), intent(out) :: Dm_mat
		real(dp), dimension(2) :: H01, H02, H11, H12
		real(dp), dimension(2) :: H01x, H02x, H11x, H12x
		real(dp), dimension(2) :: H01y, H02y, H11y, H12y
		real(dp), dimension(2) :: L01, L02
		real(dp), dimension(4) :: L_0
		real(dp), dimension(16) :: H_0, Hx, Hy

		H01(1) = (1._dp/a**3) * (2._dp*xi**3 - 3._dp*a*xi**2 + a**3)    !! changed the first (), is this right?
		H02(1) = -(1._dp/a**3) * (2._dp*xi**3 - 3._dp*a*xi**2)
		H11(1) = (1._dp/a**2) * (xi**3 - 2._dp*a*xi**2 + (a**2)*xi)
		H12(1) = (1._dp/a**2) * (xi**3 - a*xi**2)
		
		H01x(1) = (1._dp/a**3)* (6._dp*xi**2 - 6._dp*a*xi)
		H02x(1) = -(1._dp/a**3) * (6._dp*xi**2 - 6._dp*a*xi)
		H11x(1) = (1._dp/a**2) * (3._dp*xi**2 - 4._dp*a*xi + a**2)
		H12x(1) = (1._dp/a**2) * (3._dp*xi**2 - 2._dp*a*xi)

		H01(2) = (1._dp/b**3) * (2._dp*yj**3 - 3._dp*b*yj**2 + b**3)
		H02(2) = -(1._dp/b**3) * (2._dp*yj**3 - 3._dp*b*yj**2)
		H11(2) = (1._dp/b**2) * (yj**3 - 2._dp*b*yj**2 + (b**2)*yj)
		H12(2) = (1._dp/b**2) * (yj**3 - b*yj**2)

!		H01x(2) = 0._dp
!		H02x(2) = 0._dp		
!		H11x(2) = 0._dp
!		H12x(2) = 0._dp

!		H01y(1) = 0._dp
!		H02y(1) = 0._dp
!		H11y(1) = 0._dp
!		H12y(1) = 0._dp

		H01y(2) = (1._dp/b**3) * (6._dp*yj**2 - 6._dp*b*yj)
		H02y(2) = -(1._dp/b**3) * (6._dp*yj**2 - 6._dp*b*yj)
		H11y(2) = (1._dp/b**2) * (3._dp*yj**2 - 4._dp*b*yj + b**2)
		H12y(2) = (1._dp/b**2) * (3._dp*yj**2 - 2._dp*b*yj)
		
		H_0(1) = H01(1) * H01(2)
		H_0(2) = H11(1) * H01(2)
		H_0(3) = H01(1) * H11(2)
		H_0(4) = H11(1) * H11(2)

		H_0(5) = H02(1) * H01(2)
		H_0(6) = H12(1) * H01(2)
		H_0(7) = H02(1) * H11(2)
		H_0(8) = H12(1) * H11(2)

		H_0(9)  = H01(1) * H02(2)
		H_0(10) = H11(1) * H02(2)
		H_0(11) = H01(1) * H12(2)
		H_0(12) = H11(1) * H12(2)

		H_0(13) = H02(1) * H02(2)
		H_0(14) = H12(1) * H02(2)
		H_0(15) = H02(1) * H12(2)
		H_0(16) = H12(1) * H12(2)

		L01(1) = -(1._dp/a) * (xi-a)
		L02(1) = xi/a

		L01(2) = -(1._dp/b) * (yj-b)
		L02(2) = yj/b

		L_0(1) = L01(1) * L01(2)
		L_0(2) = L02(1) * L01(2)
		L_0(3) = L01(1) * L02(2)
		L_0(4) = L02(1) * L02(2)

		Nm = 0._dp
		do i = 1, 4
			do j = 1, 4
				Nm(1, i+(j-1)*14)   = H_0(i+(j-1)*4)                  !! all others are 0 ???
				Nm(2, i+(j-1)*14+4) = H_0(i+(j-1)*4)
				Nm(3, i+(j-1)*14+8) = H_0(i+(j-1)*4)
			end do
				Nm(4, 13+(i-1)*14) = L_0(i)
				Nm(5, 14+(i-1)*14) = L_0(i)
		end do
		   
		Hx(1) = H01x(1) * H01(2)
		Hx(2) = H11x(1) * H01(2)
		Hx(3) = H01x(1) * H11(2)
		Hx(4) = H11x(1) * H11(2)

		Hx(5) = H02x(1) * H01(2)
		Hx(6) = H12x(1) * H01(2)
		Hx(7) = H02x(1) * H11(2)
		Hx(8) = H12x(1) * H11(2)

		Hx(9)  = H02x(1) * H02(2)
		Hx(10) = H12x(1) * H02(2)
		Hx(11) = H02x(1) * H12(2)
		Hx(12) = H12x(1) * H12(2)

		Hx(13) = H01x(1) * H02(2)
		Hx(14) = H11x(1) * H02(2)
		Hx(15) = H01x(1) * H12(2)
		Hx(16) = H11x(1) * H12(2)

		Hy(1) = H01(1) * H01y(2)
		Hy(2) = H11(1) * H01y(2)
		Hy(3) = H01(1) * H11y(2)
		Hy(4) = H11(1) * H11y(2)

		Hy(5) = H02(1) * H01y(2)
		Hy(6) = H12(1) * H01y(2)
		Hy(7) = H02(1) * H11y(2)
		Hy(8) = H12(1) * H11y(2)

		Hy(9)  = H02(1) * H02y(2)
		Hy(10) = H12(1) * H02y(2)
		Hy(11) = H02(1) * H12y(2)
		Hy(12) = H12(1) * H12y(2)

		Hy(13) = H01(1) * H02y(2)
		Hy(14) = H11(1) * H02y(2)
		Hy(15) = H01(1) * H12y(2)
		Hy(16) = H11(1) * H12y(2)
		
		Dm_mat = 0._dp
		do i = 1, 4
			do j = 1, 4
				do k = 1, 3
				    Dm_mat(1+(k-1)*3, i+(j-1)*14+(k-1)*4) = H_0(i+(j-1)*4)
				    Dm_mat(2+(k-1)*3, i+(j-1)*14+(k-1)*4) = Hx(i+(j-1)*4)
				    Dm_mat(3+(k-1)*3, i+(j-1)*14+(k-1)*4) = Hy(i+(j-1)*4)
				end do
			end do
		
			do k = 1, 2
				Dm_mat(10+(k-1)*1,13+(i-1)*14+(k-1)) = L_0(i)
			end do
		end do

	end subroutine shapeFuncM

end module shapeFunctionM



