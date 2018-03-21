!! Xuan 02/27/2016
module shapeFunction
	use types
	use str_constants, only: a, b
	implicit none
contains
	subroutine shapeFunc2(xi, yj, Ns, D_mat)
		implicit none
		real(dp), intent(in) :: xi, yj
		real(dp), dimension(5,56),  intent(out) :: Ns
		real(dp), dimension(24,56), intent(out) :: D_mat
	
		real(dp), dimension(2)  :: H01,   H02,   H11,   H12
		real(dp), dimension(2)  :: H01x,  H02x,  H11x,  H12x
		real(dp), dimension(2)  :: H01xx, H02xx, H11xx, H12xx
		real(dp), dimension(2)  :: L01,   L02,   L01x,  L02x
		real(dp), dimension(2)  :: H01y,  H02y,  H11y,  H12y
		real(dp), dimension(2)  :: H01yy, H02yy, H11yy, H12yy
		real(dp), dimension(2)  :: L01y,  L02y
		real(dp), dimension(4)  :: L_0,   Lx,    Ly
		real(dp), dimension(16) :: H_0,   Hx,    Hy,    Hxx,    Hxy,    Hyy
		
		integer :: i, j, k

		H01(1) = (1._dp/a**3) * (2._dp*xi**3 - 3._dp*a*xi**2 + a**3)
		H02(1) = -(1._dp/a**3) * (2._dp*xi**3 - 3._dp*a*xi**2)
		H11(1) = (1._dp/a**2) * (xi**3 - 2._dp*a*xi**2 + a**2*xi)
		H12(1) = (1._dp/a**2) * (xi**3 - a*xi**2)

		H01x(1) = (1._dp/a**3) * (6._dp*xi**2 - 6._dp*a*xi)
		H02x(1) = -(1._dp/a**3) * (6._dp*xi**2 - 6._dp*a*xi)
		H11x(1) = (1._dp/a**2) * (3._dp*xi**2 - 4._dp*a*xi + a**2)
		H12x(1) = (1._dp/a**2) * (3._dp*xi**2 - 2._dp*a*xi)

		H01x(2) = 0._dp
		H02x(2) = 0._dp
		H11x(2) = 0._dp
		H12x(2) = 0._dp

		H01xx(1) = (1._dp/a**3) * (12._dp*xi - 6._dp*a)
		H02xx(1) = -(1._dp/a**3) * (12._dp*xi - 6._dp*a)
		H11xx(1) = (1._dp/a**2) * (6._dp*xi - 4._dp*a)
		H12xx(1) = (1._dp/a**2) * (6._dp*xi - 2._dp*a)

		L01(1)  = -(1._dp/a) * (xi-a)
		L02(1)  = xi/a
		L01x(1) = -1._dp/a
		L02x(1) = 1._dp/a

	!!
		H01(2) = (1._dp/b**3) * (2._dp*yj**3 - 3._dp*b*yj**2 + b**3)
		H02(2) = -(1._dp/b**3) * (2._dp*yj**3 - 3._dp*b*yj**2)
		H11(2) = (1._dp/b**2) * (yj**3 - 2._dp*b*yj**2 + b**2*yj)
		H12(2) = (1._dp/b**2) * (yj**3 - b*yj**2)

		H01y(1) = 0._dp           						        !! what is H01x(2) or H01y(1)
		H02y(1) = 0._dp
		H11y(1) = 0._dp
		H12y(1) = 0._dp

		H01y(2) = (1._dp/b**3) * (6._dp*yj**2 - 6._dp*b*yj)               !! what is H01x(2) or H01y(1)
		H02y(2) = -(1._dp/b**3) * (6._dp*yj**2 - 6._dp*b*yj)
		H11y(2) = (1._dp/b**2) * (3._dp*yj**2 - 4._dp*b*yj + b**2)
		H12y(2) = (1._dp/b**2) * (3._dp*yj**2 - 2._dp*b*yj)

		H01yy(2) = (1._dp/b**3) * (12._dp*yj - 6._dp*b)
		H02yy(2) = -(1._dp/b**3) * (12._dp*yj - 6._dp*b)
		H11yy(2) = (1._dp/b**2) * (6._dp*yj - 4._dp*b)
		H12yy(2) = (1._dp/b**2) * (6._dp*yj - 2._dp*b)

		L01(2)  = -(1._dp/b) * (yj-b)
		L02(2)  = yj/b
		L01y(2) = -1._dp/b                                      !! L01x(1) and L01y(2)? Is L01x(2) set to 0??
		L02y(2) = 1._dp/b

	!!
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

		L_0(1) = L01(1) * L01(2)
		L_0(2) = L02(1) * L01(2)
		L_0(3) = L01(1) * L02(2)
		L_0(4) = L02(1) * L02(2)
	
		Ns = 0._dp

		do i = 1, 4
			do j = 1, 4
				Ns(1, i + (j-1)*14)   = H_0(i + (j-1)*4)
				Ns(2, i + (j-1)*14+4) = H_0(i + (j-1)*4)
				Ns(3, i + (j-1)*14+8) = H_0(i + (j-1)*4)
			enddo
			Ns(4, 13 + (i-1)*14) = L_0(i)
			Ns(5, 14 + (i-1)*14) = L_0(i)
			  
		enddo
		   
		Hx(1) = H01x(1) * H01(2)
		Hx(2) = H11x(1) * H01(2)
		Hx(3) = H01x(1) * H11(2)
		Hx(4) = H11x(1) * H11(2)

		Hx(5) = H02x(1) * H01(2)
		Hx(6) = H12x(1) * H01(2)
		Hx(7) = H02x(1) * H11(2)
		Hx(8) = H12x(1) * H11(2)

		Hx(9)  = H01x(1) * H02(2)
		Hx(10) = H11x(1) * H02(2)
		Hx(11) = H01x(1) * H12(2)
		Hx(12) = H11x(1) * H12(2)

		Hx(13) = H02x(1) * H02(2)
		Hx(14) = H12x(1) * H02(2)
		Hx(15) = H02x(1) * H12(2)
		Hx(16) = H12x(1) * H12(2)

		Lx(1) = L01x(1) * L01(2)
		Lx(2) = L02x(1) * L01(2)
		Lx(3) = L01x(1) * L02(2)
		Lx(4) = L02x(1) * L02(2)

		Hy(1) = H01(1) * H01y(2)
		Hy(2) = H11(1) * H01y(2)
		Hy(3) = H01(1) * H11y(2)
		Hy(4) = H11(1) * H11y(2)

		Hy(5) = H02(1) * H01y(2)
		Hy(6) = H12(1) * H01y(2)
		Hy(7) = H02(1) * H11y(2)
		Hy(8) = H12(1) * H11y(2)

		Hy(9)  = H01(1) * H02y(2)
		Hy(10) = H11(1) * H02y(2)
		Hy(11) = H01(1) * H12y(2)
		Hy(12) = H11(1) * H12y(2)

		Hy(13) = H02(1) * H02y(2)
		Hy(14) = H12(1) * H02y(2)
		Hy(15) = H02(1) * H12y(2)
		Hy(16) = H12(1) * H12y(2)

		Ly(1) = L01(1) * L01y(2)
		Ly(2) = L02(1) * L01y(2)
		Ly(3) = L01(1) * L02y(2)
		Ly(4) = L02(1) * L02y(2)

		Hxx(1) = H01xx(1) * H01(2)
		Hxx(2) = H11xx(1) * H01(2)
		Hxx(3) = H01xx(1) * H11(2)
		Hxx(4) = H11xx(1) * H11(2)

		Hxx(5) = H02xx(1) * H01(2)
		Hxx(6) = H12xx(1) * H01(2)
		Hxx(7) = H02xx(1) * H11(2)
		Hxx(8) = H12xx(1) * H11(2)

		Hxx(9)  = H01xx(1) * H02(2)
		Hxx(10) = H11xx(1) * H02(2)
		Hxx(11) = H01xx(1) * H12(2)
		Hxx(12) = H11xx(1) * H12(2)

		Hxx(13) = H02xx(1) * H02(2)
		Hxx(14) = H12xx(1) * H02(2)
		Hxx(15) = H02xx(1) * H12(2)
		Hxx(16) = H12xx(1) * H12(2)


		Hxy(1) = H01x(1) * H01y(2)
		Hxy(2) = H11x(1) * H01y(2)
		Hxy(3) = H01x(1) * H11y(2)
		Hxy(4) = H11x(1) * H11y(2)

		Hxy(5) = H02x(1) * H01y(2)
		Hxy(6) = H12x(1) * H01y(2)
		Hxy(7) = H02x(1) * H11y(2)
		Hxy(8) = H12x(1) * H11y(2)

		Hxy(9)  = H01x(1) * H02y(2)
		Hxy(10) = H11x(1) * H02y(2)
		Hxy(11) = H01x(1) * H12y(2)
		Hxy(12) = H11x(1) * H12y(2)

		Hxy(13) = H02x(1) * H02y(2)
		Hxy(14) = H12x(1) * H02y(2)
		Hxy(15) = H02x(1) * H12y(2)
		Hxy(16) = H12x(1) * H12y(2)

		Hyy(1) = H01(1) * H01yy(2)
		Hyy(2) = H11(1) * H01yy(2)
		Hyy(3) = H01(1) * H11yy(2)
		Hyy(4) = H11(1) * H11yy(2)

		Hyy(5) = H02(1) * H01yy(2)
		Hyy(6) = H12(1) * H01yy(2)
		Hyy(7) = H02(1) * H11yy(2)
		Hyy(8) = H12(1) * H11yy(2)

		Hyy(9)  = H01(1) * H02yy(2)
		Hyy(10) = H11(1) * H02yy(2)
		Hyy(11) = H01(1) * H12yy(2)
		Hyy(12) = H11(1) * H12yy(2)

		Hyy(13) = H02(1) * H02yy(2)
		Hyy(14) = H12(1) * H02yy(2)
		Hyy(15) = H02(1) * H12yy(2)
		Hyy(16) = H12(1) * H12yy(2)

		D_mat = 0._dp
		do i = 1, 4
			do j = 1, 4
				do k = 1, 3
				    D_mat(1+(k-1)*6, i+(j-1)*14+(k-1)*4) = H_0(i+(j-1)*4)
				    D_mat(2+(k-1)*6, i+(j-1)*14+(k-1)*4) = Hx(i+(j-1)*4)
				    D_mat(3+(k-1)*6, i+(j-1)*14+(k-1)*4) = Hy(i+(j-1)*4)
				    D_mat(4+(k-1)*6, i+(j-1)*14+(k-1)*4) = Hxx(i+(j-1)*4)
				    D_mat(5+(k-1)*6, i+(j-1)*14+(k-1)*4) = Hxy(i+(j-1)*4)
				    D_mat(6+(k-1)*6, i+(j-1)*14+(k-1)*4) = Hyy(i+(j-1)*4)
				end do
			end do
		
			do k = 1, 2
				    D_mat(19+(k-1)*3, 13+(i-1)*14+(k-1)) = L_0(i)
				    D_mat(20+(k-1)*3, 13+(i-1)*14+(k-1)) = Lx(i)
				    D_mat(21+(k-1)*3, 13+(i-1)*14+(k-1)) = Ly(i)		    !! the maximum is (24, 56)
			end do
		end do

	end subroutine shapeFunc2

end module shapeFunction





