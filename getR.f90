!! Xuan 02/29/2016
module getR_mod
	use types
	use str_constants, only:k_50, k_10, k_610, k_40, k_620, k_20
	implicit none
contains

	subroutine getR(r, gamma_1d, gamma_2d, gammad, e1d, e2d, T1d, Td, R_mat)
		implicit none

		real(dp), dimension(6), intent(in) :: r
		real(dp), intent(in) :: gamma_1d, gamma_2d, gammad, e1d, e2d
		real(dp), dimension(3,3), intent(in) :: T1d, Td
	
		real(dp), dimension(24,1), intent(out) :: R_mat
	
		real(dp) :: r1, r2, r3, r4, r5, r6
		real(dp) :: p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
	
		r1 = r(1)
		r2 = r(2)
		r3 = r(3)
		r4 = r(4)
		r5 = r(5)
		r6 = r(6)

		p1 = cos(gamma_1d) / (cos(gammad)*(1._dp+e2d))
		p2 = -sin(gamma_2d) / (cos(gammad)*(1._dp+e1d))

		p3 = sin(gamma_1d) / (cos(gammad)*(1._dp+e2d))
		p4 = -cos(gamma_2d) / (cos(gammad)*(1._dp+e1d))

		p5  = (T1d(2,1)-sin(gammad)*T1d(1,1)) / (2.*cos(gammad)*(1._dp+e1d))
		p6  = (T1d(2,2)-sin(gammad)*T1d(1,2)) / (2.*cos(gammad)*(1._dp+e1d))
		p7  = (T1d(2,3)-sin(gammad)*T1d(1,3)) / (2.*cos(gammad)*(1._dp+e1d))
		p8  = -(T1d(1,1)-sin(gammad)*T1d(2,1)) / (2.*cos(gammad)*(1._dp+e2d))
		p9  = -(T1d(1,2)-sin(gammad)*T1d(2,2)) / (2.*cos(gammad)*(1._dp+e2d))
		p10 = -(T1d(1,3)-sin(gammad)*T1d(2,3)) / (2.*cos(gammad)*(1._dp+e2d))


		R_mat(1,1) = r1+r4*(p1*Td(3,2)*k_40-p1*Td(3,3)*k_620+p2*Td(3,2)*k_50-p2*Td(3,3)*k_10)+ &
						r5*(p3*Td(3,2)*k_40-p3*Td(3,3)*k_620+p4*Td(3,2)*k_50-p4*Td(3,3)*k_10)+ & 
				        r6*(p6*k_50-p7*k_10-p9*k_40-p10*k_620+1._dp/(2._dp*(1._dp+e2d))*(T1d(1,2)*k_40-k_620*T1d(1,3)- &
				        1._dp/(2._dp*(1._dp+e1d))*(T1d(2,1)*k_50-k_10*T1d(2,3))))
		R_mat(2,1) = (p5-1._dp/(2._dp*(1._dp+e1d))*(T1d(2,1)))*r6+r5*(p4*Td(3,1))+r4*p2*Td(3,1)
		R_mat(3,1) = Td(3,1)*(r4*p1+r5*p3)+(p8+1._dp/(2._dp*(1._dp+e2d))*(T1d(1,1)))*r6
		R_mat(4,1) = 0._dp

		R_mat(5,1) = 0._dp
		R_mat(6,1) = 0._dp
		R_mat(7,1) = r2+r4*(-p1*Td(3,1)*k_40-p1*Td(3,3)*k_20-p2*Td(3,1)*k_50-Td(3,3)*p2*k_610)+ &
						r5*(-p3*Td(3,1)*k_40-p3*Td(3,3)*k_20-p4*Td(3,1)*k_50-Td(3,3)*p4*k_610)+ &
				        r6*(-p5*k_50-p7*k_10-p8*k_40-p10*k_20+1._dp/(2._dp*(1._dp+e2d))*(-T1d(1,1)*k_40-T1d(1,3)*k_20)- &
				        1._dp/(2._dp*(1._dp+e1d))*(-T1d(2,1)*k_50-T1d(2,3)*k_610))
		R_mat(8,1) = (p6-1._dp/(2._dp*(1._dp+e1d))*(T1d(2,2)))*r6+r4*(p2*Td(3,2))+r5*p4*Td(3,2)

		R_mat(9,1)  = Td(3,2)*(r4*p1+r5*p3) + (p9+1./(2.*(1.+e2d))*(T1d(1,2)))*r6
		R_mat(10,1) = 0._dp
		R_mat(11,1) = 0._dp
		R_mat(12,1) = 0._dp

		R_mat(13,1) = r3+r4*(p1*Td(3,1)*k_620+p1*Td(3,2)*k_20+p2*Td(3,1)*k_10+Td(3,2)*p2*k_610)+ &
						 r5*(p3*Td(3,1)*k_620+p3*Td(3,2)*k_20+p4*Td(3,1)*k_10-Td(3,2)*p4*k_610)+ &
					     r6*(p5*k_10+p6*k_610-p8*k_620-p9*k_20+1._dp/(2._dp*(1+e2d))*(T1d(1,1)*k_620+T1d(1,2)*k_20)- &
					     1._dp/(2._dp*(1._dp+e1d))*(T1d(2,1)*k_10+T1d(2,3)*k_610))
		R_mat(14,1) = (p7-1._dp/(2._dp*(1._dp+e1d))*(T1d(2,3)))*r6+r4*(p2*Td(3,3))+r5*p4*Td(3,3)
		R_mat(15,1) = Td(3,3)*(r4*p1+r5*p3)+(p10+1._dp/(2._dp*(1._dp+e2d))*(T1d(1,3)))*r6
		R_mat(16,1) = 0._dp

		R_mat(17,1) = 0._dp
		R_mat(18,1) = 0._dp
		R_mat(19,1) = 0._dp
		R_mat(20,1) = 0._dp

		R_mat(21,1) = 0._dp
		R_mat(22,1) = 0._dp
		R_mat(23,1) = 0._dp
		R_mat(24,1) = 0._dp

	end subroutine getR

end module getR_mod
