! Xuan 02/28/2016
module getTransformation
	use types
	use str_constants, only: k_50, k_10, k_610, k_40, k_620, k_20, k_60
	use global_var 
	implicit none
contains
	subroutine getTransform(U_mat, e1, e2, T1, gamma1)
		implicit none
		real(dp), dimension(24,1), intent(in) :: U_mat                                           !! dimension //
		real(dp), intent(in) :: e1, e2
	
		real(dp), dimension(3,3), intent(out) :: T1
		real(dp), intent(out) :: gamma1

		real(dp) :: R0

		T1(1,1) = (1._dp +U_mat(2,1) -U_mat(7,1)*k_50 +U_mat(13,1)*k_10) /(e1+1._dp)
		T1(1,2) = (U_mat(8,1) +U_mat(1,1)*k_50 +U_mat(13,1)*k_610) /(e1+1._dp)
		T1(1,3) = (U_mat(14,1) -U_mat(1,1)*k_10 -U_mat(7,1)*k_610) /(e1+1._dp)

		T1(2,1) = (U_mat(3,1) -U_mat(7,1)*k_40 +U_mat(13,1)*k_620) /(e2+1._dp)
		T1(2,2) = (1._dp +U_mat(9,1) +U_mat(1,1)*k_40 +U_mat(13,1)*k_20) /(e2+1._dp)
		T1(2,3) = (U_mat(15,1) -U_mat(1,1)*k_620 -U_mat(7,1)*k_20) /(e2+1._dp)

		R0 = sqrt((T1(1,2)*T1(2,3) -T1(1,3)*T1(2,2))**2 +(T1(1,3)*T1(2,1) -T1(1,1)*T1(2,3))**2 +&
				  (T1(1,1)*T1(2,2) -T1(1,2)*T1(2,1))**2)

		T1(3,1) = (T1(1,2)*T1(2,3) -T1(1,3)*T1(2,2)) /R0
		T1(3,2) = (T1(1,3)*T1(2,1) -T1(1,1)*T1(2,3)) /R0
		T1(3,3) = (T1(1,1)*T1(2,2) -T1(1,2)*T1(2,1)) /R0

		gamma1 = asin(T1(1,1)*T1(2,1) + T1(1,2)*T1(2,2) + T1(1,3)*T1(2,3))
	end subroutine getTransform
!!----------------------------------------------------------------------------------------------------------------------
	subroutine getTransformM(U_mat, T1m, gammam)
		implicit none
		real(dp), dimension(24,1), intent(in) :: U_mat
		real(dp), dimension(3,3), intent(out) :: T1m
		real(dp), intent(out) :: gammam

		T1m(1,1) = 1._dp + U_mat(2,1) - U_mat(4,1)*k_50 + U_mat(7,1)*k_10                            !! deleted unnecessary brackets
		T1m(1,2) = U_mat(5,1) + U_mat(1,1)*k_50 + U_mat(7,1)*k_610
		T1m(1,3) = U_mat(8,1) - U_mat(1,1)*k_10 - U_mat(4,1)*k_610

		T1m(2,1) = U_mat(3,1) - U_mat(7,1)*k_40 + U_mat(7,1)*k_620
		T1m(2,2) = 1._dp + U_mat(5,1) + U_mat(1,1)*k_40 + U_mat(7,1)*k_20
		T1m(2,3) = U_mat(9,1) - U_mat(1,1)*k_620 - U_mat(4,1)*k_20

		T1m(3,1) = T1m(1,2)*T1m(2,3) - T1m(1,3)*T1m(2,2)
		T1m(3,2) = T1m(1,3)*T1m(2,1) - T1m(1,1)*T1m(2,3)
		T1m(3,3) = T1m(1,1)*T1m(2,2) - T1m(1,2)*T1m(2,1)

		gammam = 0._dp
	end subroutine getTransformM

end module getTransformation


