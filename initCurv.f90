! Xuan 02/28/2016
! 10/07/2017
module initCurvature
	use types
	use global_var, only : T_0
	implicit none
contains
	subroutine initCurv(T_0, Kcap_10, Kcap_20)  !! T_0 = initialcurvature  but not defines what is initialcurvature
		implicit none
		real(dp), dimension(3,3), intent(in)  :: T_0
		real(dp), dimension(3,3), intent(out) :: Kcap_10, Kcap_20
	
		real(dp) :: k_50, k_10, k_610, k_40, k_620, k_20, k_60
		real(dp), dimension(3,3) :: zeros
		zeros = 0._dp

		Kcap_10 = matmul(zeros, transpose(T_0))
		Kcap_20 = matmul(zeros, transpose(T_0))

		k_50  = Kcap_10(1,2)
		k_10  =-Kcap_10(1,3)
		k_610 =-Kcap_10(2,3)

		k_40  = Kcap_20(1,2)
		k_620 =-Kcap_20(1,3)
		k_20  =-Kcap_20(2,3)

		k_60  = k_610 + k_620                    !! so K_10,K_20 are matrices, K_60 is one number
	end subroutine initCurv
end module initCurvature
