!! Xuan 03/08/2016
module Curvature_mod
	use types
	use str_constants, only:k_50, k_10, k_610, k_40, k_620, k_20
	implicit none
contains
	subroutine Curvature(T_gam, T1dx, T1dy, K_cap1, K_cap2)
		implicit none
		real(dp), dimension(3,3), intent(in) :: T_gam
		real(dp), dimension(3,3), intent(in) :: T1dx, T1dy

		real(dp), dimension(3,3), intent(out) :: K_cap1, K_cap2
		real(dp) :: k1, k2, k61, k62, k5, k4, k6

		integer :: var
		real(dp) :: term1, term2, term3, term4, term5, term6
		term1 = 0._dp
		term2 = 0._dp
		term3 = 0._dp
		term4 = 0._dp
		term5 = 0._dp
		term6 = 0._dp
		do var = 1, 3
			term1 = term1 + T1dx(1,var) *T_gam(3,var)
			term2 = term2 + T1dy(2,var) *T_gam(3,var)
			term3 = term3 + T1dx(2,var) *T_gam(3,var)
			term4 = term4 + T1dy(1,var) *T_gam(3,var)
			term5 = term5 + T1dx(1,var) *T_gam(2,var)
			term6 = term6 + T1dy(2,var) *T_gam(1,var)
		end do

		k1  = -term1 - T_gam(2,1)*k_610 + T_gam(2,2)*k_10  + T_gam(2,3)*k_50
		k2  = -term2 + T_gam(1,1)*k_20  - T_gam(1,2)*k_620 - T_gam(1,3)*k_40
		k61 = -term3 + T_gam(1,1)*k_610 - T_gam(1,2)*k_10  - T_gam(1,3)*k_50
		k62 = -term4 - T_gam(2,1)*k_20  + T_gam(2,2)*k_620 + T_gam(2,3)*k_40
		k5  =  term5 - T_gam(3,1)*k_610 + T_gam(3,2)*k_10  + T_gam(3,3)*k_50
		k4  = -term6 - T_gam(3,1)*k_20  + T_gam(3,2)*k_620 + T_gam(3,3)*k_40
		k6  = k61 + k62
		
		K_cap1 = transpose( reshape( (/0._dp, k5, -k1,  -k5, 0._dp, -k61, k1,  k61, 0._dp/), (/3, 3/) ))                 !! K_cap1 and K_cap2 are not used later??
		K_cap2 = transpose( reshape( (/0._dp, k4, -k62, -k4, 0._dp, -k2,  k62, k2,  0._dp/), (/3, 3/) ))

	end subroutine Curvature

end module Curvature_mod
