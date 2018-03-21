!! Xuan 02/29/2016
module getZm_mod
	use types
	use str_constants, only: k_50, k_10, k_610, k_40, k_620, k_20
	implicit none
contains
	subroutine getZm(T1d, Ud, Zm)
		implicit none
		real(dp), dimension(3,3), intent(in) :: T1d
		real(dp), dimension(24,1), intent(in) :: Ud
	
		real(dp), dimension(18,11), intent(out) :: Zm

		Zm(1,1)  = 1._dp																!! need to compare this whether this is right?
		Zm(2,1) = T1d(2,3)*k_50 + T1d(2,2)*k_10 - T1d(1,3)*k_40 - T1d(1,2)*k_620
		Zm(2,4) = k_610*T1d(2,2) - k_20*T1d(1,2)
		Zm(2,5) = T1d(2,3)
		Zm(2,5) = -T1d(1,3)
		Zm(2,7) = k_610*T1d(2,3) - k_20*T1d(1,3)
		Zm(2,8) = -T1d(2,2)
		Zm(2,9) = T1d(1,2)

		Zm(3,2)  = Ud(10,1)
		Zm(3,4)  = -k_50*Ud(10,1)
		Zm(3,7)  = k_10*Ud(10,1)
		Zm(3,10) = T1d(1,1)

		Zm(4,2)  = Ud(11,1)
		Zm(4,4)  = -k_50*Ud(11,1)
		Zm(4,7)  = k_10*Ud(11,1)
		Zm(4,11) = T1d(1,1)

		Zm(5,3)  = Ud(10,1)
		Zm(5,4)  = -k_40*Ud(10,1)
		Zm(5,7)  = k_620*Ud(10,1)
		Zm(5,10) = T1d(2,1)

		Zm(6,3)  = Ud(11,1)
		Zm(6,4)  = -k_40*Ud(11,1)
		Zm(6,7)  = k_620*Ud(11,1)
		Zm(6,11) = T1d(2,1)

		Zm(7,4)  = 1._dp 

		Zm(8,1) = k_620*T1d(1,1) - k_10*T1d(2,1)
		Zm(8,2) = -T1d(2,3)
		Zm(8,3) = T1d(1,3)
		Zm(8,4) = T1d(1,1)*k_20 + T1d(2,3)*k_50 - T1d(1,3)*k_40-T1d(2,1)*k_610
		Zm(8,7) = k_620*T1d(1,3) - k_10*T1d(2,3)
		Zm(8,8) = T1d(2,1)
		Zm(8,9) = -T1d(1,1)

		Zm(9,1)  = k_50*Ud(10,1)
		Zm(9,5)  = Ud(10,1)
		Zm(9,7)  = k_610*Ud(10,1)
		Zm(9,10) = T1d(1,2)

		Zm(10,1) = k_50*Ud(11,1)
		Zm(10,5) = Ud(11,1)
		Zm(10,7) = k_610*Ud(11,1)
		Zm(10,11) = T1d(1,2)

		Zm(11,1)  = k_40*Ud(10,1)
		Zm(11,6)  = Ud(10,1)
		Zm(11,7)  = k_20*Ud(10,1)
		Zm(11,10) = T1d(2,2)

		Zm(12,1)  = k_40*Ud(11,1)
		Zm(12,6)  = Ud(11,1)
		Zm(12,7)  = k_20*Ud(11,1)
		Zm(12,11) = T1d(2,2)

		Zm(13,7) = 1._dp

		Zm(14,1) = k_40*T1d(1,1) - k_50*T1d(2,1)
		Zm(14,2) = T1d(2,2)
		Zm(14,3) = -T1d(1,2)
		Zm(14,4) = k_40*T1d(1,2) - k_50*T1d(2,2)
		Zm(14,5) = -T1d(2,1)
		Zm(14,6) = T1d(1,1)
		Zm(14,7) = T1d(1,1)*k_20 + T1d(2,2)*k_10 - T1d(1,2)*k_620 - T1d(2,1)*k_610

		Zm(15,1)  = -k_10*Ud(10,1)
		Zm(15,4)  = -k_610*Ud(10,1)
		Zm(15,8)  = Ud(10,1)
		Zm(15,10) = T1d(1,3)

		Zm(16,1)  = -k_10*Ud(11,1)
		Zm(16,4)  = -k_610*Ud(11,1)
		Zm(16,8)  = Ud(11,1)
		Zm(16,11) = T1d(1,3)

		Zm(17,1)  = -k_620*Ud(10,1)
		Zm(17,4)  = -k_20*Ud(10,1)
		Zm(17,9)  = Ud(10,1)
		Zm(17,10) = T1d(2,3)

		Zm(18,1)  = -k_620*Ud(11,1)
		Zm(18,4)  = -k_20*Ud(11,1)
		Zm(18,9)  = Ud(11,1)
		Zm(18,11) = T1d(2,3)

	end subroutine getZm
	
end module getZm_mod
