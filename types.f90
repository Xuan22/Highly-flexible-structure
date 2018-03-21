module types
    implicit none
    private
    public dp, hp
    integer, parameter :: dp = kind(0.d0), &              ! double precision
                          hp = selected_real_kind(15)     ! high precision
end module types
