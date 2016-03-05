module cgs_constants
    real (kind=8), parameter :: cgs_boltz = 1.38062e-16, cgs_k = 1.38062e-16,   &
        cgs_mhydr = 1.6733e-24,    &
        cgs_stef = 5.67051e-5, cgs_a = 7.5646e-15, &
        cgs_graw = 6.67259e-8, cgs_msun = 1.99d33, cgs_lsun = 3.9d33,     &
        cgs_c = 2.99792458d10, cgs_h = 6.6260755d-27,  pi = 4d0*atan(1d0),       &
        cgs_kapes = 0.34d0, cgs_hbar= 0.5*cgs_h / pi
end module
