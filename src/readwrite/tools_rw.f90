module tools_rw
    use options, only: dtfac,tinitial,c_double
    implicit none

contains

    !
    ! a lil function to calculate the iteration from the current time
    !
    integer function get_it(dx,time)
        real(c_double) :: dx, time, dt
        dt = dtfac * dx
        get_it = nint((time-tinitial)/dt) ! this can be just *below* the it we want when dt is small, so need nint()
    end function get_it


end module tools_rw
