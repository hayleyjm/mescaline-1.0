module gauge
    !
    ! A module to set the gauge used in the simulation
    !
    use options, only:alpgauge,clen,c_double,dp,harmonicF,harmonicN
    use prints, only:print_error
    implicit none

contains

    ! --------------------------------------------
    !
    ! Set the time derivative of the lapse function
    !       --> this is hard-wired for the evolution
    ! Currently have:
    !       - synchronous
    !       - general using harmonicN and harmonicF set in options.f90
    !
    ! --------------------------------------------
    subroutine get_dtalp(alp,tracek,dtalp)
      real(c_double), intent(in) :: alp,tracek
      real(c_double), intent(out) :: dtalp
      character(len=clen) :: loc
      loc = "  get_dtalp"

      select case(alpgauge)
      case("synchronous")
        !
        ! No lapse evolution
        !
        dtalp = 0._dp
      case("harmonic")
        !
        ! Set dt(alp)= -F alp^N K
        !
        dtalp = -harmonicF * alp**harmonicN * tracek
      case default
        call print_error("Please set one of the allowed gauge types.",2,loc)
      end select

    end subroutine get_dtalp


end module gauge
