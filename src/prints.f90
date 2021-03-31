module prints
  !
  ! some routines for printing things, i.e. info, warnings, logo, etc
  !
  implicit none
contains

  !
  ! A subroutine to print the main Mescaline logo at the beginning of the run
  !
  subroutine print_mescaline()
    print "(a)",'                                      |_)                 '
    print "(a)",'       __ `__ \   _ \  __|  __|  _` | | | __ \   _ \      '
    print "(a)",'       |   |   |  __/\__ \ (    (   | | | |   |  __/      '
    print "(a)",'      _|  _|  _|\___|____/\___|\__,_|_|_|_|  _|\___|      '
    print "(a)",''
    print "(a)",'   Extracting interesting things from Cactus (HDF5 files) '
    print "(a)", '                See you on the other side.               '
    print "(a)", ''
  end subroutine print_mescaline


  !
  ! A subroutine to print an error message with a particular severity.
  !
  !      severity = 0 --> not that bad, just prints error and moves on
  !      severity = 1 --> pretty bad, you should probably exit, but the code won't crash
  !      severity = 2 --> real bad, the code will (or already has) crash
  !
  subroutine print_error(message,severity,loc)
    integer, intent(in) :: severity
    character(len=*), intent(in) :: message,loc

    if (severity==0) then
       print "(a)", ""
       print "(a)"," --> BAD TRIP ("//trim(loc)//"): "//trim(message)
       print "(a)", ""
       print "(a)", "     Just some gentle nausea. Carrying on. "
       !read*
    elseif (severity==1) then
       print "(a)", ""
       print "(a)"," --> BAD TRIP ("//trim(loc)//"): "//trim(message)
       print "(a)", ""
       print "(a)", "     Maybe re-think your choices. " !Press <enter> to continue. "
       !read*
    elseif (severity==2) then
       print "(a)", ""
       print "(a)"," --> BAD TRIP ("//trim(loc)//"): "//trim(message)
       print "(a)", ""
       print "(a)", "     You've had enough. "
       stop
    endif

  end subroutine print_error


  !
  ! A subroutine to print some general info - to make sure it prints in the same format
  !
  subroutine print_info(message,loc)
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: loc ! where the message is coming from, and advance?

    if (present(loc)) then
       write(*,"(a)") trim(loc)//": "//trim(message)
    else
       write(*,"(a)") trim(message)
    endif

  end subroutine print_info


end module prints
