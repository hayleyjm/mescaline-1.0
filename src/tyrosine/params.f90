!
! these are replacements for a couple of module variables needed
! to ensure mesacaline routines are compatibile with splash
!
module params
 implicit none
 
! maximum number of columns allowed
 integer, parameter :: maxplot = 512

end module params

module labels
 implicit none
 
! maximum length of column labels
 integer, parameter :: lenlabel = 120

end module labels
