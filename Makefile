#
# Parent makefile for Mescaline -- see build/Makefile for details
#

.PHONY: mescaline clean cleanhard
mescaline:
	@cd build; ${MAKE} ${MAKECMDGOALS}

clean:
	@cd build; ${MAKE} clean

cleanhard:
	@cd build; ${MAKE} cleanhard
