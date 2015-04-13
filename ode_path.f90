MODULE ode_path
	USE nrtype
	INTEGER(I4B) :: nok,nbad,kount
	LOGICAL(LGT), SAVE :: save_steps=.true.
	REAL(SP) :: dxsav
	REAL(SP), DIMENSION(:), POINTER :: xp
	REAL(SP), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path
