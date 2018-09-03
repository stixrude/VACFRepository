	function gameq(x)
	common /gamcom/ delta

	gameq = 2.*(1. - x)**3/(2. - x) - x**(2./5.)*delta**(3./5.)

c	print*, 'in gameq',x,gameq

	return
	end
