#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Based on the original solnp by Yinyu Ye
# http://www.stanford.edu/~yyye/Col.html


#----------------------------------------------------------------------------------
# The Function SOLNP solves nonlinear programs in standard form:
#
#        minimize              J(P)
#        subject to            EC(P)  =0
#                   IB(:,1)<=  IC(P)  <=IB(:,2)
#                   PB(:,1)<=    P    <=PB(:,2).
#where
#
#  J       : Cost objective scalar function
#  EC      : Equality constraint vector function
#  IC      : Inequality constraint vector function
#  P       : Decision parameter vector
#  IB, PB  : lower and upper bounds for IC and P.
#----------------------------------------------------------------------------------

# control list
#           RHO  : penalty parameter
#           MAJIT: maximum number of major iterations
#           MINIT: maximum number of minor iterations
#           DELTA: relative step size in forward difference evaluation
#           TOL  : tolerance on feasibility and optimality
# defaults RHO=1, MAJIT=10, MINIT=10, DELTA=1.0e-5, TOL=1.0e-4

solnp = function(pars, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL, UB = NULL, control = list(), ...)
{
	# start timer
	tic = Sys.time()
	xnames = names(pars)
	# get environment
	.solnpenv <- environment()
	assign("xnames", xnames, envir = .solnpenv)
	# initiate function count
	assign(".solnp_nfn", 0, envir = .solnpenv)
	assign(".solnp_errors", 0, envir = .solnpenv)
	
	# index of function indicators
	# [1] length of pars
	# [2] has function gradient?
	# [3] has hessian?
	# [4] has ineq?
	# [5] ineq length
	# [6] has jacobian (inequality)
	# [7] has eq?
	# [8] eq length
	# [9] has jacobian (equality)
	# [10] has upper / lower bounds
	# [11] has either lower/upper bounds or ineq
	
	
	ind = rep(0, 11)
	np = ind[1]  = length(pars)
	# lower parameter bounds - indicator
	# lpb[1]=1 means lower/upper bounds present
	# lpb[2]=1 means lower/upper bounds OR inequality bounds present
	
	# do parameter and LB/UB checks
	check1 = .checkpars(pars, LB, UB, .solnpenv)
	
	# .LB and .UB assigned
	
	.LB = get(".LB", envir = .solnpenv)
	.UB = get(".UB", envir = .solnpenv)
	
	
	if( !is.null(.LB) || !is.null(.UB) ) ind[10] = 1
	
	# do function checks and return starting value
	funv = .checkfun(pars, fun, .solnpenv, ...)
	#.solnp_fun assigned
	.solnp_fun = get(".solnp_fun", envir = .solnpenv)
		
	# Analytical Gradient Functionality not yet implemented in subnp function
	
	# gradient and hessian checks
	#if(!is.null(grad)){
	#	gradv = .checkgrad(pars, grad, .solnpenv, ...)
	#	ind[2] = 1
	#} else{
	#	.solnp_gradfun = function(pars, ...) .fdgrad(pars, fun = .solnp_fun, ...)
		ind[2] = 0
	#	gradv = .solnp_gradfun(pars, ...)
	#}
	# .solnp_gradfun(pars, ...) assigned

	.solnp_hessfun = NULL
	ind[3] = 0
	#hessv = NULL
	# .solnp_hessfun(pars, ...) assigned

	# do inequality checks and return starting values
	
	if(!is.null(ineqfun)){
		ineqv 	= .checkineq(pars, ineqfun, ineqLB, ineqUB, .solnpenv, ...)
		ind[4] 	= 1
		nineq 	= length(ineqLB)
		ind[5] 	= nineq
		
		# check for infinities/nans
		.ineqLBx = .ineqLB
		.ineqUBx = .ineqUB
		.ineqLBx[!is.finite(.ineqLB)] = -1e10
		.ineqUBx[!is.finite(.ineqUB)] =  1e10
		ineqx0 	= (.ineqLBx + .ineqUBx)/2
		#if(!is.null(ineqgrad)){
		#	ineqjacv = .cheqjacineq(pars, gradineq, .ineqUB, .ineqLB, .solnpenv, ...)
		#	ind[6] = 1
		#} else{
		# .solnp_ineqjac = function(pars, ...) .fdjac(pars, fun = .solnp_ineqfun, ...)
		ind[6] = 0
		#ineqjacv = .solnp_ineqjac(pars, ...)
		#}
	} else{
		.solnp_ineqfun = function(pars, ...) .emptyfun(pars, ...)
		# .solnp_ineqjac = function(pars, ...) .emptyjac(pars, ...)
		ineqv 	= NULL
		ind[4] 	= 0
		nineq 	= 0
		ind[5] 	= 0
		ind[6] 	= 0
		ineqx0 	= NULL
		.ineqLB = NULL
		.ineqUB = NULL
	}
	# .solnp_ineqfun and .solnp_ineqjac assigned
	# .ineqLB and .ineqUB assigned
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .solnpenv)
	.ineqLB = get(".ineqLB", envir = .solnpenv)
	.ineqUB = get(".ineqUB", envir = .solnpenv)


	# equality checks
	if(!is.null(eqfun)){
		eqv 	= .checkeq(pars, eqfun, eqB, .solnpenv, ...)
		ind[7] 	= 1
		.eqB = get(".eqB", envir = .solnpenv)
		neq 	= length(.eqB)
		ind[8] 	= neq
		#if(!is.null(eqgrad)){
		#	eqjacv = .cheqjaceq(pars, gradeq, .solnpenv, ...)
		#	ind[9] = 1
		#} else{
		#	.solnp_eqjac = function(pars, ...) .fdjac(pars, fun = .solnp_eqfun, ...)
		#	eqjacv = .solnp_eqjac(pars, ...)
			ind[9] = 0
		#}
	} else {
		eqv = NULL
		#eqjacv = NULL
		.solnp_eqfun = function(pars, ...) .emptyfun(pars, ...)
		#.solnp_eqjac = function(pars, ...) .emptyjac(pars, ...)
		ind[7] 	= 0
		neq 	= 0
		ind[8] 	= 0
		ind[9] 	= 0
	}
	# .solnp_eqfun(pars, ...) and .solnp_eqjac(pars, ...) assigned
	# .solnp_eqB assigned
	.solnp_eqfun = get(".solnp_eqfun", envir = .solnpenv)

	if( ind[ 10 ] || ind [ 4 ]) ind[ 11 ] = 1
		
	# parameter bounds (pb)
	pb  = rbind( cbind(.ineqLB, .ineqUB), cbind(.LB, .UB) )
	
	# check control list
	ctrl  = .solnpctrl( control )
	rho   = ctrl[[ 1 ]]
	# maxit = outer iterations
	maxit = ctrl[[ 2 ]]
	# minit = inner iterations
	minit = ctrl[[ 3 ]]
	delta = ctrl[[ 4 ]]
	tol   = ctrl[[ 5 ]]
	trace = ctrl[[ 6 ]]
	
	# total constraints (tc) = no.inequality constraints + no.equality constraints
	tc = nineq + neq
	
	# initialize fn value and inequalities and set to NULL those not needed
	j  = jh = funv
	tt = 0 * .ones(3, 1)
	
	if( tc > 0 ) {
		# lagrange multipliers (lambda)
		lambda = 0 * .ones(tc, 1)
		# constraint vector = [1:neq 1:nineq]
		constraint = c(eqv, ineqv)
		if( ind[4] ) {
			tmpv = cbind(constraint[ (neq + 1):tc ] - .ineqLB, .ineqUB - constraint[ (neq + 1):tc ] )
			testmin = apply( tmpv, 1, FUN = function( x ) min(x[ 1 ], x[ 2 ]) )
			if( all(testmin > 0) ) ineqx0 = constraint[ (neq + 1):tc ]
			constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - ineqx0
		}
		tt[ 2 ] = .vnorm(constraint)
		if( max(tt[ 2 ] - 10 * tol, nineq, na.rm = TRUE) <= 0 ) rho = 0
	} else{
		lambda = 0
	}
	# starting augmented parameter vector
	p  = c(ineqx0, pars)
	hessv  = diag(np + nineq)
	mu = np
	.solnp_iter = 0
	ob = c(funv, eqv, ineqv)
	
	while( .solnp_iter < maxit ){
		.solnp_iter = .solnp_iter + 1
		.subnp_ctrl = c(rho, minit, delta, tol, trace)
		
		# make the scale for the cost, the equality constraints, the inequality
		# constraints, and the parameters
		if( ind[7] ) {
			# [1 neq]
			vscale = c( ob[ 1 ], rep(1, neq) * max( abs(ob[ 2:(neq + 1) ]) ) )
		} else {
			vscale = 1
		}
		
		if( !ind[ 11 ] ) {
			vscale = c(vscale, p)
		} else {
			# [ 1 neq np]
			vscale = c(vscale, rep( 1, length.out = length(p) ) )
		}
		
		vscale = apply( matrix(vscale, ncol = 1), 1, FUN = function( x ) min( max( abs(x), tol ), 1/tol ) )
		
		res   = .subnp(pars = p, yy = lambda, ob = ob, hessv = hessv, lambda = mu, vscale = vscale, 
				ctrl = .subnp_ctrl, .env = .solnpenv, ...)
		if(get(".solnp_errors", envir =  .solnpenv) == 1){
			maxit = .solnp_iter
		}
		p  = res$p
		lambda  = res$y
		hessv  = res$hessv
		mu = res$lambda
		temp = p[ (nineq + 1):(nineq + np) ]
		funv = .safefunx(temp, .solnp_fun, .env = .solnpenv, ...)
		ctmp = get(".solnp_nfn", envir =  .solnpenv)
		assign(".solnp_nfn", ctmp + 1, envir = .solnpenv)
		
		tempdf = cbind(temp, funv)
		
		if( trace ){
			.report(.solnp_iter, funv, temp)
		}
		
		eqv = .solnp_eqfun(temp, ...)		
		ineqv = .solnp_ineqfun(temp, ...)
		
		ob = c(funv, eqv, ineqv)
		
		tt[ 1 ] = (j - ob[ 1 ]) / max(abs(ob[ 1 ]), 1)
		j = ob[ 1 ]
		
		if( tc > 0 ){
			constraint = ob[ 2:(tc + 1) ]
			
			if( ind[ 4 ] ){
				tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ],
				              pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )
				              
				if( min(tempv) > 0 ) {
					p[ 1:nineq ] = constraint[ (neq + 1):tc ]
				}
				
				constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - p[ 1:nineq ]
			}
			
			tt[ 3 ] = .vnorm(constraint)
			
			if( tt[ 3 ] < 10 * tol ) { 
				rho = 0
				mu  = min(mu, tol)
			}
			
			if( tt[ 3 ] < 5 * tt[ 2 ]) {
				rho = rho/5
			}
			
			if( tt[ 3 ] > 10 * tt[ 2 ]) {
				rho = 5 * max( rho, sqrt(tol) )
			}
			
			if( max( c( tol + tt[ 1 ], tt[ 2 ] - tt[ 3 ] ) ) <= 0 ) { 
				lambda = 0
				hessv = diag( diag ( hessv ) )
			}

			tt[ 2 ] = tt[ 3 ]
		}
		
		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			maxit = .solnp_iter
		}
		
		jh = c(jh, j)
	}
	
	if( ind[ 4 ] ) {
		ineqx0 = p[ 1:nineq ]
	}
	
	p = p[ (nineq + 1):(nineq + np) ]
	
	if(get(".solnp_errors", envir =  .solnpenv) == 1){
		convergence = 2
		if( trace ) cat( paste( "\nsolnp--> Solution not reliable....Problem Inverting Hessian.\n", sep="" ) )
	} else{
		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			convergence = 0
			# if( trace ) cat( paste( "\nsolnp--> Completed in ", .solnp_iter, " iterations\n", sep="" ) )
		} else{
			convergence = 1
			# if( trace ) cat( paste( "\nsolnp--> Exiting after maximum number of iterations\n",
							# "Tolerance not achieved\n", sep="" ) )
		}
	}
	# end timer
	ctmp = get(".solnp_nfn", envir =  .solnpenv)
	toc = Sys.time() - tic
	names(p) = xnames
	ans = list(pars = p, convergence = convergence, values = as.numeric(jh), lagrange = lambda, 
			hessian = hessv, ineqx0 = ineqx0, nfuneval = ctmp, outer.iter = .solnp_iter, 
			elapsed = toc, vscale = vscale)
	return( ans )
}

#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# Based on the original subnp Yinyu Ye
# http://www.stanford.edu/~yyye/Col.html
.subnp = function(pars, yy, ob, hessv, lambda, vscale, ctrl, .env, ...)
{
	.solnp_fun = get(".solnp_fun", envir = .env)
	.solnp_eqfun = get(".solnp_eqfun", envir = .env)
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .env)
	ineqLB = get(".ineqLB", envir = .env)
	ineqUB = get(".ineqUB", envir = .env)
	LB = get(".LB", envir = .env)
	UB = get(".UB", envir = .env)
	#.solnp_gradfun = get(".solnp_gradfun", envir = .env)
	#.solnp_eqjac = get(".solnp_eqjac", envir = .env)
	#.solnp_ineqjac = get(".solnp_ineqjac", envir = .env)
	ind = get("ind", envir = .env)
	
	# pars [nineq + np]	
	rho   = ctrl[ 1 ]
	maxit = ctrl[ 2 ]
	delta = ctrl[ 3 ]
	tol   = ctrl[ 4 ]
	trace = ctrl[ 5 ]
	# [1] length of pars
	# [2] has function gradient?
	# [3] has hessian?
	# [4] has ineq?
	# [5] ineq length
	# [6] has jacobian (inequality)
	# [7] has eq?
	# [8] eq length
	# [9] has jacobian (equality)
	# [10] has upper / lower bounds
	# [11] has either lower/upper bounds or ineq
	
	
	neq   = ind[ 8 ]
	nineq = ind[ 5 ]
	np    = ind[ 1 ]
	ch    = 1
	alp   = c(0,0,0)
	nc    = neq + nineq
	npic  = np + nineq
	p0    = pars
	
	# pb [ 2 x (nineq + np) ]
	pb    = rbind( cbind(ineqLB, ineqUB), cbind(LB,UB) )
	sob   = numeric()
	ptt   = matrix()
	sc    = numeric()
	
	# scale the cost, the equality constraints, the inequality constraints, 
	# the parameters (inequality parameters AND actual parameters), 
	# and the parameter bounds if there are any
	# Also make sure the parameters are no larger than (1-tol) times their bounds
	# ob [ 1 neq nineq]
	
	ob = ob / vscale[ 1:(nc + 1) ]
	# p0 [np]
	p0 = p0 / vscale[ (neq + 2):(nc + np + 1) ]
	if( ind[ 11 ] ) {
		
		if( !ind[ 10 ] ) {
			mm = nineq
		} else {
			mm = npic
		}
		
		pb = pb / cbind(vscale[ (neq + 2):(neq + mm + 1) ], vscale[ (neq + 2):(neq + mm + 1) ])
	}

	# scale the lagrange multipliers and the Hessian
	
	if( nc > 0) {
		# yy [total constraints = nineq + neq]
		# scale here is [tc] and dot multiplied by yy
		yy = vscale[ 2:(nc + 1) ] * yy / vscale[ 1 ]
	}
	# yy = [zeros 3x1]
	
	# h is [ (np+nineq) x (np+nineq) ]
	#columnvector %*% row vector (size h) then dotproduct h then dotdivide scale[1]
	hessv = hessv * (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1)]) ) / vscale[ 1 ]
	# h[ 8x8 eye]
	j = ob[ 1 ]
	
	if( ind[4] ) {
		
		if( !ind[7] ) {
			# [nineq x (nineq+np) ]
			a = cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) ) 
		} else {
			# [ (neq+nineq) x (nineq+np)]
			a = rbind( cbind( 0 * .ones(neq, nineq), matrix(0, ncol = np, nrow = neq) ), 
					cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) ) )
		}
		
	}
	if( ind[7] && !ind[4] ) {
		a = .zeros(neq, np)
	}
	
	if( !ind[7] && !ind[4] ) {
		a = .zeros(1, np)
	}
	
	# gradient
	g= 0 * .ones(npic, 1)
	p = p0 [ 1:npic ]
	if( nc > 0 ) {
		# [ nc ]
		constraint = ob[ 2:(nc + 1) ]
		# constraint [5 0 11 3x1]
		# gradient routine
		for( i in 1:np ) {
			# scale the parameters (non ineq)
			p0[ nineq + i ] = p0[ nineq + i ] + delta
			tmpv = p0[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
			funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
			eqv 	= .solnp_eqfun(tmpv, ...)
			ineqv 	= .solnp_ineqfun(tmpv, ...)
			ctmp = get(".solnp_nfn", envir =  .env)
			assign(".solnp_nfn", ctmp + 1, envir = .env)

			ob = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
			g[ nineq + i ]   = (ob[ 1 ] - j) / delta
			a[ , nineq + i ] = (ob[ 2:(nc + 1) ] - constraint) / delta
			p0[ nineq + i ]  = p0[ nineq + i ] - delta
		}
		
		if( ind[4] ) {
			constraint[ (neq + 1):(neq + nineq) ] = constraint[ (neq + 1):(neq + nineq) ] - p0[ 1:nineq ]
		}
		
		# solver messages
		if( .solvecond(a) > 1 / .eps ) { 
			if( trace ) .subnpmsg( "m1" )
		}

		# a(matrix) x columnvector - columnvector
		# b [nc,1]
		b  = a %*% p0 - constraint
		ch = -1
		alp[ 1 ] = tol - max( abs(constraint) )
		
		if( alp[ 1 ] <= 0 ) {
			ch = 1
			
			if( !ind[11] ) {
				# a %*% t(a) gives [nc x nc]
				# t(a) %*% above gives [(np+nc) x 1]
				p0 = p0 - t(a) %*% solve(a %*% t(a), constraint)
				alp[ 1 ] = 1
			}

		}
		
		if( alp[ 1 ] <= 0 ) {
			# this expands p0 to [nc+np+1]
			p0[ npic + 1 ] = 1
			a  = cbind(a, -constraint)
			# cx is rowvector
			cx = cbind(.zeros(1,npic), 1)
			dx = .ones(npic + 1, 1)
			go = 1 
			minit = 0
			
			while( go >= tol ) {
				minit = minit + 1
				# gap [(nc + np) x 2]
				gap = cbind(p0[ 1:mm ] - pb[ , 1 ], pb[ , 2 ] - p0[ 1:mm ] )
				# this sorts every row
				gap = t( apply( gap, 1, FUN=function( x ) sort(x) ) )
				dx[ 1:mm ] = gap[ , 1 ]
				# expand dx by 1
				dx[ npic + 1 ] = p0[ npic + 1 ]
				
				if( !ind[10] ) {
					dx[ (mm + 1):npic ] = max( c(dx[ 1:mm ], 100) ) * .ones(npic - mm, 1)
				}
				# t( a %*% diag( as.numeric(dx) ) ) gives [(np+nc + 1 (or more) x nc]
				# dx * t(cx) dot product of columnvectors
				# qr.solve returns [nc x 1]
				
				# TODO: Catch errors here
				y = try( qr.solve( t( a %*% diag( as.numeric(dx) , length(dx), length(dx) ) ), dx * t(cx) ), silent = TRUE)
				if(inherits(y, "try-error")){
					p = p0 * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
					if( nc > 0 ) {
						y = 0 # unscale the lagrange multipliers
					}
					hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
					ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
					assign(".solnp_errors", 1, envir = .env)
					return(ans)
				}
				v = dx * ( dx *(t(cx) - t(a) %*% y) )
				
				if( v[ npic + 1 ] > 0 ) {
					z = p0[ npic + 1 ] / v[ npic + 1 ]
					
					for( i in 1:mm ) {
					
						if( v[ i ] < 0 ) {
							z = min(z, -(pb[ i, 2 ] - p0[ i ]) / v[ i ])
						} else if( v[ i ] > 0 ) { 
							z = min( z, (p0[ i ] - pb[ i , 1 ]) / v[ i ]) 
						}
					}
					
					if( z >= p0[ npic + 1 ] / v[ npic + 1 ] ) {
						p0 = p0 - z * v
					} else {
						p0 = p0 - 0.9 * z * v 
					}
					go = p0[ npic + 1 ]
					
					if( minit >= 10 ) {
						go = 0 
					}
					
				} else {
					go = 0
					minit = 10
				}
				
			}
			
			if( minit >= 10 ) {
				if( trace ) .subnpmsg( "m2" )
			}
			
			a = matrix(a[ , 1:npic ], ncol = npic)
			b = a %*% p0[ 1:npic ]
		}
		
	}
	
	p = p0 [ 1:npic ]
	y = 0
	
	if( ch > 0 ) {
		
		tmpv = p[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
		funv = .safefunx(tmpv, .solnp_fun, .env,...)
		eqv = .solnp_eqfun(tmpv, ...)
		ineqv = .solnp_ineqfun(tmpv, ...)
		ctmp = get(".solnp_nfn", envir =  .env)
		assign(".solnp_nfn", ctmp + 1, envir = .env)
		ob = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
	}
	
	j = ob[ 1 ]
	
	if( ind[4] ) {
		ob[ (neq + 2):(nc + 1) ] = ob[ (neq + 2):(nc + 1) ] - p[ 1:nineq ]

	}
	
	if( nc > 0 ) {
		ob[ 2:(nc + 1) ] = ob[ 2:(nc + 1) ] - a %*% p + b
		j = ob[ 1 ] - t(yy) %*% matrix(ob[ 2:(nc + 1) ], ncol=1) + rho * .vnorm(ob[ 2:(nc + 1) ]) ^ 2
	}
	
	minit = 0
	while( minit < maxit ) {
		minit = minit + 1
		
		if( ch > 0 ) {
		
			for( i in 1:np ) {
				
				p[ nineq + i ] = p[ nineq + i ] + delta
				tmpv = p[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
				funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
				eqv 	= .solnp_eqfun(tmpv, ...)
				ineqv 	= .solnp_ineqfun(tmpv, ...)
				ctmp = get(".solnp_nfn", envir =  .env)
				assign(".solnp_nfn", ctmp + 1, envir = .env)
				obm = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
				
				if( ind[4] ) {
					obm[ (neq + 2):(nc + 1)] = obm[ (neq + 2):(nc + 1) ] - p[ 1:nineq ]
				}
				
				if( nc > 0 ) {
					
					obm[ 2:(nc + 1) ] = obm[ 2:(nc + 1) ] - a %*% p + b
					obm = obm[ 1 ] - t(yy) %*% obm[ 2:(nc + 1) ] + rho * .vnorm(obm[ 2:(nc + 1 ) ]) ^ 2
				}
				
				g[ nineq + i ] = (obm - j) / delta
				p[ nineq + i ] = p[ nineq + i ] - delta
			}
			
			if( ind[4] ) {
				g[ 1:nineq ] = 0
			}
			
		}
		
		if( minit > 1 ) {
			yg = g - yg
			sx = p - sx
			sc[ 1 ] = t(sx) %*% hessv %*% sx
			sc[ 2 ] = t(sx) %*% yg
			
			if( (sc[ 1 ] * sc[ 2 ]) > 0 ) {
				sx = hessv %*% sx
				hessv  = hessv - ( sx %*% t(sx) ) / sc[ 1 ] + ( yg %*% t(yg) ) / sc[ 2 ]
			}
			
		}
		
		dx = 0.01 * .ones(npic, 1)
		if( ind[11] ) {
			
			gap = cbind(p[ 1:mm ] - pb[ , 1 ], pb[ , 2 ] - p[ 1:mm ])
			gap = t( apply( gap, 1, FUN = function( x ) sort(x) ) )
			gap = gap[ , 1 ] + sqrt(.eps) * .ones(mm, 1)
			dx[ 1:mm, 1 ] = .ones(mm, 1) / gap
			if( !ind[10] ){
				dx[ (mm + 1):npic, 1 ] = min (c( dx[ 1:mm, 1 ], 0.01) ) * .ones(npic - mm, 1)
			}
			
		}

		go = -1
		lambda = lambda / 10
		while( go <= 0 ) {
			cz = try(chol( hessv + lambda * diag( as.numeric(dx * dx), length(dx), length(dx) ) ),  silent = TRUE)
			if(inherits(cz, "try-error")){
				p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
				if( nc > 0 ) {
					y = 0 # unscale the lagrange multipliers
				}
				hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
				ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
				assign(".solnp_errors", 1, envir = .env)
				return(ans)
			}
			cz = try(solve(cz), silent = TRUE)
			if(inherits(cz, "try-error")){
				p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
				if( nc > 0 ) {
					y = 0 # unscale the lagrange multipliers
				}
				hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
				ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
				assign(".solnp_errors", 1, envir = .env)
				return(ans)
			}
			yg = t(cz) %*% g
			
			if( nc == 0 ) {
				u = -cz %*% yg
			} else{
				y = try( qr.solve(t(cz) %*% t(a), yg), silent = TRUE )
				if(inherits(y, "try-error")){
					p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
					if( nc > 0 ) {
						# y = vscale[ 1 ] * y / vscale[ 2:(nc + 1) ] # unscale the lagrange multipliers
						y = 0
					}
					hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
					ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
					assign(".solnp_errors", 1, envir = .env)
					return(ans)
				}
				
				u = -cz %*% (yg - ( t(cz) %*% t(a) ) %*% y)
			}
			
			p0 = u[ 1:npic ] + p
			if( !ind[ 11 ] ) {
				go = 1
			} else {
				go = min( c(p0[ 1:mm ] - pb[ , 1 ], pb[ , 2 ] - p0[ 1:mm ]) )
				lambda = 3 * lambda
			}
			
		}
		
		alp[ 1 ] = 0
		ob1 = ob
		ob2 = ob1
		sob[ 1 ] = j
		sob[ 2 ] = j
		ptt = cbind(p, p)
		alp[ 3 ] = 1.0
		ptt = cbind(ptt, p0)
		tmpv = ptt[ (nineq + 1):npic, 3 ] * vscale[ (nc + 2):(nc + np + 1) ]
		funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
		eqv 	= .solnp_eqfun(tmpv, ...)
		ineqv 	= .solnp_ineqfun(tmpv, ...)
		ctmp = get(".solnp_nfn", envir =  .env)
		assign(".solnp_nfn", ctmp + 1, envir = .env)
		
		ob3 = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
		sob[ 3 ] = ob3[ 1 ]
		
		if( ind[4] ) {
			ob3[ (neq + 2):(nc + 1) ] = ob3[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq, 3 ]
		}
		
		if( nc > 0 ) {
			ob3[ 2:(nc + 1) ] = ob3[ 2:(nc + 1) ] - a %*% ptt[ , 3 ] + b
			sob[ 3 ] = ob3[ 1 ] - t(yy) %*% ob3[ 2:(nc + 1) ] + rho * .vnorm(ob3[ 2:(nc + 1) ]) ^ 2
		}
		
		go = 1
		while( go > tol ) {
			alp[ 2 ] = (alp[ 1 ] + alp[ 3 ]) / 2
			ptt[ , 2 ] = (1 - alp[ 2 ]) * p + alp[ 2 ] * p0
			tmpv = ptt[ (nineq + 1):npic, 2 ] * vscale[ (nc + 2):(nc + np + 1) ]
			funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
			eqv 	= .solnp_eqfun(tmpv, ...)
			ineqv 	= .solnp_ineqfun(tmpv, ...)
			ctmp = get(".solnp_nfn", envir =  .env)
			assign(".solnp_nfn", ctmp + 1, envir = .env)
			
			ob2 = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
			
			sob[ 2 ] = ob2[ 1 ]

			if( ind[4] ) {
				ob2[ (neq + 2):(nc + 1) ] = ob2[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq , 2 ]
			}
			
			if( nc > 0 ) {
				ob2[ 2:(nc + 1) ] = ob2[ 2:(nc + 1) ] - a %*% ptt[ , 2 ] + b
				sob[ 2 ] = ob2[ 1 ] - t(yy) %*% ob2[ 2:(nc + 1) ] + rho * .vnorm(ob2[ 2:(nc + 1) ]) ^ 2
			}
			
			obm = max(sob)
			
			if( obm < j ) {
				obn = min(sob)
				go = tol * (obm - obn) / (j - obm)
			}

			condif1 = sob[ 2 ] >= sob[ 1 ]
			condif2 = sob[ 1 ] <= sob[ 3 ] && sob[ 2 ] < sob[ 1 ]
			condif3 = sob[ 2 ] <  sob[ 1 ] && sob[ 1 ] > sob[ 3 ]
			
			if( condif1 ) {
				sob[ 3 ] = sob[ 2 ]
				ob3 = ob2
				alp[ 3 ] = alp[ 2 ]
				ptt[ , 3 ] = ptt[ , 2 ]
			}
			
			if( condif2 ) {
				sob[ 3 ] = sob[ 2 ]
				ob3 = ob2
				alp[ 3 ] = alp[ 2 ]
				ptt[ , 3 ] = ptt[ , 2 ]
			}
			
			if( condif3 ) {
				sob[ 1 ] = sob[ 2 ]
				ob1 = ob2
				alp[ 1 ] = alp[ 2 ]
				ptt[ , 1 ] = ptt[ , 2 ]
			}
			
			if( go >= tol ) {
				go = alp[ 3 ] - alp[ 1 ]
			}
			
		}
		
		sx = p
		yg = g
		ch = 1
		obn = min(sob)
		if( j <= obn ) {
			maxit = minit
		}
		
		reduce = (j - obn) / ( 1 + abs(j) )
		
		if( reduce < tol ) {
			maxit = minit
		}
		
		condif1 = sob[ 1 ] <  sob[ 2 ]
		condif2 = sob[ 3 ] <  sob[ 2 ] && sob[ 1 ] >= sob[ 2 ]
		condif3 = sob[ 1 ] >= sob[ 2 ] && sob[ 3 ] >= sob[ 2 ]
		
		if( condif1 ) {
			j = sob[ 1 ]
			p = ptt[ , 1 ]
			ob = ob1
		}
		
		if( condif2 ) {
			j = sob [ 3 ]
			p = ptt[ , 3 ]
			ob = ob3
		}
		
		if( condif3 ) {
			j = sob[ 2 ]
			p = ptt[ , 2 ]
			ob = ob2
		}
		
	}
	
	p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
	
	if( nc > 0 ) {
		y = vscale[ 1 ] * y / vscale[ 2:(nc + 1) ] # unscale the lagrange multipliers
	}
	
	hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )	
	if( reduce > tol ) {
		if( trace ) .subnpmsg( "m3" )
	}

	ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
	return( ans )
}


.fdgrad = function(i, p, delta, np, vscale, constraint, j, nineq, npic, nc, .solnp_fun, 
		.solnp_eqfun, .solnp_ineqfun, .env, ...)
{
	ans = list()
	px = p
	px[ nineq + i ] = px[ nineq + i ] + delta
	tmpv = px[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
	funv = .safefunx(tmpv, .solnp_fun, .env, ...)
	eqv 	= .solnp_eqfun(tmpv, ...)
	ineqv 	= .solnp_ineqfun(tmpv, ...)
	ctmp = get(".solnp_nfn", envir =  .env)
	assign(".solnp_nfn", ctmp + 1, envir = .env)
	ob = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
	ans$p	= px
	ans$ob 	= ob
	ans$g   = (ob[ 1 ] - j) / delta
	ans$a   = (ob[ 2:(nc + 1) ] - constraint) / delta
	return( ans )
}

#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

.eps = .Machine$double.eps

.subnpmsg = function(m){
	g1 = c("solnp-->")
	m1 = paste("\n",g1, "Redundant constraints were found. Poor\n",
			g1, "intermediate results may result.Suggest that you\n",
			g1, "remove redundant constraints and re-OPTIMIZE\n", sep = "")
	m2 = paste("\n",g1, "The linearized problem has no feasible\n",
			g1, "solution.  The problem may not be feasible.\n", sep = "")
	m3 = paste("\n",g1, "Minor optimization routine did not converge in the \n",
			g1, "specified number of minor iterations.  You may need\n",
			g1, "to increase the number of minor iterations.        \n", sep = "")
	ans = switch(m,
			m1 = m1,
			m2 = m2,
			m3 = m3)
	cat(ans)
}



.checkpars = function(pars, LB, UB, .env)
{
	if(is.null(pars))
		stop("\nsolnp-->error: must supply starting parameters\n", call. = FALSE)
	if(!is.null(LB)){
		if(length(pars) != length(LB))
			stop("\nsolnp-->error: LB length not equal to parameter length\n", call. = FALSE)
		if(is.null(UB)) UB = rep(.Machine$double.xmax/2, length(LB))
	} else{
		LB = NULL
	}
	if(!is.null(UB)){
		if(length(pars) != length(UB))
			stop("\nsolnp-->error: UB length not equal to parameter length\n", call. = FALSE)
		if(is.null(LB)) LB = rep(-.Machine$double.xmax/2, length(UB))
	} else{
		UB = NULL
	}
	if(!is.null(UB) && any(LB > UB))
		stop("\nsolnp-->error: UB must be greater than LB\n", call. = FALSE)
	
	if(!is.null(UB) && any(LB == UB))
		warning("\nsolnp-->warning: Equal Lower/Upper Bounds Found. Consider\n
						excluding fixed parameters.\n", call. = FALSE)
	# deal with infinite values as these are not accepted by solve.QP
	if(!is.null(LB) && !any(is.finite(LB))){
		idx = which(!is.finite(LB))
		LB[idx] = sign(LB[idx])*.Machine$double.xmax/2
	}
	if(!is.null(UB) && !any(is.finite(UB))){
		idx = which(!is.finite(UB))
		UB[idx] = sign(UB[idx])*.Machine$double.xmax/2
	}	
	assign(".LB", LB, envir = .env)
	assign(".UB", UB, envir = .env)
	return(1)
}

.checkfun = function(pars, fun, .env, ...)
{
	if(!is.function(fun)) stop("\nsolnp-->error: fun does not appear to be a function\n", call. = FALSE)
	val = fun(pars, ...)
	if(length(val) != 1)  stop("\nsolnp-->error: objective function returns value of length greater than 1!\n", call. = FALSE)
	
	assign(".solnp_fun", fun, envir = .env)
	ctmp = get(".solnp_nfn", envir =  .env)
	assign(".solnp_nfn", ctmp + 1, envir = .env)
	return(val)
}

# Might eventually use this, but really the user must take care of such problems
# in their own function/setup
.safefunx = function(pars, fun, .env, ...){
	xnames = get("xnames", envir = .env)
	names(pars) = xnames
	v  = fun(pars, ...)
	if(is.na(v) | !is.finite(v) | is.nan(v)) {
		warning(paste("\nsolnp-->warning: ", v , " detected in function call...check your function\n", sep = ""), immediate. = FALSE)
		v = 1e24
	}
	v
}

.checkgrad = function(pars, fun, .env, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(val)!=n)
		stop("\nsolnp-->error: gradient vector length must be equal to length(pars)\n", call. = FALSE)
	assign(".solnp_gradfun", fun, envir = .env)
	return(val)
}

.checkhess = function(pars, fun, .env, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(as.vector(val)) != (n*n))
		stop("\nsolnp-->error: hessian must be of length length(pars) x length(pars)\n", call. = FALSE)
	assign(".solnp_hessfun", fun, envir = .env)
	return(val)
}

.checkineq = function(pars, fun, ineqLB, ineqUB, .env, ...)
{
	xnames = get("xnames", envir = .env)
	val = fun(pars, ...)
	n = length(val)
	if(!is.null(ineqLB)){
		if(length(ineqLB) != n)
			stop("\nsolnp-->error: inequality function returns vector of different length to
							inequality lower bounds\n", call. = FALSE)
	} else{
		stop("\nsolnp-->error: inequality function given without lower bounds\n", call. = FALSE)
	}
	if(!is.null(ineqUB)){
		if(length(ineqUB) != n)
			stop("\nsolnp-->error: inequality function returns vector of different length to
							inequality upper bounds\n", call. = FALSE)
	} else{
		stop("\nsolnp-->error: inequality function given without upper bounds\n", call. = FALSE)
	}
	if(any(ineqLB > ineqUB))
		stop("\nsolnp-->error: ineqUB must be greater than ineqLB\n", call. = FALSE)

	assign(".ineqLB", ineqLB, envir = .env)
	assign(".ineqUB", ineqUB, envir = .env)
	.solnp_ineqfun = function(x, ...){
		names(x) = xnames
		fun(x, ...)
	}
	assign(".solnp_ineqfun", .solnp_ineqfun, envir = .env)
	return(val)
}


.checkeq = function(pars, fun, eqB, .env, ...)
{
	xnames = get("xnames", envir = .env)
	n = length(eqB)
	val = fun(pars, ...) - eqB
	if(length(val)!=n)
		stop("\nsolnp-->error: equality function returns vector of different length
						to equality value\n", call. = FALSE)
	.eqB = eqB
	assign(".eqB", .eqB, envir = .env)
	.solnp_eqfun = function(x, ...){
		names(x) = xnames
		fun(x, ...) - .eqB
	}
	assign(".solnp_eqfun", .solnp_eqfun, envir = .env)
	return(val)
}


# check the jacobian of inequality
.cheqjacineq = function(pars, fun, .env,  ...)
{
	# must be a matrix -> nrows = no.inequalities, ncol = length(pars)
	val = fun(pars, ...)
	.ineqLB = get(".ineqLB", envir = .env)
	.ineqUB = get(".ineqUB", envir = .env)
	if(!is.matrix(val))
		stop("\nsolnp-->error: Jacobian of Inequality must return a matrix type object\n", call. = FALSE)
	nd = dim(val)
	if(nd[2] != length(pars))
		stop("\nsolnp-->error: Jacobian of Inequality column dimension must be equal to length
						of parameters\n", call. = FALSE)
	if(nd[1] != length(.ineqUB))
		stop("\nsolnp-->error: Jacobian of Inequality row dimension must be equal to length
						of inequality bounds vector\n", call. = FALSE)
	# as in inequality function, transforms from a 2 sided inequality to a one sided inequality
	# (for the jacobian).
	.solnp_ineqjac = function(x, ...) { retval = fun(x, ...); rbind( retval ) }
	assign(".solnp_ineqjac", .solnp_ineqjac, envir = .env)
	return(val)
}

# check the jacobian of equality
.cheqjaceq = function(pars, fun, .env, ...)
{
	# must be a matrix -> nrows = no.equalities, ncol = length(pars)
	val = fun(pars, ...)
	.eqB = get(".eqB", envir = .env)
	if(!is.matrix(val))
		stop("\nsolnp-->error: Jacobian of Equality must return a matrix type object\n", call. = FALSE)
	nd = dim(val)
	if(nd[2] != length(pars))
		stop("\nsolnp-->error: Jacobian of Equality column dimension must be equal to length of parameters\n", call. = FALSE)
	if(nd[1] != length(.eqB))
		stop("\nsolnp-->error: Jacobian of Equality row dimension must be equal to length of equality bounds vector\n", call. = FALSE)
	assign(".solnp_eqjac", fun, envir = .env)
	return(val)
}

# reporting function
.report = function(iter, funv, pars)
{
	# cat( paste( "\nIter: ", iter ," fn: ", format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE), "\t Pars: ", sep=""), 
			# format(pars, digits = 4, scientific = 6, nsmall = 5, zero.print = TRUE) )
}

# finite difference gradient
.fdgrad = function(pars, fun, ...)
{
	if(!is.null(fun)){
		
		y0 = fun(pars, ...)
		nx = length(pars)
		grd = rep(0, nx)
		deltax = sqrt(.eps)
		for(i in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			grd[i] = (fun(pars, ...) - y0) / deltax
			pars[i] = init
		}
	}
	else
	{
		grd = 0
	}
	return(grd)
}

# finite difference jacobian
.fdjac = function(pars, fun, ...)
{
	nx = length(pars)
	if(!is.null(fun))
	{
		y0 = fun(pars, ...)
		nf = length (y0)
		jac = matrix(0, nrow = nf, ncol= nx)
		deltax = sqrt (.eps)
		for(i  in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			jac[,i] = (fun(pars, ...) - y0) / deltax
			pars[i] = init
		}
	} else{
		jac = rep(0, nx)
	}
	return(jac)
}

.emptygrad = function(pars, ...)
{
	matrix(0, nrow = 0, ncol = 1)
}

.emptyjac = function(pars, ...)
{
	#matrix(0, nrow = 0, ncol = length(pars))
	NULL
}

.emptyfun = function(pars, ...)
{
	NULL
}

.ineqlbfun = function(pars, .env, ...)
{
	LB = get(".solnp_LB", envir = .env)
	UB = get(".solnp_UB", envir = .env)
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .env)
	res = c(pars - LB,  UB - pars)
	if(!is.null(.solnp_ineqfun)) res = c(.solnp_ineqfun(pars, ...), res)
	res
}

.ineqlbjac = function(pars, .env, ...)
{
	.solnp_ineqjac = get(".solnp_ineqjac", envir = .env)
	n = length(pars)
	res = rbind(diag(n), -diag(n))
	if(!is.null(.solnp_ineqjac)) res = rbind(.solnp_ineqjac(pars, ...), res)
	res
}

.solnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$rho = 1
		ans$outer.iter = 400
		ans$inner.iter = 800
		ans$delta = 1.0e-7
		ans$tol = 1.0e-8
		ans$trace = 1
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 3) == "rho")) ans$rho = as.numeric(params["rho"]) else ans$rho = 1
		if(any(substr(npar, 1, 10) == "outer.iter")) ans$outer.iter = as.numeric(params["outer.iter"]) else ans$outer.iter = 400
		if(any(substr(npar, 1, 10) == "inner.iter")) ans$inner.iter = as.numeric(params["inner.iter"]) else ans$inner.iter = 800
		if(any(substr(npar, 1, 5) == "delta")) ans$delta = as.numeric(params["delta"]) else ans$delta = 1.0e-7
		if(any(substr(npar, 1, 3) == "tol")) ans$tol = as.numeric(params["tol"]) else ans$tol = 1.0e-8
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 1
	}
	return(ans)
}

.zeros = function( n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

.ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

.vnorm = function(x)
{
	sum((x)^2)^(1/2)
}

.solvecond = function(x)
{
	z = svd(x)$d
	if(any( z == 0 )) ret = Inf else ret = max( z ) / min( z )
	return(ret)
}

#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

#---------------------------------------------------------------------------------
# optimization by randomly restarted parameters using simulated parameter strategy
# Alexios Ghalanos 2010
#---------------------------------------------------------------------------------

# allowed distributions:
# 1: uniform (no confidence in the location of the parameter...somewhere in LB-UB space)
# 2: truncnorm (high confidence in the location of the parameter)
# 3: normal (Uncertainty in Lower and Upper bounds, but some idea about the dispersion about the location)
# ...

gosolnp = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL,
		ineqUB = NULL, LB = NULL, UB = NULL, control = list(), distr = rep(1, length(LB)), distr.opt = list(),
		n.restarts = 1, n.sim = 20000, cluster = NULL, rseed = NULL, ...)
{
	if( !is.null(pars) ) gosolnp_parnames = names(pars) else gosolnp_parnames = NULL
	if(is.null(control$trace)) trace = FALSE else trace = as.logical(control$trace)
	if(is.null(control$eval.type)) parmethod = 1 else parmethod = as.integer(min(abs(control$eval.type),2))
	if(parmethod == 0) parmethod = 1
	control$eval.type = NULL
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	# function requires both upper and lower bounds
	if(is.null(LB))
		stop("\ngosolnp-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB))
		stop("\ngosolnp-->error: the function requires upper parameter bounds\n", call. = FALSE)
	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars))
		stop("\ngosolnp-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)

	np = 1:n

	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np))))
			stop("\ngosolnp-->error: fixed indices out of bounds\n", call. = FALSE)
	}
	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	# setup cluster exports:
	if( !is.null(cluster) ){
		clusterExport(cluster, c("gosolnp_parnames", "fun", "eqfun",
						"eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB"), envir = environment())
		if( !is.null(names(list(...))) ){
			# evaluate promises
			xl = names(list(...))
			for(i in 1:length(xl)){
				eval(parse(text=paste(xl[i],"=list(...)[[i]]",sep="")))
			}
			clusterExport(cluster, names(list(...)), envir = environment())
		}
		clusterEvalQ(cluster, require(Rsolnp))
	}
	# initiate random search
	gosolnp_rndpars = switch(parmethod,
			.randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
					ineqUB = ineqUB, LB = LB, UB = UB, distr = distr,
					distr.opt = distr.opt, n.restarts = n.restarts,
					n.sim = n.sim, trace = trace, rseed = rseed,
					gosolnp_parnames = gosolnp_parnames, cluster = cluster, ...),
			.randpars2(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
					ineqUB = ineqUB, LB = LB, UB = UB, distr = distr,
					distr.opt = distr.opt, n.restarts = n.restarts,
					n.sim = n.sim, rseed = rseed, trace = trace,
					gosolnp_parnames = gosolnp_parnames, cluster = cluster, ...))

	gosolnp_rndpars = gosolnp_rndpars[,1:n, drop = FALSE]
	# initiate solver restarts
	if( trace ) cat("\ngosolnp-->Starting Solver\n")
	solution = vector(mode = "list", length = n.restarts)
	if( !is.null(cluster) )
	{
		clusterExport(cluster, c("gosolnp_rndpars"), envir = environment())
		solution = parLapply(cluster, as.list(1:n.restarts), fun = function(i) {
					xx = gosolnp_rndpars[i,]
					names(xx) = gosolnp_parnames
					ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun,
									eqB = eqB, ineqfun = ineqfun,
									ineqLB = ineqLB, ineqUB = ineqUB,
									LB = LB, UB = UB,
									control = control, ...), silent = TRUE)
					if(inherits(ans, "try-error")){
						ans = list()
						ans$values = 1e10
						ans$convergence = 0
						ans$pars = rep(NA, length(xx))
					}
					return( ans )
				})
	} else {
		solution = lapply(as.list(1:n.restarts), FUN = function(i){
					xx = gosolnp_rndpars[i,]
					names(xx) = gosolnp_parnames
					ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun,
									eqB = eqB, ineqfun = ineqfun,
									ineqLB = ineqLB, ineqUB = ineqUB,
									LB = LB, UB = UB,
									control = control, ...), silent = TRUE)
					if(inherits(ans, "try-error")){
						ans = list()
						ans$values = 1e10
						ans$convergence = 0
						ans$pars = rep(NA, length(xx))
					}
					return( ans )
				})
	}
	if(n.restarts>1){
		best = sapply(solution, FUN = function(x) if(x$convergence!=0) NA else x$values[length(x$values)])
		if(all(is.na(best)))
			stop("\ngosolnp-->Could not find a feasible starting point...exiting\n", call. = FALSE)
		nb = which(best == min(best, na.rm = TRUE))[1]
		solution = solution[[nb]]
		if( trace ) cat("\ngosolnp-->Done!\n")
		solution$start.pars = gosolnp_rndpars[nb,]
		names(solution$start.pars) = gosolnp_parnames
		solution$rseed = rseed
	} else{
		solution = solution[[1]]
		solution$start.pars = gosolnp_rndpars[1,]
		names(solution$start.pars) = gosolnp_parnames
		solution$rseed = rseed
	}
	return(solution)
}



startpars = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL,
		ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL, UB = NULL,
		distr = rep(1, length(LB)), distr.opt = list(), n.sim = 20000, cluster = NULL,
		rseed = NULL, bestN = 15, eval.type = 1, trace = FALSE, ...)
{
	if( !is.null(pars) ) gosolnp_parnames = names(pars) else gosolnp_parnames = NULL
	if(is.null(eval.type)) parmethod = 1 else parmethod = as.integer(min(abs(eval.type),2))
	if(parmethod == 0) parmethod = 1
	eval.type = NULL
	#trace = FALSE
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	# function requires both upper and lower bounds
	if(is.null(LB))
		stop("\nstartpars-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB))
		stop("\nstartpars-->error: the function requires upper parameter bounds\n", call. = FALSE)

	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars))
		stop("\nstartpars-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)

	np = seq_len(n)

	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np))))
			stop("\nstartpars-->error: fixed indices out of bounds\n", call. = FALSE)
	}

	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean))
				stop(paste("\nstartpars-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd))
				stop(paste("\nstartpars-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean))
				stop(paste("\nstartpars-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd))
				stop(paste("\nstartpars-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}

	# setup cluster exports:
	if( !is.null(cluster) ){
		clusterExport(cluster, c("gosolnp_parnames", "fun", "eqfun",
						"eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB"), envir = environment())
		if( !is.null(names(list(...))) ){
			# evaluate promises
			xl = names(list(...))
			for(i in 1:length(xl)){
			  eval(parse(text = paste(xl[i], "=list(...)", "[[" , i, "]]", sep = "")))
			}
			clusterExport(cluster, names(list(...)), envir = environment())
		}
		if( !is.null(names(list(...))) ) parallel::clusterExport(cluster, names(list(...)), envir = environment())
		clusterEvalQ(cluster, require(Rsolnp))
	}

	# initiate random search
	gosolnp_rndpars = switch(parmethod,
			.randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
					LB = LB, UB = UB, distr = distr, distr.opt = distr.opt,
					n.restarts = as.integer(bestN), n.sim = n.sim, trace = trace,
					rseed = rseed, gosolnp_parnames = gosolnp_parnames,
					cluster = cluster, ...),
			.randpars2(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
					LB = LB, UB = UB, distr = distr, distr.opt = distr.opt,
					n.restarts = as.integer(bestN), n.sim = n.sim, trace = trace,
					rseed = rseed, gosolnp_parnames = gosolnp_parnames,
					cluster = cluster, ...))
	return(gosolnp_rndpars)
}


.randpars = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB,
		LB, UB, distr, distr.opt, n.restarts, n.sim, trace = TRUE, rseed,
		gosolnp_parnames, cluster, ...)
{
	if( trace ) cat("\nCalculating Random Initialization Parameters...")
	N = length(LB)
	gosolnp_rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	if(!is.null(fixed)) for(i in 1:length(fixed)) gosolnp_rndpars[,fixed[i]] = pars[fixed[i]]
	nf = 1:N
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	m = length(nf)
	set.seed(rseed)
	for(i in 1:m){
		j = nf[i]
		gosolnp_rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
		)
	}

	if( trace ) cat("ok!\n")

	if(!is.null(ineqfun)){
		if( trace ) cat("\nExcluding Inequality Violations...\n")
		ineqv = matrix(NA, ncol = length(ineqLB), nrow = n.restarts*n.sim)
		# ineqv = t(apply(rndpars, 1, FUN = function(x) ineqfun(x)))
		if(length(ineqLB) == 1){
			ineqv = apply(gosolnp_rndpars, 1, FUN = function(x){
						names(x) = gosolnp_parnames
						ineqfun(x, ...)} )
			lbviol = sum(ineqv<ineqLB)
			ubviol = sum(ineqv>ineqUB)
			if( lbviol > 0 | ubviol > 0 ){
				vidx = c(which(ineqv<ineqLB), which(ineqv>ineqUB))
				vidx = unique(vidx)
				gosolnp_rndpars = gosolnp_rndpars[-c(vidx),,drop=FALSE]
				lvx = length(vidx)
			} else{
				vidx = 0
				lvx = 0
			}
		} else{
			ineqv = t(apply(gosolnp_rndpars, 1, FUN = function(x){
								names(x) = gosolnp_parnames
								ineqfun(x, ...)} ))

			# check lower and upper violations
			lbviol = apply(ineqv, 1, FUN = function(x) sum(any(x<ineqLB)))
			ubviol = apply(ineqv, 1, FUN = function(x) sum(any(x>ineqUB)))
			if( any(lbviol > 0) | any(ubviol > 0) ){
				vidx = c(which(lbviol>0), which(ubviol>0))
				vidx = unique(vidx)
				gosolnp_rndpars = gosolnp_rndpars[-c(vidx),,drop=FALSE]
				lvx = length(vidx)

			} else{
				vidx = 0
				lvx = 0
			}
		}
		if( trace ) cat(paste("\n...Excluded ", lvx, "/",n.restarts*n.sim, " Random Sequences\n", sep = ""))
	}
	# evaluate function value
	if( trace ) cat("\nEvaluating Objective Function with Random Sampled Parameters...")
	if( !is.null(cluster) ){
		nx = dim(gosolnp_rndpars)[1]
		clusterExport(cluster, c("gosolnp_rndpars", ".safefun"), envir = environment())
		evfun = parLapply(cluster, as.list(1:nx), fun = function(i){
					.safefun(gosolnp_rndpars[i, ], fun, gosolnp_parnames, ...)
				})
		evfun = as.numeric( unlist(evfun) )
	} else{
		evfun = apply(gosolnp_rndpars, 1, FUN = function(x) .safefun(x, fun, gosolnp_parnames, ...))
	}
	if( trace ) cat("ok!\n")
	if( trace ) cat("\nSorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	ans = gosolnp_rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	if( trace ) cat("ok!\n")
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if( trace ){
		cat("\nStarting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	return(prtable)
}

# form a barrier function before passing the parameters
.randpars2 = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB, LB,
		UB, distr, distr.opt, n.restarts, n.sim, rseed, trace = TRUE,
		gosolnp_parnames, cluster, ...)
{
	if( trace ) cat("\nCalculating Random Initialization Parameters...")
	N = length(LB)
	gosolnp_idx = "a"
	gosolnp_R = NULL
	if(!is.null(ineqfun) && is.null(eqfun) ){
		gosolnp_idx = "b"
		gosolnp_R = 100
	}
	if( is.null(ineqfun) && !is.null(eqfun) ){
		gosolnp_idx = "c"
		gosolnp_R = 100
	}
	if(!is.null(ineqfun) && !is.null(eqfun) ){
		gosolnp_idx = "d"
		gosolnp_R = c(100,100)
	}
	gosolnp_rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	if(!is.null(fixed)) for(i in 1:length(fixed)) gosolnp_rndpars[,fixed[i]] = pars[fixed[i]]
	nf = 1:N
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	gosolnp_m = length(nf)
	set.seed(rseed)
	for(i in 1:gosolnp_m){
		j = nf[i]
		gosolnp_rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
		)
	}
	if( trace ) cat("ok!\n")
	# Barrier Function
	pclfn = function(x){
		z=x
		z[x<=0] = 0
		z[x>0] = (0.9+z[x>0])^2
		z
	}
	.lagrfun = function(pars, m, idx, fun, eqfun = NULL, eqB = 0, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, ...)
	{
		fn = switch(idx,
				"a" = fun(pars[1:m], ...),
				"b" = fun(pars[1:m], ...) + pars[m+1]* sum( pclfn( c(ineqLB - ineqfun(pars[1:m], ...), ineqfun(pars[1:m], ...) - ineqUB) ) ),
				"c" = fun(pars[1:m], ...) + sum( (eqfun(pars[1:m], ...) - eqB )^2 / pars[m+1]),
				"d" = fun(pars[1:m], ...) + sum( (eqfun(pars[1:m], ...) - eqB )^2 / pars[m+1]) + pars[m+2]* sum( pclfn( c(ineqLB - ineqfun(pars[1:m], ...), ineqfun(pars[1:m], ...) - ineqUB) ) ) )
		return(fn)
	}

	# evaluate function value
	if( trace ) cat("\nEvaluating Objective Function with Random Sampled Parameters...")
	if( !is.null(cluster) ){
		nx = dim(gosolnp_rndpars)[1]
		clusterExport(cluster, c("gosolnp_rndpars", "gosolnp_m", "gosolnp_idx",
						"gosolnp_R"), envir = environment())
		clusterExport(cluster, c("pclfn", ".lagrfun"), envir = environment())
		evfun = parallel::parLapply(cluster, as.list(1:nx), fun = function(i){
					.lagrfun(c(gosolnp_rndpars[i,], gosolnp_R), gosolnp_m,
							gosolnp_idx, fun, eqfun, eqB, ineqfun, ineqLB,
							ineqUB, ...)
				})
		evfun = as.numeric( unlist(evfun) )
	} else{
		evfun = apply(gosolnp_rndpars, 1, FUN = function(x){
					.lagrfun(c(x,gosolnp_R), gosolnp_m, gosolnp_idx, fun, eqfun,
							eqB, ineqfun, ineqLB, ineqUB, ...)})
	}
	if( trace ) cat("ok!\n")
	if( trace ) cat("\nSorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	#distmat = dist(evfun, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
	ans = gosolnp_rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if( trace ){
		cat("\nStarting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	return(prtable)
}


.distr1 = function(LB, UB, n)
{
	runif(n, min = LB, max = UB)
}

.distr2 = function(LB, UB, n, mean, sd)
{
	rtruncnorm(n, a = as.double(LB), b = as.double(UB), mean = as.double(mean), sd = as.double(sd))
}

.distr3 = function(n, mean, sd)
{
	rnorm(n, mean = mean, sd = sd)
}

.safefun = function(pars, fun, gosolnp_parnames, ...){
	# gosolnp_parnames = get("gosolnp_parnames", envir = .env)
	names(pars) = gosolnp_parnames
	v  = fun(pars, ...)
	if(is.na(v) | !is.finite(v) | is.nan(v)) {
		warning(paste("\ngosolnp-->warning: ", v , " detected in function call...check your function\n", sep = ""), immediate. = FALSE)
		v = 1e24
	}
	v
}

#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

benchmarkids <- function()
{
return(c("Powell", "Wright4", "Wright9", "Alkylation", "Entropy", "Box", "RosenSuzuki",
				"Electron", "Permutation"))
}


benchmark <- function( id = "Powell")
{
  if( !any(benchmarkids() == id[ 1L ]) )
    stop( "invalid benchmark id" )
 	ans = switch(id,
			Powell = .powell(),
			Wright4 = .wright4(),
			Wright9 = .wright9(),
			Alkylation = .alkylation(),
			Entropy = .entropy(),
			Box = .box(),
			RosenSuzuki = .rosensuzuki(),
			Electron = .electron(),
			Permutation = .permutation())
	return(ans)
}

.powell = function()
{
	.fn1 = function(x)
	{
		exp(x[1]*x[2]*x[3]*x[4]*x[5])
	}

	.eqn1 = function(x){
		z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
		z2=x[2]*x[3]-5*x[4]*x[5]
		z3=x[1]*x[1]*x[1]+x[2]*x[2]*x[2]
		return(c(z1,z2,z3))
	}

	.x0 = c(-2, 2, 2, -1, -1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = c(10,0,-1), control=ctrl)
	minos = list()
	minos$fn = 0.05394985
	minos$pars = c(-1.717144, 1.595710, 1.827245, 0.763643, 0.763643)
	minos$nfun = 524
	minos$iter = 12
	minos$elapsed = 0.2184

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("Powell's exponential problem is a function of five variables with
three nonlinear equality constraints on the variables.")

	return(bt)
}

.wright4 = function()
{
	.fn1 = function(x)
	{
		(x[1]-1)^2+(x[1]-x[2])^2+(x[2]-x[3])^3+(x[3]-x[4])^4+(x[4]-x[5])^4
	}

	.eqn1 = function(x){
		z1=x[1]+x[2]*x[2]+x[3]*x[3]*x[3]
		z2=x[2]-x[3]*x[3]+x[4]
		z3=x[1]*x[5]
		return(c(z1,z2,z3))
	}

	.x0 = c(1, 1, 1, 1, 1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = c(2+3*sqrt(2),-2+2*sqrt(2),2), control=ctrl)
	minos = list()
	minos$fn = 0.02931083
	minos$pars = c(1.116635, 1.220442, 1.537785, 1.972769, 1.791096)
	minos$nfun = 560
	minos$iter = 9
	minos$elapsed = 0.249

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("Wright's fourth problem is a function of five variables with
three non linear equality constraints on the variables. This popular
test problem has several local solutions and taken from Wright (1976).")
	return(bt)
}

.wright9 = function()
{
	.fn1 = function(x)
	{
		10*x[1]*x[4]-6*x[3]*x[2]*x[2]+x[2]*(x[1]*x[1]*x[1])+
				9*sin(x[5]-x[3])+x[5]^4*x[4]*x[4]*x[2]*x[2]*x[2]
	}

	.ineqn1 = function(x){
		z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
		z2=x[1]*x[1]*x[3]-x[4]*x[5]
		z3=x[2]*x[2]*x[4]+10*x[1]*x[5]
		return(c(z1,z2,z3))
	}
	ineqLB = c(-100, -2, 5)
	ineqUB = c(20, 100, 100)
	.x0 = c(1, 1, 1, 1, 1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, ineqfun = .ineqn1, ineqLB = ineqLB, ineqUB = ineqUB, control=ctrl)
	minos = list()
	minos$fn = -210.4078
	minos$pars = c(-0.08145219, 3.69237756, 2.48741102,  0.37713392, 0.17398257)
	minos$nfun = 794
	minos$iter = 11
	minos$elapsed = 0.281

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("Wright's ninth problem is a function of five variables with three
non linear inequality constraints on the variables. This popular test
problem has several local solutions and taken from Wright (1976).")
	return(bt)
}

.alkylation = function()
{
	.fn1 = function(x)
	{
		-0.63*x[4]*x[7]+50.4*x[1]+3.5*x[2]+x[3]+33.6*x[5]
	}

	.eqn1 = function(x){
		z1=98*x[3]-0.1*x[4]*x[6]*x[9]-x[3]*x[6]
		z2=1000*x[2]+100*x[5]-100*x[1]*x[8]
		z3=122*x[4]-100*x[1]-100*x[5]
		return(c(z1,z2,z3))
	}
	.ineqn1 = function(x){
		z1=(1.12*x[1]+0.13167*x[1]*x[8]-0.00667*x[1]*x[8]*x[8])/x[4]
		z2=(1.098*x[8]-0.038*x[8]*x[8]+0.325*x[6]+57.25)/x[7]
		z3=(-0.222*x[10]+35.82)/x[9]
		z4=(3*x[7]-133)/x[10]
		return(c(z1,z2,z3,z4))
	}
	ineqLB = c(0.99,0.99,0.9,0.99)
	ineqUB = c(100/99,100/99,10/9,100/99)
	eqB = c(0,0,0)
	LB = c(0,0,0,10,0,85,10,3,1,145)
	UB = c(20,16,120,50,20,93,95,12,4,162)
	.x0 = c(17.45,12,110,30,19.74,89.2,92.8,8,3.6,155)
	ctrl = list(rho = 0, trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = eqB, ineqfun = .ineqn1, ineqLB = ineqLB,
			ineqUB = ineqUB, LB = LB, UB = UB, control = ctrl)
	minos = list()
	minos$fn = -172.642
	minos$pars = c(16.996427, 16.000000, 57.685751, 30.324940, 20.000000, 90.565147, 95.000000, 10.590461, 1.561636, 153.535354)
	minos$nfun = 2587
	minos$iter = 13
	minos$elapsed = 0.811

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The Alkylation problem models a simplified alkylation process. It
is a function of ten variables with four non linear inequality and
three non linear equality constraints as well as variable bounds.
The problem is taken from Locke and Westerberg (1980).")
	return(bt)
}

.entropy = function()
{
	.fn1 = function(x)
	{
		m = length(x)
		f = 0
		for(i in 1:m){
			f = f-log(x[i])
		}
		ans = f-log(.vnorm(x-1) + 0.1)
		ans
	}

	.eqn1 = function(x){
		sum(x)
	}
	eqB = 10
	LB = rep(0,10)
	UB = rep(1000,10)
	.x0 = runif(10, 0, 1000)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = eqB, LB = LB, UB = UB, control=ctrl)
	minos = list()
	minos$fn = 0.1854782
	minos$pars = c(2.2801555, 0.8577605, 0.8577605, 0.8577605, 0.8577605, 0.8577605,
			0.8577605, 0.8577605, 0.8577605, 0.8577605)
	minos$nfun = 886
	minos$iter = 4
	minos$elapsed = 0.296

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The Entropy problem is non convex in n variables with one linear
equality constraint and variable positivity bounds.")
	return(bt)
}

.box = function()
{
	.fn1 = function(x)
	{
		-x[1]*x[2]*x[3]
	}

	.eqn1 = function(x){
		4*x[1]*x[2]+2*x[2]*x[3]+2*x[3]*x[1]
	}

	eqB = 100
	LB = rep(1, 3)
	UB = rep(10, 3)

	.x0 = c(1.1, 1.1, 9)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, eqfun = .eqn1, eqB = eqB, LB = LB, UB = UB, control=ctrl)
	minos = list()
	minos$fn = -48.11252
	minos$pars = c(2.886751, 2.886751, 5.773503)
	minos$nfun = 394
	minos$iter = 9
	minos$elapsed = 0.156

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The box problem is a function of three variables with one non
linear equality constraint and variable bounds.")
	return(bt)
}

.rosensuzuki = function()
{
	.fn1 = function(x)
	{
		x[1]*x[1]+x[2]*x[2]+2*x[3]*x[3]+x[4]*x[4]-5*x[1]-5*x[2]-21*x[3]+7*x[4]
	}

	.ineqn1 = function(x){
		z1=8-x[1]*x[1]-x[2]*x[2]-x[3]*x[3]-x[4]*x[4]-x[1]+x[2]-x[3]+x[4]
		z2=10-x[1]*x[1]-2*x[2]*x[2]-x[3]*x[3]-2*x[4]*x[4]+x[1]+x[4]
		z3=5-2*x[1]*x[1]-x[2]*x[2]-x[3]*x[3]-2*x[1]+x[2]+x[4]
		return(c(z1,z2,z3))
	}
	ineqLB = rep(0, 3)
	ineqUB = rep(1000, 3)
	.x0 = c(1, 1, 1, 1)
	ctrl=list(trace=0)
	ans = solnp(.x0, fun = .fn1, ineqfun = .ineqn1, ineqLB = ineqLB, ineqUB = ineqUB, control=ctrl)
	minos = list()
	minos$fn = -44
	minos$pars = c(2.502771e-07, 9.999997e-01, 2.000000e+00, -1.000000e+00)
	minos$nfun = 527
	minos$iter = 12
	minos$elapsed = 0.203

	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			minos =  rbind(round(minos$fn, 5L),
					round(minos$iter, 0L),
					round(0, 0L),
					round(minos$nfun, 0L),
					round(minos$elapsed, 3L),
					matrix(round(minos$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	attr(bt, "description") = paste("The Rosen-Suzuki problem is a function of four variables with
three nonlinear inequality constraints on the variables. It is taken
from Problem 43 of Hock and Schittkowski (1981).")
	return(bt)
}



#----------------------------------------------------------------------------------
# Some Problems in Global Optimization
#----------------------------------------------------------------------------------


# Distribution of Electrons on a Sphere
# Given n electrons, find the equilibrium state distribution (of minimal Coulomb potential)
# of the electrons positioned on a conducting sphere. This model is from the COPS benchmarking suite.
# See http://www-unix.mcs.anl.gov/~more/cops/.

.electron = function()
{
	gofn = function(dat, n)
	{

		x = dat[1:n]
		y = dat[(n+1):(2*n)]
		z = dat[(2*n+1):(3*n)]
		ii = matrix(1:n, ncol = n, nrow = n, byrow = TRUE)
		jj = matrix(1:n, ncol = n, nrow = n)
		ij = which(ii<jj, arr.ind = TRUE)
		i = ij[,1]
		j = ij[,2]
		#  Coulomb potential
		potential = sum(1.0/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2))
		potential
	}

	goeqfn = function(dat, n)
	{
		x = dat[1:n]
		y = dat[(n+1):(2*n)]
		z = dat[(2*n+1):(3*n)]
		apply(cbind(x^2, y^2, z^2), 1, "sum")
	}

	n = 25
	LB = rep(-1, 3*n)
	UB = rep(1, 3*n)
	eqB = rep(1, n)
	ans = gosolnp(pars  = NULL, fixed = NULL, fun = gofn, eqfun = goeqfn, eqB = eqB, LB = LB, UB = UB,
			control = list(), distr = rep(1, length(LB)), distr.opt = list(outer.iter = 10, trace = 1),
			n.restarts = 2, n.sim = 20000, rseed = 443, n = 25)

	conopt = list()
	conopt$fn  = 243.813
	conopt$iter = 33
	conopt$nfun = NA
	conopt$elapsed = 0.041
	conopt$pars = c(-0.0117133872042326,	0.627138691757704,	-0.471025867741051,	-0.164419761338935,	-0.0315460712487934,
			-0.12718981058582,	-0.540049624346613,	0.600346770449059,	0.29796281847713,	-0.740960572770077,	0.972512148478245,
			-0.870858895858346,	0.84178885636396,	-0.182471994739506,	0.603293664844919,	0.0834172554171806,	0.51317309921937,
			0.260639996237799,	-0.0972877803105543,	-0.979882559381314,	-0.64809471648373,	-0.722351411610064,	0.847184430059647,
			0.514683899757428,	-0.574607272207711, 0.114609211815613,	-0.748168886860133,	-0.379763612890494,	0.743271243797936,
			0.846784034469756,	0.220425955966718,	0.839147591392778,	-0.613810163641104,	-0.499794531840362,	-0.199680552951248,
			-0.105141937435843,	-0.434753057357539,	-0.127562956191463,	-0.895691740627038,	0.574257349984438,	-0.967631920158332,
			-0.0243647398873149,	0.959445715727407,	-0.517406241891138,	0.197677956191858,	0.503867654081605,	0.286619450971711,
			0.522509289901788,	0.474911361600545,	-0.768915699421274, -0.993341595387613,	-0.21665728244143,	-0.796187308516752,	-0.648470508368975,
			0.531000606737778,	0.967075565826839,	-0.0646353084836963,	0.512660548728168,	-0.813279524362711,	0.641174786133487,
			0.207762590603952,	-0.229335044471304,	-0.524518077390237,	-0.405512363446903,	0.55341236880543,	0.23818066376044,
			0.857939234263016,	-0.107381147942665,	0.850191665834437,	0.0274516931378701,	0.571043453369507,	-0.629331175510656,
			-0.0962423162171412,	-0.713834491989006,	0.280348229724225)
	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			conopt =  rbind(round(conopt$fn, 5L),
					round(conopt$iter, 0L),
					round(0, 0L),
					round(conopt$nfun, 0L),
					round(conopt$elapsed, 3L),
					matrix(round(conopt$pars, 5L), ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	colnames(bt) = c("solnp", "conopt")
	attr(bt, "description") = paste("The equilibrium state distribution (of minimal Coulomb potential)\n of the electrons positioned on a conducting sphere.")
	return(bt)
}

# Permutation Problem -- Unique Solution f(x) = 0 and x(i) = i

.permutation = function()
{
	.perm = function(x, n, b){
		F = 0
		for(k in 1:n){
			S = 0
			for(i in 1:n){
				S = S + ( ( (i^k) + b ) * (( x[i]/i )^k -1))
			}
			F = F + S^2
		}
		F
	}

	ans = gosolnp(pars  = NULL, fixed = NULL, fun = .perm, eqfun = NULL, eqB = NULL, LB = rep(-4, 4), UB = rep(4, 4),
			control = list(outer.iter = 25, trace = 1, tol = 1e-9), distr = rep(1, 4), distr.opt = list(),
			n.restarts = 6, n.sim = 20000, rseed = 99, n = 4, b =0.5)


	bt = data.frame( solnp = rbind(round(ans$values[length(ans$values)], 5L),
					round(ans$outer.iter, 0L),
					round(ans$convergence, 0L),
					round(ans$nfuneval, 0L),
					round(ans$elapsed, 3L),
					matrix(round(ans$pars, 5L), ncol = 1L)),
			actual =  rbind(0,
					NA,
					round(0, 0L),
					NA,
					NA,
					matrix(1:4, ncol = 1L)) )
	rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "nfunEval", "time(sec)",
			paste("par.", 1L:length(ans$pars), sep = "") )
	colnames(bt) = c("solnp", "expected")
	attr(bt, "description") = paste("Permutation Problem PERM(4,0.5).")

}