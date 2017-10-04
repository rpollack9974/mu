def hecke(f,ell,k,N):
	R = f.parent()
	q = R.gen()
	prec = f.precision_absolute()
	ans = R(0)
	coefs = f.padded_list()
	for n in range(prec):
		if n * ell < prec:
			if n % ell != 0 or N % ell == 0:
				ans += q^n * coefs[n*ell]
#				print "A) new ",n,"coef is",coefs[n*ell]
#				print ans
			else:
				ans += q^n * (coefs[n*ell] + ell^(k-1) * coefs[ZZ(n/ell)])
#				print "B) new ",n,"coef is",(coefs[n*ell] + ell^(k-1) * coefs[ZZ(n/ell)])
#				print ans
	ans += O(q^(floor(n/ell)))
	return ans