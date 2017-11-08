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

R.<q>=PolynomialRing(QQ)
E = -1/24 + sum([sigma(n)*q^n for n in range(1,100)])
E = E - 5*E.substitute(q^5)
E1 = E - 11*E.substitute(q^11)
E2 = E - E.substitute(q^11)
X11 = EllipticCurve('11a')
f = sum([X11.an(n)*q^n for n in range(1,100)])
S.<y> = PolynomialRing(pAdicField(5,10))
alpha = (y^2-y+5).roots()[0][0]
alpha = ZZ(alpha)
fa = f - QQ(5/alpha)*f.substitute(q^5)

E1 = E1 + O(q^100)
E2 = E2 + O(q^100)
fa = fa + O(q^100)

b1 = (E1-E2)/5 
b2 = E1
b3 = (E1-fa)/5