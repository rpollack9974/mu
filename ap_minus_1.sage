"""Looks at T_p acting on S_k(chi) for (p,k) an irregular pair, and checks to see if there is any 
eigenvalue defined over Z_p which is congruent to 1 mod p^2 -- same as above but computes MS over Q and then coerces to Qp"""
def ap_minus_1(chi,minp,maxp):
	ps = prime_range(minp,maxp)
	for j in range(0,len(ps)):
		p = ps[j]
		if chi.conductor() % p != 0:
			v = [(k,chi.bernoulli(k)) for k in range(1,p) if (-1)^k == chi(-1)]
			for j in range(0,len(v)):
				if v[j][1]%p==0:
					k=v[j][0]
					print "working with ",(p,k)
					M=ModularSymbols(chi,k,1).cuspidal_subspace()
					f = M.hecke_polynomial(p)
					fs = f.factor()
					S.<y>=PolynomialRing(Qp(p,10))
					for g in fs:
						R=S(g[0]).roots()
						for r in R:
							s = r[0]
							e = (s-1).valuation()
							if e>0:
								print "There is a root congruent to 1 modulo p ^",e
					del(M)
