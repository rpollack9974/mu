"""Looks at T_p acting on S_k(chi) for (p,k) an irregular pair, and checks to see if there is any 
eigenvalue defined over Z_p which is congruent to 1 mod p^2 -- same as above but computes MS over Q and then coerces to Qp"""
def ap_minus_1(chi,minp,maxp):
	filename = "val_ap_minus_1.N="+str(chi.conductor())+".txt"
	F = open(filename,'a')
	F.write("Working with character of conductor "+str(chi.conductor())+"\n"+"\n")
	F.close()
	ps = prime_range(minp,maxp)
	for j in range(0,len(ps)):
		p = ps[j]
		if chi.conductor() % p != 0:
			v = [(k,chi.bernoulli(k)) for k in range(1,p) if (-1)^k == chi(-1)]
			for j in range(0,len(v)):
				if v[j][1]%p==0:
					k=v[j][0]
					print "Working with ",(p,k)
					F = open(filename,'a')
					F.write("Working with "+str((p,k)))
					F.close()					
					M=ModularSymbols(chi,k,1).cuspidal_subspace()
					f = M.hecke_polynomial(p)
					fs = f.factor()
					S.<y>=PolynomialRing(Qp(p,10))
					for g in fs:
						pol = S(g[0]).quo_rem(y^(S(g[0]).valuation()))[0]
						if pol.degree()>0:
							R=pol.roots()
							for r in R:
								s = r[0]
								e = (s-1).valuation()
								if e>0:
									print " There is a root congruent to 1 modulo p ^",e
									F = open(filename,'a')
									F.write(" --- there is a root congruent to 1 modulo p^"+str(e)+"\n")
									F.close()
					del(M)

##Tries to see if the cuspidal Hida algebra at eisenstein ideal has rank 1
## need to enter chi a Z-valued character
def rank1(chi,minp,maxp):
	filename = "rank1_data.N="+str(chi.conductor())+".txt"
	ps=prime_range(minp,maxp)
	for j in range(0,len(ps)):
		p=ps[j]
		for k in range(1,p):
			if chi.conductor()>1:
				Bk = chi.base_extend(QQ).bernoulli(k)
			else:
				Bk = bernoulli(k)
			if Bk!=0 and Bk.valuation(p)>0:
				print "Working with ",(p,k)
				F = open(filename,'a')
				F.write("Working with "+str((p,k)))
				F.close()
				if chi.conductor()>1:
					chip = chi.base_extend(GF(p))
					M=ModularSymbols(chip,k,1,GF(p)).cuspidal_subspace()
				else:
					chip = chi
					M=ModularSymbols(1,k,1,GF(p)).cuspidal_subspace()
				q=2
				#print M.hecke_polynomial(q).substitute(q^(irr-1)+1)
				#print M.hecke_polynomial(q).derivative().substitute(q^(irr-1)+1)
				if chi.conductor()>1:
					while (M.hecke_polynomial(q).derivative().substitute(chip(q)*q^(k-1)+1)==0) and (q<20):
						print "  FAILED AT ",q
						q=next_prime(q)
					if q<20:
						F = open(filename,'a')
						F.write(" --- rank 1\n")
						F.close()
				else:
					while (M.hecke_polynomial(q).derivative().substitute(q^(k-1)+1)==0) and (q<20):
						print "  FAILED AT ",q
						q=next_prime(q)
					if q<20:
						F = open(filename,'a')
						F.write(" --- rank 1\n")
						F.close()
				del(M)
