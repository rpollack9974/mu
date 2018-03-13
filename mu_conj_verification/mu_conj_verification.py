def eisen_mult(p,k,chi=None):
	"""returns the dimension of the generalized eisenstein eigenspace 

	that is	where T_ell acts by 1 + chi(ell) ell^{k-1} and U_q acts by 1 for q | Np.
	"""
	if chi == None:
		M = ModularSymbols(1,k,1,GF(p)).cuspidal_subspace()
		N = 1
	else:
		M = ModularSymbols(chi,k,1,GF(p)).cuspidal_subspace()
		N = chi.modulus()

	d = M.dimension()
	q = 2
	sturm = M.sturm_bound()
	while d > 1 and q <= sturm:
		T = M.hecke_operator(q)
		if gcd(q,N*p) == 1:
			if chi == None:
				aq = 1 + q**(k-1)
			else:
				aq = 1 + chi(q) * q**(k-1)
		else:
			aq = 1
		aq = GF(p)(aq)
		M = ((T-aq)**d).kernel()
		d = M.dimension()
		q = next_prime(q)
	return d


def cuspidal_eisenstein_plus_eigensymbol(p,k,chi=None,acc=3):
	assert eisen_mult(p,k,chi=chi) == 1,"not unique cuspidal eisenstein symbol"
	if chi != None:
		N = chi.modulus()
	else:
		N = 1
	Phi = random_OMS(N,p,k-2,acc,char=chi)
	for j in range(acc):
		print "projecting to the ordinary subspace",(j,acc-1)
		Phi = Phi.hecke(p)
	if gcd(N,2) == 1:
		if chi == None:
			a2 = 1 + 2**(k-1)
		else:
			a2 = 1 + chi(2) * 2**(k-1)
	else:
		a2 = 1

	Phi = Phi.hecke(2) - Phi.scale(a2)

	if chi == None:
		M_full = ModularSymbols(N,k,1,GF(p)).cuspidal_subspace()
	else:
		M_full = ModularSymbols(chi,k,1,GF(p)).cuspidal_subspace()

	M = M_full
	d = M.dimension()
	q = 2
	while d > 1:
		T = M.hecke_operator(q)
		if gcd(q,N*p) == 1:
			if chi == None:
				aq = 1 + q**(k-1)
			else:
				aq = 1 + chi(q) * q**(k-1)
		else:
			aq = 1
		aq = GF(p)(aq)
		M = ((T-aq)**d).kernel()
		d = M.dimension()
		fq = M_full.hecke_polynomial(q)
		x = fq.parent().gen()
		while fq.substitute(x=aq) == 0:
			fq = fq.quo_rem(x-aq)[0]
		for j in range(acc):
			print "killing non-eisen part at q=",q,(j,acc-1)
			R = PolynomialRing(ZZ,'y')
			Phi = Phi.hecke_by_poly(q,R(fq))
		q = next_prime(q)
	Phi = Phi.plus_part()
	v = Phi.valuation()
	return Phi.change_precision(acc-v)

def mu(Phi,ap,r,lower):
	p = Phi.p()
	acc = Phi.num_moments()
	bool = false
	n = 0
	upper = Infinity
	while (not bool) and (n < acc):
		if n > 0:
			print "Checking",n,"-th derivative"
		c = pLfunction_coef(Phi,ap,n,r,1,1+p)
		v = c.valuation(p)
		assert c.precision_absolute() > v,"did not have enough accuracy for twist "+str(r)+" coefficient "+str(n)
		upper = min(v,upper)
		assert lower <= upper,"bounds messed up!!"
		bool = upper == lower
		n += 1
	return upper,bool

def test_conj_at_irregular_pair(p,k,chi=None,acc=3):
	if chi == None:
		filename = "mu_test_level_1"
	else:
		filename = "mu_test_"+str(chi)
	if chi == None:
		b = bernoulli(k)
	else:
		b = chi.bernoulli(k)
	f = open(filename,'a')
	f.write("\nWorking with irregular pair "+str((p,k))+'\n')
	f.write("  Valuation of B_k  : "+str(b.valuation(p))+'\n')
	f.close()
	Phi = cuspidal_eisenstein_plus_eigensymbol(p,k,chi=chi,acc=acc)
	print "Finding ap"
	ap,diff = Phi.is_Tq_eigen(p)
	print "Eigensymbol accurate modulo p^"+str(min(acc,diff))
	assert diff > b.valuation(p),"not enough accuracy"
	v = (ap-1).valuation(p)
	print "ord_p(ap-1) =",v
	f = open(filename,'a')
	f.write("  Valuation of a_p-1: "+str(v)+'\n')
	f.close()
	if chi == None:
		lower = bernoulli(k).valuation(p)
	else:
		lower = chi.bernoulli(k).valuation(p)
	passed = True
	for r in range(0,p-1,2):
		m,bool = mu(Phi,ap,r,lower)
		print "At twist",r,"mu = ",m
		if not bool:
			print "Failed at twist",r
		passed = passed and bool
	return passed

def test_conj(minp,maxp,chi=None,acc=3):
	if chi == None:
		filename = "mu_test_level_1"
	else:
		filename = "mu_test_"+str(chi)
	for p in primes(max(3,minp),maxp+1):
		if chi == None:
			bs = bernoulli_mod_p(p)
			irregs = [2*k for k in range(len(bs)) if bs[k]==0]
		else:
			irregs = []
			for k in range(2,p-1):
				b = chi.bernoulli(k)
				if b != 0 and b.valuation(p) > 0:
					irregs += [k]
		for k in irregs:
			print "Testing irregular pair:",(p,k)
			bool = test_conj_at_irregular_pair(p,k,chi=chi,acc=acc)
			print "Result for",(p,k),":",bool
			print "************************************************************"
			f = open(filename,'a')
			f.write("Result: "+str(bool)+"\n")
			f.close()
	return 