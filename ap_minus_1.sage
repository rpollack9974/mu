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
## HORRIBLE!!  WHY DON'T I JUST DECOMPOSE?
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
						F = open(filename,'a')
						F.write("\n  FAILED AT "+str(q))
						F.close()
						q=next_prime(q)
					F = open(filename,'a')
					if q<20:
						F.write(" --- rank 1\n")
					else:
						F.write("\n")
					F.close()
				else:
					while (M.hecke_polynomial(q).derivative().substitute(q^(k-1)+1)==0) and (q<20):
						print "  FAILED AT ",q
						F = open(filename,'a')
						F.write("\n  FAILED AT "+str(q))
						F.close()
						q=next_prime(q)
					F = open(filename,'a')
					if q<20:
						F.write(" --- rank 1\n")
					else:
						F.write("\n")
					F.close()
				del(M)

##Returns true if prod_{i=1}^{(N-1/2)} i^i mod N is not a p-th power
##This should predict that the cuspidal Hecke-algebra is rank 1
def Merel_check(p,N):
	assert ((N-1)/12).numerator()%p==0,"Not eisenstein"
	R = IntegerModRing(N)
	r = ZZ(prod([R(i)^i for i in range(1,(N+1)/2)]))
	R.<x>=PolynomialRing(ZZ)
	return len(R(x^p-r).factor_mod(N))==1

##Checks if p is a p-th power mod N
def good_check(p,N):
	R.<x>=PolynomialRing(ZZ)
	return len((x^p-p).factor_mod(N))==1

def E2sub_eigen(M):
	N=M.level()
	v = factor(N)
	Ns = [v[a][0] for a in range(len(v))]
	qs = list(primes(100)) + Ns 
	chi = trivial_character(N,GF(p))
	r=0
	while M.dimension()>1 and r<len(qs):
		q = qs[r]
		print q
		Tq = M.hecke_operator(q)
		M = ((Tq-chi(q)*q-1)).kernel()
		r = r + 1
	return M

def E2sub_gen(M,wN=None):
	N = M.level()
	v = factor(N)
	Ns = [v[a][0] for a in range(len(v))]
	qs = list(primes(M.hecke_bound())) 
	r = 0
	for q in Ns:
		if wN.count(q)==0:
			print "Using Hecke at ",q
			Uq = M.hecke_operator(q)
			M = ((Uq-1)^M.dimension()).kernel()
		else:
			print "Atkin-Lehner at ",q
			wq = M.atkin_lehner_operator(q)
			M = ((wq+1)).kernel()
		print "dim =",M.dimension()
	while M.dimension()>1 and r<len(qs):
		q = qs[r]
		if N%q!=0:
			print "Using Hecke at ",q
			Tq = M.hecke_operator(q)
			M = ((Tq-q-1)^M.dimension()).kernel()
			print "dim =",M.dimension()
		r = r + 1
	return M


"""Looks at T_p acting on S_2(Gamma_0(N)) when p | N-1 so that there is a congruence with an Eisenstein series
	This code checks to see if there is any eigenvalue defined over Z_p which is congruent to 1 mod p^2 """
def wt2_eisen(minN,maxN):
	filename = "wt2_eisen.txt"
	F = open(filename,'a')
	F.write("Running thru levels N with p | N-1")
	F.close()
	Ns = prime_range(minN,maxN+1)
	for j in range(0,len(Ns)):
		N = Ns[j]
		n = ((N-1)/12).numerator()
		ps = factor(n)
		for facts in ps:
			p = facts[0]
			if p>3:
				print "Working with ",(p,N)
				F = open(filename,'a')
				F.write("Working with "+str((p,N)))
				F.close()					

				#first in char p
				M = ModularSymbols(N*p,2,1,GF(p)).cuspidal_subspace()
				qs = [N]+[p]+list(primes(40))
				r = 0 
				chi = trivial_character(N*p,GF(p))			
				while M.dimension()>1 and r<len(qs):
					q = qs[r]
					Tq = M.hecke_operator(q)
					M = ((Tq-chi(q)*q-1)^100).kernel()
					r = r + 1
				F = open(filename,'a')
				if r < len(qs):
					print " The rank is 1"
					F.write(" --- rank 1\n")
				else:
					F.write("The rank is bounded by"+str(M.dimension()))
					print "The rank is bounded by",M.dimension()
				F.close()
				print "Merel check gives",Merel_check(p,N)
				print "Good check gives",good_check(p,N)
				del(M)

				#now char 0
				M=ModularSymbols(N,2,1).cuspidal_subspace()
				f = M.hecke_polynomial(p)
				S.<y>=PolynomialRing(Qp(p,10))
				fs = f.factor()
				for g in fs:
					pol = S(g[0]).quo_rem(y^(S(g[0]).valuation()))[0]
					if pol.degree()>0:
						R=pol.roots()
						for r in R:
							s = r[0]
							e = (s-1).valuation()
							if e>0:
								print "  There is a root congruent to 1 modulo p ^",e
								F = open(filename,'a')
								F.write(" --- there is a root congruent to 1 modulo p^"+str(e)+"\n")
								F.close()
				del(M)
				print ""

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
						F = open(filename,'a')
						F.write("\n  FAILED AT "+str(q))
						F.close()
						q=next_prime(q)
					F = open(filename,'a')
					if q<20:
						F.write(" --- rank 1\n")
					else:
						F.write("\n")
					F.close()
				else:
					while (M.hecke_polynomial(q).derivative().substitute(q^(k-1)+1)==0) and (q<20):
						print "  FAILED AT ",q
						F = open(filename,'a')
						F.write("\n  FAILED AT "+str(q))
						F.close()
						q=next_prime(q)
					F = open(filename,'a')
					if q<20:
						F.write(" --- rank 1\n")
					else:
						F.write("\n")
					F.close()
				del(M)

## this code is running through (p,N)'s with N=1 (mod p) and seeing if there is a Eisenstein p-newform
## at level Np and comparing this to whether p is a p-th power mod N
def crap(minN,maxN):
	filename = "wt2_eisen.txt"
	F = open(filename,'a')
	F.write("Running thru levels N with p | N-1")
	F.close()
	Ns = prime_range(minN,maxN+1)
	for j in range(0,len(Ns)):
		N = Ns[j]
		n = ((N-1)/12).numerator()
		ps = factor(n)
		for facts in ps:
			p = facts[0]
			if p>3:
				M = ModularSymbols(N*p,2,1,GF(p)).cuspidal_subspace()
				qs = [N]+[p]+list(primes(40))
				r = 0 
				chi = trivial_character(N*p,GF(p))			
				while M.dimension()>1 and r<len(qs):
					q = qs[r]
					Tq = M.hecke_operator(q)
					M = ((Tq-chi(q)*q-1)^100).kernel()
					r = r + 1
				dp = M.dimension()
				M = ModularSymbols(N,2,1,GF(p)).cuspidal_subspace()
				qs = [N]+[p]+list(primes(40))
				r = 0 
				chi = trivial_character(N*p,GF(p))			
				while M.dimension()>1 and r<len(qs):
					q = qs[r]
					Tq = M.hecke_operator(q)
					M = ((Tq-chi(q)*q-1)^100).kernel()
					r = r + 1
				d = M.dimension()
				if d!=dp:
					print "working with",(p,N)
					print "(dp,d)=",(dp,d)
					print "good_check=",good_check(p,N)
					print ""

########PRESTON CODE

def EisensteinRank(ells,p):
    # ells is a list of primes, and p is a prime
    # Let M=ModularSymbols(N,2,1,GF(p)) (it is non-cuspidal), and let X be the Eisenstein part of M. (with Atkin-Lehner eigenvalues all -1)
    # Returns X as a vector space
    N=prod(ells);
    M=ModularSymbols(N,2,1,GF(p));
    R.<x>=PolynomialRing(GF(p));
    d=M.dimension(); #At each step in the loop, d will be our best guess for dim(X)
    K=M.free_module(); #At each step in the loop, K will be our best guess for X
    for ell in ells:
        print "Using",ell
        if not ell==p:
            m=(M.atkin_lehner_operator(ell)).matrix()+1;
            K=m.kernel_on(K); #this is the part of K where w_ell acts as -1
        else:
            m=(M.hecke_matrix(ell))-1;
            K=m.kernel_on(K,poly=x^d); # this is the part of K where U_p acts unipotently
        print "dim =",K.dimension()
    for ell in prime_range(M.hecke_bound()+1):
        if not ell.divides(N):
            print "Using",ell
            d=K.dimension();
            m=M.hecke_matrix(ell)-ell-1; # Note that m maps K into K
            K=m.kernel_on(K,poly=x^d); # Computes the kernel of m^d as an operator on K. By definition of d, this is the part of K on which T_ell-ell-1 acts nilpotently.
        print "dim =",K.dimension()
    return K;
def EisImage(ells,p):
    # Let K=EisensteinRank(ells,p), and let I be the Eisenstein ideal (with w_ell+1 for ell in ells, except it has U_p-1 if p is in ells))
    # Returns IK as a subspace of K
    N=prod(ells);
    M=ModularSymbols(N,2,1,GF(p));
    K=EisensteinRank(ells,p);
    Im=((M.T(2)*0).matrix()).restrict(K).image() # this is a silly way to initialize Im to be the 0 subspace of K; at each step Im is our best guess for IK.
    if p.divides(N):
        m=(M.hecke_matrix(p))-1;
        mK=m.restrict(K);
        Im=Im+mK.image();
    for ell in prime_range(M.hecke_bound()+1):
        if not ell.divides(N):
            m=M.hecke_matrix(ell)-ell-1; # Note that m maps K into K
            mK=m.restrict(K);
            Im=Im+mK.image();
    return Im;
def EisKernel(ells,p):
    # Let M=ModularSymbols(N,2,1,GF(p)) (it is non-cuspidal), and let I be the Eisenstein ideal (with w_ell+1 for ell in ells, except it has U_p-1 if p is in ells)
    # Returns M[I]
    N=prod(ells);
    M=ModularSymbols(N,2,1,GF(p));
    R.<x>=PolynomialRing(GF(p));
    d=M.dimension(); 
    K=M.free_module();
    for ell in ells:
        if not ell==p:
            m=(M.atkin_lehner_operator(ell)).matrix()+1;
            K=m.kernel_on(K);
        else:
            m=(M.hecke_matrix(ell))-1;
            K=m.kernel_on(K);
    for ell in prime_range(M.hecke_bound()+1):
        if not ell.divides(N):
            m=M.hecke_matrix(ell)-ell-1;
            K=m.kernel_on(K);
    return K;