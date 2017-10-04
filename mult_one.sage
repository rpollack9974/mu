def check_mult_one(N):
	assert is_prime(N),"input must be prime"
	qs = (N-1).factor()
	for d in qs:
		p = d[0]
		if p > 3:
			print " Working with (N,p)=",(N,p)
			M = ModularSymbols(N*p,2,1)
			for q in list(primes(50))+[N,p]:
				Tq = M.hecke_operator(q)
				if gcd(q,N*p)==1:
					M = (Tq-q-1).kernel()
				else:
					M = (Tq-1).kernel()
			print "  dimension is",M.dimension()
