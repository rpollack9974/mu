def check_mult_one(N,max_prime=20,verbose=false):
	filename = "mu_data.txt"

	assert is_prime(N),"input must be prime"
	qs = (N-1).factor()
	for d in qs:
		p = d[0]
		if p > 3:
			print "Working with (N,p)=",(N,p)
			F = open(filename,'a')
			F.write("Working with (N,p)="+str((N,p))+" --- ")
			F.close()
			M = ModularSymbols(N*p,2,1)
			for q in list(primes(max_prime))+[N,p]:
				if verbose:
					print q
				Tq = M.hecke_operator(q)
				if gcd(q,N*p)==1:
					M = (Tq-q-1).kernel()
				else:
					M = (Tq-1).kernel()
			print "  dimension is",M.dimension()
			F = open(filename,'a')
			F.write("dimension is "+str(M.dimension())+'\n')
			F.close()
			if M.dimension() > 1:
				print "****************************************************"
				F = open(filename,'a')
				F.write("****************************************************"+'\n')
				F.close()

