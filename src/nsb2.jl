psi(x) = polygamma(0,x)

function find_nsb_entropy(n::Array{Int64,1}, K::Integer, precision::Real, qfun::Symbol)
	kx = StatsBase.counts(n, 1:maximum(n))
	nx = [1:maximum(n)][kx.!=0]

	edges =[precision, log(K) - precision]
	xi_lim = zeros(1,2)
end

function xi_KB(K::Int64,B::Real)
	if all(size(K)==size(B)) | length(K) == 1 | length(B) == 1
		xi = psi(B+1) - psi(1+B./K)
	else
		ArgumentError("Dimension mismatch")
	end
	return xi
end

function dxi_KB(K::Int64, B::Real)
	if all(size(K)==size(B)) | length(K) == 1 | length(B) == 1
		dxi = polygamma(1,B+1) - 1/K.*polygamma(1,1+B./K)
	else
		ArgumentError("Dimension mismatch")
	end
	return dxi
end


function max_evidence(n::Array{Int64,1}, K::Int64, precision::Real)
	kx = StatsBase.counts(n, 1:maximum(n))
	nx = [1:maximum(n)][kx.!=0]
	kx = kx[kx.!=0]
	return max_evidence(kx,nx,K,precision)
end

function get_coefficients(N::Int64)
	N2 = N*N
	N3 = N2*N
	N4 = N3*N
	N5 = N4*N
	N6 = N5*N
	N7 = N6*N
	N8 = N7*N
	N9 = N8*N
	N10 = N9*N
	N11 = N10*N
	N12 = N11*N

	Nf = float(N)
	Nm = Nf-1.0
	overN = 1/N
	Nm2 = Nm*Nm
	Nm3 = Nm2*Nm
	Nm4 = Nm3*Nm
	Nm5 = Nm4*Nm
	Nm6 = Nm5*Nm
	Nm7 = Nm6*Nm
	Nm8 = Nm7*Nm
	Nm9 = Nm8*Nm
	Nm10 = Nm9*Nm

	b = [(Nm)/(2*N), #b1
        (-2+1/N)/3, #b2
        (Nf^2 - N - 2)/(9*(Nf^2 - N)), #b3
	    (2*(2 - 3*N - 3*Nf^2 + 2*Nf^3))/(135*Nm^2*N), #b4
	    (4.0*(22.0 + 13.0*N - 12*N2 - 2.0*N3 + N4))/(405.0*Nm3*N), #b5
	    (4.0*(-40.0 + 58.0*N + 65.0*N2 - 40.0*N3 - 5.0*N4 +	
	      2.0*N5))/(1701.0*Nm4*N), #b6
	    (4*(-9496 - 6912*N + 5772*Nf^2 + 2251*Nf^3 - 
	      1053*Nf^4 - 87*Nf^5 + 29*Nf^6))/(42525*Nm^5*N), #b7
	    (16.0*(764.0 - 1030.0*N - 1434.0*Nf^2 + 757.0*Nf^3 + 295.0*Nf^4 
	       - 111.0*Nf^5 - 7.0*Nf^6 + 2.0*Nf^7))/(18225.0*Nm^6*N), #b8
	    (16*(167000 + 142516*N - 108124*Nf^2 - 66284*Nf^3 
	       + 26921*Nf^4 + 7384*Nf^5 - 2326*Nf^6 - 116*Nf^7 
	       + 29*Nf^8))/(382725*Nm^7*N), #b9
	    (16.0*(-17886224.0 + 22513608.0*N + 37376676.0*Nf^2 
	       - 17041380.0*Nf^3 - 11384883.0*Nf^4 + 3698262.0*Nf^5 
	       + 846930.0*Nf^6 - 229464.0*Nf^7 - 9387.0*Nf^8 + 2086.0*Nf^9))/
			(37889775.0*Nm^8*N), #b10
	    (16.0*(-4166651072.0 - 3997913072.0*Nf + 2783482560.0*Nf^2 +
	       2290151964.0*Nf^3 - 803439834.0*Nf^4 - 395614251.0*Nf^5 +
	       108055443.0*Nf^6 + 20215218.0*Nf^7 - 4805712.0*Nf^8 - 165395.0*Nf^9
	       + 33079.0*Nf^10))/(795685275.0*Nm^9*N), #b11
	    (32.0*(52543486208.0 - 62328059360.0*Nf - 118489458160.0*Nf^2 +
	       47185442088.0*Nf^3 + 44875379190.0*Nf^4 - 12359832987.0*Nf^5 -
	       5400540075.0*Nf^6 + 1272974916.0*Nf^7 + 200644800.0*Nf^8 - 
	       42495955.0*Nf^9 - 1255067.0*Nf^10 +
	       228194.0*Nf^11))/(14105329875.0*Nm^10*N)] #b12
   return b
end


function max_evidence(kx::Array{Int64,1}, nx::Array{Int64,1}, K::Int64, precision::Real)
	max_counter = 200
	errcode = 0

	if any(nx .<= 0)
		ArgumentError("Bins with zero count encountered")
	end
	N = sum(kx.*nx)
	if N <= 1
		ArgumentError("Too few data samples: N = $(N)")
	end
	
	ng1 = find(nx.>1)
	K1 = sum(kx)
	K2 = sum(kx[ng1])

	ep = (N-K1)/N

	if K1 >= N
		Bcl = Inf
		xicl = Inf
		errcode = 1
		dxi = NaN
		warn("MAX_EVIDENCE: No coincidence")
	elseif K1 <= 1
		Bcl = 0.0
		xicl = 0
		errcode = 2
		dxi = NaN
		warn("MAX EVIDENCE: All data coincidence")
	else
		order = 10
		b = get_coefficients(N)
		
		B0ep = N*sum(ep.^[-1:order].*b)
		println("B0ep = $(B0ep)")
		if B0ep < 0
			B0ep = precision
			warn("MAX_EVIDENCE: Series expansion for B_0 did not converge")
		end

		B0 = B0ep
		dB0 = 99999999
		counter = 0
		while ~((abs(dB0) < abs(B0*precision)) | (counter .> max_counter))
			counter += 1
			F = K1/B0 + psi(B0) - psi(B0+N)
			dF = -K1/B0^2 + polygamma(1,B0) - polygamma(1,B0+N)
			dB0 = -F/dF
			B0 += dB0
		end
		println("B0 = $(B0)")

		if counter > max_counter
			errcode = 3
			warn("MAX_EVIDENCE: Newton-Rhapson search for B_0 did not converge after $counter iterations")
		end

		Bcl = B0
		order_K = 4
		B = zeros(1,order_K)

		EG = -psi(1) 
		pg1B0 = polygamma(1,B0)
		pg1NB0 = polygamma(1,N+B0)
		denum = K1/B0^2 - pg1B0 + pg1NB0
		pg2B0 = polygamma(2,B0)
		pg2NB0 = polygamma(2,N+B0)
		pg21 = polygamma(2,1)
		pg3B0 = polygamma(3,B0)
		pg3NB0 = polygamma(4,N+B0)
		pg4B0 = polygamma(4,B0)
		pg4NB0 = polygamma(3,N+B0)

		nxng1 = nx[ng1]
		kxng1 = kx[ng1]
		skxng1 = sum(kx[ng1])
		f0 = sum(kxng1.*polygamma(0,nxng1))
		d1f0 = sum(kxng1.*polygamma(1,nxng1))
		d2f0 = sum(kxng1.*polygamma(2,nxng1))
		d3f0 = sum(kxng1.*polygamma(3,nxng1))

		#println("f0 = $f0, d1f0 = $d1f0, d2f0 = $d2f0, d3f0 = $d3f0")
		
		B[1] = (B0^2*(EG*K2 + f0)) / (B0^2*denum)
		B[2] = (K2*pi^2*B0 - (6*K1*B[1]^2)/B0^3 - 3*B[1]^2*pg2B0 +
				3*B[1]^2*pg2NB0 - 6*B0*d1f0)/(-6*denum);
		B[3] = (K2*pi^2*B[1] + (6*K1*B[1]^3)/B0^4 -(12*K1*B[1]*B[2])/B0^3 + 
				3*K2*B0^2*pg21 - 6*B[1]*B[2]*pg2B0 + 6*B[1]*B[2]*pg2NB0 - 
				B[1]^3*pg3B0 + B[1]^3*pg3NB0 - 6*B[1]*d1f0 - 3*B0^2*d2f0)/ (-6*denum)
		B[4] = -(-(K2*pi^4*B0^3)/90 + (K1*B[1]^4)/B0^5 - (K2*pi^2*B[2])/6 - 
				(3*K1*B[1]^2*B[2])/B0^4 + (K1*B[2]^2)/B0^3 + 
				(2*K1*B[1]*B[3])/B0^3 - K2*B0*B[1]*pg21 + ((B[2]^2 + 
				2*B[1]*B[3])*pg2B0)/2 
				- ((B[2]^2 + 2*B[1]*B[3])*pg2NB0)/2 + 
				(B[1]^2*B[2]*pg3B0)/2 - (B[1]^2*B[2]*pg3NB0)/2 +
				(B[1]^4*pg4B0)/ 24 - (B[1]^4*pg4NB0)/24 +  B[2]*d1f0 +
				B0*B[1]*d2f0 + (B0^3*d3f0)/6)/(-denum)

		Φ = sum(B'.*((K+0.0).^(-[1:order_K])))
		Bcl += Φ
		#println("B = $B")
		#println("Bcl = $Bcl")
		#println("Φ = $Φ")
		#ccorrect until here

		counter = 0
		dBcl = 999999999999
		F = 0.0
		dF = 0.0
		while ~((abs(dBcl) < abs(Bcl*precision)) | (counter > max_counter))
			counter += 1
			F = 1/K*sum(kx[ng1].*polygamma(0,nxng1 + Bcl/K)) - K2/K*psi(1+Bcl/K) + 
				K1/Bcl + psi(Bcl) - psi(Bcl+N)

			dF =  1/(K^2)*sum(kx[ng1].*polygamma(1,nxng1 + Bcl/K)) - K2/(K^2)*polygamma(1,1+Bcl/K)
					- K1/Bcl^2 + polygamma(1,Bcl) - polygamma(1, Bcl+N)

			dBcl = -F/dF
			if dBcl < 0.0
				dBcl = 0.0
				warn("MAX_EVIDENCE: negative (unstable) change")
				errcode = 4
			end
			Bcl = Bcl + dBcl
		end
		println("Bcl = $Bcl")

		if counter > max_counter
			errcode = 3
			warn("Newton-Raphson search for Bcl did not converge after $counter iterations")
		end

		if errcode == 3 && counter < max_counter
			println("MAX_EVIDENCE. Recovered from previous error in Bcl determination")
		end

		dBcl = 1/K^2*sum(kxng1.*polygamma(1,nxng1 + Bcl/K)) - 
				K2/K^2*polygamma(1,1+Bcl/K) - K1/Bcl^2 + polygamma(1,Bcl) - 
				polygamma(1,Bcl + N)


		xicl = xi_KB(K,Bcl)
		println(xicl)
		println(dBcl/dxi_KB(K,Bcl))
		dxi = 1/sqrt(-dBcl/dxi_KB(K,Bcl)^2)
		#FIXME: These numbers don't quite match with octave yet
		
		return Bcl, xicl, dxi, errcode
	end

end

