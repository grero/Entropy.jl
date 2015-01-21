type NSBEntropy
	N::Int64
	K::Int64
	nx::Array{Int64,1} #unique counts
	kx::Array{Int64,1} #number of each unique count
	nc::Array{Int64,1} #conidience count
	kc::Array{Int64,1} #number of each unique coincidence count
	K1::Int64 # number of bins non-zero counts
	K2::Int64 #number of bins with coincidence
	B_interp::Array{Float64} #interpolation grid
	S_nsb::Float64 #nsb entropy
	dS_nsb::Float64 #standard devation
	S_ml::Float64 #maximum likelihood estimator
end

function NSBEntropy()
	K = 0
	nx = Array(Int64,0)
	kx = Array(Int64,0)
	nc = Array(Int64,0)
	kc = Array(Int64,0)
	B_interp = Array(Float64,0)
	S_nsb = 0.0
	dS_nsb = 0.0
	S_ml = 0.0
end

function NSBEntropy(n::Array{Int64,1},K::Int64)
	k = StatsBase.counts(n, 1:maximum(n))
	nx = [1:maximum(n)][k.!=0]
	kx = k[k.>0]
	idx = nx .> 1
	nc = nx[idx]
	kc = kx[idx]
	N = sum(nx.*kx)
	K1 = sum(kx)
	K2 = sum(kc)
	return NSBEntropy(N,K,nx,kx,nc,kc,K1, K2, Float64[], 0.0, 0.0, 0.0)
end

function NSBEntropy(nx::Array{Int64,1},kx::Array{Int64,1}, K::Int64)
	idx = nx .> 1
	nc = nx[idx]
	kc = kx[idx]
	N = sum(nx.*kx)
	K1 = sum(kx)
	K2 = sum(kc)
	return NSBEntropy(N,K,nx,kx,nc,kc,K1, K2, Float64[], 0.0, 0.0, 0.0)
end


function nsb_entropy(n::Array{Int64,1},K::Int64;precision::Real=1e-5)
	S_nsb = NSBEntropy(n,K)
	find_nsb_entropy(S_nsb,precision)
	return S_nsb.S_nsb, S_nsb.dS_nsb
end

function nsb_entropy(X::Dict{Int64,Int64})
	n = collect(values(X))
	K = maximum(keys(X))
	return nsb_entropy(n,K)
end

function nsb_entropy(X::Dict{Int64,Dict{Int64,Int64}})
    E = 0
    N = 0
    for (k,v) in X
        ee,n = nsb_entropy(v)
        E += ee*n
        N += n
    end
    return E/N,N
end

function nsb_entropy(X::Dict{Int64, Dict{Int64, Dict{Int64,Int64}}})
    E = 0
    N = 0
    for (k,v) in X
        ee,n = nsb_entropy(v)
        E += ee*n
        N += n
    end
    return E/N,N
end

function get_coincidence_count(S::NSBEntropy)
	if isempty(S.nc)
		idx = S.nx .> 1
		S.nc = S.nx[idx]
		S.kc = S.kx[dx]
	end
	return S.kc,S.nc
end

psi(x) = polygamma(0,x)
psi_asymp(x) = psi(x) - log(x) #FIXME is perhaps not as accurate as we want

function find_nsb_entropy(n::Array{Int64,1}, K::Integer, precision::Real;verbose::Int64=0)
	S_nsb = NSBEntropy(n,K)
	k = StatsBase.counts(n, 1:maximum(n))
	nx = [1:maximum(n)][k.!=0]
	kx = k[k.>0]
	return find_nsb_entropy(S_nsb,precision;verbose=verbose)
end

function find_nsb_entropy(S_nsb::NSBEntropy, precision::Real;verbose::Int64=0)
	integs=  [integrand_1,integrand_S,integrand_S2]
	msgs = ["normalization", "S","S^2"]
	K = S_nsb.K
	nsb_mlog_quad = zeros(1)
	edges =[precision, log(K) - precision]
	xi_lim = zeros(2)
	B_cl, xi_cl, dS_cl,errcode = max_evidence(S_nsb, precision;verbose=verbose)
	#println("size(xi_cl) = $(size(xi_cl))")
	#println("errocode = $errcode")
	if errcode > 0
		verbose > 0 && warn("FIND_ENTROPY: Switching integration over the full range of xi due to prior errors")
		#println("B_cl = $B_cl")
		if B_cl > 0
			nsb_mlog_quad = mlog_evidence(B_cl, S_nsb)
			S_cl = meanS(B_cl, S_nsb)
		else
			nsb_mlog = 0.0
			S_cl = NaN
		end
		xi_lim[:] = edges
		delta = 0.0
		dS_cl = NaN
		xi_cl = Inf
	else
		S_cl = meanS(B_cl, S_nsb)
		verbose > 0 && warn("FIND_NSB_ENTROPY: Expect S to be near $S_cl and σ near $dS_cl")
		#println("B_cl = $B_cl")
		nsb_mlog_quad = mlog_evidence(B_cl, S_nsb) #value at the saddle
		delta = erfinv(1-precision)*sqrt(2)
		verbose > 0 && println("FIND_NSB_ENTROPY: Integrating around peak")
	end

	nsb_val_quad = zeros(length(integs),1)
	ec = zeros(nsb_val_quad)
	nfun = zeros(nsb_val_quad)
	err = zeros(nsb_val_quad)
	for i in 1:length(integs)
		if delta > 0
			#println("size(xi_cl) = $(size(xi_cl))")
			#println("size(nsb_mlog_quad) = $(size(nsb_mlog_quad))")
			cent = integs[i](xi_cl, S_nsb,nsb_mlog_quad[1])
			good_enough = 0.1*cent*exp(-(delta^2)/2.0)
			wnd = [-1.0,-1.0]
			worst = [1.0,2.0]
			best = [1.0,2.0]
			limval = good_enough.*ones(2)
			finished = false
			while !finished
				if any(xi_lim[worst].==edges)
					if any(wnd[best].>10)
						wnd[best] *= 1.2
					else
						wnd[best] += 0.5
					end
				else
					if any(wnd[worst].>10)
						wnd[worst] *= 1.2
					else
						wnd[worst] += 0.5
					end
				end
				xi_lim[1] = max(edges[1],xi_cl - (delta+wnd[1])*dS_cl)[1]
				xi_lim[2] = min(edges[2],xi_cl + (delta+wnd[2])*dS_cl)[1]
				#println("xi_lim = $xi_lim")
				#println("limval = $limval")
				#println("dS_cl = $dS_cl")
				#println("delta = $delta")
				#println("nsb_mlog_quad = $nsb_mlog_quad")

				for j in 1:length(worst)
					limval[worst[j]] = integs[i](xi_lim[worst[j]],S_nsb, nsb_mlog_quad[1])[1]
				end
				for j in 1:length(best)
					limval[best[j]] = integs[i](xi_lim[best[j]],S_nsb,nsb_mlog_quad[1])[1]
				end
				tmp,worst = findmax(limval)
				tmp,best = findmin(limval)
				#println("best = $best")
				if any(limval[worst].<good_enough)
					finished = true
				end
				if any(limval[best] .< good_enough) && any(xi_lim[worst].==edges)
					finished = true
				end
				if xi_lim == edges
					finished = true
				end
			end
		end

		verbose >0 && println("FIND_NSB_ENTROPY: Doing $(msgs[i]) integral. Limits: $(xi_lim[1]) < xi < $(xi_lim[2])")

		_func(x) = integs[i](x, S_nsb,nsb_mlog_quad[1])
		nsb_val_quad[i,:],err[i] = quadgk(_func, xi_lim[1], xi_lim[2];reltol=precision)
		ec[i] = 0
		if ec[i] > 0
			verbose > 0 && warn("FIND_NSB_ENTROPY: Problem in $(msgs[i]) integral")
		end
		if err[i] > precision
			verbose > 0 && warn("FIND_NSB_ENTROPY: Precision of $precision required by only $(err[i]) achieved")
		end

	end

	S_nb  = nsb_val_quad[2]/nsb_val_quad[1]
	#println("S_nb = $S_nb")
	S_nsb.S_nsb = S_nb
	S_nsb.dS_nsb = sqrt(abs(nsb_val_quad[3]/nsb_val_quad[1] - S_nb^2))
	return S_nsb,err
end

function meanS{T<:Real}(B::T, kx::Array{Int64,1}, nx::Array{Int64,1}, K::Int64)
	return meanS([B],kx,nx,K)
end

function meanS{T<:Real}(B::Array{T}, kx::Array{Int64,1}, nx::Array{Int64,1}, K::Int64)
	@assert size(kx) == size(nx) 
	if any(nx.<=0)
		ArgumentError("Bins with zero count encountered")
	end
	N = sum(kx.*nx)
	Nf = N + 0.0
	if N <= 1
		verbose > 0 && warn("Too few data samples: N = $N")
	end
	K1 = sum(kx)
	_B = B[:]

	ovrNB = 1./(Nf+_B)
	osB = ones(size(_B))
	osn = ones(size(nx))
	f = zeros(size(_B))

	f = psi(Nf+_B+1.0) - (ovrNB*osn').*(osB*nx' + (_B./K)*osn').*psi(osB*nx' + (_B./K)*osn'+1.0)*kx[:] -
			_B.*ovrNB*(1.0-K1/K).*psi(1+_B/K)
	f[!isfinite(_B)] = log(K)
	f = reshape(f, size(B))
	f
end

function meanS{T<:Real}(B::Array{T,1}, S_nsb::NSBEntropy)
	f = zeros(B)
	for i in 1:length(f)
		f[i] = meanS(B[i],S_nsb)
	end
	f
end

function meanS(B::Real, S_nsb::NSBEntropy)
	N = S_nsb.N
	Nf = N + 0.0
	K = S_nsb.K
	if !isfinite(B)
		return log(K)
	end
	K1 = S_nsb.K1
	nx = S_nsb.nx
	kx = S_nsb.kx

	ovrNB = 1./(Nf+B)

	#f = psi(Nf+_B+1.0) - (ovrNB*osn').*(osB*nx' + (_B./K)*osn').*psi(osB*nx' + (_B./K)*osn'+1.0)*kx[:] -
	#		_B.*ovrNB*(1.0-K1/K).*psi(1+_B/K)
	f = psi(N + B+1.0) - (ovrNB.*(nx + (B/K)).*psi(nx + B/K + 1.0))'*kx - B*ovrNB*(1.0 - K1/K).*psi(1+B/K)
	f[1]
end

function xi_KB{T<:Real}(K::Int64, B::T)
	return psi(B+1) - psi(1+B/K)
end



function xi_KB{T<:Real}(K::Int64,B::Array{T})
	xi = zeros(B)
	for i in 1:length(xi)
		xi[i] = psi(B[i]+1) - psi(1+B[i]/K)
	end
	#xi = psi(B+1) - psi(1+B./K)
	xi
end

function dxi_KB(K::Int64, B::Real)
	dxi = polygamma(1,B+1) - 1/K.*polygamma(1,1+B./K)
	return dxi
end

function dxi_KB{T<:Real}(K::Int64, B::Array{T})
	dxi = polygamma(1,B+1) - 1/K.*polygamma(1,1+B./K)
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


function max_evidence(kx::Array{Int64,1}, nx::Array{Int64,1}, K::Int64, precision::Real;verbose::Int64=0)
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
		verbose > 0 && warn("MAX_EVIDENCE: No coincidence")
	elseif K1 <= 1
		Bcl = 0.0
		xicl = 0
		errcode = 2
		dxi = NaN
		verbose >0 && warn("MAX EVIDENCE: All data coincidence")
	else
		order = 10
		b = get_coefficients(N)

		B0ep = N*sum(ep.^[-1:order].*b)
		if B0ep < 0
			B0ep = precision
			verbose > 0 && warn("MAX_EVIDENCE: Series expansion for B_0 did not converge")
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

		if counter > max_counter
			errcode = 3
			verbose > 0 && warn("MAX_EVIDENCE: Newton-Rhapson search for B_0 did not converge after $counter iterations")
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
		while !((abs(dBcl) < abs(Bcl*precision)) || (counter > max_counter))
			counter += 1
			F = 1/K*sum(kx[ng1].*psi(nxng1 + Bcl/K)) - K2/K*psi(1+Bcl/K) +
				K1/Bcl + psi(Bcl) - psi(Bcl+N)

			dF =  1/(K^2)*sum(kx[ng1].*polygamma(1,nxng1 + Bcl/K)) - K2/(K^2)*polygamma(1,1+Bcl/K) - K1/Bcl^2 + polygamma(1,Bcl) - polygamma(1, Bcl+N)
			#println("counter = $counter, F = $F, dF = $dF")

			dBcl = -F/dF
			if dBcl < 0.0
				verbose > 0 && warn("MAX_EVIDENCE: negative (unstable) change dBcl = $dBcl")
				dBcl = 0.0
				errcode = 4
			end
			Bcl = Bcl + dBcl
		end
		#println("Bcl = $Bcl")

		if counter > max_counter
			errcode = 3
			verbose > 0 && warn("Newton-Raphson search for Bcl did not converge after $counter iterations")
		end

		if errcode == 3 && counter < max_counter
			verbose > 0 && println("MAX_EVIDENCE. Recovered from previous error in Bcl determination")
		end

		dBcl = 1/K^2*sum(kxng1.*polygamma(1,nxng1 + Bcl/K)) -
				K2/K^2*polygamma(1,1+Bcl/K) - K1/Bcl^2 + polygamma(1,Bcl) -
				polygamma(1,Bcl + N)


		xicl = xi_KB(K,Bcl)
		verbose > 0 && println("xicl = $xicl")
		k = 0
		verbose > 0 && println("dBcl = $dBcl")
		#println(K)
		#println(dxi_KB(K,Bcl))
		dxi = 1./sqrt(-dBcl./dxi_KB(K,Bcl).^2)
		#FIXME: These numbers don't quite match with octave yet
	end
	return Bcl, xicl, dxi, errcode
end

function max_evidence(S_est::NSBEntropy, precision::Real;verbose::Int64=0)
	max_counter = 200
	errcode = 0

	N = S_est.N
	K = S_est.K
	K1 = S_est.K1
	K2 = S_est.K2
	kc = S_est.kc
	nc = S_est.nc

	ep = (N-K1)/N

	if K1 >= N
		Bcl = Inf
		xicl = Inf
		errcode = 1
		dxi = NaN
		verbose > 0 && warn("MAX_EVIDENCE: No coincidence")
	elseif K1 <= 1
		Bcl = 0.0
		xicl = 0
		errcode = 2
		dxi = NaN
		verbose > 0 && warn("MAX EVIDENCE: All data coincidence")
	else
		order = 10
		b = get_coefficients(N)

		B0ep = N*sum(ep.^[-1:order].*b)
		if B0ep < 0
			B0ep = precision
			verbose > 0 && warn("MAX_EVIDENCE: Series expansion for B_0 did not converge")
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

		if counter > max_counter
			errcode = 3
			verbose >0 && warn("MAX_EVIDENCE: Newton-Rhapson search for B_0 did not converge after $counter iterations")
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
		pg3NB0 = polygamma(3,N+B0)
		pg4B0 = polygamma(4,B0)
		pg4NB0 = polygamma(4,N+B0)

		nxng1 = S_est.nc
		kxng1 = S_est.kc
		skxng1 = S_est.K2 
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
		#println("B0 = $B0") #correct
		#println("Bcl = $Bcl")
		#println("Φ = $Φ")
		#ccorrect until here

		counter = 0
		dBcl = 999999999999
		F = 0.0
		dF = 0.0
		while !((abs(dBcl) < abs(Bcl*precision)) || (counter > max_counter))
			counter += 1
			F = 1/K*sum(kc.*psi(nxng1 + Bcl/K)) - K2/K*psi(1+Bcl/K) +
				K1/Bcl + psi(Bcl) - psi(Bcl+N)

			dF =  1/(K^2)*sum(kc.*polygamma(1,nxng1 + Bcl/K)) - K2/(K^2)*polygamma(1,1+Bcl/K) - K1/Bcl^2 + polygamma(1,Bcl) - polygamma(1, Bcl+N)
			#println("counter = $counter, F = $F, dF = $dF")

			dBcl = -F/dF
			if dBcl < 0.0
				verbose > 0 && warn("MAX_EVIDENCE: negative (unstable) change dBcl = $dBcl")
				dBcl = 0.0
				errcode = 4
			end
			Bcl = Bcl + dBcl
		end
		#println("Bcl = $Bcl")

		if counter > max_counter
			errcode = 3
			verbose > 0 && warn("Newton-Raphson search for Bcl did not converge after $counter iterations")
		end

		if errcode == 3 && counter < max_counter
			verbose > 0 && println("MAX_EVIDENCE. Recovered from previous error in Bcl determination")
			errcode = 0
		end

		dBcl = 1/K^2*sum(kxng1.*polygamma(1,nxng1 + Bcl/K)) -
				K2/K^2*polygamma(1,1+Bcl/K) - K1/Bcl^2 + polygamma(1,Bcl) -
				polygamma(1,Bcl + N)


		xicl = xi_KB(K,Bcl)
		verbose > 0 && println("xicl = $xicl")
		verbose > 0 && println("dBcl = $dBcl")
		#println(K)
		#println(dxi_KB(K,Bcl))
		#println("Bcl = $Bcl")
		if errcode > 0
			dxi = 0.0
		else
			dxi = 1./sqrt(-dBcl./dxi_KB(K,Bcl).^2)
		end
		#FIXME: These numbers don't quite match with octave yet
		#println("dxi = $dxi")

	end
	return Bcl, xicl, dxi, errcode
end

gammaln(x) = log(gamma(x))

mlog_evidence(B::Real, kx::Array{Int64,1}, nx::Array{Int64,1}, K::Int64) = mlog_evidence([B],kx,nx,K)

function mlog_evidence{T<:Real}(B::Array{T,1}, kx::Array{Int64,1}, nx::Array{Int64,1}, K::Int64)
	if any(nx.<=0)
		ArgumentError("Bins with zero count encountered")
	end
	N = sum(kx.*nx)
	if N <= 1
		ArgumentError("Too few data samples: N = $N")
	end
	Nf = N + 0.0
	_B = B[:]
	f = zeros(_B)
	B0 = find(_B.==0)
	if length(B0) > 0
		Bn0 = find(_B.!=0)
		f[Bn0] = mlog_evidence(B[Bn0], kx, nx, K)
		f[B0] = Inf
	else
		ng1 = find(nx.>1)
		nxng1 = nx[ng1]
		K1 = sum(kx)
		K2 = sum(kx[ng1])
		if length(ng1) > 0 #coincidences
			f1 = K2*lgamma(1+_B/K)
			f2 = -lgamma(_B*ones(1,length(ng1)))/K + ones(length(B),1)*nx[ng1]'
			f3 = f2*kx[ng1]
			f[:] = f1 .+ f3[:]
			#for i in 1:length(f)
			#	f[i] = K2*lgamma(1+_B[i]/K)
			#	f2 = -lgamma(_B[i])/K
			#	f1 = 0.0
		#		f3 = 0.0
		#		for j1 in 1:length(ng1)
		#			f1 += nxng1[j1]
		#			f3 += nxng1[j1]^2
		#		end
		#		f[i] += f2*f1 + f3
		#	end
		end
	end
	#println(f)
	large = (_B .> max(100,100*N))
	if any(large) #we have large values
		nterms = max(ceil(abs(-15 - log10(N))./log10(N./_B[large]))) + 1
		#println("nterms = $nterms")
		ifac = 1./cumprod([2:(nterms+1)])
		for i in nterms:-1:1
			f[large]  += polygamma(i,_B[large]).*ifrac[i].*Nf^(i+1)
		end
	end
	#println(f)

	other = !large
	df = -K1*log(_B[other]) + lgamma(_B[other]+N) - lgamma(_B[other])
	#println(df)
	f[other] += -K1*log(_B[other]) + lgamma(_B[other]+N) - lgamma(_B[other])
end

function mlog_evidence{T<:Real}(B::Array{T,1}, S_nsb::NSBEntropy)
	f = zeros(size(B))
	for i in 1:length(f)
		f[i] = mlog_evidence(B[i],S_nsb)
	end
	return f
end

function mlog_evidence(B::Real, S_nsb::NSBEntropy)
	N = S_nsb.N
	Nf = N + 0.0
	f = 0.0
	if B == 0.0
		f = Inf
	else
		nxng1 = S_nsb.nc
		kxng1 = S_nsb.kc
		K = S_nsb.K
		K1 = S_nsb.K1
		K2 = S_nsb.K2
		if length(nxng1) > 0 #coincidences
			f1 = K2*lgamma(1.0+B/K)
			f2 = -lgamma(B/K + nxng1)
			f3 = f2'*kxng1
			#println("f1 = $f1")
			#println("f2 = $f2")
			#println("f3 = $f3")
			f = (f1 .+ f3[:])[1]
		end
	end
	#println("f = $f")
	#println("MLOG_EVIDENCE: f = $f")
	if B > max(100,100*N)
		nterms = int(ceil(abs((-15.0 - log10(N))./log10(N./B))))
		#println("nterms = $nterms")
		ifac = 1./cumprod([2:(nterms+1)])
		#println("ifac = $ifac")
		for i in nterms:-1:1
			f += polygamma(i,B).*ifac[i].*Nf^(i+1)
			#println("f = $f")
		end
		f += (N-K1)*log(B) + psi_asymp(B)*N
		#println("f = $f")
	else
	#println(f)
		#df = -K1*log(B) + lgamma(B+N) - lgamma(B)
		f += -K1*log(B) + lgamma(B+N) - lgamma(B)
	end
	#println("f = $f")
	return f
end

integrand_1{T<:Real}(xi::T,nsb_kx_quad::Array{Int64,1}, nsb_nx_quad::Array{Int64,1}, nsb_K_quad::Int64,nsb_mlog_quad::Real) = integrand_1([xi],nsb_kx_quad, nsb_nx_quad, nsb_K_quad,nsb_mlog_quad)

function integrand_1{T<:Real}(xi::Array{T},nsb_kx_quad::Array{Int64,1}, nsb_nx_quad::Array{Int64,1}, nsb_K_quad::Int64,nsb_mlog_quad::Real)
	B = B_xiK(xi, nsb_K_quad)
	_mle = mlog_evidence(B,nsb_kx_quad, nsb_nx_quad, nsb_K_quad)
	_pxi = prior_xi(xi, nsb_K_quad)
	f = exp( -_mle + nsb_mlog_quad).*_pxi
	#println("f = $f")
	return f
end

function integrand_1(xi::Real,S_nsb::NSBEntropy,nsb_mlog_quad::Real)
	B = B_xiK(xi,S_nsb)
	_mle = mlog_evidence(B,S_nsb)
	_pxi = prior_xi(xi, S_nsb.K)
	f = exp( -_mle + nsb_mlog_quad).*_pxi
	#println("f = $f")
	#println("_mle = $_mle")
	#println("nsb_mlog_quad = $nsb_mlog_quad")
	return f
end

integrand_S{T<:Real}(xi::T,nsb_kx_quad::Array{Int64,1}, nsb_nx_quad::Array{Int64,1}, nsb_K_quad::Int64,nsb_mlog_quad::Real) = integrand_S([xi],nsb_kx_quad, nsb_nx_quad, nsb_K_quad,nsb_mlog_quad)


function integrand_S{T<:Real}(xi::Array{T},nsb_kx_quad::Array{Int64,1}, nsb_nx_quad::Array{Int64,1}, nsb_K_quad::Int64,nsb_mlog_quad::Real)
	B = B_xiK(xi, nsb_K_quad)
	_mle = mlog_evidence(B,nsb_kx_quad, nsb_nx_quad, nsb_K_quad)
	_ms = meanS(B,nsb_kx_quad, nsb_nx_quad, nsb_K_quad)
	_pxi = prior_xi(xi, nsb_K_quad)
	f = exp( -_mle + nsb_mlog_quad).*_pxi.*_ms
	return f
end

function integrand_S(xi::Real,S_nsb::NSBEntropy,nsb_mlog_quad::Real)
	B = B_xiK(xi, S_nsb)
	_mle = mlog_evidence(B,S_nsb)
	_ms = meanS(B,S_nsb)
	_pxi = prior_xi(xi, S_nsb.K)
	f = exp( -_mle + nsb_mlog_quad).*_pxi.*_ms
	return f
end

function integrand_S2(xi::Real,S_nsb::NSBEntropy,nsb_mlog_quad::Real)
	B = B_xiK(xi, S_nsb)
	_mle = mlog_evidence(B,S_nsb)
	_ms2 = meanS2(B,S_nsb)
	_pxi = prior_xi(xi, S_nsb.K)
	f = exp( -_mle + nsb_mlog_quad).*_pxi.*_ms2
	return f
end

function meanS2(B::Real, S_nsb::NSBEntropy)

	N = S_nsb.N
	K = S_nsb.K
	K1 = S_nsb.K1
	nx = S_nsb.nx
	kx = S_nsb.kx
	b = B/K
	nxb = nx + b
	pnxb1 = psi(nxb+1.0)
	p0NB2 = psi(N + B + 2.0)
	p1NB2 = polygamma(1,N+B+2.0)
	f = sum(nxb.*((pnxb1-p0NB2).*kx)*(nxb.*(pnxb1-p0NB2).*kx)' - (nxb.*kx)*(nxb.*kx)'*p1NB2)
	f += sum(2*B*(1-K1/K)*nxb.*((pnxb1-p0NB2).*(psi(b+1)-p0NB2) - p1NB2).*kx)
	f += (1-K1/K)*(1-(K1+1)/K)*B*B*((psi(b+1) - p0NB2)^2-p1NB2)
	f += -sum((nxb.*(pnxb1 - p0NB2)).^2.*kx) + sum(nxb.*nxb.*kx)*p1NB2

	f += sum((nxb.*(nxb+1).*((psi(nxb+2) - p0NB2).^2 + polygamma(1,nxb+2) - p1NB2)).*kx)

	f += B*(1-K1/K)*(1+b)*((psi(2+b) - p0NB2).^2 + polygamma(1,b+2) - p1NB2)

	#normalizing
	f = f/((N+B)*(N+B+1))
	if !isfinite(f)
		f = log(K)^2
	end
	return f
end

function prior_xi(xi::Real, K::Int64)
	return 1.0
end

function prior_xi{T<:Real}(xi::Array{T}, K::Int64)
	f = ones(size(xi))
	return f
end

function B_xiK(K::Int64)
	#create the interpolation
	K_interp_in = K
	step = 1.0e-2
	b1 = 10.0
	b2 = 100.0*K
	b = zeros(int(b1/step+1 + ceil(log(b2)/step)))
	b[1:b1/step+1] = 0:step:b1
	b[b1/step+3:end] = b1*exp([step:step:log(b2)])
	#b = [0:step:b1,b1*exp([step:step:log(b2)])]
	Bxi_interp_in = zeros(2,length(b))
	Bxi_interp_in[1,:] = xi_KB(K,b)
	Bxi_interp_in[2,:] = b
	dxi = float([1.0e-3:log(K)*1.0e-3:log(K)-1.0e-3])
	Bxi_interp_in_new = zeros(2,length(dxi))
	Bxi_interp_in_new[1,:] = dxi
	Bxi_interp_in_new[2,:]  = B_xiK(Bxi_interp_in, dxi,K)
	return Bxi_interp_in_new
end

function B_xiK{T<:Real}(xi::Array{T}, S_nsb::NSBEntropy)
	if isempty(S_nsb.B_interp)
		S_nsb.B_interp = B_xiK(S_nsb.K)
	end
	return B_xiK(S_nsb.B_interp, xi, S_nsb.K)
end

function B_xiK(xi::Real, S_nsb::NSBEntropy)
	if isempty(S_nsb.B_interp)
		S_nsb.B_interp = B_xiK(S_nsb.K)
	end
	return B_xiK(S_nsb.B_interp, xi, S_nsb.K)
end

function B_xiK{T<:Real}(xi::Array{T}, K::Int64)
	if any(xi .> log(K))
		ArgumentError("Too large xi -- bigger than log(K)")
	end
	if any(xi .< 0 )
		ArgumentError("Too small xi -- smaller than 0")
	end
	#create the interpolation
	K_interp_in = K
	step = 1.0e-2
	b1 = 10.0
	b2 = 100.0*K
	b = zeros(int(b1/step+1 + ceil(log(b2)/step)))
	b[1:b1/step+1] = 0:step:b1
	b[b1/step+3:end] = b1*exp([step:step:log(b2)])
	#b = [0:step:b1,b1*exp([step:step:log(b2)])]
	Bxi_interp_in = zeros(2,length(b))
	Bxi_interp_in[1,:] = xi_KB(K,b)
	Bxi_interp_in[2,:] = b
	dxi = float([1.0e-3:log(K)*1.0e-3:log(K)-1.0e-3])
	Bxi_interp_in_new = zeros(2,length(dxi))
	Bxi_interp_in_new[1,:] = dxi
	Bxi_interp_in_new[2,:]  = B_xiK(Bxi_interp_in, dxi,K)

	return B_xiK(Bxi_interp_in_new,xi,K)
end

function B_xiK{T<:Real}(Bxi_interp_in::Array{T,2},xi::Array{T}, K::Int64)
	max_counter = 200
	if maximum(xi) > log(K)
		ArgumentError("Too large xi -- bigger than log(K)")
	end
	if minimum(xi) < 0 
		ArgumentError("Too small xi -- smaller than 0")
	end
	_xi = xi[:]
	B = Bxi_interp_in[2,lookup(Bxi_interp_in[1,2:size(Bxi_interp_in,2)-1],_xi)+1]
	B = B[:]
	dB = zeros(B)
	counter = 0
    F   =  9999999.0*ones(_xi)
    dF   =  9999999.0*ones(_xi)
	#println(minimum(_xi))
	qq = abs(_xi*1e-13 + 1e-13)
	qv = _xi*1.0e-13
	while !(all(abs(F) .< qq) || counter > max_counter)
		counter += 1
		#F[:] = xi_KB(K,B) -_xi
		#println(maximum(F))
		#dF[:] = dxi_KB(K,B)
		#dB[:] = -F./dF
		for i in 1:length(F)
			F[i] = xi_KB(K,B[i]) -_xi[i]
			dF[i] = dxi_KB(K,B[i])
			dB[i] = abs(F[i]) >= qv[i] ? -F[i]/dF[i] : 0.0
		end
		#dB[abs(F).<qv] = 0.0
		B[:] += dB
	end
	#println(size(_xi))
	if counter >= max_counter
		error("Newton-Raphson root polishing did not converge after $counter iterations. Problems are likely.")
	end

	repl_large = ((B.>100*K)&(_xi.!=log(K)))
	B[repl_large] = (K-1.0)./(2*(log(K) - _xi[repl_large]))
	B[_xi.==log(K)] = Inf

	repl_small = _xi .<  (pi*pi*1e-6/6.0)
	B[repl_small] = (1.0 + 1/(K-1.0))*_xi[repl_small]*(6/pi^2)

	B = reshape(B,size(xi))
	return B
end

function B_xiK{T<:Real}(Bxi_interp_in::Array{T,2},xi::Real, K::Int64)
	max_counter = 200
	if xi > log(K)
		ArgumentError("Too large xi -- bigger than log(K)")
	end
	if xi < 0 
		ArgumentError("Too small xi -- smaller than 0")
	end
	B = Bxi_interp_in[2,lookup(Bxi_interp_in[1,2:size(Bxi_interp_in,2)-1],xi)+1]
	dB = 0.0
	counter = 0
    F   =  9999999.0
    dF   =  9999999.0
	#println(minimum(_xi))
	qq = abs(xi*1e-13 + 1e-13)
	qv = xi*1.0e-13
	while !((abs(F) < qq) || counter > max_counter)
		counter += 1
		#F[:] = xi_KB(K,B) -_xi
		#println(maximum(F))
		#dF[:] = dxi_KB(K,B)
		#dB[:] = -F./dF
		F = xi_KB(K,B) -xi
		dF = dxi_KB(K,B)
		dB = abs(F) >= qv ? -F/dF : 0.0
		#dB[abs(F).<qv] = 0.0
		B += dB
	end
	#println(size(_xi))
	if counter >= max_counter
		error("Newton-Raphson root polishing did not converge after $counter iterations. Problems are likely.")
	end

	if (B > 100*K)&&(xi != log(K))
		B = (K-1.0)./(2*(log(K) - xi))
	elseif xi == log(K)
		B = Inf
	elseif xi <  (pi*pi*1e-6/6.0)
		B = (1.0 + 1/(K-1.0))*xi*(6/pi^2)
	end
	return B
end

function lookup{T<:Real}(table::Array{T,2},xi::Array{T})
	return lookup(table[:],xi)
end

function lookup{T<:Real}(table::Array{T,2},xi::Real)
	return lookup(table[:],xi)
end

function lookup{T<:Real}(table::Array{T,1},xi::Array{T})
	idx = zeros(Int64,size(xi))
	for i in 1:length(xi)
		idx[i] = last(searchsorted(table,xi[i]))
	end
	return idx
end

function lookup{T<:Real}(table::Array{T,1},xi::Real)
	idx = last(searchsorted(table,xi))
	return idx
end
