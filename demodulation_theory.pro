;---------------------------------------
; Calculate the derivative
;  dD[alpha,beta]/dO[i,j]
;---------------------------------------
function derDO, invA, O, D, alpha, beta, i, j
	
	t1 = 0.d0
	if (i eq beta) then begin
		t1 = invA[alpha,j]
	endif
			
	t2 = 0.d0
	for n = 0, 3 do begin				
		t2 = t2 + invA[alpha,j]*D[n,beta]*O[i,n]
	endfor
				
	t3 = D[alpha,i]*D[j,beta]
			
	derD = t1 - t2 - t3
	return, derD
end

;---------------------------------------
; Calculate the derivative
;  deps[alpha]/dO[i,j]
;---------------------------------------
function derefficiencyO, invA, O, D, alpha, i, j
	nn = n_elements(O[*,0])
	
	den = 0.d0
	for p = 0, nn-1 do begin
		den = den + D[alpha,p]^2
	endfor
	
	denominator = den^(3.d0/2.d0)
	t = 0.d0	
	for n = 0, nn-1 do begin			
		t = t - D[alpha,n] * derDO(invA, O, D, alpha, n, i, j) / denominator			
	endfor
		
	return, t / sqrt(nn)
end

;---------------------------------------
; Calculate the influence on measurement errors on the modulation matrix
; on the demodulation matrix
;---------------------------------------
pro demodulation_theory, noise, polarimeter
	
	if (polarimeter eq 'ASP') then begin
		n = 8	
		O = transpose([[1.d0, 0.77d0, 0.41d0, -0.36d0],$
			  [1.d0, -0.06d0, 0.41d0, -0.86d0],$
		  	[1.d0, -0.06d0, -0.41d0, -0.86d0],$
		  	[1.d0, 0.77d0, -0.41d0, -0.36d0],$
		  	[1.d0, 0.77d0, 0.41d0, 0.36d0],$
		  	[1.d0, -0.06d0, 0.41d0, 0.86d0],$
		  	[1.d0, -0.06d0, -0.41d0, 0.86d0],$
		  	[1.d0, 0.77d0, -0.41d0, 0.36d0]])
	endif
	
	if (polarimeter eq 'TIP') then begin
		n = 4
		O = transpose([[1.d0, 0.47d0, -0.68d0, 0.48d0],$
		  [1.d0, -0.91d0, -0.19d0, -0.13d0],$
		  [1.d0, -0.11d0, 0.57d0, 0.72d0],$
		  [1.d0, 0.68d0, 0.27d0, -0.58d0]])
	endif
	
	if (polarimeter eq 'IDE') then begin
		n = 6
		O = transpose([[1.d0, 1.d0, 0.d0, 0.d0],$
		  [1.d0, -1.d0, 0.d0, 0.d0],$
		  [1.d0, 0.d0, 1.d0, 0.d0],$
		  [1.d0, 0.d0, -1.d0, 0.d0],$
		  [1.d0, 0.d0, 0.d0, 1.d0],$
		  [1.d0, 0.d0, 0.d0, -1.d0]])
	endif		
	
; According to del Toro Iniesta & Collados (2000), D is the pseudo-inverse that
; maximizes the efficiency of the demodulation process	
	A = O ## transpose(O)
	invA = invert(A)
	D_correct = transpose(invA ## O)
	
	sigma = replicate(noise,n,4)
	In = [1.d0,1.d-3,1.d-3,1.d-3]
	Iout = transpose(O) ## In
	
	cov_theory = dblarr(4,n,4,n)
	print, 'Calculating theoretical covariance for D...'
	for alpha = 0, 3 do begin
		for beta = 0, n-1 do begin
			for a = 0, 3 do begin
				for b = 0, n-1 do begin					
					for i = 0, n-1 do begin
						for j = 0, 3 do begin
							derD1 = derDO(invA, O, D_correct, alpha, beta, i, j)
							derD2 = derDO(invA, O, D_correct, a, b, i, j)
							cov_theory[alpha,beta,a,b] = cov_theory[alpha,beta,a,b] + $
								derD1*derD2*sigma[i,j]^2
						endfor
					endfor
				endfor
			endfor
		endfor
	endfor
	
	variance = dblarr(4,n)
	for alpha = 0, 3 do begin
		for beta = 0, n-1 do begin
			variance[alpha,beta] = cov_theory[alpha,beta,alpha,beta]
		endfor
	endfor
	
	; View what happens with an incoming Stokes vector (1,1e-3,1e-3,1e-3)	
	covariance_iin_theory = dblarr(4,4)
	print, 'Calculating theoretical covariance for Stokes vector...'
	for i = 0, 3 do begin
		for j = 0, 3 do begin
			for alpha = 0, n-1 do begin
				for beta = 0, n-1 do begin					
					covariance_iin_theory[i,j] = covariance_iin_theory[i,j] + $
						iout[alpha]*iout[beta]*cov_theory[i,alpha,j,beta]
				endfor
			endfor
		endfor
	endfor
	
	eps = dblarr(4)
	for i = 0, 3 do begin
		eps[i] = 1.d0 / sqrt(n*total(D_correct[i,*]^2))
	endfor
	
	cov_eps = dblarr(4,4)
	for alpha = 0, 3 do begin
		for beta = 0, 3 do begin
			for i = 0, n-1 do begin
				for j = 0, 3 do begin
					der1 = derefficiencyO(invA, O, D_correct, alpha, i, j)
					der2 = derefficiencyO(invA, O, D_correct, beta, i, j)
					cov_eps[alpha,beta] = cov_eps[alpha,beta] + der1*der2*sigma[i,j]^2
				endfor
			endfor
		endfor
	endfor
	
	print, 'Standard deviation D (x 1e-3)'
	print, transpose(sqrt(variance)*1.d3)
		
	print, 'Standard deviation D (x 1e-3) [LaTeX format]'
	for j = 0, 3 do begin
		res = reform(sqrt(variance[j,*])*1.d3)
		str = string(res,FORMAT='(F4.2)')
		for i = 0, n-2 do begin
			str[i] = str[i] + ' & '
		endfor
		str[n-1] = str[n-1] + '\\'
		print, str
	endfor	
	
	print, 'Covariance Iin'
	print, transpose(covariance_iin_theory)
	
	print, 'Covariance Iin (x 1e-6) [LaTeX format]'
	for j = 0, 3 do begin
		res = reform((covariance_iin_theory[j,*])*1.d6)
		str = string(res,FORMAT='(F5.2)')
		for i = 0, 2 do begin
			str[i] = str[i] + ' & '
		endfor
		str[3] = str[3] + '\\'
		print, str
	endfor
	
	print, 'Efficiencies'
	print, eps
	
	print, 'Standard deviation diagonal Iin'
	print, sqrt(diagon(covariance_iin_theory))
	
	print, 'Standard deviation efficiencies (x 1e-3)...'
	print, sqrt(diagon(cov_eps))*1.d3
	
	hes = elmhes(covariance_iin_theory)
	evals = hqr(hes,/double)
	evec = double(eigenvec(covariance_iin_theory, evals, residual = residual, /double))
	
	print, 'Eigenvalue -> Eigenvector: '
	for i = 0, 3 do begin
		print, double(evals[i]), ' -> ', evec[*,i]
	endfor
	
	print, 'Eigenvectors [LaTeX format]'
	for j = 0, 3 do begin
		res = reform((evec[*,j]))
		str = string(res,FORMAT='(F6.3)')
		for i = 0, 2 do begin
			str[i] = str[i] + ' & '
		endfor
		str[3] = str[3] + '\\'
		print, str
	endfor
	
	sigma = invert([[1.35d-6,-1.61d-6],[-1.61d-6,4.54d-6]])
	sigma_rotated = invert([[6.803d-7,0.d0],[0.d0,5.208d-6]])
	mux = 1.d0
	muy = 1.d-3
	x = findgen(50)/49.*(1.005-0.995)+0.995
	y = findgen(50)/49.*(0.01+0.01)-0.01
	res = fltarr(50,50)
	res_rotated = fltarr(50,50)
	for i = 0, 49 do begin
		for j = 0, 49 do begin
			xx = [x[i]-mux,y[j]-muy]
			res[i,j] = exp(-0.5*xx##(sigma##transpose(xx)))
			res_rotated[i,j] = exp(-0.5*xx##(sigma_rotated##transpose(xx)))
		endfor
	endfor
	contour,res,x,y,levels=[1-0.95,1-0.68]
	contour,res_rotated,x,y,levels=[1-0.95,1-0.68],/overplot
	stop
end