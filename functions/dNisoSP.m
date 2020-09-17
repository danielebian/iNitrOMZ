function iso = dNiso(varargin)

A.N = nan;
A.i14N = nan;
A.i15N = nan;
A.SP = nan;
A = parse_pv_pairs(A, varargin);

% Reference 15N/14N
Nstd1 =  0.3663./99.6337;
Nstd2 =  0.3663./100.0;

if ~isnan(A.SP) .* ~isnan(A.i15N(1)) .* ~isnan(A.N(1))
	iso.i14N=A.N-A.i15N; iso.N=A.N; iso.i15N=A.i15N;
	iso.i15N_A = 0.5.*(A.SP./1000.0.*Nstd1.*iso.i14N + iso.i15N);
	iso.i15N_B = iso.i15N - iso.i15N_A;
elseif  ~isnan(A.SP) .* ~isnan(A.i15N(1)) .* ~isnan(A.i14N(1))
        iso.N=A.i14N + A.i15N; iso.i14N=A.i14N; iso.i15N=A.i15N;
	iso.i15N_A = 0.5.*(A.SP./1000.0.*Nstd1.*iso.i14N + iso.i15N);
        iso.i15N_B = iso.i15N - iso.i15N_A;
end	

