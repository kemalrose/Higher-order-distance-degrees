-- code related to Example 5.10 in the paper
restart
m = 1 -- here we consider X = PP^m
n = 4 -- the ambient projective space of f(X) is PP^n
KK = QQ
R = KK[t_0..t_m,u_0..u_n]

-- 1) We consider the morphism f and compute DL_k(f,Q)
-- matf is the matrix whose components define the morphism f:X->PP^n
-- The one below, for m=1, defines the morphism of Example 5.10 in the paper:
matf = matrix{{t_1*(t_0^3+t_0*t_1^2+t_1^3),t_0*(t_0-t_1)^2*(t_0+t_1),-(t_0^2+t_1^2)*(t_0^2-t_0*t_1-t_1^2),-t_1*(t_0-t_1)*(t_0+t_1)^2,-(t_0+t_1)*(t_0^3-t_0^2*t_1+t_1^3)}}

jacf = diff(matrix{{t_0..t_m}},transpose matf)
rkf = rank jacf -- rkf-1 coincides with the dimension of X
singf = trim minors(rkf, jacf)
decompose singf -- this is the ideal (t_0,t_1)

k = 2
D = flatten apply(k+1, i-> first entries symmetricPower(i,matrix{{t_0..t_m}}))
Af = diff(transpose matrix{D},matf) -- this is the matrix A_p^{(k)}(f)
rank Af - 1 -- this is the dimension of the projective k-osculating space of f at some point p
kerAf = generators ker Af

-- Now we impose that the gradient of d_u^2(p) is in the column span of kerAf.
-- Here we consider the standard Euclidean inner product
grad = matf-matrix{{u_0..u_n}}

-- DCfQk is the ideal of the kth-order distance correspondence of (f,Q):
If = minors(rank(kerAf)+1, kerAf|transpose(grad));
time DCfQk = saturate(If,ideal(t_0));
time DCfQk = saturate(DCfQk,ideal(t_1));
-- DLfQk is the ideal of the kth-order distance locus of (f,Q):
time DLfQk = eliminate(toList(t_0..t_m), DCfQk);


-- 2) We consider the morphism f', the composition of f with the projection on all but the last coordinate,
-- and we compute DL_k(f',Q'), where Q' is the quadratic form induced by Q via the projection.
matf' = matf_{0..n-1}

jacf' = diff(matrix{{t_0..t_m}},transpose matf')
rkf' = rank jacf'
singf' = trim minors(rkf', jacf')
decompose singf' -- this is the ideal (t_0,t_1)

Af' = diff(transpose matrix{D},matf')
rank Af' - 1
kerAf' = generators ker Af'

-- Here the quadratic form Q' is again the standard Euclidean inner product, in this case in RR^4
grad' = matf'-matrix{{u_0..u_(n-1)}}

-- DCfQk' is the ideal of the kth-order distance correspondence of (f',Q'):
If' = minors(rank(kerAf')+1, kerAf'|transpose(grad'));
time DCfQk' = saturate(If',ideal(t_0));
time DCfQk' = saturate(DCfQk',ideal(t_1));
-- DLfQk' is the ideal of the kth-order distance locus of (f',Q'):
time DLfQk' = eliminate(toList(t_0..t_m), DCfQk');

-- we compare the two distance loci DL_k(f,Q) and DL_k(f',Q')
(codim DLfQk,degree DLfQk)
(codim DLfQk',degree DLfQk')

DLfQkW = DLfQk + ideal(u_n)
dec = decompose DLfQkW;
#dec

see = method()
see Ideal := (I) -> netList I_*

-- the following two surfaces are isomorphic in PP^3 and
-- their affine patches in {u_3 = 1} are plotted in Figure 5 of the paper.
see ideal(dec#0_1)
see DLfQk'
