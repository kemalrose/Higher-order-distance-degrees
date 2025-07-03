-- code related to the illustrative example in the introduction of the paper
restart
m = 1 -- m is the dimension of X
n = 3 -- the ambient projective space is PP^n
KK = QQ
R = KK[t_0..t_m,u_0..u_n]

-- matf is the matrix whose components define the morphism f:X->PP^n
-- The one below, for m=1, defines the rational normal curve in PP^n:
matf = symmetricPower(n,matrix{{t_0..t_m}})

jacf = diff(matrix{{t_0..t_m}},transpose matf)
rkf = rank jacf -- rkf-1 coincides with the dimension of X
singf = trim minors(rkf, jacf)

k = 2
D = flatten apply(k+1, i-> first entries symmetricPower(i,matrix{{t_0..t_m}}))
Af = diff(transpose matrix{D},matf) -- this is the matrix A_p^{(k)}(f)
rank Af - 1 -- this is the dimension of the projective k-osculating space of f at some point p
kerAf = generators ker Af

-- Now we impose that the gradient of d_u^2(p) is in the column span of kerAf.
-- First, we need to choose a quadratic form Q:
Q = diagonalMatrix(apply(n+1, i-> 1)) -- this is the Euclidean quadratic form
Q = diagonalMatrix(apply(n+1, i-> binomial(n,i))) -- this is the Bombieri-Weyl quadratic form (when m=1)

au = matf-matrix{{u_0..u_n}}
grad = au*Q

-- DCfQk is the ideal of the kth-order distance correspondence of (f,Q):
time DCfQk = saturate(minors(rank(kerAf)+1, kerAf|transpose(grad)), singf);
-- DLfQk is the ideal of the kth-order distance locus of (f,Q):
time DLfQk = eliminate(toList(t_0..t_m), DCfQk);
codim DLfQk
degree DLfQk

see = method()
see Ideal := (I) -> netList I_*

see DLfQk
