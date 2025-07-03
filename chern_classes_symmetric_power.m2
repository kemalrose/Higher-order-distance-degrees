-- this code is used to determine the Chern classes of the symmetric power of a vector bundle
restart
needsPackage "Schubert2";
m = 3 -- dimension of X
X = abstractVariety(m, QQ[{n}|apply(m, i-> c_(i+1)),Degrees=>toList(0..m)]);
r = 3 -- we consider a vector bundle of rank r with Chern classes c_0=1,c_1,...,c_r
T = abstractSheaf(X, Rank=>r, ChernClass=>1+sum(r, i-> c_(i+1)));
X.TangentBundle = T
chsymmT = chern symmetricPower(n, T)
sub(chsymmT, n=>2) -- this is the total Chern class of Sym^n(T)
