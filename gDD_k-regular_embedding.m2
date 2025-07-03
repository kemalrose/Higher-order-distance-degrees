-- we utilize the recursive formula of Eq. (6.3) to determine the Chern classes of the higher-order
-- jet bundles of a k-regular embedding f. This is used to compute the kth generic distance degree gDD(f).
restart
needsPackage "Schubert2";
chernJetBundles = (m,k) -> (
    R := QQ[{n}|{L}|apply(m, i-> c_(i+1)),Degrees=>{0,1}|toList(1..m)];
    X := abstractVariety(m,R);
    TX := abstractSheaf(X, Rank=>m, ChernClass=>1+sum(m, i-> c_(i+1)));
    X.TangentBundle = TX; -- TX is the tangent bundle of X
    CTX := cotangentBundle X; -- CTX is the cotangent bundle of X
    symCTX := symmetricPower(n, CTX); -- symCTX is the symmetric power of CTX
    LX := abstractSheaf(X, Rank=>1, ChernClass=>1+L);
    X.TautologicalLineBundle = LX; -- LX is the line bundle defining the embedding of X in projective space
    chernPolyJ_1 = sub(chern(symCTX**LX),n=>1)*chern(LX);
    for i in 2..k do chernPolyJ_i = sub(chern(symCTX**LX),n=>i)*chernPolyJ_(i-1);
    O1 := apply(m+1, i-> part(i,chernPolyJ_k));
    O2 := apply(m+1, i-> part(i,chernPolyJ_k)*L^(m-i));
    O3 := sum O2;
    return (O1,O2,O3)
    )

-- An example: the classes c_i(J_k(f)) and gDD_2(f) for m=3 and k=2
(LchernJet,Lpolar,gDD) = chernJetBundles(3,2);
LchernJet
Lpolar
gDD


-- Another example: the generic kth order distance degree of a k-regular Segre-Veronese embedding.
N = {2,1} -- N is the list of projective dimensions
D = {3,5} -- D is the list of degrees of the embedding (N and D have the same length)
S = QQ[h_1..h_(#N)]/ideal(apply(#N, i-> h_(i+1)^(N#i+1))) -- Chow ring of the Segre-Veronese variety

k = 2 -- the Segre-Veronese embedding is k-regular if k is not larger than min(D)
(LchernJet,Lpolar,gDD) = chernJetBundles(sum(N),2); 

chernProj = product(#N, i-> (1+h_(i+1))^(N#i+1)) -- total Chern class of the Segre-Veronese variety
chernClassesSegVer = apply(sum(N)+1, i-> part(i,chernProj)) -- list of Chern classes of the Segre-Veronese variety

A = S**ring(gDD)

subLpolar = apply(Lpolar, i-> sub(i, A))
-- the following is the list of kth-order polar degrees of the k-regular Segre-Veronese embedding
subPolarDegrees = apply(subLpolar, l-> contract(product(#N, i-> h_(i+1)^(N#i)), sub(l, apply(sum(N), i-> c_(i+1)=>sub(chernClassesSegVer#(i+1),A))|{L=>sum(#N, i-> (D#i)*h_(i+1))})))
sum oo -- this is the generic kth order distance degree of the given k-regular Segre-Veronese embedding
