# Vérification ligne à ligne de la version LaTeX

Comparaison de `Methodes_numeriques_LaTeX.pdf` avec le scan original
(`Methodes numeriques_compressed.pdf`, 122 pages). Les numéros `p-NNN` sont
les pages du fichier PDF scanné, en commençant à 000 (p-008 = 9ᵉ page du PDF).

Légende :

- **OK** : la transcription est fidèle.
- **CORR** : erreur de l'original corrigée dans la version LaTeX (avant → après).
- **FIDÉLITÉ** : écart que j'avais introduit sans raison lors de la première
  transcription, remis conforme à l'original lors de cette vérification.
- **DOUTE NON CORRIGÉ** : passage douteux laissé tel quel.
- **ILLISIBLE** : passage coupé ou illisible sur le scan.

Les fautes d'orthographe ou de frappe évidentes (« à calcule », « Weistrass »,
« Cranck-Nicholson », numérotation « III-3-1-1 » au lieu de « II-3-1-1 »…)
ont aussi été corrigées ; seules certaines sont listées.


- p-000 titre : OK (annotation manuscrite « CADO » non reprise)
- p-001/002 préface : OK ; « Il faut pas aussi oublier » -> « Il ne faut pas aussi oublier » (langue)
- p-003..006 table des matières originale : structure identique à la nôtre
- p-008/009 (I-1, I-2-1-1) : OK ; « à calcule » -> « à calculer »
- p-010 : CORR a'_ij = a_ij - (a_i1/a_11) a_ij  -> a_1j (dernier indice)
- p-010 : algorithme imprimé avec étape i ; corrections manuscrites (k, l) adoptées -> notre version k, l
- p-011 : « faisant varier i de 1 à n-1 » -> k (cohérent avec algorithme) ; CORR 3e ligne du système triangulaire « = b_2 » -> b_3
- p-012 : exercices, Jordan : OK
- p-013 : FIDÉLITÉ « b- Décomposition de la matrice » : le « b- » n'existe pas dans l'original -> retiré
- p-014 : relations LS=A et formules l_ij, s_ij : OK
- p-015 : Cholesky : CORR m_pj « j allant de 2 à n » -> « j allant de p+1 à n »
- p-016 fragmentation : OK
- p-017 : CORR systèmes complexes « NX - MY = d » -> « NX + MY = d » ((M+iN)(X+iY) = MX-NY + i(NX+MY))
- p-018 : OK
- p-019 : Gauss-Seidel, critères d'arrêt : OK ; CORR notation du dernier test « d_i = Σ_i ... » -> « d = ... » (somme sur i)
- p-020 : doublon de la page p-019 dans le scan
- p-021..024 (II-1, II-2-1, II-2-2) : OK
- p-025 : CORR « O(x^3) » -> O(h^3) ; CORR « h = x + h - x0 = -f/f' » (formule incohérente) -> « h = x̄ - x0 ≈ -f/f' »
- p-026 : CORR Newton complexe : y_{n+1} ≈ y_n - (GH'-HG')/(G'^2+H'^2) -> (HG'-GH') (signe) ; évaluation notée |_{x_n} -> |_{z_n} ; FIDÉLITÉ « ≈ » rétabli (j'avais mis « = »)
- p-026 : CORR dichotomie « intervalle contenant la section » -> « la solution » ; « (a-b)/2^n » -> (b-a)/2^n ; numérotation « III-3-1-1 » -> II-3-1-1
- p-027 : Horner : OK
- p-028 : CORR « P(x) = Q(x) + R » -> « P(x) = (x-r) Q(x) + R » ; numérotation « III-3-1-2 » -> II-3-1-2 ; « peut être également être utilisé » (langue)
- p-029 : CORR « Si b_n, alors » -> « Si b_n = 0, alors » ; numérotation « III-3-1-3 » -> II-3-1-3
- p-030 : CORR b_k = a_k + S b_{k-1} - P b_{k-2} « (k = 1 à n) » -> (k = 2 à n) (b_1 donné à part)
- p-031 : CORR méthode de Lin complexe : x_1 = (a_n G + b_n H)/(G^2+H^2) -> signe « - » manquant ; CORR Horner dérivée « c_n = b_n + x c_{n-1} (k=1,n-1) » -> c_k = b_k + x c_{k-1}
- p-032 : II-3-3 fin, II-3-4 : OK
- p-033 : Bairstow : CORR b_k « (k = 1,...,n) » -> (k = 2,...,n) ; CORR G(S,P) = « b_n - S b_{n-1} » -> b_n : les formules finales de α, β, Δ du livre ne sont cohérentes qu'avec G = b_n (même zéro, b_{n-1} = 0)
- p-034 : CORR relation garbled « ∂b_k/∂S - b_{k-1} = c_{k-1} » (contredit la 1re) -> « ∂b_k/∂S - b_{k-1} = S c_{k-2} - P c_{k-3} » ; FIDÉLITÉ : j'avais remplacé ces lignes par une dérivation rédigée -> remise au plus près de l'original
- p-034 : CORR suite c_k : « c_1 = a_1 + S c_0 ; c_k = a_k + S c_{k-1} - P c_{k-2} (k=1..n) » -> b_1, b_k et (k = 2..n)
- p-034 : CORR β = « b_n c_{n-3} - b_{n-1} c_{n-1} » -> b_n c_{n-2} - b_{n-1} c_{n-1} ; CORR Δ = « c_{n-2} - c_{n-1} c_{n-3} » -> c_{n-2}^2 - ... (vérifié par la résolution du système de Newton)
- p-035 : OK
- p-036/037/038 : Viète degré 3, 4, 5 : OK (relations vérifiées)
- p-039/040 : Viète général, Descartes, Gua, lacunes, Sturm : OK (exemples vérifiés)
- p-041 : CORR exemple de Sturm : (x-1)(x^3+7x^2+x+7) = x^4+6x^3-6x^2+6x « +7 » -> « -7 » (sinon 2 changements de signe et non 3)
- p-042 : Newton 2 équations : OK (signes vérifiés) ; CORR système n équations : somme sur i implicite rendue explicite (Σ_i) ; FIDÉLITÉ « , … » ajouté après f_{1,x} retiré
- p-043/044 : OK ; « Weistrass » -> Weierstrass ; CORR théorème 1 « fonction contenue » -> « continue » ; annotations manuscrites non reprises
- p-045 : CORR « f_k(x) = y_1 f_1 + ... + y_n f_n » -> « f(x) = y_1 f_1 + ... + y_{n+1} f_{n+1} » (cohérent avec la somme jusqu'à n+1)
- p-046 : CORR « (x_1,y_1) et (y_2,y_2) » -> (x_2,y_2) ; CORR « y_3 f_2(x) » -> y_3 f_3(x)
- p-048 : algorithme de Lagrange : OK
- p-049 : différences divisées, Newton : OK
- p-050 : CORR différence régressive « ∇^α y_k = ∇^{α-1} y_{k+1} - ∇^{α-1} y_k » (= définition progressive) -> ∇^{α-1} y_k - ∇^{α-1} y_{k-1} ; FIDÉLITÉ parenthèse ajoutée (∇y_k = ...) retirée
- p-050 : CORR « x_{k-1} = x_k + h » -> x_{k+1} ; FIDÉLITÉ dernier terme de Newton-Gregory remis comme l'original : Δ^{k+1} y_1 / (h^{k+1}(k+1)!) (j'avais réécrit en Δ^k/(h^k k!))
- p-051 : erreur de troncature, Spline (a)(b) : OK
- p-052 : CORR Spline (c) « i = 1,...,n+1 » -> n-1 (n points => n-1 intervalles) ; CORR S''(x) dénominateur « x_{i+1} - x_1 » -> x_i ; CORR S'(x) « S''(x_{i+1}) - S(x_i) » -> S''(x_i) ; « (d) » ajouté devant la 4e propriété mécanique (le texte parle ensuite de « quatrième propriété »)
- p-053 : CORR relation (b) « S'_{i-1} h_{i-1} + (h_i+h_{i-1}) S''_i » -> S''_{i-1} et facteur 2 (relation standard des splines cubiques) ; FIDÉLITÉ « (i = 2,...,n-1) » ajouté à (b) retiré ; libellés (a)/(b) ajoutés sur les équations car le texte y renvoie
- p-054 : NOTATION différence centrée notée « Δ » dans l'original (même symbole que la différence progressive) -> δ
- p-055 : Hermite trigonométrique, différences doubles : OK
- p-056 : CORR différence double générale « Δ^{α0} z_{1β} + β Δ^{α0} z_{1,β-1} + ... + Δ^{0β} z_11 » -> formule binomiale Δ^{α0} z_{1,β+1} - β Δ^{α0} z_{1β} + ... + (-1)^β Δ^{α0} z_11 (signes et indices)
- p-056 : CORR interpolation double : « (x-x_1)/h + Δ^{10} » -> produit ; terme d'ordre m : « (x-x_1)(x-x_2)+...+(x-x_{m-1}) » -> produit jusqu'à x_m, et indices décalés d'un cran (x_{m-1}, x_{m-2}, y_m) pour les autres termes
- p-057 : CORR « puisque n = m » -> n << m ; « f(x_k) ; f_k » -> f(x_k) ≈ f_k
- p-058/059 : moindres carrés, a_0, a_1, système normal, exercices 1-3 : OK
- p-060 : réponse de l'exercice 3 : OK (vérifiée numériquement)
- p-061 : IV intro, définitions : OK
- p-062 : OK ; « h = 10^-3.10^-1 » -> « 10^-3, 10^-1 »
- p-063 : CORR signe de l'erreur « E = 1/2 h f''(ξ) » -> -1/2 h f''(ξ) (découle de la formule de E donnée juste avant, avec x = x_1)
- p-064 : FIDÉLITÉ hypothèses remises comme l'original : C^4 et « a < x + 2h < b » (j'avais mis C^3 et a < x < x+2h < b) ; CORR « ξ ∈ [x-h ; x+h] » -> [x ; x+2h] (la formule utilise x, x+h, x+2h)
- p-064 : FIDÉLITÉ C^6 rétabli (j'avais mis C^5) ; ILLISIBLE terme d'erreur coupé par le scan « + (1/30) h^4 d/dx[… » -> forme standard (1/30) h^4 f^(5)(ξ) conservée, à vérifier sur l'original papier
- p-064 : CORR « fonction du temps dans le tableau » (verbe manquant) -> « est donnée dans le tableau »
- p-065 : CORR interpolation linéaire « f(a)(x-a)/(b-a) + f(b)(x-b)/(a-b) » (termes inversés) -> f(a)(x-b)/(a-b) + f(b)(x-a)/(b-a) ; CORR E = ∫ f''(ξ)(x-a)(x-b)dx -> facteur 1/2 manquant (sinon on n'obtient pas -1/12)
- p-066 : OK
- p-067 : CORR « E = h^2/12 (b-a) f''(η) » -> signe « - » (découle de la ligne précédente)
- p-068 : algorithme trapèzes : OK ; CORR erreur de Simpson « f''(ξ) » -> f'''(ξ) (erreur d'interpolation de degré 2 = f'''/3!) ; formule composite coupée au bord du scan, terme final f(x_n) complété
- p-069 : CORR exemple d'équation intégrale « f(x) = 5x + ∫(t+x)f(t)dt = 0 » -> « = 0 » parasite retiré ; DOUTE NON CORRIGÉ exercice 6 : c = (a+b)/4 et d = 3(a+b)/4 ne sont dans [a,b] que si a = 0 (probablement a + (b-a)/4 …) -> laissé tel quel
- p-070 : OK
- p-071 : V intro : OK
- p-072/073/074 : OK
- p-075 : CORR RK2 « f[x_k + h, y_k + f(x_k,y_k)] » -> y_k + h f(x_k,y_k) (facteur h manquant, cf. Euler) ; RK4 de Runge : approximations y_1, y_2 : OK
- p-076 : Runge ordre 4 : OK (coquille « x_0 + h/2, + y_0 » corrigée)
- p-077 : RK4 : OK ; « Kutta-Merdsen » -> Kutta-Merson (nom correct)
- p-078 : Adams-Bashforth : OK
- p-079 : OK ; « second nombre » -> « second membre »
- p-080 : prédicteur-correcteur : OK
- p-081 : méthode de tir : OK ; CORR « méthode de Newton (voir chapitre suivant) » -> « voir chapitre II » (Newton est au chapitre II)
- p-082/083 : superposition, 1re variante : OK
- p-084/085 : 2e variante (a)…(h) : OK (relations (e)(f)(g)(h) revérifiées)
- p-086 : CORR différences finies « (2 - h p_k - h^2 q_k) » -> « + h^2 q_k » (dans l'équation générale et dans tout le système ; vérifié par substitution)
- p-087 : OK
- p-088 : CORR « nœuds de subdivision de l'intégrale [a,b] » -> « de l'intervalle » ; CORR « t_i = a + ih » -> a + (i-1)h (cohérent avec t_1 = a)
- p-089 : CORR exemple « 5x/9 » -> 5x/6 (dans l'énoncé et f(x)) : seule valeur compatible avec la solution y(x) = x annoncée et avec le système qui suit
- p-090 : CORR système : « -1/9 » manquant dans les équations de y_2 et y_3 (sinon y = (0, 1/2, 1) n'est pas solution)
- p-091 : VI intro, classification : OK ; « Ce présente » -> « Ce chapitre présente »
- p-092/093/094 : OK
- p-095 : CORR « n et m les nombres réels » -> entiers ; CORR « y_i = c + jk » -> y_j ; CORR (P2) « [2(h/k)^2 + 1] U_ij » -> 2[(h/k)^2 + 1] U_ij (vérifié en multipliant (P) par -h^2)
- p-096 : (P3) « g(x_0, y_i) » -> y_j ; OK
- p-097 : (P4), théorème du maximum : OK (revérifiés)
- p-098 : OK
- p-099 : CORR « m le nombre réel défini par mh = l » -> entier (j'avais simplement supprimé « réel » : remplacé par « entier », signalé) ; CORR (C2) « j = 1, 2, …, ∞ » -> j = 0, 1, 2, … (le texte utilise (C2) avec j = 0 pour obtenir U_{i,1})
- p-100 : stabilité du schéma explicite : OK (« Le schéma est si » -> « est stable si ») ; FIDÉLITÉ « (en module) » ajouté retiré
- p-101 : CORR titre « régressives ou schéma explicite » -> implicite (le texte dit ensuite « ce schéma implicite ») ; CORR (C5) « (1+λ) U_ij » -> (1+2λ) ; « 1 ≤ i ≤ n-1 » -> m-1 (notation du paragraphe) ; facteur d'amplification : OK
- p-102 : CORR Crank-Nicolson « a^2/(2h) » -> a^2/(2h^2) ; « Cranck-Nicholson » -> Crank-Nicolson ; Richardson, Dufort-Frankel : OK
- p-103 : exercices, (O1), (O2) : OK (O2 revérifiée)
- p-104 : (O3)…(O6) : OK (O6 revérifiée)
- p-105 : OK
- p-106 : CORR « (Iak/jh) sin ph » -> (Iak/h) sin ph (« j » parasite) ; « A_{i,j+1} » -> A_{j+1} (cohérent avec V = A_j e^{Ipih})
- p-107 : CORR « tous les V_{i,j-1} » -> V_{i,j+1} ; matrice (1, Ic ; Ic, 1-c^2) et condition ak/h < 2 : OK
- p-108 : schéma implicite : OK ; CORR Lax, 2e équation « V_{i+1,j-1} - V_{i-1,j-1} » -> instant j (le texte dit : schéma explicite naturel avec moyennes)
- p-109 : CORR Lax-Wendroff : « W_{i-1,j+1} - W_{i-1,j} » -> W_{i+1,j} - W_{i-1,j} ; « a^2k^2/(2h) » -> a^2k^2/(2h^2) (2 fois) ; « V_{i+1,j+1} - V_{i-1,j+1} » -> instant j (sinon pas Lax-Wendroff explicite, et la matrice donnée ne correspond pas)
- p-110 : exercice 1 : OK
- p-111..113 : exercices 1-10 : OK ; CORR exo 7 « y''' = f(x_1, y_1, y', y'') » -> f(x, y, y', y'') (indices parasites, y_1 désigne déjà la condition y(a))
- p-114 : CORR exo 11 « d²Y/dt² + ω0 sin Y = 0 » -> ω0² (homogénéité, et T0 = 2π/ω0 donné ensuite) ; « A = sin² Y0/2 » écrit sin²(Y0/2) (notation, valeur standard de la période du pendule)
- p-115..116 : exercices 12-14 : OK
- exo 1 (p-110) : CORR loi de Planck « λ^5 exp(hc/λkT) - 1 » -> λ^5 [exp(hc/λkT) - 1] (parenthèse)
- exo 3 : CORR (vérifié sur sources : équation de De Santis) dénominateur « (1 - y^3) » -> (1 - y)^3 ; « facteur de compressibilité b » -> z (b est le paramètre cherché)
- exo 4 : « Adam - Bas Forth » -> Adams-Bashforth
- p-117 : exo 15 : OK ; CORR exo 16 « (t,x) ∈ R+ ∪ ]0,l[ » -> × ; « u(l>t) » -> u(l,t) ; « t = 2x » -> t = 2k ; liste des points « i-2 ; i-1 ; i+1 ; i+2 » -> i ajouté (f_i intervient forcément)
- p-118 : exo 17 : OK (figure redessinée : rectangle intérieur centré ; sur le croquis original il est un peu décalé vers le haut, sans cote de position)
- p-119 : CORR exo 18 « ∂²u/∂u² » -> ∂²u/∂x² ; figure : OK
- p-120 : exo 19 : « La représentation d'une grandeur physique » -> « La répartition » (comme exo 18) ; figure : OK
- p-121 : page blanche
