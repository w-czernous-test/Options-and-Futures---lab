##### TODO: przenieść przykład "control variate" z przykładu 1-wymiarowego na średnią geom. 7 cen - uzupełnić funkcję inner_control_variate_1()  ###############

TEST = FALSE
REDUCE_VARIANCE = FALSE
# TEST = TRUE

set.seed(2023) # ustalmy ziarno generatora, dla powtarzalnych eksperymentów
######################### Stałe dla rynku i opcji: #############################

m = 10  # Opcję można zrealizować w chwilach: t_0=0, t_1=0.1,...,t_{10}=1.0.
T = 1.0
dt = T / m

r = .03
delta = .05
sigma = .40
K = 100

d = 7   # Wymiar. Stan rynku jest wektorem cen: $X = S = (S_1,...,S_d)$.
S_t0 = 100 # Ta sama cena startowa dla S_1,...,S_d, w chwili t=t_0.
X_0 = rep(S_t0, d) # Wartość stanu rynku w t=t_0.

# Stałe występujące w generacji procesu X = (S_1,...,S_d) i w funkcji gęstości prawdopodobieństwa przejścia dla tego procesu:
sigsqrt = sigma * sqrt(dt)
rdelta = (r - delta - sigma^2/2) * dt

######################### Funkcja wypłaty: #############################

mean_geom = function(x) { return (exp(mean(log(x))))  } # Średnia geometryczna z wektora o dodatnich współrzędnych.
h_i = function(i,x) { return ( exp(-r*i*dt) * max(0, mean_geom(x) - K) ) } # Opcja typu "call" dla średniej geometrycznej z S_1,...,S_d.

if (TEST) print(h_i(0, X_0 * 1.2345)) # powinno być około 23.45

######################### Funkcje dla rynku i opcji: #############################

# Ponieważ w tym przykładzie numerycznym, parametry (r,delta,sigma) 
# geometrycznego ruchu Browna są takie same dla każdej współrzędnej X=(S_1,...,S_d),
# generowanie wartości wektora X dla chwili t_{i+1}, gdy dane jest X(t=t_i) = x, jest proste:

generateXnext = function(x) {
  xnext = c()  
  for (n in 1:d)
    xnext[n] = x[n] * exp(rdelta + sigsqrt * rnorm(1))
  return (xnext)
}

if (TEST) print(generateXnext(X_0))

f_i_plus_1 = function(x,y) {
  f = 1.0
  for (n in 1:d) 
    f = f * dnorm( (log(y[n]/x[n]) - rdelta)/sigsqrt ) / (y[n]*sigsqrt)
  return (f)
}

if (TEST) {
  X1 = generateXnext(X_0)
  print(X1)
  print(f_i_plus_1(X_0, X1))
}
######################### Funkcje dla redukcji wariancji (control variates): #############################
variance_reduction_with_control_variate = function(x, y, EX) { 
  # redukcja wariancji y przy pomocy "control variate" x i jej wart.oczekiwanej EX;
  v = var(x)
  beta = if (v == 0) 0 else cov(x,y) / v
  return ( mean(y - beta*(x - EX)) )
}

inner_control_variate_1_value = function(x) {return (h_i(1,x))}
inner_control_variate_1 = function(x) {
####################### Your code goes here: #########################
  
  
} # E( h_1( S_{t+1} ) | S_t = x )

######################### Generacja siatki: #############################

# Zauważmy, że siatka X ma trzy wymiary: czas (i -> t_i), numer próby (j -> X_{i,j}) oraz współrzędną wektora stanu rynku (n -> x_n, gdzie X_{i,j}=(x_1,...,x_n)).

generate_mesh_X_1 = function() {        # generujemy X[1, , ] 
  for (j in 1:b) {
    X[1, j, ] <<- generateXnext(X_0) # generuje od razu cały wektor (x_1,...,x_d)
  }
}

generate_mesh_X_i_plus_1 = function(i) { # generujemy X[i+1, , ] dla i=1,...,m-1.
  for (j in 1:b) {
    # Losujemy węzeł $X_{i\ell}$ spośród węzłów $X_{i1}$, $\ldots$, $X_{ib}$, z jednakowym prawdopodobieństwem:
    ell = sample.int(b, size=1)
    # Następnie generujemy próbkę z rozkładu o gęstości $f_{i+1}(X_{i\ell},\cdot)$:
    X[i+1, j, ] <<- generateXnext(X[i,ell, ]) # generuje od razu cały wektor (x_1,...,x_d)
  }
}

generate_mesh_X = function() {
  X <<- array(dim = c(m,b,d))   # pusta tablica X[i,j,n], dla i=1,...,m,  j=1,...,b,  n=1,...,d.
  generate_mesh_X_1()
  for (i in 1:(m-1)) {
    generate_mesh_X_i_plus_1(i)
  }
}

if(TEST) {
  b = 5
  generate_mesh_X()
  print(X)
}

############################## Funkcja wagowa: ##############################

# Przyjmijmy, że gęstość prawdopodobieństwa przejścia $f_i$ dla procesu Markowa stanu rynku,
# jest niezależna od chwili czasu i, tj. że $f_i(x,y) = f(x,y)$, dla $i=1,...,m$.

find_Wi_k_denominator = function() { 
  # Obliczenia wstępne - to obniża złożoność find_hatV_all() z O(mb^3) do O(mb^2).
  
  Wi_k_denominator <<- array(dim = c(m,b)) 
  
  for (i in 1:(m-1))
    for (k in 1:b) {
      sum = 0
      for (ell in 1:b)
        sum = sum + f_i_plus_1(X[i,ell,], X[i+1,k,])
      Wi_k_denominator[i,k] <<- sum / b
    }
}

Wi_k = function(i,k,x) {   # Dla i=1,...,m-1.
  numerator   = f_i_plus_1(x, X[i+1,k,])
  denominator = Wi_k_denominator[i,k]
  return ( numerator / denominator )
}

############################## Estymator górny: ##############################

find_hatV_mj = function(j) {    # Dla i=m.
  X_mj = X[m, j, ] # wektor (d-wymiarowy) stanu rynku
  hatV[m,j] <<- h_i(m, X_mj)
}

find_hatV_ij = function(i,j) { # Dla i=1,...,m-1.
  C = c()
  cv = c()
  X_ij = X[i, j, ] # wektor (d-wymiarowy) stanu rynku
  for (k in 1:b) {
    Wi_jk = Wi_k(i,k,X_ij)
    C[k] = Wi_jk * hatV[i+1,k]
    cv[k] = Wi_jk * inner_control_variate_1_value(X[i+1,k,])
  }
  Ecv = inner_control_variate_1(X_ij)
  if (REDUCE_VARIANCE)
    hatC_ij = variance_reduction_with_control_variate(cv, C, Ecv)
  else 
    hatC_ij = mean(C) # bez redukcji wariancji

  hatV[i,j] <<- max(h_i(i, X_ij), hatC_ij)
}

find_hatV_all = function() {
  
  hatV <<- array(dim = c(m,b)) # pusta tablica hatV[i,j]

  find_Wi_k_denominator()
  
  for (j in 1:b)
    find_hatV_mj(j)
  for (i in (m-1):1) {
    for (j in 1:b)
      find_hatV_ij(i,j)
  }
}

# Gdy mamy już gotowe $\hat V_{1,j}$, dla j=1,...,b, wówczas dopiero możemy wyliczać $\hat V_0$:

hatV_0 = function() { # Dla i=0.
  cv = c()
  for (k in 1:b) {
    cv[k] = inner_control_variate_1_value(X[1,k,])
  }
  Ecv = inner_control_variate_1(X_0)
  C = hatV[1,]
  if (REDUCE_VARIANCE)
    hatC_ij = variance_reduction_with_control_variate(cv, C, Ecv)
  else
    hatC_ij = mean(C) # bez redukcji wariancji
  
  return ( max(h_i(0, X_0), hatC_ij) ) # Dopuszczamy możliwość realizacji opcji również w chwili $t_0$.
}

if (TEST) {
  find_hatV_all()
  print(hatV_0())
}

############################## Estymator dolny: ##############################

hatC = function(i, x) { # Dla i = 0,1,...,m. Parametr x to wektor (d-wymiarowy) stanu rynku.
  if (i == m)
    return (0)

  if (i == 0)
    return (mean(hatV[1, ]))

  sum = 0
  for (k in 1:b)
    sum = sum + Wi_k(i,k,x) * hatV[i+1,k]
  return (sum / b)
}

# Generujemy losową wartość stanu rynku $X_{i+1}$,
# dla danej wartości $X_i$ oraz chwili $t=t_i$,
# i robimy tak dla każdego i=0,1,...,m-1:

generate_X_path = function() {

  X_path <<- array(dim = c(m, d)) # pusta tablica X_path[i,n], dla i=1,...,m,  n=1,...,d.
  
  # Generujemy X_path[ i+1,  ] dla i=0 :
  X_path[1, ] <<- generateXnext(X_0)

  # Generujemy X_path[ i+1,  ] dla i=1,...,m-1 :
  for (i in 1:(m-1))
    X_path[i+1, ] <<- generateXnext(X_path[i, ])
}

if (TEST) {
  generate_X_path()
  print('X_path:')
  print(X_path)
}

# Dopuszczamy możliwość realizacji opcji również w chwili $t_0$.

tau_for_X_path = function() {
  # Sprawdzamy, czy przypadkiem nie warto zrealizować opcję już w chwili 0:
  if (h_i(0, X_0)  >= hatC(0, X_0)) {
    return (0)
  }
  # Jeśli nie, szukamy pierwszej takiej chwili:
  for (i in 1:m) {
    if (h_i(i, X_path[i, ]) >= hatC(i, X_path[i, ])) {
      return (i) # powinno zostać wywołane najpóźniej dla i=m, gdyż h_i >= 0.
    }
  }
}

generate_hatv = function() {
  generate_X_path()
  tau = tau_for_X_path()
  X_tau = X_path[tau, ]
  return (h_i(tau, X_tau))
}

if (TEST) {
  print("hat v_0:")
  print(system.time(v <<- generate_hatv()))
  print(v)
}

######################## Skrypt główny: #############################
b = 50
n_paths = 500
n = 5
higher_estimator = c()  # $\hat V_0$
lower_estimator = c()   # $\hat v_0$
for (j in 1:n) {
  message("siatka nr ",j,"/",n,"...",format(Sys.time(), "%a %b %e %H:%M:%S"))
  generate_mesh_X()
  find_hatV_all()
  higher_estimator[j] = hatV_0()
  lower_estimator[j] = mean(replicate(n_paths, generate_hatv()))
}
message("mean hat V_0 = ",mean(higher_estimator), "; variance of hat V_0 = ",sd(higher_estimator))
message("mean hat v_0 = ",mean(lower_estimator), "; variance of hat v_0 = ",sd(lower_estimator))
