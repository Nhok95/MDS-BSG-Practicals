## HWE

#AA = 756; SS
#AB = 200; SB
#BB = 144; BB

AA = 0; AA
AB = 107; AB
BB = 0; BB
total = AA + AB + BB; total

p = ((2*AA) + AB)/ (2* total); p # A
q = ((2*BB) + AB)/ (2* total); q # B

p+q # should be 1

AA.exp = p*p*total; AA.exp
AB.exp = 2*p*q*total; AB.exp
BB.exp = q*q*total; BB.exp

chi = ((AA-AA.exp)^2/ AA.exp) + ((AB-AB.exp)^2/ AB.exp) + ((BB-BB.exp)^2/ BB.exp); chi

df = 3-2 

# p.value <<<< 0.01 -> 6.63 chi
# p.value < 0.05 -> reject H0 -> alleles are out of equilibrium