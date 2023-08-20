def TwoSE(a,b,
          c,d,   n,m): # 2 sim. equations
  y = ((n*c)-(m*a)) / ((b*c)-(d*a))
  x = (n - (b*y)) / a
  return [x,y]

def ThreeSE(a,b,c,
            d,e,f,
            g,h,j,   n,m,p): # 3 sim. equations
  tse = TwoSE(
    b*d-a*e, c*d-a*f,
    e*g-d*h, f*g-d*j,   n*d-m*a, m*g-p*d
  ) # calculates y and z
  x = (n - (b*tse[0] + c*tse[1])) / a
  tse.insert(0,x)
  return tse

def FourSE(a,b,c,d,
           e,f,g,h,
           i,j,k,l,
           n,m,p,q,   s,t,u,v): # 4 sim. equations
  hse = ThreeSE(
    b*e-a*f, c*e-a*g, d*e-a*h,
    j*n-m*i, k*n-p*i, l*n-q*i,
    f*i-j*e, g*i-k*e, h*i-l*e,   s*e-t*a, u*n-v*i, t*i-u*e
  ) # did you know that "hse" stands for "holy jesus this is insane"
  w = (s - (b*hse[0] + c*hse[1] + d*hse[2])) / a
  hse.insert(0,w)
  return hse

def FiveSE(a,b,c,d,e,
           f,g,h,i,j,
           k,l,m,n,o,
           p,q,r,s,t,
           A,B,C,D,E,   P,Q,R,S,T): # 5 sim. equations
  fse = FourSE(
    b*f-g*a, c*f-h*a, d*f-i*a, e*f-j*a,
    g*k-l*f, h*k-m*f, i*k-n*f, j*k-o*f,
    l*p-q*k, m*p-r*k, n*p-s*k, o*p-t*k,
    q*A-B*p, r*A-C*p, s*A-D*p, t*A-E*p,   P*f-Q*a, Q*k-R*f, R*p-S*k, S*A-T*p
  )
  v = (P - (b*fse[0] + c*fse[1] + d*fse[2] + e*fse[3])) / a
  fse.insert(0,v)
  return fse

def SixSE(a,b,c,d,e,f,
          A,B,C,D,E,F,
          g,h,i,j,k,l,
          G,H,I,J,K,L,
          m,n,o,p,q,r,
          M,N,O,P,Q,R,   S,T,U,V,W,X): # 6 sim. equations
  ise = FiveSE(
    b*A-B*a, c*A-C*a, d*A-D*a, e*A-E*a, f*A-F*a,
    B*g-h*A, C*g-i*A, D*g-j*A, E*g-k*A, F*g-l*A,
    h*G-H*g, i*G-I*g, j*G-J*g, k*G-K*g, l*G-L*g,
    H*m-n*G, I*m-o*G, J*m-p*G, K*m-q*G, L*m-r*G,
    n*M-N*m, o*M-O*m, p*M-P*m, q*M-Q*m, r*M-R*m,   S*A-T*a, T*g-U*A, U*G-V*g, V*m-W*G, W*M-X*m
  )
  u = (S - (b*ise[0] + c*ise[1] + d*ise[2] + e*ise[3] + f*ise[4])) / a
  ise.insert(0,u)
  return ise

### inputting points

print("Point form: (x1,y1)\n")
howmanypoints = int(input("How many points do you have? (2,3,4,5,6) >>> "))
print()
p = []
for points in range(howmanypoints): # input points
  point1 = float(input("Enter x" + str(points+1) + ": "))
  point2 = float(input("Enter y" + str(points+1) + ": "))
  p.append([point1,point2])
print()

### outputting polynomial coefficients

if len(p) == 2:
  coef = TwoSE(   p[0][0], 1,
                  p[1][0], 1,   p[0][1],p[1][1]
              )
  print("a =",coef[0],"\nb =",coef[1])
elif len(p) == 3:
  coef = ThreeSE( p[0][0]**2, p[0][0], 1,
                  p[1][0]**2, p[1][0], 1,
                  p[2][0]**2, p[2][0], 1,   p[0][1], p[1][1], p[2][1]
                )
  print("a =",coef[0],"\nb =",coef[1],"\nc =",coef[2])
elif len(p) == 4:
  coef = FourSE(  p[0][0]**3, p[0][0]**2, p[0][0], 1,
                  p[1][0]**3, p[1][0]**2, p[1][0], 1,
                  p[2][0]**3, p[2][0]**2, p[2][0], 1,
                  p[3][0]**3, p[3][0]**2, p[3][0], 1,   p[0][1], p[1][1], p[2][1], p[3][1]
               )
  print("a =",coef[0],"\nb =",coef[1],"\nc =",coef[2],"\nd =",coef[3])
elif len(p) == 5:
  coef = FiveSE(  p[0][0]**4, p[0][0]**3, p[0][0]**2, p[0][0], 1,
                  p[1][0]**4, p[1][0]**3, p[1][0]**2, p[1][0], 1,
                  p[2][0]**4, p[2][0]**3, p[2][0]**2, p[2][0], 1,
                  p[3][0]**4, p[3][0]**3, p[3][0]**2, p[3][0], 1,
                  p[4][0]**4, p[4][0]**3, p[4][0]**2, p[4][0], 1,   p[0][1], p[1][1], p[2][1], p[3][1], p[4][1]
               )
  print("a =",coef[0],"\nb =",coef[1],"\nc =",coef[2],"\nd =",coef[3],"\ne =",coef[4])
elif len(p) == 6:
  coef = SixSE(   
                  p[0][0]**5, p[0][0]**4, p[0][0]**3, p[0][0]**2, p[0][0], 1,
                  p[1][0]**5, p[1][0]**4, p[1][0]**3, p[1][0]**2, p[1][0], 1,
                  p[2][0]**5, p[2][0]**4, p[2][0]**3, p[2][0]**2, p[2][0], 1,
                  p[3][0]**5, p[3][0]**4, p[3][0]**3, p[3][0]**2, p[3][0], 1,
                  p[4][0]**5, p[4][0]**4, p[4][0]**3, p[4][0]**2, p[4][0], 1,
                  p[5][0]**5, p[5][0]**4, p[5][0]**3, p[5][0]**2, p[5][0], 1,   p[0][1], p[1][1], p[2][1], p[3][1], p[4][1], p[5][1]
  )
  print("a =",coef[0],"\nb =",coef[1],"\nc =",coef[2],"\nd =",coef[3],"\ne =",coef[4],"\nf =",coef[5])

### derivative -- everything after this comment was added for my MST "artifact" after reading zero. Everything before this was for my MST distinctions from trim. 3 last year.

n = howmanypoints - 1 # assumed degree (as long as a!=0)
dcoef = [(n-i)*coef[i] for i in range(n)]
dexp = [int(n-1-i) for i in range(n)]
ans=[] # this is a list because it's easier to edit for formatting
for i in range(n): # building the derivative
  ans.append( str(dcoef[i]) )
  ans.append( "x^" + str(dexp[i]) )
  ans.append( " + " )

def der(x,dcoef,dexp): # making the function
  func=""
  for i in range(n):
    func += str(dcoef[i]) + "*x**" + str(dexp[i])
    if i != n-1:
      func += "+"
  return eval(func)

# formatting!!!

def rep(a,b,c): # replaces an item in a list with something else
  i = a.index(b)
  a.remove(b)
  a.insert(i,c)
  return a

print("\nDerivative:\n")
del ans[-1] # get rid of extra " + " at the end
del ans[-1] # get rid of "x^0" because it equals 1
ans = rep(ans,"x^1","x") # turn "x^1" into simply "x"

deriv="" # smushes "ans" into a string called "deriv"
for i in ans:
  deriv += i

while "+ -" in deriv: # if there is a plus negative
  deriv = [*deriv] # split string into list of characters to edit
  
  for i in range(len(deriv)-2):
    if deriv[i] == "+":
      if deriv[i+2] == "-": # find where it is
        del deriv[i+2]
        del deriv[i]
        deriv.insert(i,"-") # replace plus negative with minus
  
  deriv2 = ""  # smush the list of characters back into a string to check again
  for i in deriv:
    deriv2 += i
  deriv2, deriv = deriv, deriv2
      
print(deriv)

# derivative at each original point

def int2(x): # get rid of ".0" on a float (looks nicer)
  if int(x)==x:
    return int(x)
  else:
    return x

print("\nSlope at point...\n")
for i in range(howmanypoints): # find slope of function at each inputted point
  print(
    "(" + str(int2(p[i][0])) + ",", str(int2(p[i][1])) + "):", # point
    der(p[i][0],dcoef,dexp) # slope
       )
print("\nRemember, if f(x) is the original function, and g(x) is the derivative, the line that touches f(x) at any given x value n is:\ng(n)(x-n) + f(n)")

# :)