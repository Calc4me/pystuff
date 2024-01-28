import constants as c
from scipy.integrate import quad
import math

#Max factorial is 998 int, and around 13.49699996 for floats

def abs(n):
  if n < 0:
    return -n
  else:
    return n
  
def circle_area(radius):
    return c.PI * float(radius**2)

def circ(r):
    return 2.0 * c.PI * float(r)

def fac(n):
  if n == 0:
    return 1
  if int(n) == n:
    return n * fac(n-1)
  if n < 0:
    n = n * -1
  int_result = quad(lambda x: (-math.log(x))**n, 0.0, 1.0)
  intmul = round(int_result[0],4)
  intpw = round(int_result[1],4)
  return round(intmul * 10.0**intpw,4)

def power(base, exponent):
  return base ** exponent

def sqrt(n):
  return n**(1/2)

def crt(n):
  return n**(1/3)

#Finds c from a and b
def hypot(a,b):
  sqrt(a^2+b^2)

def rad_to_deg(n):
  return n * 180/c.PI

def deg_to_rad(n):
  return n * c.PI/180

def convdeg(n,unit):
  if unit == "deg":
    return n * c.PI/180
  elif unit == "rad":
    return n * 180/c.PI
  else:
    print("Invalid unit")

def area(a,b):
  return a*b

def checkprime(n):
  if n <= 1:
    return False
  for i in range(2, int(n**0.5)+1):
    if n % i == 0:
      return False
  return True

def print_prime(n):
  print(checkprime(n))

def heron(a,b,c):
  s = (a+b+c)/2
  area = sqrt(s*(s-a)*(s-b)*(s-c))
  return area

def sphere(rad):
  return 4/3*c.PI*rad**3

def gcd(a,b):
  if a == 0:
    return b
  return gcd(b % a, a)

def lcm(a,b):
  return (a * b) // gcd(a, b)

def convert(n,curunit,convunit):
  if curunit == "cel":
    if convunit == "fah":
      return n * 9/5 + 32
    elif convunit == "kel":
      return n + 273.15
    else:
      print("Invalid unit")
  if curunit == "fah":
    if convunit == "cel":
      return (n - 32) * 5/9
    elif convunit == "kel":
      return (n - 32) * 5/9 + 273.15
    else:
      print("Invalid unit")
  if curunit == "kel":
    if convunit == "cel":
      return n - 273.15
    elif convunit == "fah":
      return (n - 273.15) * 9/5 + 32
    else:
      print("Invalid unit")
  else:
    print("Invalid unit")

def perm(n,k):
  return fac(n)/fac(n-k)

def comb(n,k):
  return fac(n)/fac(n-k)*fac(k)

def fibb(n):
  if n <= 1:
    return n
  else:
    return fibb(n-1) + fibb(n-2)
  
def amean(l):
  temp = 0
  for i in range(len(l)):
    temp = temp + l[i]
  return temp/len(l)

def gmean(l):
  temp = 1
  for i in range(len(l)):
    temp = temp*l[i]
  return temp**(1/len(l))

def hmean(l):
  temp = 0
  for i in range(len(l)):
    temp = temp + 1/l[i]
  return len(l)/temp

def median(l):
  l = sorted(l)
  if len(l) % 2 == 1:
    return l[int((len(l)-1)/2)]
  else:
    return ((l[int(len(l)/2)])+(l[int((len(l)/2)-1)]))/2
  
def mode(l):
  l = sorted(l)
  templist = []
  numlist = []
  testvar = 0
  maxvar = 0
  highscore = 0
  #Get list of all numbers in the list
  for i in range(len(l)):
    if l[i] not in templist:
      templist.append(l[i])
  for i in range(len(templist)):
    #Number to be measured
    testvar = templist[i]
    maxvar = l.count(testvar)
    #Is it the most measured?
    if maxvar == highscore:
      numlist.append(testvar)
    elif maxvar >= highscore:
      highscore = maxvar
      numlist.clear()
      numlist.append(testvar)
  if len(numlist) == 1:
    return numlist[0]
  elif len(numlist) > 1:
    for i in range(len(numlist)):
      return amean(numlist)
    
def variance(l,ddof=0):
  n = len(l)
  mean = sum(l) / n
  return sum((x - mean) ** 2 for x in l) / (n - ddof)

def sdev(l,ddof=0):
  return sqrt(variance(l,ddof))

def eucliddist(x1,y1,x2,y2):
  return sqrt((abs(x1-x2)**2)+(abs(y1-y2)**2))

def eucliddist(p1,p2):
  return sqrt((abs(p1[0]-p2[0])**2)+(abs(p1[1]-p2[1])**2))

def A(m,n,s ="% s"):
  print(s % ("A(% d, % d)" % (m, n)))
  if m == 0:
      return n + 1
  if n == 0:
      return A(m - 1, 1, s)
  n2 = A(m, n - 1, s % ("A(% d, %% s)" % (m - 1)))
  return A(m - 1, n2, s)

def atri(b,h):
  return (b*h)/2

def simpleint(P,r,t):
  return P*r*t

def compint(P,r,n,t):
  return P*((1 + r/n)**(n*t))

def perimiter(w,h):
  return 2*(w+h)

#Finds a or b from a&c/b&c
def pythag(n1,n2):
  if n1 > n2:
    return sqrt((n1**2)-(n2**2))
  elif n2 > n1:
    return sqrt((n2**2)-(n1**2))
  elif n1 == n2:
    return None
  
def slope(p1,p2):
  return  (p1[1]-p2[1])/(p1[0]-p2[0])

def volcube(s):
  return s^3

def percent(tot,part):
  return 100*(part/tot)

def slopeint(x1,y1,slope):
  return [1*slope, (slope*((-1)*x1))+y1]

def quadratic(a,b,c):
  determ = b**2 - 4*a*c
  if determ >= 0:
    x1 = (-b+sqrt(determ))/2*a
    x2 = (-b-sqrt(determ))/2*a
  elif determ < 0:
    x1 = str((sqrt(determ*(-1)))/2) + "i" + " + " + str((-b)/2)
    x2 = str((sqrt(determ*(-1)))/2) + "i" + " + " + str(b/2)
  return [x1,x2]

def cyl(r,h):
  return h*(circle_area(r))

def surfacearea(a,b,c):
  return 2*((a*b)+(a*c)+(b*c))

def areasec(radius, angle): 
  if angle >= 360: 
    raise Exception("Angle not possible")
  else: 
    return circle_area(radius)*(angle/360) 
  
def conearea(r,h):
  return (1/3)*circle_area(r)*h

def dist3d(p1,p2):
  return sqrt(((p1[0]-p2[0])**2)+((p1[1]+p2[1])**2)+((p1[2]+p2[2])**2))

def slopquad(a,b,c,x):
  return (2*a*x)+b+c*0

def averageslope(func,start,end):
  return slope([start,func(start)],[end,func(end)])

def divisors(n):
  divs = 0
  for i in range(1,n+1):
    if n % i == 0:
      divs += 1
  return divs

def mobius(n):
  divs = []
  if n == 1:
    return 1
  else:
    for i in range(1,n+1):
      if n % i == 0:
        divs.append(i)
    for i in range(len(divs)):
      if checkprime(divs[i]) == True:
        primedivs += 1
      if sqrt(divs[i]) == int(sqrt(divs[i])) and divs[i] != 1:
        return 0
    return (-1)**primedivs

def zeta(n, max):
  temp = 0
  for i in range(1,max+1):
    temp += 1/((i)**n)
  return temp
def collatz(n):
  values = []
  if n == 1:
    values.append(n)
    return values
  elif n > 1:
    while n > 1:
      if n % 2 == 0:
        values.append(n)
        n = int(n/2)
      elif n % 2 == 1:
        values.append(n)
        n = n*3 + 1
    values.append(n)
    return values
def aad(l):
  temp = 0
  for i in range(len(l)):
    temp = temp + (abs(l[i]-amean(l)))
  return temp/len(l)

class fraction:
  def __init__(self, numerator, denominator):
    self.n = numerator
    self.d = denominator
    self.float = self.n / self.d
    
  def __trunc__(self):
    return fraction(math.trunc(self.n/self.d), 1)
  def __ceil__(self):
    return fraction(math.ceil(self.n/self.d), 1)
  def __floor__(self):
    return fraction(math.floor(self.n/self.d), 1)
  def __round__(self, ndigits=0):
    return fraction(round(self.n/self.d, ndigits) * (10**ndigits), 10**ndigits)
  def __abs__(self):
    return fraction(abs(self.n), abs(self.d))
  def __neg__(self):
    return fraction(-self.n, self.d)
  def __add__(self, other):
    return fraction(self.n * other.d + self.d * other.n, self.d * other.d)
class circle:
  def __init__(self,radius,x,y):
    self.rad = radius
    self.posx = x
    self.posy = y
    self.area = c.PI*(self.rad)^2
    self.circum = 2*c.PI*self.rad
    self.diam = 2*self.rad
  def circ(self):
    return 2.0 * c.PI * float(self.rad)
  def circle_area(radius):
    return c.PI * float(radius^2)
  def print_circle_data(radius_from_user):
    circumference = circ(radius_from_user)
    area = circle_area(radius_from_user)
    print('The circumfrence of the circle is ' + str(circumference) + '.')
    print('The area of the circle is ' + str(area) + '.')
