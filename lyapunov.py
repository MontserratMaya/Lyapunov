import math as m

def main():
  parametros()
  Ntot = Nr*ntau

  print("Sigma = " , sigma, ", B = ", B, ", R = ", R)
  print("Delta_t = ", Delta_t, ", Ntot = ", Ntot, ", Nr = ", Nr, ", ntau = ", ntau)
  print("Imprimieros cada ", IO, " iteraciones de ortogonalización.")
  print("CI x0 = ", ci1, ", y0 = ", ci2, ", z0 = ", ci3)
  print("Daremos ", Ntrans, " pasos para acercarnos al atractor;")
  print("Es decir, integraremos de t = 0 a t = ", Ntrans*Delta_t)
  print("Cálculo transitorio: cargando... ")

  y = []
  yLorenz = [ci1,ci2,ci3]
  yout = []
  yLorenzout = []
  f = []
  cum = []

  global t
  n = int(3)
  t = 0.0

  for i in range(1,Ntrans,1):
    RK4(yLorenz,f,yLorenzout,n,Delta_t,t)
    t += Delta_t
    yLorenz = yLorenzout
    with open("LE_trans.dat", 'w') as fich:         
          print(t, yLorenz[0], yLorenz[1], yLorenz[2],  file=fich)
  
  x0 = float(yLorenzout[0])
  y0 = float(yLorenzout[1])
  z0 = float(yLorenzout[2])
  t = 0.0

  print("Fin del ciclo transitorio. Iniciaremos el cálculo en: ")
  print("t = ", t, "x0 = ", x0, "y0 = ", y0, "z0 = ", z0)

  n = int(12)

  y[0] = x0
  y[1] = y0
  y[2] = z0
  for i in range(3,11,1):
    y[i] = float(0)
  y[3] = float(1.0)
  y[7] = float(1.0)
  y[11] = float(1.0)
  cum[0] = float(0.0)
  cum[1] = float(0.0)
  cum[2] = float(0.0)
  v1 = []
  v2 = []
  v3 = []
  V1 = []
  V2 = []
  V3 = []
  Normas = []

  print("Se integrarán las ", n, "ecuaciones de t=0 a t = ", Ntot*Delta_t)

  for i in range(1,Nr,1):
    for j in range(1,ntau,1):
      RK4(y,f,yout,n,Delta_t,t)
      y = yout
      t += Delta_t
    for k in range(0,2,1):
      v1[k] = float(y[3*k + 1])
      v2[k] = float(y[3*k + 2])
      v3[k] = float(y[3*k + 3])
    OGS(v1,v2,v3,V1,V2,V3,Normas)
    for k in range(0,2,1):
      cum[k] = float(cum[k]) + m.log(float(Normas[k]))
    if i%IO == 0 | i == Nr:
      with open("LE.dat", 'w') as fich:         
          print(t-Delta_t, float(cum[0])/t, float(cum[1])/t, float(cum[2])/t,  file=fich)
    for k in range(0,2,1):
      y[3*k + 1] = float(V1[k])
      y[3*k + 2] = float(V2[k])
      y[3*k + 3] = float(V3[k])

def parametros():

  global sigma, B, R, Delta_t, Nr, ntau, IO, Ntrans,ci1,ci2,ci3


  sigmas, Bs, Rs = input("Parametros SIGMA,B,R: ").split()
  sigma = float(sigmas)
  B = float(Bs)
  R = float(Rs)
  Delta_t = float(input("Escriba el paso de tiempo Delta_t "))
  Nrs, ntaus = input("Introduzca Nr: Num de ortogonalizaciones y ntau: Num de pasos de tiempo entre ortogonalizaciones: ").split()
  Nr = int(Nrs)
  ntau = int(ntaus)
  print("A cada cuantos pasos imprimimos ")
  IO = int(input("Escriba un entero IO tal que Nr/IO sea entero "))
  if Nr%IO != 0:
    print("IO debe ser tal que Nr/IO es cero y no lo es, por favor introduzca de nuevo...")
  x00s, y00s, z00s = input("Condiciiones iniciales x0,y0,z0: ").split()
  ci1 = float(x00s)
  ci2 = float(y00s)
  ci3 = float(z00s)
  Ntrans = int(input("Escriba el numero de pasos dados para el edo transitorio: "))

def DerivadasLorentz(varModelo, funModelo):
  global X,Y,Z

  X = varModelo[0]
  Y = varModelo[1]
  Z = varModelo[2]

  funModelo[0] = sigma*(Y-X)
  funModelo[1] = R*X - Y - X*Z
  funModelo[2] = X*Y - B*Z

def Derivadas(var,func):
  func.append(sigma*(var[1]-var[0]))
  func[1] = R*var[0] - var[1] - var[0]*var[2]
  func[2] = var[0]*var[1] - B*var[2]

  for i in range(0,2):
    func[3+i] = sigma*(var[6+i] - var[3+i])
    func[6+i] = (R - var[2])*var[3+i] - var[6+i] - var[0]*var[9+i]
    func[9+i] = var[1]*var[3+i] + var[0]*var[6+i] - B*var[9+i]


    

def norma(vector):
    return m.sqrt(sum(pow(element, 2) for element in vector))

def prod(vector1, vector2):
    for j in range(0,2):
      dot = float(vector1[j])*float(vector2[j])
    return dot
    

def OGS(v1,v2,v3,V1,V2,V3,N):
  N[0] = norma(v1)
  for k in range(0,2,1):
    V1[k] = v1[k]/N[0]
  
  for k in range(0,2,1):
    V2[k] = v2[k] - prod(v2,V1)*V1[k] 
  N[1] = norma(V2)
  for k in range(0,2,1):
    V2[k] = V2[k]/N[1]
  
  for i in range(0,2,1):
    V3[i] = v3[i] - prod(v3,V2)*V2[i] - prod(v3,V1)*V1[i]
  N[2] = norma(V3)
  for i in range(0,2,1):
    V3[i] = V3[i]/N[2]
  

def RK4(y,f,yout,n,h,t):
  ytemp = [] 
  k1 = []
  k2 = []
  k3 = []
  k4 = []

  Derivadas(t,y,f)
  for j in range(0,n-1,1):
    k1[j] = h*f[j]
    ytemp[j] = y[j] + 0.5*k1[j]
  
  Derivadas(t+0.5*h,ytemp,f)
  for j in range(0,2,1):
    k2[j] = h*f[j]
    ytemp[j] = y[j] + 0.5*k2[j]

  Derivadas(t+0.5*h,ytemp,f)
  for j in range(0,2,1):
    k3[j] = h*f[j]
    ytemp[j] = y[j] + k2[j]

  Derivadas(t+h,ytemp,f)
  for j in range(0,2,1):
    k4 = h*f[j]

  for j in range(0,2,1):
    yout = y[j] + (1.0/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])
  
  del ytemp, k1, k2, k3, k4
  return yout




main()