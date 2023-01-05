import matplotlib.pyplot as plt
from sympy import *
import sys
import numpy as np
import ast


def string_to_expression(expression_string):
    x = symbols('x')
    y = symbols('y')
    return sympify(expression_string, locals={'x': x, 'y': y})


def convert_to_array(string_array):
    return ast.literal_eval(string_array)


x, y = symbols("x y")
variables = [x, y]
strtpoints = convert_to_array(sys.argv[2])


# code for obtaining hessian matrix

hessian = np.zeros((2, 2))
print(hessian)

for i in range(0, hessian.shape[0]):
    for j in range(0, hessian.shape[1]):
        hessian[i, j] = diff(diff(sys.argv[1], variables[i]), variables[j]).subs(
            x, strtpoints[0]).subs(y, strtpoints[1])
        print(hessian[i, j])
print(hessian)
#####################################


# calculate gradient at some point
def calculate_deltafunc(point, expr):
    pointx = point[0, 0]
    pointy = point[1, 0]
    x, y = symbols('x y')
    grad = np.zeros((2, 1))
    dfbydx = diff(expr, x).subs(x, pointx).subs(y, pointy)
    dfbydy = diff(expr, y).subs(y, pointy).subs(x, pointx)
    print("This function dfbydy", dfbydx)
    grad[0, 0] = dfbydx
    grad[1, 0] = dfbydy
    return grad


#######################################
hessian = np.linalg.inv(hessian)
print(hessian)


# start the algorithm and do newton raphton

def newtonrapthon(strtpoints):

    xn = np.array(strtpoints)
    print(xn)
    xn = xn.reshape(2, 1)
    print('xn is shape ', xn.shape)
    i = 0
    plotx1 = []
    while (True and (i != 30)):
        deltafunc = calculate_deltafunc(xn, sys.argv[1])
        isZero = np.array_equal(deltafunc, np.zeros((2, 1)))
        xn1 = xn.reshape(2, 1) - np.dot(hessian, deltafunc.reshape(2, 1))
        plotx1.append(deltafunc[0, 0])
        if (isZero):
            print("deltafucn is zero")
            break
        xn = xn1

        print(xn)
        i += 1
    return xn1, plotx1

xn1,plotx1 = newtonrapthon(strtpoints)
print("The extreme points are ", xn1)

if(i==30):
    print("The initial point is not valid inital point")
    
def gradentDescent(X0,iexpr):
    x, y, Y = symbols('x y Y')
    expr = string_to_expression(iexpr)   #x-y+2*x**2+2*x*y+y**2#(x-2)2 + (y-2)*2
    dx = expr.diff(x)
    Xaxis = []
    Yaxis = []
    dy =expr.diff(y)

    F = np.array([[dx.subs([(x, X0[0][0]), (y, X0[1][0])])], [dy.subs([(x, X0[0][0]), (y, X0[1][0])])]]) 

    Yaxis.append(F[0][0])

    X1 = X0+Y*F*(-1)

    Xi = solve((expr.subs([(x, X1[0][0]), (y, X1[1][0])])).diff(Y))  #GET THE VALUE OF LAMDA AFTER DIFFERENTIATION

    X1 = np.array([[X1[0][0].subs(Y,Xi[0])], [X1[1][0].subs(Y,Xi[0])]]) 

    err= X1 - X0


    i=0
    Xaxis.append(i)
    while round(F[0][0],2)!=0 and round(F[1][0],2)!=0 :
        X0 = X1
        F = np.array([[dx.subs([(x, X0[0][0]), (y, X0[1][0])])], [dy.subs([(x, X0[0][0]), (y, X0[1][0])])]]) 
        if((F[0][0])==0 and (F[1][0])==0):
            break
        else:
            Yaxis.append(F[0][0])
            X1 = X0+Y*F*(-1)

            Xi = solve((expr.subs([(x, X1[0][0]), (y, X1[1][0])])).diff(Y))  #GET THE VALUE OF LAMDA AFTER DIFFERENTIATION
     
            X1 = np.array([[X1[0][0].subs(Y,Xi[0])], [X1[1][0].subs(Y,Xi[0])]])
     
            err = X1-X0
            i=i+1
            Xaxis.append(i)

    plt.plot(Xaxis, Yaxis, color='r',label = "Gradent Descent")
  
    # naming the x axis
    plt.xlabel('Number of Iterations')
    # naming the y axis
    plt.ylabel('Convergence Criteria')
  
    # giving a title to my graph
    plt.title('Problem 1')

    # show a legend on the plot
    plt.legend()
  
    # function to show the plot
    plt.show()

    print("The extreme points are ", X1)



#X1 = np.array([[0], [0]]) 

#X2 = np.array([[0], [1]]) 
X3 = np.array([[1], [0]]) 
# gradentDescent(X1)
# gradentDescent(X2)
# gradentDescent(X3)    
def draw(plotx1,X1):    
    plt.figure()
    plt.plot(plotx1,color='b',label = "Newton Raphthon")
    gradentDescent(X1,sys.argv[1])
    

draw(plotx1,X3)    

