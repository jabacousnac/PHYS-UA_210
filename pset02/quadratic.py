import numpy as np

def x(a,b,c):
    xp = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    xm = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    print([xp, xm])
    xp = ((2*c)/(-b - np.sqrt(b**2 - 4*a*c)))
    xm = ((2*c)/(-b + np.sqrt(b**2 - 4*a*c)))
    print([xp, xm])

def x_accurate(a,b,c):
    xp = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    xm = ((2*c)/(-b - np.sqrt(b**2 - 4*a*c)))
    print([xp, xm])
    
if __name__ == "__main__":
    x(.001, 1000, .001)
    print("\n")
    x_accurate(.001, 1000, .001)
