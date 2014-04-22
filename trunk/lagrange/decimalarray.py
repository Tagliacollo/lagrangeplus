from decimal import *

class DecimalArray (list): 
    def __init__(self, listOfNumbers):
        junk = map( self.append, map(Decimal, listOfNumbers) )

    def __add__(self, b):
        ret = DecimalArray(self)
        for i, a in enumerate(self):
            ret[i] = a + b[i]

        return ret

    #def __mul__(self, b):
        #ret = DecimalArray(self)
        #for i, a in enumerate(self):
            #ret[i] = a * b[i]
        #return ret

    def __mul__(self,b):
        ret = list()
        for i,a in enumerate(self):
            ret.append( a * b[i] )
        return ret


if __name__ == "__main__":
    a = DecimalArray([1,2,3])
    b = a + a
    print b

    c = b * b
    print c
    
    print sum(c)
