import numpy as np
import matplotlib.pyplot as plt

def pos_errors(array):
    return [v.plus for v in array]

def neg_errors(array):
    return [v.minus for v in array]

class AsymmetricUncertainty:
    """
    Class for handling propagation of asymmetric uncertainties assuming a pseudo-Gaussian probability distribution
    where the plus and minus errors in each direction of the nominal value are like modified 1-sigma standard devations.
    See 2008ApJ...672..433S for more.
    
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    nominal : numeric
        the nominal value of the represented quantity
    pos_err : numeric
        the plus error on the value, as in value = nominal (+pos_err, -neg_err)
    neg_err : numeric
        the minus error on the value, as in value = nominal (+pos_err, -neg_err)

    Methods
    -------
    pdf()
        plots the probability distribution function
    
    cdf()
        plots the cumulative distribution function

    Raises
    ------
    ValueError
        if you don't pass compatible values for the arguments
    TypeError
        if you don't pass values of the correct format for the arguments
    """
    
    def __init__(self, nominal, pos_err=0, neg_err=0):
        if isinstance(nominal,str):
            stripped = nominal.replace(" ","")
            if "±" in stripped:
                self.value = float(stripped.split("±")[0])
                self.plus = self.minus = float(nominal.split("±")[1])
            elif "+" in stripped and "-" in stripped:
                self.value = float(stripped.split("(")[0])
                err_str = stripped.split("(")[1].replace(")","")
                self.plus = float(err_str.split(",")[0][1:])
                self.minus = float(err_str.split(",")[1][1:])
            else:
                raise TypeError("Failed to parse string, likely due to improper formatting.")
        else:
            self.value = float(nominal)
            self.plus = np.abs(float(pos_err))
            self.minus = np.abs(float(neg_err))
        self.maximum = self.value+self.plus
        self.minimum = self.value-self.minus
        self.sign = 1 if self.value >= 0 else -1
        self.is_symmetric = np.isclose(self.plus, self.minus)
        
    def __str__(self):
        if np.isclose(self.plus, self.minus):
            return "{} ± {}".format(self.value,self.plus)
        else:
            return "{} (+{}, -{})".format(self.value,self.plus,self.minus)
        
    def _repr_latex_(self):
        if np.isclose(self.plus, self.minus):
            return "$%f \pm %f$" %(self.value,self.plus)
        else:
            return "$%f_{-%f}^{+%f}$" %(self.value,self.minus,self.plus)
        
    def pdf(self,x):
        return np.piecewise(x, [x<self.value, x>=self.value],
                            [lambda x : np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(x-self.value)**2 / (2*self.minus**2)),
                             lambda x : np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(x-self.value)**2 / (2*self.plus**2))])
    
    def cdf(self,x):
        return np.cumsum(self.pdf(x))/np.sum(self.pdf(x))
        
    def pdfplot(self,num_sigma=5,discretization=100,**kwargs):
        neg_x = np.linspace(self.value-(num_sigma*self.minus),self.value,discretization)
        pos_x = np.linspace(self.value,self.value+(num_sigma*self.minus),discretization)
        p_neg = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(neg_x-self.value)**2 / (2*self.minus**2))
        p_pos = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(pos_x-self.value)**2 / (2*self.plus**2))
        x = np.array(list(neg_x)+list(pos_x))
        pdf = self.pdf(x)
        plt.plot(x,pdf,**kwargs)
        plt.show()
        
    def cdfplot(self,num_sigma=5,discretization=100,**kwargs):
        neg_x = np.linspace(self.value-(num_sigma*self.minus),self.value,discretization)
        pos_x = np.linspace(self.value,self.value+(num_sigma*self.minus),discretization)
        p_neg = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(neg_x-self.value)**2 / (2*self.minus**2))
        p_pos = np.sqrt(2)/np.sqrt(np.pi)/(self.plus+self.minus) * np.exp(-1*(pos_x-self.value)**2 / (2*self.plus**2))
        x = np.array(list(neg_x)+list(pos_x))
        pdf = self.pdf(x)
        cdf = np.cumsum(pdf)/np.sum(pdf)
        plt.plot(x,cdf,**kwargs)
        plt.show()
        
    def add_error(self, delta, method="quadrature", inplace=False):
        if method=="quadrature":
            new_pos = np.sqrt(self.plus**2 + delta**2)
            new_neg = np.sqrt(self.minus**2 + delta**2)
        elif method=="straight":
            new_pos = self.plus + delta
            new_neg = self.minus + delta
        elif method=="split":
            new_pos = self.plus + delta/2
            new_neg = self.minus + delta/2
        else:
            raise ValueError
        if inplace:
            self.plus = new_pos
            self.minus = new_neg
        else:
            return AsymmetricUncertainty(self.value,new_pos,new_neg)
    
    def items(self):
        return (self.value,self.plus,self.minus)
    
    def __int__(self):
        return int(self.value)
    
    def __float__(self):
        return float(self.value)
    
    def __neg__(self):
        return AsymmetricUncertainty(-self.value,self.minus,self.plus)
        
    def __add__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value + other.value
        pos = np.sqrt(self.plus**2 + other.plus**2)
        neg = np.sqrt(self.minus**2 + other.minus**2)
        #print("added",self,"+",other,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __radd__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value + other.value
        pos = np.sqrt(self.plus**2 + other.plus**2)
        neg = np.sqrt(self.minus**2 + other.minus**2)
        #print("added",other,"+",self,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __sub__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value - other.value
        pos = np.sqrt(self.plus**2 + other.minus**2)
        neg = np.sqrt(self.minus**2 + other.plus**2)
        #print("subtracted",other,"from",self,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rsub__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = other.value - self.value
        pos = np.sqrt(self.minus**2 + other.plus**2)
        neg = np.sqrt(self.plus**2 + other.minus**2)
        #print("subtracted",self,"from",other,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __mul__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        
        result = self.value * other.value
        pos = np.sqrt((self.plus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
        neg = np.sqrt((self.minus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
        #print("multiplied",self,"by",other,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rmul__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        
        result = self.value * other.value
        pos = np.sqrt((self.plus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
        neg = np.sqrt((self.minus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
        #print("multiplied",other,"by",self,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __truediv__(self,other): # self divided by something
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value / other.value
        pos = np.sqrt((self.plus/self.value)**2 + (other.minus/other.value)**2) * np.abs(result)
        neg = np.sqrt((self.minus/self.value)**2 + (other.plus/other.value)**2) * np.abs(result)
        #print("divided",self,"by",other,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rtruediv__(self,other): # something divided by self
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = other.value / self.value
        pos = np.sqrt((other.plus/other.value)**2 + (self.minus/self.value)**2) * np.abs(result)
        neg = np.sqrt((other.minus/other.value)**2 + (self.plus/self.value)**2) * np.abs(result)
        #print("divided",other,"by",self,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __pow__(self,other): # self to the something power
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = self.value**other.value
        pos = np.abs(result)*np.sqrt((self.plus*other.value/self.value)**2 + (other.plus*np.log(self.value))**2)
        neg = np.abs(result)*np.sqrt((self.minus*other.value/self.value)**2 + (other.minus*np.log(self.value))**2)
        #print("raised",self,"to",other,"=",AsymmetricUncertainty(result,pos,neg))        
        return AsymmetricUncertainty(result,pos,neg)
    
    def __rpow__(self,other): # something to the self power
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        result = other.value**self.value
        pos = np.abs(result)*np.sqrt((other.plus*self.value/other.value)**2 + (self.plus*np.log(other.value))**2)
        neg = np.abs(result)*np.sqrt((other.minus*self.value/other.value)**2 + (self.minus*np.log(other.value))**2)
        #print("raised",other,"to",self,"=",AsymmetricUncertainty(result,pos,neg))                
        return AsymmetricUncertainty(result,pos,neg)
    
    def log10(self):
        result = np.log10(self.value)
        pos = self.plus/(self.value*np.log(10))
        neg = self.minus/(self.value*np.log(10))
        #print("logged",self,"=",AsymmetricUncertainty(result,pos,neg))
        return AsymmetricUncertainty(result,pos,neg)
    
    def log(self):
        result = np.log(self.value)
        pos = self.plus/self.value
        neg = self.minus/self.value
        return AsymmetricUncertainty(result,pos,neg)
    
    def sqrt(self):
        return self**0.5
    
    def __eq__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value == other.value and self.plus == other.plus and self.minus == other.minus
    
    def __gt__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value > other.value
    
    def __lt__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value < other.value
    
    def __lshift__(self,other): # overloaded <<; definitively less than
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.maximum < other.minimum
    
    def __rshift__(self,other): # overloaded >>; definitively greater than
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.minimum > other.maximum
    
    def __le__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value <= other.value
    
    def __ge__(self,other):
        if isinstance(other,type(self)):
            pass
        else:
            other = AsymmetricUncertainty(other,0,0)
        return self.value >= other.value
    
    def conjugate(self):
        return self.value
    
    def __isfinite__(self):
        return all(np.isfinite(self.items()))
    
    def isna(self):
        return pd.isna(self.value)
    
    def notna(self):
        return pd.notna(self.value)

class UncertaintyArray(list):
    
    def refresh(self):
        for i in range(len(self)):
            try:
                self[i].value
                self[i].plus
                self[i].minus
            except AttributeError:
                self[i] = AsymmetricUncertainty(self[i],0,0)
                
        self.as_numpy = np.array(self.as_list)
        self.flattened = self.as_numpy.flatten()
        self.shape = self.as_numpy.shape
        self.ndim = self.as_numpy.ndim

        self.minus = [v.minus for v in self.as_list]
        self.plus = [v.plus for v in self.as_list]
        self.values = [v.value for v in self.as_list]
    
    def __init__(self,array=[]):

        self.as_list = list(array)
        self.refresh()

    def __len__(self):
        return len(self.as_list)

    def __iter__(self):
        return iter(self.as_list)

    def __getitem__(self, key):
        return self.as_list[key]

    def __setitem__(self, key, val):
        self.as_list[key] = val

    def __str__(self):
        return str([str(entry) for entry in self.as_list])

    def __repr__(self):
        return str(self)

    def __contains__(self, item):
        return item in self.as_list

    def append(self,entry):
        self.as_list.append(entry)
        self.refresh()
    
    def pdf(self,x):
        return np.sum([entry.pdf(x) for entry in self], axis=0)