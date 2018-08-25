__all__ = [
    'abacus',
    'replace_by_aliases',
    ]

from calla.html import *
    
def replace_by_aliases(expression, aliases):
    """
    Replace parameters in expression by their aliases.
    e.g. 'alpha*beta-gamma' -> 'α*β-γ'

    >>> replace_by_aliases('(a*3+b)/c',{'a':'A','b':'B'})
    '(A*3+B)/c'
    """
    operators = ('+','-','*','/','=','(',')','{','}','&','|','\\','·',':',' ')
    s = expression
    i = j = 0
    while i < len(s):
        if s[i].isalpha():
            j = i+1
            while j <= len(s):
                if j == len(s) or s[j] in operators:
                    attr = s[i:j]
                    if attr in aliases:
                        result = aliases[attr]
                        if isinstance(result, tuple):
                            result = result[0]
                        length = len(result)
                        s = s[:i] + result + s[j:]
                        i += length
                    else:
                        i = j
                    break
                else:
                    j += 1
        else:
            i += 1
    return s
    
class abacus:
    """
    Base class for all calculators.
    Initialize data and generate html format report.
    """
    
    def __init__(self, **kwargs):
        """
        Initialize parameters from '__inputs__'.
        parameters format: (parameter, (symbol, unit, default_value, name, description[, choices]))
        e.g. __inputs__ = OrderedDict((
            ('Es',('E<sub>s</sub>','MPa',2.0E5,'Elastic modulus of rebar')),
            ))
        'alias' is usually in html style that can be displayed better in browser.
        'choices' is optional, valid for multi-select parameters.
        """
        if hasattr(self, '__inputs__'):
            for k in self.__inputs__:
                if not hasattr(self, k):
                    infos = self.__inputs__[k]
                    if isinstance(infos, tuple):
                        v = infos[2] if len(infos)>2 else 0
                        setattr(self,k,v)
        # initialize parameters by kwargs
        for key in kwargs:
            if hasattr(self, key):
                setattr(self, key, kwargs[key])
##        if hasattr(self, '__inputs__'):
##            for key in kwargs:
##                if key in self.__inputs__:
##                    infos = self.__inputs__[key]
##                    if isinstance(infos, tuple):
##                        v = infos[2] if len(infos)>2 else 0
##                        setattr(infos,k,v)

##    def __getattr__(self, name):
##        if hasattr(self.__inputs__, name):
##            return getattr(self.__inputs__, name)
##        return super().__getattr__(name)
##
##    def __setattr__(self, name, value):
##        if hasattr(self.__inputs__, name):
##            setattr(self.__inputs__, name, value)
##        return super().__setattr__(name, value)

    @property
    def parameters(self):
        return { attr:(getattr(self, attr) if hasattr(self, attr) else 0) for attr in self.__inputs__}

    @parameters.setter
    def parameters(self, values):
        for key in values:
            if hasattr(self, key):
                setattr(self, key, values[key])

    def disableds(self):
        """
        Get parameters disabled.
        self.__disableds__ format:{parameter:{option:(disabled_parameters),...},...}
        e.g. __disableds__ = {
        'option':{0:('wmax'),1:('As')},
        'force_type':{0:('l0','Nq'),1:(),2:('l0'),3:('Mq','l0')}
        }

        >>> a=abacus()
        >>> a.option=1
        >>> a.__toggles__={'option':{0:('wmax'),1:('As')}}
        >>> a.disableds()
        ['As']
        """
        r = []
        if not hasattr(self,'__toggles__'):
            return r
        for key in self.__toggles__:
            v = getattr(self,key)
            items = self.__toggles__[key][v]
            if isinstance(items,tuple):
                for item in items:
                    r.append(item)
            else:
                r.append(items)
        return r
                        
    def formatI(self, name, digits=None, sep='：'):
        """Format Input parameters"""
        info = self.__inputs__[name]
        s = '{}{}'.format(info[3],sep) if (len(info)>3 and info[3] != '') else ''
        v = getattr(self,name)
        if digits != None:
            v = '{1:.{0}f}'.format(digits, v)
        s += '{} = {} {}'.format(info[0],v,info[1])
        return s
    
    def formatD(self, name, value=None, digits = 2, sep='：'):
        """Format Derived parameters, keep given digits after dicimal point."""
        info = self.__deriveds__[name]
        s = '{}{}'.format(info[3],sep) if (len(info)>3 and info[3] != '') else ''
        if value==None:
            value = getattr(self,name)
        if digits != None:
            v = '{1:.{0}f}'.format(digits, value)
        s += '{} = {} {}'.format(info[0],v,info[1])
        return s
    
    def formatX(self,name,digits=2, sep='：'):
        """Format parameters，choosing diferent format automatically."""
        if name in self.__inputs__:
            return self.formatI(name,digits,sep)
        elif name in self.__deriveds__:
            return self.formatD(name,None,digits,sep)
        else:
            return name
        
    def symbol(self,x):
        """get alias of x"""
        a = ''
        if x in self.__inputs__:
            a = self.__inputs__[x]
        elif x in self.__deriveds__:
            a = self.__deriveds__[x]
        else:
            raise 'Alias not exists.'
        if isinstance(a, tuple):
            a = a[0]
        return a
        
    def express(self, expression, digits = 2, use_symbols = True, include_self = True):
        """translate expression to value formula style.
        Args:
            use_aliases:使用别名
            include_self: 包含表达式本身
            
        >>> a = abacus()
        >>> a.a=1.012
        >>> a.b=2.375
        >>> a.c = 3.55
        >>> a.__inputs__ = {'a':'A','b':'B'}
        >>> exp = '(a*3+b)/c'
        >>> a.express(exp)
        '(A*3+B)/c = (1.01*3+2.38)/3.55'
        """
        result = replace_by_aliases(expression, self.__inputs__) if (hasattr(self, '__inputs__') and use_symbols) else expression
        if hasattr(self, '__deriveds__'):
            result = replace_by_aliases(result, self.__deriveds__)
        num_formula = self.replace_by_values(expression,digits)
        if include_self:
            return '{} = {}'.format(result, num_formula)
        return num_formula
    
    def replace_by_values(self, expression, digits = 2):
        """
        Replace parameters in expression by their values.
        e.g. 'a*b-c' -> '3*2-4'
        """
        operators = ('+','-','*','/','(',')','{','}',' ')
        s = expression
        i = j = 0
        while i < len(s):
            if s[i].isalpha():
                j = i+1
                while j <= len(s):
                    if j == len(s) or s[j] in operators:
                        attr = s[i:j]
                        if hasattr(self,attr):
                            v = getattr(self,attr)
                            if type(v) is float:
                                result = '{0:.{1}f}'.format(v,digits)
                            else:
                                result = str(v)
                            length = len(result)
                            s = s[:i] + result + s[j:]
                            i += length
                        else:
                            i = j
                        break
                    else:
                        j += 1
            else:
                i += 1
        return s
    
    def replace_by_symbols(self, expression):
        """
        Replace parameters in expression by their symbols.
        e.g. 'alpha*beta-gamma' -> 'α*β-γ'
        """
        if hasattr(self, '__inputs__'):
            s = replace_by_aliases(expression, self.__inputs__)
        if hasattr(self, '__deriveds__'):
            s = replace_by_aliases(s, self.__deriveds__)
        return s
    
    def _html(self, digits = 2):
        """Subclasses should implement this method in order to output html format reports."""
        yield 'Method for generating html is not implemented.'
        
    def html(self, digits = 2):
        result = ''
        gen = self._html(digits)
        for p in gen:
            result += '<p>' + p + '</p>'
        gen.close()
        return result

    def text(self, digits = 2, ignore_sub=True):
        return html2text(self.html(digits), ignore_sub)
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
