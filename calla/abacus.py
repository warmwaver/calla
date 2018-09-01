__all__ = [
    'abacus',
    'replace_by_aliases',
    'InputError',
    ]

from calla.html import *
    
def replace_by_aliases(expression, aliases):
    """
    Replace parameters in expression by their aliases.
    e.g. 'alpha*beta-gamma' -> 'α*β-γ'

    >>> replace_by_aliases('(a*3+b)/c>c',{'a':'A','b':'B','c':'C'})
    '(A*3+B)/C>C'
    """
    operators = ('+','-','*','/','=','(',')','{','}','&','|','\\','·',':',' ','<','>','^')
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
    
    def __init__(self, **inputs):
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
        # initialize parameters by inputs
        self.inputs = inputs
##        for key in inputs:
##            if hasattr(self, key):
##                setattr(self, key, inputs[key])

    @property
    def inputs(self):
        """ Get a dictionary of input parameters. """
        return { attr:(getattr(self, attr) if hasattr(self, attr) else 0) for attr in self.__inputs__} if hasattr(self, '__inputs__') else None

    @inputs.setter
    def inputs(self, values):
        """ Set value for inputs """
        for key in values:
            if hasattr(self, key):
                setattr(self, key, values[key])

    def deriveds(self):
        """ Get a dictionary of derived parameters. """
        return { attr:(getattr(self, attr) if hasattr(self, attr) else 0) for attr in self.__deriveds__} if hasattr(self, '__deriveds__') else None

    def parameters(self):
        """ Get a dictionary of all parameters. """
        paras = self.inputs
        dds = self.deriveds()
        for key in dds:
            if not key in paras:
                paras[key] = dds[key]
        return paras

    @classmethod
    def _disableds(cls,**toggles):
        """
        Get parameters disabled according to the given toggles.

        >>> a=abacus
        >>> a.__toggles__={'option1':{0:('a'),1:('b')},'option2':{True:('c'),False:('d')}}
        >>> a._disableds(option1=1,option2=False)
        ['b', 'd']
        """
        r = []
        if not hasattr(cls,'__toggles__'):
            return r
        for key in cls.__toggles__:
            v = None
            if key in toggles:
                v = toggles[key]
            elif hasattr(cls, key):
                v = getattr(cls, key)
            elif hasattr(cls, '__inputs__') and hasattr(cls.__inputs__, key)\
                 and len(cls__inputs__[key])>2:
                v = cls__inputs__[key][2]
            if v == None:
                continue
            items = cls.__toggles__[key][v]
            if isinstance(items,tuple):
                for item in items:
                    r.append(item)
            else:
                r.append(items)
        return r

    def disableds(self):
        """
        Get parameters disabled.
        self.__toggles__ format:{parameter:{option:(disabled_parameters),...},...}
        e.g. __toggles__ = {
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
                        
    def format(self, parameter, digits = 2, value=None, sep=''):
        """Format Input parameters as {name}{sep}{symbol} = {value} {unit}"""
        info = None
        if parameter in self.__inputs__:
            info = self.__inputs__[parameter]
        elif parameter in self.__deriveds__:
            info = self.__deriveds__[parameter]
        value = value or getattr(self,parameter)
        if info == None or len(info)<1:
            return '{} = {}'.format(parameter,v)
        name = info[3] if (len(info)>3 and info[3] != '') else ''
        s = '{}{}'.format(name,sep) if name != '' else ''
        if digits != None and digits >= 0:
            try:
                value = '{1:.{0}f}'.format(digits, value)
            except: # v is not decimal or numbers
                pass
        symbol = info[0]
        unit = info[1] if len(info)>2 else ''
        s += '{} = {} {}'.format(symbol,value,unit)
        return s
    
    def formatX(self, *parameters, digits=2, sep='', sep_names=', '):
        """
        Format parameters，choosing diferent format automatically.
        e.g. alpha α = value1, beta β = value2, ...
        """
        s = ''
        for parameter in parameters:
            s += self.format(parameter,digits=digits,sep=sep)
            s += sep_names
        return s[:len(s)-len(sep_names)]

    @classmethod
    def symbol(self,x):
        """get symbol of x"""
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
        """
        translate expression to value formula style.
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
        operators = ('+','-','*','/','=','(',')','{','}','&','|','\\','·',':',' ','<','>','^')
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
        """
        Subclasses should implement this method in order to output
        html format reports.
        """
        if hasattr(self, '__inputs__'):
            for attr in self.__inputs__:
                yield self.formatX(attr, digits = digits)
        if hasattr(self, '__deriveds__'):
            for attr in self.__deriveds__:
                if hasattr(self, attr):
                    yield self.formatX(attr, digits = digits)
        
    def html(self, digits = 2):
        result = ''
        gen = self._html(digits)
        for p in gen:
            result += '<p>' + p + '</p>'
        gen.close()
        return result

    def text(self, digits = 2, ignore_sub=True):
        return html2text(self.html(digits), ignore_sub)

    def none_zero_check(self, *inputs:str):
        for parameter in inputs:
            if hasattr(self, parameter) and getattr(self, parameter) == 0:
                raise InputError(self, parameter, '不能=0')

    def positive_check(self, *inputs:str):
        for parameter in inputs:
            if hasattr(self, parameter):
                value = getattr(self, parameter)
                t = type(value)
                if (t == int or t == float) and not value> 0:
                    raise InputError(self, parameter, '必须>0')

    def solve(self):
        """
        Subclass should rewrite this method to do its own calculation.
        Returns should be a dictionary with key in 'self.__deriveds__'.
        e.g. {'result1':1,'result2':2}
        """
        pass

class InputError(Exception):
    def __init__(self, calculator:abacus, parameter:str, message:str=None):
        self.calculator = calculator
        self.parameter = parameter
        self.message = '{} = {} 错误{}'.format(parameter, getattr(calculator, parameter), '' if message == None else ' ({})'.format(message))
        self.xmessage = '{} 错误{}'.format(calculator.format(parameter,digits=None), '' if message == None else ' ({})'.format(message))
        Exception.__init__(self, self.message)
    def html(self):
        return self.xmessage

def common_attrs(calculators:[abacus]):
    """ Get common attributes of calculators. """
    if calculators == None or len(calculators)<1:
        return None
    types = []
    calc0 = calculators[0]
    attrs = list(calc0.__inputs__.keys())+list(calc0.__deriveds__.keys())
    for i in range(1,len(calculators)):
        calc = calculators[i]
        t = type(calc)
        if not t in types:
            types.append(t)
            attrs_t = list(t.__inputs__.keys())+list(t.__deriveds__.keys())
            attrs = [attr for attr in attrs if attr in attrs_t]
    return attrs

def parameters_table(calculators:[abacus]):
    """ Get parameters table of calculators based on common attributes. """
    attrs = common_attrs(calculators)
    if attrs == None or len(attrs)<1:
        return None
    tbl = []
    tbl.append([calc0.symbol(attr) for attr in attrs])
    for calc in calculators:
        paras = calc.parameters()
        tbl.append([paras[attr] for attr in attrs])
    return tbl
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
