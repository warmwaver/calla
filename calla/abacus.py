__all__ = [
    'abacus',
    'replace_by_aliases',
    'InputError',
    'SolvingError',
    'InputWarning'
    ]

from calla.html import html2text
from collections import OrderedDict
    
def replace_by_aliases(expression, aliases):
    """
    Replace parameters in expression by their aliases.
    e.g. 'alpha*beta-gamma' -> 'α*β-γ'

    >>> replace_by_aliases('(a*3+b)/c>c',{'a':'A','b':'B','c':'C'})
    '(A*3+B)/C>C'
    """
    operators = (
            '+','-','*','/','=','(',')','{','}','&','|','\\',
            '·',':',' ','<','>','^')
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
        Recommended parameters format: (parameter, symbol, unit, default_value, name, description[, choices])
        e.g. __inputs__ = [
            ('Es','E<sub>s</sub>','MPa',2.0E5,'Elastic modulus of rebar'),
            ]
        Supported parameters format (obselete): (parameter, (symbol, unit, default_value, name, description[, choices]))
        e.g. __inputs__ = OrderedDict((
            ('Es',('E<sub>s</sub>','MPa',2.0E5,'Elastic modulus of rebar')),
            ))
        'alias' is usually in html style that can be displayed better in browser.
        'choices' is optional, valid for multi-select parameters.
        """
        self._inputs_ = OrderedDict()
        if hasattr(self, '__inputs__'):
            if isinstance(self.__inputs__, list):
                for item in self.__inputs__:
                    if isinstance(item, tuple) and len(item) > 0:
                        k = item[0]
                        v = item[3] if len(item) > 3 else 0
                        setattr(self, k, v)
                        # create parameter object with dynamic type
                        self._inputs_[k] = self.create_param_object(item)
            else:
                for k in self.__inputs__:
                    infos = self.__inputs__[k]
                    if isinstance(infos, tuple):
                        self._inputs_[k] = self.para_attrs(k)
                        if not hasattr(self, k):
                            v = infos[2] if len(infos) > 2 else 0
                            setattr(self, k, v)

        # initialize toggles
        if hasattr(self, '__toggles__'):
            if isinstance(self.__toggles__, list):
                self._toggles_ = OrderedDict()
                count = int(len(self.__toggles__)/2)
                for i in range(count):
                    k = self.__toggles__[2*i]
                    v = self.__toggles__[2*i+1]
                    if isinstance(k, str) and isinstance(v, dict):
                        self._toggles_[k] = v
            else:
                self._toggles_ = self.__toggles__

        # initialize parameters by inputs
        self.inputs = inputs

    def _init_deriveds_(self):
        if not hasattr(self, '_deriveds_'):
            self._deriveds_ = OrderedDict()
            if hasattr(self, '__deriveds__'):
                if isinstance(self.__deriveds__, list):
                    for item in self.__deriveds__:
                        if isinstance(item, tuple) and len(item) > 1:
                            self._deriveds_[item[0]] = self.create_param_object(item)
                else:
                    for key in self.__deriveds__:
                        self._deriveds_[key] = self.para_attrs(key)

    @property
    def order_of_inputs(self):
        order = []
        if hasattr(self, '__inputs__'):
            if isinstance(self.__inputs__, list):
                for item in self.__inputs__:
                    if isinstance(item, tuple) and len(item) > 0:
                        order.append(item[0])
            else:
                order = list(self.__inputs__.keys())
        return order

    @property
    def inputs(self):
        """ Get a dictionary of input parameters. """
        return { attr:(getattr(self, attr) if hasattr(self, attr) else 0) for attr in self._inputs_} if hasattr(self, '__inputs__') else None

    @inputs.setter
    def inputs(self, values):
        """ Set value for inputs """
        def _setvalue(cls, key, value):
            t = type(getattr(cls,key))
            if t is bool and not isinstance(value, bool):
                value = True if str(value) == 'True' else False 
            setattr(cls, key, value)

        # set values for all inputs
        for key in values:
            if hasattr(self, key):
                _setvalue(self, key, values[key])

        # toggles may reset the value of some disabled inputs
        if hasattr(self, '__toggles__'):
            for key in self._toggles_:
                if key in values:
                    _setvalue(self, key, values[key])

    def deriveds(self):
        """ Get a dictionary of derived parameters. """
        return { attr:(getattr(self, attr) if hasattr(self, attr) else 0) for attr in self._deriveds_} if hasattr(self, '__deriveds__') else None

    @property
    def parameters(self):
        """ Get a dictionary of all parameters. """
        paras = self.inputs
        dds = self.deriveds()
        for key in dds:
            if not key in paras:
                paras[key] = dds[key]
        return paras

    @classmethod
    def _disableds(cls, **toggles):
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
        
        if isinstance(cls.__toggles__, list):
            count = int(len(cls.__toggles__)/2)
            for i in range(count):
                key = cls.__toggles__[2*i]
                choices = cls.__toggles__[2*i+1]
                # if isinstance(key, str) and isinstance(v, dict):
                #     self._toggles_[key] = v
                # A latter toggle can be diabled by a previous toggle,
                # together with its' toggle function.
                if key in r:
                    continue
                v = None
                if key in toggles:
                    v = toggles[key]
                # elif hasattr(cls, key):
                #     v = getattr(cls, key)
                # elif hasattr(cls, '__inputs__') and hasattr(cls.__inputs__, key)\
                #     and len(cls.__inputs__[key])>2:
                #     v = cls.__inputs__[key][2]
                if v == None:
                    continue
                
                if v not in choices:
                    continue
                items = choices[v]
                if isinstance(items,tuple):
                    for item in items:
                        r.append(item)
                else:
                    r.append(items)
        else:
            # old format calculator
            for key in cls.__toggles__:
                # A latter toggle can be diabled by a previous toggle,
                # together with its' toggle function.
                if key in r:
                    continue
                v = None
                if key in toggles:
                    v = toggles[key]
                elif hasattr(cls, key):
                    v = getattr(cls, key)
                elif hasattr(cls, '__inputs__') and hasattr(cls.__inputs__, key)\
                    and len(cls.__inputs__[key])>2:
                    v = cls.__inputs__[key][2]
                if v == None:
                    continue
                choices = cls.__toggles__[key]
                if v not in choices:
                    continue
                items = choices[v]
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
        if not hasattr(self,'_toggles_'):
            return r
        for key in self._toggles_:
            # A latter toggle can be diabled by a previous toggle,
            # together with its' toggle function.
            if key in r: 
                continue
            choices = self._toggles_[key]
            if hasattr(self, key):
                v = getattr(self, key)
                if v not in choices:
                    continue
                items = choices[v]
                r.extend(items)
        return r
                        
    def format(self, parameter, digits = 2, value=None, sep=' ', omit_name=False, eq=None):
        """
        Format parameter as {name}{sep}{symbol} = {value} {unit}
        e.g. force N = 100 kN
        
        Arguments:
            sep: seperator between parameter's name and symbol
            sep_names: seperator between parameters' names
            omit_name: format with parameter name or not
            eq: equation (or expression) of parameter
        """
        self._init_deriveds_()
        if hasattr(self, '_inputs_'):
            info = None
            if parameter in self._inputs_:
                info = self._inputs_[parameter]
            elif parameter in self._deriveds_:
                info = self._deriveds_[parameter]
            value = value if value != None else getattr(self,parameter)
            if info == None:
                value = '{1:.{0}f}'.format(digits, value) if digits != None and digits >= 0 else value
                return '{} = {}'.format(parameter, value)
            # use choices' value to substitude parameter value
            if info.choices:
                if isinstance(info.choices, dict) and value in info.choices:
                    value = info.choices[value]
            s = ''
            if not omit_name:
                s = '{}{}'.format(info.name,sep) if info.name != '' else ''
            if digits != None and digits >= 0:
                try:
                    vabs = abs(value)
                    value = '{1:.{0}{2}}'.format(digits, value, 'e' if (vabs>1e4 or (vabs>0 and vabs<10**-digits)) else 'f')
                except: # v is not decimal or numbers
                    pass
            if info.symbol:
                s += '{} = '.format(info.symbol)
            if not eq:
                s += '{} {}'.format(value, info.unit)
            else:
                s += '{} = {} {}'.format(
                        self.replace_by_symbols(eq), value, info.unit)
            return s
        info = None
        if parameter in self.__inputs__:
            info = self.__inputs__[parameter]
        elif parameter in self.__deriveds__:
            info = self.__deriveds__[parameter]
        value = value if value != None else getattr(self,parameter)
        if info == None or len(info)<1:
            value = '{1:.{0}f}'.format(digits, value) if digits != None and digits >= 0 else value
            return '{} = {}'.format(parameter, value)
        # use choices' value to substitude parameter value
        if len(info)>5:
            choices = info[5]
            if type(choices) is dict and value in choices:
                value = choices[value]
        s = ''
        if not omit_name:
            name = info[3] if (len(info)>3 and info[3] != '') else ''
            s = '{}{}'.format(name,sep) if name != '' else ''
        if digits != None and digits >= 0:
            try:
                vabs = abs(value)
                value = '{1:.{0}{2}}'.format(digits, value, 'e' if (vabs>1e4 or (vabs>0 and vabs<10**-digits)) else 'f')
            except: # v is not decimal or numbers
                pass
        symbol = info[0]
        if symbol != '' and symbol != None:
            s += '{} = '.format(symbol)
        unit = info[1] if len(info)>2 else ''
        if eq == None or eq == '':
            s += '{} {}'.format(value, unit)
        else:
            s += '{} = {} {}'.format(
                    self.replace_by_symbols(eq), value, unit)
        return s
    
    def formatx(self, *parameters, digits=2, sep='', sep_names=', ', omit_name=True, toggled=True):
        """
        Format parameters，choosing diferent format automatically.
        e.g. force N = 100 kN, moment M=100 kN·m, ...
        
        Arguments:
            sep: seperator between parameter's name and symbol
            sep_names: seperator between parameters' names
            omit_name: format with parameter name or not
            toggled: if True, visibility is decided by toggles
        """
        s = ''
        disableds = self.disableds()
        for parameter in parameters:
            if not hasattr(self, parameter):
                continue
            if toggled and parameter in disableds:
                continue
            s += self.format(parameter,digits=digits,sep=sep,omit_name=omit_name)
            s += sep_names
        return s[:len(s)-len(sep_names)]

    def format_conclusion(self, ok: bool, eq_left, comparison_symbol, eq_right, message: str=None):
        if not message:
            message = '{}满足规范要求。'.format('' if ok else '不')
        return '<span class="conclusion_{}">{} {} {}，{}</span>'.format(
            'ok' if ok else 'ng', eq_left, comparison_symbol, eq_right, message
            )

    @staticmethod
    def create_param_object(attrs: tuple):
        """create parameter object with dynamic type"""
        if isinstance(attrs, tuple):
            n = len(attrs)
            _choices = None
            if n > 6:
                _choices_origin = attrs[6]
                if isinstance(_choices_origin, list):
                    _choices = OrderedDict()
                    for item in _choices_origin:
                        if isinstance(item, tuple):
                            _choices[item[0]] = item[1]
                    if len(_choices) < 1: # if isinstance(_choices_origin, dict):
                        _choices = _choices_origin
                else: # if isinstance(_choices_origin, dict):
                    _choices = _choices_origin

            para_obj = type(
                'param', (object,),
                dict(
                    symbol=attrs[1] if n > 1 else '',
                    unit=attrs[2] if n > 2 else '',
                    default_value=attrs[3] if n > 3 else None,
                    name=attrs[4] if n > 4 else '',
                    description=attrs[5] if n > 5 else '',
                    choices=_choices
                ))
            return para_obj()
        return None

    # TODO: 不适用于新的初始化格式，需改造
    @classmethod
    def para_attrs(self, parameter):
        """get attributes of parameter"""
        if parameter in self.__inputs__:
            attrs = self.__inputs__[parameter]
        elif parameter in self.__deriveds__:
            attrs = self.__deriveds__[parameter]
        else:
            raise Exception('parameter not exists.')
        n = len(attrs)
        cls_para = type(
                'para', (object,),
                dict(symbol=attrs[0] if n>0 else '',
                 unit=attrs[1] if n>1 else '',
                 default_value=attrs[2] if n>2 else None,
                 name=attrs[3] if n>3 else '',
                 description=attrs[4] if n>4 else '',
                 choices=attrs[5] if n>5 else None
                 ))
        return cls_para()
        
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
        if not hasattr(self, '_params_'):
            self._params_ = {}
            if hasattr(self, '_inputs_'):
                for key in self._inputs_:
                    self._params_[key] = self._inputs_[key].symbol
            self._init_deriveds_()
            if hasattr(self, '_deriveds_'):
                for key in self._deriveds_:
                    self._params_[key] = self._deriveds_[key].symbol      
        s = replace_by_aliases(expression, self._params_)
        return s
        # if hasattr(self, '__inputs__'):
        #     params = self.__inputs__.copy()
        #     if hasattr(self, '__deriveds__'):
        #         params.update(self.__deriveds__)
        #     s = replace_by_aliases(expression, params)
        #     return s
        # return expression
    
    def _html(self, digits = 2):
        """
        Subclasses should implement this method in order to output
        html format reports.
        """
        disableds = self.disableds()
        if hasattr(self, '_inputs_'):
            for attr in self._inputs_:
                if hasattr(self, attr) and (not attr in disableds):
                    yield self.format(attr, digits = None)
        if hasattr(self, '_deriveds_'):
            for attr in self._deriveds_:
                if hasattr(self, attr) and (not attr in disableds):
                    yield self.format(attr, digits = digits)
        
    def html(self, digits = 2):
        result = ''
        gen = self._html(digits)
        for p in gen:
            result += '<p>' + p + '</p>'
        gen.close()
        return result

    def text(self, digits = 2, sub='', sup=''):
        return html2text(self.html(digits), sub, sup)

    def validate(self, requirement:str, *inputs:str):
        for parameter in inputs:
            if hasattr(self, parameter):
                value = getattr(self, parameter)
                t = type(value)
                if (t == int or t == float):
                    if requirement == 'positive':
                        if value <= 0:
                            raise InputError(self, parameter, '应>0')
                    elif requirement == 'non-negative':
                        if value < 0:
                            raise InputError(self, parameter, '应≥0')
                    elif requirement == 'non-zero':
                        if value == 0:
                            raise InputError(self, parameter, '应≠0')

    def none_zero_check(self, *inputs:str):
        # obselete
        for parameter in inputs:
            if hasattr(self, parameter) and getattr(self, parameter) == 0:
                raise InputError(self, parameter, '不能=0')

    def positive_check(self, *inputs:str):
        # obselete
        for parameter in inputs:
            if hasattr(self, parameter):
                value = getattr(self, parameter)
                t = type(value)
                if (t == int or t == float) and not value> 0:
                    raise InputError(self, parameter, '应>0')

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
        self.message = '输入错误：{} = {} {}'.format(parameter, getattr(calculator, parameter), '' if message == None else ' ({})'.format(message))
        self.xmessage = '输入错误：{} {}'.format(calculator.format(parameter,digits=None), '' if message == None else ' ({})'.format(message))
        Exception.__init__(self, self.message)
    def html(self):
        return self.xmessage

class SolvingError(Exception):
    def __init__(self, message:str):
        self.message = message
    def html(self):
        return self.message

class InputWarning(Warning):
    pass

def common_attrs(calculators:list[abacus]):
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

def parameters_table(calculators:list[abacus]):
    """ Get parameters table of calculators based on common attributes. """
    attrs = common_attrs(calculators)
    if attrs == None or len(attrs)<1:
        return None
    tbl = []
    calc0 = calculators[0]
    tbl.append([calc0.symbol(attr) for attr in attrs])
    for calc in calculators:
        paras = calc.parameters()
        tbl.append([paras[attr] for attr in attrs])
    return tbl
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
