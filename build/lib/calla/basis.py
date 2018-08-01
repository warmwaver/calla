__all__ = [
    'open_html',
    'html2text',
    'abacus',
    ]

def html2text(html, ignore_sub=True):
    """Convert html to plain text

    >>> html2text('<div><a>hello,</a><p>abc</p></div><span>gogogo</span>')
    'hello,abc\\n\\ngogogo'
    """
    s = html
    n = 0
    while True:
        i = s.find('<',n)
        j = s.find('>',i)
        if i >= 0 and j > 0:
            element = s[i+1:j]
            if element != '':
                n = i
                prefix = '<{}>'.format(element)
                suffix = '</{}>'.format(element)
                if s.find(suffix)>0:
                    # Convertion
                    if element == 'p' or element == 'div':
                        s = s.replace(suffix,'\n')
                    elif element == 'sub':
                        if ignore_sub:
                            s = s.replace(prefix,'')
                        else:
                            s = s.replace(prefix,'_')
                    elif element == 'sup':
                        s = s.replace(prefix,'^')
                    s = s.replace(prefix,'')
                    s = s.replace(suffix,'')
                else:
                    n += 1
            else:
                n = i + 1
        else:
            break
    return s

def open_html(content,path = 'calla-result.html',font='12px times-new-roman'):
    """Save and open html file in default browser"""
    f = open(path, 'w')
    f.write('<html><head><style>body{{font:{}}}</style></head><body>'.format(font))
    f.write(content)
    f.write('</body></html>')
    f.close()
    import os
    os.startfile(path)
    
class abacus:
    """Base class for all calculators.
    Generate html format report.
    """
    
    def __init__(self):
        """
        Initialize parameters from '__inputs__'.
        parameters format: (name, (alias, unit, default_value, description, choices))
        e.g. __inputs__ = OrderedDict((
            ('Es',('E<sub>s</sub>','MPa',2.0E5,'Elastic modulus of rebar')),
            ))
        'alias' is usually in html style that can be displayed better in browser.
        'choices' is optional, valid for multi-select parameters.
        """
        if hasattr(self, '__inputs__'):
            for k in self.__inputs__:
                if not hasattr(self, k):
                    v = self.__inputs__[k]
                    if isinstance(v, tuple) and len(v)>2:
                        setattr(self,k,v[2])
                        
    def formatI(self, x):
        """Format Input parameters"""
        info = self.__inputs__[x]
        s = '{} = {} {}'.format(info[0],getattr(self,x),info[1])
        return s
    
    def formatD(self, x, digits = 2):
        """Format Derived parameters, keep given digits after dicimal point."""
        info = self.__deriveds__[x]
        s = '{1} = {2:.{0}f} {3}'.format(digits, info[0],getattr(self,x),info[1])
        return s
    
    def formatX(self,x,digits=2):
        """Format parameters，choosing diferent format automatically."""
        if x in self.__inputs__:
            return self.formatI(x)
        elif x in self.__deriveds__:
            return self.formatD(x,digits)
        else:
            return x
        
    def alias(self,x):
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
        
    def express(self, expression, digits = 2, use_aliases = True, include_self = True):
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
        result = abacus.replace_by_aliases(expression, self.__inputs__) if use_aliases else expression
        if hasattr(self, '__deriveds__'):
            result = abacus.replace_by_aliases(result, self.__deriveds__)
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
    
    def replace_by_aliases(expression, aliases):
        """
        Replace parameters in expression by their aliases.
        e.g. 'alpha*beta-gamma' -> 'α*β-γ'

        >>> abacus.replace_by_aliases('(a*3+b)/c',{'a':'A','b':'B'})
        '(A*3+B)/c'
        """
        operators = ('+','-','*','/','=','(',')','{','}','&','|',' ')
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
