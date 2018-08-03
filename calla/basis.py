__all__ = [
    'abacus',
    'html2text',
    'table2html',
    'table2text',
    'open_html',
    'default_html_style',
    ]

def html2text(html, ignore_sub=True):
    """Convert html to plain text

    >>> html2text('<div><a>hello,</a><p>world</p></div><span>i miss you</span>')
    'hello,world\\n\\ni miss you'
    """
    def split_block(html):
        """
        Returns:
        element,front,prefix,content,suffix,back
        """
        start = 0
        i = html.find('<',start)
        j = html.find('>',i)
        if i >= 0 and j > 0:
            content = html[i+1:j]
            element = content.split()[0]
            if element != '':
                start = i
                prefix = '<{}>'.format(content)
                suffix = '</{}>'.format(element)
                end = html.find(suffix,start+len(prefix))
                # locate the right position of suffix
                pre = '<{}'.format(element)
                n = end
                while True:
                    m = html.rfind(pre, start+len(prefix), n)
                    if m>0 and html.find('>',m+len(pre),n) > 0:
                        end = html.find(suffix,end+len(suffix))
                        n = m
                    else:
                        break
                if end>0:
                    return (element,html[:i],prefix,html[start+len(prefix):end],suffix,html[end+len(suffix):])
        return ('','','',html,'','')
    
    s = ''
    start = 0
    while html != '':
        element,front,prefix,content,suffix,back = split_block(html)
        if element != '' and prefix != '' and suffix != '':
            prefix = ''
            suffix = ''
            if element == 'p' or element == 'div'\
                or element == 'table' or element == 'tr':
                suffix = '\n'
            elif element == 'sub':
                prefix = '' if ignore_sub else '_'
            elif element == 'sup':
                prefix = '^'
            elif element == 'head' or  element == 'style':
                content = ''
            elif element == 'td' or element == 'th':
                suffix = '\t'
            elif element.startswith('h') and element[1:].isdigit():
                suffix = '\n'
            s += front + prefix + html2text(content) + suffix
            html = back
        else:
            s += content
            break
    return s

def table2html(table, digits=2):
    def f(table, digits=2):
        yield '<table border="1" cellpadding="5">'
        i = 0
        for row in table:
            element = 'td' if i>0 else 'th'
            yield '<tr>'
            j = 0
            for col in row:
                align = 'right'  if j>0 else 'left'
                if isinstance(col,float):
                    yield '<{1} align="{3}">{2:.{0}f}</{1}>'.format(digits, element, col, align)
                else:
                    yield '<{0} align="{2}">{1}</{0}>'.format(element, col, align)
                j += 1
            yield '</tr>'
            i += 1
        yield '</table>'
    # generate html
    result = ''
    gen = f(table,digits)
    for p in gen:
        result += p
    gen.close()
    return result
    
def table2text(table, digits=2):
    def f(table, digits=2):
        for row in table:
            for col in row:
                if isinstance(col,float):
                    yield '{1:.{0}f}\t'.format(digits, col)
                else:
                    yield '{0}\t'.format(col)
            yield '\n'
    # generate text
    result = ''
    gen = f(table,digits)
    for p in gen:
        result += p
    gen.close()
    return result

default_html_style = '''
body{font:12px times-new-roman}
table{border-collapse:collapse; font-size:12px;}
'''
            
def open_html(content,path = 'calla-result.html',style=default_html_style):
    """Save and open html file in default browser"""
    f = open(path, 'w')
    f.write('<html><head><style>{}</style></head><body>'.format(style))
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
