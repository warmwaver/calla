__all__ = [
    'html2text',
    'table2html',
    'table2text',
    'save',
    'save_and_open',
    'default_html_style',
    ]

default_html_style = '''
body{font:16px times new roman,宋体}
table{border-collapse:collapse; font-size:14px;}
'''

def html2text(html, sub='', sup=''):
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
    while html != '':
        element,front,prefix,content,suffix,back = split_block(html)
        if element != '' and prefix != '' and suffix != '':
            prefix = ''
            suffix = ''
            if element == 'p' or element == 'div'\
                or element == 'table' or element == 'tr':
                suffix = '\n'
            elif element == 'sub':
                prefix = sub
            elif element == 'sup':
                prefix = sup
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
    # substitude special symbols
    specials = {'&lt;':'<', '&le;':'<=', '&gt;':'>', '&ge;':'>='}
    for key in specials:
        s = s.replace(key, specials[key])
    return s

def table2html(table, digits=2, reverse=False, attr_class='result'):
    def f(table, digits=2):
        yield '<table class="{}" border="1" cellpadding="5">'.format(attr_class)
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
    def rf(table, digits=2):
        yield '<table class="{}" border="1" cellpadding="5">'.format(attr_class)
        ncol = 0
        for i in range(len(table)):
            ncol = max(ncol, len(table[i]))
        for j in range(ncol):
            yield '<tr>'
            for i in range(len(table)):
                element = 'td' if i>0 else 'th'
                align = 'right'
                col = table[i][j]
                if isinstance(col,float):
                    yield '<{1} align="{3}">{2:.{0}f}</{1}>'.format(digits, element, col, align)
                else:
                    yield '<{0} align="{2}">{1}</{0}>'.format(element, col, align)
            yield '</tr>'
        yield '</table>'
    # generate html
    result = ''
    gen = rf(table,digits) if reverse else f(table,digits)
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

def save(content,path = 'calla-result.html',style=default_html_style):
    """Save and open html file in default browser"""
    f = open(path, 'w')
    f.write('<html><head><style>{}</style></head><body>'.format(style))
    f.write(content)
    f.write('</body></html>')
    f.close()

def save_and_open(content,path = 'calla-result',style=default_html_style):
    if not (path.endswith('.html') or path.endswith('htm')):
        path += '.html'
    save(content, path,style)
    import os
    if os.path.isfile(path):
        os.startfile(path)
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
