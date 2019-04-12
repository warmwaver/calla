# calla
Open implementation of formulas in specifications or codes for structural engineering.

结构工程标准或规范中计算公式的开放实现。

访问项目主页 https://callapy.xyz 获取更多信息。

## 使用方法
所有计算类实例都有着一致的输入及输出方式：

```python
# 示例：矩形截面抗弯承载力计算
from calla.GB.flexural_capacity import fc_rect
f = fc_rect(M=180,b=250,h0=460,fc=14.3)
f.solve()
print(f.text())
```
如果需要html格式的输出，请使用html()函数替换text()函数。

## 开发者注意事项
### 如何开发一个新的计算类
为了保证一致的输入及输出方式，方便使用，所有计算类应派生自抽象基类abacus。abacus主要的类成员如下：

- **\_\_title\_\_**：计算类名称
- **\_\_inputs\_\_**: 输入参数字典
采用OrderedDict，每个条目的格式为

```python
__inputs__ = OrderedDict((
	(parameter, (symbol, unit, default_value, name, description, choices])),
	...
	))
```

- **\_\_deriveds\_\_**：结果参数字典
采用OrderedDict，格式与__inputs__相同。
- **\_\_toggles\_\_**：参数开关字典。
某个参数的取值对其它参数可用性有影响时，将其列入__toggles__字典可以实现简单的开关功能。格式为

```python
__toggles__ = {
	parameter:{option:(disabled_parameters),...},
	...
	}
```

- **solve**: 计算函数

```python
def solve(self):
	# Do some calculations here.
	return
```

- **html**：html格式的文本输出函数

```python
def html(self, digits=2):
	return '<p>Write some outputs here.</p>'
```

- **\_html\_**：一般情况下，推荐改写\_html()函数而不是html()函数。html()会自动调用\_html()函数生成html格式的输出

```python
def _html(self, digits=2):
	yield 'Write some outputs here.'
```

- **text**：plain text格式的文本输出函数
基类abacus实现了由html格式自动转换成纯文本格式。除非有特别需求，否则无需重写此函数。

```python
def text(self, digits=2):
	return 'Write some outputs here.'
```

### 计算类输入参数命名规则
字母符号尽可能的与规范标准中的一致，下标连写（除非与python关键字冲突，例如<i>a</i><sub>s</sub>，命名为a_s），上标'采用_代替。