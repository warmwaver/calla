# calla
Open implementation of formulas in specifications or codes for structural engineering.

结构工程标准或规范中计算公式的开放实现。

访问项目主页 https://callapy.org/ 获取更多信息.

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
	yield 'Write some outputs here.'
```

- **text**：plain text格式的文本输出函数
基类abacus实现了由html格式自动转换成纯文本格式。除非有特别需求，否则无需重写此函数。

### 计算类输入参数命名规则
字母符号尽可能的与规范标准中的一致，下标连写（除非与python关键字冲突，例如<i>a</i><sub>s</sub>，命名为a_s），上标'采用_代替。