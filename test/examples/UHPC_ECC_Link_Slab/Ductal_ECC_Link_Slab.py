from calla.basis import open_html
from calla.GB import flexural_capacity
from calla.JTG import crack_width

# UHPC桥面连续结构验算
flc=flexural_capacity.fc_rect()
flc.fc = 150/1.45
flc.fy = 360
flc.fcu_k = 150
flc.Es = 2E5
flc.M = 46
flc.b = 1570
flc.h0 = 80-21
flc.solve()
print(flc.text())

cw = crack_width.crack_width()
cw.b = cw.bf = flc.b
cw.h = 80
cw.h0=flc.h0
cw.cs = 15
cw.ftck = 8.5
cw.As = int(1570/50)*113.1
cw.Ml = 46
cw.Ms = 33
cw.de = 12
cw.C3=1.15
cw.solve()
print(cw.text(3))

result = '{}<br>{}'.format(flc.html(),cw.html(3))
f = open_html(result)

# 主梁（T形截面）验算
flct=flexural_capacity.fc_T()
flct.fc = 17.48
flct.fy = 360
flct.fcu_k = 38
flct.Es = 2E5
flct.M = 700
flct.b = 220
flct.h0 = 600-50
flct.bf_=1570
flct.hf_=80
flct.As=int(1570/100)*490.9
flct.solve()
print(flct.text())
