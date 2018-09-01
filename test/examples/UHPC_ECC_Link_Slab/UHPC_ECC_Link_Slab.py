from calla import html
from calla.GB import flexural_capacity
from calla.JTG import crack_width

# UHPC桥面连续结构验算
# 受弯承载力
flc=flexural_capacity.fc_rect(
    fc = 150/1.45, fy = 360, fcu_k = 150, Es = 2E5, M = 46,
    b = 1570, h0 = 80-21)
flc.solve()
print(flc.text())

# 裂缝宽度
cw = crack_width.crack_width(
    b = flc.b, bf = flc.b, h = 80, h0=flc.h0, cs = 15, ftck = 8.5,
    As = int(1570/50)*113.1, Ml = 46, Ms = 33, de = 12, C3=1.15)
cw.solve()
print(cw.text(3))

# 主梁（T形截面）验算
flct=flexural_capacity.fc_T(
    fc = 17.48, fy = 360, fcu_k = 38, Es = 2E5, M = 700,
    b = 220, h0 = 600-50, bf_=1570, hf_=80, As=int(1570/100)*490.9)
flct.solve()
print(flct.text())
