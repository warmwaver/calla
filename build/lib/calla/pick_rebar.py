"""
Given total area of rebars, pick the available combinations. 
"""

import copy
import math
def rebar_area(diameter):
    return math.pi/4*diameter**2
class pick_rebar:
    As = 1200
    # number of rebars
    nrebars = 7
    # rebar diameter that can be used
    dia_choice = [18,20,22,25,28] #mm
    # max number of rebar diameters to be combined
    nmax = 2
    combination = []
    combinations = []
    # total area limitation (by percent)
    limit = 1000
    # pick the first combination which meet the requirement
    pickfirst = False
    out = ''
    def pick(self):
        self.combinations = []
        self.out = ''
        m = self.dia_choice.__len__()
        n = 1 #number of combination
        seprater = '\n'
        while n<=self.nmax:
            combination = list(0 for x in range(n))
            i = n-1
            # 循环枚举m个元素中的n个值的组合
            while True:
                total = 0
                for choice in combination:
                    total += self.nrebars*rebar_area(self.dia_choice[choice])
                # save result
                if total>=self.As and total<=self.As*(1+self.limit):
                    result = ''
                    for j in combination:
                        result += '{}*{}+'.format(self.nrebars,self.dia_choice[j])
                    if self.pickfirst:
                        self.combination = combination
                        self.out = result[:-1]
                        return self.out
                    self.combinations.append(copy.copy(combination))
                    self.out += result[:-1] + seprater
                # 判断取值是否已到最后
                if combination[i]<m-1:
                    combination[i]+=1
                else:
                    while i>0:
                        i-=1
                        if combination[i]<m-(n-i)+1:
                            combination[i]+=1
                            while i<n-1:
                                i+=1
                                combination[i]=combination[i-1]
                            break
                if (n==1 and combination[0]==m-1) or (combination[0]==m-1 and combination[n-1]==m-1):
                    total = 0
                    for choice in combination:
                        total += self.nrebars*rebar_area(self.dia_choice[choice])
                    if total>=self.As and total<=self.As*(1+self.limit):
                        result = ''
                        for j in combination:
                            result += '{}*{}+'.format(self.nrebars,self.dia_choice[j])
                        if self.pickfirst:
                            self.combination = combination
                            self.out = result[:-1]
                            return self.out
                        self.combinations.append(copy.copy(combination))
                        self.out += result[:-1] + seprater
                    break
            n+=1
        self.out = self.out[:-1]
        # choose the combination of minimum area
        if self.combinations.__len__()>0:
            min_cb = self.combinations[0]
            min_area = 0
            for i in self.combinations[0]:
                min_area += rebar_area(self.dia_choice[i])
            for cb in self.combinations:
                area = 0
                for i in cb:
                    area += rebar_area(self.dia_choice[i])
                if area<min_area:
                    min_area = area
                    min_cb = cb
            self.combination = min_cb
            result = ''
            for i in min_cb:
                result += '{}*{}+'.format(self.nrebars,self.dia_choice[i])
            result = result[:-1]
            return result
    # chosen rebar combination's total area
    def area(self, combination=None):
        a = 0
        if combination==None:
            combination = self.combination
        for i in combination:
            a += rebar_area(self.dia_choice[i])
        return self.nrebars*a

def _test1_():
    pr = pick_rebar()
    pr.As = 1110
    pr.limit = 5
    pr.dia_choice = [18,20,22,25,28,32]
    print('---rebar picked---')
    print(pr.pick())
    print('---available combinations---')
    print(pr.out)

def _test2_():
    pr = pick_rebar()
    pr.As = 10744
    pr.nrebars = 26
    pr.nmax = 1
    pr.dia_choice = [18,20,22,25]
    print('---rebar picked---')
    print(pr.pick())
    print('---available combinations---')
    print(pr.out)

if __name__ == '__main__':
    _test1_()
    _test2_()
