#!/usr/bin/python\
import math
print('supply arguments (oldhight, old radius)')
def radhight(hight, rad):
    new_rad= hight/2
    new_hight= rad*2
    change_pos = -.5*(hight- new_hight)
    ans = 'new hight = {0:.5f}\nnew radius= {1:.5f}'.format(new_hight,new_rad)
    ans2 = 'change pos by {0:.5f}'.format(change_pos) 
    print (ans)
    print (ans2) 
