#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 10:50:50 2019

@author: peterwu
"""

import numpy as np
import timeit


start = timeit.default_timer()

times=100000
n=1000

x=[None]*times
for i in range(times):
    x[i]=np.random.standard_normal(n)
    
    
stop = timeit.default_timer()

run_time_string=str(stop - start)   

print(run_time_string)

f = open( 'performace_test.txt', 'w' )
f.write( run_time_string )
f.close()
