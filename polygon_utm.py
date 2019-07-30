#! /usr/bin/python

import csv
import os
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="hossen"
__date__ ="$01/06/2011 10:49:24 AM$"

if __name__ == "__main__":
    print "bounding polygon";

import sys
import numpy as num

#polygon1

from pyproj import Proj

zone = 54
p = Proj(proj='utm',zone=zone,ellps='WGS84')
	
point=[]
L=[[140.5, 35.0],
[140.5, 41.],
[146.5, 41],
[146.5, 35.0],
[140.5, 35.0]]
L1=num.array(L)
for i in xrange(L1.shape[0]):
	easting, northing=p(L1[i,0],L1[i,1])
        print easting, northing, L1[i,0], L1[i,1]
	point.append([easting, northing])	

ofile=open('bounding_large.csv', 'w')
wr=csv.writer(ofile, delimiter=',')
wr.writerows(point)
ofile.close()

'''
L=[[140.85, 37.05],
[140.93, 37.48],
[141.15 ,37.48],
[141.07, 37.05],
[140.85, 37.05]]

L1=num.array(L)
point=[]
for i in xrange(L1.shape[0]):
	zone, easting, northing=LLtoUTM(L1[i,1],L1[i,0])
	point.append([easting, northing])	
ofile=open('inner2.csv', 'w')
wr=csv.writer(ofile, delimiter=',')
wr.writerows(point)
ofile.close()

L=[[140.93, 37.1],
[140.97, 37.43],
[141.12, 37.43],
[141.045, 37.1],
[140.93, 37.1]]

L1=num.array(L)
point=[]
for i in xrange(L1.shape[0]):
	zone, easting, northing=LLtoUTM(L1[i,1],L1[i,0])
	point.append([easting, northing])	
ofile=open('inner3.csv', 'w')
wr=csv.writer(ofile, delimiter=',')
wr.writerows(point)
ofile.close()
'''
