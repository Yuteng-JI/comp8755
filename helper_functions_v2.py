#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math
from appscript import k
import mmh3
from bitarray import bitarray
import copy 
import random
import matplotlib.pyplot as plt
import csv
import numpy as np


# In[3]:


class WeightedBloomFilter(object):
    def __init__(self, m, n, rs, s2) -> None:
        self.m = m  # m is the bit-array size 
        self.n = n  # n is the number of input elements
        self.rs = rs  # expected re of all the elements 
        self.s2 = s2  # One parameter in hash number formula
        self.bit_array = bitarray(self.m)
        self.bit_array.setall(0)

    def add(self, item, re_index):
        # get the number of hash functions 
        hash_count = self.get_hash(re_index)
        for i in range(hash_count):
            digest = mmh3.hash(item, i) % self.m
            self.bit_array[digest] = True 

    def check(self, item, re_index):
        # check the element
        hash_count = self.get_hash(re_index)
        for i in range(hash_count):
            digest = mmh3.hash(item, i) % self.m
            if self.bit_array[digest] == False:
                return False
        return True

    def get_hash(self, re_index):
        # get the number of hash functions 
        re = self.rs[re_index]
        a = self.m*math.log(2)/self.n
        b = math.log(re) - self.s2/self.n 
        return round(a + b/math.log(2))


# In[ ]:


'generate a bloom filter of one person SNP'
def generate_bloom(in_path, p, size):
    in_file = open(in_path, 'r')
    lines = in_file.readlines()
    n = len(lines)
    BF = BloomFilter(n,p,size)
    for line in lines:
        BF.add(line)
    return BF

'Write the bloom filter into txt file'
def write_bloom(BF, out_path):
    out_file = open(out_path, "w")
    write_str = ''
    for bit in BF.bit_array:
        if bit==True:
            write_str += '1'
        else:
            write_str += '0'
    out_file.write(write_str)


'Cardinality Estimation of two bloom filter intersection'
def card_inter(bf1, bf2):
    c_1 = bf1.card_check()
    c_2 = bf2.card_check()
    bf_union = copy.copy(bf1)
    bf_union = bf_union.union(bf2)
    return c_1 + c_2 - bf_union.card_check()

'Cardinality Estimation of three bloom filter intersection'
def card_inter3(bf1, bf2, bf3):
    c_1 = bf1.card_check()
    c_2 = bf2.card_check()
    c_3 = bf3.card_check()
    c_12 = card_inter(bf1, bf2)
    c_13 = card_inter(bf1, bf3)
    c_23 = card_inter(bf2, bf3)
    bf_union = copy.deepcopy(bf1)
    bf_union = bf_union.union(bf2).union(bf3)
    c_123 = bf_union.card_check()
    return c_123 - c_1 - c_2 - c_3 + c_12 + c_13 + c_23

'Test child has how many of homozygote of fathers, output |(a1 & a2) & (b1 | b2)|'
def inherit_test(bfa1, bfa2, bfb1, bfb2):
    b_union = copy.deepcopy(bfb1)
    b_union = b_union.union(bfb2)
    return card_inter3(bfa1, bfa2, b_union)

    

'Compute actural intersection of two SNP data'
def actual_insec(file1, file2):
    '''
    The function computes the actual intersection of two snp data set
    Inputs are SNP variates of person1 and person2 
    Return the number of same SNPs
    and the same positions
    '''
    file1 = open(file1)
    file2 = open(file2)
    l1 = file1.readlines()
    l2 = file2.readlines()
    size1 = len(l1)
    size2 = len(l2)
    i = 1
    j = 1 
    X = 0
    res = []
    while i < size1 and j < size2:
        ls1 = l1[i].split('/')
        ls2 = l2[j].split('/')
        pos1 = int(ls1[1])
        pos2 = int(ls2[1])
        if pos1 == pos2:
            X += 1
            res.append(pos1)
            i += 1 
            j += 1 
        elif pos1 < pos2:
            i += 1
        else:
            j += 1 
    return X, res

