from Bio import SeqIO, Align
import time
import tkinter as tk
import csv
import ast
import random
import copy
import numpy
import gc
from tkinter import filedialog , simpledialog
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import ast
#import pandas as pd

""""this part will ask for inqury of the sequencing list, fw primer, rw primer and barcode file and UAS file"""
def inq():
    tk.messagebox.showinfo(title="select expresion file", message="select expresion file")
    exfile = filedialog.askopenfilename()
    tk.messagebox.showinfo(title="select UAS file out", message="select UAS checked file location")
    outpath = filedialog.askdirectory()
    outname = str('3X')
    return  exfile,outpath,outname

"""adding the modlule for calculating the average expression level"""
def expresnum(xchecked):
    count = 0
    expdic= {}
    expff = list()
    #exlist = list()
    for i in xchecked:
        bins = list()
        unum = len(i[0])
        o = 0
        num = 0
        exp = 0
        bins = list()
        t = 0
        std = []
        for y in i[1].keys():
            if y > 0:
                exstr = y
                """numeber of reads for one order in one bin"""
                numb = i[1][y]
                exp = exp+(numb * exstr)
                num = num + numb
        """add an condition for checking the total number of reads for one order"""
        if num > 0:
            count = count + num
            expf = round(exp/num)
            for bin in i[1]:
                #expbarc = bin[0] - 1
                #exprep = bin[1]
                #uu = [exlist[expbarc]]*exprep
                #std.extend(uu)
                #expbarc = expdic[bin]
                expbarc = bin
                exprep = i[1][bin]
                uu = [expbarc]*exprep
                std.extend(uu)
            t = round(numpy.std(std),1)
            expff.append((i[0],expf,num,unum,i[1],t))
    print ('number of express checked  is %i' %count)
    return expff

"""Here it will read the CSV files and return a list"""
def read(exfile):
    exdict = list()
    with open(exfile, 'r') as f:
        exfiled = list(csv.DictReader(f,delimiter=';', quotechar='"'))
    for i in exfiled:
        key = list(ast.literal_eval(i['order']))
        exdict.append([key,ast.literal_eval(i['bin distribution'])])
    return exdict


"""start the time"""
t0 = time.time()

"""this seqction is for loading the data:
l[0] = forward primer
l[1] = reverse primer 
l[2] = the sequencing file
l[3] = the barcode file 
l[4] = the UAS list file

s = the number of sequences to analyze

readlist read csv files and return list
readlist[0] = barcode list
readlist[1] = UAS list
readlist[2] = barcode expression list
"""
l = inq()

"""this section is for program parameters:
s = numebr of sequence reads
per = primer anealing error
ber = barcode anealing error
uer = UAS anealing error
bl = barcode length
ler = final length error"""
s = 100
per = 5
ber = 5
uer = 15
bl = 10
ler = 20

"""a function to read csv list of UAS and barcode"""
readlist = read(l[0])

"""a function to exprt the reports of UAS order expression and sequences into CSV file"""
def out(expression,outpath,outname):
    with open('%s/%s_expressionreport.csv'%(outpath,outname), 'w',newline='') as csvfile:
            filewriter = csv.DictWriter(csvfile, fieldnames=['order', 'expression', 'Number of repeats','number of UAS elements','distribution','str'],delimiter=';')
            filewriter.writeheader()
            #wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
            for i in expression:
                filewriter.writerow({'order':i[0],'expression':i[1],'Number of repeats':i[2],'number of UAS elements':i[3],'distribution':i[4],'str':i[5]})
    csvfile.close()

"""this section is for actions:
pchecked = checks the primer bindign
bchecked = checks the barcode and trim the barcode and primer sequence
uandxchecked = it aligns the sequence with UAS elements and check homany time one order is repeated in which bin
experssion = it will find the average expression of on UAS order from all the bins it was in. so retrun order,expression,total read
"""
pchecked = list() 
#pcheck = SeqIO.parse("/Users/alitafazoliyazdi/Desktop/python/python scripts/git/UASnanoporeanalyze/plasmid and primers/17-6-23/pchecked.faa","fasta")

xcheckedt = list()

xchecked = readlist
experssion = expresnum(xchecked)
out(experssion,l[1],l[2])

"""finish time"""
t1 = time.time()
total = round(t1-t0)
print('time is %i' %total)