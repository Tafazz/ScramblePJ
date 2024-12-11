from Bio import SeqIO, Align
import time
import tkinter as tk
import csv
import ast
import random
import copy
import numpy
import gc
import collections
from tkinter import filedialog , simpledialog
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import ast
#import pandas as pd

""""this part will ask for inqury of the sequencing list, fw primer, rw primer and barcode file and UAS file"""
def inq():
    tk.messagebox.showinfo(title="select UASchecked file", message="select UASchecked file")
    exfile = filedialog.askopenfilenames()
    exfile = list(exfile)
    tk.messagebox.showinfo(title="select UAS file out", message="select UAS checked file location")
    outpath = filedialog.askdirectory()
    outname = str('3X')
    return  exfile,outpath,outname


def binstat(uaschecked):
    explist = list()
    for item in uaschecked:
        o = 0
        yy = item[1]
        lisstt = item[0]
        for t in explist:
            o = 0
            if t[0] == lisstt:
                o += 1
                if yy in t[1].keys():
                    t[1][yy] += 1
                else:
                    t[1][yy] = 1

                break
        if o == 0:
            a = dict()
            a[yy] = 1
            explist.append([lisstt,a])
    for x in explist:
        x[1] = dict(sorted(x[1].items()))
    return explist

"""Here it will read the CSV files and return a list"""
def read(UASfile):
    UASfiled = list()
    exdict = list()
    for c in UASfile:
        with open(c, 'r') as f:
            UASfiled = UASfiled + (list(csv.DictReader(f,delimiter=';', quotechar='"')))
    for i in UASfiled:
        order = list(ast.literal_eval(i['UAS order']))
        exdict.append([order,ast.literal_eval(i['expression'])])
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
            filewriter = csv.DictWriter(csvfile, fieldnames=['order','Number of repeats','number of UAS elements','barcode'],delimiter=';')
            filewriter.writeheader()
            #wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
            for i in expression:
                filewriter.writerow({'order':i[0],'Number of repeats':i[2],'number of UAS elements':len(i[0]),'barcode':i[1]})
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
UASchecked = readlist
xchecked = binstat(UASchecked)
barcodenum = 4
for i in range(barcodenum):
    bcode = i+1
    ilist = list()
    for j in xchecked:
        if bcode in j[1].keys():
            ilist.append([j[0],i+1, j[1].get(bcode)])
    out(ilist,l[1],bcode)

"""finish time"""
t1 = time.time()
total = round(t1-t0)
print('time is %i' %total)