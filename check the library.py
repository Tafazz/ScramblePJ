from Bio import SeqIO, Align
import time
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
import csv
import random
import copy
import numpy
import gc
#from tkinter import filedialog , simpledialog
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import ast
#import pandas as pd

""""this part will ask for inqury of the sequencing list, fw primer, rw primer and barcode file and UAS file"""
def inq():
    root = tk.Tk()
    root.withdraw()
    seq = 0
    #tk.messagebox.showinfo(title="select primer checked file", message="select seq")
    #seq = filedialog.askopenfilename()
    seq = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/check lib/resultall.fastq'
    #tk.messagebox.showinfo(title="select barcode file", message="select barcode file")
    #bc = filedialog.askopenfilename()
    bc = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/check lib/barcode.csv'
    #tk.messagebox.showinfo(title="select barcode file out", message="select barcode file location")
    #out = filedialog.askdirectory()
    out = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/check lib'
    return  seq,bc,out

def bcheck(seq,barcode,ber,outname,outpath):
    countt = 0
    nuasc = dict()
    count = 0
    blist = dict()
    bchecked = dict()
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    reader = barcode
    row = len(reader)
    for w in range(len(reader)):
        b = int(reader[w]['num'])
        """here is the list of all barcodes that will be passed to UAScheck function"""
        blist[b]=0
        bchecked = list()
    r = list(SeqIO.parse(seq, "fastq"))
    #random.shuffle(r)
    for seq_record in r:
        sc = 0
        for i in range(row):
            sc = 0
            q = str(reader[i]['fw'])
            c = str(reader[i]['rv'])
            l = len(q)+len(c)

            x = aligner.align(q,seq_record.seq)
            y = aligner.align(c,seq_record.seq, strand="-")
            if x[0].score < len(q)-per or y[0].score < len(c)-per:
                seq_record = seq_record.reverse_complement()
                x = aligner.align(q,seq_record.seq)
                y = aligner.align(c,seq_record.seq, strand="-")
            """ it will check if the fw,rw primer alignment score is great then take the fw,rw sequence an delet it and upstream down stream seq""" 
            if x[0].score >= len(q)-per and y[0].score >= len(c)-per:
                bcheck = copy.deepcopy(seq_record)
                "find the last position of alignemtn of fw primer"
                rr = (x[0].aligned)[1][0][0]
                "find the last position of alignemtn of rw primer"
                rrr = (y[0].aligned)[1][0][0]
                #k = x[0].query
                #kk = x[0].path
                #read = seq_record.seq[rr-10:rrr+10]
                bcheck.letter_annotations = {}
                bcheck.seq = seq_record.seq[rr:rrr]
                if len(bcheck.seq) > 150:
                    #print(pcheck)
                    b = int(reader[i]['num'])
                    c = int (reader[i]['expression'])
                    bcheck.annotations["barcode number"]= b
                    bcheck.annotations["expression"] = c
                    #bcheck.seq = record.seq[fw+10:-(rv+10)]
                    #iii = bcheck.annotations["barcode number"]
                    bcheck.annotations["length"] = "%i bp"%(len(bcheck.seq))
                    #print(record.seq, bcheck.seq)
                    #print(bcheck.format("fasta"))
                    bchecked.append(bcheck)
                    blist[b] += 1
                    countt += 1 
                    break
        #print(bchecked)
    SeqIO.write(bchecked, "%s/%s.faa"%(outpath,outname), "fasta")
    print ('number of barcode checked is %i' %count)
    return nuasc,blist,bchecked
"""Here it will read the CSV files and return a list"""
def read(barcode):
    with open(barcode, 'r') as f:
        barclist = list(csv.DictReader(f,delimiter=';', quotechar='"'))
    """WO barcode"""
    #length = [216,332,448,564,680,796,912,1028,1144,1260,1376]
    """With barcode"""
    length = [286,402,518,634,750,866,982,1098,1214,1330,1446,1562,1678,1794]
    #length = [941,1057,1173,1289,1405,1521,1637,1753,1869,1985,2101]
    return barclist,length

"""checking the lengh of the seq and report the number of UAS prediction"""

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
outpath = l[2]
seq = l[0]
"""this section is for program parameters:
s = numebr of sequence reads
per = primer anealing error
ber = barcode anealing error
uer = UAS anealing error
bl = barcode length
ler = final length error"""
s = 100
per = 15
ber = 5
uer = 15
bl = 10
ler = 20

"""a function to read csv list of UAS and barcode"""
readlist = read(l[1])

def out(nuasl,blist,out,bchecked,outname):
    with open('%s/%s.csv'%(out,outname), 'w',newline='') as csvfile:
        filewriter = csv.DictWriter(csvfile, fieldnames=['seqeunce','length','barcode number','expression'],delimiter=';')
        filewriter.writeheader()
        #wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        for i in bchecked:
            filewriter.writerow({'seqeunce':i.seq,'length':i.annotations["length"],'barcode number':i.annotations["barcode number"],'expression':i.annotations["expression"]})
    #with open('%s/TOTAL LIBRARY LENGTH.csv'%out, 'w',newline='') as csvfile:
        #filewriter = csv.DictWriter(csvfile, fieldnames=['number of UAS', 'Number of reads with the same length'],delimiter=';')
        #filewriter.writeheader()
        #wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        #for i in nuasl:
            #filewriter.writerow({'number of UAS':i,'Number of reads with the same length':nuasl[i]})
            #wr.writerow(i[0])
    with open('%s/barcode list after check'%out, 'w',newline='') as csvfile:
        filewriter = csv.DictWriter(csvfile, fieldnames=['barcode number','quantity'],delimiter=';')
        filewriter.writeheader()
        for i in blist:
            filewriter.writerow({'barcode number':i,'quantity':blist[i]})

    csvfile.close()
"""this section is for actions:
pchecked = checks the primer bindign
bchecked = checks the barcode and trim the barcode and primer sequence
uandxchecked = it aligns the sequence with UAS elements and check homany time one order is repeated in which bin
experssion = it will find the average expression of on UAS order from all the bins it was in. so retrun order,expression,total read
"""
outname = str('X.bchecked')
bchecked = bcheck(seq,readlist[0],ber,outname,outpath)
out(bchecked[0],bchecked[1],l[2],bchecked[2],outname)
t1 = time.time()
total = round(t1-t0)
print('time is %i' %total)