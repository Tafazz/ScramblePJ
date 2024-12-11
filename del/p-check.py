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
    tk.messagebox.showinfo(title="select the sequencing file", message="select the sequencing file")
    seq = filedialog.askopenfilename()
    #seq = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/UASnanoporeanalyze/plasmid and primers/resulpass(fast).fastq'
   # seq = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/UASnanoporeanalyze/plasmid and primers/26.10.23/3X.nonscram.fastq'
    tk.messagebox.showinfo(title="select the sequencing file", message="select the fw seq")
    fw = filedialog.askopenfilename() 
    #fw = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/UASnanoporeanalyze/plasmid and primers/AT16.fa'
    tk.messagebox.showinfo(title="select the sequencing file", message="select the rw seq")
    rw = filedialog.askopenfilename()
   # rw = '/Users/alitafazoliyazdi/Desktop/python/python scripts/git/UASnanoporeanalyze/plasmid and primers/CR201.fa'
    tk.messagebox.showinfo(title="select primer check file out", message="select primer check file location")
    out = filedialog.askdirectory()
    return  fw, rw, seq, out

"""this part will align two primers with the fragment and"""
def al(fw,rv,seq,s,per,nn,outpath):
    count = 0
    countt = 0
    pcheckl = list()
    t = 0
    q = 0
    from Bio.pairwise2 import format_alignment
    fw = SeqIO.parse(fw,"fasta")
    fw = next(fw)

    rw = SeqIO.parse(rv,"fasta")
    rw = next(rw)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    """a loop that goes through all sequences and check them and put the correct one in a list and returns the list= max size 100"""
    r = list(SeqIO.parse(seq, "fastq"))
    random.shuffle(r)
    for seq_record in r:
        """just a counter for fasta input. s is the limit if reading file"""
        count += 1
        if count > s:
            count -= 1
            break
        xn = aligner.align(fw.seq,seq_record.seq,strand="-")
        yn = aligner.align(rw.seq,seq_record.seq)
        if xn[0].score >= len(fw)-per and yn[0].score >= len(rw)-per:
            seq_record = seq_record.reverse_complement()
        x = aligner.align(fw.seq,seq_record.seq)
        y = aligner.align(rw.seq,seq_record.seq, strand="-")
        """ it will check if the fw,rw primer alignment score is great then take the fw,rw sequence an delet it and upstream down stream seq""" 
        if x[0].score >= len(fw)-per and y[0].score >= len(rw)-per:
            pcheck = copy.deepcopy(seq_record)
            "find the last position of alignemtn of fw primer"
            rr = (x[0].aligned)[1][0][0]
            "find the last position of alignemtn of rw primer"
            rrr = (y[0].aligned)[1][0][0]
            #k = x[0].query
            #kk = x[0].path
            #read = seq_record.seq[rr-10:rrr+10]
            pcheck.letter_annotations = {}
            pcheck.seq = seq_record.seq[rr-10:rrr+10]
            if len(pcheck.seq) > 50:
                #print(pcheck)
                pcheckl.append(pcheck)
                countt += 1 
    print ('Total sequence checked is %i' %count)
    print ('number of primer checked is %i' %countt)
    SeqIO.write(pcheckl, "%s/%s.faa"%(outpath,nn), "fasta")
    return pcheckl

"""start the time"""
t0 = time.time()
"""name of the export file"""
nn = str('primer checked.23.10.23')
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
s = 100000000
per = 5


"""this section is for actions:
pchecked = checks the primer bindign
bchecked = checks the barcode and trim the barcode and primer sequence
uandxchecked = it aligns the sequence with UAS elements and check homany time one order is repeated in which bin
experssion = it will find the average expression of on UAS order from all the bins it was in. so retrun order,expression,total read
"""
pchecked = al(l[0],l[1],l[2],s,per,nn,l[3])

"""finish time"""
t1 = time.time()
total = round(t1-t0)

print('time is %i' %total)