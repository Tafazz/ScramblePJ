from Bio import SeqIO, Align
import time
import tkinter as tk
import csv
from tkinter import filedialog

""""this part will ask for inqury of the sequencing list, fw primer, rw primer and barcode file and UAS file"""
def inq():
    tk.messagebox.showinfo(title="select barcode checked file", message="select barcode checked file")
    bcheck = filedialog.askopenfilename()
    
    tk.messagebox.showinfo(title="select UAS file out", message="select UAS checked file location")
    outpath = filedialog.askdirectory()

    tk.messagebox.showinfo(title=" UAS file", message="select UAS list")
    UAS = filedialog.askdirectory()

    return UAS,outpath, outname,bcheck                                           

"""here the code checked the trimed sequences for UAS. """
def UAScheck(bchecked,UAS,length,uer,ler):
    count = 0
    uchecked = list()
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    explist = list()
    reader = UAS
    row = len(reader)
    for item in bchecked:
        seq = item['seqeunce']
        ulist = list()
        lisstt = list()
        for i in range(row):
            tt = 0
            while tt == 0:
                h = []
                b = str(reader[i]['seq'])
                x = aligner.align(b,seq)
                y = aligner.align(b,seq,strand="-")
                if x[0].score > len(b)-uer:
                    rr = (x[0].aligned)[1][0][0]
                    end = (x[0].aligned)[1][-1][1]
                    h = [rr,(int(reader[i]['number']),1)]
                    ulist.append(h)
                    seq = seq[:rr] +seq[rr:end].lower()+ seq[end:]

                elif y[0].score > len(b)-uer:
                    rrr = (y[0].aligned)[1][0][0]
                    ende = (y[0].aligned)[1][-1][1]
                    h = [ende,(int(reader[i]['number']),2)]
                    ulist.append(h)
                    seq = seq[:ende] +seq[ende:rrr].lower()+ seq[rrr:]
                else:
                    tt = 1
        ulist = sorted(ulist)
        xx = len(ulist)

        """here I check if the length of my sequence is in range of +-20 bp of the 
        expected DNA seq from this combination of UAS sequences"""
        nu = nuas(seq)
        if nu != None:
            if (len(seq) > (length[0]+length[1]*xx)-ler and len(seq) < (length[0]+length[1]*xx)+ler):
                for i in ulist:
                    lisstt.append(i[1])
                item["UAS order"] = lisstt
                uchecked.append(item)
                count += 1
                """in here a list is taken from ucheck and it will find all the UAS combinations 
                , then for each combination it will count how many is recorderd for each barcode"""
                o = 0
                yy = int(item["expression"])
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

    """here I check the bin data and delete empty bins"""
    popl = list()
    for rec in explist:
        popl = []
        for key in rec[1]:
            if rec[1][key] == 0:
                popl.append(key)
        for ii in popl:
            rec[1].pop(ii)
    print ('number of U checked is %i' %count)
    return uchecked,explist

"""Here it will read the CSV files and return a list"""
def read(UAS,bcheck):
    with open(bcheck, 'r') as f:
        bchecked = list(csv.DictReader(f,delimiter=';', quotechar='"'))
    with open(UAS, 'r') as f:
        uaslist = list(csv.DictReader(f,delimiter=';', quotechar='"'))
    """length of UAS combinations from 1 to 5"""
    length = list()
    """WO barcode"""
    """With barcode"""
    length = [83,116]

    return uaslist,length,bchecked

"""checking the lengh of the seq and report the number of UAS prediction"""
def nuas(seq):
    a = readlist[1]
    emp = a[0]
    ul = a[1]
    l = len(seq)
    if (l-emp > -uer and l-emp < uer):
        return 0
    elif l-emp > uer:
        uscount = round((l-emp)/ul)
        return uscount
    else:
        return None 
    
 
"""this seqction is for loading the data:
l[0] = List of UAS names with their corresponding sequence
l[1] = Output path
l[2] = Output name handle
l[3] = list of checked sequences with their barcodes

readlist read csv files and return list
readlist[0] = UAS list
readlist[1] = Expected length for the fragement with No UAS and length of each UAS element
readlist[2] = List of barcode checked sequences
"""
l = inq()


"""this section is for program parameters:
uer = UAS anealing error
ler = final length error"""
uer = 15
ler = 20

"""a function to read csv list of UAS and barcode"""
readlist = read(l[0],l[3])

"""a function to exprt the reports of UAS order expression and sequences into CSV file"""
def out(ucheckedl,UASinbin,outname,out):
    with open('%s/%s_UASinbin.csv'%(out,outname), 'w',newline='') as csvfile:
        filewriter = csv.DictWriter(csvfile, fieldnames=['order', 'bin distribution'],delimiter=';')
        filewriter.writeheader()
        for k in UASinbin:
            filewriter.writerow({'order':k[0],'bin distribution':k[1]})

    with open('%s/%s_UASSEQchecked.csv'%(out,outname), 'w',newline='') as csvfile:
        filewriter = csv.DictWriter(csvfile, fieldnames=['seqeunce','length','barcode number','expression','UAS order'],delimiter=';')
        filewriter.writeheader()
        for i in ucheckedl:
            filewriter.writerow({'seqeunce':i['seqeunce'],'length':i['length'],'barcode number':i['barcode number'],'expression':i['expression'],'UAS order':i['UAS order']})
    csvfile.close()


"""this section is for actions:
pchecked = checks the primer bindign
bchecked = checks the barcode and trim the barcode and primer sequence
uandxchecked = it aligns the sequence with UAS elements and check how many time one order is repeated in which bin
experssion = it will find the average expression of on UAS order from all the bins it was in. so retrun order,expression,total read
"""
pchecked = list() 

xcheckedt = list()
ucheckedt = list()
bchecked = readlist[2]
uandxchecked = UAScheck(bchecked,readlist[0],readlist[1],uer,ler)
uchecked = uandxchecked[0]
xchecked = uandxchecked[1]
outname = str('3X')
out(uchecked,xchecked,outname,l[1])