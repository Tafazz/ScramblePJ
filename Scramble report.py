import tkinter as tk
import csv
import ast
from tkinter import filedialog

""""this part will ask for inqury"""
def inq():
    tk.messagebox.showinfo(title="select UAS checked file", message="select UAS checked file")
    exfile = filedialog.askopenfilenames()
    exfile = list(exfile)

    tk.messagebox.showinfo(title="select UAS file out", message="select UAS checked file location")
    outpath = filedialog.askdirectory()

    return  exfile,outpath

"""Here it will count the number of time each UAS order is repeated"""
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
        exdict.append([order,ast.literal_eval(i['barcode number'])])
    return exdict


"""this seqction is for loading the data:
l[0] = forward primer
l[1] = reverse primer """
l = inq()

"""readlist read csv files and return a list of each seuquence UAS order and its barcode"""
readlist = read(l[0])

"""a function to exprt the reports of UAS order and sequences into CSV file"""
def out(Crep,outpath,outname):
    with open('%s/barcode%s_scramblereport.csv'%(outpath,outname), 'w',newline='') as csvfile:
            filewriter = csv.DictWriter(csvfile, fieldnames=['order','Number of repeats','number of UAS elements','barcode'],delimiter=';')
            filewriter.writeheader()
            for i in Crep:
                filewriter.writerow({'order':i[0],'Number of repeats':i[2],'number of UAS elements':len(i[0]),'barcode':i[1]})
    csvfile.close()


"""this section is for program parameters:
barcodenum = number of different barcodes for the report
"""
barcodenum = 4


"""this section is for actions:
xchecked = returns the list of each UAS order and the number of times it is repeated for each barcode
"""
xchecked = binstat(readlist)
for i in range(barcodenum):
    bcode = i+1
    ilist = list()
    for j in xchecked:
        if bcode in j[1].keys():
            ilist.append([j[0],i+1, j[1].get(bcode)])
    out(ilist,l[1],bcode)