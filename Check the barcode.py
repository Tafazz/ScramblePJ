from Bio import SeqIO, Align
import tkinter as tk
from tkinter import filedialog
import csv
import copy

""""this part will ask for inqury of the sequencing list, fw primer, rw primer and barcode file and UAS file"""
def inq():
    root = tk.Tk()
    root.withdraw()

    tk.messagebox.showinfo(title="Sequencing file", message="select the sequencing file")
    seq = filedialog.askopenfilename()
    
    tk.messagebox.showinfo(title="select barcode file", message="select barcode file")
    bc = filedialog.askopenfilename()

    tk.messagebox.showinfo(title="select output", message="select output location")
    out = filedialog.askdirectory()

    return  seq,bc,out

"""here the code checks the sequences for the barcode and put each sequence in its related dict and trim the primer and barcode """                                          
def bcheck(seq,barcode,outname,outpath):
    nuasc = dict()
    count0= 0
    countf= 0
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

    for seq_record in r:
        check = 0
        count0+= 1
        for i in range(row):
            q = str(reader[i]['fw'])
            c = str(reader[i]['rv'])
            x = aligner.align(q,seq_record.seq)
            y = aligner.align(c,seq_record.seq)

            if x[0].score < len(q)-(len(q)/10):
                seq_record = seq_record.reverse_complement()
                x = aligner.align(q,seq_record.seq)
                y = aligner.align(c,seq_record.seq)
                #y = aligner.align(c,seq_record.seq, strand="-")
            """ it will check if the fw,rw primer alignment score is great then take the fw,rw sequence an delet it and upstream down stream seq
            https://biopython.org/docs/1.76/api/Bio.Align.html?highlight=aligner#Bio.Align.PairwiseAligner""" 
            if x[0].score >= len(q)-(len(q)/10) and y[0].score >= len(c)-(len(c)/5):
                check = 1
                bcheck = copy.deepcopy(seq_record)

                "find the last position of alignemtn of fw primer.."
                rr = (x[0].aligned)[1][-1][-1]
                "find the last position of alignemtn of rw primer"
                rrr = (y[0].aligned)[1][0][0]
                ext = int(reader[i]['ext'])
                bcheck.letter_annotations = {}
                bcheck.seq = seq_record.seq[rr+ext:rrr]
                b = int(reader[i]['num'])
                bcheck.annotations["barcode number"]= b
                bcheck.annotations["length"] = "%i bp"%(len(bcheck.seq))

                if len(bcheck.seq)> 0:
                    bchecked.append(bcheck)
                    count += 1
                    nu = nuas(seq_record)
                    if nu != None:
                        if nuasc.get(nu) == None:
                            nuasc.update({nu: 1})
                        else:
                            su = nuasc.get(nu)
                            nuasc.update({nu:su+1})
                break
        if check == 0:
            countf += 1

    #SeqIO.write(bchecked, "%s/%s.faa"%(outpath,outname), "fasta")
    print ('number of reads is %i' %count0)
    print ('number of barcode checked is %i' %count)
    print ('number of failed reads is %i' %countf)

    return nuasc,blist,bchecked

"""Here it will read the CSV files and return a list"""
def read(barcode):
    with open(barcode, 'r') as f:
        barclist = list(csv.DictReader(f,delimiter=';', quotechar='"'))
    """the first number is the length of an empty cassette after PCR amlification and the second is the length of each UAS"""
    length = [83,116]
    return barclist,length

"""checking the lengh of the seq and report the number of UAS prediction"""
def nuas(seq):
    a = readlist[1]
    emp = a[0]
    ul = a[1]
    l = len(seq.seq)
    if (l-emp > -uer and l-emp < uer):
        return 0
    elif l-emp > uer:
        uscount = round((l-emp)/ul)
        return uscount
    else:
        return None

def out(out,bchecked,outname):
    with open('%s/%s.csv'%(out,outname), 'w',newline='') as csvfile:
        filewriter = csv.DictWriter(csvfile, fieldnames=['seqeunce','length','barcode number'],delimiter=';')
        filewriter.writeheader()
        for i in bchecked:
            filewriter.writerow({'seqeunce':i.seq,'length':i.annotations["length"],'barcode number':i.annotations["barcode number"]})
    csvfile.close()


"""this seqction is for loading the data:
l[0] = Sequencing file
l[1] = the barcode file
l[2] = the directory to save the output"""""
l = inq()
outpath = l[2]

""""readlist read csv files and return list
readlist[0] = barcode list
readlist[1] = defult DNA fragemnt lentgh information
"""
readlist = read(l[1])

"""this section is for actions:
pchecked = checks the primer bindign
bchecked = checks the barcode and trim the barcode and primer sequence
out =  export the result
"""

"""uer = UAS anealing error"""
uer = 15
pchecked = list() 
outname = str('5X.bchecked')
bchecked = bcheck(l[0],readlist[0],outname,outpath)
out(l[2],bchecked[2],outname)