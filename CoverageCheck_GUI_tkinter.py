from Tkinter import *
from tkFileDialog   import askopenfilename, askdirectory      
import CoverageCheck as CC

fields = 'Minimum Coverage', 'Maximum Strandbias'
bamfolder = None
bedfile = None

args = dict(bed=None,
            target_folder=None,
            min_dp=None,
            max_strand_ratio=None)

def fetch(entries):
  for entry in entries:
    field = entry[0]
    text  = entry[1].get()
    print('%s: "%s"' % (field, text)) 
  print "Bamfolder: %s" % bamfolder
  print "Bedfile: %s" % bedfile
  print "var: %s" % var

  
def makeform(root, fields):
  entries = []
  for field in fields:
    row = Frame(root)
    lab = Label(row, width=15, text=field, anchor='w')
    ent = Entry(row)
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent.pack(side=RIGHT, expand=YES, fill=X)
    entries.append((field, ent))
  return entries

def callback():
    name= askopenfilename() 
    print name

def fetch_bamfolder():
    global bamfolder
    selection = askdirectory()
    bamfolder.set(selection)
    args["bamfolder"] = selection

def fetch_bedfile():
    global bedfile
    selection = askopenfilename()
    bedfile.set(selection)
    args["bedfile"] = selection    

def run():
    bed = bedfile.get()
    target_folder = bamfolder.get()
    min_dp_value = int(min_dp.get())
    max_strand_ratio_value = int(max_strand_ratio.get())
    #min_dp = int(ents[0][1].get())
    #max_strand_ratio = int(ents[1][1].get())
  
    #print "Min DP %s" % min_dp
    #print "Bamfolder: %s" % args["bamfolder"] 
    #print "Bedfile: %s" % args["bedfile"]
    
    #root.destroy()
    CC.run(bed, target_folder, min_dp_value, max_strand_ratio_value)
    
if __name__ == '__main__':
  root = Tk()
  root.wm_title("CoverageCheck")

  #ents = makeform(root, fields)
  #root.bind('<Return>', (lambda event, e=ents: fetch(e)))   

  ##########################################
  # Select the min_dp
  
  min_dp = StringVar()
  min_dp.set(50)
  
  row = Frame(root)
  label = Label(row, width=15, text="Minimum Coverage", anchor='w')
  entry = Entry(row, textvariable=min_dp)
  row.pack(side=TOP, fill=X, padx=5, pady=5)
  label.pack(side=LEFT)
  entry.pack(side=RIGHT, expand=YES, fill=X)

  ##########################################
  # Select the max_strand_ratio
  
  max_strand_ratio = StringVar()
  max_strand_ratio.set(5)
  
  row = Frame(root)
  label = Label(row, width=15, text="Max. Strand Ratio", anchor='w')
  entry = Entry(row, textvariable=max_strand_ratio)
  row.pack(side=TOP, fill=X, padx=5, pady=5)
  label.pack(side=LEFT)
  entry.pack(side=RIGHT, expand=YES, fill=X)
  
  ##########################################
  # Show the selected Bamfolder
  
  bamfolder = StringVar()
  #bamfolder.set('Not selected')  
  bamfolder.set('/home/andreas/bioinfo/projects/nhs_coverage_tool/data/Alignment/')
  
  row = Frame(root)
  fixed_label = Label(row, width=15, text="Bamfolder", anchor='w')
  var_label = Label(row, textvariable = bamfolder)
  
  row.pack(side=TOP, fill=X, padx=5, pady=5)
  fixed_label.pack(side=LEFT)
  var_label.pack(side=RIGHT, expand=YES)

  ##########################################
  # Show the selected Bedfile

  bedfile = StringVar()
  #bedfile.set('Not selected')  
  bedfile.set('/home/andreas/bioinfo/projects/nhs_coverage_tool/data/Alignment/CLL_Dis_V1_TruSeq_CAT_Manifest_TC0031836-CAT.bed')
  
  row = Frame(root)
  fixed_label = Label(row, width=15, text="Bedfile", anchor='w')
  var_label = Label(row, textvariable = bedfile)
  
  row.pack(side=TOP, fill=X, padx=5, pady=5)
  fixed_label.pack(side=LEFT)
  var_label.pack(side=RIGHT, expand=YES)
  
  ###########################################
  # Create the buttons
  
  b1 = Button(root, text='Select bamfolder',command=fetch_bamfolder)
  b1.pack(side=LEFT, padx=5, pady=5)
  
  b2 = Button(root, text='Select region file', command=fetch_bedfile)
  b2.pack(side=LEFT, padx=5, pady=5)
  
  b2 = Button(root, text='Quit', command=root.quit)
  b2.pack(side=RIGHT, padx=5, pady=5)
  
  b1 = Button(root, text='Start CoverageCheck',command=run)
  b1.pack(side=RIGHT, padx=5, pady=5)
  
  ###########################################
  
  root.mainloop()