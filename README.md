# pgrna




#### A computational toolkit for the design of "polyvalent" guide RNAs (pgRNAs) that are optimized for CRISPR activity at multiple viral sites, simultaneously, while also avoiding interactions with the host genome or transcriptome.

## Installation instructions

### Software needed to be installed for running the code in Python:

1) Install Python 
      Go to the Anaconda Distribution Page at [Anaconda Installation](https://www.anaconda.com/products/individual). Click Download and select the latest Python version. Please ensure that you check the box that says "add to PATH" when installing on PC.

      Libraries/packages that needed to be installed in Python to run the code:
      
      1) pandas - pandas library is available as a part of the latest anaconda package. It can also be installed using conda and pip, the two main tools that install python packages.
      - Using conda – Type, conda install pandas, in the command prompt.
      - Using pip – If you are using pip, type !pip install pandas in the Jupyter Notebook App or type pip install pandas in the command prompt.
      2) biopython - It can be installed using pip and conda. Detailed instructions is available at [Biopython Installation](https://biopython.org/wiki/Packages)

      2) NCBI Standalone Blast
      - Download and install BLAST 2.8.1+ installer for your machine which is available from NCBI at [Blast Executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
      - Set up a searchable database for the Blast search of the host genome or transcriptome (e.g., from the NCBI https://www.ncbi.nlm.nih.gov/genome/guide/human/) using the command (for host genome or transcriptome in FASTA format, here named hostdb): makeblastdb -in hostdb.fna -dbtype nucl -parse_seqids  -out hostdb 
      - Add the created database to the system path, and re-direct the blastn query from the code (python) to find blast hits to the directory with the host database.


      3) Vienna RNA Package
      - Install the RNAfold installer compatible for your operating system available at https://www.tbi.univie.ac.at/RNA/. Add the path of the package installation directory to your PATH variable manually.
      - Use the subprocess call command from the script (python) to make queries in the RNAfold software.

### How the pgRNA design code works:

Input
1) To run the code, you need to provide an initial input in the form of a csv file which contains three columns mainly, potential gRNAs, their respective scores, and positions relative to a target sequence. The model input for SARS-CoV-2 DNA genome sequence is available for reference as CovidCasRxguides.csv:
2)  cas13designer, sgRNA design
3) To run the Cutting Frequency Determination (CFD) calculations, you need to input two csv files. For example for Cas13d design, Mismatch_Scores.csv and Adjacent_Mismatch.csv is needed as input. 

Output will have two main elements:

1) Initially it provides a csv containing homologous target pairs which are in the top quartile. The pairs are sorted according to the respective guide efficiency score. The model output is shown as a csv file, df_pairs_sorted_top_quartile.csv. 
2) The final output is a csv file which gives a list of potential gRNAs with following characteristics:
- High levels of predicted activity at two homeologous targets in a given viral genome.
- Low predicted relative activity at potential host “off-targets”.
- Biophysical characteristics such as high GC content, direct repeats, secondary structure free energy etc. that suggest high CRISPR activity for potential antiviral application.

The model output can be referred at df_crRNA_withhits.csv.

### Running the code:

Clone the repo by typing following command in the terminal
```
git clone https://github.com/ejosephslab/pgrna

```



The code can be run in following two ways:
1) Run in the Jupyter Notebook - The Jupyter Notebook App can be launched by typing in a terminal. Change the Jupyter Notebook startup folder by using the command cd /the_folder_name in the command prompt (Windows). Once Jupyter Notebook is launched, it will automatically open a browser window and will show the following page. Click on new and then on “Python 3”. Run each cell in the jupyter notebook to get all the outputs.
2) Running through the Windows in the command terminal (cmd) - navigate to the folder containing the python/Cas13_pgRNA.py (or python/Cas9_pgRNA.py) by typing , cd path/to/ Cas13_pgRNA.py.
Once you have navigated to the folder, run the python file in the windows terminal (cmd) by typing, 
python [Cas13pgRNA.py or Cas9pgRNA.py] [ViralcrRNA.csv] [hostdb]

### Common errors:

1) ModuleNotFoundError: No module named 'pandas' when import pandas – Make sure you have pandas installed by typing pandas –-h in the terminal(cmd).
2) ImportError: No module named RNA – Check your RNAfold installation by typing RNAfold –-h in the terminal(cmd).
3) BLAST Database error: No alias or index file found for nucleotide database – Make sure you have given the correct path to the database while making the blastn query from your python script. This has to be done at line 1022 in the Cas13_pgRNA.py file.
```
      # In[71]:


      command = 'blastn -query out_pgRNA_filtered.fsa -task blastn-short -db C:/Users/18064/Documents/Blast/db/Human_NCBI_rnadb -outfmt "7 qacc sacc qstart qend sstart send" -out pgRNA_NCBI_blast.out'
      subprocess.call(command, shell=True)


      # In[72]:
```
      
4) df_crRNA1["Free Energy"] = energy_list, length mismatch error, make sure you do not have a file named energyoutput.txt in the folder containg the python file/jupyter notebook file.

### Typical run time:

The run time depends on the length of the genome. For SARS-CoV-2 (a large RNA virus), typical run time ranges from 1.5 - 2 hours.

---------------------------------------
 Software version 0.1
 
 If you spot any bug in our code or if you have any questions or comments, please contact eajoseph@uncg.edu.

February 25 2021
