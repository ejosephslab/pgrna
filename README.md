# pgrna




#### A computational toolkit for the design of "polyvalent" guide RNAs (pgRNAs) that are optimized for CRISPR activity at multiple viral sites, simultaneously, while also avoiding interactions with the host genome or transcriptome.

This repository is intended to accompany our manuscript. For more information please refer to:

Polyvalent Guide RNAs for CRISPR Antivirals  
Rammyani Bagchi, Rachel Tinker-Kulberg, Tinku Supakar, Sydney Chamberlain, Ayalew Ligaba-Osena, and Eric A. Josephs*  
BIORXIV/2021/430352

## Installation instructions

### Software needed to be installed for running the code in Python:

1) Install Python 
      
      Go to the Anaconda Distribution Page at [Anaconda Installation](https://www.anaconda.com/products/individual). Click Download and select the latest Python version. Please ensure that you check the box that says "add to PATH" when installing on PC.

      Libraries/packages that needed to be installed in Python to run the code:
      
      1) pandas - pandas library is available as a part of the latest anaconda package. It can also be installed using conda and pip, the two main tools that install python packages.
      - Using conda – Type, conda install pandas, in the command prompt.
      - Using pip – If you are using pip, type !pip install pandas in the Jupyter Notebook App or type pip install pandas in the command prompt.
      2) biopython - It can be installed using pip and conda. Detailed instructions is available at [Biopython Installation](https://biopython.org/wiki/Packages).

      3) NCBI Standalone Blast
      - Download and install BLAST 2.8.1+ installer for your machine which is available from NCBI at [Blast Executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
      - Set up a searchable database for the Blast search of the host genome or transcriptome (e.g., from the NCBI, [Human Genome Sequences](https://www.ncbi.nlm.nih.gov/genome/guide/human/)) using the command (for host genome or transcriptome in FASTA format, here named hostdb): makeblastdb -in hostdb.fna -dbtype nucl -parse_seqids  -out hostdb. 
      - Add the created database to the system path, and re-direct the blastn query from the code (python) to find blast hits to the directory with the host database.


     4) Vienna RNA Package
      - Install the RNAfold installer compatible for your operating system available at [The ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/). Add the path of the package installation directory to your PATH variable manually.
      - Use the subprocess call command from the script (python) to make queries in the RNAfold software.

### How the pgRNA design code works:

Input
1) To run the code, you need to provide an initial input in the form of a csv file with three columns that contain: potential gRNAs, their predicted activity scores, and positions relative to a target sequence. Our code was designed to use the output files for cas13designer (https://gitlab.com/sanjanalab/cas13 or https://cas13design.nygenome.org/) (Ref. 1) and GPP sgRNA Design tool (https://portals.broadinstitute.org/gpp/public/) (Refs. 2 and 3) for RfxCas13d and SpyCas9, respectively. A model input for SARS-CoV-2 DNA genome sequence is provided as an example as [CovidCasRxguides.csv](https://github.com/ejosephslab/pgrna/blob/main/cas13pgrna/Input/CovidCasRxguides.csv).
2) You will need to edit the .py file to redirect the code to your input file by changing the line:
```
df = pd.read_csv("Input/CovidCasRxguides.csv")
```
3) To run the Cutting Frequency Determination (CFD) calculations, you also need to input two files provided in the Input folder. For example for Cas13d design, [Mismatch_Scores.csv](https://github.com/ejosephslab/pgrna/blob/main/cas13pgrna/Input/Mismatch_Scores.csv) and [Adjacent_Mismatch.csv](https://github.com/ejosephslab/pgrna/blob/main/cas13pgrna/Input/Adjacent_Mismatch.csv) is needed as input. 
```
#### CFD Scoring part

# In[22]:

df_mismatch = pd.read_csv("Input/Mismatch_Scores.csv")
#df_mismatch

# In[23]:

df_mismatch_scores = df_mismatch.set_index('Pair')
df_mismatch_scores

# In[24]:

df_adjacent_mismatch = pd.read_csv("Input/Adjacent_Mismatch.csv")
df_adjacent_mismatch
```
Output will have two main elements:

1) Initially it provides a csv containing homologous target pairs which are in the top quartile. The pairs are sorted according to the respective guide efficiency score. The model output is shown as a csv file, [df_pairs_sorted_top_quartile.csv](https://github.com/ejosephslab/pgrna/blob/main/cas13pgrna/Output/df_pairs_sorted_top_quartile.csv). 
2) The final output is a csv file which gives a list of potential gRNAs with following characteristics:
- High levels of predicted activity at two homeologous targets in a given viral genome.
- Low predicted relative activity at potential host “off-targets”.
- Biophysical characteristics such as high GC content, direct repeats, secondary structure free energy etc. that suggest high CRISPR activity for potential antiviral application.

The model output can be referred at [df_crRNA_withhits.csv](https://github.com/ejosephslab/pgrna/blob/main/cas13pgrna/Output/df_crRNA_withhits.csv).

### Running the code:

Clone the repo by typing following command in the terminal
```
git clone https://github.com/ejosephslab/pgrna

```
Navigate to the folder,cas13pgrna/cas9pgrna containing the python files/Cas13_pgRNA.py or Cas9pgRNA.py or the jupyter notebook files/Cas13_pgRNA.ipynb or Cas9_pgRNA.ipynb
```
cd path/to/cas13pgrna

For example: cd C:\Users\tinku\Documents\GitHub\pgrna\cas13pgrna
```

The code can be then run in following two ways:
1) Run in the Jupyter Notebook - The Jupyter Notebook App can be launched by typing jupyter notebook in the terminal. Once Jupyter Notebook is launched, it will automatically open a browser window and will show all the files in cas13pgrna/cas9pgrna folder. Click on Cas13_pgRNA.ipynb file. Run each cell in that notebook file to get all the outputs.
2) Running in the terminal/cmd for windows - Run the python file in the terminal (cmd) by typing, python Cas13_pgRNA.py or Cas9_pgRNA.py

```
For example: C:\Users\tinku\Documents\GitHub\pgrna\cas13pgrna>python Cas13_pgRNA.py
```
### Common errors:

1) ModuleNotFoundError: No module named 'pandas' when import pandas – Make sure you have pandas installed by typing pandas –-h in the terminal(cmd).
2) ImportError: No module named RNA – Check your RNAfold installation by typing RNAfold –-h in the terminal(cmd).
3) BLAST Database error: No alias or index file found for nucleotide database – Make sure you have given the correct path to the database while making the blastn query from your python script. This has to be done, for example, at line 1022 in the Cas13_pgRNA.py file (here for a database named 'Human_NCBI_rnadb').
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
### References for CRISPR activity prediction used for inputs:

(1) Hans-Hermann Wessels*, Alejandro Méndez-Mancilla*, Xinyi Guo, Mateusz Legut, Zharko Daniloski, Neville E. Sanjana. (2020)
Massively parallel Cas13 screens reveal principles for guide RNA design. Nature Biotechnology. 38, 722–727.

(2) Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan, K. F., … Root, D. E. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature biotechnology, 34(2), 184-191.

(3) Sanson, K. R., Hanna, R. E., Hegde, M., Donovan, K. F., Strand, C., Sullender, M. E., … Doench, J. G. (2018). Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities. Nature communications, 9(1), 5416.

---------------------------------------
 Last edit: February 25 2021  
 Software versions 0.9.0  
 
 Josephs lab  
 Department of Nanoscience  
 Joint School of Nanoscience and Nanoengineering  
 University of North Carolina at Greensboro
 
 *If you have any questions, comments, or issues, please contact Eric using the email address found in the accompanying manuscript.
