# Development of tools for scientific use

Hi there, the main aim of this repository is to deposit scripts that are developed for easing data query and analysis processes in scientific research, hopefully it can help ones who need helps. I am also looking forward to new opportunities/challenges in scientific software development as well as AI/ML-powered application development, you are always welcome to contact me via email **sicheng.xu.work@gmail.com** if you are looking for someone who may solve problems in your research through programming or have a good idea and want to build a team to achieve it :).

## Project 1: fetch needed information from online RNAseq database

### 1) File name: 
RNAseqDB_fetch_v1.0.py

### 2) Problem:
The main problem is, my girlfriend... no I am kiding, but yeah she encountered the situation where she wanted to download expression data of her choice of genes (more than 1000 in a list) but the online database seems not to support batch-download. Do it one by one manually would be a disaster (for me also since I already see myself sitting the whole day in front of PC and collect data for her...)

### 3) Goal:
This script aims to communicate with a online RNAseq database and automate the query as well as data download processes for a bunch of genes one desires

### 4) Description:

- #### Main structure: RNAInfo_query
```
function 'formatter_print_start_info': all paramters used while runing the code will be printed directly in the console as start info
function 'log_start_info': same info from 'formatter_print_start_info' will be documanted in the log file 'message.log' in the same directory as the script
function 'read_list': read the gene list provided in the program
function 'presearch': perform a 'presearch' in the database by mimicing the behavior of a manual search via browser and get data loaded and cached on the server for downstream formal data fetch (without this step data fetch for genes will be rejected as code 404 return)
function 'get_content': formal fetch for data needed
function 'write_results': receive the content of data fetch from last step and reformatted it in a more reader-friendly form, final output will be stored into csv files
function 'print_final_info': after the whole process, the number of genes whose data have been successfully fetched will be summed
function 'iter_gene_list': process will be repeated for all genes in the gene list provided
````
- #### Parameters/Tags:
```
--base-url: you can also replace the default database with the one as you wish, since the code represented here is only a prototype and applicable for the default database, problems may occur when other databases are given due to different db and data structure
--gene-list: full path of gene list file which contains all genes you want to retrieve data from
--list-format: important for code correctly recognizing all genes in the gene list, can be only csv or text file.
--list-sep: delimiters used in gene list file, can only be comma, semicolon or tab.
--tag: tag for information that you need, e.g. in the default databse, a 'tbox' tag refers to gene expression data under different treatments which are visualized in a box plot.
--pattern: specify treatment(s) you want to extract expression data from, a list of treatments can be given by connecting with underscore.
--dest: full path of directory you want to store the formatted results in.
```

### 5) Examples:
**Whole_extract_mode:**
In whole extract mode, each gene will have its own output file that contains all data retrieved from online database.
<img width="1106" alt="start_info" src="https://github.com/user-attachments/assets/06db6733-f503-4865-ad6b-437de5572e72">
above is how start info table looks like when the script runs appropriately, gene IDs that are recognized by the script will be also display on the bottom.<br/><br/>
<img width="891" alt="formatted_results" src="https://github.com/user-attachments/assets/65db8c92-4a34-4a16-8807-91f0b616b69f"><br/>
formatted results consist of in total eight columns, including average gene expression in FPKM and log2foldchange between mock and treatment.


**Specified_extract_mode:**
In this mode, pattern(s) given and only the desired data associated with the given pattern will be documented. Data from all genes will be listed in one output file with file name being the pattern itself.

<img width="144" alt="Screenshot 2024-10-28 at 22 07 21" src="https://github.com/user-attachments/assets/ad0e0791-1cbb-4e22-b3b9-c4aad99b632b">(one pattern)
<img width="1080" alt="Screenshot 2024-10-28 at 22 07 56" src="https://github.com/user-attachments/assets/5ec76576-d8ff-46cc-9546-7607d85b79e9">

<img width="122" alt="Screenshot 2024-10-28 at 22 07 33" src="https://github.com/user-attachments/assets/a0616408-e5e5-412e-b8b6-6c4080e84ce6">(multiple patterns)
<img width="1004" alt="Screenshot 2024-10-28 at 22 05 48" src="https://github.com/user-attachments/assets/fa832d5d-87a5-44d6-9406-1e668167cb3d">

