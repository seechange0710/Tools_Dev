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
- #### Main structure: RNAInfo_query = 
    - **function 'formatter_print_start_info'**: all paramters used while runing the code will be printed directly in the console as start info
    - **function 'log_start_info'**: same info from 'formatter_print_start_info' will be documanted in the log file 'message.log' in the same directory as the script
    - **function 'read_list'**: read the gene list provided in the program
    - **function 'presearch'**: perform a 'presearch' in the database by mimicing the behavior of a manual search via browser and get data loaded and cached on the server for downstream formal data fetch (without this step data fetch for genes will be rejected as code 404 return)
    - **function 'get_content'**: formal fetch for data needed
    - **function 'write_results'**: receive the content of data fetch from last step and reformatted it in a more reader-friendly form, final output will be stored into csv files
    - **function 'print_final_info'**: after the whole process, the number of genes whose data have been successfully fetched will be summed
    - **function 'iter_gene_list'**: process will be repeated for all genes in the gene list provided

- #### Parameters/Tags:
    - *--base-url*: you can also replace the default database with the one as you wish, since the code represented here is only a prototype and applicable for the default database, problems may occur when other databases are given due to different db and data structure
    - *--gene-list*: full path of gene list file which contains all genes you want to retrieve data from
    - *--list-format*: important for code correctly recognizing all genes in the gene list, can be only csv or text file.
    - *--list-sep*: delimiters used in gene list file, can only be comma, semicolon or tab.
    - *--tag*: tag for information that you need, e.g. in the default databse, a 'tbox' tag refers to gene expression data under different treatments which are visualized in a box plot.
    - *--dest*: full path of directory you want to store the formatted results in.
