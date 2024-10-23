#---- load packages ----#
import requests as rq
import os
import logging
from urllib.parse import urljoin
import pandas as pd
import argparse
#---- info ----#

# version 1.1 - update on function: add new function extract expression data of specific treatment through specifying '--pattern' tag

#---- functions ----#

class RNAseqInfo_query:
    
    def __init__(self, base_url='https://plantrnadb.com/athrdb/',gene_list='', list_format='csv', list_sep=',', tag='tbox', pattern = '', dest='/Volumes/Expansion/RNA_seq_data'):
        self.base_url = base_url
        self.gene_list = gene_list
        self.list_format = list_format
        self.list_sep = list_sep
        self.sep_name = 'comma' if self.list_sep == ',' else 'semicolon' if self.list_sep == ';' else 'tab' if self.list_sep == '\t' else 'cannot recognized'
        self.tag = tag
        self.dest = dest
        self.pattern = pattern

        self.info_dict = {'database':self.base_url,
                          'gene list': self.gene_list,
                          'format': self.list_format,
                          'delimiter': self.sep_name,
                          'data source': self.tag,
                          'results location': self.dest}
        
        self.genes = []
        self.suc_genes = []
        self.fail_genes = []
        self.target_url = ''
        self.content = None
        
        self.schema = {'Regulation':pd.StringDtype,
                'Treatment_project':pd.StringDtype,
                'Exprssion_fpkm':pd.Float32Dtype,
                'Experiment':pd.StringDtype,
                'StdErr':pd.Float32Dtype,
                'Avg_log2foldchange':pd.Float32Dtype,
                'Single_datapoints':pd.Float32Dtype}
        self.condition = ['mock','treated']
        
        self.session = rq.Session()
        self.session.headers.update({
            'Accept':'*/*',
            'Accept-Encoding':'gzip,deflate,br',
            'Accept-Language':'en-US,en;q=0.9',
            'Connection':'keep-alive',
            'Host':'plantrnadb.com',
            'Referer':'https://plantrnadb.com/athrdb/',
            'Sec-Fetch-Dest':'empty',
            'Sec-Fetch-Mode':'cors',
            'Sec-Fetch-Site':'same-origin',
            'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.4.1 Safari/605.1.15',
            'X-Requested-With':'XMLHttpRequest'
        })
        initial_cookies = {
            'newigv':'no',
            'newinfo':'no',
            'newplot':'no',
            'newtable':'no',
            'trans':'yes'
        }

        for key,value in initial_cookies.items():
            self.session.cookies.set(key,value)
        logging.info('Initial cookies set')

    def formatter_print_start_info(self):
        """
        print info table in the console
        """
        
        main_title = 'INFO TABLE'
        titles =[]
        len_titles = []
        items = []
        len_items = []
        new_formatted_infoD = {}
        for i, (k,v) in enumerate(self.info_dict.items()):
            len_title = len(k)
            len_item = len(v)
            titles.append(k)
            len_titles.append(len_title)
            items.append(v)
            len_items.append(len_item)
        len_title_max = max(len_titles)
        len_item_max = max(len_items)
        
        for i, title in enumerate(titles):
            num_space_title = len_title_max - len_titles[i]
            fill_title = ' ' * num_space_title if num_space_title > 0 else ''
            new_title = title + fill_title

            num_space_item = len_item_max - len_items[i]
            fill_item = ' ' * num_space_item if num_space_item > 0 else ''
            new_item = fill_item + items[i]

            new_formatted_infoD[new_title] = new_item

        main_title_fill_left = ((len_item_max+len_title_max+5) - len(main_title)) // 2
        main_title_fill_right = (len_item_max+len_title_max+5) - len(main_title) - main_title_fill_left
        formatted_main_title = '|' + '#' * (main_title_fill_left - 1) + main_title + '#' * (main_title_fill_right - 1) + '|'
        top_bottom = '-' * (len_item_max+len_title_max+5)
        print(top_bottom)
        print(formatted_main_title)
        for key, value in new_formatted_infoD.items():
            print(f'|{key}:  {value}|')
        print(top_bottom)

    def log_start_info(self):
        logging.info('# INFO TABLE #\ndatabase:{}\ngene list:{}\nformat:{}\ndelimiter:{}\ndata source:{}\nresults location:{}'.format(self.base_url, self.gene_list, self.list_format, self.list_sep, self.tag, self.dest))

    def read_list(self):
        """
        read in gene list for which data need to be fetched
        """

        if self.gene_list == '':
            raise ValueError('Gene path cannot be empty, please check!')
        if not os.path.exists(self.gene_list):
            raise FileNotFoundError('path cannot be found, pleae check!')
        if self.list_format not in ['csv', 'txt']:
            raise ValueError('list can only be in the form of either CSV file or TXT file!')
        if self.list_sep not in [',',';','\t']:
            raise ValueError('list delimiter can only be either ",", ";" or "\\t"!')
        
        try:
            with open(self.gene_list, 'r') as gl:
                while True:
                    line = gl.readline().strip().upper()
                    if line == '':
                        break
                    self.genes.append(line)
                print(self.genes)
        except Exception as e:
            raise e
        return self.genes

    def pre_search(self,gene):
        """
        Simulate a manual search through interaction with 'Search.php' in the browser
        This should set necessary cookies/cache on the server
        """
        search_url = urljoin(self.base_url, 'Search.php')
        query_value = f'{gene}/0/max////'
        params = {'query':query_value}

        try:
            response = self.session.get(search_url, params=params)

            if response.status_code == 200:
                logging.info(f'Successful pre-search for gene {gene}')
            else:
                logging.warning(f'Unexpected status code {response.status_code}')
        except rq.exceptions.HTTPError as httperr:
            logging.error(f'HTTP error occured while searching for gene {gene}:{httperr}')
        except Exception as e:
            logging.error(f'Error occured while searching for gene {gene}:{e}')


    
    def get_content(self, gene):
        """
        get desired data from database
        """
        self.content = None
        target_url = f'{self.base_url}user/{gene}.{self.tag}'
        logging.info(f'Trying to fetch desired data for gene {gene} from databse')
        try:
            response = self.session.get(target_url)
            if response.status_code == 200:
                logging.info(f'Successfully fetched data for gene {gene}')
                self.content = response.content.decode('utf-8')
            else:
                logging.warning(f'Requested data for gene {gene} cannot be found.')
                self.content = None
        except rq.exceptions.HTTPError as http_err:
            print(f'HTTP error occurs:{http_err}')
            self.content = None
        except Exception as e:
            print(f'Error: fecth data from server failed: {e}')
            self.content = None

    def extract_conditioned_info(self, gene):
        """
        only extract the expression data of given condition/treatment, e.g. flg22, for each gene
        """
        pattern_list = []
        if self.pattern == '':
            self.write_results(gene)
        else:
            pattern_list.append(self.pattern)
            file_name = f'{self.dest}/{self.pattern}_RNAseq_data.csv'
            if '_' in self.pattern:
                pattern_list = self.pattern.split('_')
            if self.content == None:
                raise ValueError(f'No data has been retrieved for gene {gene}')
            ##try:
            up_extract = []
            down_extract = []
            up_extract_dict = {}
            down_extract_dict = {}
            
            df = {}
            col_dtype = {}
            rep_dict = {'\nup':',up',
                        '\ndown':',down',
                        '\n':''}
            for k,v in rep_dict.items():
                self.content = self.content.replace(k,v)
            items = self.content.split(',')
            up_treatment = items[8]
            down_treatment = items[15]
            up_treat_items = up_treatment.split(';')
            down_treat_items = down_treatment.split(';')
            if len(pattern_list) == 1:
                up_extract = [self.pattern in item for item in up_treat_items]
                down_extract = [self.pattern in item for item in down_treat_items]

                # function start
                for i,title in enumerate(self.schema):
                    column = title
                    index_row_up = i+len(self.schema)
                    index_row_down = i + (2*len(self.schema))
                    item_up = items[index_row_up]
                    item_down = items[index_row_down]
                    up_extract_items = None
                    down_extract_items = None
                    if i in [1,5]:
                        elements_up = item_up.split(';')
                        elements_down = item_down.split(';')
                        if len(elements_up) != len(up_extract):
                            raise ValueError(f'Step{i}:length of up-regulated treatment is NOT equal to that of {title} while processing gene {gene}')
                        if len(elements_down) != len(down_extract):
                            raise ValueError(f'Step{i}:length of down-regulated treatment is NOT equal to that of {title} while processing gene {gene}:\n{len(elements_down)}\n{len(down_extract)}')
                        len_elements_up = len(elements_up)
                        len_elements_down = len(elements_down)
                        try:
                            up_extract_items = [item for item, flag in zip(elements_up,up_extract) if flag]
                            down_extract_items = [item for item, flag in zip(elements_down,down_extract) if flag]
                            if i == 1:
                                extract_items_len = len(up_extract_items) + len(down_extract_items)
                                gene_col = [gene] * extract_items_len
                                df['gene'] = gene_col
                            extract_items = up_extract_items + down_extract_items
                            df[column] = extract_items
                            col_dtype[column] = pd.StringDtype() if i == 1 else pd.Float32Dtype()
                        except Exception as e:
                            print(f'in step{i}:{e}')
                    elif i in [2,4,6]:
                        item_up = item_up.replace('_',';')
                        item_down = item_down.replace('_',';')
                        elements_up = item_up.split(';')
                        elements_down = item_down.split(';')
                        for n, condition in enumerate(self.condition):
                            column = f'{title}_{condition}'
                            elements_up_col = elements_up[n * len_elements_up:(n + 1) * len_elements_up]
                            elements_down_col = elements_down[n * len_elements_down:(n + 1) * len_elements_down]
                            if len(elements_up_col) != len(up_extract):
                                raise ValueError(f'Step{i}:length of up-regulated treatment in {condition} is NOT equal to that of {title} while processing gene {gene}:\n{len(elements_up_col)}\n{len(up_extract)}')
                            if len(elements_down_col) != len(down_extract):
                                raise ValueError(f'Step{i}:length of down-regulated treatment in {condition} is NOT equal to that of {title} while processing gene {gene}:\n{len(elements_down_col)}\n{len(down_extract)}')
                            try:
                                up_extract_items = [item for item, flag in zip(elements_up_col,up_extract) if flag]
                                down_extract_items = [item for item, flag in zip(elements_down_col,down_extract) if flag]
                                both_updown_items = up_extract_items + down_extract_items
                                df[column] = both_updown_items
                                col_dtype[column] = pd.StringDtype() if i == 6 else pd.Float32Dtype()
                            except Exception as e:
                                print(f'in step{i}:{e}')
                    else:
                        pass
                # function end
            else:
                print('multiple patterns mode is still under construction')
                #for treatment in pattern_list:
                #    if up_extract == [] and down_extract == []:
                #        pass
                #    else:
                #        up_extract = []
                #        down_extract = []
                #    up_extract = [treatment in item for item in up_treat_items]
                #    up_extract_dict[treatment] = up_extract
                #    down_extract = [treatment in item for item in down_treat_items]
                #    down_extract_dict[treatment] = down_extract
                #    up_extract = []
                #    down_extract = []
                #    for key, value in up_extract_dict.items():
                #        up_extract = value
                #        down_extract = down_extract_dict[key]
                        #function start

                        # function end
            
            try:
                pd_df = pd.DataFrame(df)
                df_copy = pd_df.astype(dtype=col_dtype, copy=True)
                csv_pd_df = df_copy.to_csv(index=False)
                if os.path.exists(file_name):
                    csv_pd_df = df_copy.to_csv(index=False, header=False)
                with open(file_name, 'a') as w_file:
                    w_file.write(csv_pd_df)
            except Exception as e:
                print(f'before last step:{e}')
            ##except Exception as e:
                ##print(f'last step:{e}')

    def write_results(self, gene):
        """
        parse and reformat the data,store it in csv file
        """
        file_name = f'{self.dest}/{gene}_RNAseq_data.csv'
        if self.content == None:
            raise ValueError(f'No data has been retrieved for gene {gene}')
        
        logging.info(f'Data has been retrieved from database and is being written to csv table for gene {gene}')
        try:
            rep_dict = {'\nup':',up',
                        '\ndown':',down',
                        '\n':''}
            for k,v in rep_dict.items():
                self.content = self.content.replace(k,v) 
            items = self.content.split(',')
            df = {}
            col_dtype = {}
            for i,title in enumerate(self.schema):
                final_title = title
                index_row_a = i+len(self.schema)
                index_row_b = i+(2 * len(self.schema))
                item_a = items[index_row_a]
                item_b = items[index_row_b]
                if i in [1,5]:
                    elements_a = item_a.split(';')
                    elements_b = item_b.split(';')
                    len_elements_a = len(elements_a)
                    len_elements_b = len(elements_b)
                    elements = elements_a + elements_b
                    df[final_title] = elements
                    col_dtype[final_title] = pd.StringDtype() if i == 1 else pd.Float32Dtype()
                elif i in[2,4,6]:
                    item_a = item_a.replace('_',';')
                    item_b = item_b.replace('_',';')
                    elements_a = item_a.split(';')
                    elements_b = item_b.split(';')
                    for n, condition in enumerate(self.condition):
                        final_title = f'{title}_{condition}'
                        elements_a_col = elements_a[n * len_elements_a:(n + 1) * len_elements_a]
                        elements_b_col = elements_b[n * len_elements_b:(n + 1) * len_elements_b]
                        elements = elements_a_col + elements_b_col
                        df[final_title] = elements
                        col_dtype[final_title] = pd.StringDtype() if i == 6 else pd.Float32Dtype()
                else:
                    pass
            try:
                pd_df = pd.DataFrame(df)
                df_copy = pd_df.astype(dtype=col_dtype, copy=True)
                df_treat = df_copy['Treatment_project'].str.contains(self.pattern, na=False)
                csv_pd_df = df_copy.to_csv(index=False)
                with open(file_name, 'w') as w_file:
                    w_file.write(csv_pd_df)
            except Exception as e:
                logging.error(f'Pandas dataframe failed to be established: {e}')
        except Exception as e:
            logging.error(f'Writing process for gene {gene} failed: {e}')

    def print_final_info(self):
        logging.critical(f'Successfully retieved data from in total {len(self.suc_genes)} genes:\n{self.suc_genes}')

    def iter_gene_list(self):
        for gene in self.genes:
            self.pre_search(gene)
            self.get_content(gene)
            self.extract_conditioned_info(gene)



#---- main ----#
if __name__ == '__main__': 
    # intialize parser
    parser = argparse.ArgumentParser(description='Fetch RNA expression infomation from a database.')
    # add parser arguments
    parser.add_argument('--base-url', type=str, default='https://plantrnadb.com/athrdb/',
                        help='Base URL of the RNAseq database (default: https://plantrnadb.com/athrdb/)')
    parser.add_argument('--gene-list', type=str, required=True,
                        help='Full path to the gene list file')
    parser.add_argument('--list-format', type=str, default='csv',
                        help='Format of the gene list file, either csv or txt (default: csv) ')
    parser.add_argument('--list-sep', type=str, default=',',
                        help='Delimiter used in the gene list, either ",", ";" or "\t" (default: ",")')
    parser.add_argument('--tag', type=str, default='tbox',
                        help='The type of desired information, e.g. "tbox" stands for treatment. (default: tbox)')
    parser.add_argument('--pattern', type=str, default='',
                        help='The motif for which will be searched in treatment column to extract expression data of specific conditions, e.g. flg22')
    parser.add_argument('--dest', type=str, default=os.path.dirname(__file__),
                        help='The folder full path for formatted results (default: current script directory)')
    args = parser.parse_args()

    logging.basicConfig(filename='message.log',
                        format='%(asctime)s - **%(levelname)s**: %(message)s',
                        level=logging.INFO)

    instance = RNAseqInfo_query(base_url=args.base_url,
                                gene_list=args.gene_list,
                                list_format=args.list_format,
                                list_sep=args.list_sep,
                                tag=args.tag,
                                pattern=args.pattern,
                                dest=args.dest)
    instance.formatter_print_start_info()
    instance.log_start_info()
    try:
        instance.read_list()
        instance.iter_gene_list()
        instance.print_final_info()
    except ValueError as ve:
        print(f'Error: {ve}')
    except Exception as e:
        print(f'extract RNAseq expression failed: {e}')

"""
Overview of data structure from plantrnadb (after split by delimter ',': R=row, C=cell, "|"=column sepatator, "-"=row separator)
====================================================================================================================================================
|#R1#|--C1:group--|-------C2:tagname--------|------C3:FPKM------|------C4:infotag------|-----C5:se-----|--------C6:fc---------|-----C7:contlib-----|
----------------------------------------------------------------------------------------------------------------------------------------------------
|#R2#|----C8:up---｜--C9:treatment_project--｜--C10:mFPKM_tFPKM--|--C11:mock_treatment--|--C12:mSE_tSE--|--C13:log2foldchange--|--C14:singleValue--|
----------------------------------------------------------------------------------------------------------------------------------------------------
|#R3#|--C15:down--|--C16:treatment_project--｜--C17:mFPKM_tFPKM--|--C18:mock_treatment--|--C19:mSE_tSE--|--C20:log2foldchange--|--C21:singleValue--|
====================================================================================================================================================
"""