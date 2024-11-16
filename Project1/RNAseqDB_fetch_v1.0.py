# ---- load packages ---- #

import requests as rq
import os
import logging
from urllib.parse import urljoin
import pandas as pd
import argparse
from tqdm import tqdm

from dataclasses import dataclass, field, asdict
# from typing import Any, Optional

# ---- info ---- #

# version 1.1 - update on function: add new function extract
# expression data of specific treatment through specifying '--pattern' tag
# version 1.1.1 - 15.11.2024
# reorganize code and improve readability


# ---- params ---- #

@dataclass
class QueryMetaData:
    '''
    metadate regarding data fetching parameters
    '''
    db_url: str = field(default='https://plantrnadb.com/athrdb/',
                        metadata={'help': 'database url'})
    gene_list: str = field(default='',
                           metadata={'help': 'full path of gene list'})
    list_format: str = field(default='csv',
                             metadata={'help': 'data format of gene list'})
    list_sep: str = field(default=',',
                          metadata={'help': 'separator in gene list'})
    data_tag: str = field(default='tbox',
                          metadata={'help': 'data source to be downloaded'})
    data_pattern: str = field(default='',
                              metadata={'help': 'data type to be extracted'})
    out_dir: str = field(default=os.path.dirname(__file__))


@dataclass
class ResultsDataFrameParams:
    '''
    metadata regarding results formatting parameters
    '''
    regulation: pd.StringDtype = field(default=pd.StringDtype(),
                                       metadata={'help': 'indicator of gene regulation, either up or down'})  # noqa
    treatment_project: pd.StringDtype = field(default=pd.StringDtype())
    expression_fpkm: pd.Float32Dtype = field(default=pd.Float32Dtype(),
                                             metadata={'help': 'expression value in the form of FPKM'})  # noqa
    experiment: pd.StringDtype = field(default=pd.StringDtype(),
                                       metadata={'help': 'indicator of experimental conditions, either mock or treated'})  # noqa
    stderr: pd.Float32Dtype = field(default=pd.Float32Dtype(),
                                    metadata={'help': 'standard error calculated based on all data points for given treatment & project'})  # noqa
    avg_log2fc: pd.Float32Dtype = field(default=pd.Float32Dtype(),
                                        metadata={'help': 'average log2 foldchange of all data points'})  # noqa
    single_dpts: pd.Float32Dtype = field(default=pd.Float32Dtype(),
                                         metadata={'help': 'all single data points of expression data'})  # noqa


@dataclass
class RequestSessionConfiguration:
    accept: str = field(default='*/*')
    accept_encoding: str = field(default='gzip,deflate,br')
    accept_language: str = field(default='en-US,en;q=0.9')
    connection: str = field(default='keep-alive')
    host: str = field(default='plantrnadb.com')
    referer: str = field(default='https://plantrnadb.com/athrdb/')
    sec_fetch_dest: str = field(default='empty')
    sec_fetch_mode: str = field(default='cors')
    sec_fetch_site: str = field(default='same-origin')
    user_agent: str = field(default='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.4.1 Safari/605.1.15')  # noqa
    x_requested_with: str = field(default='XMLHttpRequest')
    newigv: str = field(default='no')
    newinfo: str = field(default='no')
    newplot: str = field(default='no')
    newtable: str = field(default='no')
    trans: str = field(default='yes')

    def SessionSettingGetter(self):
        return {
                'Accept': self.accept,
                'Accept-Encoding': self.accept_encoding,
                'Accept-Language': self.accept_language,
                'Connection': self.connection,
                'Host': self.host,
                'Referer': self.referer,
                'Sec-Fetch-Dest': self.sec_fetch_dest,
                'Sec-Fetch-Mode': self.sec_fetch_mode,
                'Sec-Fetch-Site': self.sec_fetch_site,
                'User-Agent': self.user_agent,
                'X-Requested-With': self.x_requested_with
                }

    def CookiesSettingGetter(self):
        return {
                'newigv': self.newigv,
                'newinfo': self.newinfo,
                'newplot': self.newplot,
                'newtable': self.newtable,
                'trans': self.trans
                }


# ---- utils ---- #


def ReadGenesfromList(gene_list, list_format, list_sep):
    '''
    extract gene names / IDs from gene_list file provided
    store them in list[genes]
    '''
    logger.info('Reading gene list is started')
    genes = []
    if gene_list == '':
        raise ValueError('The full path for gene list cannot be empty!')  # noqa
    if not isinstance(gene_list, str):
        raise ValueError(f'The full path for gene list must be given as str, but got {type(gene_list)}!')  # noqa
    if not os.path.exists(gene_list):
        raise FileNotFoundError(f'The given path for gene list "{gene_list}" cannot be found!')  # noqa
    if list_format not in ['csv', 'txt']:
        raise ValueError(f'Gene list can only be in the form of either ".csv" file or ".txt" file, but got {list_format}!')  # noqa
    if list_sep not in [',', ';', '\t']:
        raise ValueError(f'Delimiter can only be in the form of either comma, semicolon or tab, but got {list_sep}')
    try:
        with open(gene_list, 'r') as gl:
            while True:
                line = gl.readline()
                if line == '':
                    break
                gene_id = line.split(list_sep)
                gene_id = [gene.strip().upper() for gene in gene_id]
                genes += gene_id
    except Exception as e:
        logger.error(f'Reading of gene list failed: {e}')
        raise f'Gene list cannot be read correctly, please make sure the file is not corrupted' from e
    logger.info('Reading gene list is finished')
    return genes


def ConductPreSearch(gene, session, db_url):
    '''
    simulate maunal search in the browser
    to trigger data deposit as cache / cookies
    enabling further fetching process
    '''
    logger.info(f'Pre-search for gene {gene} is started (1/4)')
    server_url = urljoin(db_url, 'Search.php')
    query_value = f'{gene}/0/max////'
    params = {'query': query_value}
    try:
        # session = rq.Session()
        response = session.get(server_url, params=params, timeout=10)
        # logger.info(f'Conducting pre-search for gene {gene} in database {db_url}...')  # noqa
        if response.status_code == 200:
            logger.info(f'Pre-search is done for gene {gene}')
        else:
            logger.error(f'Unexpected status code got when carring out pre-search for gene {gene}: {response.status_code}')  # noqa
            raise ConnectionError(f'STEP 1 - presearch for gene {gene} failed')
    except Exception as e:
        #logger.error(f'Error occured when carrying out pre-search for gene {gene} from database')  # noqa
        raise e
    logger.info(f'Pre-search for gene {gene} is finished (1/4)')
    return session


def FetchGeneData(session, gene, db_url, data_tag):
    content = None
    target_url = f'{db_url}user/{gene}.{data_tag}'
    logger.info(f'Fetching process for gene {gene} is started (2/4)')
    try:
        response = session.get(target_url)
        if response.status_code == 200:
            # logger.info(f'Successfully fetched data for gene {gene}')
            content = response.content.decode('utf-8')
        else:
            logger.error(f'Data for gene {gene} cannot be found.')
            raise RuntimeError(f'STEP2 - data fetch for gene {gene} failed')
    except Exception as e:
        #logger.error(f'Error occured when fetching data for gene {gene}: {e}')
        raise e
    logger.info(f'Fetching process for gene {gene} is finished (2/4)')
    return content


def FormatRawDataFromDb(gene, content, data_pattern):
    results = None
    pattern_list = []
    conditions = ('mock', 'treated')
    schema = asdict(ResultsDataFrameParams())
    results = {}
    results_dtype = {}
    logger.info(f'Results formatting for gene {gene} is started (3/4)')
    if content is None:
        logger.warning(f'Results formatting step for gene {gene} is skipped since no data was retrieved from db.')  # noqa
        raise ValueError(f'STEP 3 - Results formatting for gene {gene} failed')
    else:
        rep_dict = {'\nup': ',up',
                    '\ndown': ',down',
                    '\n': ''}
        for k, v in rep_dict.items():  # replace ununiform format in the dataset for further data processing # noqa
            content = content.replace(k, v)
    items = content.split(',')
    if len(items) != 21:
        logger.warning(f'Results formatting step for gene {gene} is skipped since data retrieved is incomplete and has {len(items)} cells!')  # noqa
        raise ValueError(f'STEP 3 - Results formatting for gene {gene} failed')  # noqa
    if data_pattern == '':
        pass  # pattern mode off
    else:
        if '_' in data_pattern:  # pattern mode on
            pattern_list = data_pattern.split('_')  # multi-pattern mode
        else:
            pattern_list.append(data_pattern)  # single pattern mode
    if len(pattern_list) == 0:
        for n, title in enumerate(schema):
            column = title
            idx_row_up = n + len(schema)
            idx_row_down = n + (2 * len(schema))
            items_up = items[idx_row_up]
            items_down = items[idx_row_down]
            if n in [1, 5]:
                elements_up = items_up.split(';')
                elements_down = items_down.split(';')
                len_elements_up = len(elements_up)
                len_elements_down = len(elements_down)
                elements = elements_up + elements_down
                results[column] = elements
                results_dtype[column] = pd.StringDtype() if n == 1 else pd.Float32Dtype()
            elif n in [2, 4, 6]:
                items_up = items_up.replace('_', ';')
                items_down = items_down.replace('_', ';')
                elements_up = items_up.split(';')
                elements_down = items_down.split(';')
                for m, condition in enumerate(conditions):
                    column = f'{title}_{condition}'
                    elements_up_col = elements_up[m * len_elements_up:(m + 1) * len_elements_up]
                    elements_down_col = elements_down[m * len_elements_down:(m + 1) * len_elements_down]
                    elements = elements_up_col + elements_down_col
                    results[column] = elements
                    results_dtype[column] = pd.StringDtype() if n == 6 else pd.Float32Dtype()
    else:
        up_treatment = items[8]  # cell 8 contains treatment_project info for up-regulated expression data # noqa
        up_treat_items = up_treatment.split(';')
        down_treatment = items[15]  # cell 15 for down-regulated expression data # noqa
        down_treat_items = down_treatment.split(';')
        kw_sum = []
        gene_sum = []
        for pattern in pattern_list:
            up_extract = [pattern.lower() in item.lower() for item in up_treat_items]  # noqa
            down_extract = [pattern.lower() in item.lower() for item in down_treat_items]  # noqa
            for n, title in enumerate(schema):
                column = title
                idx_row_up = n + len(schema)
                idx_row_down = n + (2 * len(schema))
                items_up = items[idx_row_up]
                items_down = items[idx_row_down]
                if n in [1, 5]:
                    elements_up = items_up.split(';')
                    elements_down = items_down.split(';')
                    if len(elements_up) != len(up_extract):
                        logger.error(f'length of up-regulated treatment is NOT equal to that of {title} for gene {gene}')  # noqa
                        raise ValueError(f'STEP 3 - Results formatting for gene {gene} failed')
                    if len(elements_down) != len(down_extract):
                        logger.error(f'length of down-regulated treatment is NOT equal to that of {title} for gene {gene}')
                        raise ValueError(f'STEP 3 - Results formatting for gene {gene} failed')  # noqa
                    if n == 1:
                        len_elements_up = len(elements_up)
                        len_elements_down = len(elements_down)
                    up_extract_items = [item for item, flag in zip(elements_up, up_extract) if flag]  # noqa
                    down_extract_items = [item for item, flag in zip(elements_down, down_extract) if flag]  # noqa
                    extract_items = up_extract_items + down_extract_items
                    if n == 1:
                        extract_items_len = len(extract_items)
                        gene_col = [gene] * extract_items_len
                        gene_sum = gene_sum + gene_col
                        kw_col = [pattern] * extract_items_len
                        kw_sum = kw_sum + kw_col
                    if column in results.keys():
                        for it in extract_items:
                            results[column].append(it)
                    else:
                        results[column] = extract_items
                        results_dtype[column] = pd.StringDtype() if n == 1 else pd.Float32Dtype()
                elif n in [2, 4, 6]:
                    items_up = items_up.replace('_', ';')
                    items_down = items_down.replace('_', ';')
                    elements_up = items_up.split(';')
                    elements_down = items_down.split(';')
                    for m, condition in enumerate(conditions):
                        column = f'{title}_{condition}'
                        elements_up_col = elements_up[m * len_elements_up:(m + 1) * len_elements_up]
                        elements_down_col = elements_down[m * len_elements_down:(m + 1) * len_elements_down]
                        if len(elements_up_col) != len(up_extract):
                            logger.error(f'length of up-regulated treatment in {condition} is NOT equal to that of {title} for gene {gene}')
                            raise ValueError(f'STEP 3 - Results formatting for gene {gene} failed')
                        if len(elements_down_col) != len(down_extract):
                            logger.error(f'length of down-regulated treatment in {condition} is NOT equal to that of {title} for gene {gene}')
                            raise ValueError(f'STEP 3 - Results formatting for gene {gene} failed')
                        up_extract_items = [item for item, flag in zip(elements_up_col,up_extract) if flag]
                        down_extract_items = [item for item, flag in zip(elements_down_col,down_extract) if flag]
                        extract_items = up_extract_items + down_extract_items
                        if column in results.keys():
                            for item in extract_items:
                                results[column].append(item)
                        else:
                            results[column] = extract_items
                            results_dtype[column] = pd.StringDtype() if n == 6 else pd.Float32Dtype()
        results['keyword'] = kw_sum
        results['gene'] = gene_sum
        up_extract = []
        down_extract = []
    logger.info(f'Results formatting for gene {gene} is finished (3/4)')
    return results, results_dtype


def WriteFormattedResults(results, dtype, data_pattern, gene):
    logger.info(f'Formatted results writing for gene {gene} is started (4/4)')
    file_name = f'{QueryMetaData.out_dir}/{data_pattern}_RNAseq_data.csv' if data_pattern != '' else f'{QueryMetaData.out_dir}/{gene}_RNAseq_data.csv'
    pd_df_results = pd.DataFrame(results)
    copy = pd_df_results.astype(dtype=dtype, copy=True)
    if os.path.exists(file_name):
        csv = copy.to_csv(index=False, header=False)
    else:
        csv = copy.to_csv(index=False)
    try:
        logger.info(f'Trying to write formatted results in csv file with path {file_name}')
        with open(file_name, 'a') as w_file:
            w_file.write(csv)
    except Exception as e:
        logger.error(f'Error occured when writing formatted results into given path: {e}]')
        raise e(f'STEP 4 - results writing for gene {gene} failed')
    logger.info(f'Formatted results writing for gene {gene} is finished (4/4)')


class InfoTablePrinter:
    main_title = 'INFO TABLE'
    sub_titles = []
    len_titles = []
    details = []
    len_details = []
    fmt_metainfo = {}

    def __init__(self, metadata: QueryMetaData):
        self.metadata = metadata
    
    @property
    def OriMetaInfo(self):
        return {
                'database': self.metadata.db_url,
                'gene list': self.metadata.gene_list,
                'format': self.metadata.list_format,
                'delimiter': self.metadata.list_sep,
                'data source': self.metadata.data_tag,
                'treatment': self.metadata.data_pattern if self.metadata.data_pattern != '' else 'none',
                'output dir': self.metadata.out_dir
                }

    def PaddingPrint(self):
        ori_metainfo = self.OriMetaInfo
        padding_char = ' '
        for k, v in ori_metainfo.items():
            self.sub_titles.append(k)
            self.details.append(v)
            len_title = len(k)
            len_detail = len(v)
            self.len_titles.append(len_title)
            self.len_details.append(len_detail)
        len_title_max = max(self.len_titles)
        len_detail_max = max(self.len_details)
        for i, title in enumerate(self.sub_titles):
            # reset padding
            padding_title = ''
            padding_detail = ''
            # padding titles
            num_padding_title = len_title_max - self.len_titles[i]
            padding_title = padding_char * num_padding_title if num_padding_title > 0 else ''  # noqa
            new_title = title + padding_title
            # padding details
            num_padding_detail = len_detail_max - self.len_details[i]
            padding_detail = padding_char * num_padding_detail if num_padding_detail > 0 else ''  # noqa
            new_detail = padding_detail + self.details[i]
            self.fmt_metainfo[new_title] = new_detail
        # padding main title
        num_padding_mtitle_left = ((len_title_max + len_detail_max + 5) - len(self.main_title)) // 2  # noqa
        num_padding_mtitle_right = (len_title_max + len_detail_max + 5) - len(self.main_title) - num_padding_mtitle_left  # noqa
        fmt_mtitle = '|' + '#' * (num_padding_mtitle_left - 1) + self.main_title + '#' * (num_padding_mtitle_right - 1) + '|'  # noqa
        # set upper and lower borders
        top_bottom = '-' * (len_title_max + len_detail_max + 5)
        # print info table
        print(top_bottom)
        print(fmt_mtitle)
        for k, v in self.fmt_metainfo.items():
            print(f'|{k}:  {v}|')
        print(top_bottom)


def SummaryPrinter(suc_genes:list, fail_genes:dict):
    if len(suc_genes) != 0:
        print(f'Successful data fetching for following genes:')
        info = ''
        for i, gene in enumerate(suc_genes):
            info = info + gene + ', ' if i < len(suc_genes) - 1 else info + gene
        print(info)
    else:
        print(f'None of genes provided succeeded')
    print('Data fetching for following genes failed (gene ID | reason):')
    for k, v in fail_genes.items():
        print(f'{k} | {v}')

def main():
    suc_genes = []
    fail_genes = {}
    # print metainfo table
    printer = InfoTablePrinter(QueryMetaData)
    printer.PaddingPrint()
    # initialize session and set headers / cookies
    print('Initializing request session settings...')
    logger.info('Request session initialization is started')
    try:
        session = rq.session()
        session_config = RequestSessionConfiguration()
        sessionHeaderArgs = session_config.SessionSettingGetter()
        session.headers.update(sessionHeaderArgs)
        sessionCookiesArgs = session_config.CookiesSettingGetter()
        for k, v in sessionCookiesArgs.items():
            session.cookies.set(k, v)
        print('Request session is successfully established')
        logger.info('Request session initialization is started is done')
    except Exception as e:
        print('Process is stopped since the initialization failed')
        logger.error(f'Error occured in request session initialization step: {e}')  # noqa
        raise e
    # read gene list
    genes = ReadGenesfromList(gene_list=QueryMetaData.gene_list,
                              list_format=QueryMetaData.list_format,
                              list_sep=QueryMetaData.list_sep)
    print(f'In total {len(genes)} genes are recognized from which data will be fetched')
    # gene-wise data fetch and results format
    for gene in tqdm(genes, desc='Iterating over all genes', unit='gene'):
        try:
            session = ConductPreSearch(gene=gene,
                                    session=session,
                                    db_url=QueryMetaData.db_url,)
            content = FetchGeneData(session=session,
                                    gene=gene,
                                    db_url=QueryMetaData.db_url,
                                    data_tag=QueryMetaData.data_tag)
            results, dtype = FormatRawDataFromDb(gene=gene,
                                                content=content,
                                                data_pattern=QueryMetaData.data_pattern)
            WriteFormattedResults(results=results,
                                  dtype=dtype,
                                  gene=gene,
                                  data_pattern=QueryMetaData.data_pattern)
            suc_genes.append(gene)
        except Exception as e:
            fail_genes[gene] = e
    SummaryPrinter(suc_genes=suc_genes,
                   fail_genes=fail_genes)


if __name__ == '__main__':
    # intialize parser
    parser = argparse.ArgumentParser(description='Fetch RNA expression data from online database.')  # noqa
    # add parser arguments
    parser.add_argument('--base-url', type=str, default=None,    # noqa
                        help='Base URL of the RNAseq database (default: https://plantrnadb.com/athrdb/)')    # noqa
    parser.add_argument('--gene-list', type=str, required=True,
                        help='Full path to the gene list file')
    parser.add_argument('--list-format', type=str, default='csv',
                        help='Format of the gene list file, either csv or txt (default: csv) ')  # noqa
    parser.add_argument('--list-sep', type=str, default=',',
                        help='Delimiter used in the gene list, either ",", ";" or "\t" (default: ",")')  # noqa
    parser.add_argument('--data-tag', type=str, default='tbox',
                        help='Type of desired information, e.g. "tbox" for expression data under different treatments (default: tbox)')  # noqa
    parser.add_argument('--data-pattern', type=str, default='',
                        help='Motif(s) specifies experimental condition from which data will be extracted, e.g. flg22')  # noqa
    parser.add_argument('--out-dir', type=str, default=os.path.dirname(__file__),  # noqa
                        help='Output directory for final results (default: current script directory)')  # noqa
    args = parser.parse_args()
    # change metainfo accordingly
    if args.base_url is None:
        pass  # no databse given, default database remains
    if args.list_format in ['csv', 'txt']:
        QueryMetaData.list_format = args.list_format  # noqa
    else:
        raise ValueError(f'Gene list format can only be either csv or txt file, but got {args.list_format}!')  # noqa
    if args.list_sep in [',', ';', '\t']:
        QueryMetaData.list_sep = args.list_sep
    else:
        raise ValueError(f'Gene list delimiter can only be either comma, semicolon or tab, but got {args.list_sep}!')  # noqa
    QueryMetaData.gene_list = args.gene_list
    QueryMetaData.data_tag = args.data_tag
    QueryMetaData.data_pattern = args.data_pattern
    QueryMetaData.out_dir = args.out_dir
    # initialize logger
    logging.basicConfig(filename=os.path.join(QueryMetaData.out_dir,'message.log'),
                        format='%(asctime)s - **%(levelname)s**: %(message)s',
                        level=logging.INFO)
    logger = logging.getLogger(__name__)
    # start main loop
    main()