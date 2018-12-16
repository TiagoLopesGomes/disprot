import re, gzip, time, itertools, io
from collections import OrderedDict
from multiprocessing import Pool
import numpy as np
import pandas as pd
import requests

# CONSTANTS

# APPRIS isoform annotations. See https://ensembl.org/info/genome/genebuild/transcript_quality_tags.html
appris_rank = OrderedDict([
    ('principal1', 1),
    ('principal2', 2),
    ('principal3', 3),
    ('principal4', 4),
    ('principal5', 5),
    ('alternative1', 6),
    ('alternative2', 7),
    (np.nan, 8)
])

predictorIdToName = OrderedDict([
    (1,  'VLXT'),
    (2,  'VSL2b'), 
    (3,  'PrDOS'),
    (4,  'IUPred-S'),
    (6,  'PV2'),
    (8,  'IUPred-L'),
    (9,  'Espritz-X'),
    (10, 'Espritz-N'),
    (11, 'Espritz-D')
])

predictorNameToId = OrderedDict([(v, k) for k, v in predictorIdToName.items()])

ensemblProteomeColumnMap = {
    'gene': 'ensembl_gene_id',
    'transcript': 'ensembl_transcript_id',
    'id': 'ensembl_peptide_id',
    'gene_symbol': 'hgnc_symbol'
}

uniprotProteomeColumnMap = {
    'id': 'uniprot_id',
    'gn': 'hgnc_symbol'
}

def createFileObject(file, mode=None):
    '''
    Create file object.
    
    Args
    - file: str
        Path to file. If extension ends with '.gz', gzip compression is assumed.
    - mode: str. default=None
        Mode to open file with.
    
    Returns: file object (io.TextIOBase)
    '''
    
    assert(type(file) is str)
    if file.endswith('.gz'):
        mode = mode if mode is not None else 'rt'
        f = gzip.open(file, mode=mode)
    else:
        mode = mode if mode is not None else 'r'
        f = open(file, mode=mode)
    return f

def readFile(file):
    '''
    Reads entire file and returns lines as a list.
    '''
    
    with createFileObject(file) as f:
        contents = f.read().splitlines()
    return contents


def parseEnsemblPepHeader(header, header_prefix='>'):
    '''
    Parse Ensembl Peptide FASTA header. See the README file in the FTP directory (e.g.,
    http://ftp.ensembl.org/pub/release-63/fasta/homo_sapiens/pep/README) for a
    detailed description of the header line format and the file naming conventions.
    See https://uswest.ensembl.org/info/data/ftp/index.html for a description of
    the LOCATION (e.g., 'chromosome:NCBI35:1:904515:910768:1') attribute.
    
    Args:
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix
    
    Returns: dict
      Map of metadata of protein sequence.
      Keys: not all may be present, e.g., the GRCh37.63 release only contains up to 'transcript'
      - id: ENSP ID
      - seqtype: sequence type (pep)
      - status: status of model
          known: can be mapped to species-specific entries UniProt or RefSeq
          novel: cannot be mapped (e.g., predicted based on homology)
      - coord_system: coordinate system (chromosome, supercontig)
      - version: coordinate system version (e.g., GRCh38)
      - name: location name (e.g., 12 - twelfth chromosome)
      - start: start coordinate
      - end: end coordinate
      - strand: strandedness (1, -1)
      - gene: ENSG ID
      - transcript: ENST ID
      - gene_biotype: see https://ensembl.org/info/genome/genebuild/biotypes.html
      - transcript_biotype: see https://ensembl.org/info/genome/genebuild/biotypes.html
      - gene_symbol: HGNC gene symbol
      - description: gene description
      - source: 
      - accession: accession id of sequence from source
    '''
    
    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]
    
    # extract attributes
    p = re.compile(r'(?P<id>[^\s]+)\s+'
                   + r'(?P<seqtype>[^:\s]+):?'
                   + r'(?P<status>[^\s]*)\s+'
                   + r'(?P<coord_system>[^:]+):'
                   + r'(?P<version>[^:]+):'
                   + r'(?P<name>[^:]+):'
                   + r'(?P<start>[^:]+):'
                   + r'(?P<end>[^:]+):'
                   + r'(?P<strand>[^\s]+)\s+'
                   + r'gene:(?P<gene>[^\s]+)\s+'
                   + r'transcript:(?P<transcript>[^\s]+)\s*'
                   + r'(gene_biotype:(?P<gene_biotype>[^\s]+)\s*)*'
                   + r'(transcript_biotype:(?P<transcript_biotype>[^\s]+)\s*)*'
                   + r'(gene_symbol:(?P<gene_symbol>[^\s]+)\s*)*'
                   + r'(description:(?P<description>[^[]+)\s*)*'
                   + r'(\[Source:(?P<source>[^;]+);)*'
                   + r'(Acc:(?P<accession>[^]]+)\])*')
#     # Old set of regex that works specifically with GRCh37.63 release
#     p = re.compile(r'(?P<id>[^\s]+)\s+'
#                    + r'(?P<seqtype>[^:]+):'
#                    + r'(?P<status>[^\s]+)\s+'
#                    + r'(?P<coord_system>[^:]+):'
#                    + r'(?P<version>[^:]+):'
#                    + r'(?P<name>[^:]+):'
#                    + r'(?P<start>[^:]+):'
#                    + r'(?P<end>[^:]+):'
#                    + r'(?P<strand>[^\s]+)\s+'
#                    + r'gene:(?P<gene>[^\s]+)\s+'
#                    + r'transcript:(?P<transcript>[^\s]+)\s*')
    m = p.match(header)
    
    # extract regex match to dict
    data = m.groupdict()
    
    # remove leading/trailing whitespace from each value in dict
    data = {key: value.strip() for key, value in data.items() if value is not None}
    
    return data

def parseUniProtHeader(header, header_prefix='>'):
    '''
    Parse UniProt FASTA header. See https://www.uniprot.org/help/fasta-headers.
    
    Args:
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix
    
    Returns: dict
      Map of metadata of protein sequence.
      Keys: db, id, uniprotName, proteinName, os, ox, gn (may be empty), pe, sv
    '''
    
    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]
    
    # extract gene name if present
    split = re.split(r'GN=(?P<gn>.*)(?=\sPE=)\s+', header)
    m_gn = ''
    if len(split) not in [1,3]:
        raise
    if len(split) > 1:
        m_gn = split[1]
        header = split[0] + split[2]
    
    # extract other variables
    p = re.compile(r'(?P<db>sp|tr)\|(?P<id>[^|]+)\|(?P<uniprotName>\S+)\s+'
                   + r'(?P<proteinName>.*(?=\sOS=)) '
                   + r'OS=(?P<os>.*(?=\sOX=))\s+'
                   + r'OX=(?P<ox>.*(?=\sPE=))\s+'
                   + r'PE=(?P<pe>.*(?=\sSV=))\s+'
                   + r'SV=(?P<sv>.*)')
    m = p.match(header)
    
    # extract regex match to dict
    data = m.groupdict()
    
    # add gene name if present
    data['gn'] = m_gn
    
    # remove leading/trailing whitespace from each value in dict
    data = {key: value.strip() for key, value in data.items()}
    
    return data

def fastaToDF(file, save='', header_prefix='>', headerParser=parseUniProtHeader, **kwargs):
    '''
    Parse FASTA file into pandas DataFrame.
    
    Args
    - file: str or io.IOBase
        Path to FASTA file, or file object. Gzip-compressed files with extension '.gz'
        are accepted.
    - save: str. default=''
        Path to save DataFrame
    - header_prefix: str. default='>'
        FASTA header line prefix
    - headerParser: function. default = parseUniProtHeader
        Function to parse header line into dict
    - **kwargs
        Additional keyword arguments to pass to pandas.DataFrame.to_csv() if saving DataFrame.
    
    Returns: pandas.DataFrame
      Rows: protein / DNA entries
      Columns: data about entries. Always includes 'seq' (sequence) column.
    '''
    
    closeAtEnd = False
    if isinstance(file, str):
        f = createFileObject(file)
        closeAtEnd = True
    elif isinstance(file, io.IOBase):
        f = file
    else:
        raise ValueError('`file` must be a string or file object')
    
    entries = []
    entry = {'seq': ''}
    while True:
        # read 1 line at a time to avoid memory overflow
        line = f.readline()

        if line == '':
            entries.append(entry)
            break

        if line.startswith(header_prefix):
            # add previous entry to running list of entries
            if entry['seq'] is not '':
                # add sequence to entry
                entries.append(entry)
            
            # parse new entry
            entry = headerParser(line)
            entry['seq'] = ''
        else:
            entry['seq'] += line.strip()
    
    if closeAtEnd:
        f.close()
    
    # construct pandas DataFrame from list of dicts
    df = pd.DataFrame(entries)
    if save is not '':
        df.to_csv(save, **kwargs)
    
    return df

def keepValidDisorder(disorder, proteome, makeCopy=False,
                      dropIfNotInProteome=False, dropIfFlippedStartEnd=False,
                      startCol='start', endCol='end', disorderIdCol='d2p2_id',
                      proteomeIdCol='d2p2_id', seqCol='seq', widthCol='width'):
    '''
    Remove invalid disorder entries.
    
    Assumptions
    - Disordered regions are 1-indexed
    
    Args
    - disorder: pandas.DataFrame
        Required columns: `startCol`, `endCol`, `disorderIdCol`
    - proteome: pandas.DataFrame
        Requied columns: `proteomeIdCol`, one of `seqCol` or `widthCol`.
        `proteomeIdCol` should contain the same type of ID as `disorderIdCol`
    - makeCopy: bool. default=False
        Make a deep copy of the `disorder` DataFrame that is then returned.
        If False, this function may modify the `disorder` DataFrame in-place.
    - dropIfNotInProteome: bool. default=False
        Only keep proteins in the proteome
    - dropIfFlippedStartEnd: bool. default=False
        Drop entries (disordered regions) where start > end
    - startCol, endCol, disorderIdCol, proteomeIdCol, seqCol, widthCol: str
        Relevant columns in `disorder` and `proteome`
    
    Returns: pandas.DataFrame
      Guarantees:
      1. startCol <= endCol
      2. startCol, endCol > 0
      3. startCol, endCol <= width
           (where width information can be be determined from the proteome)
      If `makeCopy` is False, the returned DataFrame may be a view (slice) into `disorder`.
    '''
    
    if makeCopy:
        disorder = disorder.copy()
  
    valid = (disorder[startCol].values > 0) & (disorder[endCol].values > 0)
    flipped = disorder[startCol].values > disorder[endCol].values
    nFlipped = sum(flipped)
    
    if proteome is not None:
        if dropIfNotInProteome:
            valid = valid & disorder[disorderIdCol].isin(set(proteome[proteomeIdCol]))
        
        # add width column to proteome
        if widthCol not in proteome.columns:
            proteome[widthCol] = list(map(len, proteome[seqCol]))
        
        # create temporary data frame with width information
        if sum(proteome[proteomeIdCol].duplicated() > 0):
            print('Dropping duplicate entries with the same {} value in the proteome.'.format(proteomeIdCol))
            proteome.drop_duplicates(subset=proteomeIdCol, inplace=True)
        tmp = disorder.merge(proteome[[proteomeIdCol, widthCol]], how='left',
                             left_on=disorderIdCol, right_on=proteomeIdCol)
        tmp.loc[tmp[widthCol].isnull(), widthCol] = np.Inf
        assert(tmp.shape[0] == disorder.shape[0])
        
        # disorder region must be within length of the protein
        valid = valid & \
          (disorder[startCol].values <= tmp[widthCol].values) & \
          (disorder[endCol].values <= tmp[widthCol].values)
    
    if dropIfFlippedStartEnd:
        valid = valid & ~flipped
    elif nFlipped > 0:
        start_orig = disorder.loc[flipped, 'start'].values
        end_orig = disorder.loc[flipped, 'end'].values
        disorder.loc[flipped, 'start'] = end_orig
        disorder.loc[flipped, 'end'] = start_orig
        print("Flipped {} entries: {}".format(nFlipped, np.where(flipped)[0]))
    
    nInvalid = sum(~valid)
    print("Removed {} ({:.2%}) entries: {}".format(nInvalid, nInvalid/disorder.shape[0], np.where(~valid)[0]))
    
    return disorder[valid]

def parseQuickGOJSON(results):
    '''
    Parse JSON body returned by QuickGO 'search' API into pandas.DataFrame with same columns
    as returned by the 'download' API
    
    Arg: list of dict
    
    Returns: pandas.DataFrame
      Columns will have the same names and order as returned by the QuickGO 'download' API
    '''
    
    goAspect_map = {'cellular_component': 'C', 'biological_process': 'P', 'molecular_function': 'F'}
    
    colNames_map = OrderedDict([
        # entries in order by 'search' API method (see getQuickGO())
        ('geneProductDb', 'GENE PRODUCT DB'),
        ('geneProductId', 'GENE PRODUCT ID'),
        ('symbol', 'SYMBOL'),
        ('qualifier', 'QUALIFIER'),
        ('goId', 'GO TERM'),
        ('goAspect', 'GO ASPECT'),
        ('evidenceCode', 'ECO ID'),
        ('goEvidence', 'GO EVIDENCE CODE'),
        ('reference', 'REFERENCE'),
        ('withFrom', 'WITH/FROM'),
        ('taxonId', 'TAXON ID'),
        ('assignedBy', 'ASSIGNED BY'),
        ('extensions', 'ANNOTATION EXTENSION'),
        ('date', 'DATE'),
        
        # additional entries retrieveable only via 'download' method
        ('goName', 'GO NAME'),
        ('synonyms', 'SYNONYMS')
    ])
    
    for i in range(len(results)):
        # extract 'withFrom' dictionaries into strings
        withFrom_list = []
        if results[i]['withFrom'] is not None:
            for j in range(len(results[i]['withFrom'])):
                for k in range(len(results[i]['withFrom'][j]['connectedXrefs'])):
                    withFrom_list.append(':'.join([results[i]['withFrom'][j]['connectedXrefs'][0]['db'],
                                                   results[i]['withFrom'][j]['connectedXrefs'][0]['id']]))
        results[i]['withFrom'] = '|'.join(withFrom_list)
        
        # extract 'geneProductId' into 'GENE PRODUCT DB' and 'GENE PRODUCT ID'
        results[i]['geneProductDb'], results[i]['geneProductId'] = results[i]['geneProductId'].split(':')
        
        # convert 'goAspect' to single-character symbol
        results[i]['goAspect'] = goAspect_map[results[i]['goAspect']]
    
    # rename and reorder columns
    df = pd.DataFrame(results).rename(mapper=colNames_map, axis='columns')
    df = df[list(colNames_map.values())]
    
    # convert 'DATE' to date format; convert None values to np.nan
    df['DATE'] = pd.to_datetime(df['DATE'], yearfirst=True, format='%Y%m%d')
    df.loc[df['ANNOTATION EXTENSION'].isnull(), 'ANNOTATION EXTENSION'] = np.nan
    return df

def getQuickGO(goIds, useDefaults=True, method='opt', parseResponse=True,
               attempts=5, sleep=0.5, nThreads=1, **kwargs):
    '''
    Download annotations from QuickGO. See https://www.ebi.ac.uk/QuickGO/api/index.html.
    
    Args
    - goIds: str
        Comma-separated string of GO IDs (e.g., 'GO:0016592' for mediator complex)
    - useDefaults: bool. default=True
        Apply the following 'default' filters:
        - taxonId: 9606 (Homo sapiens)
        - geneProductType: 'protein'
        - geneProductSubset: 'Swiss-Prot' (only reviewed UniProtKB entries)
        - proteome: 'gcrpCan' (Gene Centric Reference Proteome Canonical)
        - goUsage: 'descendants'
        - goUsageRelationships: 'is_a,part_of,occurs_in' (excludes 'regulates')
        - limit: 100 (maximum limit per page for 'search' API)
        - downloadLimit: 50000 (maximum limit for 'download' API)
    - method: str. default='opt'
        'search': Use QuickGO 'search' API, which returns a JSON response body.
          Required for large (>50000 returned entries) query
        'download': Use QuickGO 'download' API, which returns a text (TSV) response body.
          Maximum 50000 returned entries.
        'opt': Try 'search'. If the number of hits is < 50000, switch to 'download'.
    - parseResponse: bool. default=True
        True: parse response into pandas.DataFrame
        False: return respose bodies
    - attempts: int. default=5
        Number of attempts to retry a request upon error
    - sleep: float. default=0.5
        Seconds to sleep for in between each attempt
    - nThreads: int. default=1
        Number of threads to use. Error handling is not implemented if nThreads > 1.
    - **kwargs
        Additional parameters to pass to send in the body of the HTTP request
    
    Returns: pandas.DataFrame, list of requests.Response, or None
      `parseResponse` == True --> pandas.DataFrame
      `parseResponse` == False --> list of requests.Response
      Error in request after `attempt` attempts --> None
    
    Notes
    - The following parameters do not seem to have any effect using the download API:
      includeFields, limit
    '''
    
    # validate arguments
    assert(type(attempts) is int and attempts > 0)
    assert(not(nThreads > 1 and method == 'download'))
    assert(method in ['search', 'download', 'opt'])
    
    if method in ['search', 'opt']:
        url = 'https://www.ebi.ac.uk/QuickGO/services/annotation/search'
        headers = {'Accept': 'application/json'}
    else:
        url = 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch'
        headers = {'Accept': 'text/tsv'}
    
    # set HTTP request body parameters
    defaults = {
        'taxonId': 9606,
        'geneProductType': 'protein',
        'geneProductSubset': 'Swiss-Prot',
        'proteome': 'gcrpCan',
        'goUsage': 'descendants',
        'goUsageRelationships': 'is_a,part_of,occurs_in',
        'includeFields': 'goName,synonyms',
        'limit': 100,
        'downloadLimit': 50000
    }
    params = {'goId': goIds}
    if useDefaults:
        params.update(defaults)
    params.update(kwargs)
    
    (r, attempts_remaining) = getQuickGO_helper(url, params, headers, method, attempts, sleep)
    if attempts_remaining > 0:
        allResponses = [r]
    else:
        return None
    
    if method == 'download':
        if parseResponse:
            df = pd.read_table(io.StringIO(r.text))
            if 'DATE' in df.columns:
                df['DATE'] = pd.to_datetime(df['DATE'], yearfirst=True, format='%Y%m%d')
            return df
        else:
            return allResponses
    
    response_body = r.json()
    totalPages = response_body['pageInfo']['total']
    params.update({'page': response_body['pageInfo']['current']})
    if totalPages > 1 and params['page'] < totalPages:
        if method == 'opt' and response_body['numberOfHits'] < 50000:
            print('Switching to \'download\' API...')
            return getQuickGO(goIds, useDefaults, 'download', parseResponse, **kwargs)
        
        params['page'] += 1
        if nThreads > 1:
            allResponses.extend(getQuickGO_mt(url, params, headers, attempts, sleep, nThreads,
                                              pages=range(params['page'],totalPages+1)))
        else:
            while True:
                (r, attempts_remaining) = getQuickGO_helper(url, params, headers, attempts, sleep)
                if attempts_remaining > 0:
                    allResponses.append(r)
                else:
                    print('Skipping page {}'.format(params['page']))
                if params['page'] >= totalPages: 
                    break
                params.update({'page': params['page'] + 1})
    
    if parseResponse:
        results = list(itertools.chain.from_iterable([r.json()['results'] for r in allResponses]))
        return parseQuickGOJSON(results)
    else:
        return allResponses

def getQuickGO_helper(url, params, headers, method, attempts=5, sleep=0.5):
    '''
    Returns: (request.Response, int)
      HTTP response and unused attempts
    '''
    
    while attempts > 0:
        try:
            r = requests.get(url, params=params, headers=headers)
            if not r.ok:
                r.raise_for_status()
            if method is not 'download':
                response_body = r.json()
                print('numberOfHits: {}'.format(response_body['numberOfHits']),
                      response_body['pageInfo'], sep=', ')
            break
        except Exception as err:
            print(err)
            if (attempts > 1):
                print('Attempts remaining: {}'.format(attempts-1))
                time.sleep(sleep)
            else:
                print('All attempts exhausted')
        attempts -= 1
    
    return (r, attempts)

def getQuickGO_mt(url, params, headers, attempts, sleep, nThreads, pages):
    '''
    Multithreaded implementation of getQuickGO(method='search') without error handling
    '''
    
    pool = Pool(nThreads)
    print("Using {:d} threads...".format(pool._processes))
    sleep = max(sleep, 0.1*nThreads)
    
    allResponses = []
    for i in pages:
        params.update({'page': i})
        allResponses.append(pool.apply_async(getQuickGO_helper,
                                             (url, params, headers, 'search', attempts, sleep)))
    pool.close()
    pool.join()
    
    for i in range(len(allResponses)):
        allResponses[i] = allResponses[i].get()
    
    return(allResponses)

def processQuickGO(df, save='', **kwargs):
    '''
    Process annotations returned by QuickGO such that
    1. there is only one entry per gene name
    2. there is only one entry per UniProtKB Accession ID
    
    Args:
    - df: pandas.DataFrame
    - save: str. default=''
        Path to save DataFrame
    - **kwargs
        Additional keyword arguments to pass to pandas.DataFrame.to_csv() if saving DataFrame.
    
    Returns: pandas.DataFrame
    '''
    
    df = df[~df.duplicated(subset=['GENE PRODUCT ID', 'SYMBOL']).values]
    assert(~any(df.duplicated(subset='GENE PRODUCT ID')))
    assert(~any(df.duplicated(subset='SYMBOL')))
    
    if save is not '':
        df.to_csv(save, **kwargs)
    
    return df

def matchDups(df, ref, dup_col='seq', df_id_col='ensembl_peptide_id', ref_id_col='id', ref_seq_col='seq'):
    '''
    Find duplicated values of df[dup_col] whose corresponding (matching df[df_id_col] and ref[ref_id_col])
    values in ref[ref_seq_col] are / are not all duplicates.
    
    Args
    - df: pandas.DataFrame
    - ref: pandas.DataFrame
    - dup_col: str. default='seq'
    - df_id_col: str. default='ensembl_peptide_id'
    - ref_id_col: str. default='id'
    - ref_seq_col: str. default='seq'
    
    Returns: (list of df[dup_col] values whose corresponding ref[ref_seq_col] entries are duplicates
              list of df[dup_col] values whose corresponding ref[ref_seq_col] entries are not all duplicates)
    '''
    
    dupVals = set(df.loc[df[dup_col].duplicated(keep=False), dup_col])
    unMatched = []
    matched = []
    for dupVal in dupVals:
        if pd.isnull(dupVal):
            ids = set(df.loc[df[dup_col].isnull(), df_id_col])
        else:
            ids = set(df.loc[df[dup_col] == dupVal, df_id_col])
        if ~ref.loc[ref[ref_id_col].isin(ids), ref_seq_col].duplicated(keep=False).all():
            unMatched.append(dupVal)
        else:
            matched.append(dupVal)
    return (matched, unMatched)

def removeDupsByRank(df, dup_col='ensembl_peptide_id', rank_col='transcript_appris',
                     rank_map=None, verbose=True):
    '''
    Remove duplicates by ranking.
    
    Args
    - df: pandas.DataFrame
    - dup_cols: str or list of str. default='ensembl_peptide_id'
    - rank_col: str. default='transcript_appris'
    - rank_map: dict. default=None
    - verbose: bool. default=True
        Print out dropped rows.
    
    Returns: pandas.DataFrame
    '''
    
    assert(rank_map is None or isinstance(rank_map, dict))
    if isinstance(rank_map, dict):
        ranks = np.array([rank_map[val] for val in df[rank_col]])
    else:
        ranks = df[rank_col].values
    
    toDropAll = []
    dupVals = set(df.loc[df[dup_col].duplicated(keep=False), dup_col])
    for dupVal in dupVals:
        if pd.isnull(dupVal):
            dupRows = np.where(df[dup_col].isnull())[0]
        else:
            dupRows = np.where(df[dup_col] == dupVal)[0]
        minRank = min(ranks[dupRows])
        toDrop = [row for row in dupRows if ranks[row] > minRank]
        toDropAll.extend(toDrop)
        if verbose and len(toDrop) > 0:
            print('Dropping rows where column {} == {}: {}'.format(dup_col, dupVal, toDrop))
    
    return df.drop(df.index[toDropAll], axis='index')

def removeDups(df, ref, dup_col='uniprot_id', metadata_cols=['gn', 'hgnc_symbol'], verbose='low'):
    '''
    Form 1:1 mapping between UniProt, Ensembl, and HGNC IDs
    
    Args
    - df: pandas.DataFrame
    - dup_col: str. default='uniprot_id'
    - metadata_cols: list of str. default=['gn', 'hgnc_symbol']
    - verbose: str, bool, or None. default='low'
        3 levels of verbosity: True/'high', False/'low', None/'off'
    
    Returns: pandas.DataFrame
    '''
    
    # argument validation
    assert(verbose in [True, False, None, 'high', 'low', 'off'])
    assert(dup_col in df.columns)
    assert(all([col in df.columns for col in metadata_cols]))
    
    # process arguments
    tmp = {True: 'high', False: 'low', None: 'off'}
    if verbose in [True, False, None]:
        verbose = tmp[verbose]
    
    df = df[df['db'] == 'sp']
    df = df[~(((df['gn'].isnull()) | (df['gn'] == '')) &
              ((df['hgnc_symbol'].isnull()) | (df['hgnc_symbol'] == '')))]
    
    matched, unmatched = matchDups(df, ref, dup_col=dup_col)
    if verbose == 'high':
        print('unmatched: {}'.format(unmatched))
    
    # matching sequences
    # 1. entry with more metadata
    # 2. arbitrarily choose entry
    toDropAll = []
    for dupVal in matched:
        toDrop = []
        if pd.isnull(dupVal):
            dupRows_bool = df[dup_col].isnull()
        else:
            dupRows_bool = df[dup_col] == dupVal
        dupRows = np.where(dupRows_bool)[0]
        
        quantityMetadata = [sum(~(df.loc[row, metadata_cols].isnull())) for row in dupRows]
        lackingMetadata = [dupRows[i] for i in range(len(dupRows)) if quantityMetadata[i] < max(quantityMetadata)]
        if verbose == 'high' and len(lackingMetadata) > 0:
            print('Dropping rows where column {} == {} and lacking metadata: {}'.format(dup_col, dupVal, lackingMetadata))
        toDrop.extend(lackingMetadata)
        dupRows = list(set(dupRows) - set(toDrop))
        
        if len(dupRows) > 1:
            toDrop.extend(dupRows[1:])
        toDropAll.extend(toDrop)
        
        if verbose == 'high' and len(toDrop) > 0:
            print('Dropping rows where column {} == {}: {}'.format(dup_col, dupVal, toDrop))

    df = df.drop(df.index[toDropAll], axis='index')
    
    # unmatched sequences
    # fill in null HGNC symbols with UniProt gene names
    rowsNullSymbol = df['hgnc_symbol'].isnull()
    df.loc[rowsNullSymbol, 'hgnc_symbol'] = df.loc[rowsNullSymbol, 'gn']
    if verbose == 'high':
        print('Filling in rows with null HGNC symbols with UniProt gene names: {}'.format(np.where(rowsNullSymbol)))
    
    # check where 3 out of 4 match: uniprot_id, ensembl_peptide_id, hgnc_symbol, seq
    idCols = ['uniprot_id', 'ensembl_peptide_id', 'hgnc_symbol', 'seq']
    combinations = itertools.combinations(idCols, len(idCols)-1)
    for comb in combinations:
        df = df.drop_duplicates(subset=comb, keep='first')
    
    # manual curation
    # - keep hgnc_symbol:UGT2A1, drop hgnc_symbol:UGT2A2
    #     "According to HGNC, UGT2A1 and UGT2A2 are 2 separate genes. However, they share
    #     common exons at the C-terminus, suggesting that UGT2A1 and UGT2A2 are different
    #     isoforms encoded by the same locus." (https://www.uniprot.org/uniprot/Q9Y4X1)
    # - edit gn:NUDT4B
    #      NUDT4B is a different gene (located on 2 separate chromosomes) than NUDT4. Appears
    #      that BioMart failed to recognize the difference, even though Ensembl, UniProt, and
    #      HGNC all have separate entries. 
    # - keep ensembl_peptide_id:ENSP00000361746, drop ensembl_peptide_id:ENSP00000424176
    #      ENSP00000361746 corresponds to the canonical sequence O95925-1
    #      (https://www.uniprot.org/uniprot/O95925)
    # - keep ensembl_peptide_id: ENSP00000337450, drop ensembl_peptide_id:ENSP00000480571
    #      ENSP00000337450 corresponds to the canonical sequence P24462-1
    #      (https://www.uniprot.org/uniprot/P24462)
    # - keep ensembl_peptide_id: ENSP00000250405, drop ensembl_peptide_id:ENSP00000451320
    #      ENSP00000250405 corresponds to the canonical sequence Q92843-1
    #      (https://www.uniprot.org/uniprot/Q92843)
    # - keep ensembl_peptide_id: ENSP00000335385, drop ensembl_peptide_id:ENSP00000364777
    #      ENSP00000335385 corresponds to the canonical sequence Q96JG8-1
    #      (https://www.uniprot.org/uniprot/Q96JG8). However, HGNC and Ensembl both have 2
    #      different gene names: MAGED4 and MAGED4B
    # - keep ensembl_peptide_id: ENSP00000400870, drop ensembl_peptide_id:ENSP00000409924
    #      ENSP00000400870 corresponds to the canonical sequence Q9BTE6-1
    #      (https://www.uniprot.org/uniprot/Q9BTE6)
    df.loc[df['gn'] == 'NUDT4B', ['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'hgnc_symbol']] = \
      ['ENSG00000177144', 'ENST00000322209', 'ENSP00000492425', 'NUDT4B']
    toDrop_ensp = set(['ENSP00000424176', 'ENSP00000480571', 'ENSP00000451320', 'ENSP00000364777', 'ENSP00000409924'])
    df = df.drop(df.index[(df['hgnc_symbol'] == 'UGT2A2') | (df['ensembl_peptide_id'].isin(toDrop_ensp))], axis='index')
    
    if verbose in ['high', 'low']:
        print('Number of remaining duplicates:')
        print('UniProt IDs: {}'.format(sum(df['uniprot_id'].duplicated(keep=False))))
        print('Ensembl Peptide IDs: {}'.format(sum(df['ensembl_peptide_id'].duplicated(keep=False))))
        print('HGNC symbols: {}'.format(sum(df['hgnc_symbol'].duplicated(keep=False))))
        print('Sequences: {}'.format(sum(df['seq'].duplicated(keep=False))))
    return df

def manualMap(idMap, uniprotProteome, ensemblProteome):
    # Args
    # - idMap: pandas.DataFrame
    # - uniprotProteome: pandas.DataFrame
    # - ensemblProteome: pandas.DataFrame
    # 
    # Notes
    # - MYC: UCSC knownCanonical (GRCh38) lists ENSP00000479618 (ENST00000613283) as the canonical transcript, whereas
    #     Uniprot lists ENSP00000367207 (ENST00000377970) as the canonical sequence. Only UniProt's canonical sequence
    #     is in the D2P2 proteome.
    # - MED19: UCSC knownCanonical (GRCh38) and APPRIS (in the Ensembl online browser) list ENSP00000416227 (ENST00000431606)
    #     as the canonical transcript, whereas Uniprot lists ENSP00000337340 (ENST00000337672) as the canonical sequence.
    #     Both are in the D2P2 proteome as separate entries.
    
    # MYC: UCSC knownCanonical (GRCh38) lists ENSP00000479618 (ENST00000613283) as the canonical transcript, whereas
    #   Uniprot lists ENSP00000367207 (ENST00000377970) as the canonical sequence. Only UniProt's canonical sequence
    #   is in the D2P2 proteome.
    idMap.loc[idMap['hgnc_symbol'] == 'MYC', ['ensembl_transcript_id', 'ensembl_peptide_id']] = ['ENST00000377970', 'ENSP00000367207']
    
    # MED19: Use UniProt entry
    #   UCSC knownCanonical (GRCh38) and APPRIS (in the Ensembl online browser) list ENSP00000416227 (ENST00000431606)
    #   as the canonical transcript, whereas Uniprot lists ENSP00000337340 (ENST00000337672) as the canonical sequence.
    #   Both are in the D2P2 proteome as separate entries.
    med19 = ensemblProteome.loc[ensemblProteome['id'] == 'ENSP00000337340',
                                ['id', 'gene_symbol', 'seq', 'gene', 'transcript']] \
        .squeeze() \
        .rename(ensemblProteomeColumnMap)
    med19 = med19.append(uniprotProteome.loc[uniprotProteome['gn'] == 'MED19', ['id', 'proteinName']] \
        .squeeze() \
        .rename(uniprotProteomeColumnMap))
    idMap = idMap.append(med19, ignore_index=True)
    
    # MED29: Use Ensembl/knownCanonical/D2P2 entry
    # - UniProt: ENST00000599213, ENSP00000471802, ENSG00000063322 [Q9NX70-1]
    # - Biomart: UniProt + 1 other
    # - Ensembl / UCSC knownCanonical: ENST00000315588, ENSP00000314343, ENSG00000063322 [B4DUA7]
    # - D2P2: Ensembl / knownCanonical + 1 other
    med29 = ensemblProteome.loc[ensemblProteome['id'] == 'ENSP00000314343',
                                ['id', 'gene_symbol', 'seq', 'gene', 'transcript']] \
        .squeeze() \
        .rename(ensemblProteomeColumnMap)
    med29 = med29.append(uniprotProteome.loc[uniprotProteome['gn'] == 'MED29', ['id', 'proteinName']] \
        .squeeze() \
        .rename(uniprotProteomeColumnMap))
    idMap = idMap.append(med29, ignore_index=True)
    return idMap

def getUniProtFasta(id, parseToDF=True):
    '''
    Args
    - id: str
        UniProt Accession ID (e.g., P12345)
    - parseToDF: bool. default=True
        Parse FASTA to pandas.DataFrame using fastaToDF()
    
    Returns: str or pandas.DataFrame
    '''
    
    url = 'https://www.uniprot.org/uniprot/{}.fasta'.format(id)
    headers = {'accept': 'text/html'}
    r = requests.get(url=url, headers=headers)
    if not r.ok:
        r.raise_for_status()
    if parseToDF:
        return fastaToDF(io.StringIO(r.text), headerParser=parseUniProtHeader)
    else:
        return r.text

def generateProteome(ids, ref=None, idCol='id', getFasta=getUniProtFasta, nThreads=1):
    '''
    Create proteome from protein IDs by subsetting local reference proteome and retrieving
    additional protein data (online).
    
    Args
    - ids: list of str
        Type of IDs corresponds to `idCol`
    - ref: pandas.DataFrame. default=None
        Reference proteome for local lookup of sequences, metadata. Must contain `idCol`.
    - idCol: str. default='id'
        Column in `ref` to match `ids`
    - getFasta: function. default=getUniProtFasta
        Function to retrieve (download) protein data if not in `ref`. Should return entries (rows)
        in the same format as in `ref`
    - nThreads: int
        Number of threads to use with `getFasta` function
    
    Returns: pandas.DataFrame
      Columns will depend on `ref` and `getFasta`.
    '''
    
    ids = set(ids)
    local = set()
    df_local = None
    df_download = None
    
    if ref is not None:
        local = set([id for id in ids if id in set(ref[idCol])])
        df_local = ref[ref[idCol].isin(local)]
    toDownload = list(ids - local)
    
    if len(toDownload) > 0 and getFasta is not None:
        if nThreads > 1:
            pool = Pool(nThreads)
            print("Using {:d} threads...".format(pool._processes))
            df_download = []
            for i in range(len(toDownload)):
                df_download.append(pool.apply_async(getFasta, (toDownload[i],)))
            pool.close()
            pool.join()
            df_download = pd.concat([df_download[i].get() for i in range(len(df_download))])
        else:
            df_download = pd.concat([getFasta(id) for id in toDownload])
    return pd.concat([df_download, df_local])