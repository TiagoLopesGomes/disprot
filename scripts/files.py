import platform, os

compressionExt = {'gzip': '.gz', 'bzip2': '.bz2', 'zip': '.zip', 'xz': '.xz'}
save_kwargs = {'sep': '\t', 'index': False, 'compression': 'gzip'}

def getProjectDir():
    '''
    Returns project directory, depending on the operating system:
    - Linux: ~/projects/disprot
    - Windows: Documents/Projects/disprot
    '''
    
    if platform.system() == 'Linux':
        return os.path.join(os.environ['HOME'], 'projects', 'disprot')
    elif platform.system() == 'Windows':
        # See https://stackoverflow.com/a/30924555
        import ctypes.wintypes
        CSIDL_PERSONAL = 5       # My Documents
        SHGFP_TYPE_CURRENT = 0   # Get current, not default value
        buf = ctypes.create_unicode_buffer(ctypes.wintypes.MAX_PATH)
        ctypes.windll.shell32.SHGetFolderPathW(None, CSIDL_PERSONAL, None, SHGFP_TYPE_CURRENT, buf)
        my_documents = buf.value
        return os.path.join(my_documents, 'Projects', 'disprot')
    else:
        raise Exception('OS platform not recognized.')

def getPaths(projectDir=None, ext='.gz', makedirs=True):
    '''
    Returns dictionary mapping names to file paths.
    
    Args
    - projectDir: str. default=None
        Project directory
    - ext: str. default='.gz'
        (Compression) extension to add to processed data files
    - makedirs: bool. default=True
        Recursively create all directories.
    
    Returns: dict (str -> str)
    '''
    
    if projectDir is None:
        projectDir=getProjectDir()
    
    dirData = os.path.join(projectDir, 'data')
    dirDataRaw = os.path.join(dirData, 'ref_raw')
    dirDataProc = os.path.join(dirData, 'ref_proc')
    dirDataAux = os.path.join(dirData, 'data_aux')
    dirDataGeneLists = os.path.join(dirData, 'gene_lists')
    dirDataGeneListsRaw = os.path.join(dirData, 'gene_lists', 'raw')
    dirDataML = os.path.join(dirData, 'ml_datasets')
    dirResults = os.path.join(projectDir, 'results')
    dirResults_seqPlots = os.path.join(dirResults, 'seqPlots')
    dirResults_mlPlots = os.path.join(dirResults, 'mlPlots')
    dirScripts = os.path.join(projectDir, 'scripts')
    
    paths = {
        # directories
        'dirProject': projectDir,
        'dirDataRaw': dirDataRaw,
        'dirDataProc': dirDataProc,
        'dirDataAux': dirDataAux,
        'dirDataGeneLists': dirDataGeneLists,
        'dirDataGeneListsRaw': dirDataGeneListsRaw,
        'dirDataML': dirDataML,
        'dirResults': dirResults,
        'dirResults_seqPlots': dirResults_seqPlots,
        'dirResults_mlPlots': dirResults_mlPlots,
        'dirScripts': dirScripts,
        
        # raw data
        'uniprotProteome_raw': os.path.join(dirDataRaw, 'uniprot', 'UP000005640_9606.fasta.gz'),
        'uniprotIdMap_raw': os.path.join(dirDataRaw, 'uniprot', 'UP000005640_9606.idmapping.gz'),
        'uniprotIdMapFull_raw': os.path.join(dirDataRaw, 'uniprot', 'HUMAN_9606_idmapping.dat.gz'),
        'ensemblProteome63_raw': os.path.join(dirDataRaw, 'ensembl', 'Homo_sapiens.GRCh37.63.pep.all.fa.gz'),
        'ensemblProteome_raw': os.path.join(dirDataRaw, 'ensembl', 'Homo_sapiens.GRCh38.pep.all.fa.gz'),
        'ensemblBiomart_raw': os.path.join(dirDataRaw, 'ensembl', 'biomart.tsv.gz'),
        'd2p2Proteome_raw': os.path.join(dirDataRaw, 'd2p2', 'protein', 'genomes.protein.gz'),
        'd2p2Disorder_raw': os.path.join(dirDataRaw, 'd2p2', 'disorder', 'all.disrange.gz'),
        'd2p2Disorder_vlxt_raw': os.path.join(dirDataRaw, 'd2p2', 'disorder', 'vlxt.disrange.gz'),
        'd2p2Disorder_vsl2b_raw': os.path.join(dirDataRaw, 'd2p2', 'disorder', 'vsl2b.disrange.gz'),
        'd2p2IdMap_raw': os.path.join(dirDataRaw, 'd2p2', 'd2p2_protein_to_uniprot.tsv.bz2'),
        'hgncGeneGroups': os.path.join(dirDataRaw, 'hgnc', 'gene_groups.tsv.gz'),
        'ucscKnownCanonical37_raw': os.path.join(dirDataRaw, 'ucsc', 'knownCanonical_GRCh37.txt.gz'),
        'ucscKnownCanonical38_raw': os.path.join(dirDataRaw, 'ucsc', 'knownCanonical_GRCh38.txt.gz'),
        'disprot_raw': os.path.join(dirDataRaw, 'disprot.csv.gz'),
        'hippie_raw': os.path.join(dirDataRaw, 'hippie_current.txt.gz'),
        
        # processed data
        # - all have headers
        'uniprotProteome': os.path.join(dirDataProc, 'uniprot', 'UP000005640_9606.tsv' + ext),
        'ensemblProteome63': os.path.join(dirDataProc, 'ensembl', 'Homo_sapiens.GRCh37.63.pep.all.tsv' + ext),
        'ensemblProteome': os.path.join(dirDataProc, 'ensembl', 'Homo_sapiens.GRCh38.pep.all.tsv' + ext),
        'd2p2ProteomeHuman': os.path.join(dirDataProc, 'd2p2', 'human.protein.tsv' + ext),
        'd2p2DisorderHuman': os.path.join(dirDataProc, 'd2p2', 'human.disrange.tsv' + ext),
        'd2p2DisorderHuman_vlxt': os.path.join(dirDataProc, 'd2p2', 'human.disrange_vlxt.tsv' + ext),
        'd2p2DisorderHuman_vsl2b': os.path.join(dirDataProc, 'd2p2', 'human.disrange_vsl2b.tsv' + ext),
        'ucscKnownCanonical37': os.path.join(dirDataProc, 'ucsc', 'knownCanonical_GRCh37.tsv' + ext),
        'ucscKnownCanonical38': os.path.join(dirDataProc, 'ucsc', 'knownCanonical_GRCh38.tsv' + ext),
        'idMap': os.path.join(dirDataProc, 'idMap.tsv' + ext),
        
        # raw gene lists
        # - mediatorTFs: transcription factor-Mediator subunit interactions
        #     columns: MED1, MED12, MED14, MED15, MED16, MED17, MED19, MED21,
        #       MED23, MED25, MED26, MED29, CDK8
        #     source: Boija, A. et al. Cell 175, (2018). Table S1.
        # - QuickGO: selection of transcription-related GO term annotations
        #     columns: GENE PRODUCT DB, GENE PRODUCT ID, SYMBOL, QUALIFIER, GO TERM,
        #       GO ASPECT, ECO ID, GO EVIDENCE CODE, REFERENCE, WITH/FROM, TAXON ID,
        #       ASSIGNED BY, ANNOTATION EXTENSION, DATE, (optional) GO NAME,
        #       (optional) SYNONYMS
        #     source: processData.ipynb > Download Data section
        # - interest: genes of interest
        #     each line contains a single gene symbol
        'gl_mediatorTFs_raw': os.path.join(dirDataGeneListsRaw, 'mediatorTFs.tsv.gz'),
        'gl_QuickGO_raw': os.path.join(dirDataGeneListsRaw, 'QuickGO.tsv.gz'),
        'gl_interest_raw': os.path.join(dirDataGeneListsRaw, 'gl_interest.txt'),
        
        # gene lists: processed gene lists
        #   columns: uniprot_id, hgnc_symbol, d2p2_id
        # 
        #   Transcription factor-Mediator subunit interactions
        #   - mediatorTFs: see raw gene lists
        #       list (not table) of UniProt IDs, including human and non-human proteins
        #   - mediatorTFs_human: human subset of mediatorTFs
        #   - MED1: proteins that interact with MED1
        #   - MED15: proteins that interact with MED15
        # 
        #   GO term annotations
        #   - QuickGO: see raw gene lists
        #   - TFs: GO:003700 (DNA-binding transcription factor activity)
        #   - activator: GO:0001228 (DNA-binding transcription activator activity, RNA polymerase II-specific)
        #   - repressor: GO:0001227 (DNA-binding transcription repressor activity, RNA polymerase II-specific)
        # 
        #   HGNC gene groups subsets
        #   - GTFs: general transcription factors
        #   - medComplex: mediator complex
        #   - POLR: RNA polymerase subunits
        'gl_mediatorTFs': os.path.join(dirDataGeneLists, 'mediatorTFs.txt' + ext),
        'gl_mediatorTFs_human': os.path.join(dirDataGeneLists, 'mediatorTFs_human.tsv' + ext),
        'gl_MED1': os.path.join(dirDataGeneLists, 'gl_MED1.tsv' + ext),
        'gl_MED15': os.path.join(dirDataGeneLists, 'gl_MED15.tsv' + ext),
        'gl_QuickGO': os.path.join(dirDataGeneLists, 'QuickGO.tsv' + ext),
        'gl_TFs': os.path.join(dirDataGeneLists, 'gl_TFs.tsv' + ext),
        'gl_activator': os.path.join(dirDataGeneLists, 'gl_activator.tsv' + ext),
        'gl_repressor': os.path.join(dirDataGeneLists, 'gl_repressor.tsv' + ext),
        'gl_GTFs': os.path.join(dirDataGeneLists, 'gl_GTFs.tsv' + ext),
        'gl_medComplex': os.path.join(dirDataGeneLists, 'gl_medComplex.tsv' + ext),
        'gl_POLR': os.path.join(dirDataGeneLists, 'gl_POLR.tsv' + ext),
        'gl_interest': os.path.join(dirDataGeneLists, 'gl_interest.tsv' + ext),
        'gl_random': os.path.join(dirDataGeneLists, 'gl_random.tsv' + ext),
        
        # miscellaneous
        'standardChromosomes': os.path.join(dirDataRaw, 'standardChromosomes.txt'),
        'aaProps': os.path.join(dirDataAux, 'aaProps.tsv'),
        
        # results
        'seqPlots': os.path.join(dirResults_seqPlots, '{}.png'),
        'seqPlots_error': os.path.join(dirResults_seqPlots, 'errors.txt'),
        'mlPlots': os.path.join(dirResults_mlPlots, '{}.png'),
        'pca': os.path.join(dirResults_mlPlots, 'pca{}.png'),
        
        # IDPpi data
        'IDPpi_datasets': os.path.join(dirDataAux, '41598_2018_28815_MOESM2_ESM.xlsx'),
        'IDPpi_set': os.path.join(dirDataML, '{}_set{:02d}.tsv' + ext), # fill in with .format('train'/'test', #)
        'IDPpi_features': os.path.join(dirDataML, '{}_features{:02d}.tsv' + ext), # fill in with .format('train'/'test', #)
        'IDPpi_labels': os.path.join(dirDataML, '{}_labels{:02d}.tsv' + ext), # fill in with .format('train'/'test', #),
        'IDPpi_dataAll': os.path.join(dirDataML, 'dataAll.tsv' + ext),
        'IDPpi_dataAll_trim': os.path.join(dirDataML, 'dataAll_trim.tsv' + ext),
        'IDPpi_dataAll_aug': os.path.join(dirDataML, 'dataAll_aug.tsv' + ext),
        'IDPpi_labelsAll': os.path.join(dirDataML, 'labelsAll.tsv' + ext),
        'IDPpi_labelsAll_trim': os.path.join(dirDataML, 'labelsAll_trim.tsv' + ext),
        'IDPpi_labelsAll_aug': os.path.join(dirDataML, 'labelsAll_aug.tsv' + ext),
        'IDPpi_featuresAll_noTC_trim': os.path.join(dirDataML, 'featuresAll_noTC_trim.tsv' + ext),
        'IDPpi_featuresAll_noTC_trim_unc': os.path.join(dirDataML, 'featuresAll_noTC_trim.tsv'),
        'IDPpi_featuresPCA': os.path.join(dirDataML, 'featuresPCA.tsv' + ext),
        'IDPpi_featuresPCA_unc': os.path.join(dirDataML, 'featuresPCA.tsv'),
        
        'ml_train_features': os.path.join(dirDataML, 'train_features.tsv' + ext),
        'ml_train_labels': os.path.join(dirDataML, 'train_labels.tsv' + ext),
        'ml_val_features': os.path.join(dirDataML, 'val_features.tsv' + ext),
        'ml_val_labels': os.path.join(dirDataML, 'val_labels.tsv' + ext),
        'ml_test_features': os.path.join(dirDataML, 'test_features.tsv' + ext),
        'ml_test_labels': os.path.join(dirDataML, 'test_labels.tsv' + ext),
        
        'featuresDB_full_0': os.path.join(dirDataAux, 'featuresDB_full_0.RDS'),
        'featuresDB_full_50': os.path.join(dirDataAux, 'featuresDB_full_50.RDS'),
        'featuresDB_full_100': os.path.join(dirDataAux, 'featuresDB_full_100.RDS'),
        'featuresDB_small': os.path.join(dirDataAux, 'featuresDB_small.RDS')
    }
    
    if makedirs:
        for key, path in paths.items():
            if key.startswith('dir'):
                os.makedirs(path, exist_ok=True)
    
    return paths