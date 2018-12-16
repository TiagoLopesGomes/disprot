import matplotlib.pyplot as plt
import numpy as np
import pdb

# CONSTANTS
# - aaOrder suggested by Adrian Sanborn in email on 12/1/2018

aaOrder = ['P', 'G', 'A', 'E', 'D', 'K', 'R', 'H', 'Q', 'N', 'S', 'T', 'V', 'L', 'I', 'M', 'C', 'Y', 'F', 'W']
aaCharge = {
    'K': 1, 'R': 1, # basic positive
    'D': -1, 'E': -1, # acidic negative
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0, 'P': 0, 'V': 0, 'W': 0, # neutral nonpolar
    'H': 0, 'N': 0, 'Q': 0, 'S': 0, 'T': 0, 'Y': 0 # neutral polar
}
aaChargeColors = {
    'K': '#ffcccc', 'R': '#ffcccc', # basic positive: light red
    'E': '#cceeff', 'D': '#cceeff' # acidic negative: light blue
}

def plotSeq(seq, charOrder=aaOrder, ax=None, highlight=aaChargeColors, **kwargs):
    '''
    2-D plot of sequence.
    - Vertical axis: unique characters. Numerical range 0:len(charOrder)
    - Horizontal axis: position. Numerical range 1:len(seq)
    
    Args
    - seq: str
        Sequence of characters to visualize
    - charOrder: list or dict (str -> int). default=aaOrder
        list: the first character corresponds to lowest position on the vertical axis
        dict: specifies vertical position of each character
        None: uses alphabetical order of unique characters in `seq`
    - ax: matplotlib.axes.Axes. default=None
        Plotting axes.
    - highlight: dict (str -> matplotlib colors)
        Row highlight colors
    - **kwargs
        Arguments to pass to ax.scatter() (if ax is given) or matplotlib.pyplot.scatter()
        (if ax is None)
    
    Returns: matplotlib.collections.PathCollection
      Returned by ax.scatter() / matplotlib.pyplot.scatter()
    '''
    
    uChars = set(seq) # unique characters in sequence
    seqLen = len(seq)
    
    # ----------------------------
    # Validate / process arguments
    if charOrder is None:
        print("No character order given. Assuming alphabetical order.")
        charOrder = sorted(uChars)
    
    if any([char not in charOrder for char in uChars]):
        raise ValueError('One or more characters in `seq` is not in `charOrder`.')
    
    if type(charOrder) is not dict:
        charOrderMap = {char: index for index, char in enumerate(charOrder)}
    
    if ax is None:
        ax = plt.axes()
    # -----------------------------
    
    ## Generate points to plot
    y = [charOrderMap[char] for char in seq]
    x = list(range(1,seqLen+1))

    for char, color in highlight.items():
        if char in charOrder:
            ax.barh(y=charOrderMap[char], width=seqLen, height=1, color=color, alpha=0.5)
    paths = ax.scatter(x, y, s=200, c='black', marker='|', **kwargs)
    ax.set_ylabel('Character')
    ax.set_yticks(list(range(0,len(charOrder))))
    ax.set_yticklabels(charOrder)
    
    return paths

def avgCharge(seq, window=2, charCharge=aaCharge):
    '''
    Moving average of sequence charge.
    
    Args
    - seq: str
        Sequence of characters to analyze
    - window: int. default=2
        Moving average width parameter. Must be >=0.
    - charCharge: dict (str -> int). default=aaCharge
        Charge of each character
    
    Returns: numpy.ndarray of shape (len(seq),)
      Moving average of charges centered at each inner character. Towards the ends of the sequence,
      the averages may not be centered.
    '''
    
    # Validate arguments
    if (type(window) is not int or window < 0):
        raise ValueError('`window` must be a nonnegative integer.')
    if any([char not in charCharge for char in set(seq)]):
        raise ValueError('One or more characters in `seq` is not in `charCharge`.')
    
    seqLen = len(seq)
    avg = np.zeros(seqLen)
    for i in range(min(seqLen, window)):
        avg[i] = np.mean([charCharge[char] for char in seq[:(i+window+1)]])
    
    for i in range(min(seqLen, window), seqLen):
        avg[i] = np.mean([charCharge[char] for char in seq[(i-window):(i+window+1)]])
    
    return avg

def plotCharge(seq, window=2, charCharge=aaCharge, ax=None, **kwargs):
    '''
    Charge per amino acid analysis.
    - Vertical axis: charge Numerical range -1:1
    - Horizontal axis: position. Numerical range 1:len(seq)
    
    Args
    - seq: str
        Sequence of characters to visualize
    - window: int. default=2
        Moving average width parameter. Must be >=0.
    - charCharge: dict (str -> int). default=aaCharge
        Charge of each character
    - ax: matplotlib.axes.Axes. default=None
        Plotting axes.
    - **kwargs
        Arguments to pass to ax.bar() (if ax is given) or matplotlib.pyplot.bar()
        (if ax is None)
    
    Returns: matplotlib.collections.PathCollection
      Returned by ax.scatter() / matplotlib.pyplot.scatter()
    '''
    
    if ax is None:
        ax = plt.axes()

    x = list(range(1,len(seq)+1))
    y = avgCharge(seq, window, charCharge)
    
    color_map = {-1: 'blue', 0: 'black', 1: 'red'}
    colors = [color_map[sign] for sign in np.sign(y).astype(int)]
    ax.set_ylim(-1,1)
    ax.set_ylabel('Charge')
    return ax.bar(x, y, color=colors, edgecolor=colors)

def plotDisorder(seq, disorder, name='', ax=None, **kwargs):
    '''
    Args
    - seq: str
        Sequence of characters to visualize
    - disorder: list of 2-tuple
        List of (start, end) regions of disorder, 1-indexed
    - name: str. default=''
        Name of sequence (e.g., gene name).
    - ax: matplotlib.axes.Axes. default=None
        Plotting axes.
    - **kwargs
        Arguments to pass to ax.barh() (if ax is given) or matplotlib.pyplot.barh()
        (if ax is None)
    '''
    
    if ax is None:
        ax = plt.axes()
    ax.set_yticklabels('')
    ax.set_yticks([])
    ax.set_ylabel(name, rotation='horizontal', horizontalalignment='right')
    
    numDRs = len(disorder)
    if numDRs == 0:
        # no disordered regions
        return ax.barh(y=0, width=0)
    
    # Validate arguments
    if any([start > end for start, end in disorder]):
        raise ValueError('Disorder regions must be specified with a start position <= end position')
    if max([end for start, end in disorder]) > len(seq):
        raise ValueError('Disorder region specified past sequence length')
    
    starts, ends = zip(*disorder)
    starts = np.array(starts)
    ends = np.array(ends)
    
    return ax.barh(y=0, width=ends-starts, left=starts, color='red')

def plotSeqAndProps(seq, disorder, window=2, name='', scaleFigure=True, **kwargs):
    '''
    2-D sequence plot, disorder horizontal bar plot, and charge plot.
    
    Args
    - seq: str
        Sequence of characters to visualize
    - disorder: list of 2-tuple
        List of (start, end) regions of disorder, 1-indexed
    - window: int. default=2
        Moving average width parameter. Must be >=0.
    - name: str. default=''
        Name of sequence (e.g., gene name).
    - scaleFigure: bool. default=True
        Scale width of plot to sequence length.
    - **kwargs
        Arguments passed to plotting functions
    
    Returns: (matplotlib.figure.Figure, numpy.ndarray)
      matplotlib figure and axes for the plots. The numpy.ndarray has shape (3,) and dtype 'O' (object)
      containing matplotlib.axes.Axes objects.
    
    Based on Figure 4A from Boija, A. et al. Transcription Factors Activate Genes
    through the Phase-Separation Capacity of Their Activation Domains. Cell 175, 1–14 (2018).
    '''
    
    # Figure size constraints
    # - matplotlib limits printed figure size to 65536 x 65536 (2^16) dots
    # - assuming figure will be saved with dpi=300, the limit in inches is 65536/300 = 218.45 inches
    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True,
                            gridspec_kw={'height_ratios':[6.25,0.5,3.25]},
                            figsize=(min(len(seq)/50, 218), 8) if scaleFigure else (8,8))
    plotSeq(seq, ax=axs[0], **kwargs)
    plotDisorder(seq, disorder, name='Disorder', ax=axs[1], **kwargs)
    plotCharge(seq, ax=axs[2], **kwargs)
    axs[2].set_xlim(0, len(seq))
    axs[2].set_xlabel('Position')
    axs[0].set_title(name)
    fig.tight_layout()
    return (fig, axs)

def plotSeqAndPropsFromID(id, idType, idMap, disorderDB, proteome=None, name_col=None):
    '''
    Wrapper for plotSeqAndProps() given a protein ID.
    
    Args
    - id: str
        Protein identifier (e.g., Ensembl protein ID, UniProtKB accession ID, D2P2 ID, or gene symbol)
    - idType: str
        Type of protein identifer: 'ensembl_protein_id', 'uniprot_id', 'd2p2_id', or `name`
    - idMap: pandas.DataFrame
        1:1 map between different ID types.
    - disorder: pandas.DataFrame
        Disorder database. Must contain 'd2p2_id', 'start', and 'end'.
        Disordered regions should not overlap.
    - proteome: pandas.DataFrame. default=None
        Proteome (e.g., from Ensembl or UniProt). Must contain `idType` and 'seq'.
        If None, `idMap` is assumed to be the proteome as well.
    - name_col: str. default=None
        Column in `idMap` to use as the name of the sequence in the plot.
        If None, `id` is passed as the name.
    
    Returns: see plotSeqAndProps()
    '''
    
    if proteome is None:
        proteome = idMap
    
    seq = proteome.loc[proteome[idType] == id, 'seq'].values.item()
    d2p2_id = idMap.loc[idMap[idType] == id, 'd2p2_id'].values.item()
    disorder = disorderDB.loc[disorderDB['d2p2_id'] == d2p2_id, ['start', 'end']] \
                         .sort_values('start').to_records(index=False).tolist()
    if name_col is None:
        name = id
    else:
        name = idMap.loc[idMap[idType] == id, name_col].values.item()
    return plotSeqAndProps(seq, disorder, name=name)