#!/usr/bin/env python
# coding: utf-8

color_20 = ['#4428cb', '#ea5721', '#549b7c', '#337ab7', '#e1bab4', '#4b6846', '#bb61c3', '#214e75', '#a8fffd', '#420420', '#efe8ef', '#177abe', '#093646', '#fe22f6', '#43442d', '#d4da3f', '#c88dd1', '#0c3b72']
hatches = [ "..." , "---", "///" , "|" , "-" , "+" , "x", "o", "O", ".", "*", "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*"]
vega_20_scanpy = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
    '#9467bd', '#8c564b', '#e377c2',  # '#7f7f7f' removed grey
    '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896',
    '#c5b0d5', '#c49c94', '#f7b6d2',  # '#c7c7c7' removed grey
    '#dbdb8d', '#9edae5',
    '#ad494a', '#8c6d31']  # manual additions

#################################################################################################
def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:.1f}%\n({:.0f})".format(pct, absolute)

def lv_stat_barplot(data, splitby='key_labels', colorby='louvain_groups', hatch=True, ylab="% in each group", 
                    xlab="Clusters", figsize=(10, 3), startangle = 90, ncols=3, stacked=True, plotMethod='Barplot', 
                    orientation  = 'vertical', fontsize=11, bbox_to_anchor=(1.2, 0.34), color=None, save=None):
    """    User defined barplot function
    
    Parameters
    ----------
    adata
        Annotated data matrix
    splitby
        The key of the observations for spliting the bars. Default: key_labels
    colorby    
        The key of the obvervations for coloring the bars. Default: louvain_groups
    hatch
        True/False. whether to use hatch. Default: True
    ylab
        Label for y axis. Default: % in each group
    xlab
        Label for x axis. Default: Clusters
    figsize
        Figure of the figure. Default: (10,3)
    startangle
        Startangle for pie chart. Default: 90
    ncols
        Number of plots per row, onlly for pie chart. Default: 3
    stacked
        True/False. Whether to stack the bars. Default: True
    plotMethod
        Type of plot, options: Barplot, Stackedarea, Panel_Column_Chart and PieChart. Default: Barplot
    orientation
        Orientation of barplot. Default: vertical
    fontsize
        fontsize of labels in pie chart. Default: 11
    color
        User defined color for barplot. Default: None
    save
        File name for saving the plot
    
    Returns
    -------
    **User defined type of plot, either barplot or piechart**
    
    Examples
    --------
    >>> lv_stat_barplot(adata, xlab="Louvain Clusters", splitby="Louvain_Clusters", colorby="Donor", color='white', save='xxx.pdf')
    
    """
    
    if splitby in data.obs_keys():
        g1D = data.obs[splitby]
    elif splitby in data.var_names:
        g1D = data[:, splitby].X
        g1D = ['Positive' if x > 0 else 'Negative' for x in g1D]
        g1D = pd.Series(g1D, dtype = 'category')
        g1D = g1D.cat.reorder_categories(['Negative', 'Positive'])
        g1D = g1D.values
    else:
        raise ValueError('"' + splitby + '" is invalid!'
                                 + ' specify valid sample annotation, one of '
                                 + str(data.obs_keys()) + ' or a gene name '
                                 + str(data.var_names))
    
    if colorby in data.obs_keys():
        g2D = data.obs[colorby]
    elif colorby in data.var_names:
        g2D = data[:, colorby].X
        g2D = ['Positive' if x > 0 else 'Negative' for x in g2D]
        g2D = pd.Series(g2D, dtype = 'category')
        g2D = g2D.cat.reorder_categories(['Negative', 'Positive'])
        g2D = g2D.values
    else:
        raise ValueError('"' + colorby + '" is invalid!'
                                 + ' specify valid sample annotation, one of '
                                 + str(data.obs_keys()) + ' or a gene name '
                                 + str(data.var_names))
    
    df = pd.crosstab(g1D, g2D)
    df_new = df.div(df.sum(axis=1),axis=0)
    
    #print(df_new)
    if plotMethod=='Barplot':
        fig, ax = plt.subplots(figsize= figsize)
        if color is None:
            color = [[vega_20_scanpy[x]]*df_new.shape[0] for x in range(df_new.shape[1])]
            color = np.stack(color)
            #print(color)
        if orientation == 'horizontal':
            df_new = -df_new
            df_new.plot.barh(stacked=stacked, color=color, edgecolor="black", ax=ax)
            plt.xticks(np.arange(0, 101, 20)/-100)
            ax.set_ylabel(xlab)
            ax.set_xlabel(ylab)
            ax.set_xticklabels(np.arange(100, -1, -20),rotation=0)
            ax.grid()
        else:
            df_new = -df_new
            df_new.plot.bar(stacked=stacked, color=color, edgecolor="black", ax=ax)
            plt.yticks(np.arange(0, 101, 20)/-100, np.arange(100, -1, -20))
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.set_xticklabels(df_new.index,rotation=0)
            if len(data.obs[splitby].cat.categories) >= 5:
                plt.xticks(rotation=90)
            ax.grid()
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            plt.xticks(range(df_new.shape[0]),df.sum(axis=1),rotation=90)
            ax2.grid(False)
        if hatch is True:
            hatch1 = [[hatches[x]]*df_new.shape[0] for x in range(df_new.shape[1])]
            hatch1 = np.hstack(hatch1)
            for i, thisbar in enumerate(ax.patches):
                thisbar.set_hatch(hatch1[i])
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    elif plotMethod=='Stackedarea':
        hatch1 = hatches[0:df_new.shape[1]]
        if color is None:
            color = color_20[0:df_new.shape[1]]
        ax = df_new.plot.area(color=color)
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.grid()
    elif plotMethod=='Panel_Column_Chart':
        sns.set_style("whitegrid")
        ax = df_new.plot.bar(subplots=True, sharey=True, 
              figsize=(6, 5), legend=False,  
              grid=False, edgecolor='none', 
              fontsize=12)    
        plt.text(-1, 2, ylab, fontsize=12, rotation=0) #df_new.shape[1]/2
        sns.despine(left=True)
        for ax1 in ax:  # set the names beside the axes
            #ax1.lines[0].set_visible(False)  # remove ugly dashed line
            ax1.set_title('')
            sername = ax1.get_legend_handles_labels()[1][0]
            ax1.text(7, 0.5, sername, fontsize=12)
        #plt.suptitle("My panel chart", fontsize=18)
    elif plotMethod=='PieChart':
        import math
        nrow = math.ceil(len(df.index)/ncols)
        fig, axes = plt.subplots(nrow,ncols, figsize=figsize)
        for i in range(len(df.index)):
            if nrow==1:
                ax = axes[i % ncols]
            else:
                ax = axes[i // ncols, i % ncols]
            patches, texts, autotexts = ax.pie(df.iloc[i,:] ,startangle = startangle, counterclock=False, colors=color, autopct = lambda pct: func(pct, df.iloc[i,:]))  
            #[ _.set_fontsize(3) for _ in texts]
            ax.set_title(df.index[i], fontsize=fontsize) #,  loc='left'
            plt.setp(autotexts, size=8, weight="bold")
        plt.figlegend(patches, df.columns, bbox_to_anchor=bbox_to_anchor,  loc='right', fontsize=fontsize)
        fig.subplots_adjust(top=0.8,right=0.8)     
    
    plt.tight_layout()
    if save is not None:
        plt.savefig('./figures/Barplot_'+save, bbox_inches='tight')
    return df

def ReadinFiles_PBMC(h5dir, fn, Donor):
    adata = sc.read(h5dir+'/'+fn+'_withvec_EC10k/outs/filtered_gene_bc_matrices/hg19_MP22668_IFNB1_5LTR/matrix.mtx', cache=True).T
    adata.var = pd.read_csv(h5dir+'/'+fn+'_withvec_EC10k/outs/filtered_gene_bc_matrices/hg19_MP22668_IFNB1_5LTR/genes.tsv', sep='\t', header=None, names=['Ensembl Name', 'Gene Name'])
    adata.obs = pd.read_csv(h5dir+'/'+fn+'_withvec_EC10k/outs/filtered_gene_bc_matrices/hg19_MP22668_IFNB1_5LTR/barcodes.tsv', sep='\t', header=None, names = ['Barcode'])
    adata.var_names = adata.var['Gene Name']
    adata.obs_names = adata.obs['Barcode']
    #display(adata.var.head(3))
    #display(adata.obs.head(3))
    adata.obs['Vector'] = adata[:,'IFNB1-5LTR'].X.toarray()
    adata = adata[:, ~(adata.var_names == 'IFNB1-5LTR')].copy()
    adata.var_names_make_unique()
    adata.obs['Donor'] = Donor
    return adata

import louvain
def get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""
    import igraph as ig

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warning(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g

def _choose_graph(adata, obsp, neighbors_key):
    """Choose connectivities from neighbbors or another obsp column"""
    if obsp is not None and neighbors_key is not None:
        raise ValueError(
            'You can\'t specify both obsp, neighbors_key. ' 'Please select only one.'
        )

    if obsp is not None:
        return adata.obsp[obsp]
    else:
        neighbors = NeighborsView(adata, neighbors_key)
        if 'connectivities' not in neighbors:
            raise ValueError(
                'You need to run `pp.neighbors` first '
                'to compute a neighborhood graph.'
            )
        return neighbors['connectivities']
    
# source from online
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
np.random.seed(2016)
def repel_labels(ax, x, y, labels, c='r', k=0.01):
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi)

    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
    scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.items():
        pos[key] = (val*scale) + shift

    for label, data_str in G.edges():
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data',
                    xytext=pos[label], textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                                    shrinkA=0, shrinkB=0,
                                    connectionstyle="arc3", 
                                    color=c), fontsize=12)
    # expand limits
#     all_pos = np.vstack(pos.values())
#     x_span, y_span = np.ptp(all_pos, axis=0)
#     mins = np.min(all_pos-x_span*0.15, 0)
#     maxs = np.max(all_pos+y_span*0.15, 0)
#     ax.set_xlim([mins[0], maxs[0]])
#     ax.set_ylim([mins[1], maxs[1]])

def ReadinFiles_Tcells(h5dir, fn, Donor, Condition):
    adata = sc.read(h5dir+'/'+fn+'_withvec_EC10k/outs/filtered_gene_bc_matrices/hg19_MP22668_IFNB1_5LTR/matrix.mtx', cache=True).T
    adata.var = pd.read_csv(h5dir+'/'+fn+'_withvec_EC10k/outs/filtered_gene_bc_matrices/hg19_MP22668_IFNB1_5LTR/genes.tsv', sep='\t', header=None, names=['Ensembl Name', 'Gene Name'])
    adata.obs = pd.read_csv(h5dir+'/'+fn+'_withvec_EC10k/outs/filtered_gene_bc_matrices/hg19_MP22668_IFNB1_5LTR/barcodes.tsv', sep='\t', header=None, names = ['Barcode'])
    adata.var_names = adata.var['Gene Name']
    adata.obs_names = adata.obs['Barcode']
    #display(adata.var.head(3))
    #display(adata.obs.head(3))
    adata.obs['Vector'] = adata[:,'IFNB1-5LTR'].X.toarray()
    adata = adata[:, ~(adata.var_names == 'IFNB1-5LTR')].copy()
    adata.var_names_make_unique()
    adata.obs['Donor'] = Donor
    adata.obs['Condition'] = Condition
    return adata

def autolabel(rects, err, values):
    """
    Attach a text label above each bar displaying its height
    """
    seg = err[0].get_segments()
    for idx in range(len(rects)):
        height = seg[idx][1][1]
        ax.text(rects[idx].get_x() + rects[idx].get_width()/2., 1.02*height,
                '%d' % int(values[idx]),
                ha='center', va='bottom', fontsize=8, weight='bold')

