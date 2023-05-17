import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from ..tools import create_deg_df, get_nrows_maxcols, check_list, nth_repl_all, collect_deg_data, extract_groups, check_raw
from ..readwrite import save_and_show_figure
import seaborn as sns
import networkx as nx
from ..readwrite import save_and_show_figure
from ..utils import rename_net, check_if_adjustText
import textwrap
from scipy.stats import pearsonr
from matplotlib import ticker, colors

# @testfunc
def volcano_plot(adata, keys, groups=None, foldchanges_label='logfoldchanges', pvals_label='pvals_adj', gene_names='names', 
                fc_threshold=1, pval_threshold=0.05, show_thresholds=False, colorbar_label=None,
                annotate=False, correct_annotation=True, annotate_top_n=20, annotation_fontsize=8, max_cols=4,
                #marker_colors=None, 
                marker_size=6, ax_label_size=14, title_size=16,
                threshold_color='grey', up_color='green', down_color='red',
                title=None, xlim=None, ylim=None, figsize=(7,6),
                return_only=False, show=True,
                savepath=None, save_only=False, dpi_save=300):

    '''
    Plotting volcano_plot for a DEG analysis.
    Data is expected in `adata.uns[key]` as output of `sc.tl.rank_genes_groups`.

    Using `return_only=True` it is possible to return dataframes for each key/group instead of a plot.
    '''

    if groups is not None:
        groups = [groups] if isinstance(groups, str) else list(groups)

    
    # if marker_colors is None:
    #     marker_colors = ['black', 'green', 'red']
    # else:
    #     marker_colors = [marker_colors] if isinstance(marker_colors, str) else list(marker_colors)

    groups, keys, key_data_dict, key_updown_dict, params_dict = collect_deg_data(adata, keys=keys, groups=groups, foldchanges_label=foldchanges_label,
                                    pvals_label=pvals_label, gene_names=gene_names, 
                                    fc_threshold=fc_threshold, pval_threshold=pval_threshold, 
                                    added_cats=None)

    if return_only:
        return {k: pd.concat(key_updown_dict[k]) for k in key_updown_dict}
        #return key_updown_dict

    #else:
    multikeys = True if len(keys) > 1 else False
    multigroups = True if len(groups) > 1 else False

    # decide if 2d axes or 1d axes
    axis_2d = multikeys and multigroups

    # check if xlim or ylim were set beforehand
    xlim_was_set = False
    ylim_was_set = False
    if xlim is not None:
        xlim_was_set = True
    if ylim is not None:
        ylim_was_set = True

    if axis_2d:
        n_rows = len(keys)
        n_cols = len(groups)
        n_plots = n_rows * n_cols
    elif multigroups:
        # extract plotting configuration
        n_plots, n_rows, n_cols = get_nrows_maxcols(groups, max_cols)
    elif multikeys:
        n_plots, n_rows, n_cols = get_nrows_maxcols(keys, max_cols)
    else:
        n_plots = n_rows = n_cols = 1
        # print("Both groups and keys empty.")
        # return
        
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(figsize[0]*n_cols, figsize[1]*n_rows))

    if n_plots > 1:
        axs = axs.ravel()
    else:
        axs = [axs]

    i = 0
    for key in keys:
        data_dict = key_data_dict[key]
        updown_dict = key_updown_dict[key]
        reference = params_dict[key]['reference']

        groups = data_dict.keys()

        for group in groups:
            data = data_dict[group]
            updown = updown_dict[group]

            # extract parameters
            logpvals = data['-logpvals'].copy()
            logfold = data[foldchanges_label].copy()

            if xlim is None:
                xlim = (logfold.min(), logfold.max())
            if ylim is None:
                ylim = (logpvals.min(), logpvals.max()*1.1)

            up_exists = False
            down_exists = False
            if 'up' in updown.index:
                up = updown.loc['up']
                up_exists = True
            if 'down' in updown.index:
                down = updown.loc['down']
                down_exists = True

            #if len(marker_colors) == 3:
            # plot all genes
            axs[i].scatter(x=logfold, y=logpvals, color='k', s=marker_size, edgecolor='k')
            if up_exists:
                axs[i].scatter(x=up[foldchanges_label], y=-np.log10(up[pvals_label]), color=up_color, s=marker_size, edgecolor=up_color)
            if down_exists:
                axs[i].scatter(x=down[foldchanges_label], y=-np.log10(down[pvals_label]), color=down_color, s=marker_size, edgecolor=down_color)

            plot_colorbar = False

            
            # elif len(marker_colors) == 1:
            #     assert marker_colors[0] in data.columns, 'marker_colors {} not in dataframe'.format(marker_colors[0])

            #     colors = data[marker_colors[0]].values
            #     colors += np.nanmin(colors[np.nonzero(colors)])
            #     colors = -np.log10(colors)
            #     s = axs[i].scatter(x=logfold, y=logpvals, c=colors, s=6)
                
            #     plot_colorbar = True
            
            # else:
            #     print("Unknown length for marker_colors [{}]".format(len(marker_colors)))

            # set plot parameters
            axs[i].set_xlabel('log2(foldchanges)', fontsize=ax_label_size)
            axs[i].set_ylabel('log10(adjusted p-value)', fontsize=ax_label_size)

            if title is None:
                axs[i].set_title("{} vs. {}".format(group, reference), fontsize=title_size)
            else:
                axs[i].set_title(title, fontsize=title_size)

            if xlim is not None:
                axs[i].set_xlim(xlim[0], xlim[1])
            if ylim is not None:
                axs[i].set_ylim(ylim[0], ylim[1])

            if show_thresholds:
                axs[i].axhline(-np.log10(pval_threshold), c=threshold_color, linestyle='dashed')
                axs[i].axvline(-fc_threshold, c=threshold_color, linestyle='dashed')
                axs[i].axvline(fc_threshold, c=threshold_color, linestyle='dashed')
            
            if annotate:
                # plot annotation for upregulated genes
                #for elem in [up, down]:
                texts = []
                x_up = list(up[foldchanges_label][:annotate_top_n])
                x_down = list(down[foldchanges_label][:annotate_top_n])
                y_up = list(-np.log10(up[pvals_label][:annotate_top_n]))
                y_down = list(-np.log10(down[pvals_label][:annotate_top_n]))

                if up_exists and down_exists:
                    x_values = x_up + x_down
                    y_values = y_up + y_down
                elif up_exists:
                    x_values = x_up
                    y_values = y_up
                elif down_exists:
                    x_values = x_down
                    y_values = y_down
                else:
                    print("Neither up or down-regulated genes in group {}".format(group))
                #x_values = list(up[foldchanges_label][:annotate_top_n]) + list(down[foldchanges_label][:annotate_top_n])
                #y_values = list(-np.log10(up[pvals_label][:annotate_top_n])) + list(-np.log10(down[pvals_label][:annotate_top_n]))
                names = list(up[gene_names][:annotate_top_n]) + list(down[gene_names][:annotate_top_n])

                for x, y, g in zip(x_values, y_values, names):
                    if (x > xlim[0] and x < xlim[1]) & (y > ylim[0] and y < ylim[1]):
                        texts.append(axs[i].text(x, y, g, size=annotation_fontsize))
            
                if correct_annotation:
                    from adjustText import adjust_text
                    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=axs[i])

            # counter
            i+=1

            # check if ylim has to be reset
            if not xlim_was_set:
                xlim = None
            if not ylim_was_set:
                ylim = None
    
    # check if there are empty plots remaining
    while i < n_rows * n_cols:
        # remove empty plots
        axs[i].set_axis_off()
        i += 1

    if plot_colorbar:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.83, 0.2, 0.01, 0.5])
        clb = plt.colorbar(s, cax=cbar_ax, orientation='vertical')
        if colorbar_label is not None:
            clb.set_label(colorbar_label)

    save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, dpi_save=dpi_save)

    # if savepath is not None:
    #     print("Saving figure to file " + savepath)
    #     plt.savefig(savepath, dpi=dpi_save, bbox_inches='tight')
    #     print("Saved.")
    # if save_only:
    #     plt.close(fig)
    # elif show:
    #     return plt.show()
    # else:
    #     if isinstance(axs, list):
    #         axs = axs[0]
    #     return fig, axs

def correlate_significance_two_groups(adata, deg_key, spatialde_key, groupby, 
deg_splits='all', threshold = 0.05, increment = 10e-16, offset = 1.02, vmin = -3, vmax = 3, 
linewidth=0.5, colormap = 'PRGn', edgecolor = 'k', adjust=False, savepath=None, save_only=False, show=True, dpi_save=300):

    if adjust:
        from adjustText import adjust_text

    # extract the groups that were used in the DEG analysis
    splitby = adata.uns[deg_key]['params']['groupby']

    # check whether to use all possible groups or only a subset of it
    if deg_splits == 'all':
        deg_splits = list(adata.obs[splitby].unique())
    else:
        deg_splits = [elem for elem in deg_splits if elem in adata.obs[splitby].unique()]
    if len(deg_splits) > 2:
        print("Too many groups ({} instead of 2)".format(len(deg_splits)))
        return
    elif len(deg_splits) < 2:
        print("Found only {} of the give groups in `adata.obs[{}]`.".format(str(len(deg_splits)), splitby))

    deg_df = {deg_split: create_deg_df(adata, keys=deg_key, groups=deg_split) for deg_split in deg_splits}

    exp_conditions = {elem[0]: elem[1] for elem in list(set([tuple((well, g)) for well, g in zip(adata.obs[groupby], adata.obs[splitby])]))}
    # create results dataframe from spatialde results
    results_df = {}
    for well in adata.uns[spatialde_key].keys():
        results_df[well] = adata.uns[spatialde_key][well]["results"]
        results_df[well][groupby] = well
        results_df[well][splitby] = exp_conditions[well]

    results_df = pd.concat(results_df, ignore_index=True)
    results_df = results_df[[splitby, groupby, 'g', 'l', 'qval']]

    # find patterned genes
    patterned_genes = results_df.query('qval < 0.05')['g'].unique()

    # filter for patterned genes
    mask = [gene in patterned_genes for gene in results_df['g']]
    results_df = results_df[mask]

    # group the results dataframe
    grouped_df = results_df.groupby([splitby, "g"]).mean()
    grouped_df["std"] = results_df.groupby([splitby, "g"]).std()["qval"].values
    grouped_df["count"] = results_df.groupby([splitby, "g"]).count()["qval"].values
    grouped_df = grouped_df.dropna()

    # filter for patterned genes that are in group1
    mask = [gene in patterned_genes for gene in deg_df[deg_splits[0]].names]
    deg_group1 = deg_df[deg_splits[0]][mask]
    deg_group1 = deg_group1.sort_values('names', ascending=True)

    # sort grouped df by index
    grouped_df = grouped_df.sort_index(ascending=True)

    # make sure that the genes are in both groups. Filter unique genes out.
    in_both = grouped_df.xs(deg_splits[0]).index[grouped_df.xs(deg_splits[0]).index.get_level_values('g').isin(grouped_df.xs(deg_splits[1]).index.get_level_values('g'))]
    grouped_df = grouped_df.loc[(slice(None), in_both), :]

    deg_group1 = deg_group1[deg_group1.names.isin(in_both)]

    #return deg_group1

    # extract q-values for both groups
    qval_group1 = grouped_df.loc[(deg_splits[0])].qval.values
    qval_group2 = grouped_df.loc[(deg_splits[1])].qval.values

    # extract gene names
    gene_names = grouped_df.loc[(deg_splits[0])].index.values

    # extract logfoldchange
    logfold_group1 = deg_group1['logfoldchanges']

    # plotting
    fig, ax = plt.subplots(2,2, figsize=(10,10))
    ax = ax.ravel()

    ## Plot merged plots
    # plot non-significant points
    ns_mask = (qval_group1>=threshold) & (qval_group2>=threshold)
    xns, yns = -np.log10(qval_group1[ns_mask]+increment), -np.log10(qval_group2[ns_mask]+increment)
    s = ax[2].scatter(x=xns, y=yns, c=logfold_group1[ns_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)
    #clb = plt.colorbar(s, ax=ax[3])
    #clb.set_label('log10(Fold change) {} vs. {}'.format(groups[0], groups[1]))

    # plot all-significant points
    s_mask = (qval_group1<threshold) & (qval_group2<threshold)
    xs, ys = -np.log10(qval_group1[s_mask]+increment), -np.log10(qval_group2[s_mask]+increment)
    ax[2].scatter(x=xs, y=ys, c=logfold_group1[s_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)

    # plot group1-significant points
    y_mask = (qval_group1<threshold) & (qval_group2>=threshold)
    xy, yy = -np.log10(qval_group1[y_mask]+increment), -np.log10(qval_group2[y_mask]+increment)
    ax[2].scatter(x=xy, y=yy, c=logfold_group1[y_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)

    # plot group2-significant points
    o_mask = (qval_group1>=threshold) & (qval_group2<threshold)
    xo, yo = -np.log10(qval_group1[o_mask]+increment), -np.log10(qval_group2[o_mask]+increment)
    ax[2].scatter(x=xo, y=yo, c=logfold_group1[o_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)

    ## Plot detail plots
    # plot all-significant points as detail
    s_textmask = s_mask & (np.absolute(logfold_group1)>1) # mask for text annotation
    xs, ys = -np.log10(qval_group1[s_mask]+increment), -np.log10(qval_group2[s_mask]+increment)
    ax[1].scatter(x=xs, y=ys, c=logfold_group1[s_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)

    # plot annotation for all-significant subset
    xs_text, ys_text = -np.log10(qval_group1[s_textmask]+increment), -np.log10(qval_group2[s_textmask]+increment)

    texts_all = []
    for x, y, g in zip(xs_text, ys_text, gene_names[s_textmask]):
        texts_all.append(ax[1].text(x*offset, y*offset, g, size=8))

    if adjust:
        adjust_text(texts_all, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=ax[1])

    # plot group1-significant points
    y_mask = (qval_group1<threshold) & (qval_group2>=threshold)
    xy, yy = -np.log10(qval_group1[y_mask]+increment), -np.log10(qval_group2[y_mask]+increment)
    ax[3].scatter(x=xy, y=yy, c=logfold_group1[y_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)

    # plot annotation for group1 subset
    texts_group1 = []
    for x, y, g in zip(xy, yy, gene_names[y_mask]):
        texts_group1.append(ax[3].text(x*offset, y*offset, g, size=8))

    if adjust:
        adjust_text(texts_group1, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=ax[3])

    # plot group2-significant points
    o_mask = (qval_group1>=threshold) & (qval_group2<threshold)
    xo, yo = -np.log10(qval_group1[o_mask]+increment), -np.log10(qval_group2[o_mask]+increment)
    ax[0].scatter(x=xo, y=yo, c=logfold_group1[o_mask], vmin=vmin, vmax=vmax, edgecolor='k', linewidth=linewidth, cmap=colormap)

    # plot annotation for group2 subset
    texts_group2 = []
    for x, y, g in zip(xo, yo, gene_names[o_mask]):
        texts_group2.append(ax[0].text(x*offset, y*offset, g, size=8))
    
    if adjust:    
        adjust_text(texts_group2, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=ax[0])


    for i in range(len(ax)):
        ax[i].axhline(-np.log10(threshold), ls='dashed', color='k')
        ax[i].axvline(-np.log10(threshold), ls='dashed', color='k')

    for i in range(len(ax)):
        ax[i].set_xlabel("-log10(qval) {}".format(deg_splits[0]))
        ax[i].set_ylabel("-log10(qval) {}".format(deg_splits[1]))
        
    fig.subplots_adjust(bottom=0.25)
    cbar_ax = fig.add_axes([0.1, 0.15, 0.8, 0.02])
    clb = plt.colorbar(s, cax=cbar_ax, orientation='horizontal')
    clb.set_label('log2(Fold change) {} vs. {}'.format(deg_splits[0], deg_splits[1]))

    save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, dpi_save=dpi_save, tight=False)

def corr_heatmap(adata, keys, groupby, groups=None,
    max_cols=4, use_raw=True,
    clustermap=False, cluster_names=None,
    savepath=None, save_only=False, dpi_save=300):

    if groups is None:
        groups = list(adata.obs[groupby].unique())
    else:
        groups = [groups] if isinstance(groups, str) else list(groups)

    n_plots, n_rows, max_cols = get_nrows_maxcols(groups, max_cols)

    ## Plotting
    if not clustermap:
        fig, axs = plt.subplots(n_rows, max_cols, figsize=(8*max_cols, 7*n_rows))

        if n_plots > 1:
            axs = axs.ravel()
        else:
            axs = [axs]

    for i, group in enumerate(groups):
        # get data for one well
        subset = extract_groups(adata, groupby=groupby, groups=group)
        # check if plotting raw data
        X, var, var_names = check_raw(subset, use_raw=use_raw)

        var_ids = [var_names.get_loc(key) for key in keys]
        # subset = subset[:, keys].copy()

        # extract expression
        # expr = pd.DataFrame(subset.X, columns=subset.var_names)
        expr = pd.DataFrame(X[:, var_ids], columns=var_names[var_ids])

        # calculate correlation matrix
        excorr = expr.corr()

        if clustermap:
            assert len(groups) == 1

            #cluster_names = pd.Series(cluster_names)
            lut = dict(zip(cluster_names.unique(), "rbg"))

            colors = cluster_names.map(lut)
            sns.clustermap(excorr, cmap='RdYlGn', vmin=-1, vmax=1, 
            #row_cluster=False, col_cluster=False, 
            row_colors=colors, col_colors=colors)
        else:
            sns.set(font_scale=1)
            sns.heatmap(excorr, cmap='RdYlGn', vmin=-1, vmax=1, ax=axs[i])

            axs[i].set_title('Pearson correlation\n{}: {}'.format(groupby, group))

    plt.tight_layout()
    save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save)

def violinplot(data, x, y, ax, hue=None,
            order=None, ylabel="# of {} per spot", jitter=0.4, hue_order=None, scale='area',
            show=True, **kwargs):
    # plotting
    sns.violinplot(x=x, y=y, data=data, ax=ax, hue=hue, 
        inner=None, dodge=False, order=order, hue_order=hue_order, scale=scale, **kwargs)
    sns.stripplot(x=x, y=y, data=data, ax=ax, hue=None, size=1, color='black', dodge=True, jitter=jitter, order=order)
    
    # settings
    ax.set_ylim(0, None)
    ax.tick_params(axis='x', which='major', labelsize=16, rotation=45)
    ax.tick_params(axis='y', which='major', labelsize=16, rotation=0)
    ax.set_ylabel(ylabel.format(y), fontsize=16)
    
    if hue is not None:
        # set legend
        legend = ax.legend(title="Method",
                        #loc="upper left",
                        #bbox_to_anchor=(1, 1.03),
                        fontsize=16)
        plt.setp(legend.get_title(),fontsize=16)
    
    if show:
        plt.tight_layout()
        return plt.show()

def pearson_facet_plot(adata, groupby, 
                       savepath=None, save_only=False, dpi_save=300):
    
    df = adata.to_df()
    df = pd.concat([adata.obs, df], axis=1)
    long_df = df.melt(id_vars=adata.obs.columns)

    indices = long_df[groupby].unique()
    n = len(indices)

    fig, axs = plt.subplots(n, n, figsize=(8*n, 6*n))

    for r, rid in enumerate(indices):
        for c, cid in enumerate(indices):
            x = long_df.query('{} == "{}"'.format(groupby, cid))['value'].values
            y = long_df.query('{} == "{}"'.format(groupby, rid))['value'].values

            axs[r, c].scatter(x, y)
            axs[r, c].set_title("r = {}".format(round(pearsonr(x, y)[0], 2)), fontsize=20)
            if r == n-1:
                axs[r, c].set_xlabel(cid, fontsize=20)

            if c == 0:
                axs[r, c].set_ylabel(rid, fontsize=20)
                
    save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, dpi_save=dpi_save)
    
    
def barplot_custom(data, x, y, 
                   hue=None, 
                   order=["Periportal", "Midzone", "Pericentral"], 
                   hue_order=["young", "old"],
                   palette=["white", "grey"],
                   max_cols=4,
                   figsize=(6,4),
                   stat_annot=True,
                   test='t-test_welch',
                   comparisons_correction=None,
                   verbose=1,
                   savepath=None, fig=None, save_only=False, show=True, dpi_save=300
                   ):
    
    y = [y] if isinstance(y, str) else list(y)
    n_plots, n_rows, max_cols = get_nrows_maxcols(keys=y, max_cols=max_cols)
    
    fig, axs = plt.subplots(n_rows, max_cols, figsize=(figsize[0]*max_cols, figsize[1]*n_rows))
    
    if n_plots > 1:
        axs = axs.ravel()
    else:
        axs = [axs]
            
    for i, yy in enumerate(y):
        sns.barplot(data=data, x=x, y=yy, 
                    hue=hue, 
                    order=order,
                    hue_order=hue_order,
                    dodge=True,
                    linewidth=1, 
                    edgecolor='k',
                    palette=palette,
                    ci="sd", errcolor='k', errwidth=1.5, capsize=0.2,
                    ax=axs[i]
                )

        sns.stripplot(data=data, x=x, y=yy, 
                      hue=hue,
                    order=order,
                    hue_order=hue_order,
                    color="k",
                    dodge=True, 
                    #legend=False,
                    ax=axs[i]
                )
        
        # plot only the barplot legends
        handles, labels = axs[i].get_legend_handles_labels()
        l = axs[i].legend(handles[2:4], 
                          labels[2:4], 
                          #bbox_to_anchor=(1.05, 1), 
                          #loc=2, 
                          #borderaxespad=0.
                          )
        
        axs[i].set_title(yy)
        axs[i].set_ylabel("Expression")
        
        if stat_annot:
            #from statannot import add_stat_annotation
            from statannotations.Annotator import Annotator
            # add statistical asterisks
            box_pairs=[
            (("Pericentral", "old"), ("Pericentral", "young")),
            (("Midzone", "old"), ("Midzone", "young")),
            (("Periportal", "old"), ("Periportal", "young"))
            ]
            
            annot = Annotator(axs[i], pairs=box_pairs, data=data, 
                              x=x, y=yy, 
                              order=order,
                              hue=hue, hue_order=hue_order
                              )
            annot.configure(test=test, comparisons_correction=comparisons_correction, verbose=verbose)
            annot.apply_test()
            annot.annotate()
            
            # add_stat_annotation(axs[i], data=data, x=x, y=yy, 
            #                     order=order,
            #                     hue=hue, hue_order=hue_order,
            #                     box_pairs=box_pairs,
            #                     test='t-test_ind', comparisons_correction=None,
            #                     loc='inside', verbose=verbose)
            
    save_and_show_figure(savepath=savepath, fig=fig, save_only=save_only, show=show, dpi_save=dpi_save)
    plt.show()
    

def filter_limits(df, sign_limit=None, lFCs_limit=None):

    # Define limits if not defined
    if sign_limit is None:
        sign_limit = np.inf
    if lFCs_limit is None:
        lFCs_limit = np.inf

    # Filter by absolute value limits
    msk_sign = df['pvals'] < np.abs(sign_limit)
    msk_lFCs = np.abs(df['logFCs']) < np.abs(lFCs_limit)
    df = df.loc[msk_sign & msk_lFCs]

    return df
    
def plot_volcano(logFCs, pvals, contrast, 
                 size=None, 
                 name=None, net=None, top=5, source='source', target='target',
                 weight='weight',
                 size_scale_factor=10,
                 sign_thr=0.05, lFCs_thr=0.5, 
                 sign_limit=None, lFCs_limit=None, 
                 return_data=False,
                 figsize=(7, 5),
                 dpi=100, ax=None, return_fig=False, 
                 save=None, save_only=False, show=True, dpi_save=300
                 ):
    """
    Function adapted from https://github.com/saezlab/decoupler-py/blob/main/decoupler/plotting.py
    
    Plot logFC and p-values. If name and net are provided, it does the same for the targets of a selected source.
    Parameters
    ----------
    logFCs : DataFrame
        Data-frame of logFCs (contrasts x features).
    pvals : DataFrame
        Data-frame of p-values (contrasts x features).
    contrast : str
        Name of the contrast (row) to plot.
    name : str, None
        Name of the source to plot. If None, plot classic volcano (without subsetting targets).
    net : DataFrame, None
        Network dataframe. If None, plot classic volcano (without subsetting targets).
    top : int
        Number of top differentially expressed features.
    source : str
        Column name in net with source nodes.
    target : str
        Column name in net with target nodes.
    weight : str, None
        Column name in net with weights. If none, set to None.
    sign_thr : float
        Significance threshold for p-values.
    lFCs_thr : float
        Significance threshold for logFCs.
    sign_limit : float
        Limit of p-values to plot in -log10.
    lFCs_limit : float
        Limit of logFCs to plot in absolute value.
    figsize : tuple
        Figure size.
    dpi : int
        DPI resolution of figure.
    ax : Axes, None
        A matplotlib axes object. If None returns new figure.
    return_fig : bool
        Whether to return a Figure object or not.
    save : str, None
        Path to where to save the plot. Infer the filetype if ending on {`.pdf`, `.png`, `.svg`}.
    Returns
    -------
    fig : Figure, None
        If return_fig, returns Figure object.
    df  : DataFrame, None
        If return_data, returns DataFrame object.
    """
    
    # Load plotting packages
    at = check_if_adjustText()

    # Match dfs
    index = logFCs.index.intersection(pvals.index)
    columns = logFCs.columns.intersection(pvals.columns)
    logFCs = logFCs.loc[index, columns]
    pvals = pvals.loc[index, columns]

    # Transform sign_thr
    logsign_thr = -np.log10(sign_thr)

    # Plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    # Check for net
    if net is not None:
        # Rename nets
        net = rename_net(net, source=source, target=target, weight=weight)

        # Get max and if + and -
        max_n = np.std(np.abs(net['weight']), ddof=1)*2
        has_neg = np.any(net['weight'] < 0)

        # Filter by shared targets
        if name is None:
            raise ValueError('If net is given, name cannot be None.')
        df = net[net['source'] == name].set_index('target')
        df['logFCs'] = logFCs.loc[[contrast]].T
        df['pvals'] = -np.log10(pvals.loc[[contrast]].T)
        df = df[~np.any(pd.isnull(df), axis=1)]
        df = filter_limits(df, sign_limit=sign_limit, lFCs_limit=lFCs_limit)

        if has_neg:
            vmin = -max_n
        else:
            vmin = 0
        df.plot.scatter(x='logFCs', y='pvals', c='weight', cmap='coolwarm',
                        vmin=vmin, vmax=max_n, sharex=False, ax=ax)
        ax.set_title('{0} | {1}'.format(contrast, name))
    else:
        df = logFCs.loc[[contrast]].T.rename({contrast: 'logFCs'}, axis=1)
        df['pvals'] = pvals.loc[[contrast]].T
        df['logpvals'] = -np.log10(pvals.loc[[contrast]].T)
        
        if size is not None:
            df['means'] = size.loc[[contrast]].T
            df['means_scaled'] = df['means'] * size_scale_factor
            
        df = df[~np.any(pd.isnull(df), axis=1)]
        df = filter_limits(df, sign_limit=sign_limit, lFCs_limit=lFCs_limit)
        df['weight'] = 'gray'
        df.loc[(df['logFCs'] >= lFCs_thr) & (df['logpvals'] >= logsign_thr), 'weight'] = '#D62728'
        df.loc[(df['logFCs'] <= -lFCs_thr) & (df['logpvals'] >= logsign_thr), 'weight'] = '#1F77B4'
        
        if size is None:
            df.plot.scatter(x='logFCs', y='logpvals', c='weight', sharex=False, ax=ax)
        else:
            df.plot.scatter(x='logFCs', y='logpvals', c='weight', s='means_scaled', 
                            edgecolor='k',
                            linewidth=0.5,
                            sharex=False, ax=ax)
        
        ax.set_title('{0}'.format(contrast))

    # Draw sign lines
    ax.axhline(y=logsign_thr, linestyle='--', color="black")
    ax.axvline(x=lFCs_thr, linestyle='--', color="black")
    ax.axvline(x=-lFCs_thr, linestyle='--', color="black")

    # Plot top sign features
    signs = df[(np.abs(df['logFCs']) >= lFCs_thr) & (df['logpvals'] >= logsign_thr)].sort_values('logpvals', ascending=False)
    signs = signs.iloc[:top]
    
    # Add labels
    ax.set_ylabel('-log10(pvals)')
    texts = []
    for x, y, s in zip(signs['logFCs'], signs['logpvals'], signs.index):
        texts.append(ax.text(x, y, s))
    if len(texts) > 0:
        at.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'), ax=ax)
                    
    if size is not None:
        # produce a legend with a cross-section of sizes from the scatter
        pc = ax.collections[0] # get PathCollections
        handles, labels = pc.legend_elements(prop="sizes",
                                            alpha=1,
                                            func=lambda x: x / size_scale_factor,
                                            markerfacecolor='gray',
                                            markeredgecolor='k'
                                            )
        legend = ax.legend(handles, labels, 
                            #loc="upper right", 
                            title="Mean expression"
                            )
    
    if save is not None:
        if ax is None:
            save_and_show_figure(savepath=save, fig=fig, save_only=save_only, show=show, dpi_save=dpi_save)
        else:
            raise ValueError("ax is not None, cannot save figure.")

    if return_fig:
        if ax is None:
            return fig
    elif return_data:
        return df

def custom_dotplot(c, s, genes,
                   size_scale_factor = 100,
                   size_thresh = -np.log(0.05),
                   size_legend_title = "-log10(pval)",
                   order = None,
                   figsize=(8, 6),
                   xmargin=0.1, ymargin=0.2,
                   xticks_rotation=90,
                   savepath=None,
                   save_only=False,
                   show=True, 
                   dpi_save=300,
                   **kwargs
                  ):
    
    # select genes
    c = c.loc[:, genes]
    
    # reshape data
    lfcmelt = c.melt(var_name="gene", value_name="color", ignore_index=False)
    lfcmelt = lfcmelt.reset_index(names="zone")
    lfcmelt = lfcmelt.set_index(["zone", "gene"])

    pvalmelt = s.melt(var_name="gene", value_name="size", ignore_index=False)
    pvalmelt = pvalmelt.reset_index(names="zone")
    pvalmelt = pvalmelt.set_index(["zone", "gene"])
    
    # combine pval and lfc data
    lfcmelt["size"] = pvalmelt
    lfcmelt = lfcmelt.reset_index()
    
    # order zones
    if order is not None:
        lfcmelt = lfcmelt.set_index("zone").loc[order].reset_index()
    
    # select genes
    #subset = lfcmelt.loc[lfcmelt.gene.isin(genes), :].copy()
    
    lfcmelt["size_scaled"]= lfcmelt["size"] * size_scale_factor
    
    # calculate significance score
    #subset["-log10(pval)"] = -np.log10(subset["pval"])*size_scale_factor
        
    # plotting
    fig, ax = plt.subplots(1,1, figsize=figsize)
    lfcmelt.plot.scatter(x="gene", y="zone", c="color", s="size_scaled", ax=ax, 
                        norm=colors.CenteredNorm(), 
                        **kwargs)
    
    # identify significant measurements and draw circle around them
    sig = lfcmelt.copy()
    sig.loc[sig["size"] < size_thresh, "size_scaled"] = np.nan
    sig.plot.scatter(x="gene", y="zone", c="None", s="size_scaled", 
                        alpha=1, edgecolor='k', linewidth=1,
                        ax=ax
                        )
    
    # produce a legend with a cross-section of sizes from the scatter
    pc = ax.collections[0] # get PathCollections    
    handles, labels = pc.legend_elements(prop="sizes", 
                                         num=ticker.FixedLocator([-np.log10(0.1),
                                                                  -np.log10(0.05), 
                                                                  -np.log10(0.01)]),
                                         alpha=1,
                                         func=lambda x: x / size_scale_factor,
                                         markerfacecolor='gray',
                                         markeredgecolor='k'
                                        )
    legend = ax.legend(handles, labels, 
                       bbox_to_anchor=(1.3, 1),
                        loc="upper left", 
                        title=size_legend_title
                        )
    
    plt.xticks(rotation=xticks_rotation)
    plt.margins(x=xmargin, y=ymargin)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
