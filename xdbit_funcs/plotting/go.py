import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import textwrap
import networkx as nx
from ..tools import get_nrows_maxcols, check_list, nth_repl_all
from ..readwrite import save_and_show_figure

def go_plot(enrichment, 
            style="dot",
            color_key=None, size_key=None, groups=None,
    fig=None, axis=None,
    max_to_plot=25, max_cols=4, cmap='viridis', cmin=None, cmax=None,
    ax_label_size=16, markersize=240, 
    figsize=(8,6),
    remove_yticklabels=False,
    xtick_label_size=16, ytick_label_size=16,
    clb_label_size=16, clb_tick_label_size=16, clb_pos=None, clb_norm=False,
    title_size=16, max_line_length=25, custom_headers=None,
    max_name_length=60,
    value_to_plot='Enrichment score', 
    colorcol=None,
    sortby=None, sort=True, additional_sortby=None,
    name_key='name', libraries=None, ascending=False,
    savepath=None, save_only=False, show=True, dpi_save=300
    ):

    # get groups from index of enrichment dataframe
    extracted_groups = enrichment.index.unique(level=0).tolist()

    if groups is not None:
        groups = [groups] if isinstance(groups, str) else list(groups)
        
        groups = check_list(groups, extracted_groups)

        if len(groups) == 0:
            return
    else:
        groups = extracted_groups

    # check possible custom headers
    if custom_headers is not None:
        custom_headers = [custom_headers] if isinstance(custom_headers, str) else list(custom_headers)

    # check for libraries
    if libraries is None:
        libraries = list(enrichment['source'].unique())
    else:
        # convert libraries to list if they are string
        libraries = [libraries] if isinstance(libraries, str) else list(libraries)

        # check if all libraries are in dataframe
        notin = [elem for elem in libraries if elem not in enrichment['source'].unique()]
        assert len(notin) == 0, "Following libraries could not be found in the `source` column: {}".format(notin)

        # filter for libraries
        enrichment = enrichment[enrichment['source'].isin(libraries)].copy()

    # sort dataframe
    if sortby is None:
        sortby = value_to_plot

    if sort:
        enrichment.sort_values(sortby, ascending=ascending, inplace=True)

    # shorten to max length if necessary
    if max_to_plot is not None:
        enrichment = enrichment.groupby(level=0).head(max_to_plot)
        
    if additional_sortby is not None:
        enrichment.sort_values(additional_sortby, inplace=True)
        #enrichment[name_key] = ["{} ({})".format(a,b) for a,b in zip(enrichment[name_key], enrichment[additional_sortby])]

    # Prepare names for plotting
    # Shorten name and add GO term name if too short
    if max_name_length is not None:
        enrichment[name_key] = ["{}...({})".format(n[:max_name_length], go) if len(n)>max_name_length else n for go, n in zip(enrichment['native'], enrichment[name_key])]
    # Insert line breaks if sentence too long
    #enrichment[name_key] = [nth_repl_all(elem, ' ', '\n', 1) for elem in enrichment[name_key]]

    # get min and max for the colorbar
    if color_key is not None:
        if clb_norm:
            cmax = enrichment[color_key].max() if cmax is None else cmax
            cmin = enrichment[color_key].min() if cmin is None else cmin
        else:
            cmax = None
            cmin = None

    # Plotting
    if axis is None:
        n_plots, n_rows, n_cols = get_nrows_maxcols(groups, max_cols=max_cols)
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(figsize[0]*n_cols, figsize[1]*n_rows))
        
    else:
        assert len(groups) == 1, "Enrichment dataframe contains more than one group in first index level 0."
        axs = axis
        n_plots = 1
        show = False

    axs = axs.ravel() if n_plots > 1 else [axs]

    for i, group in enumerate(groups):
        if group in enrichment.index.unique(level=0):
            # select data
            df = enrichment.xs(group).copy()

            # if max_to_plot is not None and len(df) > max_to_plot:
            #         df = df[:max_to_plot]
            if color_key is not None:
                color = df[color_key]
            else:
                color = 'k'

            if size_key is not None:
                markersize=df[size_key]*5
            
            if max_line_length is not None:
                # introduce line breaks if the names are too long
                df[name_key] = [textwrap.fill(elem, width=max_line_length, break_long_words=True) for elem in df[name_key]]

            # plotting
            if style == "dot":
                s = axs[i].scatter(df[value_to_plot], df[name_key], 
                    c=color, cmap=cmap, 
                    s=markersize, 
                    edgecolors='k')
            elif style == "dot_fixed":
                s = axs[i].scatter(
                    x = [0] * len(df[name_key]), 
                    y = df[name_key], 
                    s=df[value_to_plot] * 50,
                    c=color, 
                    cmap=cmap, 
                    edgecolors='k')
            elif style == "bar":
                ys = df[name_key]
                ys_pos = np.arange(len(ys))
                s = axs[i].barh(
                    ys_pos,
                    width=df[value_to_plot], 
                    height=0.8,
                    color='k', 
                #cmap=cmap, 
                #s=markersize, 
                #edgecolors='k'
                )
                axs[i].set_yticks(ys_pos, labels=ys)
            else:
                raise ValueError('Invalid `style` parameter ("{}"). Possible options'.format(style))

            if custom_headers is None:
                axs[i].set_title("{}\n{}".format(group, libraries), fontsize=title_size)
            else:
                axs[i].set_title(custom_headers[i], fontsize=title_size)
            
            axs[i].invert_yaxis()
            axs[i].set_xlabel(value_to_plot, fontsize=ax_label_size)
            axs[i].tick_params(axis='x', which='major', labelsize=xtick_label_size)
            axs[i].tick_params(axis='y', which='major', labelsize=ytick_label_size)
            
            if style == "dot":
                axs[i].grid(axis='y')
            
            if style == "dot_fixed":
                axs[i].spines['top'].set_visible(False)
                axs[i].spines['right'].set_visible(False)
                axs[i].spines['bottom'].set_visible(False)
                axs[i].spines['left'].set_visible(False)
                axs[i].get_xaxis().set_ticks([])
                axs[i].set_xlabel("")

            if remove_yticklabels:
                axs[i].set_yticklabels([])

            # cbar_ax = fig.add_axes([0.95, 0.12, 0.02, 0.75])
            # clb = plt.colorbar(s, cax=cbar_ax, orientation='vertical')
            # clb.set_label('Fraction of patterned genes', loc='center')

            if style in ["dot", "dot_fixed"]:
                if color_key is not None:
                    if clb_pos is None:
                        clb = fig.colorbar(s, ax=axs[i])
                        clb.set_label(color_key, loc='center', fontsize=clb_label_size)
                        clb.ax.tick_params(labelsize=clb_tick_label_size)
                        if clb_norm:
                            clb.mappable.set_clim(cmin, cmax)
                    else:
                        if i == clb_pos:
                            clb = fig.colorbar(s, ax=axs[clb_pos])
                            clb.set_label(color_key, loc='center', fontsize=clb_label_size)
                            clb.ax.tick_params(labelsize=clb_tick_label_size)
                            if clb_norm:
                                clb.mappable.set_clim(cmin, cmax)
                                


            if size_key is not None:
                kw = dict(prop="sizes", num=5,
                    #color=s.cmap(0.7), 
                    #fmt="$ {x:.2f}",
                    #func=lambda s: np.sqrt(s/.3)/3
                    #func=lambda s: np.sqrt(s)
                    )
                size_legend = axs[i].legend(*s.legend_elements(**kw, alpha=0.6), 
                                    #markerscale=0.5,
                                    loc="lower right", title=size_key)
                
            if colorcol is not None:
                color_dict = {a: b for a, b in zip(df["name"], df[colorcol])}
                for xtick in axs[i].get_yticklabels():
                    xtick.set_color(color_dict[xtick.get_text()])
                # for xtick, color in zip(axs[i].get_xticklabels(), enrichment[colorcol]):
                #     xtick.set_color(color)


        else:
            #print('No significant results for selected libraries {} in group {}'.format(libraries, group))
            axs[i].set_title("{}\n{}".format(group, libraries), fontsize=12, fontweight='bold')
            axs[i].text(0.5,0.5, 'No significant results for selected\nlibraries {}\nin group {}'.format(libraries, group), 
                        va='center', ha='center', fontsize=20)
            axs[i].set_axis_off()

    # check if there are empty plots remaining
    if n_plots > 1:
        while i < n_rows * n_cols - 1:
            i += 1
            # remove empty plots
            axs[i].set_axis_off()

        fig.tight_layout()

    fig.tight_layout()
    if show:
        save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save, fig=fig)
    else:
        return fig, axs


def go_dotplot(enrichment, color_key=None, size_key=None, groups=None,
    fig=None, axis=None,
    max_to_plot=25, max_cols=4, cmap='viridis', cmin=None, cmax=None,
    ax_label_size=16, markersize=240, 
    figsize=(8,6),
    remove_yticklabels=False,
    xtick_label_size=16, ytick_label_size=16,
    clb_label_size=16, clb_tick_label_size=16, clb_pos=None, clb_norm=False,
    title_size=16, max_line_length=25, custom_headers=None,
    max_name_length=60,
    value_to_plot='Enrichment score', 
    colorcol=None,
    sortby=None, sort=True,
    name_key='name', libraries=None, ascending=False,
    savepath=None, save_only=False, show=True, dpi_save=300
    ):

    # get groups from index of enrichment dataframe
    extracted_groups = enrichment.index.unique(level=0).tolist()

    if groups is not None:
        groups = [groups] if isinstance(groups, str) else list(groups)
        
        groups = check_list(groups, extracted_groups)

        if len(groups) == 0:
            return
    else:
        groups = extracted_groups

    # check possible custom headers
    if custom_headers is not None:
        custom_headers = [custom_headers] if isinstance(custom_headers, str) else list(custom_headers)

    # check for libraries
    if libraries is None:
        libraries = list(enrichment['source'].unique())
    else:
        # convert libraries to list if they are string
        libraries = [libraries] if isinstance(libraries, str) else list(libraries)

        # check if all libraries are in dataframe
        notin = [elem for elem in libraries if elem not in enrichment['source'].unique()]
        assert len(notin) == 0, "Following libraries could not be found in the `source` column: {}".format(notin)

        # filter for libraries
        enrichment = enrichment[enrichment['source'].isin(libraries)].copy()

    # sort dataframe
    if sortby is None:
        sortby = value_to_plot

    if sort:
        enrichment.sort_values(sortby, ascending=ascending, inplace=True)

    # shorten to max length if necessary
    if max_to_plot is not None:
        enrichment = enrichment.groupby(level=0).head(max_to_plot)

    # Prepare names for plotting
    # Shorten name and add GO term name if too short
    if max_name_length is not None:
        enrichment[name_key] = ["{}...({})".format(n[:max_name_length], go) if len(n)>max_name_length else n for go, n in zip(enrichment['native'], enrichment[name_key])]
    # Insert line breaks if sentence too long
    #enrichment[name_key] = [nth_repl_all(elem, ' ', '\n', 1) for elem in enrichment[name_key]]

    # get min and max for the colorbar
    if color_key is not None:
        if clb_norm:
            cmax = enrichment[color_key].max() if cmax is None else cmax
            cmin = enrichment[color_key].min() if cmin is None else cmin
        else:
            cmax = None
            cmin = None

    # Plotting
    if axis is None:
        n_plots, n_rows, n_cols = get_nrows_maxcols(groups, max_cols=max_cols)
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(figsize[0]*n_cols, figsize[1]*n_rows))
        
    else:
        assert len(groups) == 1, "Enrichment dataframe contains more than one group in first index level 0."
        axs = axis
        n_plots = 1
        show = False

    axs = axs.ravel() if n_plots > 1 else [axs]

    for i, group in enumerate(groups):
        if group in enrichment.index.unique(level=0):
            # select data
            df = enrichment.xs(group).copy()

            # if max_to_plot is not None and len(df) > max_to_plot:
            #         df = df[:max_to_plot]
            if color_key is not None:
                color = df[color_key]
            else:
                color = 'k'

            if size_key is not None:
                markersize=df[size_key]*5
            
            if max_line_length is not None:
                # introduce line breaks if the names are too long
                df[name_key] = [textwrap.fill(elem, width=max_line_length, break_long_words=True) for elem in df[name_key]]

            # plotting
            s = axs[i].scatter(df[value_to_plot], df[name_key], 
                c=color, cmap=cmap, 
                s=markersize, 
                edgecolors='k')
            
            # s = axs[i].bar(df[value_to_plot], df[name_key], 
            #     color=color, 
            #     #cmap=cmap, 
            #     #s=markersize, 
            #     #edgecolors='k'
            #     )
            
            

            if custom_headers is None:
                axs[i].set_title("{}\n{}".format(group, libraries), fontsize=title_size)
            else:
                axs[i].set_title(custom_headers[i], fontsize=title_size)
            
            axs[i].invert_yaxis()
            axs[i].set_xlabel(value_to_plot, fontsize=ax_label_size)
            axs[i].tick_params(axis='x', which='major', labelsize=xtick_label_size)
            axs[i].tick_params(axis='y', which='major', labelsize=ytick_label_size)
            axs[i].grid(axis='y')

            if remove_yticklabels:
                axs[i].set_yticklabels([])

            # cbar_ax = fig.add_axes([0.95, 0.12, 0.02, 0.75])
            # clb = plt.colorbar(s, cax=cbar_ax, orientation='vertical')
            # clb.set_label('Fraction of patterned genes', loc='center')

            if color_key is not None:
                if clb_pos is None:
                    clb = fig.colorbar(s, ax=axs[i])
                    clb.set_label(color_key, loc='center', fontsize=clb_label_size)
                    clb.ax.tick_params(labelsize=clb_tick_label_size)
                    if clb_norm:
                        clb.mappable.set_clim(cmin, cmax)
                else:
                    if i == clb_pos:
                        clb = fig.colorbar(s, ax=axs[clb_pos])
                        clb.set_label(color_key, loc='center', fontsize=clb_label_size)
                        clb.ax.tick_params(labelsize=clb_tick_label_size)
                        if clb_norm:
                            clb.mappable.set_clim(cmin, cmax)

            if size_key is not None:
                kw = dict(prop="sizes", num=5,
                    #color=s.cmap(0.7), 
                    #fmt="$ {x:.2f}",
                    #func=lambda s: np.sqrt(s/.3)/3
                    #func=lambda s: np.sqrt(s)
                    )
                size_legend = axs[i].legend(*s.legend_elements(**kw, alpha=0.6), 
                                    #markerscale=0.5,
                                    loc="lower right", title=size_key)
                
            if colorcol is not None:
                color_dict = {a: b for a, b in zip(df["name"], df["color"])}
                for xtick in axs[i].get_yticklabels():
                    xtick.set_color(color_dict[xtick.get_text()])
                # for xtick, color in zip(axs[i].get_xticklabels(), enrichment[colorcol]):
                #     xtick.set_color(color)


        else:
            #print('No significant results for selected libraries {} in group {}'.format(libraries, group))
            axs[i].set_title("{}\n{}".format(group, libraries), fontsize=12, fontweight='bold')
            axs[i].text(0.5,0.5, 'No significant results for selected\nlibraries {}\nin group {}'.format(libraries, group), 
                        va='center', ha='center', fontsize=20)
            axs[i].set_axis_off()

    # check if there are empty plots remaining
    if n_plots > 1:
        while i < n_rows * n_cols - 1:
            i += 1
            # remove empty plots
            axs[i].set_axis_off()

        fig.tight_layout()

    fig.tight_layout()
    if show:
        save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save, fig=fig)
    else:
        return fig, axs

def go_barplot(enrichment, max_to_plot=25, max_cols=4, groups=None, name_dict=None, value_to_plot='Enrichment score', libraries=None,
cutoff=0.05, ascending=False, name_key='name',
savepath=None, save_only=False, show=True, dpi_save=300):

    extracted_groups = list(enrichment.index.unique(level=0))

    if groups is not None:
        groups = [groups] if isinstance(groups, str) else list(groups)
        groups = check_list(groups, extracted_groups)

        if len(groups) == 0:
            return
    else:
        groups = extracted_groups.copy()

    # check for libraries
    if libraries is None:
        libraries = list(enrichment['source'].unique())
    else:
        # convert libraries to list if they are string
        libraries = [libraries] if isinstance(libraries, str) else list(libraries)

        # filter for libraries
        enrichment = enrichment[enrichment['source'].isin(libraries)].copy()

    # sort dataframe
    enrichment.sort_values(value_to_plot, ascending=ascending, inplace=True)

    # shorten to max length if necessary
    if max_to_plot is not None:
        enrichment = enrichment.groupby(level=0).head(max_to_plot).copy()

    # Prepare names for plotting
    # Shorten name and add GO term name if too short
    enrichment.loc[:, name_key] = ["{}...({})".format(n[:60], go) if len(n)>60 else n for go, n in zip(enrichment['native'], enrichment[name_key])]
    # Insert line breaks if sentence too long
    enrichment.loc[:, name_key] = [nth_repl_all(elem, ' ', '\n', 5) for elem in enrichment[name_key]]

    n_plots, n_rows, n_cols = get_nrows_maxcols(groups, max_cols=max_cols)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(10*n_cols, 8*n_rows))
    axs = axs.ravel() if n_plots > 1 else [axs]

    for i, group in enumerate(groups):
        if group in enrichment.index.unique(level=0):
            # select data
            df = enrichment.xs(group)

            # set title
            libs = df['source'].unique()

            if name_dict is not None:
                name = name_dict[group]
                axs[i].set_title("{}\n{}".format(name, libs))
            else:
                axs[i].set_title("{}\n{}".format(group, libs))
            
            # filter by cutoff
            df = df[df[value_to_plot] > -np.log10(cutoff)]
            
            if len(df) > 0:
                df.sort_values(value_to_plot, ascending=ascending, inplace=True)

                sns.barplot(x=value_to_plot, y=name_key, data=df, ax=axs[i])
            else:
                print('No significant results for group "{}" with cutoff {}.'.format(group, cutoff))
                axs[i].set_axis_off()
        
        else:
            #print('No significant results for selected libraries {} in group {}'.format(libraries, group))
            axs[i].set_title("{}\n{}".format(group, libraries), fontsize=12, fontweight='bold')
            axs[i].text(0.5,0.5, 'No significant results for\nselected libraries {} in group {}'.format(libraries, group), 
                        va='center', ha='center', fontsize=20)
            axs[i].set_axis_off()

    while i < n_rows * n_cols - 1:
        i += 1
        # remove empty plots
        axs[i].set_axis_off()
        
    fig.tight_layout()

    save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save)




def go_network(enrichment, color_key=None, groups=None,
max_to_plot=25, max_cols=4, cmap='viridis', cmin=None, cmax=None, normalize_cmap=True, #clb_minmax=(0,1), 
edgeweight_factor=0.2, frame_factor=1.2, nodesize_factor=200, seed=0,
value_to_plot='Enrichment score', name_key='name', libraries=None, ascending=False, size_legend=False, clb_label=None,
savepath=None, save_only=False, dpi_save=300):
    
    # get groups from index of enrichment dataframe
    extracted_groups = list(enrichment.index.unique(level=0))

    if groups is not None:
        groups = [groups] if isinstance(groups, str) else list(groups)
        groups = check_list(groups, extracted_groups)

        if len(groups) == 0:
            return
    else:
        groups = extracted_groups

    # check for libraries
    if libraries is None:
        libraries = list(enrichment['source'].unique())
    else:
        # convert libraries to list if they are string
        libraries = [libraries] if isinstance(libraries, str) else list(libraries)

        # filter for libraries
        enrichment = enrichment[enrichment['source'].isin(libraries)].copy()


    # sort dataframe
    enrichment.sort_values(value_to_plot, ascending=ascending, inplace=True)

    # shorten to max length if necessary
    if max_to_plot is not None:
        enrichment = enrichment.groupby(level=0).head(max_to_plot).copy()

    # get min and max for the colorbar
    if color_key is not None:
        if normalize_cmap:
            cmax = enrichment[color_key].max() if cmax is None else cmax
            cmin = enrichment[color_key].min() if cmin is None else cmin
        else:
            cmax = None
            cmin = None
    
    # Prepare names for plotting
    # Shorten name and add GO term name if too short 
    enrichment[name_key] = ["{}...({})".format(n[:60], go) if len(n)>60 else n for go, n in zip(enrichment['native'], enrichment[name_key])]
    # Insert line breaks if sentence too long
    enrichment[name_key] = [nth_repl_all(elem, ' ', '\n', 2) for elem in enrichment[name_key]]

    # Plotting
    n_plots, n_rows, n_cols = get_nrows_maxcols(groups, max_cols=max_cols)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9*n_cols, 6*n_rows))
    axs = axs.ravel() if n_plots > 1 else [axs]

    for i, group in enumerate(groups):
        
        if group in enrichment.index.unique(level=0):
            # select data
            df = enrichment.xs(group)

            adj_mtx = pd.DataFrame([[len(set(a) & set(b)) if ida != idb else 0 
                                    for ida, a in zip(df.index, df['intersections'])] 
                                    for idb, b in zip(df.index, df['intersections'])])

            G = nx.from_pandas_adjacency(adj_mtx)

            mapping = {node: name for node, name in zip(list(G.nodes()), df['name'])}
            G = nx.relabel_nodes(G, mapping)

            pos = nx.spring_layout(G, k=4, seed=seed)  # positions for all nodes
            #pos = nx.circular_layout(G)  # positions for all nodes
            #pos = nx.random_layout(G)  # positions for all nodes
            #pos = nx.nx_agraph.graphviz_layout(G)

            edges = G.edges()
            weights = [G[u][v]['weight'] for u,v in edges]

            if color_key is not None:
                colors = df[color_key]
            else:
                colors = 'lightcoral'

            # if clb_minmax is None:
            #     clb_minmax = (colors.min(), colors.max())

            # nodes
            nodes = nx.draw_networkx_nodes(G, pos, node_size=df[value_to_plot] * nodesize_factor, ax=axs[i], node_color=colors, cmap=cmap)

            # edges
            nx.draw_networkx_edges(G, pos, edgelist=edges, width=np.array(weights) * edgeweight_factor, edge_color='lightgreen', ax=axs[i])

            # labels
            nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif", font_weight='bold', ax=axs[i])

            axs[i].set_title("{}\n{}".format(group, libraries), fontsize=12, fontweight='bold')
            axs[i].set_xlim([frame_factor*x for x in axs[i].get_xlim()])
            axs[i].set_ylim([frame_factor*y for y in axs[i].get_ylim()])
            axs[i].axis("off")

            if color_key is not None:
                # set colorbar
                clb = fig.colorbar(nodes, ax=[axs[i]], orientation='vertical', fraction=0.05, 
                #anchor=(-0.5, 0.5),
                aspect=40
                )

                if clb_label is None:
                    clb.set_label(color_key, loc='center')
                else:
                    clb.set_label(clb_label, loc='center')

                clb.mappable.set_clim(cmin, cmax)

            if size_legend:
                # good website: https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_with_legend.html
                # size legend
                x = [pos[n][0] for n in pos.keys()]
                y = [pos[n][1] for n in pos.keys()]

                s = axs[i].scatter(x, y, s=df[value_to_plot] * nodesize_factor, c='k')
                # produce a legend with a cross section of sizes from the scatter
                handles, labels = s.legend_elements(prop="sizes", alpha=0.6, num=3)
                legend2 = axs[i].legend(handles, labels, loc="upper left", title="-log(Adjusted p-value)")
        else:
            print('No significant results for selected libraries {} in group {}'.format(libraries, group))
            axs[i].set_title("{}\n{}".format(group, libraries), fontsize=12, fontweight='bold')
            axs[i].set_axis_off()

    # check if there are empty plots remaining
    while i < n_rows * n_cols - 1:
        i += 1
        # remove empty plots
        axs[i].set_axis_off()

    #plt.tight_layout()

    save_and_show_figure(savepath=savepath, save_only=save_only, dpi_save=dpi_save)