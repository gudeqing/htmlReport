import os
from functools import partial
from collections import OrderedDict
from glob import glob
import pandas as pd
import numpy as np
import json
import textwrap
from plotly import tools
import plotly.graph_objs as go
from plotly.offline import plot as plt
import plotly.io as pio
from jinja2 import Template
plt = partial(plt, auto_open=False)


def get_color_pool(n):
    import colorlover
    if n <= 12:
        return colorlover.scales['12']['qual']['Paired']

    from colorsys import hls_to_rgb
    color_pool = []
    for i in np.arange(60., 360., 300. / n):
        hue = i / 300.
        rand_num = np.random.random_sample()
        lightness = (50 + rand_num * 10) / 100.
        saturation = (90 + rand_num * 10) / 100.
        rgb = hls_to_rgb(hue, lightness, saturation)
        color_pool.append(tuple([int(x * 255) for x in rgb]))
    return colorlover.to_rgb(color_pool)


def get_marker_pool(n):
    maker_pool = [
        0, 'circle', 100, 'circle-open', 200, 'circle-dot', 300,
        'circle-open-dot', 1, 'square', 101, 'square-open', 201,
        'square-dot', 301, 'square-open-dot', 2, 'diamond', 102,
        'diamond-open', 202, 'diamond-dot', 302,
        'diamond-open-dot', 3, 'cross', 103, 'cross-open', 203,
        'cross-dot', 303, 'cross-open-dot', 4, 'x', 104, 'x-open',
        204, 'x-dot', 304, 'x-open-dot', 5, 'triangle-up', 105,
        'triangle-up-open', 205, 'triangle-up-dot', 305,
        'triangle-up-open-dot', 6, 'triangle-down', 106,
        'triangle-down-open', 206, 'triangle-down-dot', 306,
        'triangle-down-open-dot', 7, 'triangle-left', 107,
        'triangle-left-open', 207, 'triangle-left-dot', 307,
        'triangle-left-open-dot', 8, 'triangle-right', 108,
        'triangle-right-open', 208, 'triangle-right-dot', 308,
        'triangle-right-open-dot', 9, 'triangle-ne', 109,
        'triangle-ne-open', 209, 'triangle-ne-dot', 309,
        'triangle-ne-open-dot', 10, 'triangle-se', 110,
        'triangle-se-open', 210, 'triangle-se-dot', 310,
        'triangle-se-open-dot', 11, 'triangle-sw', 111,
        'triangle-sw-open', 211, 'triangle-sw-dot', 311,
        'triangle-sw-open-dot', 12, 'triangle-nw', 112,
        'triangle-nw-open', 212, 'triangle-nw-dot', 312,
        'triangle-nw-open-dot', 13, 'pentagon', 113,
        'pentagon-open', 213, 'pentagon-dot', 313,
        'pentagon-open-dot', 14, 'hexagon', 114, 'hexagon-open',
        214, 'hexagon-dot', 314, 'hexagon-open-dot', 15,
        'hexagon2', 115, 'hexagon2-open', 215, 'hexagon2-dot',
        315, 'hexagon2-open-dot', 16, 'octagon', 116,
        'octagon-open', 216, 'octagon-dot', 316,
        'octagon-open-dot', 17, 'star', 117, 'star-open', 217,
        'star-dot', 317, 'star-open-dot', 18, 'hexagram', 118,
        'hexagram-open', 218, 'hexagram-dot', 318,
        'hexagram-open-dot', 19, 'star-triangle-up', 119,
        'star-triangle-up-open', 219, 'star-triangle-up-dot', 319,
        'star-triangle-up-open-dot', 20, 'star-triangle-down',
        120, 'star-triangle-down-open', 220,
        'star-triangle-down-dot', 320,
        'star-triangle-down-open-dot', 21, 'star-square', 121,
        'star-square-open', 221, 'star-square-dot', 321,
        'star-square-open-dot', 22, 'star-diamond', 122,
        'star-diamond-open', 222, 'star-diamond-dot', 322,
        'star-diamond-open-dot', 23, 'diamond-tall', 123,
        'diamond-tall-open', 223, 'diamond-tall-dot', 323,
        'diamond-tall-open-dot', 24, 'diamond-wide', 124,
        'diamond-wide-open', 224, 'diamond-wide-dot', 324,
        'diamond-wide-open-dot', 25, 'hourglass', 125,
        'hourglass-open', 26, 'bowtie', 126, 'bowtie-open', 27,
        'circle-cross', 127, 'circle-cross-open', 28, 'circle-x',
        128, 'circle-x-open', 29, 'square-cross', 129,
        'square-cross-open', 30, 'square-x', 130, 'square-x-open',
        31, 'diamond-cross', 131, 'diamond-cross-open', 32,
        'diamond-x', 132, 'diamond-x-open', 33, 'cross-thin', 133,
        'cross-thin-open', 34, 'x-thin', 134, 'x-thin-open', 35,
        'asterisk', 135, 'asterisk-open', 36, 'hash', 136,
        'hash-open', 236, 'hash-dot', 336, 'hash-open-dot', 37,
        'y-up', 137, 'y-up-open', 38, 'y-down', 138,
        'y-down-open', 39, 'y-left', 139, 'y-left-open', 40,
        'y-right', 140, 'y-right-open', 41, 'line-ew', 141,
        'line-ew-open', 42, 'line-ns', 142, 'line-ns-open', 43,
        'line-ne', 143, 'line-ne-open', 44, 'line-nw', 144,
        'line-nw-open'
    ]
    return sorted([x for x in maker_pool if type(x) == int])[: n]


def draw(fig: go.Figure, prefix='', outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    for format in formats:
        out_name = os.path.join(outdir, '{}.{}'.format(prefix, format))
        if format == 'html':
            plt(fig, filename=out_name)
        else:
            if format in ['svg', 'pdf']:
                scale = 1
            pio.write_image(fig, format=format, file=out_name, height=height, width=width, scale=scale)


def gene_body_coverage(files:list, outdir=os.getcwd(), file_from='RSeQC',
                       formats=('html',), height:int=None, width:int=None, scale=3):
    layout = go.Layout(title="geneBodyCoverage")
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        fig = go.Figure(layout=go.Layout(title='{}'.format(sample)))
        if file_from == 'RSeQC':
            data = pd.read_table(each, header=None, index_col=0)
            normalize_y = data.iloc[1, :]/data.iloc[1, :].max()
        else:
            # assume data from Picard
            data = pd.read_table(each, header=None, skiprows=11)
            data = data.transpose()
            normalize_y = data.iloc[1, :]
        fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        all_fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        prefix = '{}.geneBodyCoverage'.format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    prefix = '{}.geneBodyCoverage'.format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def fragment_length(files:list, outdir=os.getcwd(), min_len=50, max_len=600,
                    formats=('html',), height:int=None, width:int=None, scale=3):
    layout = go.Layout(
        title="Fragment length distribution",
        xaxis=dict(title='Fragment length'),
        yaxis=dict(title='Probability'),
    )
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=0)
        data = data.loc[(data['frag_median'] >= min_len) & (data['frag_median'] <= max_len)]
        fig = go.Figure(layout=go.Layout(
            title='{}'.format(sample),
            xaxis=dict(title='Fragment length'),
            yaxis=dict(title='probability'),
        ))
        fig.add_histogram(x=data["frag_median"], histnorm='probability', name=sample)
        all_fig.add_histogram(x=data["frag_median"], histnorm='probability', name=sample)
        prefix = "{}.fragmentLengthDistribution".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    prefix = '{}.fragmentLengthDistribution'.format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def inner_distance(files:list, outdir, min_dist=-250, max_dist=250,
                   formats=('html',), height:int=None, width:int=None, scale=3):
    """
    抽样1000000得到的统计结果，第三列的值可能是：
    readPairOverlap
    sameTranscript=No,dist=genomic
    sameTranscript=Yes,nonExonic=Yes,dist=genomic
    sameTranscript=Yes,sameExon=No,dist=mRNA
    sameTranscript=Yes,sameExon=Yes,dist=mRNA
    """
    groups = [
        'sameTranscript=No,dist=genomic',
        'sameTranscript=Yes,nonExonic=Yes,dist=genomic',
        'readPairOverlap',
        'sameTranscript=Yes,sameExon=No,dist=mRNA',
        'sameTranscript=Yes,sameExon=Yes,dist=mRNA',
    ]
    layout = go.Layout(
        title="Inner distance distribution",
        xaxis=dict(title='Inner Distance'),
        yaxis=dict(title='Probability'),
    )
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0)
        values = [data[2][data[2] == x].shape[0] for x in groups]
        norm_values = [x/sum(values) for x in values]
        percents = dict(zip(groups, norm_values))
        data = data[(data[1] >= min_dist) & (data[1] <= max_dist)]
        fig = go.Figure(layout=go.Layout(
                title='{}'.format(sample),
                xaxis=dict(title='Inner Distance'),
                yaxis=dict(title='Frequency'),
        ))
        for g in groups:
            target = data[1][data[2] == g]
            name = "{}({:.2%})".format(g, percents[g])
            fig.add_histogram(x=target, name=name)
        all_fig.add_histogram(x=data[1][data[2] != "sameTranscript=No,dist=genomic"], histnorm='probability', name=sample)
        prefix = "{}.InnerDistanceDistribution".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    prefix = "{}.InnerDistanceDistribution".format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def read_distribution(files:list, outdir, formats=('html',), height:int=None, width:int=None, scale=3):
    all_data = list()
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, index_col=0, skiprows=4, skipfooter=1, sep='\s+', engine='python')
        data = data.loc[:, 'Tag_count']
        data.name = sample
        all_data.append(data)
        fig = go.Figure(layout=go.Layout(
            title='{}'.format(sample),
        ))
        fig.add_pie(labels=data.index, values=data)
        prefix = "{}.ReadDistribution".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    df = pd.concat(all_data, axis=1).T
    # print(df.head())
    data = [go.Bar(x=df.index, y=df[x], name=x) for x in df.columns]
    layout = go.Layout(
        title="Read distribution",
        xaxis=dict(title='Sample'),
        barmode='stack'
    )
    fig = go.Figure(data=data, layout=layout)
    prefix = "{}.ReadDistribution".format('samples')
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def read_duplication(files:list, outdir=os.getcwd(), max_dup=500,
                     formats=('html',), height:int=None, width:int=None, scale=3):
    traces = list()
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=None)
        data = data[data.iloc[:, 0] <= max_dup]
        trace = go.Scatter(x=data.iloc[:, 0], y=data.iloc[:, 1], name=sample, mode='markers')
        traces.append(trace)
        layout = go.Layout(
            title=sample,
            xaxis=dict(title='Occurrence'),
            yaxis=dict(title='UniqReadNumber', type='log'),
        )
        fig = go.Figure([trace], layout=layout)
        prefix = "{}.ReadDuplication".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)

    layout = go.Layout(
        title="Read duplication",
        xaxis=dict(title='Occurrence'),
        yaxis=dict(title='UniqReadNumber', type='log'),
    )
    fig = go.Figure(traces, layout=layout)
    prefix = "{}.ReadDuplication".format('samples')
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def exp_saturation(files:list, outdir=os.getcwd(), outlier_limit=5,
                   formats=('html',), height:int=None, width:int=None, scale=3):
    all_fig = tools.make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Q1',
            'Q2',
            'Q3',
            'Q4',
        ),
        shared_xaxes=False
    )
    color_pool = get_color_pool(len(files))
    for exp_file, sample_color in zip(files, color_pool):
        sample = os.path.basename(exp_file).split('.', 1)[0]
        data = pd.read_table(exp_file, header=0, index_col=0)
        # plot deviation
        describe = data['100'].describe()
        regions = [
            (describe['min'], describe['25%']),
            (describe['25%'], describe['50%']),
            (describe['50%'], describe['75%']),
            (describe['75%'], describe['max']),
        ]
        errors = data.apply(lambda column: column / data.loc[:, '100'], axis=0)
        errors = ((errors - 1)*100).abs()
        plot_data = list()
        for lower, upper in regions:
            tmp = errors[(data['100'] >= lower) & (data['100'] <= upper)]
            plot_data.append(tmp)

        fig = tools.make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Q1: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[0]),
                'Q2: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[1]),
                'Q3: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[2]),
                'Q4: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[3]),
            ),
            shared_xaxes=False
        )
        for ind, each in enumerate(plot_data):
            x = 1 if (ind+1) / 2 <= 1 else 2
            y = 1 if (ind+1) % 2 == 1 else 2
            for col in each.columns:
                upper_limit = each[col].describe()['50%']*outlier_limit
                y_data = each[col][each[col] <= upper_limit]
                # box = go.Box(y=y_data, name=col, showlegend=False, boxpoints=False)
                box = go.Box(y=y_data, name=col, showlegend=False)
                fig.append_trace(box, x, y)
            line = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines', name='median')
            fig.append_trace(line, x, y)
            if ind != 0:
                line2 = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines',
                                   legendgroup=sample, name=sample, line=dict(color=sample_color), showlegend=False)
            else:
                line2 = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines',
                                   legendgroup=sample, name=sample, line=dict(color=sample_color))
            all_fig.append_trace(line2, x, y)
        prefix = "{}.TPM.Saturation".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)

    # plot all sample together
    prefix = "{}.TPM.Saturation".format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def exp_pca(exp_table, row_sum_cutoff=0.1, exp_cutoff=0.1, cv_cutoff=0.01,
            explained_ratio=0.95, outdir=os.getcwd(), group_dict=None,
            formats=('html',), height:int=None, width:int=None, scale=3):
    from sklearn import decomposition, preprocessing
    data = pd.read_table(exp_table, header=0, index_col=0)
    data = data[data.sum(axis=1) >= row_sum_cutoff]
    pass_state = data.apply(lambda x: sum(x > exp_cutoff), axis=1)
    data = data[pass_state >= int(data.shape[1])/3]
    data = data[data.std(axis=1)/data.mean(axis=1) > cv_cutoff]
    data.to_csv(os.path.join(outdir, 'pca.filtered.data.txt'), header=True, index=True, sep='\t')
    data = np.log(data+1)
    data = data.apply(preprocessing.scale, axis=0)
    data = data.transpose()
    pca = decomposition.PCA()
    pca.fit(data)
    _ratio = list(enumerate(pca.explained_variance_ratio_, start=1))
    total_ratio, n_components = 0, 0
    for ind, each in _ratio:
        total_ratio += each
        if total_ratio >= explained_ratio:
            n_components = ind
            break
    if n_components <= 1:
        n_components = 2
    _ratio = _ratio[:n_components]
    result = pd.DataFrame(pca.transform(data), index=data.index)
    result = result.iloc[:, :n_components]
    result.index.name = 'sample'
    # result.columns = ['PC'+str(n)+'('+'{:.2f}%'.format(r*100)+')' for n, r in _ratio]
    result.columns = ['PC'+str(n) for n in range(1, result.shape[1] + 1)]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    out_dir = os.path.join(outdir, 'PCA.xls')
    result.to_csv(out_dir, sep='\t', header=True, index=True)
    pc_ratio = {'PC' + str(n): r for n, r in _ratio}
    out_dir2 = os.path.join(outdir, 'Explained_variance_ratio.xls')
    with open(out_dir2, 'w') as f:
        for each in sorted(pc_ratio.keys()):
            f.write(str(each) + '\t' + str(pc_ratio[each]) + '\n')

    # color and maker pool
    if type(group_dict) == str:
        if os.path.exists(group_dict):
            with open(group_dict) as f:
                group_dict = dict(x.strip().split()[:2] for x in f if x.strip())
        else:
            raise FileExistsError('group_dict file not existed!')

    if group_dict is None or not group_dict:
        group_dict = OrderedDict()
    for sample in result.index:
        if sample not in group_dict:
            group_dict[sample] = sample
    groups = list()
    for group in group_dict.values():
        if group not in groups:
            groups.append(group)
    colors = get_color_pool(len(group_dict))
    makers = get_marker_pool(len(groups))
    sample_colors = dict(zip(group_dict.keys(), colors))
    group_makers = dict(zip(groups, makers))
    sample_makers = dict()
    for sample in result.index:
        group = group_dict[sample]
        sample_makers[sample] = group_makers[group]

    # plot
    traces = list()
    for sample in result.index:
        trace = go.Scatter(
            x=[result.loc[sample, 'PC1']],
            y=[result.loc[sample, 'PC2']],
            mode='markers',
            name=sample,
            marker=dict(
                color=sample_colors[sample],
                line=dict(color='rgba(217, 217, 217, 0.14)', width=0.5),
                opacity=0.8,
                size=20,
                symbol=sample_makers[sample],
            ),
        )
        traces.append(trace)
    layout = dict(
        showlegend=True,
        # width=800,
        # height=600,
        xaxis=dict(title='PC1({:.2%})'.format(pc_ratio['PC1'])),
        yaxis=dict(title='PC2({:.2%})'.format(pc_ratio['PC2']))
    )
    fig = go.Figure(traces, layout=layout)
    prefix = "PC1_PC2"
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def exp_density(exp_table, outdir=os.getcwd(), row_sum_cutoff=0.1, exp_cutoff=0.1,
                formats=('html',), height:int=None, width:int=None, scale=3):

    def get_density(all_exp_pd):
        """
        sampling 800 density point for each columns of the input pandas DataFrame
        data row with log transformation failed will be removed
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        from scipy import stats
        records = list()
        target_columns = all_exp_pd.columns
        for sample in target_columns:
            exp = all_exp_pd[sample]
            exp = exp.replace([np.inf, -np.inf], np.nan).dropna()
            exp = exp[exp != 0]
            density_func = stats.gaussian_kde(exp)
            min_exp, max_exp = exp.min(), exp.max()
            x_data = np.linspace(min_exp, max_exp, num=800, endpoint=False)
            y_data = density_func(x_data)
            point_df = pd.DataFrame({'exp': x_data, 'density': y_data})
            records.append(point_df)
        return records

    data = pd.read_table(exp_table, header=0, index_col=0)
    data = data[data.sum(axis=1) >= row_sum_cutoff]
    pass_state = data.apply(lambda x: sum(x > exp_cutoff), axis=1)
    data = data[pass_state >= int(data.shape[1]) / 3]
    data.to_csv(os.path.join(outdir, 'density.filtered.data.txt'), header=True, index=True, sep='\t')
    data = np.log(data)
    traces = list()
    density_point_df_list = get_density(data)
    color_pool = get_color_pool(len(data.columns))
    for ind, (sample, color) in enumerate(zip(data.columns, color_pool)):
        # print(sample)
        trace = go.Scatter(
            x=density_point_df_list[ind]['exp'],
            y=density_point_df_list[ind]['density'],
            mode='lines',
            fill='tonexty',
            fillcolor=color,
            name=sample,
            opacity=0.7,
            line=dict(color=color)
        )
        traces.append(trace)

    layout = go.Layout(
        barmode='overlay',
        title='Expression Density',
        xaxis=dict(title='Log(Expression)', zeroline=False),
        yaxis=dict(title='Density', )
    )
    fig = go.Figure(data=traces, layout=layout)
    prefix = "Expression.density"
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def alignment_summary_table(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    summary_dict_list = [json.load(open(x), object_pairs_hook=OrderedDict) for x in files]
    sample_list = [os.path.basename(x).split('.', 1)[0] for x in files]
    df = pd.DataFrame(summary_dict_list, index=sample_list).round(2)
    df.index.name = 'sample'
    out_table = os.path.join(outdir, 'alignment_summary.txt')
    df.to_csv(out_table, index=True, header=True, sep='\t')
    # header = ['<br>'.join(textwrap.wrap(x, width=10)) for x in [df.index.name] + list(df.columns)]
    # df = df.reset_index()
    # df.columns = header
    # colorscale = [[0, '#4d004c'], [.5, '#f2e5ff'], [1, '#ffffff']]
    # table = ff.create_table(df, colorscale=colorscale)
    # fig = go.Figure(data=table)
    # out_name = os.path.join(outdir, 'AlignmentSummaryTable.html')
    # plt(fig, filename=out_name)

    df.reset_index(inplace=True)
    header_values = ['<b>' + '<br>'.join(textwrap.wrap(x, width=10)) + '</b>' for x in list(df.columns)]
    trace = go.Table(
        # columnwidth=[max(len(str(x)) for x in df[y]) for y in df.columns],
        header=dict(values=header_values,
                    fill=dict(color='#C2D4FF'),
                    align=['left'] * df.shape[1]),
        cells=dict(values=[df[x] for x in df.columns],
                   fill=dict(color='#F5F8FF'),
                   align=['left'] * df.shape[1]))
    layout = dict(
        title='Alignment Summary',
        autosize=True,
        margin=dict(t=25, l=10, r=10, b=10),
        showlegend=False,
    )
    fig = go.Figure(data=[trace], layout=layout)
    prefix = "AlignmentSummaryTable"
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def target_region_depth_distribution(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    data_dict = [json.load(open(x), object_pairs_hook=OrderedDict) for x in files]
    sample_list = [os.path.basename(x).split('.', 1)[0] for x in files]
    for distr_dict, sample in zip(data_dict, sample_list):
        layout = go.Layout(
            title="Base depth distribution",
            # xaxis=dict(title='Depth'),
            yaxis=dict(title='Base number ratio', type='log'),
        )
        total = sum(distr_dict.values())
        norm_y = [x/total for x in distr_dict.values()]
        trace = go.Bar(x=list(distr_dict.keys()), y=norm_y, )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = '{}.baseDepthDistribution'.format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def chromosome_read_distribution(files:list, outdir=os.getcwd(), top=30, formats=('html',), height:int=None, width:int=None, scale=3):
    all_data = list()
    samples = list()
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0).iloc[:-1, :]
        data = data.sort_values(by=2, ascending=False).loc[:, [2]].iloc[:top, :]
        total = data[2].sum()
        trace = go.Bar(x=data.index, y=data[2])
        layout = go.Layout(
            title='Read distribution on chromosomes/scaffolds',
            # xaxis=dict(title='Chromosomes/Scaffolds'),
            yaxis=dict(title='Mapped read number')
        )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = '{}.ChromosomeReadDistribution'.format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
        all_data.append(data[2]/total)
        samples.append(sample)
    # overall
    data = pd.concat(all_data, axis=1, sort=False).fillna(0)
    data.columns = samples
    out_table = os.path.join(outdir, 'ChromosomeReadDistribution.txt')
    data.to_csv(out_table, index=True, header=True, sep='\t')
    data = data.transpose()
    color_pool = get_color_pool(data.shape[1])
    traces = [go.Bar(x=data.index, y=data[x], name=x, marker=dict(color=c)) for x, c in zip(data.columns, color_pool)]
    layout = go.Layout(
        title='Read distribution on chromosomes/scaffolds across samples',
        # xaxis=dict(title='Sample'),
        yaxis=dict(title='Mapped read number percent'),
        barmode='stack'
    )
    fig = go.Figure(data=traces, layout=layout)
    prefix = '{}.ChromosomeReadDistribution'.format('samples')
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)


def CollectAlignmentSummaryMetrics(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    if type(files) == str:
        files = glob(files)
    data = list()
    for each in files:
        if not os.path.exists(each):
            pass
        sample = os.path.basename(each).split('.', 1)[0]
        summary = pd.read_table(each, comment='#', index_col=0, header=0)
        summary =summary.loc[["PAIR"], :]
        summary.index = [sample]
        data.append(summary)
    data = pd.concat(data, axis=0).dropna(axis=1).drop('PCT_ADAPTER', axis=1).round(4)
    data = data.transpose()
    out_table = os.path.join(outdir, 'AlignmentSummaryMetrics.xls')
    data.to_csv(out_table, index=True, header=True, sep='\t')

    for i in range(0, data.shape[1], 8):
        df = data.iloc[:, i: i+10]
        df.index.name = 'CATEGORY'
        df.reset_index(inplace=True)
        header_values = ['<b>'+'<br>'.join(textwrap.wrap(x, width=12))+'</b>' for x in list(df.columns)]
        trace = go.Table(
            columnwidth=[max(len(str(x)) for x in df[y]) for y in df.columns],
            header=dict(values=header_values,
                        fill=dict(color='#C2D4FF'),
                        align=['left'] * df.shape[1]),
            cells=dict(values=[df[x] for x in df.columns],
                       fill=dict(color='#F5F8FF'),
                       align=['left'] * df.shape[1]))
        layout = dict(
            title='Alignment Summary',
            autosize=True,
            margin=dict(t=25, l=10, r=10, b=10),
            showlegend=False,
        )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = 'AlignmentSummaryMetrics_{}'.format(i+1)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    return data


def CollectInsertSizeMetrics(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    if type(files) == str:
        files = glob(files)
    data = list()
    for each in files:
        if not os.path.exists(each):
            pass
        # histogram_line = [x[0] for x in enumerate(open(each)) if x[1].startswith('## HISTOGRAM')][0]
        sample = os.path.basename(each).split('.', 1)[0]
        summary = pd.read_table(each, comment='#', header=0, nrows=1)
        summary.index = [sample]
        data.append(summary)
    data = pd.concat(data, axis=0).dropna(axis=1).round(4)
    data = data.transpose()
    out_table = os.path.join(outdir, 'InsertSizeMetrics.xls')
    data.to_csv(out_table, index=True, header=True, sep='\t')

    for i in range(0, data.shape[1], 10):
        df = data.iloc[:, i: i + 10]
        df.index.name = 'CATEGORY'
        df.reset_index(inplace=True)
        header_values = ['<b>' + '<br>'.join(textwrap.wrap(x, width=12)) + '</b>' for x in list(df.columns)]
        trace = go.Table(
            columnwidth=[max(len(str(x)) for x in df[y]) for y in df.columns],
            header=dict(values=header_values,
                        fill=dict(color='#C2D4FF'),
                        align=['left'] * df.shape[1]),
            cells=dict(values=[df[x] for x in df.columns],
                       fill=dict(color='#F5F8FF'),
                       align=['left'] * df.shape[1]))
        layout = dict(
            title='InsertSize Summary',
            autosize=True,
            margin=dict(t=25, l=10, r=10, b=10),
            showlegend=False,
        )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = 'InsertSizeMetrics_{}'.format(i+1)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    return data


def CollectRnaSeqMetrics(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    if type(files) == str:
        files = glob(files)
        if not files:
            raise ValueError('No target file matched!')
    data = list()
    for each in files:
        if not os.path.exists(each):
            pass
        # histogram_line = [x[0] for x in enumerate(open(each)) if x[1].startswith('## HISTOGRAM')][0]
        sample = os.path.basename(each).split('.', 1)[0]
        summary = pd.read_table(each, comment='#', header=0, nrows=1)
        summary.index = [sample]
        data.append(summary)
    data = pd.concat(data, axis=0).dropna(axis=1).round(4)
    data = data.drop(['CORRECT_STRAND_READS', 'INCORRECT_STRAND_READS', 'IGNORED_READS', 'PCT_CORRECT_STRAND_READS'], axis=1)
    data = data.transpose()
    out_table = os.path.join(outdir, 'RnaSeqMetrics.xls')
    data.to_csv(out_table, index=True, header=True, sep='\t')

    for i in range(0, data.shape[1], 10):
        df = data.iloc[:, i: i + 10]
        df.index.name = 'CATEGORY'
        df.reset_index(inplace=True)
        header_values = ['<b>' + '<br>'.join(textwrap.wrap(x, width=12)) + '</b>' for x in list(df.columns)]
        trace = go.Table(
            columnwidth=[max(len(str(x)) for x in df[y]) for y in df.columns],
            header=dict(values=header_values,
                        fill=dict(color='#C2D4FF'),
                        align=['left'] * df.shape[1]),
            cells=dict(values=[df[x] for x in df.columns],
                       fill=dict(color='#F5F8FF'),
                       align=['left'] * df.shape[1]))
        layout = dict(
            title='RNA Summary',
            autosize=True,
            margin=dict(t=25, l=10, r=10, b=10),
            showlegend=False,
        )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = 'RnaSeqMetrics_{}'.format(i + 1)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    return data


def CollectTargetedPcrMetrics(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3):
    if type(files) == str:
        files = glob(files)
        if not files:
            raise ValueError('No target file matched!')
    data = list()
    for each in files:
        if not os.path.exists(each):
            pass
        # histogram_line = [x[0] for x in enumerate(open(each)) if x[1].startswith('## HISTOGRAM')][0]
        sample = os.path.basename(each).split('.', 1)[0]
        summary = pd.read_table(each, comment='#', header=0, nrows=1)
        summary.index = [sample]
        data.append(summary)
    data = pd.concat(data, axis=0).dropna(axis=1).round(4)
    data = data.transpose()
    out_table = os.path.join(outdir, 'TargetedPcrMetrics.xls')
    data.to_csv(out_table, index=True, header=True, sep='\t')

    for i in range(0, data.shape[1], 10):
        df = data.iloc[:, i: i + 10]
        df.index.name = 'CATEGORY'
        df.reset_index(inplace=True)
        header_values = ['<b>' + '<br>'.join(textwrap.wrap(x, width=12)) + '</b>' for x in list(df.columns)]
        trace = go.Table(
            columnwidth=[max(len(str(x)) for x in df[y]) for y in df.columns],
            header=dict(values=header_values,
                        fill=dict(color='#C2D4FF'),
                        align=['left'] * df.shape[1]),
            cells=dict(values=[df[x] for x in df.columns],
                       fill=dict(color='#F5F8FF'),
                       align=['left'] * df.shape[1]))
        layout = dict(
            title='Targeted Alignment Summary',
            autosize=True,
            margin=dict(t=25, l=10, r=10, b=10),
            showlegend=False,
        )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = 'TargetedPcrMetrics_{}'.format(i + 1)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale)
    return data


def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return path


def make_slider(images:list, image_ids:list=None, image_desc:list=None, template="templates/slide.jinja2",
                out='slider.html', link_images=False):
    out_dir = os.path.dirname(out)
    if not out_dir:
        out_dir = os.getcwd()
    else:
        mkdir(out_dir)
    if not images:
        return

    if link_images:
        out_dir = os.path.join(out_dir, os.path.basename(out)[:-5]+'.images')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        tmp_list = list()
        for each in images:
            target = os.path.join(out_dir, os.path.basename(each))
            if os.path.exists(target):
                continue
            os.symlink(each, target)
            tmp_list.append(target)
        images = tmp_list

    if not image_ids:
        image_ids = [os.path.basename(x).split('.')[0] for x in images]
    if len(image_ids) < len(images):
        raise Exception('number of image ids is less than the number of images')

    if not image_desc:
        image_desc = ['']*len(images)
    if len(image_desc) == 1:
        image_desc = image_desc*len(images)
    else:
        if len(image_desc) < len(images):
            raise Exception('number of image desc is less than the number of images')

    # env = Environment(loader=FileSystemLoader("templates"))
    # template = env.get_template("slide.jinja2")
    template = Template(open(template).read())
    img_info_dict = dict()
    img_cls = ["slide"]*len(images)
    img_cls[0] = "slide slide--current"
    for img, img_id, desc, tag in zip(images, image_ids, image_desc, img_cls):
        img_info_dict[img_id] = dict()
        img_info_dict[img_id]['path'] = os.path.relpath(os.path.abspath(img), start=os.path.abspath(out_dir))
        img_info_dict[img_id]['cls'] = tag
        img_info_dict[img_id]['description'] = desc
    # print(img_info_dict)
    content = template.render(img_info_dict=img_info_dict)
    with open(out, 'w', encoding='utf-8') as f:
        f.write(content)
    return out


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                arg_type = sig.parameters[arg].annotation
                if arg == 'self':
                    continue
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
                else:
                    if arg_type in [list, tuple, set] or (type(value) in [list, tuple, set]):
                        default_value = ' '.join(str(x) for x in value)
                        if type(value) in [list, tuple]:
                            one_value = value[0]
                        else:
                            one_value = value.pop()
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(one_value),
                                            metavar='Default:'+default_value, )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value),
                                            metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            try:
                with open("Argument_detail.json", 'w') as f:
                    json.dump(args, f, indent=2, sort_keys=True)
            except IOError:
                print('Current Directory is not writable, thus argument log is not written !')
            start = time.time()
            func(**args)
            print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y)}
    exclude = ['partial', 'OrderedDict', 'glob', 'plt', 'get_color_pool', 'get_marker_pool', 'draw']
    _ = [callable_dict.pop(x) for x in exclude if x in callable_dict]
    Func2Command(callable_dict)
