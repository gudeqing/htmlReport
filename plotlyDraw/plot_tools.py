import os
from functools import partial
from collections import OrderedDict
import pandas as pd
import numpy as np
import json
import textwrap
from plotly import tools
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.offline import plot as plt
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


def gene_body_coverage(files, outdir=os.getcwd()):
    layout = go.Layout(title="geneBodyCoverage")
    all_fig = go.Figure(layout=layout)
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0)
        fig = go.Figure(layout=go.Layout(title='{}'.format(sample)))
        normalize_y = data.iloc[1, :]/data.iloc[1, :].max()
        fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        all_fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        out_name = os.path.join(outdir, '{}.geneBodyCoverage.html'.format(sample))
        plt(fig, filename=out_name)
    out_name = os.path.join(outdir, 'samples.geneBodyCoverage.html')
    plt(all_fig, filename=out_name)


def fragment_length(files, outdir=os.getcwd(), min_len=50, max_len=600):
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
        out_name = os.path.join(outdir, "{}.fragmentLengthDistribution.html".format(sample))
        plt(fig, filename=out_name)
    out_name = os.path.join(outdir, "samples.fragmentLengthDistribution.html")
    plt(all_fig, filename=out_name)


def inner_distance(files, outdir, min_dist=-250, max_dist=250):
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
        out_name = os.path.join(outdir, "{}.InnerDistanceDistribution.html".format(sample))
        plt(fig, filename=out_name)
    html_file = os.path.join(outdir, 'samples.InnerDistanceDistribution.html')
    plt(all_fig, filename=html_file)


def read_distribution(files, outdir):
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
        out_name = os.path.join(outdir, "{}.ReadDistribution.html".format(sample))
        plt(fig, filename=out_name)
    df = pd.concat(all_data, axis=1).T
    # print(df.head())
    data = [go.Bar(x=df.index, y=df[x], name=x) for x in df.columns]
    layout = go.Layout(
        title="Read distribution",
        xaxis=dict(title='Sample'),
        barmode='stack'
    )
    fig = go.Figure(data=data, layout=layout)
    html_file = os.path.join(outdir, 'samples.ReadDistribution.html')
    plt(fig, filename=html_file)


def read_duplication(pos_dup_files, outdir=os.getcwd(), max_dup=500):
    traces = list()
    for each in pos_dup_files:
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
        out_name = os.path.join(outdir, "{}.ReadDuplication.html".format(sample))
        plt(fig, filename=out_name)

    layout = go.Layout(
        title="Read duplication",
        xaxis=dict(title='Occurrence'),
        yaxis=dict(title='UniqReadNumber', type='log'),
    )
    fig = go.Figure(traces, layout=layout)
    out_name = os.path.join(outdir, "samples.ReadDuplication.html")
    plt(fig, filename=out_name)


def exp_saturation(exp_files, outdir=os.getcwd(), outlier_limit=5):
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
    color_pool = get_color_pool(len(exp_files))
    for exp_file, sample_color in zip(exp_files, color_pool):
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
        out_name = os.path.join(outdir, "{}.TPM.Saturation.html".format(sample))
        plt(fig, filename=out_name)

    # plot all sample together
    out_name = os.path.join(outdir, "samples.TPM.Saturation.html")
    plt(all_fig, filename=out_name)


def exp_pca(exp_table, row_sum_cutoff=0.1, exp_cutoff=0.1, cv_cutoff=0.05,
            explained_ratio=0.95, outdir=os.getcwd(), group_dict=None):
    from sklearn import decomposition, preprocessing
    data = pd.read_table(exp_table, header=0, index_col=0)
    data = data[data.sum(axis=1) >= row_sum_cutoff]
    pass_state = data.apply(lambda x: sum(x > exp_cutoff), axis=1)
    data = data[pass_state >= int(data.shape[1])/3]
    data = data[data.std(axis=1)/data.mean(axis=1) > cv_cutoff]
    data.to_csv(os.path.join(outdir, 'pca.filtered.data.txt'), header=True, index=True, sep='\t')
    data = np.log(data+1)
    # data = data.apply(preprocessing.scale, axis=0)
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
    if group_dict is None or not group_dict:
        group_dict = dict()
    for sample in result.index:
        if sample not in group_dict:
            group_dict[sample] = sample
    groups = list(set(group_dict.values()))
    colors = get_color_pool(len(groups))
    makers = range(len(groups))
    group_colors = dict(zip(groups, colors))
    group_makers = dict(zip(groups, makers))
    sample_colors = dict()
    sample_makers = dict()
    for sample in result.index:
        group = group_dict[sample]
        sample_colors[sample] = group_colors[group]
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
        xaxis=dict(title='PC1({:.2%})'.format(pc_ratio['PC1'])),
        yaxis=dict(title='PC2({:.2%})'.format(pc_ratio['PC2']))
    )
    fig = go.Figure(traces, layout=layout)
    out_name = os.path.join(outdir, 'PC1_PC2.html')
    plt(fig, filename=out_name)


def exp_density(exp_table, outdir=os.getcwd(), row_sum_cutoff=0.1, exp_cutoff=0.1):

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
        print(sample)
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
    out_name = os.path.join(outdir, 'Expression.density.html')
    plt(fig, filename=out_name)
    return out_name


def alignment_summary_table(files, outdir=os.getcwd()):
    summary_dict_list = [json.load(open(x), object_pairs_hook=OrderedDict) for x in files]
    sample_list = [os.path.basename(x).split('.', 1)[0] for x in files]
    df = pd.DataFrame(summary_dict_list, index=sample_list).round(2)
    df.index.name = 'sample'
    out_table = os.path.join(outdir, 'alignment_summary.txt')
    df.to_csv(out_table, index=True, header=True, sep='\t')
    header = ['<br>'.join(textwrap.wrap(x, width=10)) for x in [df.index.name] + list(df.columns)]
    df = df.reset_index()
    df.columns = header
    colorscale = [[0, '#4d004c'], [.5, '#f2e5ff'], [1, '#ffffff']]
    table = ff.create_table(df, colorscale=colorscale)
    fig = go.Figure(data=table)
    out_name = os.path.join(outdir, 'AlignmentSummaryTable.html')
    plt(fig, filename=out_name)


def target_region_depth_distribution(files, outdir=os.getcwd()):
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
        out_name = os.path.join(outdir, '{}.baseDepthDistribution.html'.format(sample))
        plt(fig, filename=out_name)


def chromosome_read_distribution(files, outdir=os.getcwd(), top=30):
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
        out_name = os.path.join(outdir, '{}.ChromosomeReadDistribution.html'.format(sample))
        plt(fig, filename=out_name)
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
    out_name = os.path.join(outdir, 'samples.ChromosomeReadDistribution.html')
    plt(fig, filename=out_name)


if __name__ == '__main__':
    def introduce_command(func):
        import argparse
        import inspect
        import json
        import time
        if isinstance(func, type):
            description = func.__init__.__doc__
        else:
            description = func.__doc__
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
        func_args = inspect.getfullargspec(func)
        arg_names = func_args.args
        arg_defaults = func_args.defaults
        arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
        for arg, value in zip(arg_names, arg_defaults):
            if arg == 'self':
                continue
            if value == 'None':
                parser.add_argument('-'+arg, required=True, metavar=arg)
            elif type(value) == bool:
                if value:
                    parser.add_argument('--'+arg, action="store_false", help='default: True')
                else:
                    parser.add_argument('--'+arg, action="store_true", help='default: False')
            elif value is None:
                parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
            else:
                parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
        if func_args.varargs is not None:
            print("warning: *varargs is not supported, and will be neglected! ")
        if func_args.varkw is not None:
            print("warning: **keywords args is not supported, and will be neglected! ")
        args = parser.parse_args().__dict__
        with open("Argument_detail.json", 'w') as f:
            json.dump(args, f, indent=2, sort_keys=True)
        start = time.time()
        func(**args)
        print("total time: {}s".format(time.time() - start))

    import sys
    if sys.argv[1] == 'exp_density':
        sys.argv.remove('exp_density')
        introduce_command(exp_density)
    elif sys.argv[1] == 'exp_pca':
        sys.argv.remove('exp_pca')
        introduce_command(exp_pca)
    else:
        pass