import os
import glob
from jinja2 import Template
import importlib
import shutil
plotting = importlib.import_module('.plot_tools', package='plotlyDraw')


def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return path


def slider_report(images=None, image_ids=None, image_desc=None, template="templates/slide.jinja2",
                  exp_to_match_images=None, name='slider.html', description='', link_images=False):
    out_dir = os.path.dirname(name)
    mkdir(out_dir)
    if exp_to_match_images:
        images = glob.glob(exp_to_match_images)
        if not images:
            return
        if link_images:
            out_dir = os.path.dirname(name)
            out_dir = os.path.join(out_dir, os.path.basename(name)[:-5]+'.images')
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
        image_ids = (os.path.basename(x).split('.')[0] for x in images)
        image_desc = [description]*len(images)
    else:
        if not images:
            exit('Either images or exp_to_match_images must be provided')
        if not image_ids:
            image_ids = (os.path.basename(x).split('.')[0] for x in images)
        if not image_desc:
            image_desc = [description]*len(images)
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
    with open(name, 'w', encoding='utf-8') as f:
        f.write(content)
    return name


def qc_plot(result_dir, report_dir):
    qc_name_desc = {
        "GeneBodyCoverage": "基因覆盖度，正常为钟罩型",
        "ReadDistribution": "read的比对到各个基因区域的比例分布图",
        "InnerDistance": "pair-end read 之间的相对距离，负数表明有测序有overlap",
        "ReadDuplication": "冗余度分析图，正常时左边快速下降，尾部平滑无突起",
        "TPMSaturation": "基因测序饱和度分析，曲线越快到达最低点表示测序越饱和",
        "FragmentSize": "pair read 长度分布",
    }
    for qc_name in qc_name_desc:
        outdir = os.path.join(report_dir, qc_name)
        if os.path.exists(os.path.join(result_dir, qc_name)):
            if qc_name == 'GeneBodyCoverage':
                reg = '*.geneBodyCoverage.txt'
                stat_files = glob.glob(os.path.join(result_dir, qc_name, reg))
                if stat_files:
                    mkdir(outdir)
                    plotting.gene_body_coverage(stat_files, outdir=outdir)
            elif qc_name == 'ReadDistribution':
                reg = '*.read_distribution.txt'
                stat_files = glob.glob(os.path.join(result_dir, qc_name, reg))
                if stat_files:
                    mkdir(outdir)
                    plotting.read_distribution(stat_files, outdir=outdir)
            elif qc_name == 'InnerDistance':
                reg = '*.inner_distance.txt'
                stat_files = glob.glob(os.path.join(result_dir, qc_name, reg))
                if stat_files:
                    mkdir(outdir)
                    plotting.inner_distance(stat_files, outdir=outdir, min_dist=-250, max_dist=250)
            elif qc_name == 'ReadDuplication':
                reg = '*.pos.DupRate.xls'
                stat_files = glob.glob(os.path.join(result_dir, qc_name, reg))
                if stat_files:
                    mkdir(outdir)
                    plotting.read_duplication(stat_files, outdir=outdir, max_dup=500)
            elif qc_name == 'TPMSaturation':
                reg = '*.tpm.xls'
                stat_files = glob.glob(os.path.join(result_dir, qc_name, reg))
                if stat_files:
                    mkdir(outdir)
                    plotting.exp_saturation(stat_files, outdir=outdir, outlier_limit=5)
            elif qc_name == 'FragmentSize':
                reg = '*.fragment_size.txt'
                stat_files = glob.glob(os.path.join(result_dir, qc_name, reg))
                if stat_files:
                    mkdir(outdir)
                    plotting.fragment_length(stat_files, outdir=outdir, min_len=50, max_len=600)
    return report_dir


def exp_analysis_plot(result_dir, report_dir, level='gene'):
    if level == 'gene':
        exp_table = os.path.join(result_dir, 'MergeQuantGene', 'gene.isoform.TMM.EXPR.matrix')
    else:
        exp_table = os.path.join(result_dir, 'MergeQuantTranscript', 'transcript.isoform.TMM.EXPR.matrix')

    outdir = os.path.join(report_dir, level.capitalize()+'ExpAnalysis')
    if os.path.exists(exp_table):
        mkdir(outdir)
        plotting.exp_pca(exp_table, row_sum_cutoff=0.1, exp_cutoff=0.1, cv_cutoff=0.05,
                         explained_ratio=0.95, outdir=outdir, group_dict=None)
        plotting.exp_density(exp_table, outdir=outdir, row_sum_cutoff=0.1, exp_cutoff=0.1)
        cluster = importlib.import_module('.clusterheatmap', package='plotlyDraw')
        cluster.ClusterHeatMap(
            data_file=exp_table,
            cluster_sample=True,
            cluster_gene=True,
            do_correlation_cluster=True,
            out_name=os.path.join(outdir, 'corr_cluster_heatmap.html'),
            corr_method='pearson',
            height=600, width=600,
        )
        return outdir


def qc_slides(report_dir, slide_template="templates/slide.jinja2", reg="*.html"):
    qc_name_desc = {
        "GeneBodyCoverage": "基因覆盖度，正常为钟罩型",
        "ReadDistribution": "read的比对到各个基因区域的比例分布图",
        "InnerDistance": "pair-end read 之间的相对距离，负数表明有测序有overlap",
        "ReadDuplication": "冗余度分析图，正常时左边快速下降，尾部平滑无突起",
        "TPMSaturation": "基因测序饱和度分析，曲线越快到达最低点表示测序越饱和",
        "FragmentSize": "pair read 长度分布",
    }
    qc_result = (os.path.join(report_dir, x) for x in qc_name_desc)
    reg_list = (os.path.join(x, reg) for x in qc_result)
    result_dict = dict(QualityEvaluation=dict())
    for regexp, qc_name in zip(reg_list, qc_name_desc.keys()):
        slide = slider_report(exp_to_match_images=regexp,
                              name=os.path.join(report_dir, 'htmls', qc_name+'.html'),
                              template=slide_template, description=qc_name_desc[qc_name])
        if slide:
            result_dict['QualityEvaluation'][qc_name] = slide
    return result_dict


def exp_analysis_slides(report_dir, level='gene', slide_template="templates/slide.jinja2"):
    step_desc = {
        "ExpressionDistribution": "表达分布密度图",
        "CorrelationCluster": "基于相关性作为距离度量的样本聚类分析，用于查看样本重复性和离群样本",
        "PCA": "主成份分析，用于查看样本重复性",
    }
    name = level.capitalize()+'ExpressionAnalysis'
    result_dict = {name: dict()}

    result_dict[name]['ExpressionDistribution'] = slider_report(
        images=[os.path.join(report_dir, level.capitalize()+'ExpAnalysis', 'Expression.density.html')],
        image_ids=["Expression distribution"],
        image_desc=[step_desc['ExpressionDistribution']],
        template=slide_template,
        name=os.path.join(report_dir, 'htmls', 'ExpressionDistribution.html')
    )

    result_dict[name]['CorrelationCluster'] = slider_report(
        images=[os.path.join(report_dir, level.capitalize() + 'ExpAnalysis', 'corr_cluster_heatmap.html')],
        image_ids=["Sample correlation"],
        image_desc=[step_desc['CorrelationCluster']],
        template=slide_template,
        name=os.path.join(report_dir, 'htmls', 'CorrelationCluster.html')
    )

    result_dict[name]['PCA'] = slider_report(
        images=[os.path.join(report_dir, level.capitalize() + 'ExpAnalysis', 'PC1_PC2.html')],
        image_ids=["Principal components analysis"],
        image_desc=[step_desc['PCA']],
        template=slide_template,
        name=os.path.join(report_dir, 'htmls', 'PCA.html')
    )
    return result_dict


def alignment_summary_plot(result_dir, report_dir):
    outdir = os.path.join(report_dir, 'AlignmentSummary')
    mkdir(outdir)
    summary_files = glob.glob(os.path.join(result_dir, 'AlignmentSummary', '*.alignment_summary.json'))
    plotting.alignment_summary_table(summary_files, outdir=outdir)
    # depth distr
    depth_distr = glob.glob(os.path.join(result_dir, 'AlignmentSummary', '*.depth.distribution.json'))
    plotting.target_region_depth_distribution(depth_distr, outdir=outdir)
    # chromosome_read_distribution
    chromosome_distr = glob.glob(os.path.join(result_dir, 'AlignmentSummary', '*.chromosome.alignment.stat.txt'))
    plotting.chromosome_read_distribution(chromosome_distr, outdir=outdir, top=30)


def alignment_summary_slides(report_dir, slide_template="templates/slide.jinja2"):
    step_desc = {
        "Summary": "比对结果统计表",
        "DepthDistribution": "捕获区域的测序深度分布图",
        "ChrReadDistribution": "比对到各染色体或scaffold的read分布图",
    }
    name = 'AlignmentSummary'
    result_dict = {name: dict()}

    result_dict[name]['Summary'] = slider_report(
        images=[os.path.join(report_dir, 'AlignmentSummary', 'AlignmentSummaryTable.html')],
        image_ids=["Alignment summary table"],
        image_desc=[step_desc['Summary']],
        template=slide_template,
        name=os.path.join(report_dir, 'htmls', 'AlignmentSummaryTable.html')
    )

    result_dict[name]['DepthDistribution'] = slider_report(
        exp_to_match_images=os.path.join(report_dir, 'AlignmentSummary', '*.baseDepthDistribution.html'),
        description=step_desc['DepthDistribution'],
        template=slide_template,
        name=os.path.join(report_dir, 'htmls', 'BaseDepthDistribution.html')
    )

    result_dict[name]['ChrReadDistribution'] = slider_report(
        exp_to_match_images=os.path.join(report_dir, 'AlignmentSummary', '*.ChromosomeReadDistribution.html'),
        description=step_desc['ChrReadDistribution'],
        template=slide_template,
        name=os.path.join(report_dir, 'htmls', 'ChromosomeReadDistribution.html')
    )
    return result_dict


def mapping_summary_plot(result_dir, report_dir):
    outdir = os.path.join(report_dir, 'MappingSummary')
    mkdir(outdir)
    target_dirs = [
        'CollectAlignmentSummaryMetrics',
        'CollectInsertSizeMetrics',
        'CollectRnaSeqMetrics',
        'CollectTargetedPcrMetrics'
    ]
    for each in target_dirs:
        files = glob.glob(os.path.join(result_dir, each, '*.Collect*Metrics.xls'))
        if files:
            eval('plotting.{}(files, outdir=outdir)'.format(each))


def mapping_summary_slides(report_dir, slide_template="templates/slide.jinja2"):
    step_desc = {
        "AlignmentSummaryMetrics": "比对结果统计表",
        "InsertSizeMetrics": "Insert Size 统计表",
        "RnaSeqMetrics": "比对结果统计表，侧重基因不同区域的统计，比对到rRNA的统计",
        "TargetedSummaryMetrics": "比对到捕获区域的统计表",
    }
    name = 'MappingSummary'
    result_dict = {name: dict()}
    #
    for step, desc in step_desc.items():
        result_dict[name][step] = slider_report(
            exp_to_match_images=os.path.join(report_dir, 'MappingSummary', '*Metrics*.html'),
            description=desc,
            template=slide_template,
            name=os.path.join(report_dir, 'htmls', '{}.html'.format(step))
        )
    return result_dict


def make_report(result_dir):
    report_dir = mkdir('Report')
    index_html_path = os.path.join(report_dir, 'index.html')
    html_utils_dir = os.path.dirname(__file__)
    shutil.copytree(os.path.join(html_utils_dir, 'html.utils'), os.path.join(report_dir, 'html.utils'))
    slide_template = os.path.join(html_utils_dir, 'templates', 'slide.jinja2')
    index_template = os.path.join(html_utils_dir, 'templates', 'index.jinja2')

    qc_plot(result_dir, report_dir)
    exp_analysis_plot(result_dir, report_dir, level='gene')

    result_dict = dict()
    result_dict.update(qc_slides(report_dir, slide_template=slide_template, reg='*.html'))
    result_dict.update(exp_analysis_slides(report_dir, slide_template=slide_template))
    alignment_summary_plot(result_dir, report_dir)
    result_dict.update(alignment_summary_slides(report_dir, slide_template=slide_template))
    mapping_summary_plot(result_dir, report_dir)
    result_dict.update(mapping_summary_slides(report_dir, slide_template))

    # use relative path, this is important
    for section in result_dict:
        for each in result_dict[section]:
            # print(section)
            tmp_path = os.path.abspath(result_dict[section][each])
            result_dict[section][each] = os.path.relpath(tmp_path, start=os.path.abspath(report_dir))

    template = Template(open(index_template).read())
    content = template.render(result=result_dict)
    with open(index_html_path, 'w', encoding='utf-8') as f:
        f.write(content)


if __name__ == '__main__':
    import sys
    make_report(sys.argv[1])

