import os
import glob
from jinja2 import Template, FileSystemLoader, Environment


def slider_report(images=None, image_ids=None, image_desc=None, template="templates/slide.jinja2",
                  exp_to_match_images=None, name='slider.html', description='', link_images=True):

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
        img_info_dict[img_id]['path'] = img
        img_info_dict[img_id]['cls'] = tag
        img_info_dict[img_id]['description'] = desc
    # print(img_info_dict)
    content = template.render(img_info_dict=img_info_dict)
    with open(name, 'w', encoding='utf-8') as f:
        f.write(content)
    return name


def qc_slides(result_dir, slide_template="templates/slide.jinja2", reg="*.pdf"):
    qc_name_desc = {
        "GeneBodyCoverage": "基因覆盖度，正常为钟罩型",
        "ReadDistribution": "read的比对到各个基因区域的比例分布图",
        "InnerDistance": "pair-end read 之间的相对距离，负数表明有测序有overlap",
        "ReadDuplication": "冗余度分析图，正常时左边快速下降，尾部平滑无突起",
        "RPKMSaturation": "基因测序饱和度分析，曲线越快到达最低点表示测序约饱和",
        "FragmentSize": "pair read 长度分布",
    }
    qc_result = (os.path.join(result_dir, x) for x in qc_name_desc)
    reg_list = (os.path.join(x, reg) for x in qc_result)
    result_dict = dict(QualityEvaluation=dict())
    for regexp, qc_name in zip(reg_list, qc_name_desc.keys()):
        slide = slider_report(exp_to_match_images=regexp, name=qc_name+'.html',
                              template=slide_template, description=qc_name_desc[qc_name])
        if slide:
            result_dict['QualityEvaluation'][qc_name]=slide
    return result_dict


def make_report(result_dict, index_template="templates/index.jinja2", name='index.html'):
    template = Template(open(index_template).read())
    content = template.render(result=result_dict)
    with open(name, 'w', encoding='utf-8') as f:
        f.write(content)
    return name


if __name__ == '__main__':
    import sys
    result = qc_slides(sys.argv[1])
    make_report(result)
    # slider_report(exp_to_match_images="images/*pdf", name='test.html', link_images=False)

