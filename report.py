import os
import glob
from jinja2 import FileSystemLoader,Environment


env = Environment(loader=FileSystemLoader(os.getcwd()))    # 创建一个包加载器对象

template = env.get_template("slideTemplate.jinja2")    # 获取一个模板文件

imgs = glob.glob('images/*html')
img_info_dict = dict()
for ind, each in enumerate(imgs):
    key = os.path.basename(each).split('.')[0]
    img_info_dict[key] = dict()
    img_info_dict[key]['path'] = each
    if ind == 0:
        print(key)
        img_info_dict[key]['cls'] = "slide slide--current"
    else:
        img_info_dict[key]['cls'] = "slide"
    img_info_dict[key]['description'] = "This is desc!!!"


a = template.render(img_info_dict=img_info_dict)
with open('new.html', 'w', encoding='utf-8') as f:
    f.write(a)
