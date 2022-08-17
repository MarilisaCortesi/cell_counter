import pickle

import matplotlib.pyplot as plt
import skimage
import os
import numpy as np
import pandas
from matplotlib.backends.backend_pdf import PdfPages

def segment_nuclei(img_dapi):
    #f,ax = skimage.filters.try_all_threshold(img_dapi)
    #plt.show()
    th = skimage.filters.threshold_otsu(img_dapi)
    return img_dapi > th

def filter_gfp(lab_n, gfp_im):
    out = np.zeros_like(lab_n)
    out2 = np.zeros_like(lab_n)
    for i in range(1, np.amax(lab_n)):
        pixels = np.where(lab_n == i)
        min_x = min(pixels[0])
        max_x = max(pixels[0])
        min_y = min(pixels[1])
        max_y = max(pixels[1])
        bb = gfp_im[min_x:max_x, min_y:max_y].astype('int_')
        if max_x > min_x and max_y > min_y:
            for ip in range(len(pixels[0])):
                x_b = pixels[0][ip]-min_x
                y_b = pixels[1][ip]- min_y
                if x_b > 0:
                    x_b -=1
                if y_b > 0:
                    y_b-=1
                bb[x_b, y_b] = -1
            non_cell = bb[np.where(bb!=-1)]
            gfp_levels = gfp_im[pixels]
            if (np.mean(gfp_levels) - np.mean(non_cell))/np.mean(non_cell) > 0:
                out[pixels] = 1
            else:
                out2[pixels] = 1
    return out, out2

def filter_area(img, min_area = 50, max_area = 5000):
    out = np.zeros_like(img)
    for i in range(1, np.amax(img)):
        pixels = np.where(img == i)
        area = len(pixels[0])
        if min_area < area < max_area:
            out[pixels] = 1
    return out

def save_to_excel(ncells, file_name):
    df = pandas.DataFrame.from_dict(ncells)
    df.to_excel(file_name)


folder = r'PATH_TO_INPUT_IMAGES_FOLDER'
PDF_out = PdfPages(r'PATH_TO_OUTPUT_FILE.pdf')
file_out = r'PATH_TO_OUTPUT_FILE.xlsx'
folder_out = r'PATH_TO_FOLDER_WHERE_TO_SAVE_MASKS'
kind_experiment = '3D'
files = os.listdir(folder)
ncells = {}
for f in files:
    if 'DAPI' in f and not f.startswith('.') and 'counted' not in f:
    #if not f.startswith('.'):
        image_dapi = skimage.io.imread(folder + '\\' + f)
        name_gfp = f.split('DAPI')[0]+'GFP.tif'
        image_gfp = skimage.io.imread(folder + '\\' + name_gfp)

        nuclei = segment_nuclei(image_dapi[:, :, 2])
        label_nuclei = skimage.measure.label(nuclei)
        if kind_experiment =='3D':
            fluo, non_fluo = filter_gfp(label_nuclei, image_gfp[:, :, 1])
            label_fluo = skimage.measure.label(fluo)
            label_non_fluo = skimage.measure.label(non_fluo)
            cancer = filter_area(label_fluo)
            non_cancer = filter_area(label_non_fluo)
            label_non_cancer = skimage.measure.label(non_cancer)
            print('s')
        else:
            cancer = filter_area(label_nuclei)
        label_cancer = skimage.measure.label(cancer)
        ncells[f] = [np.amax(label_cancer)]
        image_gt = np.add(image_dapi, image_gfp).astype(np.uint8)
        image_out = skimage.color.label2rgb(label_cancer, image=image_gt)
        out_variable = {'image': image_dapi, 'mask': cancer, 'label': label_cancer}
        if kind_experiment == '3D':
            out_variable['label_non_cancer'] = label_non_cancer
            out_variable['mask_non_cancer'] = non_cancer
        with open(folder_out+'\\'+f.split('.tif')[0]+'.npy', 'wb') as F:
            pickle.dump(out_variable,F)
        #image_to_save = np.concatenate((image_gt, image_out))
        if kind_experiment == '3D':
            fg, ax = plt.subplots(2, 2)
            ax[0,0].imshow(image_dapi)
            ax[0,1].imshow((image_gfp- np.amin(image_gfp))/np.amax(image_gfp))
            ax[1, 0].imshow(label_nuclei)
            ax[1,1].imshow(image_out)
            fg.suptitle(f)

        else:
            fg, ax = plt.subplots(2,1)
            ax[0].imshow(image_dapi)
            ax[1].imshow(image_out)
        PDF_out.savefig()
        #skimage.io.imsave(folder_out + '\\' +  f.split('DAPI')[0] + 'counted.tif', 255*image_to_save)
PDF_out.close()
save_to_excel(ncells, file_out)
