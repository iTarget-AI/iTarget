import sys
from bisect import bisect
import numpy as np
import pandas as pd
from pathlib import Path
from sqlalchemy import column
from tqdm import tqdm
import scipy.io as scio
from collections import defaultdict
import time
import os
import psutil

from copy import copy
from sklearn.utils import shuffle 
import numpy as np
import pandas as pd
from model.metrics import evaluate,reshape_tf2th,to_categorical
from model.model import MultimapCNN, save_model, load_model
from feamap import basemap
from feamap import load_config
from feamap.preparation import drugllmfea_from_local, protllmfea_from_local
import torch.nn.functional as F

prj_path = Path(__file__).parent.resolve().parent.resolve()

class predict():
    def __init__(self, params):
        self.params = params

    def load_data(self, ):
        # Support both user-provided inference inputs and repository datasets with the same interface.
        if self.params.source == 'infer':
            # drugs = pd.read_csv(prj_path / 'data' / 'predict_data' / f'{self.params.source}' / f'{self.params.source}_drugs.csv', index_col=0, low_memory=False)
            # proteins = pd.read_csv(prj_path / 'data' / 'predict_data' / f'{self.params.source}' / f'{self.params.source}_prots.csv', index_col=0, low_memory=False)
            drugs = pd.read_csv(prj_path / 'data' / 'predict_data' / f'{self.params.source}_all-data-merge-drug.csv', index_col=0, low_memory=False)
            proteins = pd.read_csv(prj_path / 'data' / 'predict_data' / f'{self.params.source}_all-data-merge-prot.csv', index_col=0, low_memory=False)
            label = pd.read_csv(prj_path / 'data' / 'predict_data' / f'{self.params.source}_prot-drug-index-source.csv', low_memory=False)
        else:
            drugs = pd.read_csv(prj_path / 'data' / 'original_data' / f'{self.params.source}_all-data-merge-drug.csv', index_col=0, low_memory=False)
            proteins = pd.read_csv(prj_path / 'data' / 'original_data' / f'{self.params.source}_all-data-merge-prot.csv', index_col=0, low_memory=False)
            label = pd.read_csv(prj_path / 'data' / 'original_data' / f'{self.params.source}_prot-drug-index-source.csv', low_memory=False)

        return drugs, proteins, label

    def load_maps(self, disttype = 'uniprot+fullchembl', drugftype = ['xmol',], drugflist=[[],], protftype = ['esm2',], protflist = [[],], split_channels=True, metric = 'cosine', fitmethod = 'umap'):
        # Reload the fitted map templates so inference uses the exact same spatial layout as training.
        mp_d, mp_p = [], []
        for ftype, flist in zip(drugftype,drugflist):
            mp = basemap.Map(disttype = disttype, datatype = 'drug_fea', ftype = ftype, flist = flist, fmap_type='grid', split_channels=split_channels,  metric=metric, var_thr=1e-4)
            mp = mp.load(prj_path / 'data' / 'processed_data' / f'{ftype}_fitted.mp')
            mp_d.append(mp)
        for ftype, flist in zip(protftype,protflist):
            mp = basemap.Map(disttype = disttype, datatype = 'protein_fea', ftype = ftype, flist = flist, fmap_type='grid', split_channels=split_channels,  metric=metric, var_thr=1e-4)
            mp = mp.load(prj_path / 'data' / 'processed_data' / f'{ftype}_fitted.mp')
            mp_p.append(mp)
        
        return mp_d, mp_p

    def produce_fea(self, drugftype, protftype, mp_d, mp_p, drugs, proteins, ):
        # Transform raw LM features into the image-like tensors expected by the classifier.
        fea_drug, fea_prot = [], []
        id2idx_drug, id2idx_prot = dict(), dict()
        for channel, (ftype, mp) in enumerate(zip(drugftype,mp_d)):
            fea_o, _bitsinfo, colormaps_d = drugllmfea_from_local(drugs, feature_dict=None, type=ftype)
            fea_map, ids = mp.batch_transform(fea_o, scale=True, scale_method=self.params.scale_method)
            fea_drug.append(fea_map.astype("float32"))
        
        # feamap shape before reshape: N, W, H, C = {N, W, H, C}
        id2idx_drug = {k:v for v,k in enumerate(ids)}

        for channel, (ftype, mp) in enumerate(zip(protftype,mp_p)):
            fea_o, _bitsinfo, colormaps_p = protllmfea_from_local(proteins, feature_dict=None, type=ftype)
            fea_map, ids = mp.batch_transform(fea_o, scale=True, scale_method=self.params.scale_method)
            fea_prot.append(fea_map.astype("float32"))
            
        id2idx_prot = {k:v for v,k in enumerate(ids)}
        
        print('fea_drug.shape: ', [i.shape for i in fea_drug])
        print('fea_prot.shape: ', [i.shape for i in fea_prot])

        id2idx = dict()
        id2idx.update(id2idx_drug)
        id2idx.update(id2idx_prot)

        # print(id2idx)
        return fea_drug[0], fea_prot[0], id2idx

    def inits(self, label, id2idx, fea_drug, fea_prot):
        # Build the final inference tensor batch according to the pair index table.
        self.save_path = prj_path / 'predict_result' / f'{self.params.source}' / f"{self.params.pretrained_model_path.split('/')[-1]}"

        print(f'source: {self.params.source}')
        label = pd.concat([label.loc[label['source']==src] for src in self.params.source.split(",")]).sort_index()
        label.to_csv(prj_path / 'data' / 'predict_data' / f'{self.params.source}' / f'testdata{self.params.source}_run.csv')

        data_drug_pred = fea_drug[label['drugid'].map(id2idx).values]
        data_prot_pred = fea_prot[label['protid'].map(id2idx).values]

        # reshape for torch 
        print('reshape for torch')
        data_drug_pred = reshape_tf2th(data_drug_pred)
        data_prot_pred = reshape_tf2th(data_prot_pred)
        
        # split your data
        labelX = (data_drug_pred, data_prot_pred)
        labelY = label['label'].values
        # labelY = to_categorical(num_classes = 2, y = label['label'].values)

        return labelX, labelY

    def run(self,):
        # Full inference path: load input tables, build maps, run the model, then export scores.
        drugs, proteins, label = self.load_data()
        mp_d, mp_p = self.load_maps(disttype = self.params.disttype, drugftype = ['xmol'], drugflist=[[]], protftype = ['esm2'], protflist=[[]], split_channels=True, metric = 'cosine', fitmethod = 'umap')
        fea_drug, fea_prot, id2idx = self.produce_fea(['xmol'], ['esm2'], mp_d, mp_p, drugs, proteins, )
        labelX, labelY = self.inits(label, id2idx, fea_drug, fea_prot)


        allfold_label_data, = {},
        assess = ['tn', 'fp', 'fn', 'tp', 'acc', 'auc', 'mcc', 'precision', 'recall', 'specificity', 'sensitivity', 'f1', 'prauc', 'av_prc']
        for fold in range(self.params.kfold_num):
            kfold_label_data, = {},
            for ass in assess:
                kfold_label_data[ass] = []
            # fit your model
            print(f'>>> working on fold {fold} <<<')
            clf = load_model(params = self.params, model_path = os.path.join(prj_path, self.params.pretrained_model_path), in_channels=(fea_drug.shape[-1],fea_prot.shape[-1]), gpuid=self.params.gpu)
            uc_flag = False if self.params.uncertainty_iteration_size == 1 else True
            ROC_savepath = self.save_path / 'ROC_data' / f'{fold}th_fold'
            ROC_savepath.mkdir(parents=True, exist_ok=True)
            PRC_savepath = self.save_path / 'PRC_data' / f'{fold}th_fold'
            PRC_savepath.mkdir(parents=True, exist_ok=True)
            logits_savepath = self.save_path / 'logits_data' / f'{fold}th_fold'
            logits_savepath.mkdir(parents=True, exist_ok=True)
            for i in range(1,1+self.params.uncertainty_iteration_size):
                labelY_pred, _latent_d, _latent_p = clf.run_loop(X=labelX,batch_size=self.params.batch_size, uncertainty_flag=uc_flag)
                # fprs, tprs, thresholds_auc, pres, recs, thresholds_prc, tn, fp, fn, tp, acc, auc, mcc, precision, recall, specificity, sensitivity, f1, prauc, av_prc = evaluate(y_true=to_categorical(num_classes = 2, y=labelY), y_pred=F.softmax(labelY_pred,dim=1))
                fprs, tprs, thresholds_auc, pres, recs, thresholds_prc, tn, fp, fn, tp, acc, auc, mcc, precision, recall, specificity, sensitivity, f1, prauc, av_prc = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                # print(f'-------------------------------- finish {fold} fold cv --------------------------------')
                # print(f'INFERENCE result of {i}th attempt: acc = {acc:.4f}; auc = {auc:.4f}, mcc = {mcc:.4f}, precision = {precision:.4f}, recall = {recall:.4f}, specificity = {specificity:.4f}, sensitivity = {sensitivity:.4f}, f1 = {f1:.4f}, prauc = {prauc:.4f}, av_prc = {av_prc:.4f}')
                print(f"INFERENCE result of {i}th attempt: {label[['protid','drugid']]} -------- {F.softmax(labelY_pred,dim=1).numpy()[:,-1]}")
                for ass in assess:
                    exec(f"kfold_label_data['{ass}'].append({ass})")

                pd.DataFrame.from_dict({'fprs':fprs, 'tprs':tprs, 'thresholds':thresholds_auc}).to_csv(ROC_savepath / f'inference_ROC_for_{fold}th_fold_{i}th_attempt.csv')
                pd.DataFrame.from_dict({'pres':pres, 'recs':recs, 'thresholds':thresholds_prc}).to_csv(PRC_savepath / f'inference_PRC_for_{fold}th_fold_{i}th_attempt.csv')
                pd.DataFrame.from_dict({'y_true':labelY, 'y_pred':labelY_pred}).to_csv(logits_savepath / f'inference_logits_{fold}th_fold_{i}th_attempt.csv')
            if uc_flag:
                # Collapse repeated stochastic predictions into mean probabilities and entropy.
                from .MC_UQ_analyse import mc_uq_calculate
                uncertainty_results = mc_uq_calculate(logits_savepath, self.save_path / 'logits_data' / f'inference_logits_{fold}th_fold_mean4total{self.params.uncertainty_iteration_size}iterations.csv')
                labelY_pred = uncertainty_results['mean_probs']
                labelY_entropy = uncertainty_results['mean_entropy']
            else:
                labelY_pred = F.softmax(labelY_pred,dim=1).numpy()
                labelY_entropy = np.zeros_like(labelY_pred)

            label_inferresult_fold = pd.concat((label, pd.DataFrame({'y_pred':labelY_pred, 'y_entropy':labelY_entropy})), axis=1)
            # label_inferresult_fold.to_csv(prj_path / 'data' / 'predict_data' / f'inference_result_{fold}th_fold_{self.params.source}_prot-drug-index-source.csv')
            label_inferresult_fold.to_csv(prj_path / 'data' / 'predict_data' / f'inference_results.csv')

            allfold_label_data[fold] = kfold_label_data

        ave = self.save(assess, allfold_label_data)
        print(f'-------------------------------- {self.params.kfold_num} folds average result --------------------------------')
        # print(f'AVERAGE_LABEL result: acc = {ave.acc:.4f}; auc = {ave.auc:.4f}, mcc = {ave.mcc:.4f}, precision = {ave.precision:.4f}, recall = {ave.recall:.4f}, specificity = {ave.specificity:.4f}, sensitivity = {ave.sensitivity:.4f}, f1 = {ave.f1:.4f}, prauc = {ave.prauc:.4f}, av_prc = {ave.av_prc:.4f}')
        print(f"INFERENCE result: {label_inferresult_fold}")
        

    def save(self, assess, allfold_label_data,):
        # Persist raw inference bookkeeping even though classification metrics may be empty.
        for fold in allfold_label_data.keys():
            label_result = pd.DataFrame({})
            for ass in assess:
                label_result[ass] = allfold_label_data[fold][ass]
            label_result['Iteration'] = range(1, self.params.uncertainty_iteration_size+1)
            label_result['Foldid'] = fold

            resultdata_savepath = self.save_path / 'result_data' / f'{fold}th_fold'
            resultdata_savepath.mkdir(parents=True, exist_ok=True)

            label_result.to_csv(resultdata_savepath / f'label_result_{fold}th_fold.csv')

        result = pd.DataFrame([])
        for ass in assess:
            for k in range(self.params.kfold_num):
                result.at[k,'Foldid'] = k
                result.at[k,ass] = np.mean(allfold_label_data[k][ass])

        resultdata_savepath = self.save_path / 'result_data'
        resultdata_savepath.mkdir(parents=True, exist_ok=True)
        result.to_csv(resultdata_savepath / f'label_result_mean4total{self.params.uncertainty_iteration_size}iterations.csv')

        return result.mean()
