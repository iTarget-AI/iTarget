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
from model.model import MultimapCNN, MultimapCNN_dataset, EarlyStopping, save_model, load_model
import torch as th
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F

prj_path = Path(__file__).parent.resolve().parent.resolve()

class cross_valid():
    def __init__(self, params):
        self.params = params
        self.patience = 10

    def load_cvdata(self):
        # BindingDB uses its official split, while other benchmarks are loaded fold by fold.
        cv_data = defaultdict(list)
        if self.params.source == 'bindingdb':
            self.params.kfold_num = 1
            train_k = pd.read_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' / f'{self.params.source}' / f'{self.params.source}_train-index.csv', index_col=0, header=0, low_memory=False)
            valid_k = pd.read_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' / f'{self.params.source}' / f'{self.params.source}_dev-index.csv', index_col=0, header=0, low_memory=False)
            test_k = pd.read_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' / f'{self.params.source}' / f'{self.params.source}_test-index.csv', index_col=0, header=0, low_memory=False)
            cv_data[0] = (train_k, valid_k, test_k)
        else:
            for k in range(self.params.kfold_num):
                train_k = pd.read_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' /  f'{self.params.source}' / f'{k}th_fold' / 'train_k.csv', index_col=0, header=0, low_memory=False)
                valid_k = pd.read_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' /  f'{self.params.source}' / f'{k}th_fold' / 'valid_k.csv', index_col=0, header=0, low_memory=False)
                test_k = pd.read_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' /  f'{self.params.source}' / f'{k}th_fold' / 'test_k.csv', index_col=0, header=0, low_memory=False)
                if self.params.kfold_num==1:
                    train_k = train_k.append(valid_k,ignore_index=True).append(test_k,ignore_index=True)
                cv_data[k] = (train_k, valid_k, test_k)
        return cv_data
    
    def load_fea(self):
        # Load precomputed image-like feature maps and the ID lookup tables for both modalities.
        fea_drug = np.load(prj_path / 'data' / 'processed_data' / 'drug_fea' / 'map_transferred' / 'drug_fea.npy').astype("float32")
        fea_prot = np.load(prj_path / 'data' / 'processed_data' / 'protein_fea' / 'map_transferred' / 'protein_fea.npy').astype("float32")
        print('fea_drug.shape: ', fea_drug.shape)
        print('fea_prot.shape: ', fea_prot.shape)
        id2idx_drug = pd.read_csv(prj_path / 'data' / 'processed_data' / 'drug_fea' / 'map_transferred' / 'drug_list.csv', index_col=1).iloc[:,0].to_dict()
        id2idx_prot = pd.read_csv(prj_path / 'data' / 'processed_data' / 'protein_fea' / 'map_transferred' / 'prot_list.csv', index_col=1).iloc[:,0].to_dict()
        id2idx = dict()
        id2idx.update(id2idx_prot)
        id2idx.update(id2idx_drug)
        # print(id2idx)
        return fea_drug, fea_prot, id2idx

    def inits(self, fold, cv_data, id2idx, fea_drug, fea_prot):
        # Materialize one fold into tensors and persist the exact split used in this run.
        self.save_path = prj_path / 'pretrained' / f'{self.params.kfold_num}_fold_trainval' / f'batchsize_{self.params.batch_size}' / f'learningrate_{self.params.lr}' / f'monitor_{self.params.monitor}'
        
        if self.params.source == 'bindingdb':
            (train, valid, test) = cv_data[fold]
            print(f'source: {self.params.source}')
            train = pd.concat([train.loc[train['source']==src] for src in self.params.source.split(",")]).sort_index()
            valid = pd.concat([valid.loc[valid['source']==src] for src in self.params.source.split(",")]).sort_index()
            test = pd.concat([test.loc[test['source']==src] for src in self.params.source.split(",")]).sort_index()
            train.to_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' / f'{self.params.source}' / 'train_cv_k.csv')
            valid.to_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' / f'{self.params.source}' / 'valid_cv_k.csv')
            test.to_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' / f'{self.params.source}' / 'test_cv_k.csv')
        else:
            (train, valid, test) = cv_data[fold]
            print(f'source: {self.params.source}')
            train = pd.concat([train.loc[train['source']==src] for src in self.params.source.split(",")]).sort_index()
            valid = pd.concat([valid.loc[valid['source']==src] for src in self.params.source.split(",")]).sort_index()
            test = pd.concat([test.loc[test['source']==src] for src in self.params.source.split(",")]).sort_index()
            train.to_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' /  f'{self.params.source}' / f'{fold}th_fold' / 'train_cv_k.csv')
            valid.to_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' /  f'{self.params.source}' / f'{fold}th_fold' / 'valid_cv_k.csv')
            test.to_csv(prj_path / 'data' / 'processed_data' / 'split_cvdata' /  f'{self.params.source}' / f'{fold}th_fold' / 'test_cv_k.csv')

        data_drug_train = fea_drug[train['drugid'].map(id2idx).values]
        data_prot_train = fea_prot[train['protid'].map(id2idx).values]
        data_drug_valid = fea_drug[valid['drugid'].map(id2idx).values]
        data_prot_valid = fea_prot[valid['protid'].map(id2idx).values]
        data_drug_test = fea_drug[test['drugid'].map(id2idx).values]
        data_prot_test = fea_prot[test['protid'].map(id2idx).values]
        # reshape for torch 
        print('reshape for torch')
        data_drug_train = reshape_tf2th(data_drug_train)
        data_prot_train = reshape_tf2th(data_prot_train)
        data_drug_valid = reshape_tf2th(data_drug_valid)
        data_prot_valid = reshape_tf2th(data_prot_valid)
        data_drug_test = reshape_tf2th(data_drug_test)
        data_prot_test = reshape_tf2th(data_prot_test)
        
        # split your data
        trainX = (data_drug_train, data_prot_train)
        trainY = train['label'].values
        validX = (data_drug_valid, data_prot_valid)
        validY = valid['label'].values
        testX = (data_drug_test, data_prot_test)
        testY = test['label'].values

        return trainX, trainY, validX, validY, testX, testY

    def fit(self, model, save_path, trainX, trainY, validX, validY):
        print(model)
        # Build dataloaders once per fold, then optimize until early stopping triggers.
        training_data = MultimapCNN_dataset(trainX, trainY)
        valid_data = MultimapCNN_dataset(validX, validY)
        train_dataloader = DataLoader(training_data, batch_size=self.params.batch_size, shuffle=True)
        valid_dataloader = DataLoader(valid_data, batch_size=self.params.batch_size, shuffle=True)
        optimizer = th.optim.Adam(model.parameters(), lr=self.params.lr, weight_decay=1e-4)

        early_stopping = EarlyStopping(patience=self.patience, verbose=True, monitor = self.params.monitor)
        for t in range(self.params.n_epochs):
            time_epstart = time.time()
            print(f"-------------------------------\nEpoch {t+1}\n-------------------------------")
            train_loss, train_logits, train_label = model.train_loop(train_dataloader, optimizer)
            valid_loss, valid_logits, valid_label = model.test_loop(valid_dataloader, )
            # train_fprs, train_tprs, train_thresholds_auc, train_pres, train_recs, train_thresholds_prc, train_tn, train_fp, train_fn, train_tp, train_acc, train_auc, train_mcc, train_precision, train_recall, train_specificity, train_sensitivity, train_f1, train_prauc, train_av_prc = evaluate(y_true=to_categorical(num_classes = 2, y=train_label), y_pred=F.softmax(train_logits,dim=1))
            valid_fprs, valid_tprs, valid_thresholds_auc, valid_pres, valid_recs, valid_thresholds_prc, valid_tn, valid_fp, valid_fn, valid_tp, valid_acc, valid_auc, valid_mcc, valid_precision, valid_recall, valid_specificity, valid_sensitivity, valid_f1, valid_prauc, valid_av_prc = evaluate(y_true=to_categorical(num_classes = 2, y=valid_label), y_pred=F.softmax(valid_logits,dim=1))
            print(f'Epoch {t+1} result: valid_loss={valid_loss:.4f}, valid_acc={valid_acc:.4f}, valid_auc={valid_auc:.4f}, valid_mcc={valid_mcc:.4f}, valid_precision={valid_precision:.4f}, valid_recall={valid_recall:.4f}, valid_specificity={valid_specificity:.4f}, valid_sensitivity={valid_sensitivity:.4f}, valid_f1={valid_f1:.4f}, valid_prauc={valid_prauc:.4f}')
            print('memory used by this process', f'{(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024):.2f}', 'GB')
            print('time for running this epoch: ', f'{(time.time()-time_epstart):.4f}', 'seconds')
            early_stopping(score = {'loss_val':valid_loss, 'acc_val':valid_acc, 'auc_val':valid_auc, 'aupr_val':valid_prauc, 'mcc_val':valid_mcc, 'f1_val':valid_f1, 'recall_val':valid_recall, 'precision_val':valid_precision, 'specificity_val':valid_specificity}, model = model, model_path = save_path)
            if early_stopping.early_stop:
                print("Early stopping")
                break
        # return early_stopping._model

    def run(self,):
        # Train one model per fold, then collect repeated evaluation results if uncertainty is enabled.
        cv_data = self.load_cvdata()
        fea_drug, fea_prot, id2idx = self.load_fea()
        allfold_train_data, allfold_val_data, allfold_test_data = {},{},{}
        assess = ['tn', 'fp', 'fn', 'tp', 'acc', 'auc', 'mcc', 'precision', 'recall', 'specificity', 'sensitivity', 'f1', 'prauc', 'av_prc']
        
        for fold in cv_data.keys():
            kfold_train_data, kfold_val_data, kfold_test_data = {},{},{}
            for ass in assess:
                kfold_train_data[ass] = []
                kfold_val_data[ass] = []
                kfold_test_data[ass] = []
            
            trainX, trainY, validX, validY, testX, testY = self.inits(fold, cv_data, id2idx, fea_drug, fea_prot)
            # save model path
            save_path_model = self.save_path / 'model' / f'{fold}th_fold'
            save_path_model.mkdir(parents=True, exist_ok=True)
            # fit your model
            print(f'>>> working on fold {fold} <<<')
            clf = MultimapCNN(self.params, in_channels=(fea_drug.shape[-1],fea_prot.shape[-1]))
            self.fit(clf, save_path_model, trainX, trainY, validX, validY)
            # fit finished
            # evaluate your model
            es_model = load_model(params=self.params, model_path=save_path_model, in_channels=(fea_drug.shape[-1],fea_prot.shape[-1]) , gpuid=self.params.gpu)
            uc_flag = False if self.params.uncertainty_iteration_size == 1 else True
            ROC_savepath = self.save_path / 'ROC_data' / f'{fold}th_fold'
            ROC_savepath.mkdir(parents=True, exist_ok=True)
            PRC_savepath = self.save_path / 'PRC_data' / f'{fold}th_fold'
            PRC_savepath.mkdir(parents=True, exist_ok=True)
            logits_savepath = self.save_path / 'logits_data' / f'{fold}th_fold'
            logits_savepath.mkdir(parents=True, exist_ok=True)
            for i in range(1,1+self.params.uncertainty_iteration_size):
                trainY_pred, _latents_d, _latents_p = es_model.run_loop(X=trainX,batch_size=self.params.batch_size, uncertainty_flag=uc_flag)
                fprs, tprs, thresholds_auc, pres, recs, thresholds_prc, tn, fp, fn, tp, acc, auc, mcc, precision, recall, specificity, sensitivity, f1, prauc, av_prc = evaluate(y_true=to_categorical(num_classes = 2, y=trainY), y_pred=F.softmax(trainY_pred,dim=1))
                validY_pred, _latents_d, _latents_p = es_model.run_loop(X=validX,batch_size=self.params.batch_size, uncertainty_flag=uc_flag)
                fprs_val, tprs_val, thresholds_auc_val, pres_val, recs_val, thresholds_prc_val, tn_val, fp_val, fn_val, tp_val, acc_val, auc_val, mcc_val, precision_val, recall_val, specificity_val, sensitivity_val, f1_val, prauc_val, av_prc_val = evaluate(y_true=to_categorical(num_classes = 2, y=validY), y_pred=F.softmax(validY_pred,dim=1))
                testY_pred, _latents_d, _latents_p = es_model.run_loop(X=testX,batch_size=self.params.batch_size, uncertainty_flag=uc_flag)
                fprs_test, tprs_test, thresholds_auc_test, pres_test, recs_test, thresholds_prc_test, tn_test, fp_test, fn_test, tp_test, acc_test, auc_test, mcc_test, precision_test, recall_test, specificity_test, sensitivity_test, f1_test, prauc_test, av_prc_test = evaluate(y_true=to_categorical(num_classes = 2, y=testY), y_pred=F.softmax(testY_pred,dim=1))
                print(f'-------------------------------- finish {fold} fold cv --------------------------------')
                print(f'TRAIN result of {i}th attempt: acc = {acc:.4f}; auc = {auc:.4f}, mcc = {mcc:.4f}, precision = {precision:.4f}, recall = {recall:.4f}, specificity = {specificity:.4f}, sensitivity = {sensitivity:.4f}, f1 = {f1:.4f}, prauc = {prauc:.4f}, av_prc = {av_prc:.4f}')

                print(f'VALID result of {i}th attempt: acc = {acc_val:.4f}; auc = {auc_val:.4f}, mcc = {mcc_val:.4f}, precision = {precision_val:.4f}, recall = {recall_val:.4f}, specificity = {specificity_val:.4f}, sensitivity = {sensitivity_val:.4f}, f1 = {f1_val:.4f}, prauc = {prauc_val:.4f}, av_prc = {av_prc_val:.4f}')

                print(f'TEST result of {i}th attempt: acc = {acc_test:.4f}; auc = {auc_test:.4f}, mcc = {mcc_test:.4f}, precision = {precision_test:.4f}, recall = {recall_test:.4f}, specificity = {specificity_test:.4f}, sensitivity = {sensitivity_test:.4f}, f1 = {f1_test:.4f}, prauc = {prauc_test:.4f}, av_prc = {av_prc_test:.4f}')

                for ass in assess:
                    exec(f"kfold_train_data['{ass}'].append({ass})")
                for ass in assess:
                    exec(f"kfold_val_data['{ass}'].append({ass}_val)")
                for ass in assess:
                    exec(f"kfold_test_data['{ass}'].append({ass}_test)")

                # pd.DataFrame.from_dict({'fprs':fprs, 'tprs':tprs, 'thresholds':thresholds_auc}).to_csv(ROC_savepath / f'train_ROC_for_{fold}th_fold_{i}th_attempt.csv')
                # pd.DataFrame.from_dict({'fprs':fprs_val, 'tprs':tprs_val, 'thresholds':thresholds_auc_val}).to_csv(ROC_savepath / f'val_ROC_for_{fold}th_fold_{i}th_attempt.csv')
                pd.DataFrame.from_dict({'fprs':fprs_test, 'tprs':tprs_test, 'thresholds':thresholds_auc_test}).to_csv(ROC_savepath / f'test_ROC_for_{fold}th_fold_{i}th_attempt.csv')

                # pd.DataFrame.from_dict({'pres':pres, 'recs':recs, 'thresholds':thresholds_prc}).to_csv(PRC_savepath / f'train_PRC_for_{fold}th_fold_{i}th_attempt.csv')
                # pd.DataFrame.from_dict({'pres':pres_val, 'recs':recs_val, 'thresholds':thresholds_prc_val}).to_csv(PRC_savepath / f'val_PRC_for_{fold}th_fold_{i}th_attempt.csv')
                pd.DataFrame.from_dict({'pres':pres_test, 'recs':recs_test, 'thresholds':thresholds_prc_test}).to_csv(PRC_savepath / f'test_PRC_for_{fold}th_fold_{i}th_attempt.csv')

                # pd.DataFrame.from_dict({'y_true':trainY, 'y_pred':trainY_pred}).to_csv(logits_savepath / f'train_logits_{fold}th_fold_{i}th_attempt.csv')
                # pd.DataFrame.from_dict({'y_true':validY, 'y_pred':validY_pred}).to_csv(logits_savepath / f'val_logits_{fold}th_fold_{i}th_attempt.csv')
                pd.DataFrame.from_dict({'y_true':testY, 'y_pred':testY_pred}).to_csv(logits_savepath / f'test_logits_{fold}th_fold_{i}th_attempt.csv')
            if uc_flag:
                # Aggregate stochastic predictions when Monte Carlo uncertainty is requested.
                from MC_UQ_analyse import mc_uq_calculate
                uncertainty_results = mc_uq_calculate(logits_savepath, logits_savepath / f'test_logits_{fold}th_fold_mean4total{self.params.uncertainty_iteration_size}iterations.csv')

            allfold_train_data[fold] = kfold_train_data
            allfold_val_data[fold] = kfold_val_data
            allfold_test_data[fold] = kfold_test_data
            # evaluate finished

        ave_tr, ave_v, ave_te = self.save(assess, allfold_train_data, allfold_val_data, allfold_test_data)
        print(f'-------------------------------- {self.params.kfold_num} folds average result --------------------------------')
        print(f'AVERAGE_TRAIN result: acc = {ave_tr.acc:.4f}; auc = {ave_tr.auc:.4f}, mcc = {ave_tr.mcc:.4f}, precision = {ave_tr.precision:.4f}, recall = {ave_tr.recall:.4f}, specificity = {ave_tr.specificity:.4f}, sensitivity = {ave_tr.sensitivity:.4f}, f1 = {ave_tr.f1:.4f}, prauc = {ave_tr.prauc:.4f}, av_prc = {ave_tr.av_prc:.4f}')
        print(f'AVERAGE_VALID result: acc = {ave_v.acc:.4f}; auc = {ave_v.auc:.4f}, mcc = {ave_v.mcc:.4f}, precision = {ave_v.precision:.4f}, recall = {ave_v.recall:.4f}, specificity = {ave_v.specificity:.4f}, sensitivity = {ave_v.sensitivity:.4f}, f1 = {ave_v.f1:.4f}, prauc = {ave_v.prauc:.4f}, av_prc = {ave_v.av_prc:.4f}')
        print(f'AVERAGE_TEST result: acc = {ave_te.acc:.4f}; auc = {ave_te.auc:.4f}, mcc = {ave_te.mcc:.4f}, precision = {ave_te.precision:.4f}, recall = {ave_te.recall:.4f}, specificity = {ave_te.specificity:.4f}, sensitivity = {ave_te.sensitivity:.4f}, f1 = {ave_te.f1:.4f}, prauc = {ave_te.prauc:.4f}, av_prc = {ave_te.av_prc:.4f}')

    def save(self, assess, allfold_train_data, allfold_val_data, allfold_test_data):
        # Export per-fold summaries first, then compute fold-wise averages for the final report.
        for fold in allfold_train_data.keys():
            train_result = pd.DataFrame({})
            for ass in assess:
                train_result[ass] = allfold_train_data[fold][ass]
            train_result['Iteration'] = range(1, self.params.uncertainty_iteration_size+1)
            train_result['Foldid'] = fold
            
            val_result = pd.DataFrame({})
            for ass in assess:
                val_result[ass] = allfold_val_data[fold][ass]
            val_result['Iteration'] = range(1, self.params.uncertainty_iteration_size+1)
            val_result['Foldid'] = fold
            
            test_result = pd.DataFrame({})
            for ass in assess:
                test_result[ass] = allfold_test_data[fold][ass]
            test_result['Iteration'] = range(1, self.params.uncertainty_iteration_size+1)
            test_result['Foldid'] = fold
            
            resultdata_savepath = self.save_path / 'result_data' / f'{fold}th_fold'
            resultdata_savepath.mkdir(parents=True, exist_ok=True)
            train_result.to_csv(resultdata_savepath / f'train_result_{fold}th_fold.csv')
            val_result.to_csv(resultdata_savepath / f'val_result_{fold}th_fold.csv')
            test_result.to_csv(resultdata_savepath / f'test_result_{fold}th_fold.csv')
        
        result_tr = pd.DataFrame([])
        result_v = pd.DataFrame([])
        result_te = pd.DataFrame([])
        for ass in assess:
            for k in range(self.params.kfold_num):
                result_tr.at[k,'Foldid'] = k
                result_v.at[k,'Foldid'] = k
                result_te.at[k,'Foldid'] = k
                result_tr.at[k,ass] = np.mean(allfold_train_data[k][ass])
                result_v.at[k,ass] = np.mean(allfold_val_data[k][ass])
                result_te.at[k,ass] = np.mean(allfold_test_data[k][ass])
                
        resultdata_savepath = self.save_path / 'result_data'
        resultdata_savepath.mkdir(parents=True, exist_ok=True)
        result_tr.to_csv(resultdata_savepath / f'train_result_mean4total{self.params.uncertainty_iteration_size}iterations.csv')
        result_v.to_csv(resultdata_savepath / f'val_result_mean4total{self.params.uncertainty_iteration_size}iterations.csv')
        result_te.to_csv(resultdata_savepath / f'test_result_mean4total{self.params.uncertainty_iteration_size}iterations.csv')

        return result_tr.mean(), result_v.mean(), result_te.mean()
