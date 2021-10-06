# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 13:03:13 2020
@author: Ruibo
Class dataset
"""
import os
import numpy as np
import pandas as pd
import pickle
import Preprocessing


class DataSet(object):


    def __init__(self, training_dir):
        """
        Remember to change is_training to False,
            if this is used in SubChallenge 1 submitting or testing.
            It will then skip loading aucs.
        """
        assert isinstance(training_dir, str)
        self.training_dir = training_dir
        self.is_loaded = False
        self.is_prepared = False
        self.is_splited = False
        self.mutate_RNAs = True
        self.is_training = True
        self.is_sc2 = True

    def load_from_csv(self):
        def getCsv(fname):
            print(f'Loading {fname} data...', flush=True)
            return pd.read_csv(os.path.join(self.training_dir, fname + '.csv'))
        self.rnaseq = getCsv('rnaseq')
        self.dnaseq = getCsv('dnaseq')
        self.clinical_numerical = getCsv('clinical_numerical').\
            set_index('lab_id').astype('float64')
        self.clinical_categorical = getCsv('clinical_categorical').\
            set_index('lab_id')
        self.clinical_categorical_legend = getCsv('clinical_categorical_legend')
        # Selective load. May not have these in test cases.
        try:
            self.aucs = getCsv('aucs')
        except FileNotFoundError:
            self.is_training = False
        try:
            self.responses = getCsv('response').set_index('lab_id')
        except FileNotFoundError:
            self.is_sc2 = False

        # Replace mutated RNA expressions
        if self.mutate_RNAs:
            self._replace_muted(None)
        # Transpose and normailize RNA
        self.rnaseq = _NormSpecimens(_TransposeRnaSeqTable(self.rnaseq))
        # One-hot coding of clinical catagory
        self.clinical_categorical = Preprocessing.expand_multiclass_to_onehot(
            self.clinical_categorical, sele_cols=None)
        if self.is_training:
            # Prepare aucs, drop off too many nans, and impute the rest
            self.aucs = self.aucs.pivot(
                    index='lab_id', columns='inhibitor', values='auc')
            self.aucs = Preprocessing.drop_too_many_nans(self.aucs, 0.5)
            self.aucs = Preprocessing.Impute(self.aucs)
            # Reindex accrd to aucs
            idx = self.aucs.index
        else:
            idx = self.rnaseq.index
        self._reindex(idx)
        self.lab_id = idx
        self.is_loaded = True
        return self

    def append_set(self, leaderboard_dir):
        assert self.is_loaded
        assert isinstance(leaderboard_dir, str)
        set2 = DataSet(leaderboard_dir)
        set2.load_from_csv()
        # by default, pd.concat do axis=0, jion=outer
        self.rnaseq = pd.concat((self.rnaseq, set2.rnaseq))
        self.dnaseq = pd.concat((self.dnaseq, set2.dnaseq))
        self.clinical_numerical = pd.concat((self.clinical_numerical,
                                             set2.clinical_numerical))
        self.clinical_categorical = pd.concat((self.clinical_categorical,
                                               set2.clinical_categorical), sort=True)
        if self.is_training:
            self.aucs = pd.concat((self.aucs, set2.aucs))
            idx = self.aucs.index
        else:
            idx = self.rnaseq.index
        self._reindex(idx)
        self.lab_id = idx
        return self

    def save_to_pickle(self, fname='DataSet'):
        """ The pickle file will be saved to training folder. (ist set folder) """
        with open(os.path.join(self.training_dir, fname + '.pickle'), 'wb') as file:
            pickle.dump(self, file)
        return None
#%%
    def _return_item(self, item):
        assert isinstance(item, str)
        if 'rna' in item.lower():
            return self.rnaseq
        elif 'cat' in item.lower():
            return self.clinical_categorical
        elif 'num' in item.lower():
            return self.clinical_numerical
        elif 'auc' in item.lower():
            return self.aucs
        elif 'resp' in item.lower():
            return self.responses

    def return_temp(self, *items):
        assert self.is_loaded
        rt = []
        for each_one in items:
            rt.append(self._return_item(each_one))
        return pd.concat(rt, axis=1, join='inner')

    def return_xy(self, *xy_items):
        """ Only the last is Y, OK? """
        if len(xy_items)==0:
            x_items = ('rna', 'categorical', 'numerical')
            print("dataset.return_xy, did not specify x. Returnning all the 3.")
            y_items = ('aucs')
            print("dataset.return_xy, did not specify y. Returning aucs.")
        elif len(xy_items)==1:
            raise ValueError("Wrong, passing only one item for x and y.")
        else:
            y_items = xy_items[-1]
            x_items = xy_items[:-1]
        self.xx = self.return_temp(*x_items).reindex(self.lab_id)
        self.yy = self.return_temp(y_items).reindex(self.lab_id)
        self.is_prepared = True
        return (self.xx, self.yy)

    def kf_split(self):
        assert self.is_prepared
        from sklearn.model_selection import KFold
        kf = KFold(n_splits=5, shuffle=True, random_state=2020)
        kFolds = {}
        for i, (train_iloc, test_iloc) in enumerate(kf.split(self.lab_id)):
            train_idx = self.lab_id[train_iloc]
            test_idx = self.lab_id[test_iloc]
            kFolds['fold{}'.format(i+1)] = {'train_idx': train_idx,
                                            'test_idx': test_idx,
                                            'sele_fts': None}
        self.kfolds = kFolds
        self.is_splited = True
        return None

    def __getitem__(self, idx):
        assert self.is_splited
        if idx==len(self.kfolds):
            raise StopIteration()
        train_idx = self.kfolds['fold{}'.format(idx+1)]['train_idx']
        test_idx = self.kfolds['fold{}'.format(idx+1)]['test_idx']
        xx_train = self.xx.loc[train_idx]
        yy_train = self.yy.loc[train_idx]
        xx_test = self.xx.loc[test_idx]
        yy_test = self.yy.loc[test_idx]
        return (xx_train, yy_train, xx_test, yy_test)
        
    def _replace_muted(self, replace_value=None):
        """
        Find the mutated genes from DNA seq
        Replace the mutated RNA expressions to negative of origin or -5, or sth else.
        """
        dna_raw = self.dnaseq[['lab_id','Hugo_Symbol']]
        rna_raw = self.rnaseq.set_index('Symbol')
        for _, each_muta in dna_raw.iterrows():
            tt = rna_raw.loc[each_muta['Hugo_Symbol'], each_muta['lab_id']]
            if replace_value is None:
                if isinstance(tt, float):
                    rna_raw.loc[each_muta['Hugo_Symbol'], each_muta['lab_id']] = -tt
                else:
                    rna_raw.loc[each_muta['Hugo_Symbol'], each_muta['lab_id']] = -tt.values[0]
            else:
                rna_raw.loc[each_muta['Hugo_Symbol'], each_muta['lab_id']] = replace_value
        self.rnaseq = rna_raw.reset_index()
        return None

    def _reindex(self, index):
        """Reindex all the X and Ys. If is training, use aucs, if not, use rnas"""
        self.rnaseq = self.rnaseq.reindex(index)
        self.clinical_categorical = self.clinical_categorical.reindex(index)
        self.clinical_numerical = self.clinical_numerical.reindex(index)
        if self.is_training:
            self.aucs = self.aucs.reindex(index)
        if self.is_sc2:
            self.responses = self.responses.reindex(index)

    def _check_loaded(self):
        # to be done....
        assert(self.aucs.index == self.rnaseq.index).all()

def _TransposeRnaSeqTable(rnaseq):
        """Convert the RnaSeq table, indexed by gene, to be indexed by specimen."""
        rnaseq.index = rnaseq.Gene
        return rnaseq[rnaseq.columns[2:]].T

def _NormSpecimens(specimens):
    normed_specimens = specimens.apply(
            lambda specimen : specimen / np.linalg.norm(specimen), axis=1)
    return normed_specimens


if __name__ == '__main__':
    # Testing codes
    dataset_obj = DataSet('training')
    dataset_obj.load_from_csv()
    dataset_obj.append_set('input')
    xxx, yyy = dataset_obj.return_xy('rna','aucs')
    dataset_obj.kf_split()
    dataset_obj.save_to_pickle()

    with open('./training/DataSet.pickle', 'rb') as file:
        dataset_copy = pickle.load(file)
    for each_fold in dataset_copy:
        (x1,y1,x2,y2) = each_fold
