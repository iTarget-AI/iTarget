import pandas as pd

infer = pd.read_csv('infer_prot-drug.csv')
prots = infer[['protid', 'sequence']].drop_duplicates()
drugs = infer[['drugid', 'smiles']].drop_duplicates()
protdrug_index = infer[['source','drugid','protid','label']]

protdrug_index.to_csv('../infer_prot-drug-index-source.csv', index=False)
prots.to_csv('../infer_prots.csv', index=False)
drugs.to_csv('../infer_drugs.csv', index=False)
