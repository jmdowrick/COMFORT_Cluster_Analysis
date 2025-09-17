# Creates the data file to be used for factor analysis on the questionnaire data from the COMFORT and PSYKI cohorts
import pandas as pd
import yaml
import os

with open('config.yml') as file:
    config = yaml.safe_load(file)

# Import COMFORT questionnaire data
path = config['default']['data_raw'] + config['default']['file_PRO']
HADS_sheet = "HADS"
df_hads = pd.read_excel(io=path, sheet_name=HADS_sheet)
df_hads = df_hads[['Participant #', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E14']]
df_hads.rename(columns={'Participant #': 'id'}, inplace=True)

SAGIS_sheet = "SAGIS"
df_sagis = pd.read_excel(io=path, sheet_name=SAGIS_sheet)
df_sagis.drop(['Participant #.1', 'CASE/Control', 'CASE/CONTROL', 'Epigastric', 'IBS-D', 'Acid', 'Vomiting', 'Constipation'], axis=1, inplace=True)
dict_SAGIS = {('1. Belching with acid taste/heartburn/burning sensation in the oesophagus '
               '(tube joining mouth and stomach)'): 'sagis_1',
              '2. Dysphagia (difficulty swallowing)': 'sagis_2',
              '3. Fullness (feeling of congestion of food without relation to prior food intake)': 'sagis_3',
              ('4. Early satiety (stomach is overfilled soon after starting to eat, disproportional to the '
               'quantity of food taken, so that food cannot be finished)'): 'sagis_4',
              '5. Postprandial pain or discomfort (upper abdominal symptoms start or get worse after meals)': 'sagis_5',
              ('6. Epigastric pain/upper abdominal pain (pain between the belly button and chest/pain noticeable '
               'in the upper abdomen)'): 'sagis_6',
              '7. Retrosternal discomfort (unpleasant feeling behind the middle of the chest, painful or drawing)':
                  'sagis_7',
              '8. Pain or discomfort prior to bowel movement': 'sagis_8',
              '9. Difficulty with emptying the bowel': 'sagis_9',
              '10. Constipation': 'sagis_10',
              '11. Loose stools': 'sagis_11',
              '12. Incontinence': 'sagis_12',
              '13. Urgency to empty the bowel': 'sagis_13',
              '14. Diarrhoea': 'sagis_14',
              '15. Loss of appetite (listless for food intake)': 'sagis_15',
              '16. Abdominal cramps (spasmodic or colic like stomach pain without specified localisation)': 'sagis_16',
              '17. Sickness (discomfort combined with the impression for the need to vomit)': 'sagis_17',
              '18. Nausea (urgent feeling of the need to vomit)': 'sagis_18',
              '19. Vomiting (vomiting of mucus and gastric contents, or strong unproductive retching)': 'sagis_19',
              '20. Bloating (feeling of distension and excessive gas in the abdomen)': 'sagis_20',
              '21. Excessive gas and passing of wind': 'sagis_21',
              '22. Excessive belching': 'sagis_22',
              'Participant #': 'id'}
df_sagis.rename(columns=dict_SAGIS, inplace=True)

PROMIS_sheet = "PROMIS"
df_promis = pd.read_excel(io=path, sheet_name=PROMIS_sheet)
df_promis.drop(['CASE/CONTROL', 'CASE/Control', 'GENDER', 'AGE'], axis=1, inplace=True)
dict_PROMIS = {'Participant #': 'id',
               'P_Anxiety': 'pr_anx',
               'P_Bloating': 'pr_blo',
               'P_Constipation': 'pr_con',
               'P_DS': 'pr_ds',
               'P_Depression': 'pr_dep',
               'P_Diarrhoea': 'pr_dia',
               'P_Pain': 'pr_pain',
               'P_Reflux': 'pr_ref'}
df_promis.rename(columns=dict_PROMIS, inplace=True)

# Combine all data
df = df_hads.merge(df_promis, how='outer', on='id')
df = df.merge(df_sagis, how='outer', on='id')
              
# Replace missing values with -999
df = df.fillna(-999)

#  Save output
output_path = config['default']['data_mplus']+"comfort_PROs.csv"
df.to_csv(output_path, index=False)

output_path_names = os.path.normpath(config['default']['data_mplus']+"item_names.txt")
with open(output_path_names, 'w') as fp:
    fp.write(' '.join(df.columns.tolist()))
